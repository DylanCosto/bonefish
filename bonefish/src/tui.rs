use std::fmt;
use std::io::{self, IsTerminal};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use console::{Term, style};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressState, ProgressStyle};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Stage {
    Reference,
    Align,
    Merge,
    Bam,
    Index,
}

impl Stage {
    fn name(self) -> &'static str {
        match self {
            Self::Reference => "Prepare reference",
            Self::Align => "Align reads",
            Self::Merge => "Merge index parts",
            Self::Bam => "Sort and write BAM",
            Self::Index => "Write BAM index",
        }
    }
}

/// What the reference stage is doing, shown while it runs.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum RefWork {
    Read,
    Build,
    Load,
}

impl RefWork {
    fn verb(self) -> &'static str {
        match self {
            Self::Read => "reading",
            Self::Build => "building index",
            Self::Load => "loading index",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Outcome {
    Complete,
    Failed,
    Interrupted,
}

struct Row {
    stage: Stage,
    elapsed: Duration,
    visited: bool,
}

#[derive(Clone, Copy)]
struct Span {
    part: usize,
    end: u64,
    total: u64,
    bases: Option<(u64, u64)>,
}

struct Run {
    started: Instant,
    since: Instant,
    finished: Option<(Instant, Outcome)>,
    /// The row whose clock is running.
    active: usize,
    /// The row marked as current. It never moves back up: indexing a later reference part runs
    /// the "Prepare reference" clock while "Align reads" stays marked.
    focus: usize,
    rows: Vec<Row>,
    /// None until the reference work starts (option parsing, opening the inputs).
    ref_work: Option<RefWork>,
    /// Reference part being prepared (1-based; 0 before the first), and the last one aligned
    /// against. Reading past the last part briefly starts one more that never aligns.
    part: usize,
    aligned_part: usize,
    /// Number of reference parts, when known up front.
    parts: Option<usize>,
    /// How far into the reference (or index) file the parts so far reached, to estimate the
    /// number of parts when it is not known up front.
    span: Option<Span>,
    /// The largest fraction shown so far: a revised estimate never moves the bar back.
    shown: f64,
    /// Reads aligned against the current part, and over all parts.
    reads: usize,
    reads_total: usize,
    /// Size of all query files (None: stdin, a pipe, or no queries).
    input_total: Option<u64>,
    /// Bytes of finished query files in this part, and of the current file.
    input_done: u64,
    input_pos: u64,
    bam: Option<(u64, u64)>,
    align_start: Option<Instant>,
}

/// Where `--progress` writes how far the run has got, and when it last did.
struct ProgressFile {
    path: PathBuf,
    last: Option<Instant>,
    /// The largest fraction written: a revised estimate never moves it back.
    shown: f64,
}

/// How often the progress file is written at most, but at a stage's start or end.
const PROGRESS_INTERVAL: Duration = Duration::from_millis(500);

/// What preparing the reference is doing, for `--progress`.
fn reference_work(work: Option<RefWork>) -> String {
    match work {
        Some(RefWork::Read) => "reading the reference",
        Some(RefWork::Build) => "building the reference index",
        Some(RefWork::Load) => "loading the reference index",
        None => "preparing the reference",
    }
    .to_owned()
}

impl Run {
    /// How far the whole run has got (0 to 1), and what it is doing, for `--progress`: preparing
    /// the reference up to 2%, aligning to 85%, merging parts, writing the sorted BAM to 98% and
    /// indexing it.
    fn overall(&self, now: Instant) -> (f64, String) {
        if let Some((_, outcome)) = self.finished {
            return match outcome {
                Outcome::Complete => (1.0, format!("done · {} reads", count(self.reads_total))),
                Outcome::Failed => (0.0, format!("{} failed", self.rows[self.active].stage.name().to_lowercase())),
                Outcome::Interrupted => (0.0, "interrupted".to_owned()),
            };
        }
        match self.rows[self.focus].stage {
            Stage::Reference => (0.0, reference_work(self.ref_work)),
            Stage::Align => {
                let fraction = self.align_fraction().unwrap_or(0.0);
                let mut what = if self.active == self.focus {
                    format!("aligning · {} reads", count(self.reads_total))
                } else {
                    reference_work(self.ref_work)
                };
                if let Some(part) = self.part_label() { what.push_str(&format!(" · {part}")); }
                if let Some(started) = self.align_start {
                    let t = (now - started).as_secs_f64();
                    if (0.005..1.0).contains(&fraction) && t >= 5.0 {
                        what.push_str(&format!(" · ~{} left", short_time(t * (1.0 - fraction) / fraction)));
                    }
                }
                (0.02 + 0.83 * fraction, what)
            }
            Stage::Merge => (0.86, "merging the index parts' alignments".to_owned()),
            Stage::Bam => match self.bam.filter(|&(_, total)| total > 0) {
                Some((done, total)) => {
                    (0.87 + 0.11 * done as f64 / total as f64, format!("writing the sorted BAM · {done} of {total}"))
                }
                None => (0.87, "sorting the alignments".to_owned()),
            },
            Stage::Index => (0.98, "indexing the BAM".to_owned()),
        }
    }

    /// Fraction of all alignment work done, over all parts when their number is known.
    fn align_fraction(&self) -> Option<f64> {
        let total = self.input_total.filter(|&t| t > 0)? as f64;
        let in_part = ((self.input_done + self.input_pos) as f64 / total).min(1.0);
        // Between parts, the last aligned part is finished.
        let in_part = if self.active == self.focus { in_part } else { 1.0 };
        // Every part is a full pass over the reads, so parts weigh the same.
        let parts = self.parts_estimate().map_or(1, |(parts, _)| parts).max(self.aligned_part).max(1);
        let fraction = (self.aligned_part.saturating_sub(1) as f64 + in_part) / parts as f64;
        Some(fraction.min(1.0))
    }

    fn row_elapsed(&self, index: usize, now: Instant) -> Duration {
        let row = &self.rows[index];
        if index == self.active && self.finished.is_none() { row.elapsed + (now - self.since) } else { row.elapsed }
    }

    /// The number of reference parts, and whether it is exact: known up front, or estimated from
    /// how much of the reference file the parts so far took.
    fn parts_estimate(&self) -> Option<(usize, bool)> {
        if let Some(parts) = self.parts {
            return Some((parts, true));
        }
        let span = self.span?;
        if span.end >= span.total {
            return Some((span.part, true));
        }
        let left = (span.total - span.end) as f64;
        let more = match span.bases {
            // Each part holds at least `batch` bases, so this errs towards too many parts: the
            // bar then jumps forward rather than back.
            Some((bases, batch)) if bases > 0 && batch > 0 => {
                (left * bases as f64 / span.end.max(1) as f64 / batch as f64).ceil() as usize
            }
            _ => (left / (span.end as f64 / span.part.max(1) as f64).max(1.0)).ceil() as usize,
        };
        Some((span.part + more.max(1), false))
    }

    fn part_label(&self) -> Option<String> {
        match self.parts_estimate() {
            Some((parts, exact)) if parts > 1 => {
                let approx = if exact { "" } else { "~" };
                Some(format!("part {}/{approx}{}", self.part.min(parts), parts))
            }
            None if self.part > 1 => Some(format!("part {}", self.part)),
            _ => None,
        }
    }

    /// The live detail of the "Align reads" row, fitted to `width` characters.
    fn align_detail(&self, now: Instant, width: usize) -> String {
        let fraction = self.align_fraction().map(|f| f.max(self.shown));
        let mut pieces: Vec<String> = Vec::new();
        if self.active != self.focus {
            // Between parts: the next part's index is being prepared.
            let mut s = self.ref_work.map_or("working", RefWork::verb).to_owned();
            if let Some(part) = self.part_label() { s = format!("{s} · {part}"); }
            pieces.push(s);
        } else {
            pieces.push(format!("{} reads", count(self.reads)));
            if let Some(part) = self.part_label() { pieces.push(part); }
        }
        let align_time = self.row_elapsed(self.focus, now).as_secs_f64();
        let rate = (align_time >= 2.0 && self.reads_total > 0).then(|| format!("{}/s", count((self.reads_total as f64 / align_time) as usize)));
        // Time since alignment began, including indexing of later parts, for the estimate.
        let started = self.align_start.map(|t| (now - t).as_secs_f64());
        let eta = match (fraction, started) {
            (Some(f), Some(t)) if (0.005..1.0).contains(&f) && t >= 5.0 => Some(format!("~{} left", short_time(t * (1.0 - f) / f))),
            _ => None,
        };
        // Most useful first; keep what fits.
        let bar = fraction.map(|f| {
            let filled = (f * 10.0) as usize;
            format!("[{}{}] {:>3}%", "█".repeat(filled), "░".repeat(10 - filled.min(10)), (f * 100.0) as usize)
        });
        let mut out = String::new();
        let candidates = bar.into_iter().chain(pieces).chain(eta).chain(rate);
        for piece in candidates {
            let extra = if out.is_empty() { piece.chars().count() } else { piece.chars().count() + 3 };
            if out.chars().count() + extra > width { continue; }
            if !out.is_empty() { out.push_str(" · "); }
            out.push_str(&piece);
        }
        out
    }

    fn draw(&mut self, out: &mut dyn fmt::Write, width: usize) -> fmt::Result {
        if let Some(f) = self.align_fraction() { self.shown = self.shown.max(f); }
        let (now, status) = match self.finished {
            None => (Instant::now(), "Running"),
            Some((end, Outcome::Complete)) => (end, "Complete"),
            Some((end, Outcome::Failed)) => (end, "Failed"),
            Some((end, Outcome::Interrupted)) => (end, "Interrupted"),
        };
        let title = format!("bonefish  {status}  {}", time(now - self.started));
        let title = title
            .chars()
            .take(width.saturating_sub(1))
            .collect::<String>();
        writeln!(out, "{}", style(title).for_stderr().bold())?;

        // Keep the longest line short enough to avoid wrapping in small terminals.
        if width < 56 {
            if let Some((_, outcome)) = self.finished {
                // One summary line: what was done, or where it stopped.
                let line = match outcome {
                    Outcome::Complete => {
                        let mut s = format!("Done · {} reads", count(self.reads));
                        if let Some((done, _)) = self.bam { s.push_str(&format!(" · {done} records")); }
                        s
                    }
                    Outcome::Failed => format!("{} failed", self.rows[self.active].stage.name()),
                    Outcome::Interrupted => format!("{} interrupted", self.rows[self.active].stage.name()),
                };
                return write!(out, "{}", line.chars().take(width.saturating_sub(1)).collect::<String>());
            }
            let row = &self.rows[self.focus];
            let detail = match row.stage {
                Stage::Align => format!(" · {}", self.align_detail(now, width.saturating_sub(row.stage.name().len() + 4))),
                Stage::Bam => self
                    .bam
                    .filter(|&(_, total)| total > 0)
                    .map_or_else(String::new, |(done, total)| {
                        format!(" · {}%", done * 100 / total)
                    }),
                _ => String::new(),
            };
            return write!(
                out,
                "{}",
                format!("{}{detail}", row.stage.name())
                    .chars()
                    .take(width.saturating_sub(1))
                    .collect::<String>()
            );
        }

        let body = width.min(90) - 3;
        let details = body - 32;
        let label = "─ Alignment progress ";
        writeln!(out, "┌{label}{}┐", "─".repeat(body - label.chars().count()))?;
        let stopped = match self.finished {
            Some((_, Outcome::Failed)) => Some("failed"),
            Some((_, Outcome::Interrupted)) => Some("interrupted"),
            _ => None,
        };
        for (index, row) in self.rows.iter().enumerate() {
            let current = index == self.focus && self.finished.is_none();
            let (mark, detail) = if current {
                let detail = match row.stage {
                    Stage::Align => self.align_detail(now, details),
                    Stage::Reference => self.ref_work.map_or("working", RefWork::verb).to_owned(),
                    Stage::Bam => self.bam.filter(|&(_, total)| total > 0).map_or_else(
                        || "working".to_owned(),
                        |(done, total)| {
                            let filled = (done * 12 / total) as usize;
                            format!(
                                "[{}{}] {:>3}%",
                                "█".repeat(filled),
                                "░".repeat(12 - filled),
                                done * 100 / total
                            )
                        },
                    ),
                    _ => "working".to_owned(),
                };
                ("›", detail)
            } else if let Some(stopped) = stopped.filter(|_| index == self.active) {
                (" ", stopped.to_owned())
            } else if row.visited {
                let detail = match (row.stage, self.bam) {
                    (Stage::Align, _) if self.finished.is_some() || index < self.focus => {
                        // Drop the part count rather than cut it short.
                        let reads = format!("done · {} reads", count(self.reads));
                        let parts = format!("{reads} · {} parts", self.aligned_part);
                        if self.aligned_part > 1 && parts.chars().count() <= details { parts } else { reads }
                    }
                    (Stage::Bam, Some((done, _))) => format!("done · {done} records"),
                    _ => "done".to_owned(),
                };
                (" ", detail)
            } else {
                (" ", "queued".to_owned())
            };
            let clock = if row.visited {
                time(self.row_elapsed(index, now))
            } else {
                String::new()
            };
            let line = format!(
                "{mark} {:<19} {:<details$.details$} {:>8} ",
                row.stage.name(),
                detail,
                clock
            );
            if current {
                writeln!(out, "│{}│", style(line).for_stderr().cyan().bold())?;
            } else {
                writeln!(out, "│{line}│")?;
            }
        }
        write!(out, "└{}┘", "─".repeat(body))
    }
}

fn time(duration: Duration) -> String {
    let s = duration.as_secs();
    if s >= 3600 {
        format!("{}:{:02}:{:02}", s / 3600, s / 60 % 60, s % 60)
    } else {
        format!("{:02}:{:02}", s / 60, s % 60)
    }
}

/// A rough duration: "40s", "12m", "3h20m".
fn short_time(secs: f64) -> String {
    let s = secs.round() as u64;
    if s < 60 {
        format!("{s}s")
    } else if s < 3600 {
        format!("{}m", s.div_ceil(60))
    } else {
        format!("{}h{:02}m", s / 3600, s / 60 % 60)
    }
}

/// A count with a k/M suffix above 10,000.
fn count(n: usize) -> String {
    match n {
        0..10_000 => n.to_string(),
        10_000..10_000_000 => format!("{:.1}k", n as f64 / 1e3),
        _ => format!("{:.1}M", n as f64 / 1e6),
    }
}

#[derive(Clone)]
pub(crate) struct Tui {
    run: Arc<Mutex<Run>>,
    /// The panel on stderr, with `--tui` on a terminal.
    bar: Option<ProgressBar>,
    /// `--progress`'s file.
    progress: Option<Arc<Mutex<ProgressFile>>>,
}

impl Tui {
    /// The run's progress: drawn as a panel on stderr when `requested` and stderr is a terminal
    /// that alignment records are not printed to, and written to `progress` when given. None
    /// when neither.
    pub(crate) fn new(
        requested: bool,
        progress: Option<PathBuf>,
        stdout_data: bool,
        align: bool,
        merge: bool,
        bam: bool,
        input_total: Option<u64>,
    ) -> Option<Self> {
        let target = ProgressDrawTarget::stderr_with_hz(8);
        let display = requested
            && !(stdout_data && io::stdout().is_terminal())
            && !std::env::var_os("TERM").is_some_and(|term| term == "dumb")
            && !target.is_hidden();
        if !display && progress.is_none() {
            return None;
        }
        let mut stages = vec![Stage::Reference];
        if align {
            stages.push(Stage::Align);
        }
        if merge {
            stages.push(Stage::Merge);
        }
        if bam {
            stages.extend([Stage::Bam, Stage::Index]);
        }
        let now = Instant::now();
        let run = Arc::new(Mutex::new(Run {
            started: now,
            since: now,
            finished: None,
            active: 0,
            focus: 0,
            rows: stages
                .into_iter()
                .map(|stage| Row {
                    stage,
                    elapsed: Duration::ZERO,
                    visited: stage == Stage::Reference,
                })
                .collect(),
            ref_work: None,
            part: 0,
            aligned_part: 0,
            parts: None,
            span: None,
            shown: 0.0,
            reads: 0,
            reads_total: 0,
            input_total,
            input_done: 0,
            input_pos: 0,
            bam: None,
            align_start: None,
        }));
        let bar = display.then(|| {
            let snapshot = Arc::clone(&run);
            let theme = ProgressStyle::with_template("{panel}")
                .expect("valid panel template")
                .with_key(
                    "panel",
                    move |_: &ProgressState, out: &mut dyn fmt::Write| {
                        let (_, width) = Term::stderr().size();
                        let _ = snapshot.lock().unwrap().draw(out, usize::from(width));
                    },
                );
            ProgressBar::with_draw_target(None, target).with_style(theme)
        });
        let progress = progress.map(|path| Arc::new(Mutex::new(ProgressFile { path, last: None, shown: 0.0 })));
        let tui = Self { run, bar, progress };
        tui.report(true);
        Some(tui)
    }

    /// Whether the panel is drawn: without it, status lines go to stderr as without `--tui`.
    pub(crate) fn displayed(&self) -> bool {
        self.bar.is_some()
    }

    /// Writes the progress file, at most every [`PROGRESS_INTERVAL`] unless `now` (a stage's
    /// start or end). It is replaced by a rename, so a reader never sees half a line: one line,
    /// the fraction done and what is being done, tab-separated (`0.4130\taligning · 2.4M
    /// reads`). A file that cannot be written is left: progress is not worth failing a run for.
    fn report(&self, now: bool) {
        let Some(progress) = &self.progress else { return };
        let mut file = progress.lock().unwrap_or_else(|e| e.into_inner());
        if !now && file.last.is_some_and(|t| t.elapsed() < PROGRESS_INTERVAL) {
            return;
        }
        file.last = Some(Instant::now());
        let (fraction, what, complete) = {
            let run = self.run.lock().unwrap_or_else(|e| e.into_inner());
            let (fraction, what) = run.overall(Instant::now());
            (fraction, what, matches!(run.finished, Some((_, Outcome::Complete))))
        };
        // A failed or interrupted run keeps the fraction it reached.
        let fraction = if complete { 1.0 } else { fraction.max(file.shown) };
        file.shown = fraction;
        let line = format!("{:.4}\t{}\n", fraction.clamp(0.0, 1.0), what.replace(['\t', '\n'], " "));
        let tmp = file.path.with_extension("tmp");
        if std::fs::write(&tmp, line).is_ok() {
            let _ = std::fs::rename(&tmp, &file.path);
        }
    }

    pub(crate) fn start(&self) {
        if let Some(bar) = &self.bar {
            bar.tick();
            bar.enable_steady_tick(Duration::from_millis(200));
        }
    }

    /// Starts preparing reference part `part` (1-based).
    pub(crate) fn reference(&self, part: usize, work: RefWork) {
        {
            let mut run = self.run.lock().unwrap();
            run.ref_work = Some(work);
            run.part = part;
        }
        self.stage(Stage::Reference);
    }

    /// Reference parts 1..=`part` came from the first `end` bytes of a `total`-byte file;
    /// `bases`: (their bases, the -I batch size), when known.
    pub(crate) fn reference_span(&self, part: usize, end: u64, total: u64, bases: Option<(u64, u64)>) {
        self.run.lock().unwrap().span = Some(Span { part, end, total, bases });
    }

    /// The number of reference parts, when known before they are all read.
    pub(crate) fn parts(&self, parts: usize) {
        self.run.lock().unwrap().parts = Some(parts);
    }

    /// Starts aligning one query file against the current part.
    pub(crate) fn align_file(&self) {
        {
            let mut run = self.run.lock().unwrap();
            if run.aligned_part == run.part {
                // Another query file against the same part.
                run.input_done += run.input_pos;
            } else {
                run.aligned_part = run.part;
                run.reads = 0;
                run.input_done = 0;
            }
            run.input_pos = 0;
        }
        self.stage(Stage::Align);
    }

    pub(crate) fn stage(&self, stage: Stage) {
        let mut run = self.run.lock().unwrap();
        let next = run.rows.iter().position(|row| row.stage == stage).or_else(|| {
            // A merge that starts without --split-prefix (multi-part SAM/BAM): add its row
            // before the BAM stages.
            (stage == Stage::Merge).then(|| {
                let at = run.rows.iter().position(|row| matches!(row.stage, Stage::Bam | Stage::Index)).unwrap_or(run.rows.len());
                run.rows.insert(at, Row { stage, elapsed: Duration::ZERO, visited: false });
                if run.active >= at { run.active += 1; }
                if run.focus >= at { run.focus += 1; }
                at
            })
        });
        if let Some(next) = next {
            let previous = run.active;
            if previous != next {
                let now = Instant::now();
                let elapsed = now - run.since;
                run.rows[previous].elapsed += elapsed;
                run.since = now;
                run.active = next;
                run.rows[next].visited = true;
            }
            run.focus = run.focus.max(next);
            if stage == Stage::Align && run.align_start.is_none() {
                run.align_start = Some(Instant::now());
            }
        }
        drop(run);
        if let Some(bar) = &self.bar { bar.tick(); }
        self.report(true);
    }

    /// `reads` more reads aligned; `input_pos` bytes of the current query file consumed by them.
    pub(crate) fn progress(&self, reads: usize, input_pos: u64) {
        {
            let mut run = self.run.lock().unwrap();
            run.reads += reads;
            run.reads_total += reads;
            run.input_pos = run.input_pos.max(input_pos);
        }
        if let Some(bar) = &self.bar { bar.tick(); }
        self.report(false);
    }

    pub(crate) fn bam(&self, done: u64, total: u64) {
        self.run.lock().unwrap().bam = Some((done.min(total), total));
        self.report(false);
    }

    pub(crate) fn notice(&self, message: fmt::Arguments<'_>) {
        match &self.bar {
            Some(bar) => bar.println(message.to_string()),
            None => eprintln!("{message}"),
        }
    }

    pub(crate) fn suspend(&self, operation: impl FnOnce()) {
        match &self.bar {
            Some(bar) => bar.suspend(operation),
            None => operation(),
        }
    }

    pub(crate) fn finish(&self, ok: bool) {
        self.stop(if ok { Outcome::Complete } else { Outcome::Failed });
    }

    /// The run was stopped by a signal: draw the final panel, marked interrupted.
    pub(crate) fn interrupted(&self) {
        self.stop(Outcome::Interrupted);
    }

    fn stop(&self, outcome: Outcome) {
        let mut run = self.run.lock().unwrap_or_else(|e| e.into_inner());
        if run.finished.is_some() {
            return;
        }
        if let Some(bar) = &self.bar { bar.disable_steady_tick(); }
        let now = Instant::now();
        let elapsed = now - run.since;
        let active = run.active;
        run.rows[active].elapsed += elapsed;
        run.finished = Some((now, outcome));
        drop(run);
        self.report(true);
        if let Some(bar) = &self.bar {
            bar.finish();
            let _ = io::Write::write_all(&mut io::stderr(), b"\n");
        }
    }
}

// env_logger writes each formatted diagnostic to this sink as a complete record.
impl io::Write for Tui {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if !buf.is_empty() {
            match &self.bar {
                Some(bar) => bar.println(String::from_utf8_lossy(buf).trim_end()),
                None => io::Write::write_all(&mut io::stderr(), buf)?,
            }
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
