use std::fmt;
use std::io::{self, IsTerminal};
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

struct Row {
    stage: Stage,
    elapsed: Duration,
    visited: bool,
}

struct Run {
    started: Instant,
    since: Instant,
    finished: Option<(Instant, bool)>,
    active: usize,
    rows: Vec<Row>,
    part: usize,
    reads: usize,
    bam: Option<(u64, u64)>,
}

impl Run {
    fn draw(&self, out: &mut dyn fmt::Write, width: usize) -> fmt::Result {
        let (now, status) = match self.finished {
            None => (Instant::now(), "Running"),
            Some((end, true)) => (end, "Complete"),
            Some((end, false)) => (end, "Failed"),
        };
        let title = format!("bonefish  {status}  {}", time(now - self.started));
        let title = title
            .chars()
            .take(width.saturating_sub(1))
            .collect::<String>();
        writeln!(out, "{}", style(title).for_stderr().bold())?;

        // Keep the longest line short enough to avoid wrapping in small terminals.
        if width < 56 {
            let row = &self.rows[self.active];
            let detail = match row.stage {
                Stage::Align => format!(" · {} reads", self.reads),
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
        for (index, row) in self.rows.iter().enumerate() {
            let active = index == self.active && self.finished.is_none();
            let (mark, detail) = if active {
                let detail = match row.stage {
                    Stage::Align => {
                        if details >= 26 && self.part > 1 {
                            format!("{} reads · part {}", self.reads, self.part)
                        } else {
                            format!("{} reads", self.reads)
                        }
                    }
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
            } else if index == self.active && self.finished.is_some_and(|(_, ok)| !ok) {
                (" ", "failed".to_owned())
            } else if row.visited {
                let detail = match (row.stage, self.bam) {
                    (Stage::Align, _) if self.finished.is_some() && self.part == 1 => {
                        format!("done · {} reads", self.reads)
                    }
                    (Stage::Bam, Some((done, _))) => format!("done · {done} records"),
                    _ => "done".to_owned(),
                };
                (" ", detail)
            } else {
                (" ", "queued".to_owned())
            };
            let elapsed = if active {
                row.elapsed + (now - self.since)
            } else {
                row.elapsed
            };
            let clock = if row.visited {
                time(elapsed)
            } else {
                String::new()
            };
            let line = format!(
                "{mark} {:<19} {:<details$.details$} {:>8} ",
                row.stage.name(),
                detail,
                clock
            );
            if active {
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
    format!("{:02}:{:02}", s / 60, s % 60)
}

#[derive(Clone)]
pub(crate) struct Tui {
    run: Arc<Mutex<Run>>,
    bar: ProgressBar,
}

impl Tui {
    pub(crate) fn new(
        requested: bool,
        stdout_data: bool,
        align: bool,
        merge: bool,
        bam: bool,
    ) -> Option<Self> {
        if !requested
            || (stdout_data && io::stdout().is_terminal())
            || std::env::var_os("TERM").is_some_and(|term| term == "dumb")
        {
            return None;
        }
        let target = ProgressDrawTarget::stderr_with_hz(8);
        if target.is_hidden() {
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
            rows: stages
                .into_iter()
                .map(|stage| Row {
                    stage,
                    elapsed: Duration::ZERO,
                    visited: stage == Stage::Reference,
                })
                .collect(),
            part: 0,
            reads: 0,
            bam: None,
        }));
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
        let bar = ProgressBar::with_draw_target(None, target).with_style(theme);
        Some(Self { run, bar })
    }

    pub(crate) fn start(&self) {
        self.bar.tick();
        self.bar.enable_steady_tick(Duration::from_millis(200));
    }

    pub(crate) fn stage(&self, stage: Stage, part: usize) {
        let mut run = self.run.lock().unwrap();
        let next = run.rows.iter().position(|row| row.stage == stage).or_else(|| {
            // A merge that starts without --split-prefix (multi-part SAM/BAM): add its row
            // before the BAM stages.
            (stage == Stage::Merge).then(|| {
                let at = run.rows.iter().position(|row| matches!(row.stage, Stage::Bam | Stage::Index)).unwrap_or(run.rows.len());
                run.rows.insert(at, Row { stage, elapsed: Duration::ZERO, visited: false });
                if run.active >= at { run.active += 1; }
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
            if stage == Stage::Align && part != run.part {
                run.part = part;
                run.reads = 0;
            }
        }
        drop(run);
        self.bar.tick();
    }

    pub(crate) fn reads(&self, count: usize) {
        self.run.lock().unwrap().reads += count;
    }

    pub(crate) fn bam(&self, done: u64, total: u64) {
        self.run.lock().unwrap().bam = Some((done.min(total), total));
    }

    pub(crate) fn notice(&self, message: fmt::Arguments<'_>) {
        self.bar.println(message.to_string());
    }

    pub(crate) fn suspend(&self, operation: impl FnOnce()) {
        self.bar.suspend(operation);
    }

    pub(crate) fn finish(&self, ok: bool) {
        self.bar.disable_steady_tick();
        let mut run = self.run.lock().unwrap();
        let now = Instant::now();
        let elapsed = now - run.since;
        let active = run.active;
        run.rows[active].elapsed += elapsed;
        run.finished = Some((now, ok));
        drop(run);
        self.bar.finish();
        let _ = io::Write::write_all(&mut io::stderr(), b"\n");
    }
}

// env_logger writes each formatted diagnostic to this sink as a complete record.
impl io::Write for Tui {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if !buf.is_empty() {
            self.bar.println(String::from_utf8_lossy(buf).trim_end());
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
