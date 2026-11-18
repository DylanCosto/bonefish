//! Removes this run's temporary files when it is interrupted (Ctrl-C, SIGTERM, SIGHUP). The
//! split-index directory and the sorter's runs are otherwise removed only on a normal exit, and
//! can be as large as the output.

use std::path::{Path, PathBuf};
use std::sync::Mutex;

use crate::tui::Tui;

enum Temp {
    /// A directory created for this run.
    Dir(PathBuf),
    /// The sorter's runs: `<prefix>.NNNN.run`.
    Runs(String),
}

static TEMPS: Mutex<Vec<Temp>> = Mutex::new(Vec::new());

fn temps() -> std::sync::MutexGuard<'static, Vec<Temp>> {
    TEMPS.lock().unwrap_or_else(|e| e.into_inner())
}

pub(crate) fn register_dir(dir: &Path) {
    temps().push(Temp::Dir(dir.to_path_buf()));
}

pub(crate) fn unregister_dir(dir: &Path) {
    temps().retain(|t| !matches!(t, Temp::Dir(d) if d == dir));
}

pub(crate) fn register_runs(prefix: &str) {
    temps().push(Temp::Runs(prefix.to_owned()));
}

/// Whether `name` is `<base>.NNNN.run`.
fn is_run(name: &str, base: &str) -> bool {
    name.strip_prefix(base)
        .and_then(|rest| rest.strip_prefix('.'))
        .and_then(|rest| rest.strip_suffix(".run"))
        .is_some_and(|n| n.len() >= 4 && n.bytes().all(|b| b.is_ascii_digit()))
}

fn remove_all() {
    for temp in temps().iter() {
        match temp {
            Temp::Dir(dir) => {
                // Retried: a file created during the removal makes it fail.
                for _ in 0..20 {
                    if std::fs::remove_dir_all(dir).is_ok() || !dir.exists() { break; }
                    std::thread::sleep(std::time::Duration::from_millis(10));
                }
            }
            Temp::Runs(prefix) => {
                let prefix = Path::new(prefix);
                let Some(base) = prefix.file_name().and_then(|b| b.to_str()) else { continue };
                let dir = prefix.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
                let Ok(entries) = std::fs::read_dir(dir) else { continue };
                for entry in entries.flatten() {
                    if entry.file_name().to_str().is_some_and(|name| is_run(name, base)) {
                        let _ = std::fs::remove_file(entry.path());
                    }
                }
            }
        }
    }
}

/// On SIGINT, SIGTERM or SIGHUP: mark the progress display interrupted, remove the temporary
/// files and exit with 128 + the signal number, as the default action would.
#[cfg(unix)]
pub(crate) fn install(tui: Option<Tui>) {
    use signal_hook::consts::{SIGHUP, SIGINT, SIGTERM};
    let Ok(mut signals) = signal_hook::iterator::Signals::new([SIGINT, SIGTERM, SIGHUP]) else { return };
    std::thread::spawn(move || {
        if let Some(signal) = signals.forever().next() {
            if let Some(tui) = &tui {
                tui.interrupted();
            }
            bonefish::sorted_bam::stop_new_runs();
            remove_all();
            std::process::exit(128 + signal);
        }
    });
}

#[cfg(not(unix))]
pub(crate) fn install(_tui: Option<Tui>) {}

#[cfg(test)]
mod tests {
    use super::is_run;

    #[test]
    fn matches_only_sorter_runs() {
        assert!(is_run("out.bam.tmp.0000.run", "out.bam.tmp"));
        assert!(is_run("out.bam.tmp.12345.run", "out.bam.tmp"));
        assert!(!is_run("out.bam.tmp.run", "out.bam.tmp"));
        assert!(!is_run("out.bam.tmp.00a0.run", "out.bam.tmp"));
        assert!(!is_run("out.bam", "out.bam.tmp"));
        assert!(!is_run("other.tmp.0000.run", "out.bam.tmp"));
        assert!(!is_run("out.bam.tmp.0000.run.keep", "out.bam.tmp"));
    }
}
