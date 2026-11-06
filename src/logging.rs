//! Timing and memory instrumentation mirroring the `yak_*` helpers of the
//! original implementation (`sys.cpp`).
//!
//! Only the numbers used inside log lines are provided here; the log lines
//! themselves are emitted at their call sites so that their wording stays
//! directly comparable with the C++ sources.

use std::sync::OnceLock;
use std::time::Instant;

static START_INSTANT: OnceLock<Instant> = OnceLock::new();

/// Number of clock ticks per second reported by `/proc/self/stat`.
///
/// `sysconf(_SC_CLK_TCK)` reports 100 on every Linux configuration relevant
/// here; the value is used for CPU accounting in log lines only, never for
/// anything that influences the assembly.
const CLOCK_TICKS_PER_SECOND: f64 = 100.0;

/// Marks the beginning of the run, equivalent to `yak_reset_realtime()`.
pub fn reset_realtime() {
    let _ = START_INSTANT.set(Instant::now());
}

fn start_instant() -> Instant {
    *START_INSTANT.get_or_init(Instant::now)
}

/// Wall clock seconds elapsed since [`reset_realtime`].
pub fn realtime_seconds() -> f64 {
    start_instant().elapsed().as_secs_f64()
}

/// User plus system CPU seconds consumed by this process, like `yak_cputime()`.
pub fn cpu_seconds() -> f64 {
    let stat_content = match std::fs::read_to_string("/proc/self/stat") {
        Ok(content) => content,
        Err(_) => return 0.0,
    };

    // The process name is wrapped in parentheses and may itself contain
    // spaces, so the numeric fields are split off after the last ')'.
    let fields_after_name = match stat_content.rfind(')') {
        Some(position) => &stat_content[position + 1..],
        None => return 0.0,
    };

    let numeric_fields: Vec<&str> = fields_after_name.split_whitespace().collect();
    // utime is field 14 and stime is field 15 of the full record; the first
    // entry here is field 3 (state), so both shift by eleven positions.
    let user_ticks: f64 = numeric_fields
        .get(11)
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.0);
    let system_ticks: f64 = numeric_fields
        .get(12)
        .and_then(|value| value.parse().ok())
        .unwrap_or(0.0);

    (user_ticks + system_ticks) / CLOCK_TICKS_PER_SECOND
}

/// Peak resident set size in bytes, like `yak_peakrss()`.
pub fn peak_resident_memory_bytes() -> u64 {
    let status_content = match std::fs::read_to_string("/proc/self/status") {
        Ok(content) => content,
        Err(_) => return 0,
    };

    for line in status_content.lines() {
        if let Some(rest) = line.strip_prefix("VmHWM:") {
            let kilobytes: u64 = rest
                .split_whitespace()
                .next()
                .and_then(|value| value.parse().ok())
                .unwrap_or(0);
            return kilobytes * 1024;
        }
    }

    0
}

/// Peak resident set size in gigabytes, like `yak_peakrss_in_gb()`.
pub fn peak_resident_memory_gigabytes() -> f64 {
    peak_resident_memory_bytes() as f64 / 1_073_741_824.0
}

/// CPU usage as a fraction of the wall clock time, like `yak_cpu_usage()`.
pub fn cpu_usage() -> f64 {
    (cpu_seconds() + 1e-9) / (realtime_seconds() + 1e-9)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn realtime_grows_after_reset() {
        reset_realtime();
        std::thread::sleep(std::time::Duration::from_millis(5));
        assert!(realtime_seconds() > 0.004);
    }

    #[test]
    fn peak_resident_memory_is_positive() {
        assert!(peak_resident_memory_bytes() > 0);
    }

    #[test]
    fn cpu_usage_is_positive() {
        assert!(cpu_usage() > 0.0);
    }
}
