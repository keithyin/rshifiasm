//! Entry point, the counterpart of `main.cpp`.

use std::process::ExitCode;

use rshifiasm::assembly;
use rshifiasm::command_line::{AssemblyOptions, parser};
use rshifiasm::logging;

fn main() -> ExitCode {
    logging::reset_realtime();

    let arguments: Vec<String> = std::env::args().collect();
    let mut options = AssemblyOptions::default();

    if !parser::parse_arguments(&arguments, &mut options) {
        // The original `main` returns 0 whenever the command line is not
        // usable (help text, version, or a reported error).
        return ExitCode::SUCCESS;
    }

    let status = if options.secondary_input_files.is_some() {
        assembly::assemble_pair(&mut options)
    } else if options.debug_overlap_calculation != 0 {
        assembly::error_correction_debug(&mut options)
    } else {
        assembly::assemble(&mut options)
    };

    eprintln!(
        "[M::main] Version: {}",
        rshifiasm::command_line::options::VERSION
    );
    eprint!("[M::main] CMD:");
    for argument in &arguments {
        eprint!(" {}", argument);
    }
    eprintln!();
    eprintln!(
        "[M::main] Real time: {:.3} sec; CPU: {:.3} sec; Peak RSS: {:.3} GB",
        logging::realtime_seconds(),
        logging::cpu_seconds(),
        logging::peak_resident_memory_gigabytes()
    );

    if status == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::from(status as u8)
    }
}
