//! A Rust translation of hifiasm.
//!
//! The original project is a C++ assembler of several hundred thousand lines;
//! this crate mirrors its modules one to one so that every function has a
//! direct counterpart, which keeps the translation auditable and allows the
//! output to be compared byte for byte with the reference implementation.

pub mod alignment;
pub mod assembly;
pub mod command_line;
pub mod container;
pub mod correction;
pub mod counting;
pub mod logging;
pub mod overlap;
pub mod reads;
pub mod sequence;
pub mod utilities;

pub use command_line::{AssemblyOptions, parse_arguments};
