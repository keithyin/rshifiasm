//! Command line handling: option storage and parsing.

pub mod options;
pub mod parser;

pub use options::AssemblyOptions;
pub use parser::parse_arguments;
