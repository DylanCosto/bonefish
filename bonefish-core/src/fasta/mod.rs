pub mod record;
pub mod reader;
pub mod stream;
#[cfg(not(target_arch = "wasm32"))]
pub mod bam;

pub use record::Record;
pub use reader::{Reader, open, FastxError, parse_fasta_bytes};
pub use stream::{FastaStreamer, FastqStreamer};

pub use reader::read_fasta;
#[cfg(not(target_arch = "wasm32"))]
pub use reader::input_read_groups;
