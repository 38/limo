mod hts;
mod window;
mod bamfile;
mod depth_model;
mod scanner;

use self::bamfile::BamFile;
use self::scanner::Scanner;
/*use self::hts::{hts_open, hts_close, sam_hdr_read, sam_index_load, bam_read1, bam_init1};
use std::ffi::CString;*/
fn main() -> Result<(), ()>
{
    let bam = BamFile::new("/home/haohou/base2/lumpy-sv-2/lumpy_tests/rice/AL87.discordant.sort.bam", 1, None)?;

    let _scanner = Scanner::new(&bam); 

    Ok(())
}
