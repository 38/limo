mod hts;
mod window;
mod bamfile;
mod depth_model;
mod scanner;
mod histogram;
mod models;
mod frontend;

use self::bamfile::BamFile;
use self::scanner::Scanner;
use self::histogram::Histogram;
/*use self::hts::{hts_open, hts_close, sam_hdr_read, sam_index_load, bam_read1, bam_init1};
use std::ffi::CString;*/
fn main() -> Result<(), ()>
{
    let bam = BamFile::new("/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam", 1, None)?;

    let scanner = Scanner::new(&bam)?; 
    let mut hist = Histogram::new(1024);

    scanner.get_raw_window().iter(1).for_each(|v:i32| hist.add(v as u32));

    eprintln!("{}", hist.get_average());
    eprintln!("{}", hist.normalize(70));

    Ok(())
}
