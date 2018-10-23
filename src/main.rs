mod hts;
mod window;
mod bamfile;
mod depth_model;
mod scanner;
mod histogram;
mod models;
mod frontend;
mod event_pair;

use self::bamfile::BamFile;
use self::frontend::Frontend;
use self::models::linear::LinearModel;
use self::event_pair::EventPairProc;
use self::scanner::Scanner;

#[allow(dead_code)]
fn scan_bam() -> Result<(), ()>
{
    let bam = BamFile::new("/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam", 0, None)?;

    let scanner = Scanner::new(&bam)?;

    scanner.try_dump(&mut std::fs::File::create("/uufs/chpc.utah.edu/common/home/u0875014/limo-development/run/NA12878.limo-scan").unwrap()).expect("cannot dump");

    return Ok(());
}

fn main() -> Result<(), ()>
{
    //return scan_bam();

    let scanner = Scanner::try_load(&mut std::fs::File::open("/uufs/chpc.utah.edu/common/home/u0875014/limo-development/run/NA12878.limo-scan").unwrap()).expect("Cannot load");

    let target_copy_nums = [0,1];

    let mut frontend = Frontend::<LinearModel>::new(scanner, 300, &target_copy_nums, None)?;
    
    for (left, right) in EventPairProc::new(&mut frontend)
    {
        println!("{}\t{}\t{}\t{}\t{}", left.chrom, left.pos, right.pos, left.score, right.score); 
    }

    return Ok(()); 
}
