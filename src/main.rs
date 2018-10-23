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

const BAM_PATH:&'static str = "/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam";
const IR_PATH:&'static str = "/uufs/chpc.utah.edu/common/home/u0875014/limo-development/run/NA12878.limo-scan";

fn save_scan_result(scanner:&Scanner) -> Result<(), ()>
{
    scanner.try_dump(&mut std::fs::File::create(IR_PATH).unwrap()).expect("cannot dump");

    return Ok(());
}

fn main() -> Result<(), ()>
{
    let scanner = if !std::path::Path::new(IR_PATH).exists()
    {
        let bam = BamFile::new(BAM_PATH, 0, None)?;
        let scanner = Scanner::new(&bam)?;
        save_scan_result(&scanner)?;
        scanner
    }
    else
    {
        Scanner::try_load(&mut std::fs::File::open(IR_PATH).unwrap()).expect("Cannot load")
    };

    let target_copy_nums = [0,1];

    let mut frontend = Frontend::<LinearModel>::new(scanner, 300, &target_copy_nums, None)?;
    
    for (left, right) in EventPairProc::new(&mut frontend)
    {
        println!("{}\t{}\t{}\t{}\t{}", left.chrom, left.pos, right.pos, left.score, right.score); 
    }

    return Ok(()); 
}
