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

fn main() -> Result<(), ()>
{
    let bam = BamFile::new("/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam", 1, None)?;

    let target_copy_nums = [0,1];

    let mut frontend = Frontend::<LinearModel>::new(&bam, 300, &target_copy_nums, 300)?;
    
    for (left, right) in EventPairProc::new(&mut frontend)
    {
        println!("{}\t{}\t{}\t{}\t{}", left.chrom, left.pos, right.pos, left.score, right.score); 
    }

    return Ok(());
}
