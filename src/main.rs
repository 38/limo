mod hts;
mod window;
mod bamfile;
mod depth_model;
mod scanner;
mod histogram;
mod models;
mod frontend;

use self::bamfile::BamFile;
#[allow(unused_imports)]
use self::scanner::Scanner;
#[allow(unused_imports)]
use self::histogram::Histogram;
#[allow(unused_imports)]
use self::frontend::{Frontend, Event};
#[allow(unused_imports)]
use self::models::linear::LinearModel;

fn main() -> Result<(), ()>
{
    let bam = BamFile::new("/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam", 1, None)?;

    /*let scanner = Scanner::new(&bam)?; 
    let mut hist = Histogram::new(1024);

    scanner.get_raw_window().iter(1).for_each(|v:i32| hist.add(v as u32));

    eprintln!("{}", hist.get_average());
    eprintln!("{}", hist.normalize(70));*/


    let mut frontend = Frontend::<LinearModel>::new(&bam, 300, &[0,1], 300)?;

    //let mut recent = [Vec::<Event<LinearModel>>::new(), Vec::<Event<LinearModel>>::new()];
    let mut recent = Vec::<Event<LinearModel>>::new();

    for item in frontend.iter()
    {
        if recent.len() == 3 
        {
            recent.remove(0);
        }

        recent.push(item);

        if recent.len() == 3 && 
           recent[1].score - recent[0].score < 1e-5 &&
           recent[1].score - recent[2].score < 1e-5 
        {
            println!("{:?}", recent[1]);
        }
    }

    return Ok(());
}
