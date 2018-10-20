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
use self::frontend::{Frontend, Event, Side};
#[allow(unused_imports)]
use self::models::linear::LinearModel;

use self::depth_model::DepthModel;

fn main() -> Result<(), ()>
{
    let bam = BamFile::new("/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam", 1, None)?;

    /*let scanner = Scanner::new(&bam)?; 
    let mut hist = Histogram::new(1024);

    scanner.get_raw_window().iter(1).for_each(|v:i32| hist.add(v as u32));

    eprintln!("{}", hist.get_average());
    eprintln!("{}", hist.normalize(70));*/

    let target_copy_nums = [0,1];

    let max_copy_num = target_copy_nums.iter().fold(0, |a,b| std::cmp::max(a,*b)); 

    let mut frontend = Frontend::<LinearModel>::new(&bam, 300, &target_copy_nums, 300)?;

    let mut recent = Vec::<Event<LinearModel>>::new();

    let mut left_items : Vec<Option<Event<LinearModel>>> = vec![None; max_copy_num as usize];

    let mut last_pos = 0u32;

    for item in frontend.iter()
    {
        if item.copy_num > max_copy_num { continue; }

        if recent.len() == 3 
        {
            recent.remove(0);
        }

        recent.push(item);

        if recent.len() == 3 && 
           recent[1].score - recent[0].score < LinearModel::EPS &&
           recent[1].score - recent[2].score < LinearModel::EPS
        {
            let current = &recent[1];

            if last_pos + 1 != recent[2].pos
            {
                match current.side 
                {
                    Side::Left => {
                        /* TODO: filter on the average exclude-rate and length */
                        if left_items[current.copy_num as usize].is_none() ||
                           left_items[current.copy_num as usize].as_ref().unwrap().score > current.score
                        {
                            left_items[current.copy_num as usize] = Some(Clone::clone(current));
                        }
                    },
                    Side::Right => {
                        let mut should_reset = false;
                        if let Some(ref left_side) = left_items[current.copy_num as usize]
                        {
                            if left_side.score < 7000.0 && current.score < 7000.0
                            {
                                println!("{:?} {:?}", left_side, current); 
                                should_reset = true;
                            }
                        }

                        if should_reset
                        {
                            left_items[current.copy_num as usize] = None;
                        }
                    }
                }
            }

            last_pos = recent[2].pos;
        }
    }

    return Ok(());
}
