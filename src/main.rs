mod hts;
mod window;
mod bamfile;
mod depth_model;
mod scanner;
mod histogram;
mod models;
mod frontend;
mod event_pair;
mod edge;

use self::bamfile::BamFile;
use self::frontend::Frontend;
use self::models::linear::LinearModel;
use self::scanner::Scanner;
use self::event_pair::EventPairProc;
use self::edge::EdgeDetector;

const BAM_PATH:&'static str = "/uufs/chpc.utah.edu/common/home/u0875014/data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam";
const IR_PATH:&'static str = "/uufs/chpc.utah.edu/common/home/u0875014/limo-development/data/NA12878.chr1.limo";

fn save_scan_result(scanner:&Scanner) -> Result<(), ()>
{
    scanner.try_dump(&mut std::fs::File::create(IR_PATH).unwrap()).expect("cannot dump");

    return Ok(());
}

#[allow(dead_code)]
fn dump_frontend_events<DM:crate::depth_model::DepthModel + std::fmt::Debug>(frontend: &Frontend<DM>, fp : &mut std::fs::File)
{
    use std::io::Write;
    for event in frontend.iter() 
    { 
        write!(fp, "{:?}\n", event); 
    }
}

#[allow(dead_code)]
fn dump_event_pairs<'a, DM:crate::depth_model::DepthModel + std::fmt::Debug>(event_pair:&Vec<(crate::frontend::Event<'a, DM>, crate::frontend::Event<'a, DM>)>, fp : &mut std::fs::File)
{
    use std::io::Write;

    for ep in event_pair
    {
        write!(fp, "{}\t{}\t{}\t{}", ep.0.chrom, ep.0.pos, ep.1.pos, format!("ls:{:?};rs:{:?};cn:{}", ep.0.score, ep.1.score, ep.0.copy_num));
    }
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

    let frontend = Frontend::<LinearModel>::new(scanner, 300, &target_copy_nums, None)?;


    let event_pair:Vec<_> = EventPairProc::new(&frontend).collect();

    let mut edge_detect = EdgeDetector::new(&frontend, frontend.get_scan_size() * 2);
    
    for ref ep in event_pair
    {
        if let Some(sv) = edge_detect.detect_edge(ep) 
        {
            println!("{}\t{}\t{}", sv.chrom, sv.left_pos, sv.right_pos);
        }
    }
    
    return Ok(()); 
}
