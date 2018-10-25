/* TODO: support more chromosome && parallelize */
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
use clap::{App, load_yaml};

fn save_scan_result(scanner:&Scanner, ir_path: &str) -> Result<(), ()>
{
    scanner.try_dump(&mut std::fs::File::create(ir_path).unwrap()).expect("cannot dump");

    return Ok(());
}

fn dump_frontend_events<DM:crate::depth_model::DepthModel + std::fmt::Debug>(frontend: &Frontend<DM>, fp : &mut std::fs::File)
{
    use std::io::Write;
    for event in frontend.iter() 
    { 
        write!(fp, "{:?}\n", event); 
    }
}

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
    let options = load_yaml!("cli.yml");
    let matches = App::from_yaml(options).get_matches();

    let alignment = matches.value_of("alignment-file").unwrap();
    let scanner_dump = matches.value_of("scanner-dump-path").unwrap_or(alignment);
    let no_scanner_dump = matches.is_present("no-scanner-dump");

    let ir_path = format!("{}.limodump-{}", scanner_dump, 0);

    let scanner = if no_scanner_dump || !std::path::Path::new(&ir_path[0..]).exists()
    {
        eprintln!("Scanner dump is not available, load from the alignment file");
        let bam = BamFile::new(scanner_dump, 0, None)?;
        let scanner = Scanner::new(&bam)?;
        if !no_scanner_dump { save_scan_result(&scanner, &ir_path[0..])?; }
        scanner
    }
    else
    {
        eprintln!("Loading depth information from scanner dump");
        Scanner::try_load(&mut std::fs::File::open(ir_path).unwrap()).expect("Cannot load")
    };

    let copy_nums = matches.value_of("copy-nums").unwrap();
    let copy_nums:Vec<_> =  copy_nums.split(",").map(|s| u32::from_str_radix(s, 10).unwrap()).collect();

    let window_size = u32::from_str_radix(matches.value_of("window-size").unwrap_or("300"), 10).unwrap();
    
    let frontend = Frontend::<LinearModel>::new(scanner, window_size, &copy_nums[0..], None)?;

    if matches.is_present("dump-frontend-events")
    {
        let mut output = std::fs::File::create(matches.value_of("dump-frontend-events").unwrap()).expect("Unable to open the file");
        dump_frontend_events(&frontend, &mut output);
    }

    let event_pair: Vec<_> = EventPairProc::new(&frontend).collect();

    if matches.is_present("dump-event-pairs")
    {
        let mut output = std::fs::File::create(matches.value_of("dump-event-pairs").unwrap()).expect("Unable to open the file");
        dump_event_pairs(&event_pair, &mut output);
    }

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
