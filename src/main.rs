pub mod edge;

/* TODO: support more chromosome && parallelize */
#[cfg(with_task)]
mod task;

use frontend::prelude::run_frontend;
use self::edge::EdgeDetector;
use clap::{App, load_yaml};


fn main() -> Result<(), ()>
{
    let options = load_yaml!("cli.yml");
    let matches = App::from_yaml(options).get_matches();

    let alignment = matches.value_of("alignment-file").unwrap();
    let scanner_dump = matches.value_of("scanner-dump-path").unwrap_or(alignment);
    let no_scanner_dump = matches.is_present("no-scanner-dump");
    let dump_frontend_events = matches.value_of("dump-frontend-events");
    let dump_event_pairs = matches.value_of("dump-event-pairs");

    let copy_nums = matches.value_of("copy-nums").unwrap();
    let copy_nums:Vec<_> =  copy_nums.split(",").map(|s| u32::from_str_radix(s, 10).unwrap()).collect();

    let window_size = u32::from_str_radix(matches.value_of("window-size").unwrap_or("300"), 10).unwrap();

    let event_pair = run_frontend(alignment, scanner_dump, no_scanner_dump, 0, copy_nums, window_size, dump_frontend_events, dump_event_pairs)?;
    
    let mut edge_detect = EdgeDetector::new(&frontend, frontend.get_scan_size() * 2, &copy_nums[0..]);

    let mut event_pair_count = 0;
    
    for ref ep in event_pair
    {
        event_pair_count += 1;
        if let Some(sv) = edge_detect.detect_edge(ep) 
        {
            println!("{}\t{}\t{}\t{{\"chr\":\"{}\",\"begin\":{},\"end\":{},\"copy_num\":{},\"dep_avg\":{:.3},\"dep_var\":{:.3},\"left_model\":{:.3},\"right_model\":{:.3}}}", 
                     sv.chrom, sv.left_pos, sv.right_pos,
                     sv.chrom, sv.left_pos, sv.right_pos,
                     sv.copy_num, sv.mean, sv.sd, ep.0.score, ep.1.score);
        }
    }

    eprintln!("# of event pair poped up: {}", event_pair_count);

    return Ok(()); 
}
