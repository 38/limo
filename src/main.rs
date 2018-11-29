#[macro_use]
extern crate serde_derive;

pub mod edge;

/* TODO: support more chromosome && parallelize */
#[cfg(with_task)]
mod task;

use std::cmp::{max,min};
use frontend::prelude::*;
use self::edge::EdgeDetector;
use clap::{App, load_yaml};
use std::str::FromStr;


fn main() -> Result<(), ()>
{
    let options = load_yaml!("cli.yml");
    let matches = App::from_yaml(options).get_matches();
    
    let copy_nums = matches.value_of("copy-nums").unwrap();
    let copy_nums:Vec<_> =  copy_nums.split(",").map(|s| u32::from_str_radix(s, 10).unwrap()).collect();

    let window_size = u32::from_str_radix(matches.value_of("window-size").unwrap_or("300"), 10).unwrap();
    let alignment = matches.value_of("alignment-file").unwrap();
    let chromid = 0;

    let param = FrontendParam {
        alignment,
        scanner_dump: matches.value_of("scanner-dump-path").unwrap_or(alignment),
        no_scanner_dump: matches.is_present("no-scanner-dump"),
        chrom: chromid,
        dump_fe: matches.value_of("dump-frontend-events"),
        dump_ep: matches.value_of("dump-event-pairs"),
        copy_nums: copy_nums.clone(),
        window_size
    };

    let prob_args = match matches.value_of("prob-validate") {
        Some("off") => None,
        Some(value) => Some((alignment, None, 0, f64::from_str(value).unwrap())),
        None        => Some((alignment, None, 0, 0.2))
    };

    let frontend_ctx = run_linear_frontend(param)?;

    eprintln!("Collecting event pairs for chromosome #{}", chromid);
    let event_pair = frontend_ctx.get_result();
    
    let mut edge_detect = EdgeDetector::new(&frontend_ctx.frontend, frontend_ctx.frontend.get_scan_size() * 2, &copy_nums[0..], prob_args);
    
    eprintln!("Post processing the poped events");

    let report_unit = 10000000;
    let mut last_mb = report_unit;
    let mut event_count = 0;

    let mut events:Vec<_> = event_pair.iter().filter_map(|ep| {
        if ep.0.pos > last_mb {
            eprintln!("Postprocessed: Chromosome:{}, Offset:{}MB, Events: {}", ep.0.chrom, last_mb/1000000, event_count);
            last_mb = ((ep.0.pos + report_unit - 1) / report_unit) * report_unit;
        }
        event_count += 1;
        edge_detect.detect_edge(ep, true)
    }).collect();

    let events = {
        eprintln!("Merging the clustered events: Chromosome: #{}", chromid);
        events.sort_by(|a,b| a.left_pos.cmp(&b.left_pos));
        let mut cluster_range = (0,0);
        let mut cluster = Vec::<&crate::edge::Variant>::new();
        let mut result = Vec::<crate::edge::Variant>::new();

        let max_dist = 1000;

        let mut iter = events.iter();
        loop
        {
            let mut should_merge = false;
            let next_sv = iter.next();
            if let Some(sv) = next_sv
            {
                if max(sv.left_pos - max_dist, cluster_range.0) < min(sv.right_pos + max_dist, cluster_range.1) 
                {
                    cluster_range.0 = min(sv.left_pos - max_dist, cluster_range.0);
                    cluster_range.1 = max(sv.right_pos + max_dist, cluster_range.1);
                    cluster.push(sv);
                }
                else 
                {
                    should_merge = true;
                }
            }
            else
            {
                should_merge = true;
            }

            if should_merge
            {
                if cluster.len() > 1 
                {
                    fn update<'a, 'b>(best:&'b crate::edge::Variant<'a>, sv:&'b crate::edge::Variant<'a>) -> &'b crate::edge::Variant<'a>
                    {
                        if (best.pv_score < sv.pv_score) ||
                           (!best.boundary && sv.boundary) ||
                           (best.mean - 0.5 * best.copy_num as f64).abs() > (sv.mean - 0.5 * sv.copy_num as f64).abs() ||
                           (best.sd > sv.sd)
                        {
                            return sv;
                        }
                        return best;
                    };

                    /* Option 1: Select a best SV from the cluster */
                    let best = cluster.iter().skip(1).fold(cluster[0], |a,b| update(a, *b));

                    /* Option 2: Merge all the SV in the cluster */
                    let event_pair = make_linear_event(cluster[0].chrom, cluster[0].left_pos, cluster[cluster.len()-1].right_pos, best.copy_num);
                    let cluster_event = edge_detect.detect_edge(&event_pair, true);
                    let best = cluster_event.iter().fold(best, |a, x| update(a, &x));

                    result.push(crate::edge::Variant::clone(best));
                }
                else 
                {
                    if cluster.len() > 0 { result.push(cluster[0].clone()); }
                }

                if let Some(sv) = next_sv
                {
                    cluster.clear();
                    cluster.push(sv);
                    cluster_range = (sv.left_pos - max_dist, sv.right_pos + max_dist);
                }
                else
                {
                    break result;
                }
            }
        }
    };

    for sv in events 
    {
        println!("{}\t{}\t{}\t{}", sv.chrom, sv.left_pos, sv.right_pos, sv.json_repr());
    }
    
    return Ok(()); 
}
