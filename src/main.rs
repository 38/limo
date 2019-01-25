#![feature(core_intrinsics)]
mod edge;
mod task;


use self::task::Task;
use frontend::bamfile::BamFile;
use clap::{App, load_yaml};
use threadpool::ThreadPool;
use regex::Regex;

use std::str::FromStr;

use log::{info,debug};

fn main() -> Result<(), ()>
{
    let options = load_yaml!("cli.yml");
    let matches = App::from_yaml(options).get_matches();
    
    let copy_nums = matches.value_of("copy-nums").unwrap();
    let copy_nums:Vec<_> =  copy_nums.split(",").map(|s| u32::from_str_radix(s, 10).unwrap()).collect();

    let window_size = u32::from_str_radix(matches.value_of("window-size").unwrap_or("300"), 10).unwrap();
    let alignment = matches.value_of("alignment-file").unwrap();

    let mut nthreads = matches.value_of("threads").map(|s| usize::from_str_radix(s, 10).unwrap()).unwrap_or(1);


    let include_pattern = Regex::new(matches.value_of("include").unwrap_or(r"^([Cc]hr)?[0-9XYxy]*$")).unwrap();
    let exclude_pattern = Regex::new(matches.value_of("exclude").unwrap_or(".^")).unwrap();

    stderrlog::new()
        .module(module_path!())
        .module(frontend::get_module_path())
        .verbosity(3)
        .timestamp(stderrlog::Timestamp::Second)
        .init()
        .expect("Unable to initialize logging");

    let target_list:Vec<_> = BamFile::list_chromosomes(alignment)?.into_iter().enumerate()
        .filter(|(_, name)| include_pattern.is_match(&name) && !exclude_pattern.is_match(&name))
        .map(|(idx, name)| {
            debug!("Selected chromosome id={} name={}", idx, name);
            idx as u32
        }).collect(); 

    nthreads = nthreads.min(target_list.len());

    info!("Starting {} threads for {} chroms", nthreads, target_list.len());
    
    let tp = if nthreads > 1 { Some(ThreadPool::new(nthreads)) } else { None };

    for i in target_list.into_iter()
    {
        let task = Task {
            alignment: alignment.to_string(),
            scanner_dump: matches.value_of("scanner-dump-path").unwrap_or(alignment).to_string(),
            no_scanner_dump: matches.is_present("no-scanner-dump"),
            chrom: i,
            dump_fe: matches.value_of("dump-frontend-events").map(|x| x.to_string()),
            dump_ep: matches.value_of("dump-event-pairs").map(|x| x.to_string()),
            copy_nums: copy_nums.clone(),
            window_size: window_size,
            enable_pv: matches.value_of("prob-validate").map_or(true, |val| val != "off"),
            pv_threshold: matches.value_of("prob-validate").map_or(0.2, |val| f64::from_str(val).unwrap()),
            cluster_merge: !matches.is_present("no-cluster-merge"),
            load_events: matches.value_of("load-events").map(|x| x.to_string()),
        };

        if let Some(ref tp) = tp {
            tp.execute(move || { task.run().expect("Failed"); });
        } else {
            task.run().expect("Failed");
        }
    }

    if let Some(ref tp) = tp {
        tp.join();
    }

    return Ok(()); 
}
