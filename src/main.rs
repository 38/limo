pub mod edge;
pub mod task;

use self::task::Task;
use frontend::bamfile::BamFile;

use clap::{App, load_yaml};
use threadpool::ThreadPool;
use regex::Regex;

use std::str::FromStr;




fn main() -> Result<(), ()>
{
    let options = load_yaml!("cli.yml");
    let matches = App::from_yaml(options).get_matches();
    
    let copy_nums = matches.value_of("copy-nums").unwrap();
    let copy_nums:Vec<_> =  copy_nums.split(",").map(|s| u32::from_str_radix(s, 10).unwrap()).collect();

    let window_size = u32::from_str_radix(matches.value_of("window-size").unwrap_or("300"), 10).unwrap();
    let alignment = matches.value_of("alignment-file").unwrap();

    let nthreads = matches.value_of("threads").map(|s| usize::from_str_radix(s, 10).unwrap()).unwrap_or(1);

    let tp = ThreadPool::new(nthreads);

    let include_pattern = Regex::new(matches.value_of("include").unwrap_or(r"^([Cc]hr)?[0-9XYxy]$")).unwrap();
    let exclude_pattern = Regex::new(matches.value_of("exclude").unwrap_or(".^")).unwrap();

    let target_list:Vec<_> = BamFile::list_chromosomes(alignment)?.into_iter().enumerate()
        .filter(|(_, name)| include_pattern.is_match(&name) && !exclude_pattern.is_match(&name))
        .map(|(idx, name)| {
            eprintln!("Matched target #{}:{}", idx, name);
            idx as u32
        }).collect();

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
            enable_pv: (matches.value_of("prob-validate").unwrap_or("default") != "off"),
            pv_threshold: matches.value_of("prob-validate").iter().fold(0.2, |_,val| f64::from_str(val).unwrap())
        };

        tp.execute(move || { task.run().expect("Failed"); });
    }

    tp.join();

    return Ok(()); 
}
