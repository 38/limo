pub mod edge;
pub mod task;

use clap::{App, load_yaml};
use std::str::FromStr;
use self::task::Task;
use threadpool::ThreadPool;


fn main() -> Result<(), ()>
{
    let options = load_yaml!("cli.yml");
    let matches = App::from_yaml(options).get_matches();
    
    let copy_nums = matches.value_of("copy-nums").unwrap();
    let copy_nums:Vec<_> =  copy_nums.split(",").map(|s| u32::from_str_radix(s, 10).unwrap()).collect();

    let window_size = u32::from_str_radix(matches.value_of("window-size").unwrap_or("300"), 10).unwrap();
    let alignment = matches.value_of("alignment-file").unwrap();

    let tp = ThreadPool::new(20);

    for i in 0..23
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
            pv_threshold: matches.value_of("prob-validate").map_or(0.2, |val| f64::from_str(val).unwrap())
        };

        tp.execute(move || { task.run().expect("Failed"); });
    }

    tp.join();

    return Ok(()); 
}
