use crate::bamfile::BamFile;
use crate::depth_model::DepthModel;
use crate::frontend::{Event, Frontend};
use crate::scanner::Scanner;
use clap::ArgMatches;
use std::fmt::Debug;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;

#[derive(Clone)]
pub struct TaskDesc {
    bam_path: Arc<String>,
    copy_nums: Arc<Vec<u32>>,
    window_size: u32,
    scanner_dump_prefix: Option<Arc<String>>,
    event_dump_prefix: Option<Arc<String>>,
    paired_dump_prefix: Option<Arc<String>>,
}

pub struct Task {
    desc: TaskDesc,
    chrom_id: u32,
}

impl Task {
    pub fn new(desc: &TaskDesc, chrom: u32) -> Self {
        return Task {
            desc: desc.clone(),
            chrom_id: chrom,
        };
    }

    fn load_bam(path: &str) -> Result<Scanner, ()> {
        let bam = BamFile::new(path, 0, None)?;
        let scanner = Scanner::new(&bam)?;
        return Ok(scanner);
    }

    fn load_scanner_dump(path: &str) -> Result<Scanner, ()> {
        return Ok(
            Scanner::try_load(&mut File::open(path).unwrap()).expect("Cannot load scanner dump")
        );
    }

    fn dump_scanner(scanner: &Scanner, path: &str) -> Result<(), ()> {
        scanner
            .try_dump(&mut File::open(path).unwrap())
            .expect("Cannot dump the scanner");

        return Ok(());
    }

    fn dump_frontend_events<DM: DepthModel + Debug>(frontend: &Frontend<DM>, path: &str) {
        let fp = File::create(path).expect("Cannot open the frontend events");
        for event in frontend.iter() {
            write!(fp, "{:?}\n", event);
        }
    }

    fn dump_event_pairs<'a, DM: DepthModel + Debug>(
        event_pair: &Vec<(Event<'a, DM>, Event<'a, DM>)>,
        path: &str,
    ) {
        let fp = File::create(path).expect("Cannot open the event pairs");
        for ep in event_pair {
            write!(
                fp,
                "{}\t{}\t{}\t{}\n",
                ep.0.chrom,
                ep.0.pos,
                ep.1.pos,
                format!(
                    "ls:{:?};rs:{:?};cn:{}",
                    ep.0.score, ep.1.score, ep.0.copy_num
                )
            );
        }
    }

    pub fn run(&self) -> Result<(), ()> {
        let scanner = if let Some(scanner_dump_prefix) = self.desc.scanner_dump_prefix {
            let dump_path = format!("{}{}", scanner_dump_prefix, self.desc.chrom_id);
            if !Path::new(dump_path).exists() {
                let bam = BamFile::new(self.desc.bam_path, self.desc.chrom_id, None)?;
                let scanner = Scanner::new(&bam);
                Self::dump_scanner(&scanner, dump_path)?;
                scanner
            } else {
                Scanner::try_load(&mut File::open(dump_path).unwrap())
                    .expect("Cannot load the scanner dump")
            }
        } else {
            let bam = BamFile::new(self.desc.bam_path, self.desc.chrom_id, None)?;
            Scanner::new(&bam)
        };

        let frontend =
            Frontend::<LinearModel>::new(scanner, self.window_size, &self.copy_nums[0..], None)?;

        if let Some(event_dump_prefix) = self.desc.event_dump_prefix {}
    }
}

impl TaskDesc {
    pub fn from_arg_matches(am: &ArgMatches) -> TaskDesc {
        let bam_path = Arc::new(am.value_of("alignment-file").unwrap().to_owned());

        let copy_nums = am.value_of("copy-nums").unwrap();
        let copy_nums: Vec<_> = copy_nums
            .split(",")
            .map(|s| u32::from_str_radix(s, 10).unwrap())
            .collect();
        let copy_nums = Arc::new(copy_nums);

        let window_size =
            u32::from_str_radix(am.value_of("window-size").unwrap_or("300"), 10).unwrap();
        let scanner_dump_prefix = if am.is_present("no-scanner-dump") {
            None
        } else {
            let base = am
                .value_of("scanner-dump-path")
                .unwrap_or(am.value_of("alignment-file").unwrap());
            Some(Arc::new(format!("{}.limodump-", base)))
        };

        let event_dump_prefix = if am.is_present("dump-frontend-events") {
            let path = am.value_of("dump-frontend-events").unwrap();
            Some(Arc::new(path.to_owned()))
        } else {
            None
        };

        let paired_dump_prefix = if am.is_present("dump-event-pairs") {
            let path = am.value_of("dump-event-pairs").unwrap();
            Some(Arc::new(path.to_owned()))
        } else {
            None
        };

        return Self {
            bam_path,
            copy_nums,
            window_size,
            scanner_dump_prefix,
            event_dump_prefix,
            paired_dump_prefix,
        };
    }
}
