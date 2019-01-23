#![feature(core_intrinsics)]

pub mod hts;
pub mod window;
pub mod bamfile;
pub mod depth_model;
pub mod scanner;
pub mod histogram;
pub mod models;
pub mod frontend;
pub mod event_pair;


pub mod prelude {
    
    use crate::bamfile::BamFile;
    use crate::frontend::{Frontend,Event};
    use crate::models::linear::LinearModel;
    use crate::scanner::Scanner;
    use crate::event_pair::EventPairProc;
    use crate::depth_model::DepthModel;

    #[derive(Clone)]
    pub struct FrontendParam<'a> {
        pub alignment:&'a str, 
        pub scanner_dump:&'a str, 
        pub no_scanner_dump:bool, 
        pub chrom:u32, 
        pub copy_nums:Vec<u32>, 
        pub window_size:u32, 
        pub dump_fe: Option<&'a str>, 
        pub dump_ep: Option<&'a str>
    }

    fn save_scan_result(scanner:&Scanner, ir_path: &str) -> Result<(), ()>
    {
        scanner.try_dump(&mut std::fs::File::create(ir_path).unwrap()).expect("cannot dump");

        return Ok(());
    }

    fn dump_frontend_events<DM:DepthModel + std::fmt::Debug>(frontend: &Frontend<DM>, fp : &mut std::fs::File)
    {
        use std::io::Write;
        for event in frontend.iter() 
        { 
            write!(fp, "{:?}\n", event).expect("IO error"); 
        }
    }

    fn dump_event_pairs<'a, DM:DepthModel + std::fmt::Debug>(event_pair:&Vec<(Event<'a, DM>, Event<'a, DM>)>, fp : &mut std::fs::File)
    {
        use std::io::Write;

        for ep in event_pair
        {
            write!(fp, "{}\t{}\t{}\t{}\n", ep.0.chrom, ep.0.pos, ep.1.pos, format!("ls:{:?};rs:{:?};cn:{}", ep.0.score, ep.1.score, ep.0.copy_num)).expect("IO error");
        }
    }

    pub struct Context<DM:DepthModel> {
        pub frontend : Frontend<DM>,
        fe_path: Option<String>,
        ep_path: Option<String>
    }

    impl <DM:DepthModel + std::fmt::Debug> Context<DM> {
        pub fn get_result<'a>(&'a self) -> Vec<(Event<'a, DM>, Event<'a, DM>)> 
        {
            if self.fe_path.is_some()
            {
                let output = std::fs::File::create(self.fe_path.as_ref().unwrap().as_str());
                dump_frontend_events(&self.frontend, &mut output.unwrap());
            }

            let event_pair = EventPairProc::new(&self.frontend).collect();

            if self.ep_path.is_some()
            {
                let output = std::fs::File::create(self.ep_path.as_ref().unwrap().as_str());
                dump_event_pairs(&event_pair, &mut output.unwrap());
            }

            event_pair
        }
    }

    pub fn run_linear_frontend<'a>(param: FrontendParam<'a>) -> Result<Context<LinearModel>, ()>
    {
        let ir_path = format!("{}.limodump-{}", param.scanner_dump, param.chrom);

        let scanner = if param.no_scanner_dump || !std::path::Path::new(&ir_path[0..]).exists()
        {
            eprintln!("Scanner dump is not available, load data from the alignment file: {} chromsome: {}", param.alignment, param.chrom);
            let bam = BamFile::new(param.alignment, param.chrom, None)?;
            let scanner = Scanner::new(&bam)?;
            if !param.no_scanner_dump { save_scan_result(&scanner, &ir_path[0..])?; }
            scanner
        }
        else
        {
            eprintln!("Loading depth information from scanner dump for file: {} chromsome: {}", param.alignment, param.chrom);
            Scanner::try_load(&mut std::fs::File::open(ir_path).unwrap()).expect("Cannot load")
        };
    
        let ret = Context{ 
            frontend: Frontend::<LinearModel>::new(scanner, param.window_size, &param.copy_nums[0..], None)?,
            fe_path: if let Some(fe) = param.dump_fe { Some(format!("{}-{}", fe, param.chrom)) } else {None},
            ep_path: if let Some(ep) = param.dump_ep { Some(format!("{}-{}", ep, param.chrom)) } else {None}
        };

        return Ok(ret);
    }

    pub fn make_linear_event<'a>(chrom: &'a str, left: u32, right: u32, copy_num: u32) -> (Event<'a, LinearModel>, Event<'a, LinearModel>)
    {
        let left = Event {
            chrom: chrom,
            side : crate::frontend::Side::Left,
            score: 0.0,
            pos  : left,
            copy_num: copy_num,
            total_dep: 0,
            lowmq_dep: 0
        };
        
        let right = Event {
            chrom: chrom,
            side : crate::frontend::Side::Right,
            score: 0.0,
            pos  : right,
            copy_num: copy_num,
            total_dep: 0,
            lowmq_dep: 0
        };

        return (left, right);
    }
}
