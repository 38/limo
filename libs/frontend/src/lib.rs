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

    pub struct Frontend<'a, DM:DepthModel> {
        alignment:&'a str, 
        scanner_dump:&'a str, 
        no_scanner_dump:bool, 
        chrom:u32, 
        copy_nums:Vec<u32>, 
        window_size:u32, 
        dump_fe: Option<&'a str>, 
        dump_ep: Option<&'a str>,
        pub result: Result<Vec<(Event<'a, DM>, Event<'a, DM>), ()>>,
        pub scanner:Scanner,
        frontend: DM
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
            write!(fp, "{:?}\n", event); 
        }
    }

    fn dump_event_pairs<'a, DM:DepthModel + std::fmt::Debug>(event_pair:&Vec<(Event<'a, DM>, Event<'a, DM>)>, fp : &mut std::fs::File)
    {
        use std::io::Write;

        for ep in event_pair
        {
            write!(fp, "{}\t{}\t{}\t{}\n", ep.0.chrom, ep.0.pos, ep.1.pos, format!("ls:{:?};rs:{:?};cn:{}", ep.0.score, ep.1.score, ep.0.copy_num));
        }
    }

    impl Frontend {
        pub fn run_linear_frontend(&mut self) -> Result<Vec<(Event, Event)>, ()> 
        {
            let ir_path = format!("{}.limodump-{}", scanner_dump, 0);

            let scanner = if no_scanner_dump || !std::path::Path::new(&ir_path[0..]).exists()
            {
                eprintln!("Scanner dump is not available, load from the alignment file");
                let bam = BamFile::new(alignment, 0, None)?;
                let scanner = Scanner::new(&bam)?;
                if !no_scanner_dump { save_scan_result(&scanner, &ir_path[0..])?; }
                scanner
            }
            else
            {
                eprintln!("Loading depth information from scanner dump");
                Scanner::try_load(&mut std::fs::File::open(ir_path).unwrap()).expect("Cannot load")
            };
            
            let frontend = Frontend::<LinearModel>::new(scanner, window_size, &copy_nums[0..], None)?;

            if dump_fe.is_some()
            {
                let mut output = std::fs::File::create(dump_fe.unwrap()).expect("Unable to open the file");
                dump_frontend_events(&frontend, &mut output);
            }

            let event_pair: Vec<_> = EventPairProc::new(&frontend).collect();

            if dump_ep.is_some()
            {
                let mut output = std::fs::File::create(dump_ep.unwrap()).expect("Unable to open the file");
                dump_event_pairs(&event_pair, &mut output);
            }

            return Ok(event_pair);
        }
    }
}
