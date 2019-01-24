use std::cmp::{max,min};
use frontend::prelude::*;
use crate::edge::{EdgeDetector, Variant};

pub struct Task {
    pub alignment: String,
    pub scanner_dump: String,
    pub no_scanner_dump: bool,
    pub chrom: u32,
    pub dump_fe: Option<String>,
    pub dump_ep: Option<String>,
    pub copy_nums: Vec<u32>,
    pub window_size: u32,
    pub enable_pv: bool,
    pub pv_threshold: f64,
    pub cluster_merge: bool,
    pub load_events: Option<String>,
}

impl Task {

    pub fn run(&self) -> Result<(), ()>
    {
        let frontend_param = FrontendParam {
            alignment: self.alignment.as_str(),
            scanner_dump: self.scanner_dump.as_str(),
            no_scanner_dump: self.no_scanner_dump,
            chrom: self.chrom,
            dump_fe: self.dump_fe.iter().fold(None, |_,x| Some(x.as_str())),
            dump_ep: self.dump_ep.iter().fold(None, |_,x| Some(x.as_str())),
            copy_nums: self.copy_nums.clone(),
            window_size: self.window_size
        };

        let prob_args = if self.enable_pv {
            Some((frontend_param.alignment, None, frontend_param.chrom, self.pv_threshold))
        } else { None };


        let frontend_ctx = run_linear_frontend(frontend_param.clone())?;

        eprintln!("Constructing event detection context for chrom {}", frontend_param.chrom);

        let mut edge_detect = EdgeDetector::new(&frontend_ctx.frontend, frontend_ctx.frontend.get_scan_size() * 2, &frontend_param.copy_nums[0..], prob_args);
        
        let mut events:Vec<_> = if self.load_events.is_none() {

            eprintln!("Collecting event pairs for chromosome #{}", frontend_param.chrom);
            let event_pair = frontend_ctx.get_result();
            
            eprintln!("Post processing the poped events: Chromosome: {}", frontend_param.chrom);
            let report_unit = 10000000;
            let total_mb = event_pair.last().iter().fold(0, |_, e| e.0.pos) / 1000000;
            let mut last_mb = report_unit;
            let mut event_count = 0;
            let mut passed = 0;

            event_pair.iter().filter_map(|ep| {
                if ep.0.pos > last_mb {
                    eprintln!("Postprocessed: Chromosome:{}, Offset:{}MB/{}MB, FE_Events:{}, Passed:{}", ep.0.chrom, last_mb/1000000, total_mb, event_count, passed);
                    last_mb = ((ep.0.pos + report_unit - 1) / report_unit) * report_unit;
                }
                event_count += 1;

                edge_detect.detect_edge(ep, true).map(|x| { passed += 1; x })
            }).collect()
        } else {
            edge_detect.load_variants(std::fs::File::open(self.load_events.as_ref().unwrap()).expect("Cannot open event file"))
        };

        let events = if self.cluster_merge {
            eprintln!("Merging the clustered events: Chromosome: #{}", frontend_param.chrom);
            events.sort_by(|a,b| a.left_pos.cmp(&b.left_pos));
            let mut cluster_range = (0,0);
            let mut cluster = Vec::<&Variant>::new();
            let mut result = Vec::<Variant>::new();

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
                        fn update<'a, 'b>(best:&'b Variant<'a>, sv:&'b Variant<'a>, merged: bool) -> &'b Variant<'a>
                        {
                            let pv_diff_thres_a = if merged { 0.15 } else { 0.005 };
                            let pv_diff_thres_b = if merged { 0.005 } else { -1.0 };
                            if (best.pv_score - sv.pv_score).abs() > pv_diff_thres_a { 
                                if best.pv_score < sv.pv_score { return sv; }
                            } else if best.boundary != sv.boundary {
                                if !best.boundary { return sv; }
                            } else if (best.pv_score - sv.pv_score).abs() > pv_diff_thres_b { 
                                if best.pv_score < sv.pv_score { return sv; }
                            }  else if ((best.mean - 0.5 * best.copy_num as f64).abs() - (sv.mean - 0.5 * sv.copy_num as f64).abs()).abs() > 0.005 { 
                                if (best.mean - 0.5 * best.copy_num as f64).abs() > (sv.mean - 0.5 * sv.copy_num as f64).abs() { return sv; }
                            } else if (best.sd - sv.sd).abs() > 0.005 {
                                if best.sd > sv.sd { return sv; }
                            }
                            return best;
                        };

                        /* Option 1: Select a best SV from the cluster */
                        let best = cluster.iter().skip(1).fold(cluster[0], |a,b| update(a, *b, false));

                        /* Option 2: Merge all the SV in the cluster */
                        let event_pair = make_linear_event(cluster[0].chrom, cluster[0].left_pos, cluster[cluster.len()-1].right_pos, best.copy_num);
                        let cluster_event = edge_detect.detect_edge(&event_pair, true);
                        let mut best = cluster_event.iter().fold(best, |a, x| update(a, &x, true)).clone();
                        
                        /* Option 3: Also, it's possible we are in the middle of a huge event */
                        if cluster[cluster.len()-1].right_pos - cluster[0].left_pos > 5000 {
                            let mut event_pair = make_linear_event(cluster[0].chrom, cluster[0].left_pos, cluster[cluster.len()-1].right_pos, best.copy_num);
                            if let Some(ret) = edge_detect.extend_region(&mut event_pair, 2000) {
                                best = update(&best, &ret, true).clone();
                            }
                        }

                        result.push(best);
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
        } else {
            events
        };

        for sv in events 
        {
            println!("{}\t{}\t{}\t{}", sv.chrom, sv.left_pos, sv.right_pos, sv.json_repr());
        }

        return Ok(());
    }
}
