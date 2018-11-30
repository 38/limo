use serde_derive::Serialize;
use std::marker::PhantomData;
use frontend::depth_model::DepthModel;
use frontend::frontend::{Frontend, Event, Side};
use frontend::histogram::Histogram;
use frontend::bamfile::BamFile;
use frontend::window::Window;
use std::cmp::Ord;

pub struct EdgeDetector<'a, DM:DepthModel> {
    chrom: &'a str,
    scan_size: u32,
    read_size: u32,
    histogram: Histogram,
    raw_dep  : Box<Vec<i32>>,
    lmq_dep  : Box<Vec<i32>>,
    target_copy_num: Vec<u32>,
    bamfile : Option<BamFile>,
    pv_threshold: f64,
    phantom: PhantomData<&'a DM>
}

#[derive(Debug, Serialize, Clone)]
pub struct Variant<'a> {
    pub chrom: &'a str,   
    pub left_pos: u32,
    pub right_pos: u32,
    pub copy_num : u32,
    pub mean     : f64,
    pub sd       : f64,
    pub pv_score : f64,
    pub boundary : bool
}

impl <'a> Variant<'a> {
    pub fn json_repr(&self) -> String 
    {
        return serde_json::to_string(self).unwrap_or_else(|_| "{\"error\":1}".to_string());
    }
}

impl <'a, DM:DepthModel> EdgeDetector<'a, DM> {
    #[allow(dead_code)]
    pub fn new(frontend:&'a Frontend<DM>, scan_size: u32, copy_nums: &[u32], alignment: Option<(&str, Option<&str>, u32, f64)>) -> Self 
    {
        let raw_dep: Box<Vec<i32>> =  Box::new(frontend.get_scanner().get_raw_window().iter(1).collect());
        let lmq_dep: Box<Vec<i32>> = Box::new(frontend.get_scanner().get_low_mq_window().iter(1).collect());
        let mut histogram = Histogram::new(1024);
        raw_dep.iter().for_each(|v:&i32| histogram.add(*v as u32));
        let mut target_copy_num = vec![0u32;copy_nums.len()];
        target_copy_num[0..].clone_from_slice(copy_nums);
        target_copy_num[0..].sort();

        let ret = Self { 
            chrom: frontend.get_scanner().get_chrom(),
            scan_size,
            raw_dep,
            lmq_dep,
            histogram,
            target_copy_num,
            pv_threshold: alignment.iter().fold(0.0, |_d,v| v.3),
            phantom: PhantomData,
            bamfile: if let Some((path, refer, chrom, _)) = alignment {
                Some(BamFile::new(path, chrom, refer).unwrap())
            } else { None},
            read_size: frontend.get_scanner().get_common_read_length()
        };

        return ret;
    }

    fn compute_normalized_change_rate(&self, pos:u32, event: &Event<'a, DM>) -> i32
    {
        let raw_rate = self.raw_dep[pos as usize] - self.raw_dep[(pos - 1) as usize];

        return match (event.copy_num, event.side) {
            (copy_num, Side::Left) if copy_num < 2 => -raw_rate,
            (copy_num, Side::Left) if copy_num > 2 => raw_rate,
            (copy_num, Side::Right) if copy_num < 2 => raw_rate,
            (copy_num, Side::Right) if copy_num > 2 => -raw_rate,
            _ => 0 
        };
    }

    fn scan_edge(&self, event: &Event<'a, DM>) -> Vec<(u32, i32)>
    {
        let (left, mut right) = (if event.pos < self.scan_size { 0 } else { event.pos - self.scan_size }, event.pos + self.scan_size);
        if right as usize > self.raw_dep.len() 
        {
            right = self.raw_dep.len() as u32;
        }

        let mut ret = Vec::new();

        let mut last_scores = [0i32;3];

        for i in left..right 
        {
            last_scores[(i%3) as usize] = self.compute_normalized_change_rate(i, event);

            if i >= left + 3
            {
                let a = ((i - 2) % 3) as usize;
                let b = ((i - 1) % 3) as usize;
                let c = (i % 3) as usize;

                if last_scores[b] > last_scores[a] && last_scores[b] > last_scores[c]
                {
                    ret.push((i - 1, last_scores[b]));
                }
            }
        }

        ret.sort_unstable_by(|a,b| b.1.cmp(&a.1));

        return ret.iter().take(5).map(|a| *a).collect();
    }

    fn compute_norms(&mut self, left:u32, right:u32) -> (f64, f64, f64)
    {
        let len = (right - left) as f64;
        let mut norm1 = 0f64;
        let mut norm2 = 0f64;
        let mut lmq_avg = 0f64;
        for pos in left..right
        {
            norm1 += self.histogram.normalize(self.raw_dep[pos as usize] as u32);
            norm2 += (|x| x*x)(self.histogram.normalize(self.raw_dep[pos as usize] as u32));
            lmq_avg += self.histogram.normalize(self.lmq_dep[pos as usize] as u32);
        }

        norm1 /= len;
        norm2 /= len;

        return (norm1, norm2 - norm1 * norm1, lmq_avg / len);
    }

    /* This is the valildation based on the Read-Depth-Fragment-Depth correction
     * It based on measuring how strange the event is compared to the nearby region
     * Since the depth correction can cancel out most of the systematic bais, the average
     * of the depth is assumed to be stable.
     *
     * Problem: This needs to revisit the bamfile and slow down the program 
     */
    fn pvalue_validation(&mut self, variant: &Variant) -> Result<f64,()>
    {
        if let Some(ref mut bamfile) = self.bamfile 
        {
            let mut length = (variant.right_pos - variant.left_pos) * 100;
            
            if length < 10000 
            {
                length = 10000;
            }

            if length > 100000 
            {
                length = 100000;
            }

            if length > variant.left_pos { length = variant.left_pos; } 

            let range = ((variant.left_pos - length) as usize, (variant.right_pos + length) as usize);

            let iter = bamfile.try_iter_range(range.0, range.1)?;

            let mut window = Window::<i32>::new(range.1 - range.0);
            let mut window_r = Window::<i32>::new(range.1 - range.0);

            for read in iter 
            {
                if read.mqual() == 0 { continue; }
                if read.get_flags() & 0x80d != 1 { continue; }
                if read.begin() < range.0 as u32 || read.ref_begin() < range.0 as u32 { continue; }
                if read.get_isize() > 0 
                {
                    let beg = read.begin() - range.0 as u32;
                    let end = read.begin() + self.read_size + (read.get_isize() as u32) - range.0 as u32;

                    if end > (range.1 - range.0) as u32 { continue; }

                    if end - beg < 100000
                    {
                        window.accumulate(beg as usize, end as usize, 1);
                    }
                }

                if read.ref_end() - range.0 as u32 > (range.1 - range.0) as u32 { continue; }

                window_r.accumulate((read.ref_begin() - range.0 as u32) as usize, (read.ref_end() - range.0 as u32) as usize, 1);
            }

            let mut normal_depth_value = Vec::new();
            let mut event_depth_value = Vec::new();

            let mut pos = range.0 as u32;
    
            for (p,r) in window.iter::<i32>(1).zip(window_r.iter::<i32>(1))
            {
                if (pos as usize) < range.0 + 1000 { pos += 1; continue; }

                let correct_depth = if p == 0 { 0.0 } else { (r as f64) / (p as f64) } ;

                if pos < variant.left_pos as u32 || pos > variant.right_pos as u32
                {
                    normal_depth_value.push(correct_depth);
                }
                else
                {
                    event_depth_value.push(correct_depth);
                }

                pos += 1;
            }

            let fcmp = |a:&f64,b:&f64| a.partial_cmp(b).unwrap();

            normal_depth_value.sort_unstable_by(fcmp);
            event_depth_value.sort_unstable_by(fcmp);

            let exclude_size = event_depth_value.len() / 5;

            if exclude_size > 0 
            {

                let exclude_range = [event_depth_value[exclude_size], event_depth_value[event_depth_value.len() - exclude_size]];

                let score = exclude_range.iter().map(|x| (match normal_depth_value[..].binary_search_by(|y| fcmp(y,x)) { 
                    Err(idx) => idx, 
                    Ok(idx) => idx 
                }) as f64 / (normal_depth_value.len() as f64));
                return Ok(score.fold(4.0, |s, x| s * (x - 0.5)));
            }
            else 
            {
                return Ok(1.0);
            }
        }
        else { Ok(1.0) }
    }

    #[allow(dead_code)]
    pub fn detect_edge<'b>(&'b mut self, event: &(Event<'a, DM>, Event<'a, DM>), retry:bool) -> Option<Variant<'a>>
    {
        let left_edges = self.scan_edge(&event.0);
        let right_edges = self.scan_edge(&event.1);

        let raw_dep_sum = self.raw_dep[(event.0.pos as usize)..(event.1.pos as usize)].iter().fold(0.0, |s,v| s + (*v as f64));
        let lmq_dep_sum = self.lmq_dep[(event.0.pos as usize)..(event.1.pos as usize)].iter().fold(0.0, |s,v| s + (*v as f64));
        
        if lmq_dep_sum > raw_dep_sum * 0.2 &&  raw_dep_sum > 100000.0 { return None; }

        let mut ret:Option<Variant> = None;
        let copy_num = event.0.copy_num;

        'outer: for (left_pos, left_rate) in (&left_edges[0..]).iter().take(5)
        {
            for (right_pos, right_rate) in (&right_edges[0..]).iter().take(5)
            {
                /* If the both side changes very slowly, then we do not need to limit the change
                 * rate becuse it's meaningless */
                if (left_rate * 2 < *right_rate || right_rate * 2 < *left_rate) && (*left_rate > 10 || *right_rate > 10)  { continue; }
                if *left_pos >= *right_pos { continue; }

                let (avg, sd, lmq_avg) = self.compute_norms(*left_pos, *right_pos);

                if ret.is_none() || 
                   (avg - 0.5 * (copy_num as f64)).abs() < (ret.as_ref().unwrap().mean - 0.5 * (copy_num as f64)).abs() ||
                   (avg - lmq_avg - 0.5 * (copy_num as f64)).abs() < (ret.as_ref().unwrap().mean - 0.5 * (copy_num as f64)).abs() ||
                   ((avg - ret.as_ref().unwrap().mean).abs() < 1e-5 && 
                    (sd < ret.as_ref().unwrap().sd))
                {
                    let mut best_avg = avg - lmq_avg;

                    if (avg - 0.5 * (copy_num as f64)).abs() < (best_avg - 0.5 * (copy_num as f64)).abs()
                    {
                        best_avg = avg;
                    }
                    ret = Some(Variant {
                        chrom: self.chrom, 
                        left_pos: *left_pos,
                        right_pos: *right_pos,
                        copy_num,
                        mean: best_avg,
                        sd,
                        pv_score: 1.0,
                        boundary: true
                    });

                    break 'outer;
                }
            }
        }


        if ret.is_some() &&
           (0.5 * (copy_num as f64) - ret.as_ref().unwrap().mean).abs() > 0.2
        {
            let mut left = event.0.pos;
            let mut right = event.1.pos;

            while (self.histogram.normalize(self.raw_dep[left as usize] as u32) - (copy_num as f64) * 0.5).abs() > 0.1 && left < right { left += 1 }
            while (self.histogram.normalize(self.raw_dep[right as usize] as u32) - (copy_num as f64) * 0.5).abs() > 0.1 && left < right { right -= 1 }

            if (right - left) * 2 > event.1.pos - event.0.pos && right - left > 200
            {

                let (avg,sd, _lmq) = self.compute_norms(left, right);

                if sd < 0.2
                {
                    ret = Some(Variant {
                        chrom: self.chrom,
                        left_pos: left,
                        right_pos: right,
                        copy_num: copy_num,
                        mean: avg,
                        sd,
                        pv_score: 1.0,
                        boundary: false
                    });
                }
            }
        }
        
        let mut result = match ret {
            Some(mut variant) => if (0.5 * (copy_num as f64) - variant.mean).abs() < 0.2 { Some(variant) } else 
            {
                (|| {
                    for copy_num in self.target_copy_num.iter()
                    {
                        if (0.5 * (*copy_num as f64) - variant.mean).abs() < 0.2 
                        {
                            variant.copy_num = *copy_num;
                            return Some(variant);
                        }
                    }
                    return None;
                })()
            }
            None => None
        };

        /* P-Value validation */
        if let Some(ref mut what) = result 
        {
            if what.copy_num < 2 
            {
                what.pv_score = self.pvalue_validation(what).unwrap_or(1.0);

            }
            /* TODO: implement the PV score validation for dups as well */
        }
        
        if result.iter().fold(1.0, |_x,y| y.pv_score) < self.pv_threshold { result = None }


        if retry && result.is_some() && result.as_ref().unwrap().pv_score < 0.0
        {
            let mut new_param = event.clone();
            let cur_copy_num = copy_num;
            for copy_num in self.target_copy_num.clone().iter()
            {
                if *copy_num == cur_copy_num || *copy_num == result.as_ref().unwrap().copy_num { continue; }
                new_param.0.copy_num = *copy_num;
                new_param.1.copy_num = *copy_num;

                if let Some(result) = self.detect_edge(&new_param, false) 
                {
                    return Some(result);
                }
            }
        }

        return result;
    }
}
