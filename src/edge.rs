use std::marker::PhantomData;
use crate::depth_model::DepthModel;
use crate::frontend::{Frontend, Event, Side};
use crate::histogram::Histogram;


pub struct EdgeDetector<'a, DM:DepthModel> {
    chrom: &'a str,
    scan_size: u32, 
    histogram: Histogram,
    raw_dep  : Box<Vec<i32>>,
    lmq_dep  : Box<Vec<i32>>,
    target_copy_num: Vec<u32>,
    phantom: PhantomData<&'a DM>
}

#[derive(Debug)]
pub struct Variant<'a> {
    pub chrom: &'a str,   
    pub left_pos: u32,
    pub right_pos: u32,
    pub copy_num : u32,
    pub mean     : f64,
    pub sd       : f64
}

impl <'a, DM:DepthModel> EdgeDetector<'a, DM> {
    #[allow(dead_code)]
    pub fn new(frontend:&'a Frontend<DM>, scan_size: u32, copy_nums: &[u32]) -> Self 
    {
        let raw_dep: Box<Vec<i32>> =  Box::new(frontend.get_scanner().get_raw_window().iter(1).collect());
        let lmq_dep: Box<Vec<i32>> = Box::new(frontend.get_scanner().get_low_mq_window().iter(1).collect());
        let mut histogram = Histogram::new(1024);
        raw_dep.iter().for_each(|v:&i32| histogram.add(*v as u32));
        let mut target_copy_num = vec![0u32;copy_nums.len()];
        target_copy_num[0..].clone_from_slice(copy_nums);
        target_copy_num[0..].sort();

        return Self { 
            chrom: frontend.get_scanner().get_chrom(),
            scan_size,
            raw_dep,
            lmq_dep,
            histogram,
            target_copy_num,
            phantom: PhantomData 
        };
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

    fn compute_norms(&mut self, left:u32, right:u32) -> (f64, f64)
    {
        let len = (right - left) as f64;
        let mut norm1 = 0f64;
        let mut norm2 = 0f64;
        for pos in left..right
        {
            norm1 += self.histogram.normalize(self.raw_dep[pos as usize] as u32);
            norm2 += (|x| x*x)(self.histogram.normalize(self.raw_dep[pos as usize] as u32));
        }

        norm1 /= len;
        norm2 /= len;

        return (norm1, norm2 - norm1 * norm1);
    }

    #[allow(dead_code)]
    pub fn detect_edge(&mut self, event: &(Event<'a, DM>, Event<'a, DM>)) -> Option<Variant>
    {
        let left_edges = self.scan_edge(&event.0);
        let right_edges = self.scan_edge(&event.1);

        let raw_dep_sum = self.raw_dep[(event.0.pos as usize)..(event.1.pos as usize)].iter().fold(0.0, |s,v| s + (*v as f64));
        let lmq_dep_sum = self.lmq_dep[(event.0.pos as usize)..(event.1.pos as usize)].iter().fold(0.0, |s,v| s + (*v as f64));
        
        if lmq_dep_sum > raw_dep_sum * 0.2 &&  raw_dep_sum > 100000.0 { return None; }

        let mut ret:Option<Variant> = None;
        let copy_num = event.0.copy_num;

        'outer: for (left_pos, left_rate) in &left_edges[0..]
        {
            for (right_pos, right_rate) in &right_edges[0..]
            {
                if left_rate * 2 < *right_rate || right_rate * 2 < *left_rate { continue; }
                if *left_pos >= *right_pos { continue; }

                let (avg, sd) = self.compute_norms(*left_pos, *right_pos);

                if ret.is_none() || 
                   (avg - 0.5 * (copy_num as f64)).abs() < (ret.as_ref().unwrap().mean - 0.5 * (copy_num as f64)).abs() ||
                   ((avg - ret.as_ref().unwrap().mean).abs() < 1e-5 && 
                    (sd < ret.as_ref().unwrap().sd))
                {
                    ret = Some(Variant {
                        chrom: self.chrom, 
                        left_pos: *left_pos,
                        right_pos: *right_pos,
                        copy_num,
                        mean: avg,
                        sd
                    });
                    break 'outer;
                }
            }
        }
        
        return match ret {
            Some(ref variant) => if (0.5 * (copy_num as f64) - variant.mean).abs() < 0.2 { ret } else 
            { 
                for copy_num in self.target_copy_num.iter()
                {
                    if (0.5 * (*copy_num as f64) - variant.mean).abs() < 0.2 
                    {
                        ret.as_mut().unwrap().copy_num = *copy_num;
                        return ret;
                    }
                }
                None
            }
            None => None
        };
    }
}
