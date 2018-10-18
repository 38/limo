use crate::bamfile::BamFile;
use crate::scanner::Scanner;
use crate::depth_model::DepthModel;
use crate::histogram::Histogram;
use crate::window::WindowIter;

#[derive(Debug, Copy, Clone)]
#[allow(dead_code)]
pub enum Side {
    Left,
    Right
}

struct SVModel<DM:DepthModel + Sized> {
    copy_num : u32,
    side     : Side,
    model    : DM
}

pub struct Frontend<DM:DepthModel + Sized> {
    scanner : Scanner,
    window_size: u32,
    copy_nums: Vec<u32>,
    dmp      : DM::ParamType
}

pub struct FrontendIter<'a, DM:DepthModel + Sized> {
    pos : u32,
    correct_iter: WindowIter<'a, i32, i32>,
    exclude_iter: WindowIter<'a, i32, i32>,
    hist    : Histogram,
    left_mod: Vec<SVModel<DM>>,
    right_mod: Vec<SVModel<DM>>
}

#[derive(Debug, Copy, Clone)]
pub struct Event<DM:DepthModel> {
    side : Side,
    score: DM::Output,
    pos  : u32,
    copy_num: u32
}

impl <DM:DepthModel + Sized> Frontend<DM> {
    #[allow(dead_code)]
    pub fn new(fp:&BamFile, window_size: u32, copy_nums:&[u32], dmp : DM::ParamType) -> Result<Self,()>
    {
        let scanner = Scanner::new(fp)?;

        let ret = Self {
            scanner,
            window_size,
            copy_nums:copy_nums.to_vec(),
            dmp
        };

        return Ok(ret);
    }
}

impl <'a, DM:DepthModel + Sized> FrontendIter<'a, DM> where DM::Input : From<f64>
{
    fn get_normalized_depth(&mut self) -> Option<DM::Input>
    {
        if let Some(correct) = self.correct_iter.next()
        {
            if let Some(excluded) = self.exclude_iter.next()
            {
                let normalized = self.hist.normalize((correct - excluded) as u32);
                return Some(From::from(normalized));
            }
        }
        return None;
    }

    fn feed(&mut self, dep : DM::Input)
    {
        self.left_mod[0..].iter_mut().for_each(|m| m.model.put(dep));
        self.right_mod[0..].iter_mut().for_each(|m| m.model.put(dep));
    }

    #[allow(dead_code)]
    fn new(obj:&'a mut Frontend<DM>) -> Self
    {
        let mut hist = Histogram::new(1024);
        obj.scanner.get_raw_window().iter(1).for_each(|v:i32| hist.add(v as u32));

        let size = obj.window_size + obj.scanner.get_common_read_length();
        let correct_iter = obj.scanner.get_corrected().iter(obj.window_size as usize);
        let exclude_iter = obj.scanner.get_low_mq_window().iter(obj.window_size as usize);
        
        let mut left_mod = Vec::<SVModel<DM>>::new();
        let mut right_mod = Vec::<SVModel<DM>>::new();
        

        for copy_num in obj.copy_nums[0..].iter()
        {
            let left_side = DM::create_model(*copy_num, true, obj.dmp);
            let right_side = DM::create_model(*copy_num, false, obj.dmp);

            left_mod.push(SVModel {
                copy_num : *copy_num,
                side:      Side::Left,
                model: left_side
            });

            right_mod.push(SVModel {
                copy_num: *copy_num,
                side: Side::Right,
                model: right_side
            });
        }
        
        let mut ret = Self {
            hist,
            pos: size,
            left_mod,
            right_mod,
            correct_iter,
            exclude_iter
        };

        for _ in 0..ret.pos 
        {
            if let Some(dep) = ret.get_normalized_depth()
            {
                ret.feed(dep);
            }
        }

        return ret;
    }
}

impl <'a, DM : DepthModel + Sized> Iterator for FrontendIter<'a, DM> where DM::Input : From<f64>
{
    type Item = Event<DM>;
    fn next(&mut self) -> Option<Self::Item>
    {
        if let Some(dep) = self.get_normalized_depth()
        {
            self.feed(dep);
            self.pos += 1;

            let mut best_score = Default::default(); 
            let mut best_side: Option<Side> = None;
            let mut best_copy_num: Option<u32> = None;

            let mut scan_models = |model : &SVModel<DM>| {
                if best_side.is_some()
                {
                    let cur_score = model.model.get_score();
                    if best_score > cur_score 
                    {
                        best_score = cur_score;
                        best_copy_num = Some(model.copy_num);
                        best_side = Some(model.side);
                    }
                }
                else
                {
                    best_score = model.model.get_score();
                    best_copy_num = Some(model.copy_num);
                    best_side = Some(model.side);
                }
            };

            self.left_mod[0..].iter().for_each(&mut scan_models);
            self.right_mod[0..].iter().for_each(&mut scan_models);

            if best_side.is_some()
            {
                return Some(Event { 
                    score: best_score,
                    pos  : self.pos,
                    side : best_side.unwrap(),
                    copy_num: best_copy_num.unwrap()
                });
            }
        }
        return None;
    }
}

#[cfg(test)]
mod test {
    use crate::scanner::mock_bam::*;
}
