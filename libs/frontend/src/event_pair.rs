use crate::depth_model::DepthModel;
use crate::frontend::{Event, FrontendIter, Frontend, Side};


pub struct EventPairProc<'a, DM : DepthModel> {
    left_side: Vec<Option<Event<'a, DM>>>,
    last_pos : u32,
    recent   : Vec<Event<'a, DM>>,
    max_copy_num : u32,
    window_size: u32,
    fe_iter  : FrontendIter<'a, DM>,
    chrom_size: u32,
    last_mb: u32,
}

impl <'a, DM : DepthModel> EventPairProc<'a, DM> 
{
    pub fn new(fe:&'a Frontend<DM>) -> Self
    {
        let target_copy_nums = fe.get_copy_nums();
        let max_copy_num = target_copy_nums.iter().fold(0, |a,b| std::cmp::max(a,*b));
        return Self {
            left_side : vec![None; (max_copy_num + 1) as usize],
            last_pos  : 0,
            recent    : Vec::new(),
            max_copy_num,
            fe_iter   : fe.iter(),
            window_size: fe.get_window_size(),
            chrom_size: fe.get_chrom_size(),
            last_mb: 0,
        };
    }
}

impl <'a, DM : DepthModel> Iterator for EventPairProc<'a, DM>
{
    type Item = (Event<'a, DM>, Event<'a, DM>);

    fn next(&mut self) -> Option<Self::Item>
    {
        return loop {
            let mut ret = None as Option<Self::Item>;
            let next_event = self.fe_iter.next();
            if let Some(next_event) = next_event 
            {
                if next_event.pos / 1000_0000 != self.last_mb {
                    self.last_mb = next_event.pos / 1000_0000;
                    eprintln!("Pairing model events: Chrom {}, {}MB/{}MB", next_event.chrom, self.last_mb * 10, self.chrom_size / 100_0000);
                }

                if next_event.copy_num > self.max_copy_num { continue; }
                if self.recent.len() == 3 { self.recent.remove(0); }
                self.recent.push(next_event);

                if self.recent.len() == 3 &&
                   DM::score_cmp(Clone::clone(&self.recent[1].score), Clone::clone(&self.recent[0].score)) < 0 &&
                   DM::score_cmp(Clone::clone(&self.recent[1].score), Clone::clone(&self.recent[2].score)) < 0
                {
                    let current = &self.recent[1];
                    let cur_cn = current.copy_num as usize;
                    if self.last_pos + 1 != self.recent[2].pos 
                    {
                        match current.side 
                        {
                           Side::Left => {
                               if self.left_side[cur_cn].is_none() ||
                                  self.left_side[cur_cn].as_ref().unwrap().score > current.score
                               {
                                   self.left_side[cur_cn] = Some(Clone::clone(current));
                               }
                           },
                           Side::Right => {
                               if let Some(ref left_side) = self.left_side[cur_cn]
                               {
                                   if DM::score_threshold(self.window_size, Clone::clone(&left_side.score)) && 
                                      DM::score_threshold(self.window_size, Clone::clone(&current.score))
                                   {
                                       ret = Some((Clone::clone(left_side), Clone::clone(current)));
                                   }
                               }
                               if ret.is_some()
                               {
                                   self.left_side[cur_cn] = None;

                                   /* We assume that the SV doesn't overlap, so once we found a
                                    * good candidate, just remove everything that doesn't that good
                                    * and would overlap with current we have */
                                   for i in 0..self.left_side.len()
                                   {
                                       if self.left_side[i].is_some() && self.left_side[i].as_ref().unwrap().score > ret.as_ref().unwrap().0.score
                                       {
                                           self.left_side[i] = None;
                                       }
                                   }
                               }
                           }
                        }
                    }
                    self.last_pos = self.recent[2].pos;
                }

                if ret.is_some() { break ret; }
            }
            else 
            {
                break None;
            }
        };
    }
}
