use crate::depth_model::DepthModel;
use crate::frontend::{Event, FrontendIter, Frontend, Side};

pub struct EventPairProc<'a, DM : DepthModel> {
    left_side: Vec<Option<Event<'a, DM>>>,
    last_pos : u32,
    recent   : Vec<Event<'a, DM>>,
    max_copy_num : u32,
    fe_iter  : FrontendIter<'a, DM>
}

impl <'a, DM : DepthModel> EventPairProc<'a, DM> 
{
    pub fn new(fe:&'a mut Frontend<DM>) -> Self
    {
        let target_copy_nums = fe.get_copy_nums();
        let max_copy_num = target_copy_nums.iter().fold(0, |a,b| std::cmp::max(a,*b));
        return Self {
            left_side : vec![None; (max_copy_num + 1) as usize],
            last_pos  : 0,
            recent    : Vec::new(),
            max_copy_num,
            fe_iter   : fe.iter()
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
                               /* TODO: filter on the average exclude-rate and len */
                               if self.left_side[cur_cn].is_none() ||
                                  self.left_side[cur_cn].as_ref().unwrap().score > current.score
                               {
                                   self.left_side[cur_cn] = Some(Clone::clone(current));
                               }
                           },
                           Side::Right => {
                               if let Some(ref left_side) = self.left_side[cur_cn]
                               {
                                   if DM::score_threshold(Clone::clone(&left_side.score)) && 
                                      DM::score_threshold(Clone::clone(&current.score))
                                   {
                                       ret = Some((Clone::clone(left_side), Clone::clone(current)));
                                   }
                               }
                               if ret.is_some()
                               {
                                   self.left_side[cur_cn] = None;
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