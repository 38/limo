use std::fmt::Debug;
use crate::scanner::Scanner;

pub trait DepthModel where Self:Clone {
    type Input : Copy + From<f64>;
    type Output: PartialOrd + Default + Debug + Clone;
    type ParamType : Copy;
    fn determine_default_param(scanner:&Scanner, window_size: u32, copy_nums: &[u32]) -> Self::ParamType;
    fn score_cmp(left : Self::Output, right : Self::Output) -> i32;
    fn score_threshold(win_size:u32, score : Self::Output) -> bool;
    fn create_model(copy_num:u32, left_side:bool, param:Self::ParamType) -> Self;
    fn put(&mut self, next : Self::Input) -> ();
    fn get_score(&self) -> Self::Output;
}
