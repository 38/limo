use std::fmt::Debug;

pub trait DepthModel {
    type Input : Copy + From<f64>;
    type Output: PartialOrd + Default + Debug;
    type ParamType : Copy;
    fn create_model(copy_num:u32, left_side:bool, param:Self::ParamType) -> Self;
    fn put(&mut self, next : Self::Input) -> ();
    fn get_score(&self) -> Self::Output;
}
