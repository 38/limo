use std::fmt::Debug;

pub trait DepthModel where Self:Clone {
    type Input : Copy + From<f64>;
    type Output: PartialOrd + Default + Debug;
    type ParamType : Copy;
    const EPS : Self::Output;
    fn create_model(copy_num:u32, left_side:bool, param:Self::ParamType) -> Self;
    fn put(&mut self, next : Self::Input) -> ();
    fn get_score(&self) -> Self::Output;
}
