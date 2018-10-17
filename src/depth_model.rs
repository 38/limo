pub trait DepthModel {
    type Input;
    type Output;
    fn put(&mut self, next : Self::Input) -> ();
    fn get_score(&self) -> Self::Output;
}

