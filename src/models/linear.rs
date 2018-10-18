use crate::depth_model::DepthModel;

pub struct LinearModel;

impl LinearModel {
    #[allow(dead_code)]
    pub fn new() -> LinearModel
    {
        LinearModel
    }
}

impl DepthModel for LinearModel {
    type Input = f64;
    type Output = f64;
    fn put(&mut self, _next: Self::Input) 
    {
    }

    fn get_score(&self) -> Self::Output 
    {
        return 0f64;
    }
}
