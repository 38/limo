use crate::depth_model::DepthModel;

#[derive(Debug, Clone)]
pub struct LinearModel {
    k : f64,
    p : f64,
    pe: f64,
    buf: Vec<f64>,
    sum: f64,
    square: f64,
    left : u32,
    right: u32,
    left_idx: u32,
    len : u32
}

impl LinearModel {
    #[allow(dead_code)]
    pub fn new(len:u32, from:f64, to:f64) -> LinearModel
    {
        return LinearModel {
            k: (to - from) / ((len - 1) as f64),
            p: from,
            pe: to,
            buf: vec![0f64;len as usize],
            sum: 0f64,
            square: 0f64,
            left: 0,
            right: 0,
            left_idx: 0,
            len: len
        };
    }
}

fn square<T:std::ops::Mul + Copy>(what:T) -> T::Output
{
    return what * what;
}

impl DepthModel for LinearModel {
    type Input = f64;
    type Output = f64;
    type ParamType = u32;
    const EPS:f64 = 1e-5;
    fn create_model(copy_num:u32, left: bool, p:u32) -> Self
    {
        let target_depth = 0.5 * (copy_num as f64);
        if left
        {
            return LinearModel::new(p + 1, 1.0, target_depth);
        }

        return LinearModel::new(p + 1, target_depth, 1.0);
    }
    fn put(&mut self, next: Self::Input) 
    {
        if self.right - self.left >= self.len
        {
            self.square += 2.0 * self.k * self.sum + square(next - self.pe) - square(self.buf[self.left_idx as usize] - self.p + self.k) + (self.len as f64 * square(self.k));
            if self.square < 0.0 { self.square = 0f64 }
            
            self.sum += next - self.buf[self.left_idx as usize];

            self.buf[self.left_idx as usize] = next;

            self.left_idx = if self.left_idx + 1 == self.len { 0 } else { self.left_idx + 1};

            self.left += 1;
            self.right += 1;
        }
        else
        {
            self.buf[self.right as usize] = next;
            self.right += 1;

            self.square += square(next - self.p - self.k * (self.right as f64) + self.k);
            self.sum += next - self.p - self.k * (self.right as f64) + self.k;
        }
    }

    fn get_score(&self) -> Self::Output 
    {
        return square(self.sum) + self.square * (self.len as f64);
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_linear_function() -> Result<(), ()>
    {
        let mut lm = LinearModel::new(11, 0.0, 10.0); 

        assert_eq!(lm.k, 1.0);
        assert_eq!(lm.p, 0.0);
        assert_eq!(lm.pe, 10.0);

        (0..=10).for_each(|_| {
            lm.put(0.0)
        });

        assert_eq!(lm.sum, -55.0);
        assert_eq!(lm.square, 385.0);

        lm.put(10.0);
        
        assert_eq!(lm.sum, -45.0);
        assert_eq!(lm.square, 285.0);

        lm.put(10.0);

        assert_eq!(lm.sum, -35.0);
        assert_eq!(lm.square, 205.0);

        return Ok(());
    }

}
