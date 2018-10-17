pub struct Histogram {
    freq : Vec<u32>,
    average: Option<f64>,
    count: u32
}

impl Histogram {
    pub fn new(max:u32) -> Histogram
    {
        return Histogram {
            freq : vec![0;max as usize + 1usize],
            count: 0,
            average: None
        };
    }

    pub fn add(&mut self, val:u32) -> ()
    {
        self.average = None;
        let actual_val = if val as usize >= self.freq.len() { self.freq.len() - 1 } else { val as usize};
        self.freq[actual_val] += 1;
        self.count += 1;
    }

    pub fn get_average(&mut self) -> f64 
    {
        if self.average.is_none()
        {
            self.average = Some(self.freq.iter().zip(0..).fold(0f64, |s,(a,v)| s + (*a as f64) * (v as f64)) / (self.count as f64));
        }
        return *self.average.as_mut().unwrap();
    }

    pub fn normalize(&mut self, val:u32) -> f64
    {
        if val as usize > self.freq.len()
        {
            return 1f64;
        }
        return (val as f64) / self.get_average();
    }

}
