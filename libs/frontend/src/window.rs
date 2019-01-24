use std::ops::{Add, Sub};
use std::io::{Read, Write};
use std::slice;

#[allow(dead_code)]
pub struct Window<T> where
    T : Sized + Default + Clone,
    T : Add<Output = T>,
    T : Sub<Output = T>
{
    acc: Vec<T>,
    ext: Vec<T>
}


#[allow(dead_code)]
impl <T> Window<T> where
    T : Sized + Default + Clone,
    T : Add<Output = T>,
    T : Sub<Output = T>
{
    pub fn try_dump<TOut:Write>(&self, fp:&mut TOut) -> Result<(), std::io::Error>
    {
        fp.write_all(unsafe { slice::from_raw_parts((&self.acc.len() as *const usize) as *const u8, std::mem::size_of::<usize>()) })?;
        fp.write_all(unsafe { slice::from_raw_parts(self.acc.as_ptr() as *const u8, std::mem::size_of::<T>() * self.acc.len()) })?;
        fp.write_all(unsafe { slice::from_raw_parts(self.ext.as_ptr() as *const u8, std::mem::size_of::<T>() * self.acc.len()) })?;

        return Ok(());
    }

    pub fn try_load<TIn:Read>(fp:&mut TIn) -> Result<Self, std::io::Error>
    {
        let mut buffer = 0usize;
        let buffer_ptr = &mut buffer as *mut usize;

        let size = fp.read(unsafe{ slice::from_raw_parts_mut(buffer_ptr as *mut u8, std::mem::size_of::<usize>())})?;

        if size != std::mem::size_of::<usize>()
        { 
            return Err(std::io::Error::new(std::io::ErrorKind::Other, "unexpected header size")); 
        }

        let mut ret = Self {
            acc: vec![T::default(); buffer],
            ext: vec![T::default(); buffer]
        };

        let acc = fp.read(unsafe{slice::from_raw_parts_mut(ret.acc.as_mut_ptr() as *mut u8, std::mem::size_of::<T>() * buffer)})?;
        let ext = fp.read(unsafe{slice::from_raw_parts_mut(ret.ext.as_mut_ptr() as *mut u8, std::mem::size_of::<T>() * buffer)})?;

        if acc != buffer * std::mem::size_of::<T>() ||
           ext != buffer * std::mem::size_of::<T>()
        {
            return Err(std::io::Error::new(std::io::ErrorKind::Other, "cannot read the data block"));
        }

        return Ok(ret);
    }

    pub fn new(range:usize) -> Self 
    {
        Self {
            acc: vec![T::default(); range + 1],
            ext: vec![T::default(); range + 1]
        }
    }

    pub fn accumulate(&mut self, beg:usize, mut end:usize, weight:T) 
    {
        end = (self.acc.len() - 1).min(end);
        self.acc[beg] = T::clone(&self.acc[beg]) + T::clone(&weight);
        self.acc[end] = T::clone(&self.acc[end]) - T::clone(&weight);

        if end - beg > 1 
        {
            self.ext[beg + 1] = T::clone(&self.ext[beg + 1]) + T::clone(&weight);
            self.ext[end] = T::clone(&self.ext[end]) - T::clone(&weight);
        }
    }

    pub fn iter<'a, R>(&'a self, win_size : usize) -> WindowIter<'a, T, R> where
        R : Default + Clone + Add<Output = R> + Sub<Output = R> + From<T> + PartialOrd
    {
        WindowIter::new(self, win_size)
    }

    pub fn size(&self) -> usize {
        self.acc.len() - 1
    }
}

pub struct WindowIter<'a, T, R : From<T>> where
    T : Sized + Default + Clone,
    T : Add<Output = T>,
    T : Sub<Output = T>,
    R : Default + Clone + Add<Output = R> + Sub<Output = R>
{
    win_size : usize,
    win_obj  : &'a Window<T>,
    i        : usize,
    a_window : Vec<R>,
    e_window : Vec<R>,
    a_sum    : R,
    e_sum    : R,
    result   : R
}

impl <'a, T, R : From<T>> WindowIter<'a, T, R> where
    T : Sized + Default + Clone,
    T : Add<Output = T>,
    T : Sub<Output = T>,
    R : Default + Clone + Add<Output = R> + Sub<Output = R> + PartialOrd
{
    pub fn new(win : &'a Window<T>, size: usize) -> Self 
    {
        Self {
            win_size : size,
            win_obj  : win,
            i        : 0,
            a_window : vec![R::default(); size],
            e_window : vec![R::default(); if size > 0 { size - 1 } else { 0 } ],
            a_sum    : R::default(),
            e_sum    : R::default(),
            result   : R::default()
        }
    }

    pub fn get_next(&mut self) -> Option<R> 
    {
        if self.i >= self.win_obj.acc.len() - 1
        {
            return None;
        }

        if self.win_size == 1 
        {
            self.result = R::clone(&self.result) + R::from(T::clone(&self.win_obj.acc[self.i]));
            self.i += 1;
            assert!(self.result >= R::default());
            return Some(R::clone(&self.result));
        }

        while self.i <= self.win_obj.acc.len() - 1
        {
            self.a_sum = R::clone(&self.a_sum) + R::from(T::clone(&self.win_obj.acc[self.i]));
            self.e_sum = R::clone(&self.e_sum) + R::from(T::clone(&self.win_obj.ext[self.i]));

            let mut should_return = false;
            let mut ret = R::default();

            if self.i >= self.win_size
            {
                should_return = true;
                ret =  R::clone(&self.result);
                self.result = R::clone(&self.result) - R::clone(&self.a_window[self.i % self.win_size]);

            }

            if self.i >= self.win_size - 1 
            {
                self.result = R::clone(&self.result) + R::clone(&self.e_window[self.i % (self.win_size - 1)]);
            }

            self.a_window[self.i % self.win_size] = R::clone(&self.a_sum);
            self.e_window[self.i % (self.win_size - 1)] = R::clone(&self.e_sum);

            self.result = R::clone(&self.result) + R::clone(&self.a_sum) - R::clone(&self.e_sum);
            
            self.i += 1;

            if should_return 
            {
                assert!(ret >= R::default());
                return Some(ret);
            }
        }

        return None;
    }
}

impl <'a, T, R : From<T>> Iterator for WindowIter<'a, T, R> where
    T : Sized + Default + Clone,
    T : Add<Output = T>,
    T : Sub<Output = T>,
    R : Default + Clone + Add<Output = R> + Sub<Output = R> + PartialOrd
{
    type Item = R;
    fn next(&mut self) -> Option<Self::Item>
    {
        self.get_next()
    }
}



#[cfg(test)]
mod window_test {
    use super::*;
    #[test]
    fn window_iter() {
        let mut win = Window::<i32>::new(15);
        win.accumulate(1,5,1);
        win.accumulate(2,6,1);
        win.accumulate(3,9,1);
        win.accumulate(2,9,1);
        win.accumulate(4,6,1);

        let result : Vec<i32> = win.iter(3).collect();

        let exp = [3, 4, 5, 5, 5, 4, 2, 2, 2, 0, 0, 0];

        for (actual, expected) in result.into_iter().zip(exp.into_iter())
        {
            assert_eq!(actual, *expected);
        }
        
        let result1 : Vec<i64> = win.iter(1).collect();

        let exp1 = [1, 3, 4, 5, 4, 2, 2, 2, 0, 0, 0];

        for (actual, expected) in result1.into_iter().zip(exp1.into_iter())
        {
            assert_eq!(actual, *expected);
        }

    }
}
