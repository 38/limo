use std::ops::{Add, Sub};

pub struct Window<T : Add + Sub + Sized> {
    acc: Vec<T>,
    ext: Vec<T>
}
