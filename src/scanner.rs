use super::window::Window;
use super::bamfile::{BamFile, Alignment, BamFileIter};

pub trait AlignmentType {
    fn get_begin(&self) -> u32;
    fn get_end(&self) -> u32;
    fn get_length(&self) -> u32;
    fn check_is_split_read(&self) -> bool;
    fn get_mqual(&self) -> u32;
}

impl <'a> AlignmentType for Alignment<'a> {
    fn get_begin(&self) -> u32 { self.begin() }
    fn get_end(&self) -> u32 { self.end() }
    fn get_length(&self) -> u32 { self.length() }
    fn check_is_split_read(&self) -> bool {self.is_split_read() }
    fn get_mqual(&self) -> u32 {self.mqual() }
}

pub trait Input<'a, T:AlignmentType> {
    type IterType : Iterator<Item = T>;
    fn size(&self) -> usize;
    fn try_iter(&'a self) -> Result<Self::IterType, ()>;
}

impl <'a> Input<'a, Alignment<'a>> for BamFile {
    type IterType = BamFileIter<'a>;
    fn size(&self) -> usize { self.size() }
    fn try_iter(&'a self) -> Result<Self::IterType, ()> { self.try_iter() }
}

pub struct Scanner {
    corrected_window : Window<i32>,
    low_mq_window    : Window<i32>,
    raw_window       : Window<i32>,
    common_read_len     : u32,
    common_read_len_cnt : u32
}

impl Scanner {
    #[allow(dead_code)]
    pub fn get_common_read_length(&self) -> u32 
    {
        self.common_read_len
    }

    #[allow(dead_code)]
    pub fn get_corrected(&self) -> &Window<i32>
    {
        &self.corrected_window
    }

    #[allow(dead_code)]
    pub fn get_low_mq_window(&self) -> &Window<i32>
    {
        &self.low_mq_window
    }

    #[allow(dead_code)]
    pub fn get_raw_window(&self) -> &Window<i32>
    {
        &self.raw_window
    }

    pub fn new<'a, IType, AType>(bam:&'a IType) -> Result<Scanner, ()>
        where AType : AlignmentType,
              IType : Input<'a, AType>
    {
        let size = bam.size();
        
        let mut ret = Scanner {
            corrected_window : Window::<i32>::new(size),
            low_mq_window    : Window::<i32>::new(size),
            raw_window       : Window::<i32>::new(size),
            common_read_len  : 0,
            common_read_len_cnt: 0
        };

        for read in bam.try_iter()?
        {
           if ret.common_read_len != read.get_length() 
           {
               if ret.common_read_len_cnt > 0 
               {
                   ret.common_read_len_cnt -= 1
               }
               else
               {
                   ret.common_read_len_cnt += 1;
               }

               if ret.common_read_len_cnt == 0
               {
                   ret.common_read_len = read.get_length();
               }
           }

           let (begin, end) = (read.get_begin() as usize, read.get_end() as usize);

           ret.raw_window.accumulate(begin, end, 1);

           if read.check_is_split_read() { continue; }

           if read.get_mqual() == 0 
           {
               ret.low_mq_window.accumulate(begin, end, 1);
           }

           ret.corrected_window.accumulate(begin, end, 1);
        }

        return Ok(ret);
    }
}

#[cfg(test)]
mod scanner_test {
    use super::*;

    struct TestAlignment {
        begin : u32,
        end   : u32,
        split : bool,
        qual  : u32
    }

    impl <'a> AlignmentType for &'a TestAlignment {
        fn get_begin(&self) -> u32 { self.begin }
        fn get_end(&self) -> u32 { self.end }
        fn get_length(&self) -> u32 { (self.end - self.begin) as u32 }
        fn check_is_split_read(&self) -> bool { self.split }
        fn get_mqual(&self) -> u32 { self.qual }
    }

    impl <'a> Input<'a, &'a TestAlignment> for (usize, Vec<TestAlignment>) {
        type IterType = std::slice::Iter<'a, TestAlignment>;
        fn size(&self) -> usize { self.0 }
        fn try_iter(&'a self) -> Result<Self::IterType, ()> { Ok(self.1[0..].iter()) }
    }

    #[test]
    fn test_scanner() -> Result<(), ()>
    {
        let my_bam = (10, vec![ 
            TestAlignment{begin: 1, end: 6, split: false, qual: 100 },  // ==> 5
            TestAlignment{begin: 1, end: 5, split: false, qual: 100 },  // ==> 4
            TestAlignment{begin: 2, end: 4, split: false, qual: 100 },  // ==> 2
            TestAlignment{begin: 2, end: 7, split: false, qual: 100 },  // ==> 5
            TestAlignment{begin: 3, end: 8, split: false, qual: 100 },  // ==> 5
            TestAlignment{begin: 0, end: 10, split: true, qual: 100},   // ==> 10
            TestAlignment{begin: 0, end: 5, split: false, qual:0}       // ==> 5
        ]);

        // 01234567890
        //  aaaaa
        //  bbbb
        //   cc
        //   ddddd
        //    eeeee
        // 245543210000

        let scanner = Scanner::new(&my_bam)?;

        assert_eq!(scanner.get_common_read_length(), 5);

        eprintln!("{:?}", scanner.get_corrected().iter::<i32>(2).collect::<Vec<i32>>());

        assert_eq!(scanner.get_corrected().iter::<i32>(2).collect::<Vec<i32>>(),     vec![3,5,6,6,5,3,2,1]);
        assert_eq!(scanner.get_raw_window().iter::<i32>(2).collect::<Vec<i32>>(),    vec![4,6,7,7,6,4,3,2]);
        assert_eq!(scanner.get_low_mq_window().iter::<i32>(2).collect::<Vec<i32>>(), vec![1,1,1,1,1,0,0,0]);

        return Ok(());
    }
}
