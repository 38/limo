use super::window::Window;
use super::bamfile::BamFile;

pub struct Scanner {
    corrected_window : Window<u32>,
    low_mq_window    : Window<u32>,
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
    pub fn get_corrected(&self) -> &Window<u32>
    {
        &self.corrected_window
    }

    #[allow(dead_code)]
    pub fn get_low_mq_window(&self) -> &Window<u32>
    {
        &self.low_mq_window
    }

    pub fn new(bam:&BamFile) -> Result<Scanner, ()>
    {
        let size = bam.size();
        
        let mut ret = Scanner {
            corrected_window : Window::<u32>::new(size),
            low_mq_window    : Window::<u32>::new(size),
            common_read_len  : 0,
            common_read_len_cnt: 0
        };

        for read in bam.try_iter()?
        {
           if ret.common_read_len != read.length() 
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
                   ret.common_read_len = read.length();
               }
           }

           let (begin, end) = (read.begin(), read.end());

           if read.is_split_read() { continue; }

           if read.mqual() == 0 
           {
               ret.low_mq_window.accumulate(begin as usize, end as usize, 1);
           }

           ret.corrected_window.accumulate(begin as usize, end as usize, 1);
        }

        return Ok(ret);
    }
}
