use super::hts::*; 
use std::ffi::CString;
use std::ptr::null_mut;

pub struct Scanner {
    chrom : String,
    fp    : *mut htsFile,
    hdr   : *mut bam_hdr_t,
    idx   : *mut hts_idx_t,
    length: usize,
    iter  : *mut hts_itr_t,
    common_rl_cnt : u32,
    common_rl : u32
}

impl Scanner {
    pub fn new<'c,'b>(path:&'c str, chrom:u32, reference:Option<&'b str>) -> Result<Self, ()>
    {
        let fp = unsafe {
            hts_open(CString::new(path).unwrap().as_ptr(), 
                     CString::new("rb").unwrap().as_ptr()) 
        };

        if fp == null_mut() 
        {
            return Err(());
        }

        if reference.is_some() 
        {
            if unsafe{ hts_set_fai_filename(fp, CString::new(reference.unwrap()).unwrap().as_ptr()) } < 0
            {
                return Err(());
            }
        }

        let hdr = unsafe { sam_hdr_read(fp) };

        if hdr == null_mut()
        {
            return Err(());
        }

        let idx = unsafe { sam_index_load(fp, CString::new(path).unwrap().as_ptr()) };

        if idx == null_mut()
        {
            return Err(());
        }

        let len = unsafe{ (*(*hdr).target_len.offset(chrom as isize)) };

        let iter = unsafe{ sam_itr_queryi(idx, chrom as i32, 0,  len as i32) };
        
        return Ok(Scanner{
            chrom : unsafe{ CString::from_raw(*(*hdr).target_name.offset(chrom as isize)).into_string().unwrap() },
            fp    : fp,
            hdr   : hdr,
            idx   : idx,
            length: len as usize,
            iter  : iter,
            common_rl_cnt: 0,
            common_rl: 0
        });
    }
}
