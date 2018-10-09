use super::hts::*; 
use std::ffi::{CString, c_void};
use std::ptr::{null_mut, null};

#[allow(dead_code)]
pub struct Alignment<'a> {
    data : &'a bam1_t
}

impl <'a> Alignment<'a> {
    pub fn begin(&self) -> u32 { self.data.core.pos as u32 }
    pub fn end(&self) -> u32 { (self.data.core.pos + self.data.core.l_qseq) as u32 }
    pub fn length(&self) -> u32 {  self.data.core.l_qseq as u32 }
    pub fn mqual(&self) -> u32 { self.data.core.qual as u32 }
    pub unsafe fn read_tag(&self, tag : &str) -> *const u8 
    {
        let tag_array = CString::new(tag).unwrap().as_ptr();
        bam_aux_get(self.data as *const bam1_t, tag_array) 
    }
    pub fn is_split_read(&self) -> bool { return unsafe{ self.read_tag("SA") } != null(); }
}

#[allow(dead_code)]
pub struct BamFile {
    chrom         : String,
    chrom_id      : u32,
    fp            : *mut htsFile,
    hdr           : *mut bam_hdr_t,
    idx           : *mut hts_idx_t,
    length        : usize
}


#[allow(dead_code)]
pub struct BamFileIter<'a> {
    iter : *mut hts_itr_t,
    file : &'a BamFile,
    buffer : *mut bam1_t
}

impl <'a> Drop for BamFileIter<'a> {
    fn drop(&mut self)
    {
        if self.iter != null_mut() 
        {
            unsafe { hts_itr_destroy(self.iter) };
            self.iter = null_mut();
        }

        if self.buffer != null_mut()
        {
            unsafe { bam_destroy1(self.buffer) };
            self.buffer = null_mut();
        }
    }
}

impl Drop for  BamFile {
    fn drop(&mut self) 
    {

        if self.idx != null_mut()
        {
            unsafe { hts_idx_destroy(self.idx) };
            self.idx = null_mut();
        }

        if self.fp != null_mut()
        {
            unsafe { hts_close(self.fp) };
            self.fp = null_mut();
        }
    }
}

impl BamFile {
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

        return Ok(BamFile{
            chrom : unsafe{ CString::from_raw(*(*hdr).target_name.offset(chrom as isize)).into_string().unwrap() },
            chrom_id : chrom as u32,
            fp    : fp,
            hdr   : hdr,
            idx   : idx,
            length: len as usize
        });
    }

    #[allow(dead_code)]
    pub fn try_iter(&self) -> Result<BamFileIter, ()>
    {
        let iter = unsafe{ sam_itr_queryi(self.idx, self.chrom_id as i32, 0,  self.length as i32) };

        if iter == null_mut() 
        {
            return Err(());
        }

        let buffer = unsafe{ bam_init1() };

        if buffer == null_mut()
        {
            unsafe { hts_itr_destroy(iter) };
            return Err(());
        }

        return Ok(BamFileIter {
            iter : iter,
            file : self,
            buffer: buffer
        });
    }
}

impl <'a> Iterator for BamFileIter<'a> {
    type Item = Alignment<'a>;

    fn next(&mut self)->Option<Self::Item> 
    {
        let rc = unsafe{ hts_itr_next((*self.file.fp).fp.bgzf, self.iter, self.buffer as *mut c_void, self.file.fp as *mut c_void) };
        if rc < 0
        {
            return None;
        }

        return Some(Alignment { data : unsafe { self.buffer.as_ref().unwrap() } });
    }

}
