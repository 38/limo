use super::hts::*; 
use std::ffi::{CString, c_void};
use std::ptr::{null_mut, null};
use std::ops::Index;

#[allow(dead_code)]
pub struct Alignment<'a> {
    data : &'a bam1_t
}

#[derive(Debug)]
pub enum Nucleotide { A, C, G, T, N }

static A : Nucleotide = Nucleotide::A;
static C : Nucleotide = Nucleotide::C;
static G : Nucleotide = Nucleotide::G;
static T : Nucleotide = Nucleotide::T;
static N : Nucleotide = Nucleotide::N;

#[derive(Debug)]
pub enum Cigar {
    Match(u32),
    Insert(u32),
    Delete(u32),
    Skip(u32),
    Soft(u32),
    Hard(u32),
    Pad(u32),
    Equal(u32),
    Diff(u32),
    Back(u32)
}

impl Cigar {
    fn from_alignment(al:&Alignment, idx:usize) -> Option<Cigar> 
    {
        if idx > (al.data.core.n_cigar as usize) { return None } 
        let cigar_num : u32 = unsafe { *(al.data.data.offset(al.data.core.l_qname as isize + idx as isize) as *const u32)};

        let cigar_op = cigar_num & BAM_CIGAR_MASK;
        let cigar_len = cigar_num >> BAM_CIGAR_SHIFT;

        match cigar_op {
            BAM_CMATCH => Some(Cigar::Match(cigar_len)),
            BAM_CINS => Some(Cigar::Insert(cigar_len)),
            BAM_CDEL => Some(Cigar::Delete(cigar_len)),
            BAM_CREF_SKIP => Some(Cigar::Skip(cigar_len)),
            BAM_CSOFT_CLIP => Some(Cigar::Soft(cigar_len)),
            BAM_CHARD_CLIP => Some(Cigar::Hard(cigar_len)),
            BAM_CPAD => Some(Cigar::Pad(cigar_len)),
            BAM_CEQUAL => Some(Cigar::Equal(cigar_len)),
            BAM_CDIFF => Some(Cigar::Diff(cigar_len)),
            BAM_CBACK => Some(Cigar::Back(cigar_len)),
            _ => None
        }
    }

    fn in_alignment(&self) -> bool 
    {
        match self {
            Cigar::Match(_) | Cigar::Insert(_) | Cigar::Soft(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
            _ => false
        }
    }
    
    fn in_reference(&self) -> bool 
    {
        match self {
            Cigar::Match(_) | Cigar::Delete(_) | Cigar::Skip(_) | Cigar::Equal(_) | Cigar::Diff(_) => true,
            _ => false
        }
    }
}

impl From<u32> for &'static Nucleotide {
    fn from(what:u32) -> Self 
    {
        match what {
            1 => &A,
            2 => &C,
            4 => &G,
            8 => &T,
            _ => &N
        }
    }
}

impl Into<char> for Nucleotide {
    fn into(self) -> char 
    {
        match self {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'T',
            Nucleotide::N => 'N'
        }
    }
}

pub struct Sequence<'a> {
    alignment : &'a Alignment<'a>
}

impl <'a> Sequence<'a> {
    pub fn size(&self) -> usize { self.alignment.length() as usize }
}

impl <'a> Index<usize> for Sequence<'a> {
    type Output = Nucleotide;
    fn index(&self, offset : usize) -> &Nucleotide 
    {
        let hts_obj = self.alignment.data;
        let seq = unsafe { hts_obj.data.offset(hts_obj.core.n_cigar as isize * 4 + (hts_obj.core.l_qname as isize) + offset as isize / 2) };
        let numeric : u32 = if offset % 2 == 0 { (unsafe { *seq } as u32) >> 4 }  else  { (unsafe { * seq } as u32) & 0xf };
        return From::from(numeric);
    }
}

#[allow(dead_code)]
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
    pub fn sequence(&self) -> Sequence { Sequence { alignment : self } }
    pub fn cigar(&self, idx:usize) -> Option<Cigar> { Cigar::from_alignment(self, idx) }
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
    pub fn try_iter(&self) -> Result<BamFileIter, ()> { self.try_iter_range(0, self.length) }

    #[allow(dead_code)]
    pub fn try_iter_range(&self, begin:usize, end:usize) -> Result<BamFileIter, () >
    {
        let left = begin;

        let right = if end > self.length { self.length } else { end };

        let iter = unsafe{ sam_itr_queryi(self.idx, self.chrom_id as i32, left as i32,  right as i32) };

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
