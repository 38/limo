mod hts;
mod window;
mod scanner;

use self::scanner::Scanner;
/*use self::hts::{hts_open, hts_close, sam_hdr_read, sam_index_load, bam_read1, bam_init1};
use std::ffi::CString;*/
fn main() 
{
    let scanner = Scanner::new("/home/haohou/base2/lumpy-sv-2/lumpy_tests/rice/AL87.discordant.sort.bam", 1, None);

    println!("{:?}", scanner.is_ok());
    /*
    let fp = unsafe{hts_open(CString::new("/home/haohou/base2/lumpy-sv-2/lumpy_tests/rice/AL87.discordant.sort.bam").unwrap().as_ptr(),
                             CString::new("rb").unwrap().as_ptr())};

    println!("{:?}", fp);

    let hdr  = unsafe{sam_hdr_read(fp)};

    let read = unsafe{bam_init1()};

    println!("{}", unsafe{bam_read1((*fp).fp.bgzf, read)});

    println!("Pos = {}", unsafe{(*read).core.pos});
    
    let read = unsafe{bam_init1()};

    println!("{}", unsafe{bam_read1((*fp).fp.bgzf, read)});

    println!("Pos = {}", unsafe{(*read).core.pos});

    unsafe{ hts_close(fp) };
    */
}
