mod hts;
mod window;
mod bamfile;

use self::bamfile::BamFile;
/*use self::hts::{hts_open, hts_close, sam_hdr_read, sam_index_load, bam_read1, bam_init1};
use std::ffi::CString;*/
fn main() -> Result<(), ()>
{
    let bam = BamFile::new("/home/haohou/base2/lumpy-sv-2/lumpy_tests/rice/AL87.discordant.sort.bam", 1, None)?;

    for al in bam.try_iter()? 
    {
       //println!("{} {} {} {} {}", al.begin(), al.end(), al.length(), al.mqual(), al.is_split_read()); 
       //println!("{:?}", al.cigar(0));
       let alignment = al.alignment();

       for map in alignment 
       {
           println!("{:?}", map);
       }
    }
    Ok(())
}
