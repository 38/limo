use std::process::Command;
use std::path::Path;

use bindgen::Builder as BG;

use num_cpus::get as get_num_cpus;

use std::env;

fn create_hts_bindings(base:&str) -> Result<(), ()>
{
    let include_param = format!("-I{}/htslib/htslib", base);
    if !Path::new("generated/htslib.rs").exists()
    {
        BG::default()
            .header("include/htslib_wrap.h")
            .clang_arg(include_param.as_str())
            .layout_tests(false)
            .generate()?
            .write_to_file("generated/htslib.rs")
            .expect("Unable to write the generated file");
    }
    if !Path::new(&format!("{}/htslib/libhts-static.a", base)).exists() {
        std::os::unix::fs::symlink(format!("{}/htslib/libhts.a", base), format!("{}/htslib/libhts-static.a", base))
        .expect("Unable to create the symlink to libhts.a");
    }
    Ok(())
}
fn main() -> Result<(), std::io::Error>
{
    let base = format!("{}/../", env::var("CARGO_MANIFEST_DIR").unwrap());
    let hts_bin_path = format!("{}/htslib/libhts.a", base);
    if let Err(_) = create_hts_bindings(base.as_str()) { return Err(std::io::Error::new(std::io::ErrorKind::Other, "Bindgen failed")); }
    if !Path::new(hts_bin_path.as_str()).exists()
    {
        Command::new("make")
            .arg(format!("-j{}", get_num_cpus()))
            .current_dir(format!("{}/htslib", base))
            .spawn()
            .expect("Unable to call makefile for htslib");
    }
    println!("cargo:rustc-link-search={}/htslib/", base);
    println!("cargo:rustc-link-lib=hts-static");
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=bz2");
    println!("cargo:rustc-link-lib=lzma");
    println!("cargo:rustc-link-lib=pthread");
    Ok(())
}
