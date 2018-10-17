extern crate bindgen;

use std::process::Command;
use std::path::Path;

fn create_hts_bindings() -> Result<(), ()>
{
    if !Path::new("generated/htslib.rs").exists()
    {
        bindgen::Builder::default()
            .header("include/htslib_wrap.h")
            .clang_arg("-Ihtslib/htslib")
            .layout_tests(false)
            .generate()?
            .write_to_file("generated/htslib.rs")
            .expect("Unable to write the generated file");
    }
    Ok(())
}
fn main() -> Result<(), ()>
{
    create_hts_bindings()?;
    if !Path::new("htslib/libhts.a").exists()
    {
        Command::new("make").current_dir("htslib").spawn().expect("Unable to call makefile for htslib");
    }
    println!("cargo:rustc-link-search=htslib/");
    println!("cargo:rustc-link-lib=hts");
    println!("cargo:rustc-link-lib=z");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:rustc-link-lib=bz2");
    println!("cargo:rustc-link-lib=lzma");
    println!("cargo:rustc-link-lib=pthread");
    Ok(())
}
