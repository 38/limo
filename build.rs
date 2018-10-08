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
            .generate()?
            .write_to_file("generated/htslib.rs").
            expect("Unable to write the generated file");
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
    /*if let Ok(search_path) = env::var("PSTD_LIB_PATH")
    {
        println!("cargo:rustc-link-search={}", search_path);
    }
    else if let Some(lib_path) = search_for_pstd_lib() 
    {
        println!("cargo:rustc-link-search={}", lib_path);
    }
    else
    {
        eprintln!("Warn: Cannot find libpstd.so, plumber-rs crate could not be built.");
        eprintln!("Hint: make sure you have plumber intalled or try to set PSTD_LIB_PATH environment variable");
    }
    println!("cargo:rustc-link-lib=pstd");*/
}
