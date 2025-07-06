fn main() {
    // // Link to libdeflate
    // println!("cargo:rustc-link-lib=deflate");

    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        println!("cargo:rustc-link-search=native={}/lib", prefix);
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}/lib", prefix);
    }
}
