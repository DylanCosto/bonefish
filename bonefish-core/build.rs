// With `--features gpu`, compile cuda/gpu_dp.cu with nvcc (from $CUDA_HOME, default
// ~/micromamba/envs/cuda) for sm_120 and link cudart. Without the feature this does nothing,
// so the default and wasm builds need no CUDA toolchain.
fn main() {
    println!("cargo:rerun-if-changed=cuda/gpu_dp.cu");
    println!("cargo:rerun-if-env-changed=CUDA_HOME");
    if std::env::var_os("CARGO_FEATURE_GPU").is_none() {
        return;
    }
    let cuda = std::env::var("CUDA_HOME")
        .unwrap_or_else(|_| format!("{}/micromamba/envs/cuda", std::env::var("HOME").unwrap()));
    let out = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let obj = out.join("gpu_dp.o");
    let st = std::process::Command::new(format!("{cuda}/bin/nvcc"))
        .args(["-arch=sm_120", "-O3", "-Xcompiler", "-fPIC", "-c", "cuda/gpu_dp.cu", "-o"])
        .arg(&obj)
        .status()
        .expect("run nvcc (set CUDA_HOME)");
    assert!(st.success(), "nvcc failed");
    let st = std::process::Command::new("ar").arg("crs").arg(out.join("libgpu_dp.a")).arg(&obj).status().expect("run ar");
    assert!(st.success(), "ar failed");
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=gpu_dp");
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-link-search=native={cuda}/lib");
    println!("cargo:rustc-link-lib=dylib=cudart");
}
