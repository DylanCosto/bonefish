// With `--features gpu`, embed an rpath to the CUDA runtime (link args from a dependency's
// build script do not propagate to this binary).
fn main() {
    println!("cargo:rerun-if-env-changed=CUDA_HOME");
    if std::env::var_os("CARGO_FEATURE_GPU").is_none() {
        return;
    }
    let cuda = std::env::var("CUDA_HOME")
        .unwrap_or_else(|_| format!("{}/micromamba/envs/cuda", std::env::var("HOME").unwrap()));
    println!("cargo:rustc-link-arg=-Wl,-rpath,{cuda}/lib");
}
