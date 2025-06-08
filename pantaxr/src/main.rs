
use pantaxr::cmdline::*;
use pantaxr::fastixe;
use pantaxr::rcls;
use pantaxr::profile;
use pantaxr::stat;
use clap::Parser;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.mode {
        Mode::Profile(profile_args) => profile::profile(profile_args),
        Mode::Fastixe(fastixe_args) => fastixe::fastixe(fastixe_args),
        Mode::Rcls(rcls_args) => rcls::rcls(rcls_args),
        Mode::Stat(stat_args) => stat::stat(stat_args),
    }

}
