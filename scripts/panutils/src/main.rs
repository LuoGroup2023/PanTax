
use panutils::cmdline::*;
use panutils::fastixe;
use panutils::rcls;
use panutils::profile;
use clap::Parser;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.mode {
        Mode::Profile(profile_args) => profile::profile(profile_args),
        Mode::Fastixe(fastixe_args) => fastixe::fastixe(fastixe_args),
        Mode::Rcls(rcls_args) => rcls::rcls(rcls_args),
        
    }

}
