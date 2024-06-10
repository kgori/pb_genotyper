use clap::error::ErrorKind;
use clap::{CommandFactory, Parser};
use std::path::{Path, PathBuf};

#[derive(Parser, Debug, Clone)]
#[command(version, about, long_about = None)]
pub struct ProgramOptions {
    #[arg(short, long)]
    pub variants: PathBuf,
    
    #[arg(short, long)]
    pub regions: PathBuf,

    #[arg(short, long)]
    pub bamfile: PathBuf,

    #[arg(short, long)]
    pub output: PathBuf,
}

fn validate_file(file: &Path) {
    if !file.exists() {
        let mut cmd = ProgramOptions::command();
        cmd.error(
            ErrorKind::ValueValidation,
            format!("file `{}` not found", file.display()),
        )
        .exit();
    }
}

pub fn parse_cli() -> ProgramOptions {
    let args = ProgramOptions::parse();
    validate_file(&args.variants);
    validate_file(&args.regions);
    validate_file(&args.bamfile);
    args
}
