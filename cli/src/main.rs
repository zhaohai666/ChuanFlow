// ChuanFlow CLI (README 第6层): thin automation entry point that locates and
// launches the C++ engine binaries. Deliberately dependency-free so it builds
// offline in CI and on HPC logins.

use std::env;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::{self, Command};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug, PartialEq)]
enum Command2 {
    Version,
    Run { engine: String, case: String },
}

fn parse(args: &[String]) -> Result<Command2, String> {
    match args.first().map(String::as_str) {
        None | Some("help") | Some("--help") => Err(
            "usage: chuanflow version | chuanflow run steady <case.dict>".to_string(),
        ),
        Some("version") => Ok(Command2::Version),
        Some("run") => {
            let engine = match args.get(1).map(String::as_str) {
                Some("steady") => "chflow-1d-steady",
                Some("transient") => return Err("1D transient solver is on the roadmap (M2)".into()),
                other => {
                    return Err(format!("unknown engine '{}': expected steady|transient", other.unwrap_or("<none>")))
                }
            };
            let case = args
                .get(2)
                .ok_or_else(|| "missing <case.dict>".to_string())?;
            Ok(Command2::Run { engine: engine.to_string(), case: case.clone() })
        }
        other => Err(format!("unknown command '{}'", other.unwrap_or("<empty>"))),
    }
}

/// Engine discovery: $CHUANFLOW_BIN_DIR/<engine>, else PATH lookup.
fn engine_path(engine: &str) -> PathBuf {
    match env::var("CHUANFLOW_BIN_DIR") {
        Ok(dir) => PathBuf::from(dir).join(engine),
        Err(_) => PathBuf::from(engine),
    }
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let cmd = match parse(&args) {
        Ok(cmd) => cmd,
        Err(msg) => {
            eprintln!("chuanflow: {msg}");
            process::exit(2);
        }
    };
    let code = match cmd {
        Command2::Version => {
            println!("chuanflow {VERSION}");
            0
        }
        Command2::Run { engine, case } => match Command::new(engine_path(&engine)).arg(&case).status()
        {
            Ok(status) => status.code().unwrap_or(1),
            Err(err) => {
                let _ = writeln!(io::stderr(), "chuanflow: cannot run {engine}: {err}");
                127
            }
        },
    };
    process::exit(code);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_parses() {
        assert_eq!(parse(&["version".to_string()]), Ok(Command2::Version));
    }

    #[test]
    fn run_steady_parses() {
        let got = parse(&[
            "run".to_string(),
            "steady".to_string(),
            "case.dict".to_string(),
        ]);
        assert_eq!(
            got,
            Ok(Command2::Run { engine: "chflow-1d-steady".into(), case: "case.dict".into() })
        );
    }

    #[test]
    fn rejects_bad_usage() {
        assert!(parse(&[]).is_err());
        assert!(parse(&["run".to_string()]).is_err());
        assert!(parse(&["run".to_string(), "mesh".to_string()]).is_err());
        assert!(parse(&["run".to_string(), "steady".to_string()]).is_err());
        assert!(parse(&["transient".to_string(), "case".to_string()]).is_err());
        assert!(parse(&["run".to_string(), "transient".to_string(), "c".to_string()]).is_err());
    }

    #[test]
    fn cli_version_runs() {
        let Some(exe) = option_env!("CARGO_BIN_EXE_chuanflow") else { return };
        let out = Command::new(exe)
            .arg("version")
            .output()
            .expect("chuanflow binary");
        assert!(out.status.success());
        assert!(String::from_utf8_lossy(&out.stdout).contains("chuanflow"));
    }
}
