// ChuanFlow CLI (README 第6层): thin automation entry point that locates and
// launches the C++ engine binaries. Deliberately dependency-free so it builds
// offline in CI and on HPC logins.

use std::env;
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::{self, Command};

const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Engine registry: (CLI key, CMake target / binary name, domain group).
///
/// The binary name must stay identical to the `add_executable(...)` target in
/// the corresponding CMakeLists.txt; the test `table_is_drift_free` pins the
/// whole table so a rename has to be updated here deliberately.
///
/// Group order follows the usage listing:
/// 1D 稳态 / 1D 瞬态 / nodal / 两相 / PVT / 3D / 耦合.
const ENGINES: &[(&str, &str, &str)] = &[
    // 1D 稳态
    ("steady", "chflow-1d-steady", "1D 稳态"),
    // 1D 瞬态
    ("transient", "chflow-1d-transient", "1D 瞬态"),
    // nodal 节点分析
    ("nodal", "chflow-1d-nodal", "nodal 节点分析"),
    // 两相
    ("twofluid", "chflow-1d-twofluid", "两相"),
    // PVT
    ("pvt", "chflow-pvt", "PVT 流体相态"),
    // 3D 族（src/3d/apps/CMakeLists.txt 及其可选 hook 中的全部 add_executable 目标）
    ("3d-scalar-transport", "chflow-scalar-transport", "3D"),
    ("3d-simple", "chflow-simple-foam", "3D"),
    ("3d-piso", "chflow-piso-foam", "3D"),
    ("3d-vof", "chflow-vof", "3D"),
    ("3d-compressible", "chflow-3d-compressible", "3D"),
    ("3d-dpm", "chflow-3d-dpm", "3D"),
    ("3d-p1", "chflow-3d-p1", "3D"),
    ("3d-combustion", "chflow-3d-combustion", "3D"),
    ("3d-combustion-detailed", "chflow-3d-combustion-detailed", "3D"),
    ("3d-qmom", "chflow-3d-qmom", "3D"),
    ("3d-mixture-pbm", "chflow-3d-mixture-pbm", "3D"),
    ("3d-cht", "chflow-3d-cht", "3D"),
    ("3d-ale", "chflow-3d-ale", "3D"),
    ("3d-les", "chflow-3d-les", "3D"),
    ("3d-rsm", "chflow-3d-rsm", "3D"),
    ("3d-melting", "chflow-3d-melting", "3D"),
    ("3d-cavitation", "chflow-3d-cavitation", "3D"),
    ("3d-radiation-do", "chflow-3d-radiation-do", "3D"),
    ("3d-reconstruction", "chflow-reconstruction", "3D"),
    ("3d-parallel", "chflow-3d-parallel", "3D"),
    ("3d-dynamic-mesh", "chflow-3d-dynamic-mesh", "3D"),
    ("3d-adaptive-refine", "chflow-3d-adaptive-refine", "3D"),
    ("3d-mesh", "chflow-openfoam-mesh", "3D"),
    ("3d-backend", "chflow-3d-backend", "3D"),
    // MPI / PETSc 条件编译目标：仅在工具链具备对应后端时构建；缺二进制时 CLI 如实报 127。
    ("3d-mpi-linear-solver", "chflow-mpi-linear-solver", "3D (MPI 条件编译)"),
    ("3d-petsc-distributed", "chflow-3d-petsc-distributed", "3D (MPI+PETSc 条件编译)"),
    // 耦合 1D/3D
    ("coupled", "chflow-coupled-1d3d", "耦合"),
    ("coupled-transient", "chflow-coupled-1d3d-transient", "耦合"),
    ("coupled-multi", "chflow-coupled-1d3d-multi", "耦合"),
    ("coupled-transient-multi", "chflow-coupled-1d3d-transient-multi", "耦合"),
];

#[derive(Debug, PartialEq)]
enum Command2 {
    Version,
    Run { engine: String, case: String },
}

/// Parsed outcome. `Help` is its own path: usage on stdout, exit 0, never the
/// same code as an argument error (exit 2).
#[derive(Debug, PartialEq)]
enum ParseResult {
    Ok(Command2),
    Help,
    Err(String),
}

fn lookup_engine(key: &str) -> Option<&'static str> {
    ENGINES
        .iter()
        .find(|(k, _, _)| *k == key)
        .map(|(_, exe, _)| *exe)
}

fn available_keys() -> String {
    ENGINES
        .iter()
        .map(|(k, _, _)| *k)
        .collect::<Vec<_>>()
        .join("|")
}

/// Full grouped usage text, also used by `--help`.
fn usage() -> String {
    let mut out = String::new();
    out.push_str("usage: chuanflow version\n");
    out.push_str("       chuanflow run <engine> <case.dict>\n\n");
    out.push_str("<case.dict> is the case dictionary path, passed verbatim to the engine.\n");
    out.push_str("Binary lookup: $CHUANFLOW_BIN_DIR/<exe> first, otherwise a PATH lookup.\n\n");
    out.push_str("Engines by domain:\n");
    let mut last_group: &str = "";
    for (key, exe, group) in ENGINES {
        if *group != last_group {
            out.push_str(&format!("\n  {group}:\n"));
            last_group = group;
        }
        out.push_str(&format!("    {:26} -> {exe}\n", key));
    }
    out
}

fn parse(args: &[String]) -> ParseResult {
    match args.first().map(String::as_str) {
        None => ParseResult::Err("missing command".to_string()),
        Some("help") | Some("--help") | Some("-h") => ParseResult::Help,
        Some("version") => ParseResult::Ok(Command2::Version),
        Some("run") => {
            let engine = match args.get(1).map(String::as_str) {
                Some(key) => match lookup_engine(key) {
                    Some(exe) => exe,
                    None => {
                        return ParseResult::Err(format!(
                            "unknown engine '{}': expected one of {}",
                            key,
                            available_keys()
                        ))
                    }
                },
                None => return ParseResult::Err("missing <engine>".to_string()),
            };
            let case = match args.get(2) {
                Some(case) => case,
                None => return ParseResult::Err("missing <case.dict>".to_string()),
            };
            ParseResult::Ok(Command2::Run { engine: engine.to_string(), case: case.clone() })
        }
        other => ParseResult::Err(format!("unknown command '{}'", other.unwrap_or("<empty>"))),
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
        ParseResult::Help => {
            print!("{}", usage());
            process::exit(0);
        }
        ParseResult::Ok(cmd) => cmd,
        ParseResult::Err(msg) => {
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
                let where_ = match env::var("CHUANFLOW_BIN_DIR") {
                    Ok(dir) => format!("$CHUANFLOW_BIN_DIR={dir} ({}/{engine})", dir),
                    Err(_) => "PATH lookup".to_string(),
                };
                let _ = writeln!(io::stderr(), "chuanflow: cannot run {engine}: searched {where_}; {err}");
                127
            }
        },
    };
    process::exit(code);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Expected registry, pinned here so any CMake target rename shows up as a
    /// test failure instead of a silent 127 at runtime.
    fn expected_table() -> Vec<(&'static str, &'static str, &'static str)> {
        vec![
            ("steady", "chflow-1d-steady", "1D 稳态"),
            ("transient", "chflow-1d-transient", "1D 瞬态"),
            ("nodal", "chflow-1d-nodal", "nodal 节点分析"),
            ("twofluid", "chflow-1d-twofluid", "两相"),
            ("pvt", "chflow-pvt", "PVT 流体相态"),
            ("3d-scalar-transport", "chflow-scalar-transport", "3D"),
            ("3d-simple", "chflow-simple-foam", "3D"),
            ("3d-piso", "chflow-piso-foam", "3D"),
            ("3d-vof", "chflow-vof", "3D"),
            ("3d-compressible", "chflow-3d-compressible", "3D"),
            ("3d-dpm", "chflow-3d-dpm", "3D"),
            ("3d-p1", "chflow-3d-p1", "3D"),
            ("3d-combustion", "chflow-3d-combustion", "3D"),
            ("3d-combustion-detailed", "chflow-3d-combustion-detailed", "3D"),
            ("3d-qmom", "chflow-3d-qmom", "3D"),
            ("3d-mixture-pbm", "chflow-3d-mixture-pbm", "3D"),
            ("3d-cht", "chflow-3d-cht", "3D"),
            ("3d-ale", "chflow-3d-ale", "3D"),
            ("3d-les", "chflow-3d-les", "3D"),
            ("3d-rsm", "chflow-3d-rsm", "3D"),
            ("3d-melting", "chflow-3d-melting", "3D"),
            ("3d-cavitation", "chflow-3d-cavitation", "3D"),
            ("3d-radiation-do", "chflow-3d-radiation-do", "3D"),
            ("3d-reconstruction", "chflow-reconstruction", "3D"),
            ("3d-parallel", "chflow-3d-parallel", "3D"),
            ("3d-dynamic-mesh", "chflow-3d-dynamic-mesh", "3D"),
            ("3d-adaptive-refine", "chflow-3d-adaptive-refine", "3D"),
            ("3d-mesh", "chflow-openfoam-mesh", "3D"),
            ("3d-backend", "chflow-3d-backend", "3D"),
            ("3d-mpi-linear-solver", "chflow-mpi-linear-solver", "3D (MPI 条件编译)"),
            ("3d-petsc-distributed", "chflow-3d-petsc-distributed", "3D (MPI+PETSc 条件编译)"),
            ("coupled", "chflow-coupled-1d3d", "耦合"),
            ("coupled-transient", "chflow-coupled-1d3d-transient", "耦合"),
            ("coupled-multi", "chflow-coupled-1d3d-multi", "耦合"),
            ("coupled-transient-multi", "chflow-coupled-1d3d-transient-multi", "耦合"),
        ]
    }

    #[test]
    fn version_parses() {
        assert_eq!(parse(&["version".to_string()]), ParseResult::Ok(Command2::Version));
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
            ParseResult::Ok(Command2::Run { engine: "chflow-1d-steady".into(), case: "case.dict".into() })
        );
    }

    #[test]
    fn run_transient_parses() {
        let got = parse(&[
            "run".to_string(),
            "transient".to_string(),
            "case.dict".to_string(),
        ]);
        assert_eq!(
            got,
            ParseResult::Ok(Command2::Run { engine: "chflow-1d-transient".into(), case: "case.dict".into() })
        );
    }

    #[test]
    fn run_pvt_parses() {
        let got = parse(&[
            "run".to_string(),
            "pvt".to_string(),
            "case.dict".to_string(),
        ]);
        assert_eq!(
            got,
            ParseResult::Ok(Command2::Run { engine: "chflow-pvt".into(), case: "case.dict".into() })
        );
    }

    #[test]
    fn table_is_drift_free() {
        let expected = expected_table();
        assert_eq!(ENGINES.len(), expected.len());
        for (got, exp) in ENGINES.iter().zip(expected.iter()) {
            assert_eq!(got, exp);
        }
    }

    #[test]
    fn table_keys_and_exes_are_unique() {
        let mut keys: Vec<&str> = ENGINES.iter().map(|(k, _, _)| *k).collect();
        let mut exes: Vec<&str> = ENGINES.iter().map(|(_, e, _)| *e).collect();
        keys.sort();
        exes.sort();
        assert!(keys.windows(2).all(|w| w[0] != w[1]), "duplicate CLI key");
        assert!(exes.windows(2).all(|w| w[0] != w[1]), "duplicate binary");
    }

    #[test]
    fn every_registered_engine_parses_to_its_exe() {
        for (key, exe, _) in ENGINES {
            let got = parse(&[
                "run".to_string(),
                key.to_string(),
                "case.dict".to_string(),
            ]);
            assert_eq!(
                got,
                ParseResult::Ok(Command2::Run { engine: exe.to_string(), case: "case.dict".to_string() }),
                "key {key} did not map to {exe}"
            );
        }
    }

    #[test]
    fn unknown_engine_is_rejected_with_key_hint() {
        let ParseResult::Err(err) =
            parse(&["run".to_string(), "nonsense".to_string(), "x".to_string()])
        else {
            panic!("unknown engine must be an error outcome");
        };
        assert!(err.contains("unknown engine 'nonsense'"));
        assert!(err.contains("steady"), "hint should list existing keys");
        assert!(err.contains("3d-vof"), "hint should list 3D keys");
        assert!(err.contains("coupled-transient-multi"), "hint should list coupled keys");
    }

    #[test]
    fn rejects_bad_usage() {
        assert!(matches!(parse(&[]), ParseResult::Err(_)));
        assert!(matches!(parse(&["run".to_string()]), ParseResult::Err(_)));
        assert!(matches!(
            parse(&["run".to_string(), "mesh".to_string()]),
            ParseResult::Err(_)
        ));
        assert!(matches!(
            parse(&["run".to_string(), "steady".to_string()]),
            ParseResult::Err(_)
        ));
        assert!(matches!(
            parse(&["transient".to_string(), "case".to_string()]),
            ParseResult::Err(_)
        ));
    }

    /// Help is its own outcome with exit 0; unknown engine and missing argument
    /// stay exit 2. The three codes must not collapse onto one return value.
    #[test]
    fn help_unknown_engine_missing_arg_have_distinct_codes() {
        let help = parse(&["--help".to_string()]);
        let help_short = parse(&["-h".to_string()]);
        let unknown = parse(&["run".to_string(), "bogus".to_string(), "x.dict".to_string()]);
        let missing_arg = parse(&["run".to_string(), "steady".to_string()]);

        assert_eq!(help, ParseResult::Help);
        assert_eq!(help_short, ParseResult::Help);
        assert!(matches!(unknown, ParseResult::Err(_)));
        assert!(matches!(missing_arg, ParseResult::Err(_)));

        let code = |r: &ParseResult| match r {
            ParseResult::Help => 0,
            ParseResult::Ok(_) => 0,
            ParseResult::Err(_) => 2,
        };
        assert_ne!(code(&help), code(&unknown));
        assert_ne!(code(&help), code(&missing_arg));
        assert_eq!(code(&unknown), code(&missing_arg));
    }

    #[test]
    fn help_lists_usage_on_stdout_path() {
        // Help must carry the usage text, not an error message.
        assert_eq!(parse(&["help".to_string()]), ParseResult::Help);
        assert!(usage().starts_with("usage:"));
    }

    #[test]
    fn usage_lists_every_group_and_key() {
        let text = usage();
        for group in ["1D 稳态", "1D 瞬态", "nodal 节点分析", "两相", "PVT 流体相态", "3D", "耦合"] {
            assert!(text.contains(group), "usage missing group {group}");
        }
        for (key, exe, _) in ENGINES {
            assert!(text.contains(key), "usage missing key {key}");
            assert!(text.contains(exe), "usage missing exe {exe}");
        }
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

    #[test]
    fn cli_help_exits_zero() {
        let Some(exe) = option_env!("CARGO_BIN_EXE_chuanflow") else { return };
        for arg in ["--help", "-h"] {
            let out = Command::new(exe)
                .arg(arg)
                .output()
                .expect("chuanflow binary");
            assert_eq!(out.status.code(), Some(0), "{arg} must exit 0");
            assert!(
                String::from_utf8_lossy(&out.stdout).contains("usage:"),
                "{arg} usage text on stdout"
            );
        }
    }

    #[test]
    fn cli_unknown_engine_and_missing_arg_exit_two() {
        let Some(exe) = option_env!("CARGO_BIN_EXE_chuanflow") else { return };
        let unknown = Command::new(exe)
            .args(["run", "bogus", "x.dict"])
            .output()
            .expect("chuanflow binary");
        assert_eq!(unknown.status.code(), Some(2), "unknown engine must exit 2");
        let missing = Command::new(exe)
            .args(["run", "steady"])
            .output()
            .expect("chuanflow binary");
        assert_eq!(missing.status.code(), Some(2), "missing <case.dict> must exit 2");
    }
}
