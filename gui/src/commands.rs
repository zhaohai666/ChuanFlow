//! ChuanFlow M5 command layer.
//!
//! Everything here is plain Rust: the parsing helpers are pure functions and
//! are covered by `#[test]` without ever starting a WebView. The two
//! `#[tauri::command]` wrappers are thin enough that the tests exercise the
//! exact code path the GUI calls.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Case model + parser
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct NodeDef {
    /// Node name, e.g. `A`.
    pub name: String,
    /// `Some(p)` for a fixed-pressure node, `None` for a free node.
    pub pressure: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PipeDef {
    pub name: String,
    pub from: String,
    pub to: String,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SourceDef {
    /// Node the source is attached to.
    pub node: String,
    /// Injection rate (m^3/s); negative means withdrawal.
    pub rate: f64,
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct NetworkCase {
    pub nodes: Vec<NodeDef>,
    pub pipes: Vec<PipeDef>,
    pub sources: Vec<SourceDef>,
}

#[derive(Debug)]
struct Entry {
    name: String,
    kv: HashMap<String, String>,
}

fn strip_comment(line: &str) -> &str {
    let mut end = line.len();
    for (i, _) in line.char_indices() {
        let rest = &line[i..];
        if rest.starts_with("//") || rest.starts_with('#') {
            end = i;
            break;
        }
    }
    &line[..end]
}

fn section_kind(name: &str) -> &'static str {
    match name {
        "nodes" => "nodes",
        "pipes" => "pipes",
        "sources" => "sources",
        // solver/fluid/pvt/... are valid dictionary blocks we simply ignore.
        _ => "__skip__",
    }
}

fn parse_entry(line: &str) -> Option<Entry> {
    let open = line.find('{')?;
    let close = line.rfind('}')?;
    if close < open {
        return None;
    }
    let name = line[..open].trim();
    if name.is_empty() || !name.chars().next()?.is_ascii_alphabetic() {
        return None;
    }
    if !name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
        return None;
    }

    let mut kv = HashMap::new();
    for stmt in line[open + 1..close].split(';') {
        let stmt = stmt.trim();
        if stmt.is_empty() {
            continue;
        }
        let mut it = stmt.splitn(2, char::is_whitespace);
        let key = it.next()?.trim();
        if key.is_empty() {
            continue;
        }
        let value = it.next().unwrap_or("").trim();
        kv.insert(key.to_string(), value.to_string());
    }
    Some(Entry { name: name.to_string(), kv })
}

/// Parse a simplified ChuanFlow dictionary (OpenFOAM style) into a `NetworkCase`.
///
/// Recognised blocks are `nodes`, `pipes` and `sources`; other blocks such as
/// `solver` / `fluid` are accepted but ignored. Entries are expected on a
/// single line, e.g. `A { pressure 3.0e4; }`. Comments (`//` and `#`) are
/// stripped. This is a pure function so it is cheap to unit-test and to call
/// from the GUI for live validation.
pub fn parse_case_text(text: &str) -> Result<NetworkCase, String> {
    let mut case = NetworkCase::default();
    let mut section: Option<String> = None;
    let mut pending: Option<String> = None;

    for (idx, raw) in text.lines().enumerate() {
        let line_no = idx + 1;
        let line = strip_comment(raw).trim();
        if line.is_empty() {
            continue;
        }

        if line == "}" {
            section = None;
            pending = None;
            continue;
        }
        if line == "{" {
            section = Some(pending.take().unwrap_or_else(|| "__skip__".to_string()));
            continue;
        }

        let tokens: Vec<&str> = line.split_whitespace().collect();
        let starts_with_ident = tokens.first().map(|t| is_identifier(t)).unwrap_or(false);

        // A bare identifier opens a block whose `{` sits on the next line.
        if tokens.len() == 1 && starts_with_ident {
            pending = Some(section_kind(tokens[0]).to_string());
            continue;
        }
        // `nodes {` on a single line.
        if tokens.len() == 2 && tokens[1] == "{" && starts_with_ident {
            section = Some(section_kind(tokens[0]).to_string());
            continue;
        }

        let kind = section.clone();
        if matches!(kind.as_deref(), Some("nodes" | "pipes" | "sources")) {
            let entry = parse_entry(line)
                .ok_or_else(|| format!("line {line_no}: cannot parse entry: {line}"))?;
            match kind.as_deref().unwrap() {
                "nodes" => {
                    let pressure = entry
                        .kv
                        .get("pressure")
                        .and_then(|v| v.parse::<f64>().ok());
                    case.nodes.push(NodeDef { name: entry.name, pressure });
                }
                "pipes" => {
                    let from = entry.kv.get("from").cloned().ok_or_else(|| {
                        format!("line {line_no}: pipe {} missing 'from'", entry.name)
                    })?;
                    let to = entry.kv.get("to").cloned().ok_or_else(|| {
                        format!("line {line_no}: pipe {} missing 'to'", entry.name)
                    })?;
                    case.pipes.push(PipeDef { name: entry.name, from, to });
                }
                "sources" => {
                    let rate = entry
                        .kv
                        .get("rate")
                        .and_then(|v| v.parse::<f64>().ok())
                        .unwrap_or(0.0);
                    case.sources.push(SourceDef { node: entry.name, rate });
                }
                _ => {}
            }
            continue;
        }

        // One-line block we do not care about, e.g. `solver { type steadyNetwork; }`.
        if line.contains('{') && line.contains('}') && starts_with_ident {
            continue;
        }
        if kind.as_deref() == Some("__skip__") {
            continue;
        }
        return Err(format!("line {line_no}: content outside of a block: {line}"));
    }

    if case.nodes.is_empty() && case.pipes.is_empty() && case.sources.is_empty() {
        return Err("empty case: no nodes, pipes or sources found".to_string());
    }
    for pipe in &case.pipes {
        if !case.nodes.iter().any(|n| n.name == pipe.from) {
            return Err(format!(
                "pipe {}: unknown from-node {}",
                pipe.name, pipe.from
            ));
        }
        if !case.nodes.iter().any(|n| n.name == pipe.to) {
            return Err(format!("pipe {}: unknown to-node {}", pipe.name, pipe.to));
        }
    }
    Ok(case)
}

fn is_identifier(s: &str) -> bool {
    let mut chars = s.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
}

// ---------------------------------------------------------------------------
// Solver output model + parser
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NodeResult {
    pub name: String,
    pub pressure: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PipeResult {
    pub name: String,
    pub q: f64,
    pub v: f64,
    pub dp: f64,
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct SolverSummary {
    pub nodes: Vec<NodeResult>,
    pub pipes: Vec<PipeResult>,
    pub iterations: u64,
    pub residual: f64,
    pub converged: bool,
}

fn value_after<'a>(s: &'a str, key: &str) -> Option<&'a str> {
    let pos = s.find(key)? + key.len();
    s[pos..].split_whitespace().next()
}

/// Parse the stdout of `chflow-1d-steady` into a structured summary.
///
/// Tolerates the `[chflow] ` log prefix and ignores lines it does not
/// understand, so extra fields on pipe rows (`Re`, `f`, `regime`, ...) are
/// harmless. Values are read as token/value pairs, which is exactly how the
/// C++ logger emits them (steadyNetworkApp.cpp).
pub fn parse_solver_output(out: &str) -> SolverSummary {
    let mut summary = SolverSummary::default();

    for raw in out.lines() {
        let line = raw.trim().trim_start_matches("[chflow]").trim();
        if line.is_empty() {
            continue;
        }

        if let Some(rest) = line.strip_prefix("node ") {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.len() >= 3 && parts[1] == "p" {
                if let Ok(p) = parts[2].parse::<f64>() {
                    summary.nodes.push(NodeResult {
                        name: parts[0].to_string(),
                        pressure: p,
                    });
                }
            }
            continue;
        }

        if let Some(rest) = line.strip_prefix("pipe ") {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }
            let name = parts[0].to_string();
            let mut q = 0.0f64;
            let mut v = 0.0f64;
            let mut dp = 0.0f64;
            // Scan every token rather than assuming a fixed layout: regime
            // flags and optional fields vary from case to case.
            for i in 0..parts.len().saturating_sub(1) {
                match parts[i] {
                    "Q" => q = parts[i + 1].parse().unwrap_or(q),
                    "v" => v = parts[i + 1].parse().unwrap_or(v),
                    "dp" => dp = parts[i + 1].parse().unwrap_or(dp),
                    _ => {}
                }
            }
            summary.pipes.push(PipeResult { name, q, v, dp });
            continue;
        }

        if line.contains("steady solved:") {
            if let Some(it) = value_after(line, "iterations=").and_then(|s| s.parse::<u64>().ok()) {
                summary.iterations = it;
            }
            if let Some(r) = value_after(line, "residual=").and_then(|s| s.parse::<f64>().ok()) {
                summary.residual = r;
            }
            summary.converged = line.contains("(converged)");
        }
    }

    summary
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

/// Spawn an engine executable with a case file and parse its stdout.
#[tauri::command]
pub fn run_case(case_path: String, exe_path: String) -> Result<SolverSummary, String> {
    let output = std::process::Command::new(&exe_path)
        .arg(&case_path)
        .output()
        .map_err(|e| format!("failed to spawn '{exe_path}': {e}"))?;

    if !output.status.success() {
        return Err(format!(
            "'{exe_path}' exited with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    Ok(parse_solver_output(&stdout))
}

/// Convenience command for the editor: validate the text, write it to a
/// temporary dictionary, run the engine, then clean up.
#[tauri::command]
pub fn run_case_text(case_text: String, exe_path: String) -> Result<SolverSummary, String> {
    let _ = parse_case_text(&case_text)?;

    let path = std::env::temp_dir().join(format!("chflow_gui_case_{}.dict", std::process::id()));
    std::fs::write(&path, &case_text)
        .map_err(|e| format!("failed to write temporary case file: {e}"))?;
    let result = run_case(path.to_string_lossy().into_owned(), exe_path);
    let _ = std::fs::remove_file(&path);
    result
}

/// Validate case text from the frontend without running the engine.
#[tauri::command]
pub fn validate_case(case_text: String) -> Result<NetworkCase, String> {
    parse_case_text(&case_text)
}

// ---------------------------------------------------------------------------
// Tests (no WebView required)
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    const SAMPLE_CASE: &str = r#"
        // minimal three-node network
        solver { type steadyNetwork; }
        nodes
        {
            A { pressure 3.0e4; }
            B { pressure 2.0e4; }
            E {};
        }
        pipes
        {
            P1 { from A; to E; length 20; diameter 0.02; }
            P2 { from E; to B; length 20; diameter 0.02; }
        }
        sources
        {
            E { rate 1e-6; }
        }
    "#;

    const SAMPLE_OUT: &str = "\
[chflow] network: 3 nodes, 2 pipes
[chflow] node A p 30000
[chflow] node E p 25000
[chflow] node B p 20000
[chflow] pipe P1 Q 1e-6 v 0.003 Re 50 f 0.5 dp 5000
[chflow] pipe P2 Q -1e-6 v -0.003 Re 50 f 0.5 dp 5000
[chflow] steady solved: iterations=7 residual=1e-9 (converged)
";

    #[test]
    fn parse_case_text_valid() {
        let case = parse_case_text(SAMPLE_CASE).expect("sample case parses");
        assert_eq!(case.nodes.len(), 3);
        assert_eq!(case.nodes[0].name, "A");
        assert_eq!(case.nodes[0].pressure, Some(30_000.0));
        assert_eq!(case.nodes[2].name, "E");
        assert_eq!(case.nodes[2].pressure, None);
        assert_eq!(case.pipes.len(), 2);
        assert_eq!((case.pipes[0].from.as_str(), case.pipes[0].to.as_str()), ("A", "E"));
        assert_eq!(case.sources.len(), 1);
        assert_eq!(case.sources[0].node, "E");
        assert_eq!(case.sources[0].rate, 1e-6);
    }

    #[test]
    fn parse_case_text_empty_is_error() {
        assert!(parse_case_text("").is_err());
        assert!(parse_case_text("// just a comment\n# and another\n").is_err());
    }

    #[test]
    fn parse_case_text_bad_line_is_error() {
        let bad = "nodes\n{\n    A { pressure 3e4; }\n    this is not an entry\n}\n";
        let err = parse_case_text(bad).unwrap_err();
        assert!(err.contains("cannot parse entry"), "unexpected error: {err}");
    }

    #[test]
    fn parse_case_text_rejects_unknown_node() {
        let bad = "nodes\n{\n A { pressure 1; }\n}\npipes\n{\n P { from A; to MISSING; }\n}\n";
        let err = parse_case_text(bad).unwrap_err();
        assert!(err.contains("unknown to-node MISSING"), "unexpected error: {err}");
    }

    #[test]
    fn parse_solver_output_full() {
        let s = parse_solver_output(SAMPLE_OUT);
        assert_eq!(s.nodes.len(), 3);
        assert_eq!(s.nodes[1].pressure, 25_000.0);
        assert_eq!(s.pipes.len(), 2);
        assert_eq!(s.pipes[0].q, 1e-6);
        assert_eq!(s.pipes[0].v, 0.003);
        assert_eq!(s.pipes[0].dp, 5000.0);
        assert_eq!(s.pipes[1].q, -1e-6);
        assert_eq!(s.iterations, 7);
        assert_eq!(s.residual, 1e-9);
        assert!(s.converged);
    }

    #[test]
    fn parse_solver_output_not_converged_and_noise() {
        let out = "\
random banner text
[chflow] node X p 1.5
[chflow] pipe R Q 2 v 3 dp 4 Re 10 f 0.2 dpValve 7 CHOKED
[chflow] steady solved: iterations=200 residual=0.03 (NOT converged)
";
        let s = parse_solver_output(out);
        assert_eq!(s.nodes.len(), 1);
        assert_eq!(s.nodes[0].name, "X");
        assert_eq!(s.pipes[0].q, 2.0);
        assert_eq!(s.pipes[0].dp, 4.0);
        assert_eq!(s.iterations, 200);
        assert!(!s.converged);
    }

    #[cfg(unix)]
    #[test]
    fn run_case_with_fake_exe() {
        use std::fs;
        use std::os::unix::fs::PermissionsExt;

        let script = std::env::temp_dir().join(format!("chflow_gui_fake_{}.sh", std::process::id()));
        fs::write(
            &script,
            "#!/bin/sh\necho '[chflow] node A p 30000'\n\
             echo '[chflow] pipe P1 Q 1.5 v 2 Re 100 f 0.03 dp 400'\n\
             echo '[chflow] steady solved: iterations=7 residual=1e-9 (converged)'\n",
        )
        .unwrap();
        let mut perms = fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).unwrap();

        let exe = script.to_string_lossy().into_owned();
        let s = run_case("ignored.dict".into(), exe).expect("fake exe runs");
        assert_eq!(s.nodes[0].pressure, 30_000.0);
        assert_eq!(s.pipes[0].q, 1.5);
        assert_eq!(s.pipes[0].dp, 400.0);
        assert!(s.converged);
        let _ = fs::remove_file(&script);
    }

    #[test]
    fn run_case_missing_exe_is_error() {
        let err = run_case(
            "case.dict".into(),
            "/definitely/not/a/real/chflow-binary".into(),
        )
        .unwrap_err();
        assert!(err.contains("failed to spawn"), "unexpected error: {err}");
    }

    #[test]
    fn run_case_with_real_engine_if_available() {
        // End-to-end smoke test against the actual C++ engine when it has
        // already been built by the top-level CMake project.
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let exe = root.join("build/src/1d/apps/chflow-1d-steady");
        let case = root.join("tests/cases/steadyAllBelowThreshold/case.dict");
        if !exe.exists() || !case.exists() {
            return;
        }
        let s = run_case(case.to_string_lossy().into_owned(), exe.to_string_lossy().into_owned())
            .expect("real engine runs");
        assert_eq!(s.nodes.len(), 3);
        assert_eq!(s.pipes.len(), 5);
        assert!(s.converged);
    }
}
