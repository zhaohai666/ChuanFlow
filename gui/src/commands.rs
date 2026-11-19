//! ChuanFlow M5 command layer.
//!
//! Everything here is plain Rust: the parsing helpers are pure functions and
//! are covered by `#[test]` without ever starting a WebView. The two
//! `#[tauri::command]` wrappers are thin enough that the tests exercise the
//! exact code path the GUI calls.

use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use serde::{Deserialize, Serialize};

/// Monotonic id for temp files. Rust tests run in parallel threads within one
/// process, so `std::process::id()` alone would collide and clobber another
/// test's case script; this keeps every concurrent call's path unique.
static NEXT_TMP_ID: AtomicU64 = AtomicU64::new(0);
fn unique_tmp_id() -> u64 {
    NEXT_TMP_ID.fetch_add(1, Ordering::Relaxed)
}

// ---------------------------------------------------------------------------
// Case model + parser
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct NodeDef {
    /// Node name, e.g. `A`.
    pub name: String,
    /// `Some(p)` for a fixed-pressure node, `None` for a free node.
    pub pressure: Option<f64>,
    /// Transient-only fields. Absent in steady cases; kept optional so a
    /// steady parse -> serialize round-trip stays byte-symmetric.
    pub volume: Option<f64>,
    pub bulk_modulus: Option<f64>,
    pub liquid_fraction: Option<f64>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PipeDef {
    pub name: String,
    pub from: String,
    pub to: String,
    /// Pipe geometry, kept so drag-edited networks serialize into a case the
    /// C++ engine can actually run. Defaults match `network.cpp::loadNetwork`.
    pub length: f64,
    pub diameter: f64,
    pub roughness: f64,
    pub elevation: f64,
}

impl Default for PipeDef {
    fn default() -> Self {
        Self {
            name: String::new(),
            from: String::new(),
            to: String::new(),
            length: 1000.0,
            diameter: 0.05,
            roughness: 0.0,
            elevation: 0.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SourceDef {
    /// Dictionary entry name (the label before `{ ... }`).
    pub name: String,
    /// Node the source is attached to (the `node` key; may differ from name).
    pub node: String,
    /// Injection rate (m^3/s); negative means withdrawal.
    pub rate: f64,
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct NetworkCase {
    pub nodes: Vec<NodeDef>,
    pub pipes: Vec<PipeDef>,
    pub sources: Vec<SourceDef>,
    /// Optional `time{}` block; present for transient cases, absent for steady.
    pub time: Option<TimeDef>,
}

/// Time controls kept so a transient case can be edited and serialized.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TimeDef {
    pub dt: f64,
    pub end_time: f64,
    pub report_every: i64,
}

impl Default for TimeDef {
    fn default() -> Self {
        Self { dt: 0.01, end_time: 1.0, report_every: 1 }
    }
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
        "time" => "time",
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
        // Single-line `time { dt ...; endTime ...; }` — time has no named
        // sub-entries, so parse its whole body as key/value pairs.
        if tokens.first().copied() == Some("time") && line.contains('{') && line.contains('}') {
            let open = line.find('{').unwrap();
            let close = line.rfind('}').unwrap();
            update_time(&mut case.time, &parse_kv_body(&line[open + 1..close]));
            continue;
        }

        let kind = section.clone();
        // Multi-line `time { ... }`: each line is a bare `key value;` pair.
        if kind.as_deref() == Some("time") {
            let kv = parse_kv_body(line);
            if !kv.is_empty() {
                update_time(&mut case.time, &kv);
            }
            continue;
        }
        if matches!(kind.as_deref(), Some("nodes" | "pipes" | "sources")) {
            let entry = parse_entry(line)
                .ok_or_else(|| format!("line {line_no}: cannot parse entry: {line}"))?;
            match kind.as_deref().unwrap() {
                "nodes" => {
                    let get = |k: &str| entry.kv.get(k).and_then(|v| v.parse::<f64>().ok());
                    case.nodes.push(NodeDef {
                        name: entry.name,
                        pressure: get("pressure"),
                        volume: get("volume"),
                        bulk_modulus: get("bulkModulus"),
                        liquid_fraction: get("liquidFraction"),
                    });
                }
                "pipes" => {
                    let from = entry.kv.get("from").cloned().ok_or_else(|| {
                        format!("line {line_no}: pipe {} missing 'from'", entry.name)
                    })?;
                    let to = entry.kv.get("to").cloned().ok_or_else(|| {
                        format!("line {line_no}: pipe {} missing 'to'", entry.name)
                    })?;
                    let d = PipeDef::default();
                    let num = |k: &str, fallback: f64| {
                        entry.kv.get(k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(fallback)
                    };
                    case.pipes.push(PipeDef {
                        name: entry.name,
                        from,
                        to,
                        length: num("length", d.length),
                        diameter: num("diameter", d.diameter),
                        roughness: num("roughness", d.roughness),
                        elevation: num("elevation", d.elevation),
                    });
                }
                "sources" => {
                    let rate = entry
                        .kv
                        .get("rate")
                        .and_then(|v| v.parse::<f64>().ok())
                        .unwrap_or(0.0);
                    // Real cases: `well1 { node W1; rate ...; }`; the GUI's
                    // simplified form omits `node` and uses the entry name.
                    let node = entry.kv.get("node").cloned().unwrap_or_else(|| entry.name.clone());
                    case.sources.push(SourceDef { name: entry.name, node, rate });
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

/// Parse a bare `key value; key value;` body (no entry name) into a map. Used
/// for `time{}` which has no named sub-entries, unlike nodes/pipes/sources.
fn parse_kv_body(body: &str) -> HashMap<String, String> {
    let mut kv = HashMap::new();
    for stmt in body.split(';') {
        let stmt = stmt.trim();
        if stmt.is_empty() {
            continue;
        }
        let mut it = stmt.splitn(2, char::is_whitespace);
        let key = it.next().unwrap_or("").trim();
        let value = it.next().unwrap_or("").trim();
        if !key.is_empty() {
            kv.insert(key.to_string(), value.to_string());
        }
    }
    kv
}

/// Fold one `time{}` key/value set into the optional time block. Used by both
/// the single-line (`time { dt ...; }`) and multi-line forms.
fn update_time(time: &mut Option<TimeDef>, kv: &HashMap<String, String>) {
    let dt = kv.get("dt").and_then(|v| v.parse::<f64>().ok());
    let end = kv.get("endTime").and_then(|v| v.parse::<f64>().ok());
    if let (Some(dt), Some(end)) = (dt, end) {
        *time = Some(TimeDef {
            dt,
            end_time: end,
            report_every: kv
                .get("reportEvery")
                .and_then(|v| v.parse::<f64>().ok())
                .unwrap_or(1.0) as i64,
        });
    }
}

/// Format an `f64` the way the C++ `num()` helper does (`setprecision(12)`),
/// so a parse -> serialize -> parse round-trip keeps every value bit-identical.
#[allow(dead_code)] // used by serialize_network_case (round-trip authority + tests)
fn fmt_num(v: f64) -> String {
    let s = format!("{:.12}", v);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() || trimmed == "-0" { "0".to_string() } else { trimmed.to_string() }
}

/// Serialize a `NetworkCase` back into a dictionary the C++ engine accepts.
///
/// This is the inverse of `parse_case_text`: every field the parser keeps is
/// written back (pipe geometry, source `node=` key, transient `time{}`), and
/// `parse_case_text(serialize_network_case(&c)) == c` holds for every case the
/// parser accepts. It is a pure function so round-trip correctness is cheap to
/// unit-test and the drag editor can call it before running the solver.
#[allow(dead_code)] // authority for the JS editor serializer; exercised by round-trip tests
pub fn serialize_network_case(case: &NetworkCase) -> String {
    let mut out = String::new();
    out.push_str("// generated by ChuanFlow GUI drag editor\n");
    out.push_str("solver { type steadyNetwork; }\n");
    out.push_str("fluid { rho 850; mu 0.1; }\n\n");

    out.push_str("nodes\n{\n");
    for n in &case.nodes {
        let mut fields = Vec::new();
        if let Some(p) = n.pressure {
            fields.push(format!("pressure {}", fmt_num(p)));
        }
        if let Some(v) = n.volume {
            fields.push(format!("volume {}", fmt_num(v)));
        }
        if let Some(b) = n.bulk_modulus {
            fields.push(format!("bulkModulus {}", fmt_num(b)));
        }
        if let Some(a) = n.liquid_fraction {
            fields.push(format!("liquidFraction {}", fmt_num(a)));
        }
        if fields.is_empty() {
            out.push_str(&format!("    {} {{}}\n", n.name));
        } else {
            out.push_str(&format!("    {} {{ {} }}\n", n.name, fields.join("; ") + ";"));
        }
    }
    out.push_str("}\n\npipes\n{\n");
    for p in &case.pipes {
        out.push_str(&format!(
            "    {} {{ from {}; to {}; length {}; diameter {}; roughness {}; elevation {}; }}\n",
            p.name, p.from, p.to,
            fmt_num(p.length), fmt_num(p.diameter), fmt_num(p.roughness), fmt_num(p.elevation),
        ));
    }
    out.push_str("}\n\nsources\n{\n");
    for s in &case.sources {
        // Always emit the explicit `node` key: real C++ loader reads it, and
        // the parser round-trips it regardless of whether it equals the name.
        out.push_str(&format!(
            "    {} {{ node {}; rate {}; }}\n",
            s.name, s.node, fmt_num(s.rate),
        ));
    }
    out.push_str("}\n");

    if let Some(t) = &case.time {
        out.push_str(&format!(
            "\ntime {{ dt {}; endTime {}; reportEvery {}; }}\n",
            fmt_num(t.dt), fmt_num(t.end_time), t.report_every,
        ));
    }
    out
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

// ---------------------------------------------------------------------------
// Transient output model + parser
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TransientSample {
    pub t: f64,
    /// node name -> pressure (Pa)
    pub pressure: HashMap<String, f64>,
    /// node name -> liquid fraction alpha in [0,1]
    pub alpha: HashMap<String, f64>,
    /// node name -> hydrate fraction in [0,1]
    pub hydrate: HashMap<String, f64>,
    /// pipe name -> volumetric flow rate Q (m^3/s)
    pub flow: HashMap<String, f64>,
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct TransientSeries {
    pub node_names: Vec<String>,
    pub pipe_names: Vec<String>,
    pub dt: f64,
    pub end_time: f64,
    pub steps: u64,
    pub samples: Vec<TransientSample>,
}

/// Parse the stdout of `chflow-1d-transient` into a time-ordered series.
///
/// The C++ `transientApp.cpp` emits one `t <value>` line followed by one
/// `node <name> p <p> alpha <a> hydrate <h>` per node and one
/// `pipe <name> Q <q>` per pipe, ending with `transient solved: steps=N`.
/// Tolerates the `[chflow] ` prefix and unknown lines; this is a pure function
/// so it is unit-testable without starting the engine.
pub fn parse_transient_output(out: &str) -> Result<TransientSeries, String> {
    let mut series = TransientSeries::default();
    let mut current: Option<TransientSample> = None;

    for raw in out.lines() {
        let line = raw.trim().trim_start_matches("[chflow]").trim();
        if line.is_empty() {
            continue;
        }

        if let Some(rest) = line.strip_prefix("transient:") {
            for tok in rest.split_whitespace() {
                if let Some(v) = tok.strip_prefix("dt=") {
                    series.dt = v.parse().unwrap_or(0.0);
                } else if let Some(v) = tok.strip_prefix("endTime=") {
                    series.end_time = v.parse().unwrap_or(0.0);
                }
            }
            continue;
        }

        if let Some(rest) = line.strip_prefix("t ") {
            if let Some(prev) = current.take() {
                series.samples.push(prev);
            }
            let t = rest
                .split_whitespace()
                .next()
                .and_then(|v| v.parse::<f64>().ok())
                .ok_or_else(|| format!("cannot parse time sample: {line}"))?;
            current = Some(TransientSample {
                t,
                pressure: HashMap::new(),
                alpha: HashMap::new(),
                hydrate: HashMap::new(),
                flow: HashMap::new(),
            });
            continue;
        }

        if let Some(rest) = line.strip_prefix("node ") {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            // node <name> p <p> alpha <a> hydrate <h>  =  7 tokens
            if parts.len() >= 7 && parts[1] == "p" && parts[3] == "alpha" && parts[5] == "hydrate" {
                if let (Ok(p), Ok(a), Ok(h)) = (
                    parts[2].parse::<f64>(),
                    parts[4].parse::<f64>(),
                    parts[6].parse::<f64>(),
                ) {
                    if let Some(sample) = current.as_mut() {
                        let name = parts[0].to_string();
                        if !series.node_names.contains(&name) {
                            series.node_names.push(name.clone());
                        }
                        sample.pressure.insert(name.clone(), p);
                        sample.alpha.insert(name.clone(), a);
                        sample.hydrate.insert(name, h);
                    }
                }
            }
            continue;
        }

        if let Some(rest) = line.strip_prefix("pipe ") {
            let parts: Vec<&str> = rest.split_whitespace().collect();
            if parts.len() >= 3 && parts[1] == "Q" {
                if let Ok(q) = parts[2].parse::<f64>() {
                    if let Some(sample) = current.as_mut() {
                        let name = parts[0].to_string();
                        if !series.pipe_names.contains(&name) {
                            series.pipe_names.push(name.clone());
                        }
                        sample.flow.insert(name, q);
                    }
                }
            }
            continue;
        }

        if let Some(rest) = line.strip_prefix("transient solved:") {
            if let Some(s) = rest.split_whitespace().find(|t| t.starts_with("steps=")) {
                series.steps = s["steps=".len()..].parse().unwrap_or(0);
            }
        }
    }

    if let Some(prev) = current.take() {
        series.samples.push(prev);
    }
    if series.samples.is_empty() {
        return Err("no transient time samples found in solver output".to_string());
    }
    Ok(series)
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

/// Resolve the default engine executables for this installation.
///
/// A packaged .app carries the solvers in `Contents/Resources/bin/`; the
/// lookup walks up from the running binary (so it works wherever the .app
/// lives, regardless of the Finder-launched CWD being `/`) and falls back to
/// the developer tree relative to this source file. Pure path resolution - no
/// processes are spawned - so it is unit-testable.
pub fn resolve_default_engines() -> (String, String) {
    let steady_name = "chflow-1d-steady";
    let transient_name = "chflow-1d-transient";

    // packaged layout: <...>/ChuanFlow GUI.app/Contents/MacOS/chflow-gui
    //              -> <...>/ChuanFlow GUI.app/Contents/Resources/bin/<exe>
    if let Ok(exe) = std::env::current_exe() {
        for ancestor in exe.ancestors().skip(1) {
            // the Contents directory: its parent is the .app bundle
            if ancestor.file_name().map(|n| n == "Contents").unwrap_or(false) {
                let bin = ancestor.join("Resources/bin");
                let steady = bin.join(steady_name);
                let transient = bin.join(transient_name);
                if steady.exists() && transient.exists() {
                    return (
                        steady.to_string_lossy().into_owned(),
                        transient.to_string_lossy().into_owned(),
                    );
                }
                break;
            }
        }
    }

    // developer tree: gui/src/commands.rs -> <root>/build/src/1d/apps
    let dev = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("build/src/1d/apps");
    (
        dev.join(steady_name).to_string_lossy().into_owned(),
        dev.join(transient_name).to_string_lossy().into_owned(),
    )
}

/// Frontend hook: the default engine paths for this installation (packaged
/// Resources/bin first, developer tree second).
#[tauri::command]
pub fn default_engines() -> Vec<String> {
    let (steady, transient) = resolve_default_engines();
    vec![steady, transient]
}

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

    let path = std::env::temp_dir().join(format!(
        "chflow_gui_case_{}_{}.dict", std::process::id(), unique_tmp_id()));
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

/// Run a transient case, writing its text to a temp dictionary and parsing the
/// time series from stdout. The case must carry a `time{}` block.
#[tauri::command]
pub fn run_transient_text(case_text: String, exe_path: String) -> Result<TransientSeries, String> {
    let case = parse_case_text(&case_text)?;
    if case.time.is_none() {
        return Err("transient case needs a time{} block (dt, endTime)".to_string());
    }
    let path = std::env::temp_dir().join(format!(
        "chflow_gui_transient_{}_{}.dict", std::process::id(), unique_tmp_id()));
    std::fs::write(&path, &case_text)
        .map_err(|e| format!("failed to write temporary case file: {e}"))?;
    let output = std::process::Command::new(&exe_path)
        .arg(&path)
        .output()
        .map_err(|e| format!("failed to spawn '{exe_path}': {e}"));
    let _ = std::fs::remove_file(&path);
    let output = output?;
    if !output.status.success() {
        return Err(format!(
            "'{exe_path}' exited with {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr).trim()
        ));
    }
    parse_transient_output(&String::from_utf8_lossy(&output.stdout))
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
    fn resolve_default_engines_returns_existing_paths() {
        let (steady, transient) = resolve_default_engines();
        // At least one of the two layouts must resolve: the packaged
        // Resources/bin (only inside a bundle) or the developer tree
        // (always present after a CMake build in CI/dev).
        let dev_root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().join("build");
        if dev_root.exists() {
            assert!(steady.ends_with("chflow-1d-steady"));
            assert!(transient.ends_with("chflow-1d-transient"));
            assert!(Path::new(&steady).exists(), "dev fallback {steady} must exist");
            assert!(Path::new(&transient).exists(), "dev fallback {transient} must exist");
        }
    }

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

        let script = std::env::temp_dir().join(format!(
            "chflow_gui_fake_{}_{}.sh", std::process::id(), unique_tmp_id()));
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

    const SAMPLE_TRANSIENT_OUT: &str = "\
[chflow] transient: 2 nodes, 1 pipes, dt=0.0001 endTime=0.3
[chflow] t 0
[chflow] node A p 200000 alpha 0.5 hydrate 0
[chflow] node B p 100000 alpha 0.5 hydrate 0
[chflow] pipe CROSS Q 0
[chflow] t 0.0001
[chflow] node A p 199994.5 alpha 0.4999 hydrate 0.01
[chflow] node B p 100005.5 alpha 0.5001 hydrate 0
[chflow] pipe CROSS Q 9.98e-05
[chflow] transient solved: steps=3000
";

    #[test]
    fn parse_transient_output_full() {
        let s = parse_transient_output(SAMPLE_TRANSIENT_OUT).expect("transient parses");
        assert_eq!(s.dt, 0.0001);
        assert_eq!(s.end_time, 0.3);
        assert_eq!(s.steps, 3000);
        assert_eq!(s.samples.len(), 2);
        assert_eq!(s.node_names, vec!["A", "B"]);
        assert_eq!(s.pipe_names, vec!["CROSS"]);
        assert_eq!(s.samples[0].t, 0.0);
        assert_eq!(s.samples[0].pressure["A"], 200_000.0);
        assert_eq!(s.samples[0].alpha["B"], 0.5);
        assert_eq!(s.samples[0].hydrate["A"], 0.0);
        assert_eq!(s.samples[0].flow["CROSS"], 0.0);
        assert_eq!(s.samples[1].t, 0.0001);
        assert_eq!(s.samples[1].pressure["B"], 100_005.5);
        assert_eq!(s.samples[1].hydrate["A"], 0.01);
        assert_eq!(s.samples[1].flow["CROSS"], 9.98e-5);
    }

    #[test]
    fn parse_transient_output_empty_is_error() {
        assert!(parse_transient_output("banner only\nno samples\n").is_err());
    }

    #[cfg(unix)]
    #[test]
    fn run_transient_text_with_fake_exe() {
        use std::fs;
        use std::os::unix::fs::PermissionsExt;

        let script = std::env::temp_dir().join(format!(
            "chflow_gui_fake_transient_{}_{}.sh", std::process::id(), unique_tmp_id()));
        fs::write(&script, format!("#!/bin/sh\ncat <<'EOF'\n{}\nEOF\n", SAMPLE_TRANSIENT_OUT)).unwrap();
        let mut perms = fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).unwrap();

        let case = "nodes\n{\n A { pressure 1; }\n B {}\n}\npipes\n{\n P { from A; to B; }\n}\ntime { dt 0.0001; endTime 0.3; }\n";
        let s = run_transient_text(case.into(), script.to_string_lossy().into_owned())
            .expect("fake transient runs");
        assert_eq!(s.samples.len(), 2);
        assert_eq!(s.steps, 3000);
        let _ = fs::remove_file(&script);
    }

    #[test]
    fn run_transient_text_requires_time_block() {
        let case = "nodes\n{\n A {}\n}\npipes\n{\n P { from A; to A; }\n}\n";
        let err = run_transient_text(case.into(), "/no/exe".into()).unwrap_err();
        assert!(err.contains("time{}"), "unexpected error: {err}");
    }

    #[test]
    fn run_transient_text_with_real_engine_if_available() {
        // Round-trip authority: a model serialized by the GUI must be a case
        // the real C++ transient engine accepts and produces a parseable series.
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let exe = root.join("build/src/1d/apps/chflow-1d-transient");
        if !exe.exists() {
            return;
        }
        let mut case = parse_case_text(
            "nodes\n{\n A { volume 1; bulkModulus 1e8; pressure 2e5; liquidFraction 0.5; }\n\
             B { volume 1; bulkModulus 1e8; pressure 1e5; liquidFraction 0.5; }\n}\n\
             pipes\n{\n CROSS { from A; to B; length 2.49; diameter 0.0519; }\n}\n\
             time { dt 0.0001; endTime 0.01; reportEvery 1; }\n",
        )
        .expect("case parses");
        // Inject a zero source so both the steady and transient loaders agree.
        case.sources.clear();
        let text = serialize_network_case(&case);
        let series = run_transient_text(text, exe.to_string_lossy().into_owned())
            .expect("serialized case runs in real engine");
        assert!(series.samples.len() >= 2, "expected multiple time samples");
        assert_eq!(series.node_names, vec!["A", "B"]);
        assert_eq!(series.pipe_names, vec!["CROSS"]);
        assert_eq!(series.samples[0].t, 0.0);
        assert_eq!(series.samples[0].pressure["A"], 200_000.0);
        assert!(series.steps >= 100);
    }

    #[test]
    fn serialize_network_case_round_trip() {
        let case = parse_case_text(SAMPLE_CASE).expect("sample parses");
        let text = serialize_network_case(&case);
        let back = parse_case_text(&text).expect("serialized text re-parses");
        assert_eq!(back, case, "round-trip must preserve every field");
    }

    #[test]
    fn serialize_network_case_round_trip_with_time_and_named_source() {
        let text = "\
nodes
{
    A { pressure 5e5; }
    B {}
}
pipes
{
    LINE { from A; to B; length 2.49; diameter 0.0519; roughness 0; elevation 0; }
}
sources
{
    well1 { node A; rate 0.002; }
}
time { dt 0.01; endTime 1; reportEvery 10; }
";
        let case = parse_case_text(text).expect("parse");
        assert_eq!(case.sources[0].name, "well1");
        assert_eq!(case.sources[0].node, "A");
        assert!(case.time.is_some());
        let back = parse_case_text(&serialize_network_case(&case)).expect("re-parse");
        assert_eq!(back, case);
        // Serialized text must be runnable: real loader reads `node` and time.
        assert!(serialize_network_case(&case).contains("node A;"));
        assert!(serialize_network_case(&case).contains("endTime 1"));
    }
}
