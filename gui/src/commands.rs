//! ChuanFlow M5 command layer.
//!
//! Everything here is plain Rust: the parsing helpers are pure functions and
//! are covered by `#[test]` without ever starting a WebView. The two
//! `#[tauri::command]` wrappers are thin enough that the tests exercise the
//! exact code path the GUI calls.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

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
    /// Inline equipment the steadyNetwork kernel reads off the same pipe
    /// (`network.cpp::loadNetwork`): `valve` is the valve admittance K
    /// (0 = no valve) and the pump family is the affinity-law curve. 0 keys
    /// stay absent from the serialization so a plain pipe round-trips
    /// byte-identically.
    #[serde(default)]
    pub valve: f64,
    #[serde(default)]
    pub pump_head: f64,
    #[serde(default)]
    pub pump_rate: f64,
    #[serde(default)]
    pub pump_shutoff_head: f64,
    #[serde(default)]
    pub pump_stages: f64,
    #[serde(default)]
    pub pump_speed: f64,
    /// ESP advanced intake semantics (`network.cpp` reads both keys; the
    /// affinity-law pump is steady-network, transient rejects pumps).
    #[serde(default)]
    pub pump_gas_fraction: f64,
    #[serde(default)]
    pub pump_intake_pressure: f64,
    /// Gas-lift valve (§5.89): the dome pressure is the switch, the rest are
    /// its port/spring parameters. Steady-network element.
    #[serde(default)]
    pub lift_dome_pressure: f64,
    #[serde(default)]
    pub lift_area_ratio: f64,
    #[serde(default)]
    pub lift_spring_pressure: f64,
    #[serde(default)]
    pub lift_discharge_coeff: f64,
    #[serde(default)]
    pub lift_kappa: f64,
    /// Rod pump (§5.90): the barrel bore is the switch. Steady-network
    /// element (transient explicitly refuses it: no polished rod exists).
    #[serde(default)]
    pub rod_plunger_diameter: f64,
    #[serde(default)]
    pub rod_stroke_length: f64,
    #[serde(default)]
    pub rod_strokes_per_minute: f64,
    #[serde(default)]
    pub rod_diameter: f64,
    #[serde(default)]
    pub rod_length: f64,
    #[serde(default)]
    pub rod_young_modulus: f64,
    #[serde(default)]
    pub rod_density: f64,
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
            valve: 0.0,
            pump_head: 0.0,
            pump_rate: 0.0,
            pump_shutoff_head: 0.0,
            pump_stages: 0.0,
            pump_speed: 0.0,
            pump_gas_fraction: 0.0,
            pump_intake_pressure: 0.0,
            lift_dome_pressure: 0.0,
            lift_area_ratio: 0.0,
            lift_spring_pressure: 0.0,
            lift_discharge_coeff: 0.0,
            lift_kappa: 0.0,
            rod_plunger_diameter: 0.0,
            rod_stroke_length: 0.0,
            rod_strokes_per_minute: 0.0,
            rod_diameter: 0.0,
            rod_length: 0.0,
            rod_young_modulus: 0.0,
            rod_density: 0.0,
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

/// The `fluid{}` block the steadyNetwork engine reads. Defaults match
/// `network.cpp::loadNetwork` so a model without an explicit block stays
/// byte-compatible with the M5 serializer.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FluidDef {
    pub rho: f64,
    pub mu: f64,
}

impl Default for FluidDef {
    fn default() -> Self {
        Self { rho: 850.0, mu: 0.1 }
    }
}

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct NetworkCase {
    pub nodes: Vec<NodeDef>,
    pub pipes: Vec<PipeDef>,
    pub sources: Vec<SourceDef>,
    #[serde(default)]
    pub fluid: FluidDef,
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
        "fluid" => "fluid",
        // solver/pvt/... are valid dictionary blocks we simply ignore.
        _ => "__skip__",
    }
}

/// Fold one `fluid{}` key/value set into the case's fluid block. Keys other
/// than rho/mu are ignored (the steadyNetwork loader reads exactly these).
fn update_fluid(fluid: &mut FluidDef, kv: &HashMap<String, String>) {
    if let Some(v) = kv.get("rho").and_then(|v| v.parse::<f64>().ok()) {
        fluid.rho = v;
    }
    if let Some(v) = kv.get("mu").and_then(|v| v.parse::<f64>().ok()) {
        fluid.mu = v;
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

    'lines: for (idx, raw) in text.lines().enumerate() {
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
        // Single-line `time { dt ...; endTime ...; }` and the same one-line
        // form of `fluid { rho ...; mu ...; }` — both have no named
        // sub-entries, so parse their whole body as key/value pairs.
        for flat in ["time", "fluid"] {
            if tokens.first().copied() == Some(flat)
                && line.contains('{')
                && line.contains('}')
            {
                let open = line.find('{').unwrap();
                let close = line.rfind('}').unwrap();
                let kv = parse_kv_body(&line[open + 1..close]);
                if flat == "time" {
                    update_time(&mut case.time, &kv);
                } else {
                    update_fluid(&mut case.fluid, &kv);
                }
                continue 'lines;
            }
        }

        let kind = section.clone();
        // Multi-line `time{}`/`fluid{}`: each line is a bare `key value;` pair.
        if matches!(kind.as_deref(), Some("time" | "fluid")) {
            let kv = parse_kv_body(line);
            if !kv.is_empty() {
                if kind.as_deref() == Some("time") {
                    update_time(&mut case.time, &kv);
                } else {
                    update_fluid(&mut case.fluid, &kv);
                }
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
                        valve: num("valve", 0.0),
                        pump_head: num("pumpHead", 0.0),
                        pump_rate: num("pumpRate", 0.0),
                        pump_shutoff_head: num("pumpShutoffHead", 0.0),
                        pump_stages: num("pumpStages", 0.0),
                        pump_speed: num("pumpSpeed", 0.0),
                        pump_gas_fraction: num("pumpGasFraction", 0.0),
                        pump_intake_pressure: num("pumpIntakePressure", 0.0),
                        lift_dome_pressure: num("liftDomePressure", 0.0),
                        lift_area_ratio: num("liftAreaRatio", 0.0),
                        lift_spring_pressure: num("liftSpringPressure", 0.0),
                        lift_discharge_coeff: num("liftDischargeCoeff", 0.0),
                        lift_kappa: num("liftKappa", 0.0),
                        rod_plunger_diameter: num("rodPlungerDiameter", 0.0),
                        rod_stroke_length: num("rodStrokeLength", 0.0),
                        rod_strokes_per_minute: num("rodStrokesPerMinute", 0.0),
                        rod_diameter: num("rodDiameter", 0.0),
                        rod_length: num("rodLength", 0.0),
                        rod_young_modulus: num("rodYoungModulus", 0.0),
                        rod_density: num("rodDensity", 0.0),
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
    out.push_str(&format!(
        "fluid {{ rho {}; mu {}; }}\n\n",
        fmt_num(case.fluid.rho),
        fmt_num(case.fluid.mu)
    ));

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
        // inline equipment keys are emitted only when present (0 = absent),
        // so a plain pipe round-trips byte-identically through the parser
        let mut fields = format!(
            "from {}; to {}; length {}; diameter {}; roughness {}; elevation {};",
            p.from, p.to,
            fmt_num(p.length), fmt_num(p.diameter), fmt_num(p.roughness), fmt_num(p.elevation)
        );
        if p.valve != 0.0 {
            fields.push_str(&format!(" valve {};", fmt_num(p.valve)));
        }
        if p.pump_head != 0.0 {
            fields.push_str(&format!(" pumpHead {};", fmt_num(p.pump_head)));
        }
        if p.pump_rate != 0.0 {
            fields.push_str(&format!(" pumpRate {};", fmt_num(p.pump_rate)));
        }
        if p.pump_shutoff_head != 0.0 {
            fields.push_str(&format!(" pumpShutoffHead {};", fmt_num(p.pump_shutoff_head)));
        }
        if p.pump_stages != 0.0 {
            fields.push_str(&format!(" pumpStages {};", fmt_num(p.pump_stages)));
        }
        if p.pump_speed != 0.0 {
            fields.push_str(&format!(" pumpSpeed {};", fmt_num(p.pump_speed)));
        }
        if p.pump_gas_fraction != 0.0 {
            fields.push_str(&format!(" pumpGasFraction {};", fmt_num(p.pump_gas_fraction)));
        }
        if p.pump_intake_pressure != 0.0 {
            fields.push_str(&format!(" pumpIntakePressure {};", fmt_num(p.pump_intake_pressure)));
        }
        if p.lift_dome_pressure != 0.0 {
            fields.push_str(&format!(" liftDomePressure {};", fmt_num(p.lift_dome_pressure)));
        }
        if p.lift_area_ratio != 0.0 {
            fields.push_str(&format!(" liftAreaRatio {};", fmt_num(p.lift_area_ratio)));
        }
        if p.lift_spring_pressure != 0.0 {
            fields.push_str(&format!(" liftSpringPressure {};", fmt_num(p.lift_spring_pressure)));
        }
        if p.lift_discharge_coeff != 0.0 {
            fields.push_str(&format!(" liftDischargeCoeff {};", fmt_num(p.lift_discharge_coeff)));
        }
        if p.lift_kappa != 0.0 {
            fields.push_str(&format!(" liftKappa {};", fmt_num(p.lift_kappa)));
        }
        if p.rod_plunger_diameter != 0.0 {
            fields.push_str(&format!(" rodPlungerDiameter {};", fmt_num(p.rod_plunger_diameter)));
        }
        if p.rod_stroke_length != 0.0 {
            fields.push_str(&format!(" rodStrokeLength {};", fmt_num(p.rod_stroke_length)));
        }
        if p.rod_strokes_per_minute != 0.0 {
            fields.push_str(&format!(" rodStrokesPerMinute {};", fmt_num(p.rod_strokes_per_minute)));
        }
        if p.rod_diameter != 0.0 {
            fields.push_str(&format!(" rodDiameter {};", fmt_num(p.rod_diameter)));
        }
        if p.rod_length != 0.0 {
            fields.push_str(&format!(" rodLength {};", fmt_num(p.rod_length)));
        }
        if p.rod_young_modulus != 0.0 {
            fields.push_str(&format!(" rodYoungModulus {};", fmt_num(p.rod_young_modulus)));
        }
        if p.rod_density != 0.0 {
            fields.push_str(&format!(" rodDensity {};", fmt_num(p.rod_density)));
        }
        out.push_str(&format!("    {} {{ {} }}\n", p.name, fields));
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
// Solver job registry: start / poll / stop
// ---------------------------------------------------------------------------

/// One running (or finished) engine process. Output is accumulated by reader
/// threads so the UI can poll a tail without blocking the solver thread.
pub struct Job {
    /// Kept for the log / future registry introspection.
    #[allow(dead_code)]
    pub label: String,
    pub child: Arc<Mutex<std::process::Child>>,
    pub output: Arc<Mutex<String>>,
    pub finished: Arc<AtomicBool>,
    pub exit: Arc<Mutex<Option<i32>>>,
}

/// Tauri-managed registry of engine jobs keyed by a monotonic id.
#[derive(Default)]
pub struct JobRegistry {
    next: AtomicU64,
    jobs: Mutex<HashMap<u64, Job>>,
}

/// Keep the tail bounded: GUI log only ever shows the most recent output.
const OUTPUT_CAP: usize = 512 * 1024;

fn append_output(buf: &Mutex<String>, chunk: &str) {
    let mut s = buf.lock().unwrap();
    s.push_str(chunk);
    if s.len() > OUTPUT_CAP {
        let cut = s.len() - OUTPUT_CAP;
        // cut at a line boundary so a tail never starts mid-token
        let drop_to = s[cut..].find('\n').map(|p| cut + p + 1).unwrap_or(s.len());
        s.drain(..drop_to);
    }
}

/// Spawn an engine job. The inner logic is Tauri-free so it is unit-testable.
pub fn spawn_job(
    registry: &JobRegistry,
    label: &str,
    exe_path: &str,
    case_path: &str,
) -> Result<u64, String> {
    use std::process::{Command, Stdio};

    let mut child = Command::new(exe_path)
        .arg(case_path)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("failed to spawn '{exe_path}': {e}"))?;

    let out_buf: Arc<Mutex<String>> = Arc::new(Mutex::new(String::new()));
    let finished = Arc::new(AtomicBool::new(false));
    let exit: Arc<Mutex<Option<i32>>> = Arc::new(Mutex::new(None));

    let stdout = child
        .stdout
        .take()
        .ok_or_else(|| "child stdout unavailable".to_string())?;
    let stderr = child
        .stderr
        .take()
        .ok_or_else(|| "child stderr unavailable".to_string())?;

    let sink = out_buf.clone();
    std::thread::spawn(move || {
        use std::io::Read;
        let mut r = stdout;
        let mut chunk = [0u8; 4096];
        loop {
            match r.read(&mut chunk) {
                Ok(0) | Err(_) => break,
                Ok(n) => {
                    let text = String::from_utf8_lossy(&chunk[..n]);
                    append_output(&sink, &text);
                }
            }
        }
    });
    let sink = out_buf.clone();
    std::thread::spawn(move || {
        use std::io::Read;
        let mut r = stderr;
        let mut chunk = [0u8; 4096];
        loop {
            match r.read(&mut chunk) {
                Ok(0) | Err(_) => break,
                Ok(n) => {
                    let text = String::from_utf8_lossy(&chunk[..n]);
                    append_output(&sink, &text);
                }
            }
        }
    });

    let child_arc = Arc::new(Mutex::new(child));

    let finished2 = finished.clone();
    let exit2 = exit.clone();
    let wait_child = child_arc.clone();
    std::thread::spawn(move || {
        // wait() needs &mut Child: lock, wait, release before publishing the
        // exit code so stop_job can still lock during the wait.
        let code = {
            let mut c = wait_child.lock().unwrap();
            c.wait().map(|st| st.code().unwrap_or(-1)).unwrap_or(-1)
        };
        *exit2.lock().unwrap() = Some(code);
        finished2.store(true, Ordering::SeqCst);
    });

    let id = registry.next.fetch_add(1, Ordering::SeqCst);
    registry.jobs.lock().unwrap().insert(
        id,
        Job {
            label: label.to_string(),
            child: child_arc,
            output: out_buf,
            finished,
            exit,
        },
    );
    Ok(id)
}

/// One poll result: whether the job finished, its exit code (once finished)
/// and the recent output tail.
#[derive(Serialize)]
pub struct JobStatus {
    pub finished: bool,
    pub exit: Option<i32>,
    pub output: String,
}

const OUTPUT_TAIL: usize = 8000;

fn poll_job_inner(registry: &JobRegistry, id: u64) -> Option<JobStatus> {
    let jobs = registry.jobs.lock().unwrap();
    let job = jobs.get(&id)?;
    let finished = job.finished.load(Ordering::SeqCst);
    let output = {
        let s = job.output.lock().unwrap();
        let start = s.len().saturating_sub(OUTPUT_TAIL);
        // when the cap truncated, land the tail on a line boundary; when the
        // whole output fits, keep it from the very beginning
        let start = if start > 0 {
            s[start..].find('\n').map(|p| start + p + 1).unwrap_or(start)
        } else {
            0
        };
        s[start..].to_string()
    };
    let exit = if finished {
        *job.exit.lock().unwrap()
    } else {
        None
    };
    Some(JobStatus { finished, exit, output })
}

fn stop_job_inner(registry: &JobRegistry, id: u64) -> Result<bool, String> {
    let jobs = registry.jobs.lock().unwrap();
    match jobs.get(&id) {
        Some(job) => {
            // Locking child's mutex: the waiter thread only waits, never
            // locks; kill is safe against the reader threads.
            let mut c = job.child.lock().unwrap();
            match c.kill() {
                Ok(()) => Ok(true),
                // already exited
                Err(_) => Ok(false),
            }
        }
        None => Err(format!("unknown job id {id}")),
    }
}

fn drop_job_inner(registry: &JobRegistry, id: u64) {
    registry.jobs.lock().unwrap().remove(&id);
}

/// Start an engine executable as a background job. Returns the job id the
/// frontend polls with `poll_job` and can terminate with `stop_job`.
#[tauri::command]
pub fn start_case(
    label: String,
    exe_path: String,
    case_path: String,
    state: tauri::State<JobRegistry>,
) -> Result<u64, String> {
    spawn_job(&state, &label, &exe_path, &case_path)
}

/// Poll a job: recent output tail, finished flag and exit code (when done).
/// A finished job is removed from the registry by this call.
#[tauri::command]
pub fn poll_job(id: u64, state: tauri::State<JobRegistry>) -> Result<JobStatus, String> {
    let status = poll_job_inner(&state, id).ok_or_else(|| format!("unknown job id {id}"))?;
    if status.finished {
        drop_job_inner(&state, id);
    }
    Ok(status)
}

/// Terminate a running job (the stop button). Safe to call repeatedly.
#[tauri::command]
pub fn stop_job(id: u64, state: tauri::State<JobRegistry>) -> Result<bool, String> {
    let stopped = stop_job_inner(&state, id)?;
    if stopped {
        drop_job_inner(&state, id);
    }
    Ok(stopped)
}

// ---------------------------------------------------------------------------
// Project / artifact file IO
// ---------------------------------------------------------------------------

/// Read a UTF-8 text file (project JSON, .vtu field output, case dict).
#[tauri::command]
pub fn read_text_file(path: String) -> Result<String, String> {
    std::fs::read_to_string(&path).map_err(|e| format!("failed to read '{path}': {e}"))
}

/// Write a UTF-8 text file (project JSON, exported case dict).
#[tauri::command]
pub fn write_text_file(path: String, content: String) -> Result<(), String> {
    if let Some(parent) = Path::new(&path).parent() {
        if !parent.as_os_str().is_empty() {
            let _ = std::fs::create_dir_all(parent);
        }
    }
    std::fs::write(&path, content).map_err(|e| format!("failed to write '{path}': {e}"))
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

/// Resolve the default engine executables for this installation.
///
/// A packaged .app carries the solvers in `Contents/Resources/bin/` and the
/// coupling template case in `Contents/Resources/cases/`; the lookup walks up
/// from the running binary (so it works wherever the .app lives, regardless
/// of the Finder-launched CWD being `/`) and falls back to the developer tree
/// relative to this source file. Pure path resolution - no processes are
/// spawned - so it is unit-testable.
pub fn resolve_default_engines() -> (String, String, String, String, String, String, String) {
    let steady_name = "chflow-1d-steady";
    let transient_name = "chflow-1d-transient";
    let coupled_name = "chflow-coupled-1d3d";
    let twofluid_name = "chflow-1d-twofluid";
    let nodal_name = "chflow-1d-nodal";
    let pvt_name = "chflow-pvt";
    let template_rel = "tests/cases/coupled1d3d/case.dict";

    // packaged layout: <...>/ChuanFlow GUI.app/Contents/MacOS/chflow-gui
    //              -> <...>/ChuanFlow GUI.app/Contents/Resources/bin/<exe>
    //              -> <...>/ChuanFlow GUI.app/Contents/Resources/cases/<case>
    let packaged: Option<(PathBuf, PathBuf, PathBuf, PathBuf, PathBuf, PathBuf, PathBuf)> =
        std::env::current_exe().ok().and_then(|exe| {
            for ancestor in exe.ancestors().skip(1) {
                // the Contents directory: its parent is the .app bundle
                if ancestor.file_name().map(|n| n == "Contents").unwrap_or(false) {
                    let bin = ancestor.join("Resources/bin");
                    let cases = ancestor.join("Resources/cases");
                    let steady = bin.join(steady_name);
                    let transient = bin.join(transient_name);
                    let coupled = bin.join(coupled_name);
                    let twofluid = bin.join(twofluid_name);
                    let nodal = bin.join(nodal_name);
                    let pvt = bin.join(pvt_name);
                    if steady.exists() && transient.exists() {
                        let template = if cases.join("coupled1d3d.dict").exists() {
                            cases.join("coupled1d3d.dict")
                        } else {
                            coupled.clone()
                        };
                        return Some((steady, transient, coupled, template, twofluid, nodal, pvt));
                    }
                    break;
                }
            }
            None
        });

    if let Some((steady, transient, coupled, template, twofluid, nodal, pvt)) = packaged {
        return (
            steady.to_string_lossy().into_owned(),
            transient.to_string_lossy().into_owned(),
            coupled.to_string_lossy().into_owned(),
            template.to_string_lossy().into_owned(),
            twofluid.to_string_lossy().into_owned(),
            nodal.to_string_lossy().into_owned(),
            pvt.to_string_lossy().into_owned(),
        );
    }

    // developer tree: gui/src/commands.rs -> <root>/build/src/...
    let dev = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("build");
    (
        dev.join("src/1d/apps").join(steady_name).to_string_lossy().into_owned(),
        dev.join("src/1d/apps").join(transient_name).to_string_lossy().into_owned(),
        dev.join("src/coupling/apps").join(coupled_name).to_string_lossy().into_owned(),
        dev.parent().unwrap().join(template_rel).to_string_lossy().into_owned(),
        dev.join("src/1d/apps").join(twofluid_name).to_string_lossy().into_owned(),
        dev.join("src/1d/apps").join(nodal_name).to_string_lossy().into_owned(),
        dev.join("src/pvt/apps").join(pvt_name).to_string_lossy().into_owned(),
    )
}

/// E1 (§5.181): resolve a named 3D engine binary (packaged Resources/bin
/// first, developer tree second). The 3D family registry below keys every
/// member by its CMake target name.
pub fn resolve_three_d_engine_named(bin_name: &str) -> String {
    let packaged = std::env::current_exe().ok().and_then(|exe| {
        exe.ancestors()
            .find(|a| a.file_name().map(|n| n == "Contents").unwrap_or(false))
            .map(|c| c.join("Resources/bin").join(bin_name))
            .filter(|p| p.exists())
    });
    match packaged {
        Some(p) => p.to_string_lossy().into_owned(),
        None => {
            let dev = Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .join("build/src/3d/apps")
                .join(bin_name);
            dev.to_string_lossy().into_owned()
        }
    }
}

/// E1: resolve the default 3D solver. piso-foam is the most general
/// pressure-based 3D app and the startedChannel template (a verified PISO
/// case) is verified against it.
pub fn resolve_default_3d_engine() -> String {
    resolve_three_d_engine_named("chflow-piso-foam")
}

/// E1 (§5.181): the 3D engine family registry. One row per 3D app that has a
/// real, engine-verified case (tests/cases/<case_id>/case.dict) AND honours
/// the `output { vtk ...; }` redirect the GUI's 3D flow injects, so the field
/// view can auto-load the VTU after the run. Members without a VTU-producing
/// path (cht, ale, adaptiveRefine, dynamicMesh, mixturePbm, qmom,
/// combustionDetailed) stay out of the GUI registry on purpose - wiring them
/// would promise a field view the app cannot deliver.
/// Rows: (frontend engine key, binary name, case id, human label).
pub const THREE_D_FAMILY: &[(&str, &str, &str, &str)] = &[
    ("threeD", "chflow-piso-foam", "startedChannel", "PISO 瞬态通道启动"),
    ("threeDSimple", "chflow-simple-foam", "poiseuilleChannel", "SIMPLE 稳态泊肃叶通道"),
    ("threeDVof", "chflow-vof", "vofDroplet", "VOF 液滴界面捕捉"),
    ("threeDCompressible", "chflow-3d-compressible", "compressibleNsSod", "可压缩 NS（Sod 激波管）"),
    ("threeDScalar", "chflow-scalar-transport", "advectionDiffusion", "标量对流-扩散"),
    ("threeDDpm", "chflow-3d-dpm", "dpmChannel", "DPM 颗粒通道（单向耦合）"),
    ("threeDP1", "chflow-3d-p1", "radiationP1Slab", "P1 辐射灰色介质板"),
    ("threeDCombustion", "chflow-3d-combustion", "combustionFixedOxidiser", "一步化学燃烧"),
    ("threeDLes", "chflow-3d-les", "lesSgs", "LES 亚格子模型"),
    ("threeDMelting", "chflow-3d-melting", "meltingFront", "熔化前沿（相变）"),
    ("threeDCavitation", "chflow-3d-cavitation", "cavitationVenturi", "Venturi 空化"),
];

/// One 3D family entry handed to the frontend: the persistent engine key,
/// the resolved binary path and the case template that drives it.
#[derive(Debug, Clone, Serialize)]
pub struct ThreeDEngineEntry {
    pub key: String,
    pub label: String,
    pub exe: String,
    pub case_id: String,
    pub case_path: String,
}

/// Resolve a CTest case dictionary: packaged Resources/cases/<id>.dict wins
/// when present, developer tree tests/cases/<id>/case.dict otherwise.
fn resolve_case_path(id: &str) -> PathBuf {
    let dev = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("tests/cases")
        .join(id)
        .join("case.dict");
    let packaged = std::env::current_exe().ok().and_then(|exe| {
        exe.ancestors()
            .find(|a| a.file_name().map(|n| n == "Contents").unwrap_or(false))
            .map(|contents| contents.join("Resources/cases").join(format!("{id}.dict")))
            .filter(|p| p.exists())
    });
    packaged.unwrap_or(dev)
}

/// E1 (§5.181): the full 3D family with resolved binary and case paths.
pub fn resolve_three_d_family() -> Vec<ThreeDEngineEntry> {
    THREE_D_FAMILY
        .iter()
        .map(|(key, bin, case_id, label)| ThreeDEngineEntry {
            key: (*key).to_string(),
            label: (*label).to_string(),
            exe: resolve_three_d_engine_named(bin),
            case_id: (*case_id).to_string(),
            case_path: resolve_case_path(case_id).to_string_lossy().into_owned(),
        })
        .collect()
}

/// Frontend hook: the 3D engine family (key / label / exe / case path).
#[tauri::command]
pub fn three_d_engines() -> Vec<ThreeDEngineEntry> {
    resolve_three_d_family()
}

/// Frontend hook: the default engine paths for this installation
/// (packaged Resources/bin first, developer tree second):
/// [steady, transient, coupled-1d3d, coupled template case, two-fluid, nodal, pvt, 3d piso].
#[tauri::command]
pub fn default_engines() -> Vec<String> {
    let (steady, transient, coupled, template, twofluid, nodal, pvt) = resolve_default_engines();
    vec![steady, transient, coupled, template, twofluid, nodal, pvt, resolve_default_3d_engine()]
}

/// One built-in case template: a kernel-verified case dictionary shipped with
/// the app (Resources/cases/<id>.dict) or the developer tree
/// (tests/cases/<id>/case.dict). Loading one loads a case the engine has
/// already been verified against byte-for-byte.
#[derive(Serialize)]
pub struct TemplateEntry {
    pub id: String,
    pub name: String,
    pub path: String,
}

/// The template library. Every id names a case that is part of the CTest
/// regression suite, so a template is by construction a case the kernel
/// accepts and produces checked physics for.
pub fn resolve_case_templates() -> Vec<TemplateEntry> {
    const TEMPLATES: &[(&str, &str)] = &[
        ("steadyNetwork", "稳态集输管网（基础）"),
        ("steadyEquipmentValve", "稳态节流阀（ laminar 阀律）"),
        ("steadyEquipmentPump", "稳态离心泵（亲和律闭环）"),
        ("steadyNetworkMultiphase", "稳态黑油多相集输（Beggs-Brill）"),
        ("transientBlowdown", "瞬态放空（罐排空 2x2 解析）"),
        ("twoFluidTerrainSlug", "地形段塞（捕集-行进-卸载弛豫振荡）"),
        ("coupled1d3d", "1D↔3D 耦合（网络嵌 SIMPLE 块）"),
        ("nodalVogel", "节点分析（Vogel IPR × Beggs-Brill 垂直举升）"),
        ("pvtFlash", "PVT 闪蒸（烃组分 + 黏度/表面张力）"),
        ("startedChannel", "3D 通道启动（PISO，体积力驱动 + VTU 场）"),
        // E1 (§5.181): 3D family templates - each pairs a registered 3D
        // engine with its engine-verified case (see THREE_D_FAMILY).
        ("poiseuilleChannel", "3D 泊肃叶通道（SIMPLE 稳态解析剖面）"),
        ("vofDroplet", "3D VOF 液滴（表面张力界面捕捉）"),
        ("compressibleNsSod", "3D 可压缩 NS（Sod 激波管）"),
        ("advectionDiffusion", "3D 标量输运（对流-扩散）"),
        ("dpmChannel", "3D DPM 颗粒通道（单向耦合）"),
        ("radiationP1Slab", "3D P1 辐射（灰色介质平板）"),
        ("combustionFixedOxidiser", "3D 燃烧（固定氧化剂一步化学）"),
        ("lesSgs", "3D LES 亚格子输运"),
        ("meltingFront", "3D 熔化前沿（相变）"),
        ("cavitationVenturi", "3D Venturi 空化"),
    ];
    let mut out = Vec::new();
    for (id, name) in TEMPLATES {
        // packaged path wins when it exists (the .app carries the same dicts)
        let path = std::env::current_exe().ok().and_then(|exe| {
            exe.ancestors()
                .find(|a| a.file_name().map(|n| n == "Contents").unwrap_or(false))
                .map(|contents| {
                    contents
                        .join("Resources/cases")
                        .join(format!("{id}.dict"))
                })
        });
        let chosen = match path {
            Some(p) if p.exists() => p,
            _ => resolve_case_path(id),
        };
        out.push(TemplateEntry {
            id: (*id).to_string(),
            name: (*name).to_string(),
            path: chosen.to_string_lossy().into_owned(),
        });
    }
    out
}

/// E5: the Beggs-Brill four regime boundaries, ported 1:1 from the kernel's
/// single source of truth `beggsBrillPattern` (src/1d/src/beggsBrill.cpp,
/// D4 / §5.85). The frontend's regimemap.js carries a hardcoded copy of these
/// closed forms; this command lets it reconcile against the same formulas the
/// Rust layer owns, and a unit test reconciles *this* against the C++ engine's
/// printed L1..L4 report so a kernel-side formula change cannot drift silently.
pub(crate) fn beggs_brill_boundaries_impl(lambda: f64) -> [f64; 4] {
    [
        316.0 * lambda.powf(0.302),
        0.0009252 * lambda.powf(-2.4684),
        0.1 * lambda.powf(-1.4516),
        0.5 * lambda.powf(-6.738),
    ]
}

#[tauri::command]
pub fn beggs_brill_boundaries(lambda: f64) -> [f64; 4] {
    beggs_brill_boundaries_impl(lambda)
}

/// Frontend hook: the built-in template library entries.
#[tauri::command]
pub fn case_templates() -> Vec<TemplateEntry> {
    resolve_case_templates()
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

/// One row of the parallel sweep: per-index result so order is preserved
/// regardless of which worker finished when.
#[derive(Debug, Clone, Serialize)]
pub struct RunRow {
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub summary: Option<SolverSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// Run many case texts concurrently through the same engine binary — the
/// parameter sweep's batch primitive. Each row keeps its input index (the
/// frontend maps values to rows positionally), workers are scoped threads
/// capped by both `available_parallelism` and 8 (process storms are worse
/// than queueing), and every temporary dict is uniquely named so workers
/// never share a file.
#[tauri::command]
pub fn run_case_texts(texts: Vec<String>, exe_path: String) -> Vec<RunRow> {
    if texts.is_empty() {
        return Vec::new();
    }
    let workers = std::cmp::min(
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
        texts.len(),
    )
    .min(8);
    let rows: std::sync::Mutex<Vec<Option<RunRow>>> =
        std::sync::Mutex::new((0..texts.len()).map(|_| None).collect());
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..workers {
            let next = &next;
            let rows = &rows;
            let texts = &texts;
            let exe_path = &exe_path;
            scope.spawn(move || loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if i >= texts.len() {
                    break;
                }
                let row = match run_case_text(texts[i].clone(), exe_path.clone()) {
                    Ok(summary) => RunRow {
                        ok: true,
                        summary: Some(summary),
                        error: None,
                    },
                    Err(e) => RunRow {
                        ok: false,
                        summary: None,
                        error: Some(e),
                    },
                };
                rows.lock().expect("sweep row mutex")[i] = Some(row);
            });
        }
    });
    rows.into_inner()
        .expect("sweep rows mutex")
        .into_iter()
        .enumerate()
        .map(|(i, r)| {
            r.unwrap_or_else(|| RunRow {
                ok: false,
                summary: None,
                error: Some(format!("sweep row {i} never finished")),
            })
        })
        .collect()
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

/// One row of the transient what-if batch. `series` is positionally aligned
/// with the frontend's parameter values, just like `RunRow` for steady sweeps.
#[derive(Debug, Clone, Serialize)]
pub struct TransientRunRow {
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub series: Option<TransientSeries>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// Run many transient case texts through the same engine binary concurrently.
/// The scoped worker pool mirrors `run_case_texts` (capped at 8 workers, unique
/// temp dictionary per row); the row index is preserved regardless of which
/// worker finished first. This is the batch primitive for transient what-if
/// sweeps: one an offline case text per parameter value, all real engine runs.
#[tauri::command]
pub fn run_transient_texts(texts: Vec<String>, exe_path: String) -> Vec<TransientRunRow> {
    if texts.is_empty() {
        return Vec::new();
    }
    let workers = std::cmp::min(
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
        texts.len(),
    )
    .min(8);
    let rows: std::sync::Mutex<Vec<Option<TransientRunRow>>> =
        std::sync::Mutex::new((0..texts.len()).map(|_| None).collect());
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..workers {
            let next = &next;
            let rows = &rows;
            let texts = &texts;
            let exe_path = &exe_path;
            scope.spawn(move || loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if i >= texts.len() {
                    break;
                }
                let row = match run_transient_text(texts[i].clone(), exe_path.clone()) {
                    Ok(series) => TransientRunRow {
                        ok: true,
                        series: Some(series),
                        error: None,
                    },
                    Err(e) => TransientRunRow {
                        ok: false,
                        series: None,
                        error: Some(e),
                    },
                };
                rows.lock().expect("transient sweep row mutex")[i] = Some(row);
            });
        }
    });
    rows.into_inner()
        .expect("transient sweep rows mutex")
        .into_iter()
        .enumerate()
        .map(|(i, r)| {
            r.unwrap_or_else(|| TransientRunRow {
                ok: false,
                series: None,
                error: Some(format!("transient sweep row {i} never finished")),
            })
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Nodal analysis (U2): IPR-VLP working point + batched VLP sampling
// ---------------------------------------------------------------------------

/// One `chflow-1d-nodal` result: the IPR declaration echo, the static column,
/// the working point (rate / bottom-hole pressure) and the solve status.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct NodalSummary {
    pub ipr_form: String,
    pub p_res: f64,
    pub q_max: f64,
    pub static_p: f64,
    pub rate: f64,
    pub pwf: f64,
    pub holdup: f64,
    pub regime: String,
    pub iterations: u64,
    pub flows: bool,
    pub converged: bool,
}

/// Parse the stdout of `chflow-1d-nodal`. The app emits (12 significant
/// digits, [chflow] prefix tolerated):
///   `nodal: ipr <form> p_res <p> qmax <q>, well ...`
///   `column p <p> (the rate-0 static column)`
///   `well q <q> pwf <p> pth <p>`
///   `well holdup <hl> regime <name>`
///   `nodal solved: iterations=N flows=true|false (converged|NOT converged)`
/// Pure function, unit-testable without the engine.
pub fn parse_nodal_output(out: &str) -> Result<NodalSummary, String> {
    let mut s = NodalSummary {
        ipr_form: String::new(),
        p_res: f64::NAN,
        q_max: f64::NAN,
        static_p: f64::NAN,
        rate: f64::NAN,
        pwf: f64::NAN,
        holdup: f64::NAN,
        regime: String::new(),
        iterations: 0,
        flows: false,
        converged: false,
    };
    let mut saw_well = false;
    let mut saw_solved = false;
    for raw in out.lines() {
        let line = raw.trim().replace("[chflow] ", "");
        let line = line.trim();
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() >= 7 && toks[0] == "nodal:" && toks[1] == "ipr" {
            s.ipr_form = toks[2].to_string();
            if toks[3] == "p_res" {
                s.p_res = toks[4].parse().unwrap_or(f64::NAN);
            }
            if toks[5] == "qmax" {
                s.q_max = toks[6].trim_end_matches(',').parse().unwrap_or(f64::NAN);
            }
        } else if toks.len() >= 2 && toks[0] == "column" && toks[1] == "p" {
            s.static_p = toks[2].parse().unwrap_or(f64::NAN);
        } else if toks.len() >= 4 && toks[0] == "well" && toks[1] == "q" {
            s.rate = toks[2].parse().unwrap_or(f64::NAN);
            if toks[3] == "pwf" {
                s.pwf = toks.get(4).and_then(|v| v.parse().ok()).unwrap_or(f64::NAN);
            }
            saw_well = true;
        } else if toks.len() >= 4 && toks[0] == "well" && toks[1] == "holdup" {
            s.holdup = toks[2].parse().unwrap_or(f64::NAN);
            if toks[3] == "regime" {
                s.regime = toks.get(4).copied().unwrap_or("").to_string();
            }
        } else if line.starts_with("nodal solved:") {
            s.iterations = after_key_u64(line, "iterations=");
            s.flows = line.contains("flows=true");
            s.converged = line.contains("(converged)") && !line.contains("(NOT converged)");
            saw_solved = true;
        }
    }
    if !saw_well || !saw_solved {
        return Err("not chflow-1d-nodal output (missing well/nodal solved line)".to_string());
    }
    Ok(s)
}

fn after_key_u64(line: &str, key: &str) -> u64 {
    let pos = match line.find(key) {
        Some(p) => p + key.len(),
        None => return 0,
    };
    line[pos..]
        .split_whitespace()
        .next()
        .unwrap_or("0")
        .parse()
        .unwrap_or(0)
}

/// Run one nodal case text through `chflow-1d-nodal` and parse the working
/// point. Unlike `run_case_text` there is no network validation: a nodal
/// dictionary is not a network case.
pub fn run_nodal_text(case_text: String, exe_path: String) -> Result<NodalSummary, String> {
    let path = std::env::temp_dir().join(format!(
        "chflow_gui_nodal_{}_{}.dict", std::process::id(), unique_tmp_id()));
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
    parse_nodal_output(&String::from_utf8_lossy(&output.stdout))
}

/// One row of the nodal batch (VLP sampling): positionally aligned with the
/// input texts, like `RunRow` / `TransientRunRow`.
#[derive(Debug, Clone, Serialize)]
pub struct NodalRunRow {
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub summary: Option<NodalSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// Run many nodal case texts concurrently through the same engine binary —
/// the VLP-curve sampler's batch primitive (one linear-IPR variant per
/// productivity index; every working point lies on the wellbore's VLP curve).
/// Same scoped worker pool as `run_case_texts` / `run_transient_texts`.
#[tauri::command]
pub fn run_nodal_texts(texts: Vec<String>, exe_path: String) -> Vec<NodalRunRow> {
    if texts.is_empty() {
        return Vec::new();
    }
    let workers = std::cmp::min(
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
        texts.len(),
    )
    .min(8);
    let rows: std::sync::Mutex<Vec<Option<NodalRunRow>>> =
        std::sync::Mutex::new((0..texts.len()).map(|_| None).collect());
    let next = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..workers {
            let next = &next;
            let rows = &rows;
            let texts = &texts;
            let exe_path = &exe_path;
            scope.spawn(move || loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                if i >= texts.len() {
                    break;
                }
                let row = match run_nodal_text(texts[i].clone(), exe_path.clone()) {
                    Ok(summary) => NodalRunRow { ok: true, summary: Some(summary), error: None },
                    Err(e) => NodalRunRow { ok: false, summary: None, error: Some(e) },
                };
                rows.lock().expect("nodal row mutex")[i] = Some(row);
            });
        }
    });
    rows.into_inner()
        .expect("nodal rows mutex")
        .into_iter()
        .enumerate()
        .map(|(i, r)| {
            r.unwrap_or_else(|| NodalRunRow {
                ok: false,
                summary: None,
                error: Some(format!("nodal row {i} never finished")),
            })
        })
        .collect()
}

// ---------------------------------------------------------------------------
// PVT flash report (U4): run chflow-pvt and hand the raw report to the
// frontend parser (the flash case dictionary is not a network case, so there
// is no Rust-side validation; parsing is JS-side like the two-fluid frames).
// ---------------------------------------------------------------------------

/// Write a flash case text to a temporary dictionary, run `chflow-pvt` on it
/// and return its raw stdout (the JS `parsePvtOutput` turns it into curves).
#[tauri::command]
pub fn run_pvt_text(case_text: String, exe_path: String) -> Result<String, String> {
    let path = std::env::temp_dir().join(format!(
        "chflow_gui_pvt_{}_{}.dict", std::process::id(), unique_tmp_id()));
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
    Ok(String::from_utf8_lossy(&output.stdout).into_owned())
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
        let (steady, transient, coupled, template, twofluid, nodal, pvt) = resolve_default_engines();
        // At least one of the two layouts must resolve: the packaged
        // Resources/bin (only inside a bundle) or the developer tree
        // (always present after a CMake build in CI/dev).
        let dev_root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().join("build");
        if dev_root.exists() {
            assert!(steady.ends_with("chflow-1d-steady"));
            assert!(transient.ends_with("chflow-1d-transient"));
            assert!(coupled.ends_with("chflow-coupled-1d3d"));
            assert!(Path::new(&steady).exists(), "dev fallback {steady} must exist");
            assert!(Path::new(&transient).exists(), "dev fallback {transient} must exist");
            assert!(Path::new(&coupled).exists(), "dev fallback {coupled} must exist");
            assert!(Path::new(&twofluid).exists(), "dev fallback {twofluid} must exist");
            assert!(Path::new(&nodal).exists(), "dev fallback {nodal} must exist");
            assert!(Path::new(&pvt).exists(), "dev fallback {pvt} must exist");
            assert!(
                Path::new(&template).exists() && template.ends_with("coupled1d3d/case.dict"),
                "coupling template {template} must exist"
            );
        }
    }

    #[test]
    fn parse_nodal_output_reads_working_point() {
        let out = "[chflow] nodal: ipr vogel p_res 2.2e7 qmax 0.006, well 2400 m x 0.062 m, wellhead 2.5e6 Pa\n\
                   [chflow] column p 18642604.4027 (the rate-0 static column)\n\
                   [chflow] well q 0.00290785579167 pwf 15120468.4282 pth 2.5e6\n\
                   [chflow] well holdup 0.401234567891 regime distributed\n\
                   [chflow] nodal solved: iterations=9 flows=true (converged)\n";
        let s = parse_nodal_output(out).unwrap();
        assert_eq!(s.ipr_form, "vogel");
        assert_eq!(s.p_res, 2.2e7);
        assert_eq!(s.q_max, 0.006);
        assert_eq!(s.static_p, 18642604.4027);
        assert_eq!(s.rate, 0.00290785579167);
        assert_eq!(s.pwf, 15120468.4282);
        assert_eq!(s.holdup, 0.401234567891);
        assert_eq!(s.regime, "distributed");
        assert_eq!(s.iterations, 9);
        assert!(s.flows);
        assert!(s.converged);
    }

    #[test]
    fn parse_nodal_output_no_flow_and_rejects_noise() {
        let s = parse_nodal_output(
            "nodal: ipr linear p_res 2e7 qmax 0.002, well 1000 m x 0.1 m, wellhead 9e6 Pa\n\
             column p 10500000 (the rate-0 static column)\n\
             well q 0 pwf 10500000 pth 9e6\n\
             nodal solved: iterations=12 flows=false (converged)\n",
        )
        .unwrap();
        assert_eq!(s.rate, 0.0);
        assert!(!s.flows);
        assert!(s.converged);
        assert!(parse_nodal_output("steady solved: iterations=3 residual=1e-9 (converged)").is_err());
    }

    #[test]
    fn run_nodal_texts_keeps_order_and_reports_errors() {
        let tmp = std::env::temp_dir().join(format!("chflow_gui_nodal_fake_{}.sh", std::process::id()));
        std::fs::write(
            &tmp,
            "#!/bin/sh\ncat \"$1\" | grep -o 'productivityIndex [0-9.eE+-]*' | awk \
             '{printf \"nodal: ipr linear p_res 2e7 qmax 1, well 1 m x 1 m, wellhead 1 Pa\\n\"; \
             printf \"well q %s pwf 1234 pth 1\\n\", $2*100; \
             printf \"nodal solved: iterations=1 flows=true (converged)\\n\"}'\n",
        )
        .unwrap();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&tmp, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        let mk = |j: &str| format!("ipr {{ form linear; reservoirPressure 2e7; productivityIndex {j}; }}\n");
        let rows = run_nodal_texts(
            vec![mk("0.001"), "not a case".to_string(), mk("0.002")],
            tmp.to_string_lossy().into_owned(),
        );
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(rows.len(), 3);
        assert!(rows[0].ok);
        assert_eq!(rows[0].summary.as_ref().unwrap().rate, 0.1);
        assert!(!rows[1].ok);
        assert!(rows[1].error.as_ref().unwrap().contains("missing well"));
        assert!(rows[2].ok);
        assert_eq!(rows[2].summary.as_ref().unwrap().rate, 0.2);
    }

    #[test]
    fn case_templates_resolve_to_existing_verified_cases() {
        let templates = resolve_case_templates();
        // E1 (§5.181): 10 pre-existing + 10 new 3D family templates.
        assert_eq!(templates.len(), 20);
        for t in &templates {
            assert!(
                Path::new(&t.path).exists(),
                "template {} ({}) must resolve: {}",
                t.id,
                t.name,
                t.path
            );
        }
        let ids: Vec<&str> = templates.iter().map(|t| t.id.as_str()).collect();
        assert!(ids.contains(&"twoFluidTerrainSlug"));
        assert!(ids.contains(&"coupled1d3d"));
        assert!(ids.contains(&"nodalVogel"));
        assert!(ids.contains(&"pvtFlash"));
        assert!(ids.contains(&"startedChannel"));
        // the slug template must end at the verified case file
        let slug = templates.iter().find(|t| t.id == "twoFluidTerrainSlug").unwrap();
        assert!(slug.path.ends_with("twoFluidTerrainSlug/case.dict") || slug.path.ends_with("twoFluidTerrainSlug.dict"));
    }

    /// E1 (§5.181): the 3D engine family registry must resolve every member's
    /// binary (developer tree layout) and its engine-verified case dict.
    #[test]
    fn three_d_family_resolves_binaries_and_verified_cases() {
        let family = resolve_three_d_family();
        assert_eq!(family.len(), THREE_D_FAMILY.len());
        assert_eq!(family.len(), 11);
        // keys are unique and the PISO entry leads the family
        let mut keys: Vec<&str> = family.iter().map(|e| e.key.as_str()).collect();
        assert_eq!(keys[0], "threeD");
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), THREE_D_FAMILY.len(), "family keys must be unique");
        let bins: Vec<&str> = THREE_D_FAMILY.iter().map(|(_, bin, ..)| *bin).collect();
        let mut sorted_bins = bins.to_vec();
        sorted_bins.sort_unstable();
        sorted_bins.dedup();
        assert_eq!(sorted_bins.len(), bins.len(), "family binaries must be unique");
        for e in &family {
            assert!(e.exe.ends_with(THREE_D_FAMILY.iter()
                .find(|(k, ..)| *k == e.key).map(|(_, b, ..)| *b).unwrap()));
            assert!(
                Path::new(&e.case_path).exists(),
                "family case {} must exist: {}",
                e.case_id,
                e.case_path
            );
            assert!(e.case_path.ends_with(&format!("{}/case.dict", e.case_id)));
            assert!(!e.label.is_empty());
        }
        let piso = family.iter().find(|e| e.key == "threeD").unwrap();
        assert!(piso.exe.ends_with("chflow-piso-foam"));
        assert!(piso.case_path.ends_with("startedChannel/case.dict"));
    }

    /// E1 (§5.181): every 3D family case id must be registered as a template
    /// so the template menu and the engine family cannot drift apart.
    #[test]
    fn three_d_family_cases_are_all_registered_templates() {
        let templates = resolve_case_templates();
        let ids: Vec<&str> = templates.iter().map(|t| t.id.as_str()).collect();
        for (key, _bin, case_id, _label) in THREE_D_FAMILY {
            assert!(
                ids.contains(case_id),
                "family engine {key} case {case_id} missing from template registry"
            );
        }
        // the default 3D resolve still points at PISO (frontends persist it)
        assert!(resolve_default_3d_engine().ends_with("chflow-piso-foam"));
        assert!(resolve_three_d_engine_named("chflow-simple-foam").ends_with("chflow-simple-foam"));
    }

    #[test]
    fn resolve_default_3d_engine_points_at_piso() {
        let p = resolve_default_3d_engine();
        assert!(p.ends_with("chflow-piso-foam"), "got {p}");
        // developer tree fallback must exist when this checkout has a CMake
        // build; an isolated worktree (no build/) still resolves the right
        // name, mirroring the resolve_default_engines test convention.
        let dev_root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap().join("build");
        if dev_root.exists() {
            assert!(Path::new(&p).exists(), "dev 3D engine {p} must exist");
        }
    }

    /// E5: the Rust port of the Beggs-Brill boundaries must agree with the
    /// L1..L4 the C++ engine prints for every pipe of the shipped flow-pattern
    /// case — this is the runtime reconciliation the frontend copy lacks.
    #[test]
    fn beggs_brill_boundaries_match_engine_report() {
        let (steady, _t, _c, _tmpl, _tf, _n, _p) = resolve_default_engines();
        if !Path::new(&steady).exists() {
            return; // real-engine-if-available convention (packaged layout w/o 1D exe)
        }
        let case = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .join("tests/cases/steadyFlowPatternMap/case.dict");
        if !case.exists() {
            return;
        }
        let output = std::process::Command::new(&steady)
            .arg(&case)
            .output()
            .expect("spawn chflow-1d-steady");
        assert!(output.status.success());
        let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
        let mut compared = 0usize;
        for line in stdout.lines() {
            // [chflow] pipe <id> ... lambda <l> Nfr <nfr> L1 <a> L2 <b> L3 <c> L4 <d> dp ...
            let grab = |key: &str| -> Option<f64> {
                let idx = line.find(&format!(" {key} "))? + key.len() + 2;
                line[idx..]
                    .split(|c: char| c == ' ')
                    .next()
                    .and_then(|w| w.parse::<f64>().ok())
            };
            let (Some(lam), Some(l1), Some(l2), Some(l3), Some(l4)) = (
                grab("lambda"), grab("L1"), grab("L2"), grab("L3"), grab("L4"),
            ) else {
                continue;
            };
            let rust = beggs_brill_boundaries_impl(lam);
            for (got, want) in [(l1, rust[0]), (l2, rust[1]), (l3, rust[2]), (l4, rust[3])] {
                let scale = want.abs().max(1.0);
                assert!(
                    (got - want).abs() <= 1e-9 * scale,
                    "pipe line '{line}': engine L={got} vs rust {want}"
                );
            }
            compared += 1;
        }
        assert!(
            compared >= 5,
            "expected at least 5 regime report lines, got {compared}"
        );
    }

    // -- job registry -------------------------------------------------------

    #[cfg(unix)]
    fn fake_script(name: &str, body: &str) -> std::path::PathBuf {
        use std::os::unix::fs::PermissionsExt;
        let path = std::env::temp_dir().join(format!(
            "chflow_gui_{}_{}.sh",
            name,
            unique_tmp_id()
        ));
        std::fs::write(&path, body).unwrap();
        let mut perms = std::fs::metadata(&path).unwrap().permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(&path, perms).unwrap();
        path
    }

    #[cfg(unix)]
    #[test]
    fn job_registry_start_poll_finish() {
        let script = fake_script(
            "jobok",
            "#!/bin/sh\necho '[chflow] node A p 30000'\necho done-line\n",
        );
        let reg = JobRegistry::default();
        let id = spawn_job(&reg, "steady", script.to_str().unwrap(), "case.dict")
            .expect("job starts");
        assert_eq!(id, 0);
        // poll until finished (bounded so a broken pipe cannot hang the test)
        let mut status = None;
        for _ in 0..200 {
            status = poll_job_inner(&reg, id);
            if status.as_ref().map(|s| s.finished).unwrap_or(false) {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(20));
        }
        let status = status.expect("job exists");
        assert!(status.finished, "job must finish");
        assert_eq!(status.exit, Some(0));
        assert!(status.output.contains("node A p 30000"), "tail: {}", status.output);
        assert!(status.output.contains("done-line"));
        let _ = std::fs::remove_file(&script);
    }

    #[cfg(unix)]
    #[test]
    fn job_registry_stop_kills_running_child() {
        let script = fake_script("jobsleep", "#!/bin/sh\nsleep 30\n");
        let reg = JobRegistry::default();
        let id = spawn_job(&reg, "transient", script.to_str().unwrap(), "case.dict")
            .expect("job starts");
        std::thread::sleep(std::time::Duration::from_millis(150));
        let stopped = stop_job_inner(&reg, id).expect("stop works");
        assert!(stopped, "running job is stopped");
        // the stop command wrapper drops the registry entry; emulate it here
        drop_job_inner(&reg, id);
        // the inner poll would now return None
        assert!(poll_job_inner(&reg, id).is_none());
        let _ = std::fs::remove_file(&script);
    }

    #[cfg(unix)]
    #[test]
    fn job_registry_poll_unknown_id_is_none() {
        let reg = JobRegistry::default();
        assert!(poll_job_inner(&reg, 999).is_none());
        assert!(stop_job_inner(&reg, 999).is_err());
    }

    // -- project file io ----------------------------------------------------

    #[test]
    fn project_file_round_trip() {
        let path = std::env::temp_dir().join(format!(
            "chflow_gui_project_{}.json",
            unique_tmp_id()
        ));
        write_text_file(path.to_string_lossy().into_owned(), "{\"model\":1}".into())
            .expect("write works");
        let back = read_text_file(path.to_string_lossy().into_owned()).expect("read works");
        assert_eq!(back, "{\"model\":1}");
        let _ = std::fs::remove_file(&path);
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

    #[cfg(unix)]
    #[test]
    fn run_transient_texts_keeps_order_and_reports_errors() {
        use std::fs;
        use std::os::unix::fs::PermissionsExt;

        // The fake engine echoes each case's endTime as its last sample time,
        // so a positional mix-up is visible even when workers finish in a
        // different order than the input rows.
        let script = std::env::temp_dir().join(format!(
            "chflow_gui_fake_transient_batch_{}_{}.sh", std::process::id(), unique_tmp_id()));
        fs::write(&script, "#!/bin/sh\n\
            v=$(grep -o 'endTime [0-9.]*' \"$1\" | head -1 | cut -d' ' -f2)\n\
            cat <<EOF\n\
            transient solved: steps=10\n\
            t 0\n\
            node A p 1 alpha 1 hydrate 0\n\
            t $v\n\
            node A p 2 alpha 1 hydrate 0\n\
            pipe P Q 3\n\
            EOF\n").unwrap();
        let mut perms = fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).unwrap();

        let case = |end: f64| format!(
            "nodes\n{{\n A {{ pressure 1; }}\n}}\npipes\n{{\n P {{ from A; to A; }}\n}}\n\
             time {{ dt 0.1; endTime {end}; }}\n");
        let texts = vec![case(1.0), case(2.0), case(3.0), case(4.0), case(5.0)];
        let rows = run_transient_texts(texts, script.to_string_lossy().into_owned());
        assert_eq!(rows.len(), 5);
        for (i, row) in rows.iter().enumerate() {
            assert!(row.ok, "row {i} failed: {:?}", row.error);
            let s = row.series.as_ref().unwrap();
            assert_eq!(s.samples.len(), 2);
            assert_eq!(s.samples[1].t, (i + 1) as f64, "row {i} out of order");
            assert_eq!(s.samples[1].pressure["A"], 2.0);
        }
        // One invalid row must not poison the others; its slot carries the
        // engine-side rejection verbatim.
        let bad = "nodes\n{\n A {}\n}\npipes\n{\n P { from A; to A; }\n}\n";
        let mut texts = vec![case(9.0), bad.to_string(), case(8.0)];
        let rows = run_transient_texts(std::mem::take(&mut texts), script.to_string_lossy().into_owned());
        assert!(rows[0].ok);
        assert!(!rows[1].ok);
        assert!(rows[1].error.as_ref().unwrap().contains("time{}"));
        assert!(rows[2].ok);
        assert_eq!(rows[2].series.as_ref().unwrap().samples[1].t, 8.0);
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

    /// Real-engine what-if: the SAME case with two valve-close schedules must
    /// give two genuinely different transient series. The case text is read
    /// from the repository's CTest-verified transientValve regression case,
    /// not a hand-written stand-in.
    #[test]
    fn run_transient_texts_real_engine_what_if() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
        let exe = root.join("build/src/1d/apps/chflow-1d-transient");
        if !exe.exists() {
            return;
        }
        let base = std::fs::read_to_string(
            root.join("tests/cases/transientValve/case.dict")).expect("case file");
        let case = |start: f64| base.replace("valveClose { start 2;", &format!(
            "valveClose {{ start {start};"));
        let texts = vec![case(2.0), case(5.0)];
        let rows = run_transient_texts(texts, exe.to_string_lossy().into_owned());
        assert_eq!(rows.len(), 2);
        for (i, row) in rows.iter().enumerate() {
            assert!(row.ok, "row {i} failed: {:?}", row.error);
            assert!(row.series.as_ref().unwrap().samples.len() > 2);
        }
        let a = rows[0].series.as_ref().unwrap();
        let b = rows[1].series.as_ref().unwrap();
        let last_a = a.samples.last().unwrap();
        let last_b = b.samples.last().unwrap();
        let node = &a.node_names[0];
        assert!((last_a.pressure[node] - last_b.pressure[node]).abs() > 1.0,
            "the two valve schedules must change the tank pressure");
    }

    #[test]
    fn serialize_network_case_round_trip() {
        let case = parse_case_text(SAMPLE_CASE).expect("sample parses");
        let text = serialize_network_case(&case);
        let back = parse_case_text(&text).expect("serialized text re-parses");
        assert_eq!(back, case, "round-trip must preserve every field");
    }

    #[test]
    fn serialize_network_case_round_trip_with_inline_equipment() {
        // valve + pump keys the steadyNetwork kernel reads must survive a
        // parse -> serialize -> parse round trip, and a plain pipe must stay
        // free of device tokens
        let text = "nodes\n{\n A { pressure 1e5; }\n B { pressure 1e5; }\n}\npipes\n{\n P1 { from A; to B; length 1000; diameter 0.05; pumpHead 25; pumpRate 3e-4; pumpShutoffHead 40; }\n P2 { from A; to B; length 800; diameter 0.05; valve 0.012; }\n}\n";
        let case = parse_case_text(text).expect("equipment case parses");
        assert_eq!(case.pipes[0].pump_head, 25.0);
        assert_eq!(case.pipes[0].pump_rate, 3.0e-4);
        assert_eq!(case.pipes[0].pump_shutoff_head, 40.0);
        assert_eq!(case.pipes[0].valve, 0.0);
        assert_eq!(case.pipes[1].valve, 0.012);
        assert_eq!(case.pipes[1].pump_head, 0.0);
        let back = parse_case_text(&serialize_network_case(&case)).expect("re-parse");
        assert_eq!(back, case, "设备键必须无损往返");
        let plain = parse_case_text(
            "nodes\n{\n A {};\n}\npipes\n{\n P { from A; to A; }\n}\n",
        )
        .unwrap();
        assert!(
            !serialize_network_case(&plain).contains("valve"),
            "无设备管段不得输出 valve 键"
        );
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

    #[test]
    fn device_families_roundtrip_and_stay_absent_when_zero() {
        // Gas-lift valve (§5.89), rod pump (§5.90) and ESP advanced intake —
        // every key the steadyNetwork kernel reads must survive the GUI
        // round-trip, and a plain pipe must stay byte-identical.
        let text = "nodes\n{\n A {}\n B {}\n}\npipes\n{\n \
GL { from A; to B; length 100; diameter 0.05; liftDomePressure 1.6e7; liftAreaRatio 0.31; liftSpringPressure 9e6; liftDischargeCoeff 0.85; liftKappa 1.4; }\n \
RP { from A; to B; length 800; diameter 0.062; rodPlungerDiameter 0.038; rodStrokeLength 2.5; rodStrokesPerMinute 6; rodDiameter 0.0222; rodLength 1200; rodYoungModulus 2.1e11; rodDensity 7850; }\n \
ESP { from A; to B; length 600; diameter 0.05; pumpHead 800; pumpRate 0.01; pumpGasFraction 0.25; pumpIntakePressure 3.5e6; }\n \
PLAIN { from A; to B; length 100; diameter 0.05; }\n}\n";
        let case = parse_case_text(text).expect("parse");
        let gl = &case.pipes[0];
        assert_eq!(gl.lift_dome_pressure, 1.6e7);
        assert_eq!(gl.lift_area_ratio, 0.31);
        assert_eq!(gl.lift_spring_pressure, 9e6);
        assert_eq!(gl.lift_discharge_coeff, 0.85);
        assert_eq!(gl.lift_kappa, 1.4);
        let rp = &case.pipes[1];
        assert_eq!(rp.rod_plunger_diameter, 0.038);
        assert_eq!(rp.rod_strokes_per_minute, 6.0);
        assert_eq!(rp.rod_young_modulus, 2.1e11);
        assert_eq!(rp.rod_density, 7850.0);
        let esp = &case.pipes[2];
        assert_eq!(esp.pump_gas_fraction, 0.25);
        assert_eq!(esp.pump_intake_pressure, 3.5e6);
        assert_eq!(esp.pump_head, 800.0);

        let round = parse_case_text(&serialize_network_case(&case)).expect("re-parse");
        assert_eq!(round, case);
        let text2 = serialize_network_case(&round);
        assert!(text2.contains("liftDomePressure"));
        assert!(text2.contains("rodPlungerDiameter"));
        assert!(text2.contains("pumpGasFraction"));
        assert!(text2.contains("pumpIntakePressure"));
        // Plain pipe: no device keys leak into its serialization.
        assert!(text2.contains(
            "PLAIN { from A; to B; length 100; diameter 0.05; roughness 0; elevation 0; }"
        ));
    }

    #[cfg(unix)]
    #[test]
    fn run_case_texts_runs_every_row_and_keeps_order() {
        use std::fs;
        use std::os::unix::fs::PermissionsExt;

        let script = std::env::temp_dir().join(format!(
            "chflow_gui_fake_batch_{}_{}.sh", std::process::id(), unique_tmp_id()));
        // The fake engine echoes the case name it was handed so the test can
        // assert that result row i came from text i (order preservation under
        // concurrency), and still prints a parseable summary line.
        fs::write(
            &script,
            "#!/bin/sh\ncase_file=\"$1\"\n\
             nodes=$(grep -c 'pressure' \"$case_file\" || true)\n\
             echo \"[chflow] node A p 30000\"\n\
             echo \"[chflow] steady solved: iterations=$nodes residual=1e-9 (converged)\"\n",
        )
        .unwrap();
        let mut perms = fs::metadata(&script).unwrap().permissions();
        perms.set_mode(0o755);
        fs::set_permissions(&script, perms).unwrap();

        // Five texts with 1..5 fixed-pressure nodes: row i must report
        // iterations = i (the fake engine counts `pressure` occurrences),
        // proving results map back to their inputs, not to completion order.
        let texts: Vec<String> = (1..=5)
            .map(|k| {
                // One node per line: the fake engine counts *lines* containing
                // `pressure` (grep -c), so row i must report iterations = i.
                let mut nodes = String::new();
                for j in 0..k {
                    nodes.push_str(&format!(" F{j} {{ pressure 1e5; }}\n"));
                }
                format!(
                    "nodes\n{{\n{nodes}}}\npipes\n{{\n P {{ from F0; to F0; }}\n}}\n"
                )
            })
            .collect();
        let rows = run_case_texts(texts, script.to_string_lossy().into_owned());
        assert_eq!(rows.len(), 5);
        for (i, row) in rows.iter().enumerate() {
            assert!(row.ok, "row {i} failed: {:?}", row.error);
            let s = row.summary.as_ref().expect("summary present");
            assert_eq!(s.iterations, (i + 1) as u64, "row {i} order");
        }
        let _ = fs::remove_file(&script);

        // Missing engine: every row fails loudly, none panics.
        let rows = run_case_texts(
            vec!["nodes\n{\n A {}\n}\n".into(), "nodes\n{\n B {}\n}\n".into()],
            "/definitely/not/a/real/chflow-binary".into(),
        );
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| !r.ok && r.error.is_some()));
    }
}
