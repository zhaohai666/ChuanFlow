/**
 * Workspace-level pure helpers: the project tree model, the 1D<->3D coupled
 * solver output parser and the ASCII .vtu field reader/renderer. No DOM and
 * no Tauri imports, so node --test covers every branch.
 */

import { fmtDisplay } from './units.js';

// ---------------------------------------------------------------------------
// Project tree
// ---------------------------------------------------------------------------

/**
 * Build the hierarchical project tree shown in the left dock. Every leaf
 * carries a selection target { kind, ref } the property panel understands.
 * `extras` may carry { hasSteady, hasTransient, hasCoupled, coupledEdge }.
 */
export function buildProjectTree(model, extras = {}) {
  const tree = [];

  const fluid = model.fluid || { rho: 850, mu: 0.1 };
  tree.push({
    id: 'fluid',
    label: `流体物性  (ρ ${fluid.rho} · μ ${fluid.mu})`,
    kind: 'group',
    children: [
      { id: 'fluid.rho', label: `密度 ρ = ${fluid.rho} kg/m³`, kind: 'fluid', ref: 'rho' },
      { id: 'fluid.mu', label: `黏度 μ = ${fluid.mu} Pa·s`, kind: 'fluid', ref: 'mu' },
    ],
  });

  const nodeLeaves = model.nodes.map((n) => ({
    id: `node:${n.id}`,
    label:
      n.name +
      (n.pressure !== null && n.pressure !== undefined
        ? `  · 定压`
        : '  · 自由'),
    kind: 'node',
    ref: n.id,
  }));
  const pipeLeaves = model.pipes.map((p) => ({
    id: `pipe:${p.id}`,
    label: `${p.name}  ${p.from} → ${p.to}`,
    kind: 'pipe',
    ref: p.id,
  }));
  const sourceLeaves = model.sources.map((s) => ({
    id: `source:${s.name}`,
    label: `${s.name}  @ ${s.node}  ${s.rate >= 0 ? '+' : ''}${s.rate} m³/s`,
    kind: 'source',
    ref: s.name,
  }));
  const t = model.time || {};
  tree.push({
    id: 'net1d',
    label: `1D 生产管网  (${model.nodes.length} 节点 / ${model.pipes.length} 管段)`,
    kind: 'group',
    children: [
      { id: 'group.nodes', label: `井口/节点 (${model.nodes.length})`, kind: 'group', children: nodeLeaves },
      { id: 'group.pipes', label: `集输管线 (${model.pipes.length})`, kind: 'group', children: pipeLeaves },
      { id: 'group.sources', label: `源/汇 (${model.sources.length})`, kind: 'group', children: sourceLeaves },
      {
        id: 'time',
        label: `瞬态控制  dt ${t.dt} / endTime ${t.endTime}`,
        kind: 'time',
      },
    ],
  });

  if (extras.coupledEdge) {
    tree.push({
      id: 'coupled',
      label: `耦合接口 (1D ↔ 3D)`,
      kind: 'group',
      children: [
        {
          id: 'coupled.edge',
          label: `耦合边 ${extras.coupledEdge}：1D 管段 → 3D 详细的域`,
          kind: 'coupled',
        },
      ],
    });
  }

  const monitorLeaves = [];
  for (const n of model.nodes) {
    monitorLeaves.push({
      id: `mon.node:${n.id}`,
      label: `${n.name} · 压力 p`,
      kind: 'monitor-node',
      ref: n.id,
    });
  }
  for (const p of model.pipes) {
    monitorLeaves.push({
      id: `mon.pipe:${p.id}`,
      label: `${p.name} · 流量 Q`,
      kind: 'monitor-pipe',
      ref: p.id,
    });
  }
  tree.push({
    id: 'monitor',
    label: `监测变量  (${monitorLeaves.length})`,
    kind: 'group',
    children: monitorLeaves,
  });

  const resultLeaves = [];
  resultLeaves.push({
    id: 'res.steady',
    label: extras.hasSteady ? '稳态结果 ✓' : '稳态结果（未求解）',
    kind: 'result-steady',
    ok: !!extras.hasSteady,
  });
  resultLeaves.push({
    id: 'res.transient',
    label: extras.hasTransient ? '瞬态时序结果 ✓' : '瞬态时序结果（未求解）',
    kind: 'result-transient',
    ok: !!extras.hasTransient,
  });
  resultLeaves.push({
    id: 'res.coupled',
    label: extras.hasCoupled ? '耦合结果 ✓' : '耦合结果（未求解）',
    kind: 'result-coupled',
    ok: !!extras.hasCoupled,
  });
  tree.push({ id: 'results', label: '仿真结果', kind: 'group', children: resultLeaves });

  return tree;
}

// ---------------------------------------------------------------------------
// 1D<->3D coupled solver output
// ---------------------------------------------------------------------------

function after(text, key) {
  const pos = text.indexOf(key);
  if (pos < 0) return null;
  return text.slice(pos + key.length).trim();
}

/**
 * Parse chflow-coupled-1d3d stdout: "node final p", "pipe final ... Q ... dp",
 * "coupled solved: edge ... q ... dp ... R ... iterations ... residual ..." and
 * the per-seed outer-loop lines ("coupled seed s: iterations ... converged/
 * not converged", "interface seed s iter k ..."). Tolerates the [chflow]
 * prefix; unknown lines are ignored.
 */
export function parseCoupledOutput(out) {
  const clean = out
    .split('\n')
    .map((l) => l.replace('[chflow] ', ''))
    .join('\n');
  const result = {
    nodes: [], // { name, p }
    pipes: [], // { name, q, v, re, dp }
    solved: null, // { edge, q, dp, r, iterations, residual }
    seeds: [], // { seed, iterations, converged, dp, q }
    interfaces: [], // { seed, iteration, q, dp, resistance }
  };

  for (const raw of clean.split('\n')) {
    const line = raw.trim();
    let m = /^node final p (\S+)$/.exec(line);
    if (m) {
      result.nodes.push({ name: `#${result.nodes.length + 1}`, p: Number(m[1]) });
      continue;
    }
    m = /^node final p (\S+) (\S+)$/.exec(line);
    if (m) {
      // variant with explicit node name
      result.nodes[result.nodes.length - 1] = { name: m[1], p: Number(m[2]) };
      continue;
    }
    m = /^pipe final (\S+) Q (\S+) v (\S+) Re (\S+) dp (\S+)$/.exec(line);
    if (m) {
      result.pipes.push({
        name: m[1], q: Number(m[2]), v: Number(m[3]), re: Number(m[4]), dp: Number(m[5]),
      });
      continue;
    }
    m = /^coupled solved: edge (\S+) q (\S+) dp (\S+) R (\S+) iterations (\d+) residual (\S+)$/.exec(line);
    if (m) {
      result.solved = {
        edge: m[1], q: Number(m[2]), dp: Number(m[3]), r: Number(m[4]),
        iterations: Number(m[5]), residual: Number(m[6]),
      };
      continue;
    }
    m = /^coupled seed (\S+): iterations (\d+)(.*)$/.exec(line);
    if (m) {
      result.seeds.push({
        seed: Number(m[1]),
        iterations: Number(m[2]),
        converged: m[3].includes('converged') && !m[3].includes('not converged'),
      });
      continue;
    }
    m = /^interface seed (\S+) iter (\d+) (.*)$/.exec(line);
    if (m) {
      const rest = m[3];
      const grab = (key) => {
        const a = after(rest, `${key} `);
        return a ? Number(a.split(/\s+/)[0]) : NaN;
      };
      result.interfaces.push({
        seed: Number(m[1]),
        iteration: Number(m[2]),
        q: grab('q'),
        dp: grab('dp'),
        resistance: grab('R'),
      });
    }
  }
  return result;
}

// ---------------------------------------------------------------------------
// ASCII .vtu field reader (cf::fvm::VtkWriter output)
// ---------------------------------------------------------------------------

/**
 * Parse an ASCII .vtu written by VtkWriter: points (x y z), polygonal
 * connectivity + offsets, and the cell data arrays U (3 comps) and p.
 * Returns null when the text is not a readable VTU.
 */
export function parseVtu(text) {
  try {
    const ptsMatch = /<DataArray type="Float64" NumberOfComponents="3" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    if (!ptsMatch) return null;
    const coords = ptsMatch[1].trim().split(/\s+/).map(Number);
    const points = [];
    for (let i = 0; i + 2 < coords.length; i += 3) {
      points.push({ x: coords[i], y: coords[i + 1], z: coords[i + 2] });
    }

    const connMatch = /Name="connectivity" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    const offMatch = /Name="offsets" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    if (!connMatch || !offMatch) return null;
    const conn = connMatch[1].trim().split(/\s+/).map(Number);
    const offsets = offMatch[1].trim().split(/\s+/).map(Number);

    const uMatch = /<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    const pMatch = /<DataArray type="Float64" Name="p" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    const U = uMatch ? uMatch[1].trim().split(/\s+/).map(Number) : [];
    const p = pMatch ? pMatch[1].trim().split(/\s+/).map(Number) : [];

    const cells = [];
    let prev = 0;
    let k = 0;
    const cellU = [];
    const cellP = [];
    for (const o of offsets) {
      const vids = conn.slice(prev, o);
      prev = o;
      cells.push(vids);
      cellU.push(U.length >= 3 * (k + 1) ? [U[3 * k], U[3 * k + 1], U[3 * k + 2]] : [0, 0, 0]);
      cellP.push(p.length > k ? p[k] : 0);
      k += 1;
    }
    return { points, cells, cellU, cellP };
  } catch (_) {
    return null;
  }
}

/** Cell centroid from vertex list. */
export function cellCentroid(points, vids) {
  let x = 0, y = 0, z = 0;
  for (const v of vids) {
    x += points[v].x;
    y += points[v].y;
    z += points[v].z;
  }
  const n = vids.length || 1;
  return { x: x / n, y: y / n, z: z / n };
}

/** [min, max] of a numeric array (NaN-safe). */
export function fieldRange(values) {
  let lo = Infinity;
  let hi = -Infinity;
  for (const v of values) {
    if (!Number.isFinite(v)) continue;
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  if (lo === Infinity) return [0, 1];
  return [lo, hi];
}

/**
 * Classic CFD rainbow colormap: value in [min,max] -> css `rgb(...)`.
 * blue (low) -> cyan -> green -> yellow -> red (high).
 */
export function rainbowColor(value, min, max) {
  let t = max > min ? (value - min) / (max - min) : 0;
  t = Math.max(0, Math.min(1, t));
  const stops = [
    [0, 0, 114, 189],
    [0.25, 36, 192, 148],
    [0.5, 117, 205, 75],
    [0.75, 245, 200, 47],
    [1, 214, 47, 39],
  ];
  for (let i = 1; i < stops.length; i += 1) {
    if (t <= stops[i][0]) {
      const [t0, r0, g0, b0] = stops[i - 1];
      const [t1, r1, g1, b1] = stops[i];
      const f = (t - t0) / (t1 - t0 || 1);
      const r = Math.round(r0 + f * (r1 - r0));
      const g = Math.round(g0 + f * (g1 - g0));
      const b = Math.round(b0 + f * (b1 - b0));
      return `rgb(${r},${g},${b})`;
    }
  }
  return 'rgb(214,47,39)';
}

/**
 * Log panel line classifier: returns 'error' | 'warn' | 'ok' | 'info'.
 * Kernel-facing keywords only - no invented vocabulary.
 */
export function classifyLogLine(line) {
  const l = line.toLowerCase();
  if (l.includes('error') || l.includes('fatal') || l.includes('failed') || l.includes('nan')) {
    return 'error';
  }
  if (l.includes('warn') || l.includes('not converged')) return 'warn';
  if (l.includes('converged') || l.includes('solved') || l.includes('ok:')) return 'ok';
  return 'info';
}

// ---------------------------------------------------------------------------
// Engine stdout parsers (JS mirrors of the Rust run_case / run_transient_text
// parsers) - the job-based runs stream raw text, so the workspace parses it.
// ---------------------------------------------------------------------------

/** Parse chflow-1d-steady stdout: "node NAME p V", "pipe NAME Q v Re f dp",
 * "steady solved: iterations=N residual=R (converged)". */
export function parseSteadyOutput(out) {
  const summary = { nodes: [], pipes: [], iterations: 0, residual: 0, converged: false };
  for (const raw of out.split('\n')) {
    const line = raw.trim().replace('[chflow] ', '').trim();
    if (!line) continue;
    let m = /^node (\S+) p (\S+)$/.exec(line);
    if (m) {
      summary.nodes.push({ name: m[1], p: Number(m[2]) });
      continue;
    }
    m = /^pipe (\S+) Q (\S+) v (\S+) Re (\S+)(.*)$/.exec(line);
    if (m) {
      const rest = m[5];
      const dpMatch = /\bdp (\S+)/.exec(rest);
      summary.pipes.push({
        name: m[1], q: Number(m[2]), v: Number(m[3]), re: Number(m[4]),
        dp: dpMatch ? Number(dpMatch[1]) : 0,
      });
      continue;
    }
    if (line.includes('steady solved:')) {
      const it = /iterations=(\d+)/.exec(line);
      const res = /residual=(\S+)/.exec(line);
      summary.iterations = it ? Number(it[1]) : 0;
      summary.residual = res ? Number(res[1]) : 0;
      summary.converged = line.includes('(converged)');
    }
  }
  return summary;
}

/** Parse chflow-1d-transient stdout into a time-ordered series (same shape
 * as the Rust parse_transient_output: per sample name->value maps). */
export function parseTransientOutput(out) {
  const series = {
    nodeNames: [], pipeNames: [], dt: 0, endTime: 0, steps: 0, samples: [],
  };
  let cur = null;
  for (const raw of out.split('\n')) {
    const line = raw.trim().replace('[chflow] ', '').trim();
    if (!line) continue;
    if (line.startsWith('transient:')) {
      const dt = /dt=(\S+)/.exec(line);
      const end = /endTime=(\S+)/.exec(line);
      series.dt = dt ? Number(dt[1]) : 0;
      series.endTime = end ? Number(end[1]) : 0;
      continue;
    }
    if (line.startsWith('t ')) {
      if (cur) series.samples.push(cur);
      const t = Number(line.slice(2).trim().split(/\s+/)[0]);
      cur = Number.isFinite(t)
        ? { t, pressure: {}, alpha: {}, hydrate: {}, flow: {} }
        : null;
      continue;
    }
    let m = /^node (\S+) p (\S+) alpha (\S+) hydrate (\S+)$/.exec(line);
    if (m && cur) {
      if (!series.nodeNames.includes(m[1])) series.nodeNames.push(m[1]);
      cur.pressure[m[1]] = Number(m[2]);
      cur.alpha[m[1]] = Number(m[3]);
      cur.hydrate[m[1]] = Number(m[4]);
      continue;
    }
    m = /^pipe (\S+) Q (\S+)$/.exec(line);
    if (m && cur) {
      if (!series.pipeNames.includes(m[1])) series.pipeNames.push(m[1]);
      cur.flow[m[1]] = Number(m[2]);
      continue;
    }
    if (line.startsWith('transient solved:')) {
      const st = /steps=(\d+)/.exec(line);
      series.steps = st ? Number(st[1]) : 0;
    }
  }
  if (cur) series.samples.push(cur);
  return series;
}

/**
 * Extract the fields the drag editor model keeps (fluid, time controls and
 * per-pipe geometry) from a case dictionary, so applying a dict rebuilds a
 * FAITHFUL model instead of silently defaulting geometry. Pure.
 */
export function extractCaseDetails(text) {
  const details = {
    fluid: null,
    time: null,
    pipes: {}, // name -> { length, diameter, roughness, elevation }
  };

  const grabKv = (body) => {
    const kv = {};
    for (const stmt of body.split(';')) {
      const t = stmt.trim();
      if (!t) continue;
      const parts = t.split(/\s+/);
      kv[parts[0]] = parts.slice(1).join(' ');
    }
    return kv;
  };

  const fluid = /\bfluid\s*\{([^}]*)\}/.exec(text);
  if (fluid) {
    const kv = grabKv(fluid[1]);
    details.fluid = {
      rho: Number(kv.rho),
      mu: Number(kv.mu),
    };
  }
  const time = /\btime\s*\{([^}]*)\}/.exec(text);
  if (time) {
    const kv = grabKv(time[1]);
    details.time = {
      dt: Number(kv.dt),
      endTime: Number(kv.endTime),
      reportEvery: Math.round(Number(kv.reportEvery || 1)),
    };
  }

  // pipe entries wherever they appear: NAME { from X; to Y; length ...; }.
  // Single-brace-depth bodies only: the section header `pipes {` has nested
  // entry braces, so [^{}]* naturally skips it and lands on the entries.
  const pipeRe = /([A-Za-z_]\w*)\s*\{([^{}]*)\}/g;
  for (const m of text.matchAll(pipeRe)) {
    if (!/\bfrom\b/.test(m[2])) continue;
    const kv = grabKv(m[2]);
    if (kv.from === undefined) continue;
    details.pipes[m[1]] = {
      length: kv.length !== undefined ? Number(kv.length) : 1000,
      diameter: kv.diameter !== undefined ? Number(kv.diameter) : 0.05,
      roughness: kv.roughness !== undefined ? Number(kv.roughness) : 0,
      elevation: kv.elevation !== undefined ? Number(kv.elevation) : 0,
    };
  }
  return details;
}
