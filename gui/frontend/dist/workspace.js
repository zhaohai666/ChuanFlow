/**
 * Workspace-level pure helpers: the project tree model, the 1D<->3D coupled
 * solver output parser and the ASCII .vtu field reader/renderer. No DOM and
 * no Tauri imports, so node --test covers every branch.
 */

import { fmtDisplay } from './units.js';
import { ALL_DEVICE_FIELDS } from './editor.js';

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

  if (extras.scenarios && extras.scenarios.length) {
    tree.push({
      id: 'group.scenarios',
      label: `工况 (${extras.scenarios.length})`,
      kind: 'group',
      children: extras.scenarios.map((s) => ({
        id: `scen:${s.id}`,
        label: `${s.name} · ${s.engine === 'transient' ? '瞬态' : '稳态'} · ${s.nOverrides} 项`,
        kind: 'scenario',
        ref: s.id,
      })),
    });
  }

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
    // VtkWriter's faces array is nested: per cell [nFaces, nPts, ids..., ...].
    // It is optional for field display but required for report definitions
    // (surface/volume integrals), so parse it when present.
    const faceMatch = /Name="faces" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    const faceData = faceMatch ? faceMatch[1].trim().split(/\s+/).map(Number) : [];

    const uMatch = /<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    const pMatch = /<DataArray type="Float64" Name="p" format="ascii">\s*([\s\S]*?)<\/DataArray>/.exec(text);
    const U = uMatch ? uMatch[1].trim().split(/\s+/).map(Number) : [];
    const p = pMatch ? pMatch[1].trim().split(/\s+/).map(Number) : [];

    const cells = [];
    const cellFaces = [];
    let prev = 0;
    let k = 0;
    let fi = 0;
    const cellU = [];
    const cellP = [];
    for (const o of offsets) {
      const vids = conn.slice(prev, o);
      prev = o;
      cells.push(vids);
      if (faceData.length) {
        const rings = [];
        const nf = faceData[fi++] || 0;
        for (let f = 0; f < nf; f += 1) {
          const nv = faceData[fi++] || 0;
          rings.push(faceData.slice(fi, fi + nv));
          fi += nv;
        }
        cellFaces.push(rings);
      }
      cellU.push(U.length >= 3 * (k + 1) ? [U[3 * k], U[3 * k + 1], U[3 * k + 2]] : [0, 0, 0]);
      cellP.push(p.length > k ? p[k] : 0);
      k += 1;
    }
    return { points, cells, cellFaces, cellU, cellP };
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

/**
 * Streamlines through a cell-centred velocity field (U7 / §5.138). The cell
 * centres are bucketed into a uniform grid so each Euler step samples the
 * nearest cell in O(1) instead of scanning every cell; integration uses a
 * normalised direction and a fraction of the local cell size, and stops on
 * leaving the bounding box (one cell of slack), a vanishing velocity, a step
 * budget, or returning near the seed (closed recirculation). Pure so node
 * --test can assert straight-line and circular-field behaviour.
 */
export function traceStreamlines(points, cells, cellU, seeds, opts = {}) {
  if (!cells.length || !seeds.length) return [];
  const maxSteps = opts.maxSteps ?? 120;
  const stepScale = opts.stepScale ?? 0.2;
  const speedEps = opts.speedEps ?? 1e-12;
  const gridN = opts.gridN ?? 12;
  const centres = cells.map((vids) => cellCentroid(points, vids));
  let xlo = Infinity, ylo = Infinity, zlo = Infinity;
  let xhi = -Infinity, yhi = -Infinity, zhi = -Infinity;
  for (const c of centres) {
    xlo = Math.min(xlo, c.x); xhi = Math.max(xhi, c.x);
    ylo = Math.min(ylo, c.y); yhi = Math.max(yhi, c.y);
    zlo = Math.min(zlo, c.z); zhi = Math.max(zhi, c.z);
  }
  const hx = (xhi - xlo) / gridN || 1;
  const hy = (yhi - ylo) / gridN || 1;
  const hz = (zhi - zlo) / gridN || 1;
  const ds = stepScale * Math.min(hx, hy, hz);
  const ix = (x) => Math.max(0, Math.min(gridN - 1, Math.floor((x - xlo) / hx)));
  const iy = (y) => Math.max(0, Math.min(gridN - 1, Math.floor((y - ylo) / hy)));
  const iz = (z) => Math.max(0, Math.min(gridN - 1, Math.floor((z - zlo) / hz)));
  const buckets = new Map();
  centres.forEach((c, k) => {
    const key = ix(c.x) + gridN * (iy(c.y) + gridN * iz(c.z));
    if (!buckets.has(key)) buckets.set(key, []);
    buckets.get(key).push(k);
  });
  const nearest = (p) => {
    const bx = ix(p.x), by = iy(p.y), bz = iz(p.z);
    let best = -1, bestD = Infinity;
    for (let r = 0; r <= 3; r += 1) {
      for (let a = -r; a <= r; a += 1) for (let b = -r; b <= r; b += 1) for (let c = -r; c <= r; c += 1) {
        if (Math.max(Math.abs(a), Math.abs(b), Math.abs(c)) !== r) continue;
        const gx = bx + a, gy = by + b, gz = bz + c;
        if (gx < 0 || gy < 0 || gz < 0 || gx >= gridN || gy >= gridN || gz >= gridN) continue;
        for (const k of buckets.get(gx + gridN * (gy + gridN * gz)) || []) {
          const c0 = centres[k];
          const d = (c0.x - p.x) ** 2 + (c0.y - p.y) ** 2 + (c0.z - p.z) ** 2;
          if (d < bestD) { bestD = d; best = k; }
        }
      }
      if (best >= 0) break; // nearest bucket ring found; a nearer cell cannot be further out
    }
    return best;
  };

  const lines = [];
  for (const seed of seeds) {
    let p = { x: seed[0], y: seed[1], z: seed[2] };
    const line = [[p.x, p.y, p.z]];
    line.seedCell = -1;
    for (let step = 0; step < maxSteps; step += 1) {
      const k = nearest(p);
      if (k < 0) break;
      if (line.seedCell < 0) line.seedCell = k;
      const u = cellU[k] || [0, 0, 0];
      const mag = Math.hypot(u[0], u[1], u[2]);
      if (!Number.isFinite(mag) || mag <= speedEps) break;
      const f = ds / mag;
      const q = { x: p.x + u[0] * f, y: p.y + u[1] * f, z: p.z + u[2] * f };
      if (q.x < xlo - hx || q.x > xhi + hx || q.y < ylo - hy || q.y > yhi + hy
          || q.z < zlo - hz || q.z > zhi + hz) break;
      p = q;
      line.push([p.x, p.y, p.z]);
      if (step > 4 && Math.hypot(p.x - seed[0], p.y - seed[1], p.z - seed[2]) <= ds) break;
    }
    if (line.length > 1) lines.push(line);
  }
  return lines;
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

/** Newell's area vector of a planar vertex ring; |S| is the face area. */
export function areaVector(points, ring) {
  let nx = 0, ny = 0, nz = 0;
  for (let i = 0; i < ring.length; i += 1) {
    const a = points[ring[i]];
    const b = points[ring[(i + 1) % ring.length]];
    nx += (a.y - b.y) * (a.z + b.z);
    ny += (a.z - b.z) * (a.x + b.x);
    nz += (a.x - b.x) * (a.y + b.y);
  }
  return { x: nx / 2, y: ny / 2, z: nz / 2 };
}

/** Divergence-theorem cell volumes from the VTU face rings (U8 / §5.139). */
export function cellVolumes(points, cells, cellFaces) {
  if (!cellFaces || cellFaces.length !== cells.length) return null;
  return cells.map((vids, c) => {
    const centre = cellCentroid(points, vids);
    let vol = 0;
    for (const ring of cellFaces[c]) {
      if (ring.length < 3) continue;
      const s = areaVector(points, ring);
      const f = cellCentroid(points, ring);
      vol += Math.abs(s.x * (f.x - centre.x) + s.y * (f.y - centre.y)
        + s.z * (f.z - centre.z));
    }
    return vol / 3;
  });
}

/**
 * Evaluate report definitions against a parsed VTU field. Supported fields:
 * p, |U|, Ux/Uy/Uz; operations: volumeIntegral, volumeAverage, max, min,
 * flux (through a unit normal) and surfaceAverage (on the plane described by
 * `normal` and the fractional position along the field's projected extent).
 * All weights come from the mesh's own face geometry — no hardcoded cell size.
 */
export function evaluateReports(vtu, defs) {
  const volumes = cellVolumes(vtu.points, vtu.cells, vtu.cellFaces);
  const n = vtu.cells.length;
  const field = (name, c) => {
    if (name === 'p') return vtu.cellP[c];
    const u = vtu.cellU[c] || [0, 0, 0];
    if (name === '|U|') return Math.hypot(u[0], u[1], u[2]);
    if (name === 'Ux') return u[0];
    if (name === 'Uy') return u[1];
    if (name === 'Uz') return u[2];
    return NaN;
  };
  return defs.map((def) => {
    let value = NaN;
    if (def.op === 'volumeIntegral' || def.op === 'volumeAverage') {
      if (!volumes) return { name: def.name, value: NaN, error: 'VTU has no faces array' };
      let num = 0, den = 0;
      for (let c = 0; c < n; c += 1) {
        const v = field(def.field, c);
        if (!Number.isFinite(v)) continue;
        num += v * volumes[c];
        den += volumes[c];
      }
      value = def.op === 'volumeIntegral' ? num : (den > 0 ? num / den : NaN);
    } else if (def.op === 'max' || def.op === 'min') {
      const vals = [];
      for (let c = 0; c < n; c += 1) {
        const v = field(def.field, c);
        if (Number.isFinite(v)) vals.push(v);
      }
      value = vals.length ? (def.op === 'max' ? Math.max(...vals) : Math.min(...vals)) : NaN;
    } else if (def.op === 'flux' || def.op === 'surfaceAverage') {
      if (!volumes || !vtu.cellFaces) {
        return { name: def.name, value: NaN, error: 'VTU has no faces array' };
      }
      const nv = def.normal || [1, 0, 0];
      const len = Math.hypot(nv[0], nv[1], nv[2]) || 1;
      const nx = nv[0] / len, ny = nv[1] / len, nz = nv[2] / len;
      // plane position from the FACE-centre projection extent, so frac=0/1
      // select the true boundary face layer even in a one-cell-thick mesh
      let lo = Infinity, hi = -Infinity;
      for (const rings of vtu.cellFaces) for (const ring of rings) {
        if (ring.length < 3) continue;
        const f = cellCentroid(vtu.points, ring);
        const p = f.x * nx + f.y * ny + f.z * nz;
        lo = Math.min(lo, p); hi = Math.max(hi, p);
      }
      const target = lo + (def.frac ?? 0.5) * (hi - lo || 1);
      let num = 0, den = 0;
      for (let c = 0; c < n; c += 1) {
        for (const ring of vtu.cellFaces[c]) {
          if (ring.length < 3) continue;
          const s = areaVector(vtu.points, ring);
          const projected = s.x * nx + s.y * ny + s.z * nz;
          if (Math.abs(projected) < 1e-14) continue;
          const f = cellCentroid(vtu.points, ring);
          const fp = f.x * nx + f.y * ny + f.z * nz;
          const tol = Math.cbrt(volumes[c]) * 0.3;
          if (Math.abs(fp - target) > tol) continue;
          const w = Math.abs(projected);
          if (def.op === 'flux') {
            const u = vtu.cellU[c] || [0, 0, 0];
            num += (u[0] * nx + u[1] * ny + u[2] * nz) * w;
            den += 1;
          } else {
            const v = field(def.field, c);
            if (!Number.isFinite(v)) continue;
            num += v * w;
            den += w;
          }
        }
      }
      value = def.op === 'flux' ? num : (den > 0 ? num / den : NaN);
    }
    return { name: def.name, value, error: null };
  });
}

/** Report history -> CSV text (one row per recorded field update). */
export function reportsToCsv(history, defs) {
  const header = ['step', ...defs.map((d) => d.name)].join(',');
  const rows = history.map((h, i) => [i + 1, ...defs.map((d) => {
    const v = h[d.name];
    return Number.isFinite(v) ? String(v) : '';
  })].join(','));
  return [header, ...rows].join('\n') + '\n';
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
      // Inline equipment (network.cpp pipe keys + E4 device families).
      // The registry is shared with editor.js serialization so all three
      // sides (parse / serialize / extract) stay in lockstep. 0 = absent.
      ...Object.fromEntries(
        ALL_DEVICE_FIELDS.map((k) => [k, kv[k] !== undefined ? Number(kv[k]) : 0]),
      ),
    };
  }
  return details;
}


/**
 * Parse chflow-1d-twofluid per-frame output into an animation-ready timeline.
 * A frame starts at the "[chflow] t <t> vm <vm> dp <dp>" line and gathers the
 * per-cell lines that follow; each cell line is "cell <i> alpha <a> water <w>
 * vsl <v> vsg <g> [T <T>] [p <p>] [s <s>]" — the trailing keys depend on the
 * case's features (thermal / EOS / terrain slugs), so every pair is stored by
 * name instead of by position.
 */
export function parseTwoFluidFrames(out) {
  const frames = [];
  let cur = null;
  const tRe = /^t\s+(\S+)\s+vm\s+(\S+)\s+dp\s+(\S+)/;
  const cellRe = /^cell\s+(\d+)\s+(.*)$/;
  const nodeRe = /^node\s+(\d+)\s+p\s+(\S+)$/;
  for (const raw of out.split('\n')) {
    const line = raw.trim().replace('[chflow] ', '').trim();
    if (!line) continue;
    const tm = tRe.exec(line);
    if (tm) {
      if (cur) frames.push(cur);
      cur = {
        t: Number(tm[1]),
        vm: Number(tm[2]),
        dp: Number(tm[3]),
        nodeP: [],
        cells: [],
      };
      continue;
    }
    const nm = nodeRe.exec(line);
    if (nm && cur) {
      cur.nodeP[Number(nm[1])] = Number(nm[2]);
      continue;
    }
    const cm = cellRe.exec(line);
    if (cm && cur) {
      const cell = { index: Number(cm[1]) };
      for (const stmt of cm[2].split(/\s+(?=[A-Za-z])/)) {
        const parts = stmt.trim().split(/\s+/);
        if (parts.length >= 2) {
          const num = Number(parts[1]);
          cell[parts[0]] = Number.isFinite(num) ? num : parts[1];
        }
      }
      cur.cells.push(cell);
    }
  }
  if (cur) frames.push(cur);
  return frames;
}

// ---------------------------------------------------------------------------
// Nodal analysis (U2): IPR declaration, closed-form IPR curve, VLP sampling
// ---------------------------------------------------------------------------

const NODAL_NUM_RE = '[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?';

/** Body + inner span of a top-level `name { ... }` block; null if absent. */
function nodalBlockBody(text, name) {
  const re = new RegExp(`\\b${name}\\s*\\{`);
  const m = re.exec(text);
  if (!m) return null;
  const open = m.index + m[0].length - 1;
  let depth = 0;
  for (let i = open; i < text.length; i += 1) {
    if (text[i] === '{') depth += 1;
    else if (text[i] === '}') {
      depth -= 1;
      if (depth === 0) return { body: text.slice(open + 1, i), start: open + 1, end: i, kwStart: m.index };
    }
  }
  return null;
}

function nodalNum(body, key) {
  const m = new RegExp(`\\b${key}\\s+(${NODAL_NUM_RE})`).exec(body);
  return m ? Number(m[1]) : NaN;
}

/**
 * Parse the ipr{} declaration of a nodal case dictionary into the closed-form
 * spec the kernel's Nodal.h uses (mirrored here only for drawing the curve;
 * all physics stays in the engine). Throws without an ipr{} block.
 */
export function parseNodalIpr(text) {
  const blk = nodalBlockBody(text, 'ipr');
  if (!blk) throw new Error('nodal 字典缺少 ipr{} 块');
  const form = (/\bform\s+(\w+)/.exec(blk.body) || [])[1] || 'linear';
  const pRes = nodalNum(blk.body, 'reservoirPressure');
  if (!Number.isFinite(pRes) || pRes <= 0) throw new Error('ipr{} 缺少有效 reservoirPressure');
  const spec = { form, pRes };
  if (form === 'linear') {
    spec.J = nodalNum(blk.body, 'productivityIndex');
    if (!Number.isFinite(spec.J) || spec.J <= 0) throw new Error('linear ipr 缺少 productivityIndex');
    spec.qAof = pRes * spec.J;
  } else if (form === 'vogel') {
    spec.qmax = nodalNum(blk.body, 'maxRate');
    if (!Number.isFinite(spec.qmax) || spec.qmax <= 0) throw new Error('vogel ipr 缺少 maxRate');
    spec.qAof = spec.qmax;
  } else if (form === 'jones') {
    spec.C = nodalNum(blk.body, 'deliverability');
    spec.n = nodalNum(blk.body, 'exponent');
    if (!Number.isFinite(spec.C) || !Number.isFinite(spec.n) || spec.n <= 0) {
      throw new Error('jones ipr 缺少 deliverability/exponent');
    }
    spec.qAof = spec.C * pRes ** (2 * spec.n);
  } else {
    throw new Error(`未知 ipr form: ${form}`);
  }
  spec.jRef = spec.qAof / pRes;
  return spec;
}

/**
 * Closed-form IPR curve p(q) for the parsed spec (the same three forms as
 * Nodal.h: linear p = pRes - q/J; Vogel q/qmax = 1 - 0.2x - 0.8x^2 with
 * x = p/pRes; Jones q = C (pRes^2 - p^2)^n). n+1 points, q ascending.
 */
export function iprCurvePoints(spec, n = 48) {
  const pts = [];
  const { pRes, form } = spec;
  if (form === 'linear') {
    for (let i = 0; i <= n; i += 1) {
      const q = (spec.qAof * i) / n;
      pts.push({ q, p: pRes - q / spec.J });
    }
  } else if (form === 'vogel') {
    for (let i = 0; i <= n; i += 1) {
      const x = 1 - i / n; // p/pRes descending -> q ascending
      const p = x * pRes;
      pts.push({ q: spec.qmax * (1 - 0.2 * x - 0.8 * x * x), p });
    }
  } else { // jones
    for (let i = 0; i <= n; i += 1) {
      const x = 1 - i / n;
      const p = x * pRes;
      pts.push({ q: spec.C * (pRes ** 2 - p * p) ** spec.n, p });
    }
  }
  return pts;
}

// ---------------------------------------------------------------------------
// Template side-by-side comparison (U13)
// ---------------------------------------------------------------------------

function mergeRunRows(la, lb, field) {
  const map = new Map();
  for (const it of la || []) map.set(it.name, { name: it.name, a: it[field], b: NaN });
  for (const it of lb || []) {
    const e = map.get(it.name);
    if (e) e.b = it[field];
    else map.set(it.name, { name: it.name, a: NaN, b: it[field] });
  }
  return [...map.values()].map((r) => ({ ...r, d: r.a - r.b }));
}

/**
 * Align two steady SolverSummary batches (template vs current model) by
 * entity name: rows carry a = template value, b = current value, d = a - b.
 * Names present in only one run keep NaN on the other side (shown as —).
 */
export function compareSteadyRuns(template, current) {
  return {
    nodes: mergeRunRows(template?.nodes, current?.nodes, 'pressure'),
    pipes: mergeRunRows(template?.pipes, current?.pipes, 'q'),
  };
}

// ---------------------------------------------------------------------------
// PVT fluid wizard (U4): flash case serialization + chflow-pvt report parser
// ---------------------------------------------------------------------------

/** Component names the PR feed accepts (ComponentData built-in table). */
export const PVT_COMPONENTS = [
  'methane', 'ethane', 'propane', 'nButane', 'nPentane', 'nHexane',
  'nHeptane', 'nOctane', 'nNonane', 'nDecane', 'CO2', 'N2', 'H2S',
];

/**
 * Serialize the wizard's spec into a `chflow-pvt` flash dictionary. The
 * composition is declared as-is (the kernel normalizes); viscosity /
 * surfaceTension blocks are gated exactly as the app reads them.
 */
export function pvtFlashCaseText(spec) {
  const comps = (spec.components || []).filter(
    (c) => c.name && Number.isFinite(c.z) && c.z > 0,
  );
  if (!comps.length) throw new Error('PVT：至少需要一个组分（名称 + 摩尔分数 > 0）');
  const ps = (spec.pressures || []).filter(Number.isFinite);
  if (!ps.length) throw new Error('PVT：需要压力网格');
  const T = Number(spec.temperature);
  if (!Number.isFinite(T) || T <= 0) throw new Error('PVT：需要有效温度 (K)');
  const lines = [
    `temperature ${T};`,
    '',
    `pressures (${ps.join(' ')});`,
    '',
    'composition',
    '{',
    ...comps.map((c) => `    ${c.name} ${c.z};`),
    '}',
    '',
    'solve',
    '{',
    '    tol 1.0e-12;',
    '    maxIter 500;',
    `    bubblePoint ${spec.bubble ? 'yes' : 'no'};`,
    '}',
  ];
  if (spec.viscosity) lines.push('', 'viscosity', '{', '    model lbc;', '}');
  if (spec.surfaceTension) {
    lines.push('', 'surfaceTension', '{', '    model weinaugKatz;',
      '    parachors fromTable;', '}');
  }
  return `${lines.join('\n')}\n`;
}

/**
 * Parse the full `chflow-pvt` text report into per-(T,p) flash rows plus the
 * bubble point and any tabulated property. Numbers stay SI; the GUI converts
 * at display time. Throws on output without a single flash row.
 */
export function parsePvtOutput(out) {
  const res = {
    nComponents: 0,
    components: [],
    molarMass: NaN,
    flashes: [],
    bubble: null,
    tables: [],
  };
  let cur = null;
  let m;
  for (const raw of out.split('\n')) {
    const line = raw.trim().replace('[chflow] ', '').trim();
    if ((m = /^pvt:\s+(\d+)\s+components,\s+feed\s+(.*?)\s+molarMass\s+(\S+)$/.exec(line))) {
      res.nComponents = Number(m[1]);
      res.components = m[2].trim().split(/\s+/).map((t) => {
        const [name, z] = t.split('=');
        return { name, z: Number(z) };
      });
      res.molarMass = Number(m[3]);
    } else if ((m = /^pvt flash T (\S+) p (\S+) beta (\S+) single (\S+) iterations (\S+) converged (\S+) residual (\S+)$/.exec(line))) {
      cur = {
        T: Number(m[1]), p: Number(m[2]), beta: Number(m[3]),
        single: m[4] === '1', iterations: Number(m[5]),
        converged: m[6] === '1', residual: Number(m[7]),
        x: {}, y: {},
      };
      res.flashes.push(cur);
    } else if (cur && (m = /^pvt liquid rho (\S+) M (\S+)$/.exec(line))) {
      cur.rhoL = Number(m[1]);
      cur.ML = Number(m[2]);
    } else if (cur && (m = /^pvt vapour rho (\S+) M (\S+)$/.exec(line))) {
      cur.rhoV = Number(m[1]);
      cur.MV = Number(m[2]);
    } else if (cur && (m = /^pvt (x|y) (\S+) (\S+)$/.exec(line))) {
      cur[m[1]][m[2]] = Number(m[3]);
    } else if (cur && (m = /^pvt viscosity liquid (\S+) vapour (\S+)$/.exec(line))) {
      cur.muL = Number(m[1]);
      cur.muV = Number(m[2]);
    } else if (cur && (m = /^pvt surfaceTension (\S+)$/.exec(line))) {
      cur.sigma = Number(m[1]);
    } else if ((m = /^pvt bubble T (\S+) p (\S+)$/.exec(line))) {
      res.bubble = { T: Number(m[1]), p: Number(m[2]) };
    } else if ((m = /^pvt table (\S+) T (\S+) p (\S+) value (\S+)$/.exec(line))) {
      res.tables.push({ prop: m[1], T: Number(m[2]), p: Number(m[3]), value: Number(m[4]) });
    }
  }
  if (!res.flashes.length) throw new Error('非 chflow-pvt 输出（无 pvt flash 行）');
  return res;
}

/**
 * Read an existing flash dictionary back into the wizard spec (template ->
 * 向导回填). Composition entries are `name z;` pairs inside composition{}.
 */
export function parsePvtCaseDict(text) {
  const tm = new RegExp(`\\btemperature\\s+(${NODAL_NUM_RE})`).exec(text);
  const pm = /pressures\s*\(([^)]*)\)/.exec(text);
  const comp = nodalBlockBody(text, 'composition');
  const components = [];
  if (comp) {
    const re = new RegExp(`([A-Za-z_]\\w*)\\s+(${NODAL_NUM_RE})\\s*;`, 'g');
    let m;
    while ((m = re.exec(comp.body)) !== null) components.push({ name: m[1], z: Number(m[2]) });
  }
  return {
    temperature: tm ? Number(tm[1]) : NaN,
    pressures: pm ? pm[1].trim().split(/[\s,]+/).map(Number).filter(Number.isFinite) : [],
    components,
    bubble: /bubblePoint\s+yes/.test(text),
    viscosity: /\bviscosity\s*\{/.test(text),
    surfaceTension: /\bsurfaceTension\s*\{/.test(text),
  };
}

/**
 * Group parsed flashes by temperature for the rho/mu/sigma-vs-p charts:
 * one entry per T, ascending, each carrying the pressure-ascending rows.
 */
export function pvtCurvesByT(parsed) {
  const byT = new Map();
  for (const f of parsed.flashes) {
    if (!byT.has(f.T)) byT.set(f.T, []);
    byT.get(f.T).push(f);
  }
  return [...byT.entries()]
    .sort((a, b) => a[0] - b[0])
    .map(([T, rows]) => ({ T, rows: rows.sort((x, y) => x.p - y.p) }));
}

// ---------------------------------------------------------------------------
// Scenario manager (U1): batch-run row normalization for both engines
// ---------------------------------------------------------------------------

/**
 * Normalize a `run_transient_texts` batch: per-scenario convergence, frame
 * count, end-of-run p/Q and peak p for the first declared node/pipe, plus the
 * raw series for the overlay chart. Failed rows keep their error, no fake 0s.
 */
export function scenarioRowsTransient(batch) {
  return (batch || []).map((r) => {
    if (!r || !r.ok) return { ok: false, error: String((r && r.error) || '未知错误'), series: null };
    const s = r.series;
    const samples = s.samples || [];
    const last = samples[samples.length - 1] || {};
    const n0 = s.node_names[0];
    const p0 = s.pipe_names[0];
    return {
      ok: true,
      steps: s.steps,
      tEnd: last.t ?? NaN,
      pEnd: n0 !== undefined ? last.pressure?.[n0] : NaN,
      qEnd: p0 !== undefined ? last.flow?.[p0] : NaN,
      pPeak: n0 !== undefined
        ? Math.max(...samples.map((x) => x.pressure?.[n0] ?? -Infinity))
        : NaN,
      series: s,
    };
  });
}

/** Normalize a `run_case_texts` batch for the scenario table. */
export function scenarioRowsSteady(batch) {
  return (batch || []).map((r) => {
    if (!r || !r.ok) return { ok: false, error: String((r && r.error) || '未知错误') };
    const sum = r.summary || {};
    const pipes = sum.pipes || [];
    return {
      ok: true,
      converged: !!sum.converged,
      iterations: sum.iterations || 0,
      residual: sum.residual ?? NaN,
      qEnd: pipes.length ? pipes[pipes.length - 1].q : NaN,
      maxDp: pipes.length ? Math.max(...pipes.map((p) => Math.abs(p.dp || 0))) : NaN,
      summary: sum,
    };
  });
}

// ---------------------------------------------------------------------------
// Monitor variable tree (U5): layered checkbox selection over declared series
// ---------------------------------------------------------------------------

/**
 * Build the 节点 / 管段 / 瞬态帧 checkbox tree from a transient result. Only
 * quantities the kernel actually declared per frame appear (hydrate vanishes
 * without the hydrate model, alpha without a regime), so the tree is the
 * sample's own key set - never a fixed guess list.
 */
export function buildMonitorTree(t) {
  const s0 = t.samples[0] || {};
  const has = (obj, name) => Number.isFinite((s0[obj] || {})[name]);
  const groups = [];
  const nodeItems = t.nodeNames
    .map((n) => ({
      label: n,
      series: [
        ...(has('pressure', n) ? [`p:${n}`] : []),
        ...(has('alpha', n) ? [`alpha:${n}`] : []),
        ...(has('hydrate', n) ? [`hydrate:${n}`] : []),
      ],
    }))
    .filter((it) => it.series.length);
  if (nodeItems.length) groups.push({ id: 'nodes', label: '节点', items: nodeItems });
  const pipeItems = t.pipeNames
    .map((n) => ({ label: n, series: has('flow', n) ? [`Q:${n}`] : [] }))
    .filter((it) => it.series.length);
  if (pipeItems.length) groups.push({ id: 'pipes', label: '管段', items: pipeItems });
  const frameItems = [];
  if (s0.pressure && Object.keys(s0.pressure).length) {
    frameItems.push({ label: '压力包络 min/max', series: ['env:p'] });
  }
  if (s0.alpha && Object.keys(s0.alpha).length) {
    frameItems.push({ label: '全节点平均持液率', series: ['avg:alpha'] });
  }
  frameItems.push({ label: '帧游标（当前帧标记）', series: ['marker'] });
  groups.push({ id: 'frames', label: '瞬态帧', items: frameItems });
  return groups;
}

/** Per-sample values for one tree series id; null for unknown ids. */
export function monitorSeriesValues(samples, series) {
  const [kind, name] = series.split(':');
  if (kind === 'p') return samples.map((s) => s.pressure?.[name]);
  if (kind === 'alpha') return samples.map((s) => s.alpha?.[name]);
  if (kind === 'hydrate') return samples.map((s) => s.hydrate?.[name]);
  if (kind === 'Q') return samples.map((s) => s.flow?.[name]);
  if (kind === 'avg' && name === 'alpha') {
    return samples.map((s) => {
      const vals = Object.values(s.alpha || {}).filter(Number.isFinite);
      return vals.length ? vals.reduce((a, b) => a + b, 0) / vals.length : NaN;
    });
  }
  return null;
}

/** Frame-wide min/max envelope of a declared quantity (瞬态帧 band). */
export function monitorEnvBand(samples, qty) {
  const key = qty === 'p' ? 'pressure' : qty;
  return samples.map((s) => {
    const vals = Object.values(s[key] || {}).filter(Number.isFinite);
    if (!vals.length) return { lo: NaN, hi: NaN };
    return { lo: Math.min(...vals), hi: Math.max(...vals) };
  });
}

// ---------------------------------------------------------------------------
// Slug timeline (U3): multi-quantity overlay series + flow-assurance bands
// ---------------------------------------------------------------------------

/** Cell keys the two-fluid frame carries numerically (order of first sight,
 *  `index` excluded); drives the slug tab's series checkboxes. */
export function slugSeriesKeys(frames) {
  const seen = new Set();
  const keys = [];
  for (const fr of frames) {
    for (const c of fr.cells) {
      for (const [k, v] of Object.entries(c)) {
        if (k === 'index' || seen.has(k)) continue;
        if (typeof v === 'number' && Number.isFinite(v)) {
          seen.add(k);
          keys.push(k);
        }
      }
    }
  }
  // the terrain-slug output carries pressure on `node` lines instead of cells
  if (!seen.has('p') && frames.some((fr) => (fr.nodeP || []).some(Number.isFinite))) {
    keys.push('p');
  }
  return keys;
}

/** Band definitions over the kernel's own per-cell flags: a cell belongs to
 *  the band when its declared value is > 0 (s = slug slip marker, h = hydrate,
 *  wax = wax deposition, sand = sand/bed load, inh = inhibitor present). */
export const SLUG_BAND_DEFS = [
  { key: 's', label: '段塞', color: '#ff6b6b' },
  { key: 'h', label: '水合物', color: '#7db86a' },
  { key: 'wax', label: '蜡', color: '#c9a06a' },
  { key: 'sand', label: '砂', color: '#d9c26a' },
  { key: 'inh', label: '抑制剂', color: '#4bb3c4' },
];

/** Contiguous [from, to] cell-index runs per band key for one frame; only
 *  keys the frame actually carries produce entries. */
export function slugBands(frame) {
  const out = [];
  for (const def of SLUG_BAND_DEFS) {
    const has = frame.cells.some((c) => Number.isFinite(c[def.key]) && c[def.key] > 0);
    if (!has) continue;
    const runs = [];
    let run = null;
    frame.cells.forEach((c, i) => {
      const on = Number.isFinite(c[def.key]) && c[def.key] > 0;
      if (on && !run) run = { from: i, to: i };
      else if (on) run.to = i;
      else if (run) {
        runs.push(run);
        run = null;
      }
    });
    if (run) runs.push(run);
    out.push({ ...def, runs });
  }
  return out;
}

/**
 * Productivity indices for VLP sampling: a geometric grid around the declared
 * IPR's equivalent J (qAof/pRes), spanning five decades. Running the nodal
 * solver once per J (with the ipr replaced by that linear form) yields working
 * points that all lie on the wellbore's VLP curve at the fixed wellhead.
 */
export function vlpSampleJs(spec, count = 13) {
  const js = [];
  for (let k = 0; k < count; k += 1) {
    js.push(spec.jRef * 10 ** (-3 + (6 * k) / (count - 1)));
  }
  return js;
}

// ---------------------------------------------------------------------------
// U6 (§5.152): post-processing derived fields — Green-Gauss velocity gradient
// (vorticity |ω| and Q criterion) and cell-threshold isosurface patches. Both
// read only the kernel's own VTU geometry (points / cells / face rings); no
// physics is recomputed, these are CFD-Post-style derived quantities.
// ---------------------------------------------------------------------------

/**
 * Shared-face adjacency from the VTU face rings: key a face by its sorted
 * vertex ids, pair the (cell, ringIndex) occurrences. Returns
 * `{ pairs: [{a, ia, b, ib}], nbrs }` where `nbrs[c][i]` is the cell sharing
 * face i of cell c, or -1 for boundary faces. Null when faces are absent.
 */
export function faceAdjacency(cells, cellFaces) {
  if (!cellFaces || cellFaces.length !== cells.length) return null;
  const byId = new Map();
  const nbrs = cells.map((_, c) => cellFaces[c].map(() => -1));
  const pairs = [];
  cellFaces.forEach((rings, c) => {
    rings.forEach((ring, i) => {
      if (ring.length < 3) return;
      const key = [...ring].sort((x, y) => x - y).join(',');
      const seen = byId.get(key);
      if (seen === undefined) {
        byId.set(key, { c, i });
      } else {
        nbrs[seen.c][seen.i] = c;
        nbrs[c][i] = seen.c;
        pairs.push({ a: seen.c, ia: seen.i, b: c, ib: i });
        byId.delete(key);
      }
    });
  });
  return { pairs, nbrs };
}

/**
 * Green-Gauss gradient of a cell-centred vector field: G_c[i][j] = dU_j/dx_i
 * = (1/V_c) Σ_f A_f ⊗ U_f, with the face value U_f averaged over the two
 * adjacent cells (own value on boundary faces). Face area vectors are
 * re-oriented outward from each cell by the centroid test, so the writer's
 * ring order never matters. Exact for linear fields on symmetric stencils.
 */
export function greenGaussVelocityGradient(points, cells, cellFaces, cellU) {
  const volumes = cellVolumes(points, cells, cellFaces);
  const adj = faceAdjacency(cells, cellFaces);
  if (!volumes || !adj) return null;
  return cells.map((vids, c) => {
    const centre = cellCentroid(points, vids);
    const g = [0, 0, 0, 0, 0, 0, 0, 0, 0]; // row-major dU_j/dx_i
    cellFaces[c].forEach((ring, i) => {
      if (ring.length < 3) return;
      const s = areaVector(points, ring);
      const f = cellCentroid(points, ring);
      let ax = s.x, ay = s.y, az = s.z;
      if (ax * (f.x - centre.x) + ay * (f.y - centre.y) + az * (f.z - centre.z) < 0) {
        ax = -ax; ay = -ay; az = -az; // force outward orientation
      }
      const n = adj.nbrs[c][i];
      const own = cellU[c] || [0, 0, 0];
      const nb = n >= 0 ? (cellU[n] || [0, 0, 0]) : own;
      const uf = [(own[0] + nb[0]) / 2, (own[1] + nb[1]) / 2, (own[2] + nb[2]) / 2];
      const a = [ax, ay, az];
      for (let p = 0; p < 3; p += 1)
        for (let q = 0; q < 3; q += 1)
          g[3 * p + q] += a[p] * uf[q];
    });
    const v = volumes[c] || 1;
    return g.map((x) => x / v);
  });
}

/**
 * Vorticity magnitude and Q criterion per cell from the Green-Gauss tensor.
 * ω = curl U (g[i][j] = dU_j/dx_i), Q = (|Ω|² − |S|²)/2 with S/Ω the
 * symmetric / antisymmetric parts. Returns null when the VTU carries no
 * faces array (nothing to integrate over).
 */
export function vorticityAndQ(points, cells, cellFaces, cellU) {
  const grads = greenGaussVelocityGradient(points, cells, cellFaces, cellU);
  if (!grads) return null;
  const omega = [];
  const mag = [];
  const q = [];
  for (const g of grads) {
    const wx = g[3 * 1 + 2] - g[3 * 2 + 1];
    const wy = g[3 * 2 + 0] - g[3 * 0 + 2];
    const wz = g[3 * 0 + 1] - g[3 * 1 + 0];
    omega.push([wx, wy, wz]);
    mag.push(Math.hypot(wx, wy, wz));
    let om2 = 0;
    let s2 = 0;
    for (let i = 0; i < 3; i += 1) {
      for (let j = 0; j < 3; j += 1) {
        const antisym = (g[3 * i + j] - g[3 * j + i]) / 2;
        const sym = (g[3 * i + j] + g[3 * j + i]) / 2;
        om2 += antisym * antisym;
        s2 += sym * sym;
      }
    }
    q.push((om2 - s2) / 2);
  }
  return { omega, mag, q };
}

/**
 * Isosurface patches of a cell-centred scalar at `threshold`: every internal
 * face whose two adjacent values straddle the level becomes a quad taken from
 * the face ring, shifted along the cell-centre line to the linear crossing
 * fraction f = (t − a)/(b − a). Returns { patches: [ring-of-xyz], f: [..],
 * src: [[a, b], ..] } (src lets the view drop patches killed by clipping).
 */
export function isoSurfacePatches(points, cells, cellFaces, values, threshold) {
  const adj = faceAdjacency(cells, cellFaces);
  const out = { patches: [], f: [], src: [] };
  if (!adj) return out;
  for (const { a, ia, b } of adj.pairs) {
    const va = values[a];
    const vb = values[b];
    if (!Number.isFinite(va) || !Number.isFinite(vb)) continue;
    if ((va - threshold) * (vb - threshold) >= 0) continue;
    if (Math.abs(vb - va) < 1e-300) continue;
    const f = (threshold - va) / (vb - va);
    const ca = cellCentroid(points, cells[a]);
    const cb = cellCentroid(points, cells[b]);
    const dx = (f - 0.5) * (cb.x - ca.x);
    const dy = (f - 0.5) * (cb.y - ca.y);
    const dz = (f - 0.5) * (cb.z - ca.z);
    out.patches.push(cellFaces[a][ia].map((v) => [
      points[v].x + dx, points[v].y + dy, points[v].z + dz,
    ]));
    out.f.push(f);
    out.src.push([a, b]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// U9 (§5.153): CFD-Post-style export — every-cell field CSV and a planar
// profile CSV. Both read only the kernel's VTU payload; the plane selection
// mirrors evaluateReports' face-centre projection rule exactly.
// ---------------------------------------------------------------------------

/** CSV of all cell-centred fields: centroid, p, Ux/Uy/Uz, |U|. */
export function fieldsToCsv(points, cells, cellU, cellP) {
  const rows = cells.map((vids, c) => {
    const m = cellCentroid(points, vids);
    const u = cellU[c] || [0, 0, 0];
    return [c, m.x, m.y, m.z,
      cellP[c] === undefined ? '' : cellP[c],
      u[0], u[1], u[2], Math.hypot(u[0], u[1], u[2])].join(',');
  });
  return ['cell,cx,cy,cz,p,Ux,Uy,Uz,magU', ...rows].join('\n') + '\n';
}

/**
 * CSV of the cells cut by the plane `normal @ frac` (frac spans the face
 * centre projections, so 0/1 select the true boundary layers). Null when the
 * VTU carries no faces — the plane cannot be located without them.
 */
export function planeProfileCsv(points, cells, cellFaces, cellU, cellP, normal, frac) {
  const volumes = cellVolumes(points, cells, cellFaces);
  if (!volumes) return null;
  const len = Math.hypot(normal[0], normal[1], normal[2]) || 1;
  const nx = normal[0] / len, ny = normal[1] / len, nz = normal[2] / len;
  let lo = Infinity, hi = -Infinity;
  for (const rings of cellFaces) for (const ring of rings) {
    if (ring.length < 3) continue;
    const f = cellCentroid(points, ring);
    const p = f.x * nx + f.y * ny + f.z * nz;
    lo = Math.min(lo, p); hi = Math.max(hi, p);
  }
  const target = lo + (frac ?? 0.5) * (hi - lo || 1);
  const rows = [];
  cells.forEach((vids, c) => {
    const onPlane = cellFaces[c].some((ring) => {
      if (ring.length < 3) return false;
      const s = areaVector(points, ring);
      if (Math.abs(s.x * nx + s.y * ny + s.z * nz) < 1e-14) return false;
      const f = cellCentroid(points, ring);
      return Math.abs(f.x * nx + f.y * ny + f.z * nz - target)
        <= Math.cbrt(volumes[c]) * 0.3;
    });
    if (!onPlane) return;
    const m = cellCentroid(points, vids);
    const u = cellU[c] || [0, 0, 0];
    rows.push([m.x, m.y, m.z,
      cellP[c] === undefined ? '' : cellP[c],
      u[0], u[1], u[2], Math.hypot(u[0], u[1], u[2])].join(','));
  });
  return ['cx,cy,cz,p,Ux,Uy,Uz,magU', ...rows].join('\n') + '\n';
}
