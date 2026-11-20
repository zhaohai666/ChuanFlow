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
      // inline equipment (network.cpp pipe keys), 0 = absent
      valve: kv.valve !== undefined ? Number(kv.valve) : 0,
      pumpHead: kv.pumpHead !== undefined ? Number(kv.pumpHead) : 0,
      pumpRate: kv.pumpRate !== undefined ? Number(kv.pumpRate) : 0,
      pumpShutoffHead: kv.pumpShutoffHead !== undefined ? Number(kv.pumpShutoffHead) : 0,
      pumpStages: kv.pumpStages !== undefined ? Number(kv.pumpStages) : 0,
      pumpSpeed: kv.pumpSpeed !== undefined ? Number(kv.pumpSpeed) : 0,
      pumpGasFraction: kv.pumpGasFraction !== undefined ? Number(kv.pumpGasFraction) : 0,
      pumpIntakePressure: kv.pumpIntakePressure !== undefined ? Number(kv.pumpIntakePressure) : 0,
      liftDomePressure: kv.liftDomePressure !== undefined ? Number(kv.liftDomePressure) : 0,
      liftAreaRatio: kv.liftAreaRatio !== undefined ? Number(kv.liftAreaRatio) : 0,
      liftSpringPressure: kv.liftSpringPressure !== undefined ? Number(kv.liftSpringPressure) : 0,
      liftDischargeCoeff: kv.liftDischargeCoeff !== undefined ? Number(kv.liftDischargeCoeff) : 0,
      liftKappa: kv.liftKappa !== undefined ? Number(kv.liftKappa) : 0,
      rodPlungerDiameter: kv.rodPlungerDiameter !== undefined ? Number(kv.rodPlungerDiameter) : 0,
      rodStrokeLength: kv.rodStrokeLength !== undefined ? Number(kv.rodStrokeLength) : 0,
      rodStrokesPerMinute: kv.rodStrokesPerMinute !== undefined ? Number(kv.rodStrokesPerMinute) : 0,
      rodDiameter: kv.rodDiameter !== undefined ? Number(kv.rodDiameter) : 0,
      rodLength: kv.rodLength !== undefined ? Number(kv.rodLength) : 0,
      rodYoungModulus: kv.rodYoungModulus !== undefined ? Number(kv.rodYoungModulus) : 0,
      rodDensity: kv.rodDensity !== undefined ? Number(kv.rodDensity) : 0,
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
        cells: [],
      };
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
