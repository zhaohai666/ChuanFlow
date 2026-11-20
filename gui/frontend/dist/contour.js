/**
 * F16 contour extraction pure helpers: classic marching squares on a regular
 * 2-D scalar grid with linear-interpolated crossings and a fixed saddle rule,
 * segment chaining, a slice-plane rasterizer that feeds the 3-D field view,
 * and 1-D profile bands that reuse the same level list. No DOM, no Tauri.
 */

import { rainbowColor } from './workspace.js';

/**
 * Marching squares over `field` (row-major, index j*nx + i, node (i, j) at
 * (x0 + i*dx, y0 + j*dy)) for every threshold in `levels`.
 * Returns [{ level, segments: [[[x1, y1], [x2, y2]], ...] }].
 * Saddle rule (cells 5 and 10, the two diagonal cases): the cell centre is
 * the mean of the four corners; centre >= level connects the branches that
 * isolate the two LOW corners, centre < level the branches that isolate the
 * two HIGH corners. The rule is deterministic and applied identically to
 * both saddle cases.
 */
export function marchingSquares(field, nx, ny, levels, opts = {}) {
  const dx = opts.dx ?? 1;
  const dy = opts.dy ?? 1;
  const x0 = opts.x0 ?? 0;
  const y0 = opts.y0 ?? 0;
  const out = [];
  for (const level of levels) {
    const segments = [];
    for (let j = 0; j + 1 < ny; j += 1) {
      for (let i = 0; i + 1 < nx; i += 1) {
        const va = field[j * nx + i];           // bottom-left
        const vb = field[j * nx + i + 1];       // bottom-right
        const vc = field[(j + 1) * nx + i + 1]; // top-right
        const vd = field[(j + 1) * nx + i];     // top-left
        if (![va, vb, vc, vd].every(Number.isFinite)) continue;
        let cse = 0;
        if (va >= level) cse |= 1;
        if (vb >= level) cse |= 2;
        if (vc >= level) cse |= 4;
        if (vd >= level) cse |= 8;
        if (cse === 0 || cse === 15) continue;
        const t = (p, q) => (Math.abs(q - p) < 1e-300 ? 0 : (level - p) / (q - p));
        const bottom = [x0 + (i + t(va, vb)) * dx, y0 + j * dy];
        const right = [x0 + (i + 1) * dx, y0 + (j + t(vb, vc)) * dy];
        const top = [x0 + (i + t(vd, vc)) * dx, y0 + (j + 1) * dy];
        const left = [x0 + i * dx, y0 + (j + t(va, vd)) * dy];
        const push = (p, q) => segments.push([p, q]);
        switch (cse) {
          case 1: case 14: push(left, bottom); break;
          case 2: case 13: push(bottom, right); break;
          case 3: case 12: push(left, right); break;
          case 4: case 11: push(right, top); break;
          case 6: case 9: push(bottom, top); break;
          case 7: case 8: push(left, top); break;
          case 5: {
            const centre = (va + vb + vc + vd) / 4;
            if (centre < level) { push(left, bottom); push(right, top); } // 隔离高角 a、c
            else { push(left, top); push(bottom, right); } // 隔离低角 b、d
            break;
          }
          case 10: {
            const centre = (va + vb + vc + vd) / 4;
            if (centre < level) { push(bottom, right); push(left, top); } // 隔离高角 b、d
            else { push(left, bottom); push(right, top); } // 隔离低角 a、c
            break;
          }
          default: break;
        }
      }
    }
    out.push({ level, segments });
  }
  return out;
}

/**
 * Chain segments into polylines by matching endpoints (tolerance eps).
 * Returns arrays of point lists; open chains and closed loops are both
 * emitted as-is (a loop simply ends where it started).
 */
export function joinSegments(segments, eps = 1e-9) {
  const key = (p) => `${Math.round(p[0] / eps)},${Math.round(p[1] / eps)}`;
  const used = segments.map(() => false);
  const chains = [];
  const findNext = (p, from) => {
    const k = key(p);
    for (let i = from; i < segments.length; i += 1) {
      if (used[i]) continue;
      if (key(segments[i][0]) === k) return { i, at: 0 };
      if (key(segments[i][1]) === k) return { i, at: 1 };
    }
    return null;
  };
  for (let s = 0; s < segments.length; s += 1) {
    if (used[s]) continue;
    used[s] = true;
    const chain = [segments[s][0], segments[s][1]];
    for (;;) {
      const fwd = findNext(chain[chain.length - 1], 0);
      if (fwd) {
        used[fwd.i] = true;
        chain.push(segments[fwd.i][fwd.at === 0 ? 1 : 0]);
        continue;
      }
      const bwd = findNext(chain[0], 0);
      if (!bwd) break;
      used[bwd.i] = true;
      chain.unshift(segments[bwd.i][bwd.at === 0 ? 1 : 0]);
    }
    chains.push(chain);
  }
  return chains;
}

/**
 * Rasterize a cell-centred 3-D field onto the plane through the KEPT cells
 * (pass the clip's `keep` mask so the grid lives on the cut face; `null`
 * keeps everything) with unit-ish `normal`, into a res x res
 * marching-squares grid. Bins take the nearest kept centroid's value.
 * Returns { field, nx, ny, x0, y0, dx, dy, u, v } with u/v the in-plane
 * basis vectors, so a renderer maps grid coords back to world points:
 * world = (x, y) expressed in the basis, i.e. p = u*x + v*y (+ any fixed
 * normal offset, which an orthographic projection drops anyway).
 */
export function rasterizeSliceField(points, cells, values, keep, normal, res = 32) {
  const len = Math.hypot(normal[0], normal[1], normal[2]);
  if (!Number.isFinite(len) || len < 1e-12) return null;
  const n = [normal[0] / len, normal[1] / len, normal[2] / len];
  const ref = Math.abs(n[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  let u = [n[1] * ref[2] - n[2] * ref[1], n[2] * ref[0] - n[0] * ref[2], n[0] * ref[1] - n[1] * ref[0]];
  const ul = Math.hypot(...u) || 1;
  u = u.map((x) => x / ul);
  const v = [n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2], n[0] * u[1] - n[1] * u[0]];
  const items = [];
  cells.forEach((vids, c) => {
    if (keep && !keep[c]) return;
    let cx = 0, cy = 0, cz = 0;
    for (const vid of vids) { cx += points[vid].x; cy += points[vid].y; cz += points[vid].z; }
    const m = vids.length || 1;
    const p = [cx / m, cy / m, cz / m];
    items.push({
      p,
      pu: u[0] * p[0] + u[1] * p[1] + u[2] * p[2],
      pv: v[0] * p[0] + v[1] * p[1] + v[2] * p[2],
      value: values[c],
    });
  });
  if (!items.length) return null;
  let uLo = Infinity, uHi = -Infinity, vLo = Infinity, vHi = -Infinity;
  let mu = 0, mv = 0;
  const mean = [0, 0, 0];
  for (const it of items) {
    uLo = Math.min(uLo, it.pu); uHi = Math.max(uHi, it.pu);
    vLo = Math.min(vLo, it.pv); vHi = Math.max(vHi, it.pv);
    mu += it.pu; mv += it.pv;
    mean[0] += it.p[0]; mean[1] += it.p[1]; mean[2] += it.p[2];
  }
  mu /= items.length;
  mv /= items.length;
  mean[0] /= items.length;
  mean[1] /= items.length;
  mean[2] /= items.length;
  // a point ON the mean plane: world(x, y) = origin + u*x + v*y
  const origin = [
    mean[0] - u[0] * mu - v[0] * mv,
    mean[1] - u[1] * mu - v[1] * mv,
    mean[2] - u[2] * mu - v[2] * mv,
  ];
  const R = Math.max(2, Math.round(res));
  const du = (uHi - uLo || 1) / (R - 1);
  const dv = (vHi - vLo || 1) / (R - 1);
  const field = new Array(R * R).fill(NaN);
  for (let j = 0; j < R; j += 1) {
    for (let i = 0; i < R; i += 1) {
      const gu = uLo + i * du;
      const gv = vLo + j * dv;
      let best = null;
      let bestD = Infinity;
      for (const it of items) {
        const dd = (it.pu - gu) ** 2 + (it.pv - gv) ** 2;
        if (dd < bestD) { bestD = dd; best = it; }
      }
      if (best && Number.isFinite(best.value)) field[j * R + i] = best.value;
    }
  }
  return { field, nx: R, ny: R, x0: uLo, y0: vLo, dx: du, dy: dv, origin, u, v };
}

/** Colour for one contour level on the field's [lo, hi] rainbow scale. */
export function contourColor(level, lo, hi) {
  return rainbowColor(level, lo, hi);
}

/**
 * 1-D profile bands (沿程剖面等值带): for each level, the contiguous
 * [from, to] index runs where values[i] >= level. Non-finite entries count
 * as below every level. Reuses the same level list as the 2-D contours.
 */
export function profileBands(values, levels) {
  return levels.map((level) => {
    const runs = [];
    let run = null;
    (values || []).forEach((v, i) => {
      const on = Number.isFinite(v) && v >= level;
      if (on && !run) run = { from: i, to: i };
      else if (on) run.to = i;
      else if (run) { runs.push(run); run = null; }
    });
    if (run) runs.push(run);
    return { level, runs };
  });
}
