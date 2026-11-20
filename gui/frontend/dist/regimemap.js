/**
 * P4 flow-pattern map (流型图) pure helpers. The regime dispatch is a direct
 * port of the kernel's single source of truth `beggsBrillPattern`
 * (src/1d/src/beggsBrill.cpp, 5.85/D4): the four published boundaries
 *   L1 = 316 lambda^0.302,   L2 = 0.0009252 lambda^-2.4684,
 *   L3 = 0.1 lambda^-1.4516, L4 = 0.5 lambda^-6.738,
 * with the published precedence (last matching rule wins). No DOM, no Tauri:
 * node --test covers the geometry, placement and hit logic.
 */

export const REGIME_NAMES = ['segregated 分层', 'transition 过渡', 'intermittent 间歇', 'distributed 分散'];

/** One fill colour per regime index (0..3), matching the GUI palette family. */
export const REGIME_COLORS = ['#2b6cb0', '#4bb3c4', '#e0a33e', '#c96a9b'];

export const GRAVITY = 9.81;

/**
 * Port of beggsBrillPattern(lambda, nfr): boundaries + regime. The kernel's
 * precedence is reproduced exactly (regime 0 segregated, 1 transition,
 * 2 intermittent, 3 distributed; the last matching rule wins).
 */
export function beggsBrillPattern(lambda, nfr) {
  const pat = {
    l1: 316.0 * Math.pow(lambda, 0.302),
    l2: 0.0009252 * Math.pow(lambda, -2.4684),
    l3: 0.1 * Math.pow(lambda, -1.4516),
    l4: 0.5 * Math.pow(lambda, -6.738),
    regime: 3,
    transitionWeight: 0,
  };
  if ((lambda < 0.01 && nfr < pat.l1) || (lambda >= 0.01 && nfr < pat.l2)) {
    pat.regime = 0;
  } else if (lambda >= 0.01 && nfr >= pat.l2 && nfr <= pat.l3) {
    pat.regime = 1;
  } else if ((0.01 <= lambda && lambda < 0.4 && nfr > pat.l3 && nfr <= pat.l1)
    || (lambda >= 0.4 && nfr > pat.l3 && nfr <= pat.l4)) {
    pat.regime = 2;
  } else {
    pat.regime = 3;
  }
  if (pat.regime === 1) pat.transitionWeight = (pat.l3 - nfr) / (pat.l3 - pat.l2);
  return pat;
}

/**
 * Map superficial velocities (m/s) + pipe diameter (m) onto map coordinates:
 * lambda = vsl/vm (no-slip liquid content), Nfr = vm^2/(g D). vm <= 0 or a
 * non-positive diameter yields lambda = 0 / nfr = 0 (the map's degenerate
 * no-flow corner), never NaN.
 */
export function mapPoint(vsl, vsg, diameter, g = GRAVITY) {
  const a = Number(vsl);
  const b = Number(vsg);
  const d = Number(diameter);
  const vm = (Number.isFinite(a) ? a : 0) + (Number.isFinite(b) ? b : 0);
  const lambda = vm > 1e-12 ? (Number.isFinite(a) ? a : 0) / vm : 0;
  const nfr = (vm > 0 && Number.isFinite(d) && d > 0 && g > 0) ? (vm * vm) / (g * d) : 0;
  return { lambda, nfr, vm };
}

/** Region hit test: the regime index at one map coordinate. */
export function regimeHit(lambda, nfr) {
  return beggsBrillPattern(lambda, nfr).regime;
}

/**
 * Sample the four boundary curves L1..L4 over lambda in [lamLo, lamHi]
 * (count+1 points each, ascending lambda). Non-finite / negative / infinite
 * samples are dropped (lambda = 0 makes L2..L4 blow up), so every returned
 * point is a drawable finite coordinate.
 */
export function sampleBoundaries(lamLo, lamHi, count = 64) {
  const out = { l1: [], l2: [], l3: [], l4: [] };
  const n = Math.max(1, Math.round(count));
  for (let i = 0; i <= n; i += 1) {
    const lam = lamLo + ((lamHi - lamLo) * i) / n;
    const pat = beggsBrillPattern(lam, 1);
    for (const key of ['l1', 'l2', 'l3', 'l4']) {
      const v = pat[key];
      if (Number.isFinite(v) && v > 0) out[key].push({ lambda: lam, nfr: v });
    }
  }
  return out;
}

/**
 * Raster regime field for the colour bands: one entry per (i, j) of the
 * axis grids; field[j][i] = regimeHit(lambdaAxis[i], nfrAxis[j]).
 */
export function regimeField(lambdaAxis, nfrAxis) {
  return nfrAxis.map((nfr) => lambdaAxis.map((lam) => regimeHit(lam, nfr)));
}

/**
 * Data -> canvas placement. `axes` = { lamLo, lamHi, nfrLo, nfrHi } with the
 * nfr axis on log10 (the published boundaries span ~13 decades). Returns
 * { x, y } in [0, width] x [0, height] with y growing downward; points
 * outside the range are clamped to the border so markers stay visible.
 */
export function projectToCanvas(lambda, nfr, axes, width, height) {
  const x = width * ((lambda - axes.lamLo) / (axes.lamHi - axes.lamLo || 1));
  const lLo = Math.log10(Math.max(axes.nfrLo, 1e-12));
  const lHi = Math.log10(Math.max(axes.nfrHi, 1e-11));
  const lv = Math.log10(Math.max(nfr, 1e-12));
  const y = height - height * ((lv - lLo) / (lHi - lLo || 1));
  return {
    x: Math.max(0, Math.min(width, x)),
    y: Math.max(0, Math.min(height, y)),
  };
}

/**
 * Draw-command list for the whole map: colour-band rects (raster of
 * regimeField), boundary polylines, axis labels and the overlay markers.
 * The renderer (app.js) only executes primitives; every decision (which
 * regime, where, which colour) is made here and is unit-testable.
 */
export function regimeMapDrawCommands(opts) {
  const {
    width = 640, height = 420,
    lamLo = 0, lamHi = 1, nfrLo = 1e-3, nfrHi = 1e8,
    resX = 48, resY = 40,
    points = [],
  } = opts || {};
  const axes = { lamLo, lamHi, nfrLo, nfrHi };
  const lambdaAxis = Array.from({ length: resX }, (_, i) => lamLo + ((lamHi - lamLo) * (i + 0.5)) / resX);
  const nfrAxis = Array.from({ length: resY }, (_, j) =>
    10 ** (Math.log10(Math.max(nfrLo, 1e-12))
      + ((Math.log10(Math.max(nfrHi, 1e-11)) - Math.log10(Math.max(nfrLo, 1e-12))) * (j + 0.5)) / resY));
  const field = regimeField(lambdaAxis, nfrAxis);
  const commands = [];
  const bw = width / resX;
  const bh = height / resY;
  for (let j = 0; j < resY; j += 1) {
    for (let i = 0; i < resX; i += 1) {
      commands.push({
        type: 'rect',
        x: i * bw,
        y: height - (j + 1) * bh,
        w: bw,
        h: bh,
        color: REGIME_COLORS[field[j][i]],
        regime: field[j][i],
      });
    }
  }
  const bounds = sampleBoundaries(lamLo, lamHi, 96);
  for (const key of ['l1', 'l2', 'l3', 'l4']) {
    const pts = bounds[key].map((p) => projectToCanvas(p.lambda, p.nfr, axes, width, height));
    if (pts.length > 1) commands.push({ type: 'polyline', points: pts, boundary: key });
  }
  // one label per regime, anchored at a coordinate guaranteed to sit inside
  // that region for the default (0..1, 1e-3..1e8) map window
  const anchors = [[0.5, 1e-3], [0.5, 0.05], [0.5, 10], [0.5, 1e6]];
  for (let r = 0; r < 4; r += 1) {
    const a = projectToCanvas(anchors[r][0], anchors[r][1], axes, width, height);
    commands.push({
      type: 'label',
      x: a.x,
      y: a.y,
      text: REGIME_NAMES[r],
      color: REGIME_COLORS[r],
      regime: r,
    });
  }
  for (const p of points) {
    const q = projectToCanvas(p.lambda, p.nfr, axes, width, height);
    commands.push({
      type: 'marker',
      x: q.x,
      y: q.y,
      regime: regimeHit(p.lambda, p.nfr),
      label: p.label || '',
    });
  }
  return { commands, axes };
}
