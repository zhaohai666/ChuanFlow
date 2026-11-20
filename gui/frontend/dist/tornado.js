/**
 * P20 sensitivity tornado chart pure helpers, layered on the existing batch
 * sweep scheduling (steady run_case_texts / transient run_transient_texts).
 * The caller supplies one base metric value plus, per parameter, the metric
 * evaluated at the low and high level; everything here is sorting, scaling
 * and draw-command generation. No DOM, no Tauri.
 */

/**
 * Data model. `cases` = [{ label, low, high }] with the three numbers being
 * metric values (low/base already solved, `base` passed separately).
 * Returns rows sorted by swing = |high - low| descending:
 *   { label, low, high, base, lowDelta: low - base, highDelta: high - base,
 *     swing }
 * Non-finite levels collapse to the base value (zero-width bar) instead of
 * poisoning the sort, and ties keep input order (stable).
 */
export function tornadoData(base, cases) {
  const b = Number.isFinite(base) ? base : 0;
  const rows = (cases || []).map((c, i) => {
    const low = Number.isFinite(c.low) ? c.low : b;
    const high = Number.isFinite(c.high) ? c.high : b;
    return {
      label: String(c.label ?? `#${i + 1}`),
      low,
      high,
      base: b,
      lowDelta: low - b,
      highDelta: high - b,
      swing: Math.abs(high - low),
    };
  });
  return rows
    .map((r, i) => ({ ...r, _i: i }))
    .sort((x, y) => (y.swing - x.swing) || (x._i - y._i))
    .map(({ _i, ...r }) => r);
}

/**
 * SVG/Canvas draw commands for the rows. The base value maps exactly to
 * `baseX` (the centred vertical baseline); every row contributes a low bar
 * (baseX -> x(low), always leftward-or-zero width) and a high bar (baseX ->
 * x(high)). Degenerate rows (low == high == base) produce zero-width rects.
 * Returns { commands, baseX, domain: [lo, hi], rowHeight }.
 */
export function tornadoCommands(rows, opts = {}) {
  const { width = 560, height = 220, marginL = 130, marginR = 20, marginTop = 10 } = opts;
  const plotW = width - marginL - marginR;
  let lo = Number.isFinite(opts.domainLo) ? opts.domainLo : Infinity;
  let hi = Number.isFinite(opts.domainHi) ? opts.domainHi : -Infinity;
  for (const r of rows) {
    lo = Math.min(lo, r.low, r.base);
    hi = Math.max(hi, r.high, r.base);
  }
  if (!(hi > lo)) { lo = (Number.isFinite(lo) ? lo : 0) - 1; hi = lo + 2; }
  const base = rows.length ? rows[0].base : 0;
  const x = (v) => marginL + ((v - lo) / (hi - lo)) * plotW;
  const baseX = x(base);
  const rowH = rows.length ? (height - 2 * marginTop) / rows.length : 0;
  const commands = [];
  commands.push({
    type: 'line', x1: baseX, y1: marginTop, x2: baseX, y2: marginTop + rowH * rows.length,
    role: 'baseline',
  });
  rows.forEach((r, i) => {
    const yc = marginTop + i * rowH;
    const h = Math.max(2, rowH * 0.55);
    const y = yc + (rowH - h) / 2;
    const xl = x(r.low);
    const xh = x(r.high);
    commands.push({
      type: 'rect', x: Math.min(xl, baseX), y, w: Math.abs(baseX - xl), h,
      color: 'var(--accent)', role: 'low', label: r.label,
    });
    commands.push({
      type: 'rect', x: Math.min(xh, baseX), y, w: Math.abs(xh - baseX), h,
      color: '#e0a33e', role: 'high', label: r.label,
    });
    commands.push({
      type: 'text', x: marginL - 6, y: y + h / 2, anchor: 'end',
      text: `${r.label} (±${r.swing.toPrecision(3)})`, role: 'name',
    });
  });
  return { commands, baseX, domain: [lo, hi], rowHeight: rowH };
}
