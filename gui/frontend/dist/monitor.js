/**
 * F17 online monitor: incremental parser for the kernel's agreed per-
 * iteration stream lines
 *   ITER <n> <residU> <residV> <residW> <residP> ...
 *   MONITOR <name> <iter> <value>
 * fed as arbitrary chunks of a streaming stdout buffer. Half lines stay in
 * the buffer until their newline arrives; out-of-order iteration numbers are
 * rejected (never merged). The log-scale mapping and rolling window for the
 * live curves are pure too, so node --test covers everything. No DOM.
 */

/** Fresh stream state. `iters` are {n, resid[]}, `monitors` name -> [{iter, value}]. */
export function createMonitorState(maxPoints = 2000) {
  return {
    buffer: '',
    lastIter: 0,
    iters: [],
    monitors: {},
    rejected: 0,
    malformed: 0,
    maxPoints,
  };
}

function pushCapped(list, item, max) {
  list.push(item);
  if (list.length > max) list.splice(0, list.length - max);
}

/**
 * Feed one streamed chunk. Returns { state, events } where events lists
 * { kind: 'iter' | 'monitor' | 'rejected' | 'malformed', ... } for THIS
 * call; state is the same (mutated) object passed in. Lines may arrive split
 * anywhere: only complete lines (terminated by \n) are ever parsed, the
 * trailing fragment is kept verbatim.
 */
export function feedMonitor(state, chunk) {
  const events = [];
  state.buffer += String(chunk ?? '');
  const nl = state.buffer.lastIndexOf('\n');
  if (nl < 0) return { state, events };
  const lines = state.buffer.slice(0, nl).split('\n');
  state.buffer = state.buffer.slice(nl + 1);
  for (const raw of lines) {
    const line = raw.trim().replace('[chflow] ', '').trim();
    if (!line) continue;
    const parts = line.split(/\s+/);
    if (parts[0] === 'ITER') {
      const n = Number(parts[1]);
      const resid = parts.slice(2).map(Number);
      if (!Number.isInteger(n) || n < 0 || resid.length === 0 || resid.some((r) => !Number.isFinite(r))) {
        state.malformed += 1;
        events.push({ kind: 'malformed', line });
        continue;
      }
      if (n <= state.lastIter) {
        state.rejected += 1;
        events.push({ kind: 'rejected', line, reason: 'iter<=lastIter' });
        continue;
      }
      state.lastIter = n;
      pushCapped(state.iters, { n, resid }, state.maxPoints);
      events.push({ kind: 'iter', n, resid });
    } else if (parts[0] === 'MONITOR') {
      const name = parts[1];
      const iter = Number(parts[2]);
      const value = Number(parts[3]);
      if (!name || !Number.isInteger(iter) || iter < 0 || !Number.isFinite(value)) {
        state.malformed += 1;
        events.push({ kind: 'malformed', line });
        continue;
      }
      // per-name monotonicity: a sample older than the last accepted iter
      // for this monitor is out of order and rejected
      if (!state.monitors[name]) state.monitors[name] = [];
      const list = state.monitors[name];
      const last = list.length ? list[list.length - 1].iter : -1;
      if (iter < last) {
        state.rejected += 1;
        events.push({ kind: 'rejected', line, reason: 'monitor iter<last' });
        continue;
      }
      pushCapped(list, { iter, value }, state.maxPoints);
      events.push({ kind: 'monitor', name, iter, value });
    } else {
      // unknown lines are ignored silently: the log panel already shows them
    }
  }
  return { state, events };
}

/**
 * log10 vertical mapping for the residual panel: value v in [lo, hi]
 * (both > 0) -> y in [0, height] with y = 0 at the TOP (hi) — the usual
 * screen convention. Values outside clamp. Non-positive values map to the
 * bottom edge so a zero residual still draws.
 */
export function logY(v, lo, hi, height) {
  const l = Math.log10(Math.max(1e-300, lo));
  const h = Math.log10(Math.max(1e-299, hi));
  if (!(h > l)) return 0;
  if (!(v > 0)) return height;
  const t = (Math.log10(v) - l) / (h - l);
  return height - Math.max(0, Math.min(1, t)) * height;
}

/** Decade ticks (10^k) inside [lo, hi], ascending. */
export function logTicks(lo, hi) {
  const out = [];
  if (!(lo > 0) || !(hi >= lo)) return out;
  for (let k = Math.ceil(Math.log10(lo)); k <= Math.floor(Math.log10(hi) + 1e-9); k += 1) {
    out.push(10 ** k);
  }
  return out;
}

/** Rolling window: the last maxN points of a {t, v} series. */
export function rollingWindow(points, maxN) {
  const n = Math.max(1, Math.round(maxN));
  return (points || []).slice(Math.max(0, (points || []).length - n));
}

/**
 * Draw commands for the live panel: one polyline per named series (residual
 * components resid:0..k and monitor:NAME), log-scaled y, iteration x.
 * `state` is a monitor state; returns { commands, seriesNames, xLo, xHi }.
 */
export function monitorDrawCommands(state, width, height) {
  const commands = [];
  const xs = state.iters.map((it) => it.n);
  const xLo = xs.length ? xs[0] : 0;
  const xHi = xs.length ? xs[xs.length - 1] : 1;
  const x = (n) => ((n - xLo) / (xHi - xLo || 1)) * width;
  let lo = Infinity;
  let hi = -Infinity;
  for (const it of state.iters) {
    for (const r of it.resid) {
      if (r > 0) { lo = Math.min(lo, r); hi = Math.max(hi, r); }
    }
  }
  if (!(hi > lo)) { lo = 1e-8; hi = 1; }
  lo /= 10;
  hi *= 10;
  for (const t of logTicks(lo, hi)) {
    commands.push({
      type: 'line', x1: 0, y1: logY(t, lo, hi, height), x2: width, y2: logY(t, lo, hi, height),
      role: 'tick', value: t,
    });
  }
  const nComp = state.iters.reduce((m, it) => Math.max(m, it.resid.length), 0);
  const seriesNames = [];
  for (let c = 0; c < nComp; c += 1) {
    const pts = state.iters
      .filter((it) => Number.isFinite(it.resid[c]))
      .map((it) => `${x(it.n).toFixed(1)},${logY(it.resid[c], lo, hi, height).toFixed(1)}`)
      .join(' ');
    if (pts) {
      seriesNames.push(`resid:${c}`);
      commands.push({ type: 'polyline', points: pts, role: 'resid', series: c });
    }
  }
  let mLo = Infinity;
  let mHi = -Infinity;
  for (const list of Object.values(state.monitors)) {
    for (const p of list) {
      mLo = Math.min(mLo, p.value);
      mHi = Math.max(mHi, p.value);
    }
  }
  if (Number.isFinite(mLo) && Number.isFinite(mHi)) {
    const span = mHi - mLo || Math.abs(mHi) || 1;
    for (const [name, list] of Object.entries(state.monitors)) {
      const pts = list
        .map((p) => `${x(p.iter).toFixed(1)},${(height - ((p.value - mLo) / span) * height).toFixed(1)}`)
        .join(' ');
      if (pts) {
        seriesNames.push(`monitor:${name}`);
        commands.push({ type: 'polyline', points: pts, role: 'monitor', name });
      }
    }
  }
  return { commands, seriesNames, xLo, xHi };
}
