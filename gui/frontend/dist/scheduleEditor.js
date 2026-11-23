/**
 * E3 — structured editing of the 1D transient multi-segment event tables
 * (`valveSchedule{}` inside a pipes entry, top-level `sourceSchedule{}`).
 * All keys and guards mirror src/1d/src/transient.cpp (§5.170): the value key
 * is `k` for valve events and `rate` for source events, each event is an
 * anonymous block `{ t <t>; <key> <v>; }`, tables are piecewise linear,
 * EQUAL ADJACENT TIMES ARE A STEP, values hold before the first / after the
 * last event, and the march samples the table at t + dt. Validation mirrors
 * the loader's thirteen cfFatal rejections: non-empty table, finite values,
 * t >= 0, k/rate >= 0, ascending time (equal = step). Pure functions only —
 * no DOM, no Tauri — so node --test covers parse / round-trip / guards.
 */

import { upsertDictBlock, removeDictBlock } from './editor.js';

/** The kernel's own value key per schedule kind. */
export const SCHEDULE_KINDS = {
  valve: { block: 'valveSchedule', valueKey: 'k' },
  source: { block: 'sourceSchedule', valueKey: 'rate' },
};

// ---------------------------------------------------------------------------
// text scanning (brace matching with // comment skipping, anchored per brace
// so nested `}` can never truncate a block early — the pitfall the E3 kernel
// checker documented)
// ---------------------------------------------------------------------------

/** Skip a `// ...` comment starting at i; returns the index after it. */
function skipComment(text, i) {
  const eol = text.indexOf('\n', i);
  return eol < 0 ? text.length : eol;
}

/** Match the `}` closing the `{` at `open`; -1 when unbalanced. */
function braceEnd(text, open) {
  let depth = 0;
  for (let i = open; i < text.length; i++) {
    if (text[i] === '/' && text[i + 1] === '/') {
      i = skipComment(text, i);
      continue;
    }
    if (text[i] === '{') depth++;
    else if (text[i] === '}' && --depth === 0) return i;
  }
  return -1;
}

/** Locate a named block `name { ... }` from `from`; null when absent. */
function findBlock(text, name, from = 0) {
  const m = new RegExp(`\\b${name}\\s*\\{`).exec(text.slice(from));
  if (!m) return null;
  const nameStart = from + m.index;
  const open = nameStart + m[0].length - 1;
  const end = braceEnd(text, open);
  return end < 0 ? null : { body: text.slice(open + 1, end), nameStart, open, end };
}

/** Iterate the named entries (`NAME { ... }`) at one nesting level. */
function namedEntries(body) {
  const out = [];
  let cursor = 0;
  while (cursor < body.length) {
    const m = /\b(\w+)\s*\{/.exec(body.slice(cursor));
    if (!m) break;
    const start = cursor + m.index;
    const open = start + m[0].length - 1;
    const end = braceEnd(body, open);
    if (end < 0) break;
    // skip matches that sit inside an entry already collected (nested blocks)
    const prev = out[out.length - 1];
    if (!prev || start > prev.end) {
      out.push({ name: m[1], body: body.slice(open + 1, end), start, open, end });
    }
    cursor = end + 1;
  }
  return out;
}

/** Parse the anonymous event blocks `{ t <t>; <key> <v>; }` in order. */
function parseEventBlocks(body, valueKey) {
  const events = [];
  let cursor = 0;
  while (cursor < body.length) {
    const open = body.indexOf('{', cursor);
    if (open < 0) break;
    const end = braceEnd(body, open);
    if (end < 0) break;
    const inner = body.slice(open + 1, end);
    const num = (key) => {
      const m = new RegExp(`\\b${key}\\s+([-+0-9.eE]+)\\s*;`).exec(inner);
      return m ? Number(m[1]) : undefined;
    };
    const t = num('t');
    const v = num(valueKey);
    // lenient parse: malformed events the kernel would reject are skipped
    // here and re-caught loudly by scheduleBlockLines on generation
    if (Number.isFinite(t) && Number.isFinite(v)) events.push({ t, [valueKey]: v });
    cursor = end + 1;
  }
  return events;
}

// ---------------------------------------------------------------------------
// parse
// ---------------------------------------------------------------------------

/** Parse the case dictionary into `{ valves: {pipe: [{t,k}]}, source }`. */
export function parseSchedules(caseText) {
  const valves = {};
  const pipesBlock = findBlock(caseText, 'pipes');
  if (pipesBlock) {
    for (const entry of namedEntries(pipesBlock.body)) {
      const vs = findBlock(entry.body, 'valveSchedule');
      if (!vs) continue;
      const ev = findBlock(vs.body, 'events');
      if (!ev) continue;
      valves[entry.name] = parseEventBlocks(ev.body, 'k');
    }
  }
  let source = null;
  const ss = findBlock(caseText, 'sourceSchedule');
  if (ss) {
    const nodeM = /\bnode\s+(\w+)\s*;/.exec(ss.body);
    const ev = findBlock(ss.body, 'events');
    if (nodeM && ev) {
      source = { node: nodeM[1], events: parseEventBlocks(ev.body, 'rate') };
    }
  }
  return { valves, source };
}

/** All pipe entry names declared in the case's pipes{} block. */
export function pipeNames(caseText) {
  const pipesBlock = findBlock(caseText, 'pipes');
  return pipesBlock ? namedEntries(pipesBlock.body).map((e) => e.name) : [];
}

// ---------------------------------------------------------------------------
// validate + generate
// ---------------------------------------------------------------------------

/**
 * Kernel-aligned event-table validation (transient.cpp loader semantics).
 * Throws with a loud message on: empty table, non-finite t / value, negative
 * t / value, descending time (equal times are a STEP and stay legal).
 */
function validateEvents(events, valueKey, label) {
  if (!Array.isArray(events) || events.length === 0) {
    throw new Error(`${label}: events 表不能为空`);
  }
  let prevT = -Infinity;
  for (let i = 0; i < events.length; i++) {
    const t = Number(events[i].t);
    const v = Number(events[i][valueKey]);
    if (!Number.isFinite(t)) throw new Error(`${label}: 事件 ${i + 1} 的 t 必须是有限数值`);
    if (!Number.isFinite(v)) throw new Error(`${label}: 事件 ${i + 1} 的 ${valueKey} 必须是有限数值`);
    if (t < 0) throw new Error(`${label}: 事件 ${i + 1} 的 t 必须 >= 0（步进从 t = 0 开始）`);
    if (v < 0) throw new Error(`${label}: 事件 ${i + 1} 的 ${valueKey} 必须 >= 0`);
    if (i > 0 && t < prevT) {
      throw new Error(`${label}: 事件必须按时间升序（相等 t 为 STEP）：#${i} ${prevT} -> #${i + 1} ${t}`);
    }
    prevT = t;
  }
}

/** Compact number text, mirroring the engine's dictionary style. */
function fmt(v) {
  return String(Number(v));
}

/**
 * Validate one schedule and return its exact inner dictionary lines (the
 * content that sits inside the `valveSchedule{}` / `sourceSchedule{}` block,
 * unindented — callers add their own nesting). kind is 'valve' or 'source'.
 */
export function scheduleBlockLines(kind, events, node) {
  const spec = SCHEDULE_KINDS[kind];
  if (!spec) throw new Error(`scheduleEditor: unknown kind '${kind}'`);
  validateEvents(events, spec.valueKey, spec.block);
  const lines = [];
  if (kind === 'source') {
    if (!node || !/^\w+$/.test(String(node))) {
      throw new Error('sourceSchedule: node 必须是已声明的节点名');
    }
    lines.push(`node ${node};`);
  }
  lines.push('events', '{');
  for (const ev of events) {
    lines.push(`{ t ${fmt(ev.t)}; ${spec.valueKey} ${fmt(ev[spec.valueKey])}; }`);
  }
  lines.push('}');
  return lines;
}

// ---------------------------------------------------------------------------
// apply (atomic block upsert / removal)
// ---------------------------------------------------------------------------

/** Replace the pipe entry's valveSchedule block, or insert one when absent. */
function upsertValveSchedule(caseText, pipe, events) {
  const pipesBlock = findBlock(caseText, 'pipes');
  if (!pipesBlock) throw new Error('valveSchedule: case 字典没有 pipes{} 块');
  const entry = namedEntries(pipesBlock.body).find((e) => e.name === pipe);
  if (!entry) throw new Error(`valveSchedule: pipes{} 里没有管段 '${pipe}'`);
  // mutual exclusion (kernel cfFatal §5.170): the event table prices K(t)
  // itself, so the single-cycle keys would be read by nothing
  if (/\bvalveClose\b|\bvalveOpen\b/.test(entry.body)
      || /(^|[\s;])valve\s+[-+0-9.eE]/.test(entry.body.replace(/\bvalveSchedule\b/g, ' '))) {
    throw new Error(`valveSchedule: 管段 '${pipe}' 已声明 valve/valveClose/valveOpen，单一日程定价一个阀`);
  }
  const lines = scheduleBlockLines('valve', events);
  const base = pipesBlock.open + 1;
  // entry.body starts right after the entry's `{`, so every offset inside it
  // maps into the full text from base + entry.open + 1
  const bodyStart = base + entry.open + 1;
  const entryEnd = base + entry.end;
  const vs = findBlock(entry.body, 'valveSchedule');
  const ind = '        '; // entry content sits at 4, valveSchedule at 8
  if (vs) {
    const vsStart = bodyStart + vs.nameStart;
    const vsEnd = bodyStart + vs.end;
    const block = [`${ind}valveSchedule`, `${ind}{`,
      ...lines.map((l) => `${ind}    ${l}`), `${ind}}`].join('\n');
    return caseText.slice(0, vsStart) + block + caseText.slice(vsEnd + 1);
  }
  const block = `\n${ind}valveSchedule\n${ind}{\n${
    lines.map((l) => `${ind}    ${l}`).join('\n')}\n${ind}}\n    `;
  return caseText.slice(0, entryEnd) + block + caseText.slice(entryEnd);
}

/** Remove the pipe entry's valveSchedule block entirely (absent -> unchanged). */
function removeValveSchedule(caseText, pipe) {
  const pipesBlock = findBlock(caseText, 'pipes');
  if (!pipesBlock) return caseText;
  const entry = namedEntries(pipesBlock.body).find((e) => e.name === pipe);
  if (!entry) return caseText;
  const vs = findBlock(entry.body, 'valveSchedule');
  if (!vs) return caseText;
  const bodyStart = pipesBlock.open + 1 + entry.open + 1;
  const lineStart = caseText.lastIndexOf('\n', bodyStart + vs.nameStart) + 1;
  return caseText.slice(0, lineStart) + caseText.slice(bodyStart + vs.end + 1);
}

/** Node names declared in the case's nodes{} block. */
function nodeNames(caseText) {
  const nodesBlock = findBlock(caseText, 'nodes');
  return nodesBlock ? namedEntries(nodesBlock.body).map((e) => e.name) : [];
}

/**
 * Apply the whole schedule model to a case dictionary. Valve schedules are
 * replaced/inserted inside their pipe entries (an empty event list removes
 * the pipe's valveSchedule block entirely); the sourceSchedule block is
 * upserted / removed at the top level via upsertDictBlock / removeDictBlock.
 * Throws on every guard the kernel's loader would cfFatal on.
 */
export function applySchedules(caseText, model) {
  let out = caseText;
  const valves = (model && model.valves) || {};
  for (const [pipe, events] of Object.entries(valves)) {
    out = events.length
      ? upsertValveSchedule(out, pipe, events)
      : removeValveSchedule(out, pipe);
  }
  const source = model && model.source;
  if (source && source.enabled !== false
      && Array.isArray(source.events) && source.events.length) {
    // kernel guards: the node must exist and must not also be priced by a
    // sources{} entry (one rate account per node)
    const node = String(source.node || '');
    if (!nodeNames(out).includes(node)) {
      throw new Error(`sourceSchedule: 节点 '${node}' 未在 nodes{} 中声明`);
    }
    const sourcesBlock = findBlock(out, 'sources');
    if (sourcesBlock) {
      for (const entry of namedEntries(sourcesBlock.body)) {
        const m = /\bnode\s+(\w+)\s*;/.exec(entry.body);
        if (m && m[1] === node) {
          throw new Error(`sourceSchedule: 节点 '${node}' 已被 sources{} 条目定价（每节点一个 rate 账户）`);
        }
      }
    }
    const lines = scheduleBlockLines('source', source.events, node);
    out = upsertDictBlock(out, 'sourceSchedule', lines);
  } else {
    out = removeDictBlock(out, 'sourceSchedule');
  }
  return out;
}

// ---------------------------------------------------------------------------
// ESD preset
// ---------------------------------------------------------------------------

/**
 * ESD (emergency shutdown) preset: every DECLARED valve schedule becomes a
 * linkage table that holds its own first-table K up to `esdTime`, then STEPs
 * to `closedK` at that instant (the equal-t event pair — the kernel's step
 * semantics fire the shutdown at t + dt). Source schedules pass through
 * untouched. Pure: returns a new model, throws on a non-finite / negative
 * esdTime or closedK.
 */
export function esdPreset(model, esdTime = 5, closedK = 1e12) {
  const t = Number(esdTime);
  const k = Number(closedK);
  if (!Number.isFinite(t) || t < 0) throw new Error('esdPreset: ESD 时刻必须是 >= 0 的有限数值');
  if (!Number.isFinite(k) || k < 0) throw new Error('esdPreset: closedK 必须是 >= 0 的有限数值');
  const valves = {};
  for (const [pipe, events] of Object.entries((model && model.valves) || {})) {
    if (!Array.isArray(events) || events.length === 0) { valves[pipe] = events; continue; }
    const kOpen = Number(events[0].k);
    valves[pipe] = [
      { t: 0, k: kOpen },
      { t, k: kOpen },
      { t, k },
    ];
  }
  return { valves, source: model ? model.source : null };
}

// ---------------------------------------------------------------------------
// timeline visualization (pure SVG string, no DOM)
// ---------------------------------------------------------------------------

/**
 * Piecewise-linear event-table preview as an inline SVG string: x = time,
 * y = log10(value) (a closed K of 1e12 must not flatten a K = 20 baseline),
 * equal-t adjacent events draw the vertical STEP edge naturally. Event points
 * get marker circles. Returns '' for an empty table.
 */
export function timelineSvg(events, valueKey, width = 480, height = 90) {
  if (!Array.isArray(events) || events.length === 0) return '';
  const pad = { l: 34, r: 10, t: 8, b: 14 };
  const iw = Math.max(10, width - pad.l - pad.r);
  const ih = Math.max(10, height - pad.t - pad.b);
  const ts = events.map((e) => Number(e.t));
  const vs = events.map((e) => Math.log10(Math.max(Number(e[valueKey]), 1e-30)));
  const t0 = Math.min(0, ...ts);
  const t1 = ts[ts.length - 1] > t0 ? ts[ts.length - 1] : t0 + 1;
  const v0 = Math.min(...vs);
  const v1 = Math.max(...vs) > v0 ? Math.max(...vs) : v0 + 1;
  const x = (t) => pad.l + iw * ((t - t0) / (t1 - t0));
  const y = (v) => pad.t + ih * (1 - (v - v0) / (v1 - v0));
  const pts = events.map((e, i) => `${x(ts[i]).toFixed(1)},${y(vs[i]).toFixed(1)}`).join(' ');
  const dots = events.map((e, i) =>
    `<circle class="sched-dot" cx="${x(ts[i]).toFixed(1)}" cy="${y(vs[i]).toFixed(1)}" r="2.6"><title>t ${fmt(e.t)} · ${valueKey} ${fmt(e[valueKey])}</title></circle>`,
  ).join('');
  return `<svg class="sched-svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`
    + `<line class="sched-axis" x1="${pad.l}" y1="${height - pad.b}" x2="${width - pad.r}" y2="${height - pad.b}" />`
    + `<polyline class="sched-line" points="${pts}" />${dots}</svg>`;
}
