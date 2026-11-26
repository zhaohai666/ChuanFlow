/**
 * E2b — structured editing of the steady-network declaration blocks that the
 * flow-assurance module does not cover:
 *
 *   correlations{}  per-pipe holdup law / forced flow pattern / HB Griffith
 *                   (network.cpp:3431-3469; scope guards 1228-1269)
 *   rheology{}     per-pipe power-law / Herschel-Bulkley liquid law
 *                   (network.cpp:3813-3836)
 *   blackoil{}     nested under fluid{}: black-oil PVT declaration
 *                   (network.cpp:3042-3056; BlackOil.h:82-89)
 *
 * correlations and rheology are NOT top-level blocks: their keys live inside
 * `pipes { <name> { ... } }`, so this module upserts/removes keys inside the
 * named pipe entry rather than replacing a whole top-level block. blackoil is
 * a nested sub-block of `fluid {}`.
 *
 * Pure functions only — no DOM, no Tauri. The schema, guards and renderers
 * mirror the kernel's own dictionary reads line for line; node --test covers
 * key alignment, round-trip through the GUI parser and loud rejections.
 */

// ---------------------------------------------------------------------------
// Schemas. kind: 'num' finite scalar | 'word' bare token with options.
// required: the GUI refuses to emit the block without it. `check` mirrors the
// kernel guard as (values) -> error string | null.
// ---------------------------------------------------------------------------

export const PIPE_BLOCKS = {
  correlations: {
    label: '持液率/流型关联式（per pipe）',
    fields: [
      { key: 'correlation', kind: 'word', required: false,
        options: ['beggsBrill', 'hagedornBrown', 'dunsRos', 'olgasSteady'],
        hint: '持液率定律；缺省 Beggs-Brill' },
      { key: 'flowPattern', kind: 'num', required: false,
        hint: '流型索引 -1(地图)/0 分离/2 间歇/3 分散；1 过渡型内核拒绝' },
      { key: 'hbGriffith', kind: 'word', required: false,
        options: ['true', 'false'], hint: 'HB 分支 Grifth 图版，仅可配 hagedornBrown' },
    ],
    check: (v) => {
      if (v.correlation !== undefined && v.flowPattern !== undefined) {
        return 'correlation 与 flowPattern 不能同时声明（内核二者互斥）';
      }
      if (v.hbGriffith === 'true' && v.correlation !== 'hagedornBrown') {
        return 'hbGriffith 仅可与 correlation hagedornBrown 同用';
      }
      if (v.flowPattern !== undefined) {
        const fp = Number(v.flowPattern);
        if (!Number.isInteger(fp) || ![-1, 0, 2, 3].includes(fp)) {
          return 'flowPattern 只能是 -1/0/2/3（1 为过渡型，内核拒绝）';
        }
      }
      return null;
    },
  },
  rheology: {
    label: '非牛顿流变（per pipe）',
    fields: [
      { key: 'flowIndex', kind: 'num', required: false,
        hint: '幂律指数 n ≥0，0=牛顿流体' },
      { key: 'consistency', kind: 'num', required: false,
        hint: '稠度系数 K (Pa·s^n)；n>0 时必须 >0，n=0 时不得声明' },
      { key: 'yieldStress', kind: 'num', required: false,
        hint: '屈服应力 τ0 ≥0，0=无屈服' },
    ],
    check: (v) => {
      const n = v.flowIndex === undefined ? 0 : Number(v.flowIndex);
      if (!(n >= 0)) return 'rheology.flowIndex 不可为负';
      const K = v.consistency === undefined ? 0 : Number(v.consistency);
      if (n > 0 && !(K > 0)) return 'flowIndex>0 时 consistency 必须为正';
      if (n === 0 && K > 0) return '未声明 flowIndex 时 consistency 不被读取（内核拒绝）';
      const t = v.yieldStress === undefined ? 0 : Number(v.yieldStress);
      if (!(t >= 0)) return 'rheology.yieldStress 不可为负';
      return null;
    },
  },
};

export const PIPE_BLOCK_NAMES = Object.keys(PIPE_BLOCKS);

export const BLACKOIL_SPEC = {
  label: '黑油 PVT（fluid.blackoil{}）',
  // network.cpp:3048-3053 reads these six keys; the only kernel-fatal range
  // check is oilDensity > 0 (3054-3055). The remaining GUI guards are
  // conservative (positive/non-negative, fraction bounds); flagged in the
  // delivery report for the author to relax if the kernel later allows it.
  fields: [
    { key: 'gasGravity', kind: 'num', required: false,
      hint: '天然气相对密度（空气=1），默认 0.65' },
    { key: 'oilDensity', kind: 'num', required: true,
      hint: '储罐油密度 (kg/m³)，必须 >0' },
    { key: 'gor', kind: 'num', required: false,
      hint: '生产气油比 (sm³/sm³)，默认 0' },
    { key: 'waterCut', kind: 'num', required: false,
      hint: '含水率 [0,1)，默认 0' },
    { key: 'temperature', kind: 'num', required: false,
      hint: '等温温度 (K)，默认 300' },
    { key: 'surfaceTension', kind: 'num', required: false,
      hint: '气液表面张力 (N/m)，默认 0.03' },
  ],
  check: (v) => {
    if (!(Number(v.oilDensity) > 0)) return 'blackoil.oilDensity 必须为正';
    if (v.gasGravity !== undefined && !(Number(v.gasGravity) > 0)) {
      return 'blackoil.gasGravity 必须为正';
    }
    if (v.gor !== undefined && !(Number(v.gor) >= 0)) return 'blackoil.gor 不可为负';
    const wc = v.waterCut === undefined ? 0 : Number(v.waterCut);
    if (!(wc >= 0 && wc < 1)) return 'blackoil.waterCut 必须在 [0,1)';
    if (v.temperature !== undefined && !(Number(v.temperature) > 0)) {
      return 'blackoil.temperature 必须为正';
    }
    if (v.surfaceTension !== undefined && !(Number(v.surfaceTension) >= 0)) {
      return 'blackoil.surfaceTension 不可为负';
    }
    return null;
  },
};

// ---------------------------------------------------------------------------
// Brace-aware dictionary scanning (skips // comments), mirroring editor.js.
// ---------------------------------------------------------------------------

function stripComment(line) {
  return line.replace(/\/\/.*$/, '');
}

/** Index of the matching close brace for the opener at `open`; -1 on mismatch. */
function braceEnd(text, open) {
  let depth = 0;
  for (let i = open; i < text.length; i += 1) {
    if (text[i] === '/' && text[i + 1] === '/') {
      const eol = text.indexOf('\n', i);
      i = eol < 0 ? text.length : eol;
      continue;
    }
    if (text[i] === '{') depth += 1;
    else if (text[i] === '}' && (depth -= 1) === 0) return i;
  }
  return -1;
}

/** Locate a top-level `name { ... }` block; returns body span or null. */
function sectionSpan(text, name) {
  const m = new RegExp(`\\b${name}\\s*\\{`).exec(text);
  if (!m) return null;
  const open = m.index + m[0].length - 1;
  const end = braceEnd(text, open);
  if (end < 0) return null;
  return { body: text.slice(open + 1, end), start: open + 1, end };
}

/** Locate a nested `entryName { ... }` inside a section body span. */
function entrySpan(text, section, entryName) {
  const m = new RegExp(`\\b${entryName}\\s*\\{`).exec(section.body);
  if (!m) return null;
  const open = section.start + m.index + m[0].length - 1;
  const end = braceEnd(text, open);
  if (end < 0) return null;
  return { body: text.slice(open + 1, end), start: open + 1, end, opener: section.start + m.index };
}

/** Parse a `key value;` body into an ordered key/value map. */
function parseBodyKv(body) {
  const kv = {};
  const order = [];
  for (const stmt of body.split(';')) {
    const t = stripComment(stmt).trim();
    if (!t) continue;
    const sp = t.split(/\s+/);
    kv[sp[0]] = sp.slice(1).join(' ');
    order.push(sp[0]);
  }
  return { kv, order };
}

// ---------------------------------------------------------------------------
// Schema-driven validation + line rendering (shared by pipe and blackoil).
// ---------------------------------------------------------------------------

/** Validate one schema's values and render its inner dictionary lines. */
function renderSpecLines(spec, values, blockLabel) {
  const v = { ...(values || {}) };
  if (spec.check) {
    const err = spec.check(v);
    if (err) throw new Error(err);
  }
  const lines = [];
  for (const f of spec.fields) {
    const raw = v[f.key];
    if (raw === undefined || raw === null || raw === '') {
      if (f.required) throw new Error(`${blockLabel}.${f.key} 是必需字段`);
      continue;
    }
    const s = String(raw).trim();
    if (f.kind === 'num') {
      const n = Number(s);
      if (!Number.isFinite(n)) throw new Error(`${blockLabel}.${f.key} 必须是有限数值（got ${s}）`);
      lines.push(`${f.key} ${s};`);
    } else {
      if (f.options && !f.options.includes(s)) {
        throw new Error(`${blockLabel}.${f.key} 只能取 ${f.options.join('/')}`);
      }
      lines.push(`${f.key} ${s};`);
    }
  }
  return lines;
}

/** Public renderer for a per-pipe block. */
export function renderPipeBlockLines(blockId, values) {
  const spec = PIPE_BLOCKS[blockId];
  if (!spec) throw new Error(`dictBlocks: unknown pipe block '${blockId}'`);
  return renderSpecLines(spec, values, blockId);
}

/** Public renderer for the fluid.blackoil block. */
export function renderBlackoilLines(values) {
  return renderSpecLines(BLACKOIL_SPEC, values, 'blackoil');
}

// ---------------------------------------------------------------------------
// Per-pipe upsert / removal. Keys of the target block are reconciled inside
// the named entry; every other statement (geometry, other device keys) keeps
// its original text verbatim.
// ---------------------------------------------------------------------------

/**
 * Apply (or remove) one per-pipe declaration block.
 * @param {string} caseText dictionary text
 * @param {string} pipeName entry name under `pipes {}`
 * @param {string} blockId `correlations` | `rheology`
 * @param {object|null} values validated form values; null/undefined removes
 *   every schema key of that block from the entry.
 */
export function applyPipeBlock(caseText, pipeName, blockId, values) {
  const spec = PIPE_BLOCKS[blockId];
  if (!spec) throw new Error(`dictBlocks: unknown pipe block '${blockId}'`);
  const section = sectionSpan(caseText, 'pipes');
  if (!section) throw new Error('字典缺少 pipes{} 块，无法定位管段');
  const entry = entrySpan(caseText, section, pipeName);
  if (!entry) throw new Error(`pipes{} 中找不到管段 '${pipeName}'`);

  const schemaKeys = spec.fields.map((f) => f.key);
  const { kv, order } = parseBodyKv(entry.body);
  const next = { ...kv };

  if (values === null || values === undefined) {
    for (const k of schemaKeys) delete next[k];
  } else {
    const lines = renderSpecLines(spec, values, blockId);
    const emitted = new Set();
    for (const line of lines) {
      const sp = line.split(/\s+/);
      next[sp[0]] = sp.slice(1).join(' ').replace(/;$/, '');
      emitted.add(sp[0]);
    }
    for (const k of schemaKeys) if (!emitted.has(k)) delete next[k];
  }

  const outOrder = [...order.filter((k) => next[k] !== undefined)];
  for (const k of Object.keys(next)) {
    if (!outOrder.includes(k)) outOrder.push(k);
  }
  const newBody = outOrder.map((k) => `${k} ${next[k]};`).join(' ');
  const newEntry = `${pipeName} { ${newBody} }`;
  return caseText.slice(0, entry.opener) + newEntry + caseText.slice(entry.end + 1);
}

/** Parse a named pipe entry's key/value map (GUI's own dict parser). */
export function parsePipeEntry(caseText, pipeName) {
  const section = sectionSpan(caseText, 'pipes');
  if (!section) return null;
  const entry = entrySpan(caseText, section, pipeName);
  if (!entry) return null;
  return parseBodyKv(entry.body).kv;
}

// ---------------------------------------------------------------------------
// fluid.blackoil nested upsert / removal.
// ---------------------------------------------------------------------------

/**
 * Apply (or remove) the `fluid { blackoil { ... } }` sub-block. When the
 * values are null/undefined the blackoil sub-block is dropped, leaving fluid's
 * other keys (rho/mu/pvt/etc.) intact. Throws on the schema guards.
 */
export function applyBlackoil(caseText, values) {
  const fluid = sectionSpan(caseText, 'fluid');

  if (values === null || values === undefined) {
    if (!fluid) return caseText;
    const bo = entrySpan(caseText, fluid, 'blackoil');
    if (!bo) return caseText;
    // drop the sub-block span, then collapse a leftover blank line
    let out = caseText.slice(0, bo.opener) + caseText.slice(bo.end + 1);
    out = out.replace(/\n[ \t]*\n([ \t]*\})/, '\n$1');
    return out;
  }

  const lines = renderSpecLines(BLACKOIL_SPEC, values, 'blackoil');
  const inner = lines.map((s) => `        ${s}`).join('\n');
  const block = `blackoil\n    {\n${inner}\n    }`;

  if (!fluid) {
    const sep = caseText.endsWith('\n') ? '' : '\n';
    return caseText + sep + 'fluid\n{\n    ' + block + '\n}\n';
  }

  const bo = entrySpan(caseText, fluid, 'blackoil');
  if (bo) {
    const newBo = `blackoil\n    {\n${inner}\n    }`;
    return caseText.slice(0, bo.opener) + newBo + caseText.slice(bo.end + 1);
  }
  // append the sub-block at the end of the fluid body
  const newBody = fluid.body.replace(/\s*$/, '') + `\n    ${block}\n`;
  return caseText.slice(0, fluid.start) + newBody + caseText.slice(fluid.end);
}

/** Parse the fluid.blackoil sub-block into a key/value map; null if absent. */
export function parseBlackoil(caseText) {
  const fluid = sectionSpan(caseText, 'fluid');
  if (!fluid) return null;
  const bo = entrySpan(caseText, fluid, 'blackoil');
  if (!bo) return null;
  return parseBodyKv(bo.body).kv;
}

/**
 * Whether `fluid {}` already declares a `pvt {}` sub-block. The kernel refuses
 * blackoil and pvt together (network.cpp:3040-3041); the panel uses this to
 * warn before the user fills the blackoil form.
 */
export function fluidHasPvt(caseText) {
  const fluid = sectionSpan(caseText, 'fluid');
  return fluid ? /\bpvt\s*\{/.test(fluid.body) : false;
}
