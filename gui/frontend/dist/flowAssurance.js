/**
 * E2 — structured editing of the 1D two-fluid flow-assurance declaration
 * blocks: hydrate{}, inhibitor{}, inhibitorEquilibrium{}, slugTracking{},
 * purgeWave{}, wax{}, pigging{}, sand{} and acoustics{}. All kernel keys and
 * guards mirror src/1d/src/twoFluid.cpp (§5.127/§5.141/§5.161/§5.162/§5.163
 * /§5.96 wax /§5.176 pigging /G10 sand /§5.183 acoustics). Pure functions
 * only: parse a case dictionary into a form model, render the form back into
 * exact dictionary lines (with loud validation), apply/remove blocks via
 * upsertDictBlock / removeDictBlock, and parse the kernel's report rows for
 * the slug timeline. No DOM, no Tauri — node --test covers round-trip, guards
 * and report parsing.
 */

import { upsertDictBlock, removeDictBlock } from './editor.js';

/**
 * Per-block field schema. kind: 'num' (finite scalar) | 'word' (bare token).
 * required fields must be present in the generated block; optional fields are
 * emitted only when the form carries a non-empty value. `check` carries the
 * kernel guard as (value, form) -> error string | null.
 */
export const FA_BLOCKS = {
  hydrate: {
    label: '水合物生成/分解（hydrate{}）',
    fields: [
      { key: 'curveA', kind: 'num', required: true, hint: 'ln(p) 平衡线斜率' },
      { key: 'curveB', kind: 'num', required: true, hint: 'ln(p) 平衡线截距 (K)' },
      { key: 'growthRate', kind: 'num', required: true, hint: '生成速率 1/s' },
      { key: 'dissociationRate', kind: 'num', required: true, hint: '分解速率 1/s' },
      { key: 'gasContent', kind: 'num', required: true, hint: '水合物气体含量 (mol/m3)' },
      { key: 'latentHeat', kind: 'num', required: true, hint: '生成潜热 (J/m3)' },
      { key: 'initial', kind: 'num', required: false, hint: '初始水合物分数 [0,1]' },
    ],
    check: (v) => (v.initial !== undefined && (v.initial < 0 || v.initial > 1)
      ? 'hydrate.initial 必须在 [0,1]' : null),
  },
  inhibitor: {
    label: '抑制剂温降（inhibitor{}）',
    fields: [
      { key: 'kDepression', kind: 'num', required: true, hint: 'Hammerschmidt 温降常数' },
      { key: 'initial', kind: 'num', required: false, hint: '初始抑制剂浓度' },
    ],
  },
  inhibitorEquilibrium: {
    label: '抑制剂严格活度相平衡（inhibitorEquilibrium{}）',
    fields: [
      { key: 'enable', kind: 'word', required: false, options: ['true', 'false'], hint: 'false 退 Hammerschmidt 路径' },
      { key: 'activityK', kind: 'num', required: true, hint: '水活度系数 aW=1/(1+K·w)' },
      { key: 'basis', kind: 'word', required: false, options: ['liquid', 'water'], hint: 'w 的取值基准' },
    ],
    check: (v) => {
      if (!(v.activityK > 0)) return 'inhibitorEquilibrium.activityK 必须为正';
      if (v.enable && !['true', 'false'].includes(String(v.enable))) return 'enable 只能 true/false';
      if (v.basis && !['liquid', 'water'].includes(String(v.basis))) return 'basis 只能 liquid/water';
      return null;
    },
  },
  slugTracking: {
    label: '段塞追踪（slugTracking{}）',
    fields: [
      { key: 'enable', kind: 'word', required: false, options: ['true', 'false'], hint: 'false 关闭追踪' },
      { key: 'alphaThreshold', kind: 'num', required: false, hint: '持液率阈值 (0,1)，默认 0.5' },
      { key: 'minLength', kind: 'num', required: false, hint: '最小段塞长度 (m)，默认 0' },
      { key: 'velocityMode', kind: 'word', required: false, options: ['constant', 'gasCap'], hint: '前缘速度模式，默认 constant' },
      { key: 'slugVelocity', kind: 'num', required: false, hint: 'constant 模式段塞速度 (m/s)' },
      { key: 'gasCapTime', kind: 'num', required: false, hint: 'gasCap 模式气帽松弛时间 (s)' },
      { key: 'maxSlugVelocity', kind: 'num', required: false, hint: 'gasCap 模式速度上限，默认 0' },
    ],
    check: (v) => {
      const mode = String(v.velocityMode || 'constant');
      const th = v.alphaThreshold === undefined ? 0.5 : Number(v.alphaThreshold);
      if (!(th > 0 && th < 1)) return 'slugTracking.alphaThreshold 必须在 (0,1)';
      if (v.minLength !== undefined && Number(v.minLength) < 0) return 'minLength 不可为负';
      if (!['constant', 'gasCap'].includes(mode)) return 'velocityMode 只能 constant/gasCap';
      if (mode === 'constant' && !(Number(v.slugVelocity) > 0)) return 'constant 模式需正的 slugVelocity';
      if (mode === 'gasCap' && !(Number(v.gasCapTime) > 0)) return 'gasCap 模式需正的 gasCapTime';
      if (v.gasCapTime !== undefined && !(Number(v.gasCapTime) > 0) && mode === 'gasCap')
        return 'gasCap 模式需正的 gasCapTime';
      return null;
    },
  },
  purgeWave: {
    label: '立管吹扫波（purgeWave{}）',
    fields: [
      { key: 'enable', kind: 'word', required: false, options: ['true', 'false'], hint: 'false 关闭吹扫' },
      { key: 'triggerFactor', kind: 'num', required: false, hint: '触发倍率，默认 1' },
      { key: 'continuityAlpha', kind: 'num', required: false, hint: '连续持液阈值 (0,1)，默认 0.5' },
      { key: 'alphaResidual', kind: 'num', required: false, hint: '吹干残余持液 [0,1)，默认 0.05' },
      { key: 'waveMode', kind: 'word', required: false, options: ['closed', 'constant'], hint: '波前速度模式，默认 closed' },
      { key: 'uWave', kind: 'num', required: false, hint: 'constant 模式波速 (m/s)' },
      { key: 'gasCapTime', kind: 'num', required: true, hint: '气帽压力松弛时间 (s)，必须为正' },
      { key: 'maxCycles', kind: 'num', required: false, hint: '最大周期数，默认 0（不限）' },
      { key: 'riserStart', kind: 'num', required: false, hint: '立管窗口起点 (m)' },
      { key: 'riserEnd', kind: 'num', required: false, hint: '立管窗口终点 (m)' },
    ],
    check: (v) => {
      const ca = v.continuityAlpha === undefined ? 0.5 : Number(v.continuityAlpha);
      const ar = v.alphaResidual === undefined ? 0.05 : Number(v.alphaResidual);
      if (!(ca > 0 && ca < 1)) return 'purgeWave.continuityAlpha 必须在 (0,1)';
      if (!(ar >= 0 && ar < 1)) return 'purgeWave.alphaResidual 必须在 [0,1)';
      if (ar >= ca) return 'alphaResidual 必须小于 continuityAlpha';
      const mode = String(v.waveMode || 'closed');
      if (!['closed', 'constant'].includes(mode)) return 'waveMode 只能 closed/constant';
      if (mode === 'constant' && !(Number(v.uWave) > 0)) return 'constant 模式需正的 uWave';
      if (!(Number(v.gasCapTime) > 0)) return 'purgeWave 需正的 gasCapTime';
      if (v.riserEnd !== undefined && Number(v.riserEnd) > 0
          && Number(v.riserEnd) <= Number(v.riserStart || 0))
        return 'riserEnd 必须大于 riserStart';
      return null;
    },
  },
  wax: {
    label: '蜡沉积/溶解（wax{}）',
    // twoFluid.cpp:3481-3498 (keys), guards 396-413. Kernel defaults in
    // parentheses; the block also requires a thermal{} block to be declared.
    fields: [
      { key: 'meltingTemperature', kind: 'num', required: false, hint: '析蜡温度 Tf (K)，默认 323.75，必须 >0' },
      { key: 'fusionEnthalpy', kind: 'num', required: false, hint: '溶解焓 dHf (J/mol)，默认 54890，必须 >0' },
      { key: 'waxFraction', kind: 'num', required: false, hint: '蜡含量 z ∈ [0,1]，默认 0' },
      { key: 'growthRate', kind: 'num', required: false, hint: '沉积速率 ≥0 1/s，默认 0' },
      { key: 'dissolutionRate', kind: 'num', required: false, hint: '溶解速率 ≥0 1/s，默认 0' },
      { key: 'latentHeat', kind: 'num', required: false, hint: '析蜡潜热 ≥0，默认 0' },
      { key: 'initial', kind: 'num', required: false, hint: '初始蜡分数 [0,1]，默认 0' },
      { key: 'feedback', kind: 'word', required: false, options: ['true', 'false'], hint: '是否按沉积粗糙度反馈水动力' },
      { key: 'depositRoughness', kind: 'num', required: false, hint: '沉积层粗糙度 ≥0，默认 0' },
    ],
    check: (v) => {
      const Tf = v.meltingTemperature === undefined ? 323.75 : Number(v.meltingTemperature);
      const dHf = v.fusionEnthalpy === undefined ? 54890 : Number(v.fusionEnthalpy);
      if (!(Tf > 0)) return 'wax.meltingTemperature 必须为正';
      if (!(dHf > 0)) return 'wax.fusionEnthalpy 必须为正';
      const z = v.waxFraction === undefined ? 0 : Number(v.waxFraction);
      if (!(z >= 0 && z <= 1)) return 'wax.waxFraction 必须在 [0,1]';
      const g = v.growthRate === undefined ? 0 : Number(v.growthRate);
      const d = v.dissolutionRate === undefined ? 0 : Number(v.dissolutionRate);
      if (!(g >= 0)) return 'wax.growthRate 不可为负';
      if (!(d >= 0)) return 'wax.dissolutionRate 不可为负';
      const i = v.initial === undefined ? 0 : Number(v.initial);
      if (!(i >= 0 && i <= 1)) return 'wax.initial 必须在 [0,1]';
      const lh = v.latentHeat === undefined ? 0 : Number(v.latentHeat);
      if (!(lh >= 0)) return 'wax.latentHeat 不可为负';
      const dr = v.depositRoughness === undefined ? 0 : Number(v.depositRoughness);
      if (!(dr >= 0)) return 'wax.depositRoughness 不可为负';
      return null;
    },
  },
  pigging: {
    label: '清管器追踪（pigging{}）',
    // twoFluid.cpp:3866-3902. frontCapacity is kernel-fatal unless >0;
    // pigVelocity must be >0 only when the tracker is enabled.
    fields: [
      { key: 'enable', kind: 'word', required: false, options: ['true', 'false'], hint: 'false 关闭追踪（缺省 true）' },
      { key: 'startX', kind: 'num', required: false, hint: '发射位置 x (m)，默认 0' },
      { key: 'endX', kind: 'num', required: false, hint: '到达位置 x (m)，默认 0' },
      { key: 'pigVelocity', kind: 'num', required: false, hint: '清管器速度 (m/s)，启用时必须 >0' },
      { key: 'blockage', kind: 'num', required: false, hint: '堵塞截面比 [0,1)，默认 0' },
      { key: 'pushEfficiency', kind: 'num', required: false, hint: '堆液效率 [0,1]，默认 0' },
      { key: 'frontCapacity', kind: 'num', required: true, hint: '前端容积 C (>0)' },
    ],
    check: (v) => {
      const enable = v.enable === undefined ? true : v.enable === 'true';
      const sx = Number(v.startX || 0);
      const ex = Number(v.endX || 0);
      if (sx < 0) return 'pigging.startX 不可为负';
      if (ex < 0) return 'pigging.endX 不可为负';
      if (ex > 0 && ex <= sx) return 'pigging.endX 必须大于 startX';
      const vel = Number(v.pigVelocity || 0);
      if (vel < 0) return 'pigging.pigVelocity 不可为负';
      const b = Number(v.blockage || 0);
      if (!(b >= 0 && b < 1)) return 'pigging.blockage 必须在 [0,1)';
      const pe = Number(v.pushEfficiency || 0);
      if (!(pe >= 0 && pe <= 1)) return 'pigging.pushEfficiency 必须在 [0,1]';
      if (!(Number(v.frontCapacity) > 0)) return 'pigging.frontCapacity 必须为正';
      if (enable && !(vel > 0)) return '启用 pigging 需正的 pigVelocity';
      return null;
    },
  },
  sand: {
    label: '砂/岩屑输运（sand{}）',
    // twoFluid.cpp:3589-3605 (keys), guards 454-464. terminalVelocity and
    // criticalVelocity are loud getScalar reads: required and >0.
    fields: [
      { key: 'terminalVelocity', kind: 'num', required: true, hint: '沉降末速 vt (m/s)，必须 >0' },
      { key: 'criticalVelocity', kind: 'num', required: true, hint: '临界流速 uCrit (m/s)，必须 >0' },
      { key: 'erodeRate', kind: 'num', required: false, hint: '床层侵蚀速率 ≥0，默认 0' },
      { key: 'initial', kind: 'num', required: false, hint: '初始悬浮砂分数 [0,1]，默认 0' },
      { key: 'bedInitial', kind: 'num', required: false, hint: '初始床层分数 [0,1)，默认 0' },
      { key: 'feedback', kind: 'word', required: false, options: ['true', 'false'], hint: '床层是否反馈内径' },
    ],
    check: (v) => {
      if (!(Number(v.terminalVelocity) > 0)) return 'sand.terminalVelocity 必须为正';
      if (!(Number(v.criticalVelocity) > 0)) return 'sand.criticalVelocity 必须为正';
      const e = v.erodeRate === undefined ? 0 : Number(v.erodeRate);
      if (!(e >= 0)) return 'sand.erodeRate 不可为负';
      const i = v.initial === undefined ? 0 : Number(v.initial);
      if (!(i >= 0 && i <= 1)) return 'sand.initial 必须在 [0,1]';
      const b = v.bedInitial === undefined ? 0 : Number(v.bedInitial);
      if (!(b >= 0 && b < 1)) return 'sand.bedInitial 必须在 [0,1)';
      return null;
    },
  },
  acoustics: {
    label: '管线声波（acoustics{}）',
    // twoFluid.cpp:3904-4004. soundSpeed (wholesale override) and aModel
    // wood (derived from phase speeds) are mutually exclusive; wood needs a
    // declared positive liquidSoundSpeed aL.
    fields: [
      { key: 'enable', kind: 'word', required: false, options: ['true', 'false'], hint: 'false 关闭声波场（缺省 true）' },
      { key: 'initialGauge', kind: 'num', required: false, hint: '初始表压超压 (Pa)，启用时必须 >0' },
      { key: 'soundSpeed', kind: 'num', required: false, hint: '整体波速覆盖 (m/s)，0=按 aModel 推算；与 wood 互斥' },
      { key: 'aModel', kind: 'word', required: false, options: ['twoPhase', 'eos', 'wood'], hint: '波速模型，默认 twoPhase' },
      { key: 'liquidSoundSpeed', kind: 'num', required: false, hint: '液相声速 aL (m/s)，wood 模型必须 >0' },
      { key: 'cfl', kind: 'num', required: false, hint: 'CFL 数 (0,0.5]，默认 0.5' },
      { key: 'openEnd', kind: 'word', required: false, options: ['outlet', 'inlet'], hint: '压力释放端' },
    ],
    check: (v) => {
      const enable = v.enable === undefined ? true : v.enable === 'true';
      const g = v.initialGauge === undefined ? 0 : Number(v.initialGauge);
      if (g < 0) return 'acoustics.initialGauge 不可为负';
      const a = v.soundSpeed === undefined ? 0 : Number(v.soundSpeed);
      if (a < 0) return 'acoustics.soundSpeed 不可为负';
      const model = String(v.aModel || 'twoPhase');
      if (!['twoPhase', 'eos', 'wood'].includes(model)) return 'aModel 只能 twoPhase/eos/wood';
      const aL = v.liquidSoundSpeed === undefined ? 0 : Number(v.liquidSoundSpeed);
      if (v.liquidSoundSpeed !== undefined && !(aL > 0)) return 'acoustics.liquidSoundSpeed 必须为正';
      if (model === 'wood') {
        if (a > 0) return 'acoustics: soundSpeed 与 aModel wood 互斥';
        if (!(aL > 0)) return 'aModel wood 需正的 liquidSoundSpeed';
      }
      const cfl = v.cfl === undefined ? 0.5 : Number(v.cfl);
      if (!(cfl > 0 && cfl <= 0.5)) return 'acoustics.cfl 必须在 (0,0.5]';
      if (enable && !(g > 0)) return '启用 acoustics 需正的 initialGauge';
      return null;
    },
  },
};

export const FA_BLOCK_NAMES = Object.keys(FA_BLOCKS);

/** Empty form model: every block declared false (no block emitted). */
export function emptyFlowAssurance() {
  const out = {};
  for (const name of FA_BLOCK_NAMES) out[name] = { enabled: false, values: {} };
  return out;
}

/** Parse one top-level block body into a flat key/value map. */
function parseBlock(caseText, name) {
  const lines = caseText.split('\n');
  const nameRe = new RegExp(`^\\s*${name}\\b`);
  for (let i = 0; i < lines.length; i++) {
    const code = lines[i].replace(/\/\/.*$/, '');
    if (!nameRe.test(code)) continue;
    let openIdx = code.includes('{') ? i : -1;
    if (openIdx < 0) {
      for (let j = i + 1; j < lines.length; j++) {
        const nxt = lines[j].replace(/\/\/.*$/, '').trim();
        if (!nxt) continue;
        if (nxt === '{') openIdx = j;
        break;
      }
    }
    if (openIdx < 0) continue;
    const values = {};
    let depth = 0;
    for (let k = openIdx; k < lines.length; k++) {
      const c = lines[k].replace(/\/\/.*$/, '');
      for (const ch of c) {
        if (ch === '{') depth++;
        else if (ch === '}') depth--;
      }
      // Match keys anywhere in the line so inline blocks
      // (`hydrate { curveA 35.2; ... }`) parse the same as multi-line ones.
      // `[^;{}]+` cannot swallow the opener brace or the block name's `{`.
      for (const m of c.matchAll(/(\w+)\s+([^;{}]+);/g)) {
        values[m[1]] = m[2].trim();
      }
      if (depth <= 0) break;
    }
    return values;
  }
  return null;
}

/** Parse the case dictionary into the structured form model. */
export function parseFlowAssurance(caseText) {
  const out = emptyFlowAssurance();
  for (const name of FA_BLOCK_NAMES) {
    const values = parseBlock(caseText, name);
    if (values === null) continue;
    out[name] = { enabled: true, values };
  }
  return out;
}

/** Validate and render one block's inner dictionary lines. Throws on a guard. */
export function flowAssuranceBlockLines(name, values) {
  const spec = FA_BLOCKS[name];
  if (!spec) throw new Error(`flowAssurance: unknown block '${name}'`);
  const v = { ...(values || {}) };
  if (spec.check) {
    const err = spec.check(v);
    if (err) throw new Error(err);
  }
  const lines = [];
  for (const f of spec.fields) {
    let raw = v[f.key];
    if (raw === undefined || raw === null || raw === '') {
      if (f.required) throw new Error(`${name}.${f.key} 是必需字段`);
      continue;
    }
    const s = String(raw).trim();
    if (f.kind === 'num') {
      const n = Number(s);
      if (!Number.isFinite(n)) throw new Error(`${name}.${f.key} 必须是有限数值（got ${s}）`);
      lines.push(`${f.key} ${s};`);
    } else {
      if (f.options && !f.options.includes(s)) throw new Error(`${name}.${f.key} 只能取 ${f.options.join('/')}`);
      lines.push(`${f.key} ${s};`);
    }
  }
  return lines;
}

/**
 * Apply the whole form model to a case dictionary. Enabled blocks are
 * replaced/appended atomically; disabled blocks are removed entirely so the
 * kernel falls back to the byte-identical undeclared path.
 */
export function applyFlowAssurance(caseText, form) {
  // Cross-block guard: the strict activity equilibrium consumes the
  // inhibitor concentration field, so inhibitor{} must also be declared
  // (kernel cfFatal "inhibitorEquilibrium{} ... needs inhibitor{}" §5.163).
  if (form?.inhibitorEquilibrium?.enabled && !form?.inhibitor?.enabled) {
    const declared = parseBlock(caseText, 'inhibitor') !== null;
    if (!declared) {
      throw new Error('inhibitorEquilibrium 需要同时启用 inhibitor{}（抑制剂浓度场）');
    }
  }
  let out = caseText;
  for (const name of FA_BLOCK_NAMES) {
    const entry = form?.[name] || { enabled: false };
    if (entry.enabled) {
      const lines = flowAssuranceBlockLines(name, entry.values);
      out = upsertDictBlock(out, name, lines);
    } else {
      out = removeDictBlock(out, name);
    }
  }
  return out;
}

/**
 * Parse the kernel's `slugtrack` report rows into a timeline of events.
 * NONE rows are skipped; GEN/MERGE/OUT carry id, time, position and (for
 * OUT) out-volume/duration. This drives the slug page's event markers.
 */
export function parseSlugEvents(out) {
  const events = [];
  if (!out) return events;
  for (const line of out.split('\n')) {
    const tok = line.replace(/^[^[]*\[chflow\]\s*/, '').trim().split(/\s+/);
    if (tok[0] !== 'slugtrack') continue;
    const kindIdx = tok.indexOf('event');
    if (kindIdx < 0) continue;
    const kind = tok[kindIdx + 1];
    if (!kind || kind === 'NONE' || kind === 'count') continue;
    const num = (key) => {
      const i = tok.indexOf(key);
      return i >= 0 ? Number(tok[i + 1]) : undefined;
    };
    const ev = {
      kind,
      id: Number(tok[1]),
      t: num('t'),
      xR: num('xR'),
      xF: num('xF'),
      L: num('L'),
      inv: num('inv'),
      uF: num('uF'),
      uR: num('uR'),
    };
    if (kind === 'OUT') {
      ev.out = num('out');
      ev.duration = num('duration');
    }
    events.push(ev);
  }
  return events;
}

/**
 * Parse the kernel's `purge` rows into a phase timeline. Two row shapes are
 * recognised (src/1d/src/twoFluid.cpp §5.162):
 *   state row  : purge t <t> phase <phase> x_wave <x> u_wave <u> cap .. cycles ..
 *   event row  : purge t <t> x_wave <x> head .. cap .. phase COMPLETE
 *                cycle <n> out .. duration .. period ..
 * Both are useful: state rows colour the timeline, COMPLETE rows mark a
 * completed cycle. Each row is tagged { t, phase, ..., cycle?, period? }.
 */
export function parsePurgePhases(out) {
  const phases = [];
  if (!out) return phases;
  for (const line of out.split('\n')) {
    const tok = line.replace(/^[^[]*\[chflow\]\s*/, '').trim().split(/\s+/);
    if (tok[0] !== 'purge') continue;
    const get = (key) => {
      const i = tok.indexOf(key);
      return i >= 0 ? tok[i + 1] : undefined;
    };
    const t = Number(get('t'));
    if (!Number.isFinite(t)) continue;
    const isEvent = get('cycle') !== undefined;
    const isState = get('u_wave') !== undefined;
    if (!isEvent && !isState) continue;
    phases.push({
      t,
      phase: get('phase'),
      xWave: Number(get('x_wave')),
      head: get('head') === undefined ? undefined : Number(get('head')),
      cap: get('cap') === undefined ? undefined : Number(get('cap')),
      cycles: get('cycles') === undefined ? undefined : Number(get('cycles')),
      out: get('out') === undefined ? undefined : Number(get('out')),
      cycle: isEvent ? Number(get('cycle')) : undefined,
      period: get('period') === undefined ? undefined : Number(get('period')),
      kind: isEvent ? 'event' : 'state',
    });
  }
  return phases;
}
