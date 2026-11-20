/**
 * F18 Workbench-style parameter sets: named groups of parameter overrides
 * (+ description) that activate against the current model / case text and
 * batch-solve through the existing sweep scheduling. Serialization, override
 * merging and conflict detection are pure; node --test covers the merge
 * semantics and the exactness of activation. No DOM, no Tauri.
 */

import { networkToCaseText, applyWhatIfOverrides } from './editor.js';

export const PARAMSET_KIND = 'chuanflow-paramsets';
export const PARAMSET_VERSION = 1;

let psSeq = 0;
export function nextParamSetId() {
  psSeq += 1;
  return `ps${psSeq}`;
}

/** Identity of one override: parameter + target selects the same slot. */
export function overrideKey(o) {
  return `${o.param}@${o.target || ''}`;
}

/**
 * Merge two override lists: `patch` wins over `base` on identical keys
 * (same param + target); base order is preserved, new patch entries append.
 * Returns a fresh array; inputs are never mutated.
 */
export function mergeOverrides(base, patch) {
  const out = (base || []).map((o) => ({ ...o }));
  for (const p of patch || []) {
    const k = overrideKey(p);
    const i = out.findIndex((o) => overrideKey(o) === k);
    if (i >= 0) out[i] = { ...out[i], ...p };
    else out.push({ ...p });
  }
  return out;
}

/**
 * Conflict detection across parameter sets: every pair of sets that write
 * the same override key with different values. Returns
 * [{ a, b, key, valueA, valueB }] using set ids.
 */
export function findConflicts(sets) {
  const conflicts = [];
  const list = sets || [];
  for (let x = 0; x < list.length; x += 1) {
    for (let y = x + 1; y < list.length; y += 1) {
      const mapA = new Map((list[x].overrides || []).map((o) => [overrideKey(o), o]));
      for (const o of list[y].overrides || []) {
        const a = mapA.get(overrideKey(o));
        if (a && Number(a.value) !== Number(o.value)) {
          conflicts.push({
            a: list[x].id, b: list[y].id, key: overrideKey(o),
            valueA: Number(a.value), valueB: Number(o.value),
          });
        }
      }
    }
  }
  return conflicts;
}

/**
 * Serialize the workspace structure: { version, kind, activeId, sets }.
 * Sets keep id / name / description / engine / overrides only.
 */
export function paramSetsToJson(sets, activeId = null) {
  return JSON.stringify({
    version: PARAMSET_VERSION,
    kind: PARAMSET_KIND,
    activeId,
    sets: (sets || []).map((s) => ({
      id: s.id,
      name: s.name,
      description: s.description || '',
      engine: s.engine === 'transient' ? 'transient' : 'steady',
      overrides: (s.overrides || []).map((o) => ({
        param: o.param, target: o.target || '', value: Number(o.value),
      })),
    })),
  }, null, 2);
}

/**
 * Read the structure back; throws on a foreign kind, a bad version or a
 * non-array sets payload. Unknown fields on entries are dropped, values are
 * coerced to numbers.
 */
export function parseParamSets(text) {
  const data = typeof text === 'string' ? JSON.parse(text) : text;
  if (!data || data.kind !== PARAMSET_KIND) throw new Error('不是 ChuanFlow 参数集文件');
  if (data.version !== PARAMSET_VERSION) throw new Error(`参数集版本不支持: ${data.version}`);
  if (!Array.isArray(data.sets)) throw new Error('参数集缺少 sets 数组');
  return {
    activeId: data.activeId ?? null,
    sets: data.sets.map((s) => ({
      id: String(s.id),
      name: String(s.name ?? s.id),
      description: String(s.description ?? ''),
      engine: s.engine === 'transient' ? 'transient' : 'steady',
      overrides: (Array.isArray(s.overrides) ? s.overrides : []).map((o) => ({
        param: String(o.param), target: String(o.target ?? ''), value: Number(o.value),
      })),
    })),
  };
}

/**
 * Steady activation: apply a set's overrides to a COPY of the editor model
 * (same three parameters the steady sweep knows), throwing on an unknown
 * target so a stale set can never silently solve the wrong case.
 */
export function applySetToModel(model, set) {
  const m = JSON.parse(JSON.stringify(model));
  for (const o of set.overrides || []) {
    if (o.param === 'sourceRate') {
      const s = m.sources.find((x) => x.name === o.target);
      if (!s) throw new Error(`参数集 ${set.name}：源汇 ${o.target} 不存在`);
      s.rate = Number(o.value);
    } else if (o.param === 'pipeDiameter' || o.param === 'pipeLength') {
      const p = m.pipes.find((x) => x.name === o.target);
      if (!p) throw new Error(`参数集 ${set.name}：管段 ${o.target} 不存在`);
      p[o.param === 'pipeDiameter' ? 'diameter' : 'length'] = Number(o.value);
    } else {
      throw new Error(`参数集 ${set.name}：稳态不支持参数 ${o.param}（边界/阀门时序请切瞬态引擎）`);
    }
  }
  return m;
}

/**
 * Activation text: the exact case dictionary the set produces. Steady
 * serializes the overridden model; transient composes the what-if edits on
 * the raw dictionary (byte-exact single-number replacements).
 */
export function activateSet(set, model, caseText) {
  if (set.engine === 'transient') {
    return (set.overrides || []).reduce(
      (text, o) => applyWhatIfOverrides(text, o.param, o.target, o.value),
      caseText,
    );
  }
  return networkToCaseText(applySetToModel(model, set));
}
