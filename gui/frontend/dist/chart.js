/**
 * Pure helpers shared by the GUI and `node --test`. No DOM and no Tauri
 * imports here, so this file can be imported both by the WebView and by Node.
 */

/**
 * Lightweight validation of a simplified ChuanFlow dictionary. Mirrors the
 * Rust `parse_case_text` closely enough for live editor feedback. Supports
 * both `nodes\n{ ... }` and `nodes { A {...} }` layouts, and trailing `;` on
 * entries (the real dictionary style).
 */
export function validateCaseText(text) {
  const result = {
    valid: false,
    errors: [],
    nodes: [],
    pipes: [],
    sources: [],
  };
  if (typeof text !== 'string' || text.trim() === '') {
    result.errors.push('case 文本为空');
    return result;
  }

  function addEntry(section, name, body) {
    const kv = {};
    for (const stmt of body.split(';')) {
      const t = stmt.trim();
      if (!t) continue;
      const parts = t.split(/\s+/);
      kv[parts[0]] = parts.slice(1).join(' ');
    }
    if (section === 'nodes') {
      result.nodes.push({ name, pressure: kv.pressure ? Number(kv.pressure) : null });
    } else if (section === 'pipes') {
      result.pipes.push({ name, from: kv.from, to: kv.to });
    } else {
      // Real case style: `well1 { node A; rate ...; }`; simplified GUI style
      // omits `node` and uses the entry name as the node.
      result.sources.push({ node: kv.node || name, rate: Number(kv.rate || 0) });
    }
  }

  let section = null;
  for (const raw of text.split('\n')) {
    const line = raw.replace(/\/\/.*$/, '').replace(/#.*$/, '').trim();
    if (!line || line === '{') continue;
    if (line === '}') {
      section = null;
      continue;
    }

    const head = line.match(/^(nodes|pipes|sources)\b\s*\{?(.*)$/);
    if (head) {
      section = head[1];
      const rest = head[2].trim();
      if (rest.includes('}')) {
        // Single-line block: every `NAME { ... }` inside is one entry.
        for (const em of rest.matchAll(/([A-Za-z_]\w*)\s*\{([^}]*)\}\s*;?/g)) {
          addEntry(section, em[1], em[2]);
        }
        section = null;
      }
      continue;
    }
    if (!section) continue;

    const em = line.match(/^([A-Za-z_]\w*)\s*\{([^}]*)\}\s*;?\s*$/);
    if (!em) {
      result.errors.push(`无法解析条目: ${line}`);
      continue;
    }
    addEntry(section, em[1], em[2]);
  }

  if (result.nodes.length === 0) result.errors.push('缺少 nodes 块');
  if (result.pipes.length === 0) result.errors.push('缺少 pipes 块');
  const names = new Set(result.nodes.map((n) => n.name));
  for (const p of result.pipes) {
    if (!names.has(p.from)) result.errors.push(`管 ${p.name} 起点 ${p.from} 不存在`);
    if (!names.has(p.to)) result.errors.push(`管 ${p.name} 终点 ${p.to} 不存在`);
  }
  result.valid = result.errors.length === 0;
  return result;
}

/** Map values to bar heights (0..maxHeight), normalised by the largest |value|. */
export function toBarHeights(values, maxHeight = 100) {
  if (!Array.isArray(values) || values.length === 0) return [];
  const peak = Math.max(0, ...values.map((v) => Math.abs(Number(v) || 0)));
  if (peak === 0) return values.map(() => 0);
  return values.map((v) => (Math.abs(Number(v) || 0) / peak) * maxHeight);
}

/** Map values to an SVG `x,y x,y ...` polyline string inside width x height. */
export function toLinePoints(values, width, height, maxValue) {
  if (!Array.isArray(values) || values.length === 0) return '';
  const nums = values.map((v) => Number(v) || 0);
  const hi = maxValue !== undefined ? Number(maxValue) : Math.max(...nums);
  const lo = Math.min(0, ...nums);
  const span = hi - lo || 1;
  const step = nums.length === 1 ? 0 : width / (nums.length - 1);
  return nums
    .map((v, i) => {
      const x = i * step;
      const y = height - ((v - lo) / span) * height;
      return `${x.toFixed(2)},${y.toFixed(2)}`;
    })
    .join(' ');
}

/**
 * Map transient samples onto an SVG polyline using the real time `t` as the
 * horizontal axis (samples may be non-uniformly spaced). `valueOf(sample)`
 * selects the series (pressure / alpha / hydrate / flow for a given name).
 * Returns `{ points, min, max, tMin, tMax }` so the caller can draw axes.
 */
export function toTimeSeriesPoints(samples, valueOf, width, height) {
  const empty = { points: '', min: 0, max: 0, tMin: 0, tMax: 0 };
  if (!Array.isArray(samples) || samples.length === 0) return empty;
  const pairs = samples
    .map((s) => ({ t: Number(s.t), v: Number(valueOf(s)) }))
    .filter((p) => Number.isFinite(p.t) && Number.isFinite(p.v));
  if (pairs.length === 0) return empty;
  const tMin = Math.min(...pairs.map((p) => p.t));
  const tMax = Math.max(...pairs.map((p) => p.t));
  const min = Math.min(0, ...pairs.map((p) => p.v));
  const max = Math.max(...pairs.map((p) => p.v));
  const tSpan = tMax - tMin || 1;
  const vSpan = max - min || 1;
  const points = pairs
    .map((p) => {
      const x = ((p.t - tMin) / tSpan) * width;
      const y = height - ((p.v - min) / vSpan) * height;
      return `${x.toFixed(2)},${y.toFixed(2)}`;
    })
    .join(' ');
  return { points, min, max, tMin, tMax };
}
