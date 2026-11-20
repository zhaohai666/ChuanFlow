/**
 * Pure helpers for the drag-and-drop network editor. No DOM and no Tauri
 * imports here, so the layout / hit-test / serialization logic is shared by
 * the WebView canvas and by `node --test`.
 *
 * The model is intentionally small: nodes carry a canvas position plus the
 * fields the case needs; pipes reference node ids and carry geometry. Nothing
 * here mutates its input — every operation returns a new model.
 */

let uid = 0;
function nextId(prefix) {
  uid += 1;
  return `${prefix}${uid}`;
}

/** Create an empty editor model. */
export function emptyModel() {
  return {
    nodes: [], // { id, name, x, y, pressure: number|null }
    pipes: [], // { id, name, from, to, length, diameter, roughness, elevation }
    sources: [], // { name, node, rate }
    fluid: { rho: 850, mu: 0.1 }, // matches network.cpp defaults
    time: { dt: 0.01, endTime: 1, reportEvery: 1 },
    selection: null, // { kind: 'node'|'pipe', id }
    pendingLink: null, // node id awaiting a second click to complete a pipe
  };
}

function uniqueName(nodes, base) {
  const names = new Set(nodes.map((n) => n.name));
  if (!names.has(base)) return base;
  let i = 2;
  while (names.has(`${base}${i}`)) i += 1;
  return `${base}${i}`;
}

/** Add a node at (x, y). Returns { model, id }. */
export function addNode(model, x, y, name) {
  const id = nextId('n');
  const label = name || uniqueName(model.nodes, 'N');
  return {
    model: {
      ...model,
      nodes: [...model.nodes, { id, name: label, x, y, pressure: null }],
    },
    id,
  };
}

/** Move a node to (x, y); clamped to the canvas box. */
export function moveNode(model, id, x, y, width = 640, height = 420) {
  return {
    ...model,
    nodes: model.nodes.map((n) =>
      n.id === id
        ? {
            ...n,
            x: Math.max(12, Math.min(width - 12, x)),
            y: Math.max(12, Math.min(height - 12, y)),
          }
        : n,
    ),
  };
}

/** Remove a node and every pipe/source attached to it. */
export function removeNode(model, id) {
  const node = model.nodes.find((n) => n.id === id);
  if (!node) return model;
  return {
    ...model,
    nodes: model.nodes.filter((n) => n.id !== id),
    pipes: model.pipes.filter((p) => p.from !== node.name && p.to !== node.name),
    sources: model.sources.filter((s) => s.node !== node.name),
    selection: null,
    pendingLink: null,
  };
}

/** Set a node's fixed pressure; null makes it a free node. */
export function setNodePressure(model, id, pressure) {
  return {
    ...model,
    nodes: model.nodes.map((n) =>
      n.id === id
        ? { ...n, pressure: pressure === null || pressure === '' ? null : Number(pressure) }
        : n,
    ),
  };
}

/**
 * Connect two node ids with a pipe. Defaults match `network.cpp::loadNetwork`
 * so the generated case is runnable. Self-links and duplicates are ignored.
 */
export function addPipe(model, fromId, toId) {
  if (!fromId || !toId || fromId === toId) return model;
  const from = model.nodes.find((n) => n.id === fromId);
  const to = model.nodes.find((n) => n.id === toId);
  if (!from || !to) return model;
  if (model.pipes.some((p) => p.from === from.name && p.to === to.name)) return model;
  const id = nextId('p');
  const base = `P${model.pipes.length + 1}`;
  const name = uniqueName(model.pipes, base);
  return {
    ...model,
    pipes: [
      ...model.pipes,
      {
        id,
        name,
        from: from.name,
        to: to.name,
        length: 1000,
        diameter: 0.05,
        roughness: 0,
        elevation: 0,
      },
    ],
    pendingLink: null,
  };
}

/** Remove a pipe by id. */
export function removePipe(model, id) {
  return { ...model, pipes: model.pipes.filter((p) => p.id !== id), selection: null };
}

/** Add / remove a source attached to a node name. */
export function addSource(model, nodeName, rate = 1e-6) {
  if (!model.nodes.some((n) => n.name === nodeName)) return model;
  const name = uniqueName(model.sources, 'S');
  return { ...model, sources: [...model.sources, { name, node: nodeName, rate: Number(rate) }] };
}

export function removeSource(model, name) {
  return { ...model, sources: model.sources.filter((s) => s.name !== name) };
}

/** Update transient time controls (subset of fields). */
export function setTime(model, patch) {
  return { ...model, time: { ...model.time, ...patch } };
}

/** Update the fluid block (subset of fields: rho, mu). */
export function setFluid(model, patch) {
  const base = model.fluid || { rho: 850, mu: 0.1 };
  return { ...model, fluid: { ...base, ...patch } };
}

/** Edit one pipe's geometry fields (length/diameter/roughness/elevation). */
export function setPipeFields(model, pipeId, patch) {
  return {
    ...model,
    pipes: model.pipes.map((p) => (p.id === pipeId ? { ...p, ...patch } : p)),
  };
}

/** Edit a source's rate by its dictionary name. */
export function updateSourceRate(model, name, rate) {
  return {
    ...model,
    sources: model.sources.map((s) =>
      s.name === name ? { ...s, rate: Number(rate) } : s,
    ),
  };
}

/**
 * Canvas hit test: nearest node whose centre is within `radius`, else null.
 * Pure so it can be tested without a DOM.
 */
export function hitTestNode(model, x, y, radius = 14) {
  let best = null;
  let bestDist = radius * radius;
  for (const n of model.nodes) {
    const dx = n.x - x;
    const dy = n.y - y;
    const d2 = dx * dx + dy * dy;
    if (d2 <= bestDist) {
      bestDist = d2;
      best = n;
    }
  }
  return best;
}

/** Format a number compactly, mirroring the engine's 12-significant style. */
function fmt(v) {
  if (!Number.isFinite(Number(v))) return '0';
  return String(Number(v));
}

/**
 * Serialize the editor model into a ChuanFlow dictionary. Emits the explicit
 * `node=` key on sources and the full pipe geometry, so the text is directly
 * runnable by the C++ solver (steady) and transient engine (when time is set).
 */
export function networkToCaseText(model) {
  const f = model.fluid || { rho: 850, mu: 0.1 };
  const lines = [
    '// generated by ChuanFlow GUI drag editor',
    'solver { type steadyNetwork; }',
    `fluid { rho ${fmt(f.rho)}; mu ${fmt(f.mu)}; }`,
    '',
    'nodes',
    '{',
  ];
  for (const n of model.nodes) {
    lines.push(
      n.pressure === null || n.pressure === undefined
        ? `    ${n.name} {}`
        : `    ${n.name} { pressure ${fmt(n.pressure)}; }`,
    );
  }
  lines.push('}', '', 'pipes', '{');
  for (const p of model.pipes) {
    // inline equipment (valve / pump curve) only when present (0 = absent),
    // so a plain pipe serializes byte-identically to the M5 output
    let fields =
      `from ${p.from}; to ${p.to}; length ${fmt(p.length)}; ` +
      `diameter ${fmt(p.diameter)}; roughness ${fmt(p.roughness)}; elevation ${fmt(p.elevation)};`;
    if (p.valve) fields += ` valve ${fmt(p.valve)};`;
    if (p.pumpHead) fields += ` pumpHead ${fmt(p.pumpHead)};`;
    if (p.pumpRate) fields += ` pumpRate ${fmt(p.pumpRate)};`;
    if (p.pumpShutoffHead) fields += ` pumpShutoffHead ${fmt(p.pumpShutoffHead)};`;
    if (p.pumpStages) fields += ` pumpStages ${fmt(p.pumpStages)};`;
    if (p.pumpSpeed) fields += ` pumpSpeed ${fmt(p.pumpSpeed)};`;
    if (p.pumpGasFraction) fields += ` pumpGasFraction ${fmt(p.pumpGasFraction)};`;
    if (p.pumpIntakePressure) fields += ` pumpIntakePressure ${fmt(p.pumpIntakePressure)};`;
    if (p.liftDomePressure) fields += ` liftDomePressure ${fmt(p.liftDomePressure)};`;
    if (p.liftAreaRatio) fields += ` liftAreaRatio ${fmt(p.liftAreaRatio)};`;
    if (p.liftSpringPressure) fields += ` liftSpringPressure ${fmt(p.liftSpringPressure)};`;
    if (p.liftDischargeCoeff) fields += ` liftDischargeCoeff ${fmt(p.liftDischargeCoeff)};`;
    if (p.liftKappa) fields += ` liftKappa ${fmt(p.liftKappa)};`;
    if (p.rodPlungerDiameter) fields += ` rodPlungerDiameter ${fmt(p.rodPlungerDiameter)};`;
    if (p.rodStrokeLength) fields += ` rodStrokeLength ${fmt(p.rodStrokeLength)};`;
    if (p.rodStrokesPerMinute) fields += ` rodStrokesPerMinute ${fmt(p.rodStrokesPerMinute)};`;
    if (p.rodDiameter) fields += ` rodDiameter ${fmt(p.rodDiameter)};`;
    if (p.rodLength) fields += ` rodLength ${fmt(p.rodLength)};`;
    if (p.rodYoungModulus) fields += ` rodYoungModulus ${fmt(p.rodYoungModulus)};`;
    if (p.rodDensity) fields += ` rodDensity ${fmt(p.rodDensity)};`;
    lines.push(`    ${p.name} { ${fields} }`);
  }
  lines.push('}', '', 'sources', '{');
  for (const s of model.sources) {
    lines.push(`    ${s.name} { node ${s.node}; rate ${fmt(s.rate)}; }`);
  }
  lines.push('}');

  const t = model.time || {};
  if (t.dt !== undefined && t.endTime !== undefined) {
    lines.push(
      '',
      `time { dt ${fmt(t.dt)}; endTime ${fmt(t.endTime)}; reportEvery ${Math.round(t.reportEvery || 1)}; }`,
    );
  }
  return lines.join('\n') + '\n';
}

/**
 * Arbitrary-plane clip: keep the cells whose centroid sits on the low side of
 * the plane with unit normal `normal`. `frac` (0..1) slides the plane between
 * the lowest and highest centroid projections, so axis-aligned clipping is the
 * special case normal = [1,0,0]. Degenerate normals keep every cell.
 * Pure so node --test can assert against the axis clip it replaces.
 */
export function clipKeepsByPlane(cells, points, normal, frac) {
  const len = Math.hypot(normal[0], normal[1], normal[2]);
  if (!Number.isFinite(len) || len < 1e-12) return cells.map(() => true);
  const nx = normal[0] / len;
  const ny = normal[1] / len;
  const nz = normal[2] / len;
  const proj = cells.map((vids) => {
    let cx = 0, cy = 0, cz = 0;
    for (const v of vids) {
      cx += points[v].x;
      cy += points[v].y;
      cz += points[v].z;
    }
    const m = vids.length || 1;
    return nx * (cx / m) + ny * (cy / m) + nz * (cz / m);
  });
  let lo = Infinity;
  let hi = -Infinity;
  for (const s of proj) {
    lo = Math.min(lo, s);
    hi = Math.max(hi, s);
  }
  const d = lo + frac * (hi - lo || 1);
  return proj.map((s) => s <= d);
}

/**
 * Magnitude-adaptive vector thinning: keep ceil(max/2) of the strongest cells
 * (ties broken by lower index) plus floor(max/2) evenly spaced indices from
 * the remainder, so hotspots stay dense without losing coverage. Deterministic
 * and capped at maxArrows. Pure so node --test can verify the cap and the
 * always-keep-strongest contract.
 */
export function adaptiveVectorIndices(mags, maxArrows) {
  const idx = mags.map((_, i) => i);
  if (!Number.isFinite(maxArrows) || maxArrows <= 0) return [];
  if (idx.length <= maxArrows) return idx;
  const hot = Math.ceil(maxArrows / 2);
  const rest = Math.floor(maxArrows / 2);
  const byMag = [...idx].sort((a, b) => (mags[b] - mags[a]) || (a - b));
  const keep = new Set(byMag.slice(0, hot));
  const remaining = idx.filter((i) => !keep.has(i));
  if (remaining.length > 0 && rest > 0) {
    const step = remaining.length / rest;
    for (let k = 0; k < rest; k += 1) {
      keep.add(remaining[Math.min(remaining.length - 1, Math.floor(k * step))]);
    }
  }
  return [...keep].sort((a, b) => a - b);
}
