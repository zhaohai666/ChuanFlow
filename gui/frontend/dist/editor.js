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

// ---------------------------------------------------------------------------
// Steady equipment field registry (P6-P8 fields + E4 families). One registry
// drives serialize, extractCaseDetails and the model rebuild so the three can
// never drift. Key names are the kernel's own (network.cpp pipe keys for
// ESP/lift/rod/compressor, pointElementsApp.cpp keys for the four standalone
// point elements, prefixed with the family tag to avoid cross-device clashes).
// ---------------------------------------------------------------------------
export const DEVICE_FIELDS = [
  // ESP advanced intake semantics (steady kernel reads both keys)
  'pumpGasFraction', 'pumpIntakePressure',
  // gas-lift valve (§5.89)
  'liftDomePressure', 'liftAreaRatio', 'liftSpringPressure',
  'liftDischargeCoeff', 'liftKappa',
  // rod pump (§5.90)
  'rodPlungerDiameter', 'rodStrokeLength', 'rodStrokesPerMinute', 'rodDiameter',
  'rodLength', 'rodYoungModulus', 'rodDensity',
  // E4: jet pump (§5.140)
  'jetNozzleArea', 'jetThroatArea', 'jetThroatLength', 'jetFrictionFactor',
  'jetMotiveRate', 'jetSuctionPressure', 'jetSuctionRate', 'jetDischargePressure',
  // E4: progressive cavity pump (§5.140)
  'pcpCavitiesPerRev', 'pcpCavityVolume', 'pcpRpm', 'pcpSlip',
  'pcpPressureRise', 'pcpDischargePressure', 'pcpIntakePressure',
  // E4: plunger gas lift (§5.140)
  'pglGasVolumePerCycle', 'pglReferencePressure', 'pglGasDensityAtReference',
  'pglCyclePeriod', 'pglOpenWindow', 'pglReleaseDepth',
  'pglCasingSurfacePressure', 'pglTubingSurfacePressure', 'pglLiquidDensity',
  'pglDomePressure', 'pglAreaRatio', 'pglSpringPressure', 'pglPortDiameter',
  'pglDischargeCoeff', 'pglKappa',
  // E4: separator (§5.115 / pointElementsApp.cpp)
  'sepFeedRate', 'sepPressure', 'sepTemperature', 'sepEfficiency', 'sepDP',
  'sepGasGravity', 'sepOilDensity', 'sepGor', 'sepWaterCut',
  // E4: compressor — real network.cpp pipe keys (§5.156), no prefix needed
  'pressureRatio', 'shutoffRatio', 'ratedRate', 'polytropicEfficiency',
  'kappa', 'specificGasConstant', 'compressibility',
];

/** Base (M5) pipe keys that also carry equipment but predate the registry. */
export const BASIC_DEVICE_FIELDS = [
  'valve', 'pumpHead', 'pumpRate', 'pumpShutoffHead', 'pumpStages', 'pumpSpeed',
];

/** Every equipment key the GUI can round-trip through a case dictionary. */
export const ALL_DEVICE_FIELDS = [...BASIC_DEVICE_FIELDS, ...DEVICE_FIELDS];

/**
 * E4 device families: `switch` is the > 0 trigger, `fields` the editable
 * properties, `block` the top-level point-element dictionary block consumed
 * by `chflow-1d-point-elements` (null for compressor, whose keys are native
 * network.cpp pipe keys solved by `chflow-1d-steady`).
 */
export const DEVICE_FAMILIES = [
  {
    id: 'jetPump', label: '射流泵', switch: 'jetNozzleArea', block: 'jetPumps',
    fields: ['jetNozzleArea', 'jetThroatArea', 'jetThroatLength', 'jetFrictionFactor',
      'jetMotiveRate', 'jetSuctionPressure', 'jetSuctionRate', 'jetDischargePressure'],
    defaults: {
      jetNozzleArea: 1e-4, jetThroatArea: 4e-4, jetThroatLength: 0.2,
      jetFrictionFactor: 0.02, jetMotiveRate: 2e-3, jetSuctionPressure: 2e6,
      jetSuctionRate: 0, jetDischargePressure: 2.05e6,
    },
  },
  {
    id: 'pcp', label: '螺杆泵 PCP', switch: 'pcpCavityVolume', block: 'pcpPumps',
    fields: ['pcpCavitiesPerRev', 'pcpCavityVolume', 'pcpRpm', 'pcpSlip',
      'pcpPressureRise', 'pcpDischargePressure', 'pcpIntakePressure'],
    defaults: {
      pcpCavitiesPerRev: 2, pcpCavityVolume: 1.2e-5, pcpRpm: 120, pcpSlip: 2e-12,
      pcpPressureRise: 8e6, pcpDischargePressure: 0, pcpIntakePressure: 0,
    },
  },
  {
    id: 'plungerGasLift', label: '柱塞气举', switch: 'pglGasVolumePerCycle',
    block: 'plungerGasLifts',
    fields: ['pglGasVolumePerCycle', 'pglReferencePressure',
      'pglGasDensityAtReference', 'pglCyclePeriod', 'pglOpenWindow',
      'pglReleaseDepth', 'pglCasingSurfacePressure', 'pglTubingSurfacePressure',
      'pglLiquidDensity', 'pglDomePressure', 'pglAreaRatio', 'pglSpringPressure',
      'pglPortDiameter', 'pglDischargeCoeff', 'pglKappa'],
    defaults: {
      pglGasVolumePerCycle: 0.5, pglReferencePressure: 1e5,
      pglGasDensityAtReference: 1.2, pglCyclePeriod: 600, pglOpenWindow: 30,
      pglReleaseDepth: 1500, pglCasingSurfacePressure: 8e6,
      pglTubingSurfacePressure: 4e6, pglLiquidDensity: 900, pglDomePressure: 9e6,
      pglAreaRatio: 0.1, pglSpringPressure: 0, pglPortDiameter: 0.01,
      pglDischargeCoeff: 0.665, pglKappa: 1.3,
    },
  },
  {
    id: 'separator', label: '分离器', switch: 'sepFeedRate', block: 'separators',
    fields: ['sepFeedRate', 'sepPressure', 'sepTemperature', 'sepEfficiency',
      'sepDP', 'sepGasGravity', 'sepOilDensity', 'sepGor', 'sepWaterCut'],
    defaults: {
      sepFeedRate: 0.01, sepPressure: 2e6, sepTemperature: 340, sepEfficiency: 1,
      sepDP: 0, sepGasGravity: 0.75, sepOilDensity: 850, sepGor: 180,
      sepWaterCut: 0.05,
    },
  },
  {
    id: 'compressor', label: '压缩机', switch: 'pressureRatio', block: null,
    fields: ['pressureRatio', 'shutoffRatio', 'ratedRate',
      'polytropicEfficiency', 'kappa', 'specificGasConstant', 'compressibility'],
    defaults: {
      pressureRatio: 2, shutoffRatio: 2.4, ratedRate: 0.05,
      polytropicEfficiency: 0.9, kappa: 1.3, specificGasConstant: 518.28,
      compressibility: 0.95,
    },
  },
];

/**
 * The five E4 families (excludes the three P6-P8 devices already wired in).
 * Exported so the toolbar, property panel, removal action and block emitter
 * all iterate the same list.
 */
export const E4_DEVICE_FAMILIES = DEVICE_FAMILIES;

/** Pull every equipment key out of an `extractCaseDetails` pipe record. */
export function pipeDeviceFieldsFrom(geo = {}) {
  const out = {};
  for (const k of ALL_DEVICE_FIELDS) out[k] = geo[k] ?? 0;
  return out;
}

/**
 * Emit the top-level point-element blocks for the four devices whose physics
 * lives in `chflow-1d-point-elements` (jet pump / PCP / PGL / separator).
 * Keys are the kernel's exact names from `pointElementsApp.cpp`; each block
 * entry is named after its source pipe, so one GUI pipe maps to one element.
 * No device fields => no blocks (plain-pipe byte identity preserved).
 */
export function deviceBlocksFromModel(model) {
  const rho = fmt((model.fluid || {}).rho ?? 850);
  const lines = [];
  const active = (key) => model.pipes.filter((p) => Number(p[key]) > 0);
  const emit = (block, entries) => {
    if (!entries.length) return;
    lines.push('', block, '{', ...entries, '}');
  };

  emit('jetPumps', active('jetNozzleArea').map((p) => {
    let b = `rho ${rho}; nozzleArea ${fmt(p.jetNozzleArea)}; ` +
      `throatArea ${fmt(p.jetThroatArea)}; throatLength ${fmt(p.jetThroatLength)}; ` +
      `frictionFactor ${fmt(p.jetFrictionFactor)}; motiveRate ${fmt(p.jetMotiveRate)}; ` +
      `dischargePressure ${fmt(p.jetDischargePressure)};`;
    b += Number(p.jetSuctionRate) > 0
      ? ` suctionRate ${fmt(p.jetSuctionRate)};`
      : ` suctionPressure ${fmt(p.jetSuctionPressure)};`;
    return `    ${p.name} { ${b} }`;
  }));

  emit('pcpPumps', active('pcpCavityVolume').map((p) => {
    let b = `rho ${rho}; cavitiesPerRev ${fmt(p.pcpCavitiesPerRev)}; ` +
      `cavityVolume ${fmt(p.pcpCavityVolume)}; rpm ${fmt(p.pcpRpm)}; ` +
      `slip ${fmt(p.pcpSlip)};`;
    b += Number(p.pcpDischargePressure) > 0 && Number(p.pcpIntakePressure) > 0
      ? ` dischargePressure ${fmt(p.pcpDischargePressure)}; intakePressure ${fmt(p.pcpIntakePressure)};`
      : ` pressureRise ${fmt(p.pcpPressureRise)};`;
    return `    ${p.name} { ${b} }`;
  }));

  emit('plungerGasLifts', active('pglGasVolumePerCycle').map((p) => {
    const b = `gasVolumePerCycle ${fmt(p.pglGasVolumePerCycle)}; ` +
      `referencePressure ${fmt(p.pglReferencePressure)}; ` +
      `gasDensityAtReference ${fmt(p.pglGasDensityAtReference)}; ` +
      `cyclePeriod ${fmt(p.pglCyclePeriod)}; openWindow ${fmt(p.pglOpenWindow)}; ` +
      `releaseDepth ${fmt(p.pglReleaseDepth)}; ` +
      `casingSurfacePressure ${fmt(p.pglCasingSurfacePressure)}; ` +
      `tubingSurfacePressure ${fmt(p.pglTubingSurfacePressure)}; ` +
      `liquidDensity ${fmt(p.pglLiquidDensity)}; ` +
      `liftDomePressure ${fmt(p.pglDomePressure)}; liftAreaRatio ${fmt(p.pglAreaRatio)}; ` +
      `liftSpringPressure ${fmt(p.pglSpringPressure)}; ` +
      `liftPortDiameter ${fmt(p.pglPortDiameter)}; ` +
      `liftDischargeCoeff ${fmt(p.pglDischargeCoeff)}; liftKappa ${fmt(p.pglKappa)};`;
    return `    ${p.name} { ${b} }`;
  }));

  emit('separators', active('sepFeedRate').map((p) =>
    `    ${p.name} { feedRate ${fmt(p.sepFeedRate)}; pressure ${fmt(p.sepPressure)}; ` +
    `temperature ${fmt(p.sepTemperature)}; efficiency ${fmt(p.sepEfficiency)}; ` +
    `dP ${fmt(p.sepDP)}; ` +
    `fluid { blackoil { gasGravity ${fmt(p.sepGasGravity)}; ` +
    `oilDensity ${fmt(p.sepOilDensity)}; gor ${fmt(p.sepGor)}; ` +
    `waterCut ${fmt(p.sepWaterCut)}; } } }`,
  ));

  return lines.length ? lines.join('\n') + '\n' : '';
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
    for (const key of DEVICE_FIELDS) {
      if (p[key]) fields += ` ${key} ${fmt(p[key])};`;
    }
    lines.push(`    ${p.name} { ${fields} }`);
  }
  lines.push('}', '', 'sources', '{');
  for (const s of model.sources) {
    lines.push(`    ${s.name} { node ${s.node}; rate ${fmt(s.rate)}; }`);
  }
  lines.push('}');

  // E4 point-element blocks (jet pump / PCP / PGL / separator). Absent when
  // no such device is present, so plain-pipe output stays byte-identical.
  const blocks = deviceBlocksFromModel(model);
  if (blocks) lines.push(blocks.trimEnd());

  const t = model.time || {};
  if (t.dt !== undefined && t.endTime !== undefined) {
    lines.push(
      '',
      `time { dt ${fmt(t.dt)}; endTime ${fmt(t.endTime)}; reportEvery ${Math.round(t.reportEvery || 1)}; }`,
    );
  }
  return lines.join('\n') + '\n';
}

// Default compliance for a free node when the drag editor runs a transient:
// the editor model carries only a fixed pressure (or none), never a tank volume
// or bulk modulus, so a free node serialized as `{}` has no storage term and the
// transient pressure matrix is singular. Supplying a declared default compliance
// makes the default model runnable; the user can override it in the raw dict.
export const TRANSIENT_NODE_DEFAULTS = {
  volume: 1,        // m³
  bulkModulus: 1e7, // Pa
  pressure: 1e5,    // Pa, free-node initial pressure
  dt: 0.1,          // s
  endTime: 10,      // s
};

/**
 * Serialize the editor model as a TRANSIENT network case. Unlike the steady
 * serialization, a fixed-pressure node becomes `fixedPressure` (the transient
 * Dirichlet port) and a free node gets a declared default volume + bulkModulus
 * so its pressure row carries a storage term (no singular matrix). A time{}
 * block is always emitted (model.time or the defaults). Point-element device
 * blocks are omitted: the transient row does not consume them.
 */
export function modelToTransientCaseText(model) {
  const f = model.fluid || { rho: 850, mu: 0.1 };
  const D = TRANSIENT_NODE_DEFAULTS;
  const lines = [
    '// generated by ChuanFlow GUI drag editor (transient)',
    'solver { type transientNetwork; }',
    `fluid { rho ${fmt(f.rho)}; mu ${fmt(f.mu)}; }`,
    '',
    'nodes',
    '{',
  ];
  for (const n of model.nodes) {
    lines.push(
      n.pressure === null || n.pressure === undefined
        ? `    ${n.name} { pressure ${fmt(D.pressure)}; volume ${fmt(D.volume)}; `
          + `bulkModulus ${fmt(D.bulkModulus)}; liquidFraction 1; }`
        : `    ${n.name} { fixedPressure ${fmt(n.pressure)}; }`,
    );
  }
  lines.push('}', '', 'pipes', '{');
  for (const p of model.pipes) {
    let fields =
      `from ${p.from}; to ${p.to}; length ${fmt(p.length)}; ` +
      `diameter ${fmt(p.diameter)}; roughness ${fmt(p.roughness)}; elevation ${fmt(p.elevation)};`;
    if (p.valve) fields += ` valve ${fmt(p.valve)};`;
    lines.push(`    ${p.name} { ${fields} }`);
  }
  lines.push('}', '', 'sources', '{');
  for (const s of model.sources) {
    lines.push(`    ${s.name} { node ${s.node}; rate ${fmt(s.rate)}; }`);
  }
  lines.push('}');
  const t = model.time || {};
  const dt = t.dt !== undefined ? t.dt : D.dt;
  const endTime = t.endTime !== undefined ? t.endTime : D.endTime;
  lines.push('', `time { dt ${fmt(dt)}; endTime ${fmt(endTime)}; reportEvery ${Math.round(t.reportEvery || 1)}; }`);
  return lines.join('\n') + '\n';
}

/**
 * Transient what-if parameters (U10 / §5.137). Each entry maps a UI label to
 * a precise dictionary location. `section` decides which target list the GUI
 * offers; `global` parameters need no target at all.
 */
export const TRANSIENT_WHAT_IF_PARAMS = [
  { value: 'valveCloseStart', label: '关阀开始 (s)', section: 'pipes', block: 'valveClose', key: 'start' },
  { value: 'valveCloseEnd', label: '关阀结束 (s)', section: 'pipes', block: 'valveClose', key: 'end' },
  { value: 'closedK', label: '关闭阻力 K', section: 'pipes', block: 'valveClose', key: 'closedK' },
  { value: 'valveOpenStart', label: '开阀开始 (s)', section: 'pipes', block: 'valveOpen', key: 'start' },
  { value: 'valveOpenEnd', label: '开阀结束 (s)', section: 'pipes', block: 'valveOpen', key: 'end' },
  { value: 'scheduleOpenRate', label: '源注入率 (m³/s)', section: 'sources', block: 'schedule', key: 'scheduleOpenRate' },
  { value: 'scheduleCloseStart', label: '源关闭开始 (s)', section: 'sources', block: 'schedule', key: 'scheduleCloseStart' },
  { value: 'scheduleCloseEnd', label: '源关闭结束 (s)', section: 'sources', block: 'schedule', key: 'scheduleCloseEnd' },
  { value: 'scheduleOpenStart', label: '源重启开始 (s)', section: 'sources', block: 'schedule', key: 'scheduleOpenStart' },
  { value: 'scheduleOpenEnd', label: '源重启结束 (s)', section: 'sources', block: 'schedule', key: 'scheduleOpenEnd' },
  { value: 'endTime', label: '终止时间 (s)', section: 'time', block: 'time', key: 'endTime' },
  { value: 'dt', label: '时间步 (s)', section: 'time', block: 'time', key: 'dt' },
  { value: 'reportEvery', label: '报告间隔 (步)', section: 'time', block: 'time', key: 'reportEvery' },
];

const NUMBER_RE = '[-+]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][-+]?\\d+)?';

/** Locate `{ ... }` matching the opening brace at `open`, skipping // comments. */
function braceEnd(text, open) {
  let depth = 0;
  for (let i = open; i < text.length; i++) {
    if (text[i] === '/' && text[i + 1] === '/') {
      const eol = text.indexOf('\n', i);
      i = eol < 0 ? text.length : eol;
      continue;
    }
    if (text[i] === '{') depth++;
    else if (text[i] === '}' && --depth === 0) return i;
  }
  return -1;
}

function blockBody(text, name, from = 0) {
  const m = new RegExp(`\\b${name}\\s*\\{`).exec(text.slice(from));
  if (!m) return null;
  const open = from + m.index + m[0].length - 1;
  const end = braceEnd(text, open);
  return end < 0 ? null : { body: text.slice(open + 1, end), start: open + 1, end };
}

function namedEntry(text, section, name) {
  const sec = blockBody(text, section);
  if (!sec) return null;
  const m = new RegExp(`\\b${name}\\s*\\{`).exec(sec.body);
  if (!m) return null;
  const open = sec.start + m.index + m[0].length - 1;
  const end = braceEnd(text, open);
  return end < 0 ? null : { body: text.slice(open + 1, end), start: open + 1, end };
}

function replaceNumber(body, key, value) {
  const re = new RegExp(`(\\b${key}\\s+)(${NUMBER_RE})(\\s*;)`);
  if (!re.test(body)) return null;
  return body.replace(re, `$1${value}$3`);
}

/**
 * Apply one transient what-if override to a case-dictionary TEXT, changing
 * exactly one number and leaving every other byte untouched. The transient
 * sweeps operate on the raw dictionary rather than the editor model because
 * tank volume, bulk modulus and valve/source schedules cannot be faithfully
 * reconstructed from the model — this is the honest counterpart of the
 * steady sweep's model serialization. Throws on an unsupported parameter,
 * missing target entry, non-finite value or a key that is not present.
 */
export function applyWhatIfOverrides(caseText, param, target, value) {
  const spec = TRANSIENT_WHAT_IF_PARAMS.find((p) => p.value === param);
  if (!spec) throw new Error(`unsupported what-if parameter: ${param}`);
  const v = Number(value);
  if (!Number.isFinite(v)) throw new Error(`what-if value is not finite: ${value}`);
  if (spec.section === 'time') {
    const timeBlock = blockBody(caseText, 'time');
    if (!timeBlock) throw new Error('case has no time{} block');
    const changed = replaceNumber(timeBlock.body, spec.key, v);
    if (changed === null) throw new Error(`time{} has no ${spec.key} key`);
    return caseText.slice(0, timeBlock.start) + changed + caseText.slice(timeBlock.end);
  }
  if (!target) throw new Error(`what-if ${param} needs a target entry`);
  const entry = namedEntry(caseText, spec.section, target);
  if (!entry) throw new Error(`case has no ${spec.section} entry named ${target}`);
  let block = entry;
  const isStartEnd = spec.key.endsWith('Start') || spec.key.endsWith('End');
  if (spec.section === 'sources') {
    // sources.<name>.schedule { openRate; close { start; end; }; open { ... }; }
    const sched = blockBody(entry.body, 'schedule');
    if (!sched) throw new Error(`source ${target} has no schedule{} block`);
    block = { body: sched.body, start: entry.start + sched.start, end: entry.start + sched.end };
    if (isStartEnd) {
      const inner = spec.key.startsWith('scheduleClose') ? 'close' : 'open';
      const sub = blockBody(block.body, inner);
      if (!sub) throw new Error(`source ${target} has no schedule.${inner}{} block`);
      block = { body: sub.body, start: block.start + sub.start, end: block.start + sub.end };
    }
  } else if (spec.section === 'pipes') {
    // pipes.<name>.valveClose { start; end; closedK; } / valveOpen { start; end; }
    const sub = blockBody(entry.body, spec.block);
    if (!sub) throw new Error(`pipe ${target} has no ${spec.block}{} block`);
    block = { body: sub.body, start: entry.start + sub.start, end: entry.start + sub.end };
  }
  const lookupKey = isStartEnd ? (spec.key.endsWith('Start') ? 'start' : 'end')
    : spec.key === 'scheduleOpenRate' ? 'openRate' : spec.key;
  const changed = replaceNumber(block.body, lookupKey, v);
  if (changed === null) throw new Error(`entry ${target} has no ${lookupKey} key`);
  return caseText.slice(0, block.start) + changed + caseText.slice(block.end);
}

/**
 * Apply a scenario's whole override list to a transient dictionary, one
 * what-if at a time (each edit re-scans the previous result, so multi-key
 * scenarios compose exactly like the kernel reads them). Throws with the
 * offending override index.
 */
export function applyScenarioOverrides(caseText, overrides) {
  let text = caseText;
  (overrides || []).forEach((o, i) => {
    try {
      text = applyWhatIfOverrides(text, o.param, o.target, o.value);
    } catch (e) {
      throw new Error(`工况第 ${i + 1} 项 (${o.param}${o.target ? ` @${o.target}` : ''}=${o.value}) 失败: ${e.message || e}`);
    }
  });
  return text;
}

/**
 * Replace a nodal case's whole top-level ipr{} block with a linear inflow of
 * productivity index `j` at reservoir pressure `pRes`, leaving every other
 * byte (fluid/well/solve/wellhead) untouched. Each such run's working point
 * lies on the wellbore's VLP curve - this is the nodal sampler's text edit.
 * Throws when the case has no balanced ipr{} block.
 */
export function iprLinearVariant(caseText, pRes, j) {
  const m = /\bipr\s*\{/.exec(caseText);
  if (!m) throw new Error('nodal 字典缺少 ipr{} 块');
  const open = m.index + m[0].length - 1;
  const close = braceEnd(caseText, open);
  if (close < 0) throw new Error('ipr{} 花括号不平衡');
  const block = `ipr {\n` +
    `    form              linear;\n` +
    `    reservoirPressure ${pRes};\n` +
    `    productivityIndex ${j};\n` +
    `}`;
  return caseText.slice(0, m.index) + block + caseText.slice(close + 1);
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

/**
 * E1: point a case's `output{}` block at a fixed VTU path so the GUI can
 * render the field after the run. If an output block exists, only the `vtk`
 * word inside it is rewritten (other keys — e.g. a field list — survive);
 * otherwise a minimal output block is appended. Comment lines are ignored
 * during the scan, mirroring the engine's dictionary semantics closely
 * enough for the app-authored cases this is used on.
 */
export function upsertOutputVtk(caseText, vtuPath) {
  // find the last uncommented `output` block opener: the name and the `{`
  // may sit on the same line or on consecutive lines (dict convention)
  const lines = caseText.split('\n');
  let openIdx = -1;
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i].replace(/\/\/.*$/, '').replace(/#.*$/, '');
    if (!/\boutput\b/.test(line)) continue;
    if (line.includes('{')) {
      openIdx = i;
    } else {
      for (let j = i + 1; j < lines.length; j++) {
        const nxt = lines[j].replace(/\/\/.*$/, '').replace(/#.*$/, '').trim();
        if (!nxt) continue;
        if (nxt.includes('{')) openIdx = j;
        break;
      }
    }
  }
  if (openIdx >= 0) {
    // bracket-match to the block's closing brace
    let depth = 0;
    let endIdx = -1;
    for (let i = openIdx; i < lines.length; i++) {
      const code = lines[i].replace(/\/\/.*$/, '').replace(/#.*$/, '');
      for (const ch of code) {
        if (ch === '{') depth++;
        else if (ch === '}') depth--;
      }
      if (depth <= 0 && i >= openIdx && code.includes('}')) {
        endIdx = i;
        break;
      }
    }
    if (endIdx < 0) return caseText + `\noutput\n{\n    vtk ${vtuPath};\n}\n`;
    const block = lines.slice(openIdx, endIdx + 1);
    let replaced = false;
    const rewritten = block.map((line) => {
      if (replaced) return line;
      const code = line.replace(/\/\/.*$/, '').replace(/#.*$/, '');
      const m = code.match(/(\bvtk\s+)(\S+?)(\s*;)/);
      if (m) {
        replaced = true;
        return line.replace(/(\bvtk\s+)(\S+?)(\s*;)/, `$1${vtuPath}$3`);
      }
      return line;
    });
    if (replaced) {
      lines.splice(openIdx, endIdx - openIdx + 1, ...rewritten);
      return lines.join('\n');
    }
    // output block without a vtk word: inject one after the opening line
    lines.splice(openIdx + 1, 0, `    vtk ${vtuPath};`);
    return lines.join('\n');
  }
  return caseText + `\noutput\n{\n    vtk ${vtuPath};\n}\n`;
}

/**
 * E2: replace or append a top-level dictionary block. If a block with the
 * same name exists (dict convention: the name and its `{` may share a line
 * or sit on consecutive lines), the whole block is replaced; otherwise the
 * block is appended at the end. Pure so node --test can verify the three
 * paths (replace / append / name-with-inline-brace).
 */
export function upsertDictBlock(caseText, name, innerLines) {
  const lines = caseText.split('\n');
  const blockRe = new RegExp(`\\b${name}\\b`);
  let openIdx = -1;
  for (let i = 0; i < lines.length; i++) {
    const code = lines[i].replace(/\/\/.*$/, '').replace(/#.*$/, '');
    if (!blockRe.test(code)) continue;
    if (code.includes('{')) {
      openIdx = i;
      break;
    }
    // the name may carry trailing words (e.g. `node X;`) before its own
    // brace on the next line — only treat it as a block opener when the
    // next non-empty line is exactly a `{` (top-level block convention)
    for (let j = i + 1; j < lines.length; j++) {
      const nxt = lines[j].replace(/\/\/.*$/, '').replace(/#.*$/, '').trim();
      if (!nxt) continue;
      if (nxt === '{') openIdx = j;
      break;
    }
    if (openIdx >= 0) break;
  }
  const newBlock = [`${name}`, '{', ...innerLines.map((s) => `    ${s}`), '}'];
  if (openIdx < 0) {
    const sep = caseText.endsWith('\n') ? '' : '\n';
    return caseText + sep + newBlock.join('\n') + '\n';
  }
  // bracket-match to the block's closing brace (start depth from the opener)
  let depth = 0;
  let endIdx = openIdx;
  for (let i = openIdx; i < lines.length; i++) {
    const code = lines[i].replace(/\/\/.*$/, '').replace(/#.*$/, '');
    for (const ch of code) {
      if (ch === '{') depth++;
      else if (ch === '}') depth--;
    }
    if (depth <= 0) {
      endIdx = i;
      break;
    }
  }
  lines.splice(openIdx - (lines[openIdx].trim() === '{' ? 1 : 0), endIdx - openIdx + 1 + (lines[openIdx].trim() === '{' ? 1 : 0), ...newBlock);
  return lines.join('\n');
}

/**
 * E2: remove a top-level dictionary block entirely (the inverse of
 * upsertDictBlock). Used when a flow-assurance declaration is switched off:
 * the kernel treats an absent block as byte-identical to the pre-block path,
 * so removing it (rather than emitting `enable false`) keeps the case clean.
 * Anchors the name at line start (after whitespace) to avoid matching a key
 * inside another block. Pure; returning the input unchanged when absent.
 */
export function removeDictBlock(caseText, name) {
  const lines = caseText.split('\n');
  const nameRe = new RegExp(`^\\s*${name}\\b`);
  for (let i = 0; i < lines.length; i++) {
    const code = lines[i].replace(/\/\/.*$/, '').replace(/#.*$/, '');
    if (!nameRe.test(code)) continue;
    const braceOnSame = code.includes('{');
    let openIdx = -1;
    if (braceOnSame) openIdx = i;
    else {
      for (let j = i + 1; j < lines.length; j++) {
        const nxt = lines[j].replace(/\/\/.*$/, '').replace(/#.*$/, '').trim();
        if (!nxt) continue;
        if (nxt === '{') openIdx = j;
        break;
      }
    }
    if (openIdx < 0) continue;
    let depth = 0;
    let endIdx = openIdx;
    for (let k = openIdx; k < lines.length; k++) {
      const c = lines[k].replace(/\/\/.*$/, '').replace(/#.*$/, '');
      for (const ch of c) {
        if (ch === '{') depth++;
        else if (ch === '}') depth--;
      }
      if (depth <= 0) { endIdx = k; break; }
    }
    const start = braceOnSame ? i : i;
    lines.splice(start, endIdx - start + 1);
    // drop a single blank line left where the block stood
    if (lines[start] && !lines[start].trim() && lines[start - 1] && lines[start - 1].trim()) {
      lines.splice(start, 1);
    }
    return lines.join('\n');
  }
  return caseText;
}

/**
 * E4: pure SVG glyph command generator. Returns `{ className, title, commands }`
 * where each command is `{ type, ...attrs }`; the canvas turns commands into
 * real SVG nodes and `node --test` can assert the geometry without a DOM.
 * The three P6-P8 families keep their exact pre-E4 shapes and ordering.
 */
export function deviceGlyphCommands(p, mx, my, ax, ay, bx, by) {
  const el = (type, attrs) => ({ type, ...attrs });
  const bits = (label, pairs) => {
    const out = [`${label} ${p.name}:`];
    for (const [k, v] of pairs) if (v) out.push(`${k}=${v}`);
    return out.join(' ');
  };

  if (p.liftDomePressure) {
    return {
      className: 'dev-gaslift',
      title: bits('气举阀', [['dome', p.liftDomePressure], ['Ar', p.liftAreaRatio],
        ['spring', p.liftSpringPressure], ['Cd', p.liftDischargeCoeff],
        ['kappa', p.liftKappa]]),
      commands: [
        el('path', {
          d: `M ${mx - 8},${my + 6} Q ${mx - 8},${my - 9} ${mx},${my - 9} ` +
             `Q ${mx + 8},${my - 9} ${mx + 8},${my + 6} Z`,
        }),
        el('circle', { cx: mx, cy: my + 8, r: 2.5 }),
      ],
    };
  }
  if (p.rodPlungerDiameter) {
    return {
      className: 'dev-rodpump',
      title: bits('有杆泵', [['bore', p.rodPlungerDiameter], ['stroke', p.rodStrokeLength],
        ['spm', p.rodStrokesPerMinute], ['rodD', p.rodDiameter], ['rodL', p.rodLength]]),
      commands: [
        el('rect', { x: mx - 10, y: my - 5, width: 20, height: 10 }),
        el('path', {
          d: `M ${mx},${my - 5} L ${mx},${my - 12} L ${mx + 7},${my - 12}`,
        }),
      ],
    };
  }
  if (p.valve) {
    return {
      className: 'dev-valve',
      title: `节流阀 ${p.name}: valve K=${p.valve}`,
      commands: [
        el('polygon', { points: `${mx - 7},${my - 7} ${mx - 7},${my + 7} ${mx},${my}` }),
        el('polygon', { points: `${mx + 7},${my - 7} ${mx + 7},${my + 7} ${mx},${my}` }),
        el('line', { x1: mx, y1: my - 7, x2: mx, y2: my - 12 }),
      ],
    };
  }
  if (p.pumpHead || p.pumpRate) {
    const dx = bx - ax;
    const dy = by - ay;
    const len = Math.hypot(dx, dy) || 1;
    const ux = dx / len;
    const uy = dy / len;
    return {
      className: 'dev-pump',
      title: bits(`泵`, [['pumpHead', p.pumpHead], ['pumpRate', p.pumpRate],
        ['shutoff', p.pumpShutoffHead], ['stages', p.pumpStages],
        ['speed', p.pumpSpeed], ['gasFrac', p.pumpGasFraction],
        ['intake', p.pumpIntakePressure]]),
      commands: [
        el('circle', { cx: mx, cy: my, r: 9 }),
        el('polygon', {
          points: `${mx + ux * 6},${my + uy * 6} ` +
            `${mx - ux * 4 - uy * 4},${my - uy * 4 + ux * 4} ` +
            `${mx - ux * 4 + uy * 4},${my - uy * 4 - ux * 4}`,
        }),
      ],
    };
  }
  // ---- E4: five new device families -------------------------------------
  if (p.jetNozzleArea) {
    return {
      className: 'dev-jetpump',
      title: bits('射流泵', [['An', p.jetNozzleArea], ['At', p.jetThroatArea],
        ['motive', p.jetMotiveRate], ['pd', p.jetDischargePressure]]),
      commands: [
        el('polygon', { points: `${mx - 13},${my - 5} ${mx - 13},${my + 5} ${mx - 3},${my}` }),
        el('rect', { x: mx - 3, y: my - 4, width: 11, height: 8 }),
        el('path', {
          d: `M ${mx + 8},${my - 4} L ${mx + 14},${my - 8} ` +
             `M ${mx + 8},${my + 4} L ${mx + 14},${my + 8}`,
        }),
      ],
    };
  }
  if (p.pcpCavityVolume) {
    return {
      className: 'dev-pcp',
      title: bits('螺杆泵 PCP', [['cav', p.pcpCavityVolume], ['rpm', p.pcpRpm],
        ['slip', p.pcpSlip], ['dp', p.pcpPressureRise]]),
      commands: [
        el('rect', { x: mx - 12, y: my - 7, width: 24, height: 14, rx: 6 }),
        el('path', {
          d: `M ${mx - 8},${my} Q ${mx - 4},${my - 7} ${mx},${my} ` +
             `T ${mx + 8},${my}`,
        }),
      ],
    };
  }
  if (p.pglGasVolumePerCycle) {
    return {
      className: 'dev-pgl',
      title: bits('柱塞气举', [['vol', p.pglGasVolumePerCycle],
        ['T', p.pglCyclePeriod], ['w', p.pglOpenWindow],
        ['depth', p.pglReleaseDepth], ['dome', p.pglDomePressure]]),
      commands: [
        el('rect', { x: mx - 6, y: my - 11, width: 12, height: 22 }),
        el('rect', { x: mx - 5, y: my, width: 10, height: 5 }),
        el('circle', { cx: mx, cy: my + 8, r: 2 }),
        el('circle', { cx: mx - 10, cy: my + 5, r: 1.6 }),
        el('circle', { cx: mx + 10, cy: my - 4, r: 1.6 }),
      ],
    };
  }
  if (p.sepFeedRate) {
    return {
      className: 'dev-separator',
      title: bits('分离器', [['feed', p.sepFeedRate], ['p', p.sepPressure],
        ['T', p.sepTemperature], ['eff', p.sepEfficiency]]),
      commands: [
        el('rect', { x: mx - 12, y: my - 8, width: 24, height: 16, rx: 8 }),
        el('line', { x1: mx - 18, y1: my, x2: mx - 12, y2: my }),
        el('line', { x1: mx + 12, y1: my - 4, x2: mx + 18, y2: my - 4 }),
        el('line', { x1: mx + 12, y1: my + 4, x2: mx + 18, y2: my + 4 }),
        el('line', { x1: mx, y1: my - 8, x2: mx, y2: my + 2 }),
      ],
    };
  }
  if (p.pressureRatio) {
    return {
      className: 'dev-compressor',
      title: bits('压缩机', [['PR', p.pressureRatio], ['PR0', p.shutoffRatio],
        ['rate', p.ratedRate], ['eta', p.polytropicEfficiency],
        ['kappa', p.kappa]]),
      commands: [
        el('circle', { cx: mx, cy: my, r: 9 }),
        el('polygon', {
          points: `${mx},${my} ${mx + 5},${my - 6} ${mx + 7},${my + 2}`,
        }),
        el('path', { d: `M ${mx + 9},${my} L ${mx + 15},${my} L ${mx + 12},${my - 3}` }),
      ],
    };
  }
  return null;
}
