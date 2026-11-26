/**
 * ChuanFlow professional workspace (M7): four-column CAE layout on the
 * existing Tauri commands. Everything wired here is real: the steady /
 * transient / coupled engines run as stoppable background jobs, the 3D field
 * view renders the kernel's own ASCII .vtu output, and the unit system /
 * theme / project save-load are display-layer state over the SI dictionary.
 */

import { invoke } from './tauri-api.js';
import {
  emptyModel, addNode, moveNode, removeNode, setNodePressure, addPipe,
  removePipe, addSource, removeSource, setTime, hitTestNode,
  networkToCaseText, setPipeFields, updateSourceRate, setFluid,
  clipKeepsByPlane, adaptiveVectorIndices, upsertOutputVtk,
  applyWhatIfOverrides, TRANSIENT_WHAT_IF_PARAMS, iprLinearVariant,
  applyScenarioOverrides, deviceGlyphCommands, pipeDeviceFieldsFrom,
  E4_DEVICE_FAMILIES, ALL_DEVICE_FIELDS,
} from './editor.js';
import {
  validateCaseText, toBarHeights, toTimeSeriesPoints,
} from './chart.js';
import {
  UNIT_SYSTEMS, fmtDisplay, toDisplay, fromDisplay, unitLabel, parseInput,
} from './units.js';
import {
  buildProjectTree, parseCoupledOutput, parseVtu, cellCentroid, fieldRange,
  rainbowColor, classifyLogLine, parseSteadyOutput, parseTransientOutput,
  extractCaseDetails, parseTwoFluidFrames, traceStreamlines,
  evaluateReports, reportsToCsv,
  parseNodalIpr, iprCurvePoints, vlpSampleJs,
  vorticityAndQ, isoSurfacePatches, fieldsToCsv, planeProfileCsv,
  slugSeriesKeys, slugBands,
  buildMonitorTree, monitorSeriesValues, monitorEnvBand, hydrateSeriesNote,
  compareSteadyRuns,
  scenarioRowsTransient, scenarioRowsSteady,
  pvtFlashCaseText, parsePvtOutput, pvtCurvesByT, parsePvtCaseDict, PVT_COMPONENTS,
} from './workspace.js';
import { wellControlPanel, riskPanel } from './safety.js';
import {
  templateTarget, isPartialFidelity, isSteadyComparable, THREE_D_FAMILY, COUPLED_FAMILY,
} from './templates.js';
import {
  REGIME_NAMES, REGIME_COLORS, mapPoint, regimeHit, regimeMapDrawCommands,
  beggsBrillPattern,
} from './regimemap.js';
import {
  tornadoData, tornadoCommands,
} from './tornado.js';
import {
  marchingSquares, rasterizeSliceField, contourColor, profileBands,
} from './contour.js';
import {
  nextParamSetId, paramSetsToJson, parseParamSets, activateSet, findConflicts,
  applySetToModel, mergeOverrides,
} from './paramsets.js';
import {
  createMonitorState, feedMonitor, monitorDrawCommands,
} from './monitor.js';
import {
  FA_BLOCKS, FA_BLOCK_NAMES, parseFlowAssurance, applyFlowAssurance,
  parseSlugEvents, parsePurgePhases,
} from './flowAssurance.js';
import {
  parseSchedules, pipeNames, applySchedules, esdPreset, timelineSvg,
} from './scheduleEditor.js';

// ---------------------------------------------------------------------------
// state
// ---------------------------------------------------------------------------

const DEFAULT_EXE = '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-steady';
// Audit fix (P1-2): the point-element engine (chflow-1d-point-elements) is the
// only consumer of the E4 device top-level blocks (jetPumps/pcpPumps/
// plungerGasLifts/separators — the steady network app ignores unknown blocks,
// so device physics silently dropped out of a plain ▶稳态 run). Lives in the
// same build directory as the steady binary; E9 persistence covers it.
const DEFAULT_POINT_ELEMENTS_EXE = DEFAULT_EXE.replace('chflow-1d-steady', 'chflow-1d-point-elements');
const POINT_ELEMENT_BLOCK_RE = /^\s*(jetPumps|pcpPumps|plungerGasLifts|separators)\s*\{/m;
const DEFAULT_TRANSIENT_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-transient';
const DEFAULT_COUPLED_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/coupling/apps/chflow-coupled-1d3d';
const DEFAULT_TWOFLUID_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-twofluid';
const DEFAULT_NODAL_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-nodal';
const DEFAULT_PVT_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/pvt/apps/chflow-pvt';
const DEFAULT_3D_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/3d/apps/chflow-piso-foam';
// E1 (§5.181): the 3D engine family's developer-tree defaults. The Rust
// `three_d_engines` resolve beats these at startup, exactly like the 8
// pre-existing engine keys.
const DEFAULT_3D_FAMILY_DIR = '/Users/haizai/project/ChuanFlow/build/src/3d/apps';
// Coupling engine family developer-tree defaults; the Rust `coupled_engines`
// resolve beats these at startup, same contract as the 3D family.
const DEFAULT_COUPLED_FAMILY_DIR = '/Users/haizai/project/ChuanFlow/build/src/coupling/apps';
const THREE_D_VTU = '/tmp/chflow_gui_3d_field.vtu';
const DEFAULT_COUPLED_CASE =
  '/Users/haizai/project/ChuanFlow/tests/cases/coupled1d3d/case.dict';
const COUPLED_VTU = '/tmp/chflow_gui_coupled_block.vtu';

const DEFAULT_CASE = `// ChuanFlow 简化稳态管网
solver { type steadyNetwork; }
fluid  { rho 850; mu 0.1; }

nodes
{
    A { pressure 3.0e4; }
    B { pressure 2.0e4; }
    E {};
}

pipes
{
    P1 { from A; to E; length 20; diameter 0.02; yieldStress 3.0; }
    P2 { from E; to B; length 20; diameter 0.02; yieldStress 3.0; }
}

sources
{
    E { rate 1e-6; }
}

time { dt 0.01; endTime 1; reportEvery 1; }
`;

const state = {
  model: emptyModel(),
  units: UNIT_SYSTEMS.SI,
  theme: 'dark',
  sel: null, // { kind, ref } from the tree / canvas
  steady: null,
  transient: null,
  monitors: [],
  compare: null,
  scenarios: [],
  scenarioSel: null,
  scenarioResults: null,
  pvt: { comps: [{ name: 'methane', z: 0.5 }, { name: 'nPentane', z: 0.5 }], parsed: null },
  coupled: null,
  coupledEdge: null,
  vtu: null,
  colorMode: 'none',
  fieldMode: 'p',
  view3d: { yaw: 0, pitch: 0, clipAxis: 'none', clipFrac: 0.5, clipNormal: [1, 0, 0], vectors: false, streamlines: false, iso: { on: false, frac: 0.5 }, contours: { on: false, levels: 5 } },
  streamlineCache: null,
  vortCache: null,
  reports: {
    defs: [
      { name: 'pAvg', field: 'p', op: 'volumeAverage' },
      { name: 'fluxX', field: 'Ux', op: 'flux', normal: [1, 0, 0], frac: 1 },
    ],
    history: [], key: null,
  },
  slug: {
    frames: [], idx: 0, playing: false, timer: null, series: ['p'], bands: false,
    // E2: slugTracking events + purgeWave phase timeline from the kernel log.
    events: [], purge: [], tMin: 0, tMax: 1,
  },
  flowAssurance: null,
  regime: { diameter: 0.05, vsl: 0.1, vsg: 0.3 },
  tornado: { rows: [], running: false, data: null },
  paramSets: [],
  activeParamSet: null,
  monitorStream: createMonitorState(1500),
  nodal: null,
  sweep: { rows: [], running: false, out: '', engine: 'steady' },
  loadedTemplate: null,
  job: { id: null, label: null, running: false, offset: 0, full: '' },
  collapsedLeft: false,
  collapsedRight: false,
};

const ENGINES_STORE_KEY = 'chflow.engines.v1';
// E1 (§5.181): the 3D engine family rides the same persisted store — one
// key per family member, aligned with the pre-existing 8-key style.
const THREE_D_FAMILY_KEYS = THREE_D_FAMILY.map((e) => e.key);
const COUPLED_FAMILY_KEYS = COUPLED_FAMILY.map((e) => e.key);
const ENGINES_KEYS = [
  'steady', 'pointElements', 'transient', 'coupled', 'coupledCase', 'twofluid', 'nodal', 'pvt', 'threeD',
  ...THREE_D_FAMILY_KEYS.filter((k) => k !== 'threeD'),
  ...COUPLED_FAMILY_KEYS.filter((k) => k !== 'coupled'),
];
const ENGINE_INPUT_IDS = {
  steady: 'exePath',
  pointElements: 'pointElementsExePath',
  transient: 'transientExePath',
  coupled: 'coupledExePath',
  coupledCase: 'coupledCasePath',
  twofluid: 'twofluidExePath',
  nodal: 'nodalExePath',
  pvt: 'pvtExePath',
  threeD: 'threeDExePath',
};

const engines = {
  steady: DEFAULT_EXE,
  pointElements: DEFAULT_POINT_ELEMENTS_EXE,
  transient: DEFAULT_TRANSIENT_EXE,
  coupled: DEFAULT_COUPLED_EXE,
  coupledCase: DEFAULT_COUPLED_CASE,
  twofluid: DEFAULT_TWOFLUID_EXE,
  nodal: DEFAULT_NODAL_EXE,
  pvt: DEFAULT_PVT_EXE,
  threeD: DEFAULT_3D_EXE,
};
// E1 (§5.181): per-member developer-tree defaults; threeD (PISO) is declared
// above so the original constant keeps its single meaning.
for (const e of THREE_D_FAMILY) {
  if (e.key !== 'threeD') engines[e.key] = `${DEFAULT_3D_FAMILY_DIR}/${e.bin}`;
}
// coupling family per-member defaults; the steady single-edge `coupled` key
// keeps the original DEFAULT_COUPLED_EXE above.
for (const e of COUPLED_FAMILY) {
  if (e.key !== 'coupled') engines[e.key] = `${DEFAULT_COUPLED_FAMILY_DIR}/${e.bin}`;
}
// which family member the next ▶ 3D run submits to (set by template load /
// the dictionary-page family selector; defaults to PISO)
let active3DEngine = 'threeD';
// which coupling member the next ▶ 耦合 run submits to (set by template load;
// defaults to the steady single-edge coupled engine)
let activeCoupledEngine = 'coupled';

// E9: restore engine paths persisted from a previous session (settings beat
// hardcoded developer defaults; the backend resolve beats both below).
(function restoreEngines() {
  try {
    const raw = localStorage.getItem(ENGINES_STORE_KEY);
    if (!raw) return;
    const saved = JSON.parse(raw);
    for (const k of ENGINES_KEYS) {
      if (typeof saved[k] === 'string' && saved[k]) engines[k] = saved[k];
    }
  } catch (_) {
    /* corrupted store -> fall back to defaults */
  }
})();

function persistEngines() {
  try {
    const out = {};
    for (const k of ENGINES_KEYS) out[k] = engines[k];
    localStorage.setItem(ENGINES_STORE_KEY, JSON.stringify(out));
  } catch (_) {
    /* storage unavailable (e.g. file:// without quota) -> session-only */
  }
}

function reflectEngines() {
  for (const k of ENGINES_KEYS) {
    const el = document.getElementById(ENGINE_INPUT_IDS[k]);
    if (el) el.value = engines[k];
  }
}

let uidTimer = Date.now() % 100000;
function tmpPath(base) {
  uidTimer += 1;
  return `/tmp/chflow_gui_${base}_${uidTimer}.dict`;
}

// ---------------------------------------------------------------------------
// dom shortcuts
// ---------------------------------------------------------------------------

const $ = (id) => document.getElementById(id);
const svg = $('networkCanvas');
const treeEl = $('projectTree');
const propsEl = $('propPanel');
const logEl = $('logPanel');
const convEl = $('convergence');

// ---------------------------------------------------------------------------
// logging
// ---------------------------------------------------------------------------

function appendLog(text) {
  const lines = text.split('\n').filter((l) => l.trim() !== '');
  for (const line of lines) {
    const div = document.createElement('div');
    div.className = `log-line ${classifyLogLine(line)}`;
    div.textContent = line;
    logEl.appendChild(div);
  }
  while (logEl.childElementCount > 800) logEl.removeChild(logEl.firstChild);
  logEl.scrollTop = logEl.scrollHeight;
}

function setStatus(msg, cls = '') {
  const el = $('status');
  el.textContent = msg;
  el.className = `status ${cls}`;
}

function setJobChip(running, label = '') {
  const chip = $('jobChip');
  chip.textContent = running ? `运行中 · ${label}` : '空闲';
  chip.className = `chip ${running ? 'running' : 'idle'}`;
  $('tbStop').disabled = !running;
}

// ---------------------------------------------------------------------------
// tabs
// ---------------------------------------------------------------------------

function switchTab(tabId) {
  for (const b of $('mainTabs').querySelectorAll('.tab')) {
    b.classList.toggle('active', b.dataset.tab === tabId);
  }
  for (const sec of document.querySelectorAll('.tab-page')) {
    sec.classList.toggle('active', sec.id === `tab-${tabId}`);
  }
  if (tabId === 'field') renderField();
  if (tabId === 'slug' && state.slug.frames.length) renderSlugFrame();
  // E2: the timeline canvas has zero clientWidth while hidden; redraw on entry.
  if (tabId === 'slug') renderFlowAssuranceTimeline();
  if (tabId === 'regime') renderRegimeMap();
  if (tabId === 'paramsets') renderParamSets();
}

$('mainTabs').addEventListener('click', (e) => {
  const b = e.target.closest('.tab');
  if (b) switchTab(b.dataset.tab);
});

$('bottomTabs').addEventListener('click', (e) => {
  const b = e.target.closest('.tab');
  if (!b) return;
  for (const t of $('bottomTabs').querySelectorAll('.tab')) {
    t.classList.toggle('active', t === b);
  }
  for (const page of document.querySelectorAll('.btab-page')) {
    page.classList.toggle('active', page.id === `btab-${b.dataset.btab}`);
  }
  if (b.dataset.btab === 'curves') renderCharts();
  if (b.dataset.btab === 'monitor') renderMonitorPanel();
});

// ---------------------------------------------------------------------------
// project tree
// ---------------------------------------------------------------------------

function renderTree() {
  const tree = buildProjectTree(state.model, {
    hasSteady: !!state.steady,
    hasTransient: !!(state.transient && state.transient.samples.length),
    hasCoupled: !!state.coupled,
    coupledEdge: state.coupledEdge,
    scenarios: state.scenarios.map((s) => ({
      id: s.id, name: s.name, engine: s.engine, nOverrides: s.overrides.length,
    })),
  });
  treeEl.innerHTML = '';

  const renderLevel = (items, container) => {
    for (const item of items) {
      const row = document.createElement('div');
      row.className = 'row' + (isItemSelected(item) ? ' selected' : '');
      row.dataset.treeId = item.id;
      if (item.children) {
        row.innerHTML = `<span class="caret">▾</span><span>${item.label}</span>`;
        const kids = document.createElement('div');
        kids.className = 'group-children';
        renderLevel(item.children, kids);
        row.addEventListener('click', () => kids.classList.toggle('collapsed'));
        container.appendChild(row);
        container.appendChild(kids);
      } else {
        const mark =
          item.ok === true ? ' <span class="leaf-ok">✓</span>'
          : item.ok === false ? ''
          : '';
        row.innerHTML = `<span class="caret">·</span><span>${item.label}${mark}</span>`;
        row.addEventListener('click', () => selectTreeItem(item));
        container.appendChild(row);
      }
    }
  };
  renderLevel(tree, treeEl);
}

function isItemSelected(item) {
  if (!state.sel) return false;
  if (state.sel.treeId) return state.sel.treeId === item.id;
  if (item.kind === 'node' && state.sel.kind === 'node') return state.sel.ref === item.ref;
  if (item.kind === 'pipe' && state.sel.kind === 'pipe') return state.sel.ref === item.ref;
  if (item.kind === 'source' && state.sel.kind === 'source') return state.sel.ref === item.ref;
  return false;
}

function selectTreeItem(item) {
  state.sel = { kind: item.kind, ref: item.ref, treeId: item.id };
  if (item.kind === 'scenario') {
    state.scenarioSel = item.ref;
    switchTab('scenarios');
    renderTree();
    renderScenarios();
    return;
  }
  if (item.kind === 'node' || item.kind === 'pipe') {
    state.model.selection = { kind: item.kind, id: item.ref };
  }
  switchTab('canvas');
  renderTree();
  renderProps();
  renderCanvas();
}

// ---------------------------------------------------------------------------
// property panel
// ---------------------------------------------------------------------------

function numInput(label, value, family, onApply) {
  const wrap = document.createElement('label');
  wrap.className = 'pfield';
  const unit = family ? unitLabel(family, state.units) : '';
  wrap.innerHTML = `<span>${label}${unit ? ` (${unit})` : ''}</span>`;
  const input = document.createElement('input');
  input.type = 'text';
  input.spellcheck = false;
  input.value = family && Number.isFinite(value)
    ? String(Number(toDisplay(value, family, state.units).toPrecision(8)))
    : String(value);
  input.addEventListener('change', () => {
    const v = family ? parseInput(input.value, family, state.units) : Number(input.value);
    if (!Number.isFinite(v)) {
      setStatus(`无效数值：${label} = ${input.value}`, 'error');
      return;
    }
    onApply(v);
  });
  wrap.appendChild(input);
  return wrap;
}

// E4 device field labels + display unit families (null = raw SI number).
const DEVICE_FIELD_META = {
  jetNozzleArea: ['喷嘴面积 An (m²)', null],
  jetThroatArea: ['喉道面积 At (m²)', null],
  jetThroatLength: ['喉道长度 Lt (m)', 'length'],
  jetFrictionFactor: ['喉道摩擦系数 f', null],
  jetMotiveRate: ['动力液流量 (m³/s)', 'rate'],
  jetSuctionPressure: ['吸入压力 (Pa)', 'pressure'],
  jetSuctionRate: ['吸入流量 (m³/s，0 = 按压力解)', 'rate'],
  jetDischargePressure: ['排出压力 (Pa)', 'pressure'],
  pcpCavitiesPerRev: ['每转密封腔数', null],
  pcpCavityVolume: ['单腔体积 (m³)', null],
  pcpRpm: ['转速 (rpm)', null],
  pcpSlip: ['滑移系数 (m³/s/Pa)', null],
  pcpPressureRise: ['压差 (Pa)', 'pressure'],
  pcpDischargePressure: ['排出压力 (Pa，与吸入同时填 = 端口模式)', 'pressure'],
  pcpIntakePressure: ['吸入压力 (Pa，与排出同时填 = 端口模式)', 'pressure'],
  pglGasVolumePerCycle: ['每周期气量 Vg (m³)', null],
  pglReferencePressure: ['参考压力 (Pa)', 'pressure'],
  pglGasDensityAtReference: ['参考气密度 (kg/m³)', 'density'],
  pglCyclePeriod: ['周期 T (s)', null],
  pglOpenWindow: ['开阀窗口 w (s)', null],
  pglReleaseDepth: ['释放深度 (m)', 'length'],
  pglCasingSurfacePressure: ['套管井口压力 (Pa)', 'pressure'],
  pglTubingSurfacePressure: ['油管井口压力 (Pa)', 'pressure'],
  pglLiquidDensity: ['液相密度 (kg/m³)', 'density'],
  pglDomePressure: ['阀充压 Pd (Pa)', 'pressure'],
  pglAreaRatio: ['阀面积比 R', null],
  pglSpringPressure: ['弹簧压力 St (Pa)', 'pressure'],
  pglPortDiameter: ['阀孔径 (m)', 'length'],
  pglDischargeCoeff: ['流量系数 Cd', null],
  pglKappa: ['绝热指数 k', null],
  sepFeedRate: ['进料量 (m³/s)', 'rate'],
  sepPressure: ['分离压力 (Pa)', 'pressure'],
  sepTemperature: ['分离温度 (K)', 'temperature'],
  sepEfficiency: ['液相回收效率 (0..1)', null],
  sepDP: ['级间压降 dP (Pa)', 'pressure'],
  sepGasGravity: ['气体相对密度', null],
  sepOilDensity: ['标况油密度 (kg/m³)', 'density'],
  sepGor: ['气油比 GOR (sm³/sm³)', null],
  sepWaterCut: ['含水率 (0..1)', null],
  pressureRatio: ['压比 PR (rated)', null],
  shutoffRatio: ['关死压比 PR0', null],
  ratedRate: ['额定流量 (m³/s)', 'rate'],
  polytropicEfficiency: ['多变效率 (0..1)', null],
  kappa: ['绝热指数 k (>1)', null],
  specificGasConstant: ['气体常数 R (J/kg·K)', null],
  compressibility: ['压缩因子 Z', null],
};

function renderProps() {
  propsEl.innerHTML = '';
  const sel = state.sel;
  const head = document.createElement('h3');

  if (!sel) {
    head.textContent = '属性';
    propsEl.appendChild(head);
    const note = document.createElement('p');
    note.className = 'note';
    note.textContent = '在项目树或画布中选择节点 / 管段 / 源 / 流体，参数在此编辑。';
    propsEl.appendChild(note);
    return;
  }

  if (sel.kind === 'node') {
    const n = state.model.nodes.find((x) => x.id === sel.ref);
    if (!n) return;
    head.textContent = `节点 ${n.name}`;
    propsEl.appendChild(head);
    propsEl.appendChild(numInput('边界压力（留空 = 自由节点）', n.pressure ?? '', 'pressure', (v) => {
      state.model = setNodePressure(state.model, n.id, v);
      afterModelChange();
    }));
    const prow = document.createElement('div');
    prow.className = 'prow';
    const del = document.createElement('button');
    del.textContent = '删除节点';
    del.addEventListener('click', () => {
      state.model = removeNode(state.model, n.id);
      state.sel = null;
      afterModelChange();
    });
    prow.appendChild(del);
    propsEl.appendChild(prow);
    return;
  }

  if (sel.kind === 'pipe') {
    const p = state.model.pipes.find((x) => x.id === sel.ref);
    if (!p) return;
    head.textContent = `管段 ${p.name}`;
    propsEl.appendChild(head);
    const kv = document.createElement('div');
    kv.className = 'kv';
    kv.innerHTML = `<span>连接</span><span>${p.from} → ${p.to}</span>`;
    propsEl.appendChild(kv);
    propsEl.appendChild(numInput('管长', p.length, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { length: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('内径', p.diameter, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { diameter: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('粗糙度', p.roughness, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { roughness: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('高程差', p.elevation, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { elevation: v });
      afterModelChange();
    }));

    // inline equipment (network.cpp pipe keys), 0 = absent
    const devHead = document.createElement('p');
    devHead.className = 'note';
    devHead.textContent = '内嵌设备（0 = 无）：';
    propsEl.appendChild(devHead);
    propsEl.appendChild(numInput('阀门 valve K', p.valve ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { valve: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('泵扬程 pumpHead (m)', p.pumpHead ?? 0, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpHead: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('泵工况流量 pumpRate (m³/s)', p.pumpRate ?? 0, 'rate', (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpRate: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('关死扬程 pumpShutoffHead (m)', p.pumpShutoffHead ?? 0, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpShutoffHead: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('泵级数 pumpStages', p.pumpStages ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpStages: Math.max(0, Math.round(v)) });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('转速 pumpSpeed (fraction)', p.pumpSpeed ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpSpeed: v });
      afterModelChange();
    }));
    // ESP advanced intake semantics (steady kernel reads both keys)
    propsEl.appendChild(numInput('气锁体积分数 pumpGasFraction', p.pumpGasFraction ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpGasFraction: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('吸入口压力 pumpIntakePressure (Pa)', p.pumpIntakePressure ?? 0, 'pressure', (v) => {
      state.model = setPipeFields(state.model, p.id, { pumpIntakePressure: v });
      afterModelChange();
    }));
    // gas-lift valve (§5.89, steady-only): dome pressure is the switch
    propsEl.appendChild(numInput('气举阀充压 liftDomePressure (Pa)', p.liftDomePressure ?? 0, 'pressure', (v) => {
      state.model = setPipeFields(state.model, p.id, { liftDomePressure: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('气举阀面积比 liftAreaRatio', p.liftAreaRatio ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { liftAreaRatio: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('弹簧压力 liftSpringPressure (Pa)', p.liftSpringPressure ?? 0, 'pressure', (v) => {
      state.model = setPipeFields(state.model, p.id, { liftSpringPressure: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('流量系数 liftDischargeCoeff', p.liftDischargeCoeff ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { liftDischargeCoeff: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('流道指数 liftKappa', p.liftKappa ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { liftKappa: v });
      afterModelChange();
    }));
    // rod pump (§5.90, steady-only; transient refuses it)
    propsEl.appendChild(numInput('柱塞直径 rodPlungerDiameter (m)', p.rodPlungerDiameter ?? 0, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { rodPlungerDiameter: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('冲程长度 rodStrokeLength (m)', p.rodStrokeLength ?? 0, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { rodStrokeLength: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('冲次 rodStrokesPerMinute', p.rodStrokesPerMinute ?? 0, null, (v) => {
      state.model = setPipeFields(state.model, p.id, { rodStrokesPerMinute: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('抽油杆直径 rodDiameter (m)', p.rodDiameter ?? 0, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { rodDiameter: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('抽油杆长度 rodLength (m)', p.rodLength ?? 0, 'length', (v) => {
      state.model = setPipeFields(state.model, p.id, { rodLength: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('杆材杨氏模量 rodYoungModulus (Pa)', p.rodYoungModulus ?? 0, 'pressure', (v) => {
      state.model = setPipeFields(state.model, p.id, { rodYoungModulus: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('杆材密度 rodDensity (kg/m³)', p.rodDensity ?? 0, 'density', (v) => {
      state.model = setPipeFields(state.model, p.id, { rodDensity: v });
      afterModelChange();
    }));
    // E4: five new device families (§5.176). Grouped by family; the switch
    // field (> 0) is listed first in each group.
    for (const fam of E4_DEVICE_FAMILIES) {
      const grp = document.createElement('p');
      grp.className = 'note';
      grp.textContent = `${fam.label}（${fam.switch} > 0 生效）：`;
      propsEl.appendChild(grp);
      for (const key of fam.fields) {
        const [label, unitFamily] = DEVICE_FIELD_META[key] || [key, null];
        propsEl.appendChild(numInput(label, p[key] ?? 0, unitFamily, (v) => {
          state.model = setPipeFields(state.model, p.id, { [key]: v });
          afterModelChange();
        }));
      }
    }
    if (state.steady) {
      const r = state.steady.pipes.find((x) => x.name === p.name);
      if (r) {
        const res = document.createElement('div');
        res.className = 'kv';
        res.innerHTML = `<span>结果</span><span>Q ${fmtDisplay(r.q, 'rate', state.units)} · dp ${fmtDisplay(r.dp, 'pressure', state.units)}</span>`;
        propsEl.appendChild(res);
      }
    }
    const prow = document.createElement('div');
    prow.className = 'prow';
    const del = document.createElement('button');
    del.textContent = '删除管段';
    del.addEventListener('click', () => {
      state.model = removePipe(state.model, p.id);
      state.sel = null;
      afterModelChange();
    });
    prow.appendChild(del);
    propsEl.appendChild(prow);
    return;
  }

  if (sel.kind === 'source') {
    const s = state.model.sources.find((x) => x.name === sel.ref);
    if (!s) return;
    head.textContent = `源/汇 ${s.name}`;
    propsEl.appendChild(head);
    const kv = document.createElement('div');
    kv.className = 'kv';
    kv.innerHTML = `<span>节点</span><span>${s.node}</span>`;
    propsEl.appendChild(kv);
    propsEl.appendChild(numInput('注入率（负 = 抽出）', s.rate, 'rate', (v) => {
      state.model = updateSourceRate(state.model, s.name, v);
      afterModelChange();
    }));
    return;
  }

  if (sel.kind === 'fluid') {
    head.textContent = '流体物性';
    propsEl.appendChild(head);
    const f = state.model.fluid || { rho: 850, mu: 0.1 };
    propsEl.appendChild(numInput('密度 ρ', f.rho, 'density', (v) => {
      state.model = setFluid(state.model, { rho: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('黏度 μ', f.mu, 'viscosity', (v) => {
      state.model = setFluid(state.model, { mu: v });
      afterModelChange();
    }));
    return;
  }

  if (sel.kind === 'time') {
    head.textContent = '瞬态控制';
    propsEl.appendChild(head);
    const t = state.model.time;
    propsEl.appendChild(numInput('时间步 dt (s)', t.dt, null, (v) => {
      state.model = setTime(state.model, { dt: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('结束时间 endTime (s)', t.endTime, null, (v) => {
      state.model = setTime(state.model, { endTime: v });
      afterModelChange();
    }));
    propsEl.appendChild(numInput('reportEvery', t.reportEvery, null, (v) => {
      state.model = setTime(state.model, { reportEvery: Math.max(1, Math.round(v)) });
      afterModelChange();
    }));
    return;
  }

  if (sel.kind === 'coupled') {
    head.textContent = '耦合接口 (1D ↔ 3D)';
    propsEl.appendChild(head);
    const note = document.createElement('p');
    note.className = 'note';
    note.textContent = '耦合边把 1D 管段的阻力交给 3D SIMPLE 块求解，交界面按 rate-for-drop 交换直到收敛。';
    propsEl.appendChild(note);
    if (state.coupled) {
      const ex = state.coupled.solved;
      if (ex) {
        const kv = document.createElement('div');
        kv.className = 'kv';
        kv.innerHTML =
          `<span>边 ${ex.edge}</span><span>q ${fmtDisplay(ex.q, 'rate', state.units)}</span>`;
        propsEl.appendChild(kv);
        const kv2 = document.createElement('div');
        kv2.className = 'kv';
        kv2.innerHTML =
          `<span>压降</span><span>${fmtDisplay(ex.dp, 'pressure', state.units)} · R ${ex.r.toPrecision(6)}</span>`;
        propsEl.appendChild(kv2);
      }
    }
    return;
  }

  head.textContent = '属性';
  propsEl.appendChild(head);
}

// ---------------------------------------------------------------------------
// 1D canvas (SVG drag editor)
// ---------------------------------------------------------------------------

function appendDeviceGlyph(p, ax, ay, bx, by) {
  const mx = (ax + bx) / 2;
  const my = (ay + by) / 2;
  const glyph = deviceGlyphCommands(p, mx, my, ax, ay, bx, by);
  if (!glyph) return;
  const NS = 'http://www.w3.org/2000/svg';
  const g = document.createElementNS(NS, 'g');
  g.dataset.pipeId = p.id;
  g.setAttribute('class', glyph.className);
  for (const cmd of glyph.commands) {
    const { type, ...attrs } = cmd;
    const node = document.createElementNS(NS, type);
    for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
    g.appendChild(node);
  }
  const title = document.createElementNS(NS, 'title');
  title.textContent = glyph.title;
  g.appendChild(title);
  g.addEventListener('click', () => {
    state.sel = { kind: 'pipe', ref: p.id };
    state.model.selection = { kind: 'pipe', id: p.id };
    renderCanvas();
    renderProps();
    renderTree();
  });
  svg.appendChild(g);
}

function renderCanvas() {
  const m = state.model;
  svg.innerHTML = '';
  const W = 900;
  const H = 560;

  // pipes
  const vals = state.steady
    ? state.steady.pipes.map((x) => Math.abs(state.colorMode === 'dp' ? x.dp : x.q))
    : [];
  const [lo, hi] = vals.length ? fieldRange(vals) : [0, 1];
  for (const p of m.pipes) {
    const a = m.nodes.find((n) => n.name === p.from);
    const b = m.nodes.find((n) => n.name === p.to);
    if (!a || !b) continue;
    const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
    line.setAttribute('x1', a.x);
    line.setAttribute('y1', a.y);
    line.setAttribute('x2', b.x);
    line.setAttribute('y2', b.y);
    let cls = 'edge';
    let stroke = null;
    const isCoupled = state.coupled && state.coupled.solved && state.coupled.solved.edge === p.name;
    if (isCoupled) cls += ' coupled-edge';
    if (state.colorMode !== 'none' && state.steady) {
      const r = state.steady.pipes.find((x) => x.name === p.name);
      if (r) {
        const mine = Math.abs(state.colorMode === 'dp' ? r.dp : r.q);
        stroke = rainbowColor(mine, lo, hi);
      }
    }
    if (stroke) line.setAttribute('stroke', stroke);
    line.setAttribute('class', cls);
    svg.appendChild(line);

    const res = state.steady && state.steady.pipes.find((x) => x.name === p.name);
    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    label.setAttribute('x', (a.x + b.x) / 2);
    label.setAttribute('y', (a.y + b.y) / 2 - 6);
    label.setAttribute('text-anchor', 'middle');
    label.setAttribute('class', 'edge-label');
    label.textContent = res
      ? `${p.name}  Q ${fmtDisplay(res.q, 'rate', state.units, 4)}  dp ${fmtDisplay(res.dp, 'pressure', state.units, 4)}`
      :       `${p.name}`;
    svg.appendChild(label);

    // Audit fix (P1-3): the gate must cover every device family the glyph
    // generator supports, not just the three M5-era inline keys — otherwise a
    // pipe carrying only a jet pump / PCP / PGL / separator / compressor gets
    // serialized correctly but draws no glyph on the canvas.
    if (p.valve || p.pumpHead || p.pumpRate || ALL_DEVICE_FIELDS.some((k) => p[k])) {
      appendDeviceGlyph(p, a.x, a.y, b.x, b.y);
    }
  }

  // nodes
  for (const n of m.nodes) {
    const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    g.setAttribute('class', 'node' + (state.sel && state.sel.kind === 'node' && state.sel.ref === n.id ? ' selected' : ''));
    g.dataset.nodeId = n.id;
    const c = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    c.setAttribute('cx', n.x);
    c.setAttribute('cy', n.y);
    c.setAttribute('r', 10);
    c.setAttribute('fill', n.pressure !== null && n.pressure !== undefined ? 'var(--accent)' : 'var(--panel)');
    c.setAttribute('stroke', 'var(--text)');
    g.appendChild(c);
    const t = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    t.setAttribute('x', n.x);
    t.setAttribute('y', n.y - 16);
    t.setAttribute('text-anchor', 'middle');
    t.textContent = n.name + (n.pressure !== null && n.pressure !== undefined
      ? ` ${fmtDisplay(n.pressure, 'pressure', state.units, 4)}`
      : '');
    g.appendChild(t);
    svg.appendChild(g);
  }
}

function canvasPoint(evt) {
  const pt = svg.createSVGPoint();
  pt.x = evt.clientX;
  pt.y = evt.clientY;
  return pt.matrixTransform(svg.getScreenCTM().inverse());
}

let dragging = null;
svg.addEventListener('pointerdown', (e) => {
  const pt = canvasPoint(e);
  const hit = hitTestNode(state.model, pt.x, pt.y, 16);
  if (hit) {
    if (state.model.pendingLink) {
      state.model = addPipe(state.model, state.model.pendingLink, hit.id);
      state.sel = { kind: 'pipe', ref: state.model.pipes[state.model.pipes.length - 1]?.id };
      afterModelChange();
      return;
    }
    dragging = { id: hit.id, dx: hit.x - pt.x, dy: hit.y - pt.y };
    state.sel = { kind: 'node', ref: hit.id };
    state.model.selection = { kind: 'node', id: hit.id };
    svg.setPointerCapture(e.pointerId);
    renderCanvas();
    renderProps();
    renderTree();
  }
});
svg.addEventListener('pointermove', (e) => {
  if (!dragging) return;
  const pt = canvasPoint(e);
  state.model = moveNode(state.model, dragging.id, pt.x + dragging.dx, pt.y + dragging.dy, 900, 560);
  renderCanvas();
});
svg.addEventListener('pointerup', () => {
  if (dragging) {
    dragging = null;
    afterModelChange(false);
  }
});
svg.addEventListener('click', (e) => {
  if (e.target.closest('.node')) return;
  const pt = canvasPoint(e);
  if (state.model.pendingLink) {
    state.model.pendingLink = null;
    renderCanvas();
    return;
  }
  const res = addNode(state.model, pt.x, pt.y);
  state.model = res.model;
  state.sel = { kind: 'node', ref: res.id };
  state.model.selection = { kind: 'node', id: res.id };
  afterModelChange();
});

document.addEventListener('keydown', (e) => {
  if (e.key !== 'Delete' && e.key !== 'Backspace') return;
  if (['INPUT', 'TEXTAREA', 'SELECT'].includes(document.activeElement?.tagName)) return;
  if (!state.sel) return;
  if (state.sel.kind === 'node') {
    state.model = removeNode(state.model, state.sel.ref);
  } else if (state.sel.kind === 'pipe') {
    state.model = removePipe(state.model, state.sel.ref);
  } else {
    return;
  }
  state.sel = null;
  afterModelChange();
});

$('addSourceBtn').addEventListener('click', () => {
  if (!state.sel || state.sel.kind !== 'node') {
    setStatus('先在画布或项目树中选中一个节点，再添加源。', 'error');
    return;
  }
  const n = state.model.nodes.find((x) => x.id === state.sel.ref);
  if (!n) return;
  state.model = addSource(state.model, n.name, 1e-6);
  state.sel = { kind: 'source', ref: state.model.sources.at(-1).name };
  afterModelChange();
});
$('deleteBtn').addEventListener('click', () => {
  if (!state.sel) return;
  if (state.sel.kind === 'node') state.model = removeNode(state.model, state.sel.ref);
  else if (state.sel.kind === 'pipe') state.model = removePipe(state.model, state.sel.ref);
  else if (state.sel.kind === 'source') {
    state.model = removeSource(state.model, state.sel.ref);
  } else return;
  state.sel = null;
  afterModelChange();
});
$('colorMode').addEventListener('change', (e) => {
  state.colorMode = e.target.value;
  renderCanvas();
});

// ---------------------------------------------------------------------------
// inline equipment (valve / pump) on the selected pipe
// ---------------------------------------------------------------------------

function selectedPipe() {
  if (!state.sel || state.sel.kind !== 'pipe') return null;
  return state.model.pipes.find((x) => x.id === state.sel.ref) || null;
}

function withSelectedPipe(fn) {
  const p = selectedPipe();
  if (!p) {
    setStatus('先选中一条管段（画布或项目树），再操作设备。', 'error');
    return;
  }
  fn(p);
  afterModelChange();
}

$('addValveBtn').addEventListener('click', () => {
  withSelectedPipe((p) => {
    state.model = setPipeFields(state.model, p.id, { valve: p.valve || 20 });
    setStatus(`管段 ${p.name} 已加节流阀（valve K=${p.valve || 20}）。`, 'ok');
  });
});
$('addPumpBtn').addEventListener('click', () => {
  withSelectedPipe((p) => {
    const patch = {
      pumpHead: p.pumpHead || 25,
      pumpRate: p.pumpRate || 3.0e-4,
    };
    state.model = setPipeFields(state.model, p.id, patch);
    setStatus(`管段 ${p.name} 已加离心泵（pumpHead=${patch.pumpHead} m, pumpRate=${patch.pumpRate} m³/s）。`, 'ok');
  });
});
$('addGasLiftBtn').addEventListener('click', () => {
  withSelectedPipe((p) => {
    const patch = { liftDomePressure: p.liftDomePressure || 1.6e7 };
    state.model = setPipeFields(state.model, p.id, patch);
    setStatus(
      `管段 ${p.name} 已加气举阀（liftDomePressure=${patch.liftDomePressure} Pa，§5.89 稳态元件）。`,
      'ok',
    );
  });
});
$('addRodPumpBtn').addEventListener('click', () => {
  withSelectedPipe((p) => {
    const patch = {
      rodPlungerDiameter: p.rodPlungerDiameter || 0.038,
      rodStrokeLength: p.rodStrokeLength || 2.5,
      rodStrokesPerMinute: p.rodStrokesPerMinute || 6,
      // point element: the kernel refuses a barrel bore on a marched pipe
      // (length) and reads no pipe bore ("the bore that matters is the
      // barrel's"), so both geometry fields go to zero
      length: 0,
      diameter: 0,
    };
    state.model = setPipeFields(state.model, p.id, patch);
    setStatus(
      `管段 ${p.name} 已加有杆泵（bore=${patch.rodPlungerDiameter} m, spm=${patch.rodStrokesPerMinute}，§5.90 稳态点元件——管长已置 0）。`,
      'ok',
    );
  });
});
// E4: one "add" button per new family; all fields get kernel-valid defaults
// so the generated dictionary snippet is runnable without further edits.
for (const fam of E4_DEVICE_FAMILIES) {
  const btnId = `add${fam.id.charAt(0).toUpperCase()}${fam.id.slice(1)}Btn`;
  const btn = $(btnId);
  if (!btn) continue;
  btn.addEventListener('click', () => {
    withSelectedPipe((p) => {
      const patch = {};
      for (const key of fam.fields) patch[key] = p[key] || fam.defaults[key];
      state.model = setPipeFields(state.model, p.id, patch);
      setStatus(
        `管段 ${p.name} 已加${fam.label}（${fam.switch}=${patch[fam.switch]}）。`,
        'ok',
      );
    });
  });
}
$('removeDeviceBtn').addEventListener('click', () => {
  withSelectedPipe((p) => {
    const patch = Object.fromEntries(ALL_DEVICE_FIELDS.map((k) => [k, 0]));
    state.model = setPipeFields(state.model, p.id, patch);
    setStatus(`管段 ${p.name} 的内嵌设备已移除。`, 'ok');
  });
});

// ---------------------------------------------------------------------------
// dict tab
// ---------------------------------------------------------------------------

const caseText = $('caseText');

function syncDictText() {
  caseText.value = networkToCaseText(state.model);
}

$('applyDictBtn').addEventListener('click', () => applyDict(false));
$('runDictBtn').addEventListener('click', () => applyDict(true));
$('validateBtn').addEventListener('click', () => {
  const local = validateCaseText(caseText.value);
  let msg = local.valid
    ? `本地校验通过：${local.nodes.length} 节点 / ${local.pipes.length} 管道 / ${local.sources.length} 源`
    : '本地校验失败:\n' + local.errors.join('\n');
  invoke('validate_case', { caseText: caseText.value })
    .then((data) => {
      setStatus(`${msg}\nRust 校验通过：${data.nodes.length} 节点 / ${data.pipes.length} 管道`, 'ok');
    })
    .catch((e) => setStatus(`${msg}\nRust 校验失败: ${e}`, 'error'));
});
$('resetBtn').addEventListener('click', () => {
  freshModel(DEFAULT_CASE);
  setStatus('已恢复示例算例。', 'ok');
});

function freshModel(text) {
  const r = validateCaseText(text);
  const details = extractCaseDetails(text);
  state.loadedTemplate = null;
  const m = emptyModel();
  if (details.fluid && Number.isFinite(details.fluid.rho)) m.fluid = details.fluid;
  if (details.time && Number.isFinite(details.time.dt)) m.time = details.time;
  if (r.valid) {
    let id = 0;
    for (const n of r.nodes) {
      id += 1;
      m.nodes.push({
        id: `n${id}`,
        name: n.name,
        x: 90 + ((id - 1) % 4) * 200,
        y: 120 + Math.floor((id - 1) / 4) * 150,
        pressure: n.pressure,
      });
    }
    for (const p of r.pipes) {
      id += 1;
      const geo = details.pipes[p.name] || {};
      m.pipes.push({
        id: `p${id}`,
        name: p.name,
        from: p.from,
        to: p.to,
        length: geo.length ?? 1000,
        diameter: geo.diameter ?? 0.05,
        roughness: geo.roughness ?? 0,
        elevation: geo.elevation ?? 0,
        // E4: every equipment field round-trips into the model (not just the
        // six M5-era ones) — the third leg of the three-sync.
        ...pipeDeviceFieldsFrom(geo),
      });
    }
    for (const s of r.sources) {
      m.sources.push({ name: `S${m.sources.length + 1}`, node: s.node, rate: s.rate });
    }
  }
  state.model = m;
  state.sel = null;
  state.steady = null;
  state.transient = null;
  state.coupled = null;
  state.coupledEdge = null;
  afterModelChange();
}

// apply dict text -> model; then optionally run steady
async function applyDict(andRun) {
  const local = validateCaseText(caseText.value);
  if (!local.valid) {
    setStatus('字典校验失败:\n' + local.errors.join('\n'), 'error');
    appendLog('字典校验失败');
    return;
  }
  const text = caseText.value;
  try {
    await invoke('validate_case', { caseText: text });
  } catch (e) {
    setStatus(`Rust 校验失败: ${e}`, 'error');
    appendLog(`Rust 校验失败: ${e}`);
    return;
  }
  freshModel(text);
  // the raw dict stays the run authority until the user edits the model
  // (the drag editor does not model every key, e.g. yieldStress)
  state.rawDictText = text;
  syncDictText();
  setStatus('字典已应用到模型。', 'ok');
  if (andRun) startJob('steady', text);
}

// ---------------------------------------------------------------------------
// job runner (start / poll / stop)
// ---------------------------------------------------------------------------

async function writeCaseFile(text, base) {
  const path = tmpPath(base);
  await invoke('write_text_file', { path, content: text });
  return path;
}

async function startJob(label, caseText_, isPath = false) {
  if (state.job.running) {
    setStatus('已有作业在运行，请先停止。', 'error');
    return;
  }
  const exe =
    label === 'steady' ? engines.steady
    : label === 'pointElements' ? engines.pointElements
    : label === 'transient' ? engines.transient
    : label === 'slug' ? engines.twofluid
    : label === '3d' ? (engines[active3DEngine] || engines.threeD)
    : label === 'coupled' ? (engines[activeCoupledEngine] || engines.coupled)
    : engines.coupled;
  let casePath;
  if (label === 'coupled') {
    casePath = await prepareCoupledCase();
    if (!casePath) return;
  } else if (label === 'slug') {
    casePath = isPath ? caseText_ : await writeCaseFile(caseText_, label);
  } else if (label === '3d') {
    // E1: redirect the case's output{vtk} at a fixed path the field view
    // loads after the run; other output keys survive the rewrite.
    casePath = await writeCaseFile(upsertOutputVtk(caseText_, THREE_D_VTU), label);
  } else {
    casePath = await writeCaseFile(caseText_, label);
  }
  try {
    const id = await invoke('start_case', { label, exePath: exe, casePath });
    state.job = { id, label, running: true, offset: 0, full: '' };
    state.monitorStream = createMonitorState(1500); // F17: fresh iteration clock per job
    setJobChip(true, label);
    appendLog(`[job ${id}] 启动 ${label}: ${exe} ${casePath}`);
    setStatus(`作业 ${label} 已启动。`);
    switchTab(label === 'coupled' ? 'coupled' : label === 'slug' ? 'slug' : label === '3d' ? 'field' : 'canvas');
  } catch (e) {
    appendLog(`启动失败: ${e}`);
    setStatus(`启动失败: ${e}`, 'error');
  }
}

$('tbStop').addEventListener('click', async () => {
  if (!state.job.running) return;
  try {
    await invoke('stop_job', { id: state.job.id });
    appendLog(`[job ${state.job.id}] 已请求停止`);
    setJobChip(false);
    state.job = { id: null, label: null, running: false, offset: 0, full: '' };
    setStatus('作业已停止。', 'error');
  } catch (e) {
    setStatus(`停止失败: ${e}`, 'error');
  }
});

setInterval(async () => {
  if (!state.job.running) return;
  let st;
  try {
    st = await invoke('poll_job', { id: state.job.id });
  } catch (_) {
    return;
  }
  if (st.output && st.output.length > state.job.offset) {
    const chunk = st.output.slice(state.job.offset);
    appendLog(chunk);
    state.job.full += chunk;
    state.job.offset = st.output.length;
    // F17: incremental ITER/MONITOR parsing off the same stream
    feedMonitor(state.monitorStream, chunk);
    if (document.querySelector('#btab-monitor.active')) renderMonitorPanel();
  }
  if (st.finished) {
    const label = state.job.label;
    const full = state.job.full;
    setJobChip(false);
    state.job = { id: null, label: null, running: false, offset: 0, full: '' };
    if (st.exit === 0) {
      appendLog(`[job] ${label} 完成`);
      if (label === 'steady') onSteadyDone(full);
      else if (label === 'pointElements') onPointElementsDone(full);
      else if (label === 'transient') onTransientDone(full);
      else if (label === 'coupled') onCoupledDone(full);
      else if (label === 'slug') onSlugDone(full);
      else if (label === '3d') on3DDone(full);
    } else {
      appendLog(`[job] ${label} 退出码 ${st.exit}`);
      setStatus(`${label} 求解失败（退出码 ${st.exit}），详见求解日志。`, 'error');
    }
  }
}, 300);

// Audit fix (P1-2): the point-element engine prints one verdict line per
// device (PUMPING / NO-SUCTION / OPEN / CLOSED / closure=...); surface them in
// the status area and the run log — no chart plumbing, the physics is closed-form.
function onPointElementsDone(fullOutput) {
  const lines = fullOutput.split('\n').map((l) => l.replace('[chflow] ', '').trim())
    .filter((l) => /(PUMPING|NO-SUCTION|OPEN|CLOSED|SEPARAT|closure=|point elements)/.test(l));
  if (!lines.length) {
    setStatus('点元求解完成（无设备报告行）。');
    return;
  }
  appendLog(`[point-elements] ${lines.join(' | ')}`);
  setStatus(`点元设备求解完成：${lines.length} 条报告行，详见求解日志。`, 'ok');
}

function onSteadyDone(fullOutput) {
  const s = parseSteadyOutput(fullOutput);
  if (!s.nodes.length) {
    setStatus('稳态输出解析失败（日志无 node/pipe 行）。', 'error');
    return;
  }
  state.steady = s;
  setStatus(`稳态完成：iterations=${s.iterations} residual=${s.residual} ${s.converged ? '(converged)' : '(NOT converged)'}`, s.converged ? 'ok' : 'error');
  renderCanvas();
  renderTree();
  renderProps();
  renderCharts();
  renderConvergence();
}

// ---------------------------------------------------------------------------
// §5.185 / §5.187 readout panels (safety.js folds the kernel lines; this only
// paints them). Both stay hidden when the case declares neither block.
// ---------------------------------------------------------------------------

function safetyRowHtml(row) {
  const flag = row.flag ? ` ${row.flag}` : '';
  const note = row.note ? `<span class="note">${escAttr(row.note)}</span>` : '';
  return `<div class="row${flag}"><span class="label">${escAttr(row.label)}</span>` +
    `<span class="value">${escAttr(row.value)}</span>${note}</div>`;
}

function renderWellControlPanel(series) {
  const panel = $('wellControlPanel');
  if (!panel) return;
  let p = null;
  try {
    p = wellControlPanel(series);
  } catch (e) {
    p = null;
  }
  if (!p) {
    panel.classList.add('hidden');
    panel.innerHTML = '';
    return;
  }
  panel.classList.remove('hidden');
  panel.innerHTML = '<h4>井控读数（储层流入 · §5.185）</h4>' +
    p.rows.map(safetyRowHtml).join('') +
    `<div class="summary">${escAttr(p.summary)}</div>`;
}

function renderRiskPanel(frames) {
  const panel = $('riskPanel');
  if (!panel) return;
  let p = null;
  try {
    p = riskPanel(frames);
  } catch (e) {
    p = null;
  }
  if (!p) {
    panel.classList.add('hidden');
    panel.innerHTML = '';
    return;
  }
  panel.classList.remove('hidden');
  panel.innerHTML = '<h4>流动保障读数（只读 · §5.187）</h4>' +
    p.lines.map((l) => safetyRowHtml({
      label: l.label, value: l.value, note: l.detail, flag: l.flag,
    })).join('') +
    `<div class="summary">${escAttr(p.summary)}</div>`;
}

function onTransientDone(fullOutput) {
  const t = parseTransientOutput(fullOutput);
  if (!t.samples.length) {
    setStatus('瞬态输出解析失败（日志无时间帧）。', 'error');
    return;
  }
  state.transient = t;
  renderCharts();
  renderConvergence();
  renderWellControlPanel(t);
  setStatus(`瞬态完成：${t.samples.length} 帧，steps=${t.steps}。`, 'ok');
}

// ---------------------------------------------------------------------------
// coupled view
// ---------------------------------------------------------------------------

async function prepareCoupledCase() {
  const srcPath = $('coupledCasePath').value.trim() || engines.coupledCase;
  let text;
  try {
    text = await invoke('read_text_file', { path: srcPath });
  } catch (e) {
    setStatus(`读取耦合算例失败: ${e}`, 'error');
    appendLog(`读取耦合算例失败: ${e}`);
    return null;
  }
  // redirect the block's vtk output to a fixed path the 3D tab renders
  const patched = /vtk\s+[^;\s]+;/.test(text)
    ? text.replace(/vtk\s+[^;\s]+;/, `vtk ${COUPLED_VTU};`)
    : text.replace(/\n(\s*)solve\s*\{/, `\n$1output { vtk ${COUPLED_VTU}; }\n\n$1solve {`);
  const casePath = tmpPath('coupled');
  await invoke('write_text_file', { path: casePath, content: patched });
  $('vtuPath').value = COUPLED_VTU;
  return casePath;
}

function onCoupledDone(fullOutput) {
  const c = parseCoupledOutput(fullOutput);
  if (!c.solved && !c.pipes.length) {
    setStatus('耦合输出解析失败（日志无 coupled solved 行）。', 'error');
    return;
  }
  state.coupled = c;
  state.coupledEdge = c.solved ? c.solved.edge : null;
  renderCoupled();
  renderTree();
  renderProps();
  renderConvergence();
  switchTab('field');
  loadVtu();
  setStatus('耦合求解完成，交界面已收敛，3D 场已载入。', 'ok');
}

/** E1: after a 3D run, load the fixed-path VTU into the field view. */
function on3DDone() {
  $('vtuPath').value = THREE_D_VTU;
  switchTab('field');
  loadVtu();
  setStatus('3D 求解完成，场已载入 3D 视图。', 'ok');
}

function renderCoupled() {
  const c = state.coupled;
  if (!c) return;
  const u = state.units;
  const nodeRows = c.nodes
    .map((n) => `<tr><td>${n.name}</td><td>${fmtDisplay(n.p, 'pressure', u)}</td></tr>`)
    .join('');
  $('coupledNodes').querySelector('tbody').innerHTML = nodeRows || '<tr><td>—</td></tr>';
  const pipeRows = c.pipes
    .map((p) => `<tr><td>${p.name}</td><td>${fmtDisplay(p.q, 'rate', u)}</td><td>${fmtDisplay(p.dp, 'pressure', u)}</td></tr>`)
    .join('');
  $('coupledPipes').querySelector('tbody').innerHTML = pipeRows || '<tr><td>—</td></tr>';

  const ex = c.solved;
  if (ex) {
    $('coupledBlock').innerHTML =
      `<div class="kv"><span>耦合边</span><span>${ex.edge}</span></div>` +
      `<div class="kv"><span>交换流量 q</span><span>${fmtDisplay(ex.q, 'rate', u)}</span></div>` +
      `<div class="kv"><span>块压降 dp</span><span>${fmtDisplay(ex.dp, 'pressure', u)}</span></div>` +
      `<div class="kv"><span>阻力 R</span><span>${ex.r.toPrecision(8)} Pa·s/m³</span></div>` +
      `<div class="kv"><span>外层迭代</span><span>${ex.iterations}</span></div>` +
      `<div class="kv"><span>残差</span><span>${ex.residual.toExponential(3)}</span></div>`;
  }

  // interface exchange convergence: dp per outer iteration (last seed)
  const seeds = [...new Set(c.interfaces.map((x) => x.seed))];
  const lastSeed = seeds.at(-1);
  const iters = c.interfaces.filter((x) => x.seed === lastSeed);
  if (iters.length >= 2) {
    const pts = toTimeSeriesPoints(
      iters.map((x) => ({ t: x.iteration, v: x.dp })),
      (s) => s.v, 560, 90,
    );
    $('coupledChart').innerHTML =
      `<svg class="chart" viewBox="0 0 600 100" preserveAspectRatio="none">` +
      `<polyline points="${pts.points}" fill="none" stroke="var(--coupled)" stroke-width="2"/></svg>` +
      `<p class="hint">seed ${lastSeed}: dp ${fmtDisplay(iters[0].dp, 'pressure', u)} → ${fmtDisplay(iters.at(-1).dp, 'pressure', u)}（外层迭代收敛）</p>`;
  } else {
    $('coupledChart').innerHTML = '<p class="hint">交界面迭代数据不足（单迭代收敛）。</p>';
  }
}

$('runCoupledBtn').addEventListener('click', async () => {
  const casePath = await prepareCoupledCase();
  if (!casePath) return;
  state.lastCoupledCasePath = casePath;
  startJob('coupled');
});
$('viewFieldBtn').addEventListener('click', () => {
  loadVtu();
  switchTab('field');
});

// ---------------------------------------------------------------------------
// 3D field view (VTU)
// ---------------------------------------------------------------------------

async function loadVtu() {
  const path = $('vtuPath').value.trim();
  if (!path) {
    setStatus('VTU 路径为空。', 'error');
    return;
  }
  try {
    const text = await invoke('read_text_file', { path });
    const vtu = parseVtu(text);
    if (!vtu || vtu.cells.length === 0) {
      setStatus(`VTU 解析失败（不是 VtkWriter 输出？）: ${path}`, 'error');
      return;
    }
    state.vtu = vtu;
    state.streamlineCache = null;
    state.vortCache = null;
    renderField();
    setStatus(`已载入 ${vtu.cells.length} 个单元的场数据。`, 'ok');
  } catch (e) {
    setStatus(`读取 VTU 失败: ${e}`, 'error');
  }
}

function projectPoint(p, c, yaw, pitch) {
  // rotate around the block centroid (yaw about z), then tilt with pitch and
  // take an orthographic x-right / y-up projection
  const dx = p.x - c.x;
  const dy = p.y - c.y;
  const dz = p.z - c.z;
  const cy = Math.cos(yaw);
  const sy = Math.sin(yaw);
  const rx = dx * cy - dy * sy;
  const ry = dx * sy + dy * cy;
  const cp = Math.cos(pitch);
  const sp = Math.sin(pitch);
  return { x: rx, y: ry * cp - dz * sp };
}

function renderField() {
  const cv = $('fieldCanvas');
  const vtu = state.vtu;
  if (!vtu) {
    $('fieldInfo').textContent = '尚无场数据：运行耦合算例后点「查看 3D 场」，或载入 VTU 路径。';
    $('legendLo').textContent = '—';
    $('legendHi').textContent = '—';
    return;
  }
  const { points, cells, cellU, cellP, cellFaces } = vtu;
  // U6: vorticity / Q are Green-Gauss derived quantities over the kernel's
  // own face rings; a VTU without faces cannot support them (note it, fall
  // back to |U| so the view never blanks).
  let isoNote = '';
  let values;
  if (state.fieldMode === 'vort' || state.fieldMode === 'q') {
    if (!state.vortCache) state.vortCache = vorticityAndQ(points, cells, cellFaces, cellU);
    if (!state.vortCache) {
      isoNote = ' · 该 VTU 无 faces 数组，涡量/Q 不可用（退回 |U|）';
      values = cellU.map((u) => Math.hypot(u[0], u[1], u[2]));
    } else {
      values = state.fieldMode === 'vort' ? state.vortCache.mag : state.vortCache.q;
    }
  } else {
    values = state.fieldMode === 'p'
      ? cellP
      : cellU.map((u) => Math.hypot(u[0], u[1], u[2]));
  }
  const [lo, hi] = fieldRange(values);
  const v3 = state.view3d;
  const yaw = (v3.yaw * Math.PI) / 180;
  const pitch = (v3.pitch * Math.PI) / 180;

  // clip: arbitrary plane through cell centroids. `clipAxis` picks a unit
  // axis (x/y/z) or 'custom' to use the free normal `clipNormal`; both are
  // the same plane equation, so the axis form is the n=[1,0,0] special case.
  const clipNormal = v3.clipAxis === 'none'
    ? null
    : v3.clipAxis === 'x' ? [1, 0, 0]
      : v3.clipAxis === 'y' ? [0, 1, 0]
        : v3.clipAxis === 'z' ? [0, 0, 1]
          : v3.clipNormal;
  const keep = clipNormal
    ? clipKeepsByPlane(cells, points, clipNormal, v3.clipFrac)
    : cells.map(() => true);

  const dpr = window.devicePixelRatio || 1;
  cv.width = cv.clientWidth * dpr;
  cv.height = cv.clientHeight * dpr;
  const ctx = cv.getContext('2d');
  ctx.clearRect(0, 0, cv.width, cv.height);

  // projection bounds over visible vertices only
  let cx = 0, cy2 = 0, cz = 0, n = 0;
  for (const c of cells) {
    for (const v of c) { cx += points[v].x; cy2 += points[v].y; cz += points[v].z; n += 1; }
  }
  const centroid = { x: cx / (n || 1), y: cy2 / (n || 1), z: cz / (n || 1) };
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (let k = 0; k < cells.length; k += 1) {
    if (!keep[k]) continue;
    for (const v of cells[k]) {
      const q = projectPoint(points[v], centroid, yaw, pitch);
      xmin = Math.min(xmin, q.x);
      xmax = Math.max(xmax, q.x);
      ymin = Math.min(ymin, q.y);
      ymax = Math.max(ymax, q.y);
    }
  }
  const mx = 60 * dpr;
  const my = 40 * dpr;
  const scale = Math.min(
    (cv.width - 2 * mx) / (xmax - xmin || 1),
    (cv.height - 2 * my) / (ymax - ymin || 1),
  );
  const px = (x) => mx + (x - xmin) * scale;
  const py = (y) => cv.height - my - (y - ymin) * scale;

  values.forEach((val, k) => {
    if (!keep[k]) return;
    const vids = cells[k];
    if (vids.length < 3) return;
    ctx.beginPath();
    const q0 = projectPoint(points[vids[0]], centroid, yaw, pitch);
    ctx.moveTo(px(q0.x), py(q0.y));
    for (let i = 1; i < vids.length; i += 1) {
      const q = projectPoint(points[vids[i]], centroid, yaw, pitch);
      ctx.lineTo(px(q.x), py(q.y));
    }
    ctx.closePath();
    ctx.fillStyle = rainbowColor(val, lo, hi);
    ctx.fill();
    ctx.strokeStyle = 'rgba(0,0,0,0.15)';
    ctx.lineWidth = Math.max(0.5, 0.5 * dpr);
    ctx.stroke();
  });

  // isosurface layer (U6): threshold-crossing face patches of the rendered
  // scalar, drawn as translucent shells over the cell field; patches whose
  // adjacent cells are clipped out are dropped with them.
  let isoCount = 0;
  let isoLevel = NaN;
  if (v3.iso.on && cellFaces && cellFaces.length === cells.length) {
    isoLevel = lo + v3.iso.frac * (hi - lo || 1);
    const iso = isoSurfacePatches(points, cells, cellFaces, values, isoLevel);
    ctx.lineWidth = Math.max(0.8, 1 * dpr);
    iso.patches.forEach((ring, pi) => {
      const [a, b] = iso.src[pi];
      if (!keep[a] || !keep[b]) return;
      ctx.beginPath();
      ring.forEach((p, i) => {
        const q = projectPoint({ x: p[0], y: p[1], z: p[2] }, centroid, yaw, pitch);
        if (i === 0) ctx.moveTo(px(q.x), py(q.y));
        else ctx.lineTo(px(q.x), py(q.y));
      });
      ctx.closePath();
      ctx.fillStyle = 'rgba(255,255,255,0.30)';
      ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,0.75)';
      ctx.stroke();
      isoCount += 1;
    });
  }

  // contour layer (F16): marching squares over the clip plane (or a z-slice
  // through the kept cells when no clip is active), levels evenly spaced
  // inside [lo, hi]; colours reuse the field rainbow so the legend applies.
  let contourLines = 0;
  if (v3.contours.on) {
    const grid = rasterizeSliceField(points, cells, values, keep,
      clipNormal || [0, 0, 1], 28);
    if (grid) {
      const nL = Math.max(2, Math.min(15, Number(v3.contours.levels) || 5));
      const levels = Array.from({ length: nL },
        (_, k) => lo + ((k + 1) / (nL + 1)) * (hi - lo || 1));
      const cs = marchingSquares(grid.field, grid.nx, grid.ny, levels,
        { dx: grid.dx, dy: grid.dy, x0: grid.x0, y0: grid.y0 });
      const wp = (p) => ({
        x: grid.origin[0] + grid.u[0] * p[0] + grid.v[0] * p[1],
        y: grid.origin[1] + grid.u[1] * p[0] + grid.v[1] * p[1],
        z: grid.origin[2] + grid.u[2] * p[0] + grid.v[2] * p[1],
      });
      ctx.lineWidth = Math.max(1, 1.2 * dpr);
      for (const lv of cs) {
        ctx.strokeStyle = contourColor(lv.level, lo, hi);
        for (const [p, q] of lv.segments) {
          const a = projectPoint(wp(p), centroid, yaw, pitch);
          const b = projectPoint(wp(q), centroid, yaw, pitch);
          ctx.beginPath();
          ctx.moveTo(px(a.x), py(a.y));
          ctx.lineTo(px(b.x), py(b.y));
          ctx.stroke();
          contourLines += 1;
        }
      }
    }
  }

  // streamline layer (U7): cell-centred velocity integration with a 64-line
  // budget; colour follows the seed-cell speed on the same rainbow scale as
  // the |U| field, and the geometry is cached per loaded VTU so rotating the
  // view never re-integrates.
  let streamCount = 0;
  const magsAll = cellU.map((u) => Math.hypot(u[0], u[1], u[2]));
  const [magLo, magHi] = fieldRange(magsAll);
  if (v3.streamlines) {
    if (!state.streamlineCache) {
      const centres = cells.map((cids) => cellCentroid(points, cids));
      const stride = Math.max(1, Math.ceil(centres.length / 64));
      const seeds = centres.filter((_, k) => k % stride === 0 && magsAll[k] > 1e-12)
        .map((c) => [c.x, c.y, c.z]);
      state.streamlineCache = traceStreamlines(points, cells, cellU, seeds,
        { maxSteps: 120 });
    }
    ctx.lineWidth = Math.max(0.8, 1.0 * dpr);
    for (const line of state.streamlineCache) {
      ctx.strokeStyle = rainbowColor(magsAll[line.seedCell] ?? magLo, magLo, magHi);
      ctx.beginPath();
      line.forEach((p, i) => {
        const q = projectPoint({ x: p[0], y: p[1], z: p[2] }, centroid, yaw, pitch);
        if (i === 0) ctx.moveTo(px(q.x), py(q.y));
        else ctx.lineTo(px(q.x), py(q.y));
      });
      ctx.stroke();
      streamCount += 1;
    }
  }

  // velocity vectors at cell centroids (magnitude-adaptive thinning to
  // <= 400 arrows: hotspots stay dense, coverage is sampled evenly)
  let drawn = 0;
  if (v3.vectors) {
    const visible = cells.map((_, k) => k).filter((k) => keep[k]);
    const mags = cellU.map((u) => Math.hypot(u[0], u[1], u[2]));
    const drawnIdx = new Set(adaptiveVectorIndices(
      visible.map((k) => mags[k]), 400,
    ).map((local) => visible[local]));
    const maxMag = Math.max(...mags, 1e-12);
    const arrowLen = Math.min(cv.width, cv.height) * 0.09;
    ctx.lineWidth = Math.max(1, 1.2 * dpr);
    for (const k of visible) {
      if (!drawnIdx.has(k)) continue;
      const u = cellU[k];
      const mag = Math.hypot(u[0], u[1], u[2]);
      if (mag <= 0) continue;
      const c = cellCentroid(points, cells[k]);
      const q = projectPoint(c, centroid, yaw, pitch);
      // rotate the vector with the view yaw so arrows match the picture
      const cyaw = Math.cos(yaw);
      const syaw = Math.sin(yaw);
      const ux = u[0] * cyaw - u[1] * syaw;
      const uy = u[0] * syaw + u[1] * cyaw;
      const vp = { x: ux, y: uy * Math.cos(pitch) - u[2] * Math.sin(pitch) };
      // magnitude-based colour on the same scale as the |U| field legend
      const vecColor = rainbowColor(mag, magLo, magHi);
      ctx.strokeStyle = vecColor;
      const s = arrowLen / maxMag;
      const x1 = px(q.x);
      const y1 = py(q.y);
      const x2 = px(q.x + vp.x * s / scale);
      const y2 = py(q.y + vp.y * s / scale);
      ctx.beginPath();
      ctx.moveTo(x1, y1);
      ctx.lineTo(x2, y2);
      ctx.stroke();
      const ang = Math.atan2(y2 - y1, x2 - x1);
      ctx.beginPath();
      ctx.moveTo(x2, y2);
      ctx.lineTo(x2 - 6 * dpr * Math.cos(ang - 0.4), y2 - 6 * dpr * Math.sin(ang - 0.4));
      ctx.lineTo(x2 - 6 * dpr * Math.cos(ang + 0.4), y2 - 6 * dpr * Math.sin(ang + 0.4));
      ctx.closePath();
      ctx.fillStyle = ctx.strokeStyle;
      ctx.fill();
      drawn += 1;
    }
  }

  const FIELD_NAMES = { p: '压力 p', u: '流速 |U|', vort: '涡量 |ω|', q: 'Q 判据' };
  const unit = state.fieldMode === 'p' ? unitLabel('pressure', state.units)
    : state.fieldMode === 'vort' ? '1/s'
      : state.fieldMode === 'q' ? '1/s²' : 'm/s';
  const conv = state.fieldMode === 'p'
    ? (v) => toDisplay(v, 'pressure', state.units)
    : (v) => v;
  $('legendLo').textContent = `${Number(conv(lo).toPrecision(5))} ${unit}`;
  $('legendHi').textContent = `${Number(conv(hi).toPrecision(5))} ${unit}`;
  const clipTxt = v3.clipAxis === 'none'
    ? '无剖切'
    : v3.clipAxis === 'custom'
      ? `剖切 n=[${v3.clipNormal.map((x) => x.toFixed(2)).join(',')}]\u2009@${v3.clipFrac.toFixed(2)}`
      : `剖切 ${v3.clipAxis}@${v3.clipFrac.toFixed(2)}`;
  $('fieldInfo').textContent =
    `块网格 ${cells.length} 单元 / ${points.length} 顶点 · ` +
    `渲染场 ${FIELD_NAMES[state.fieldMode] || state.fieldMode} · ` +
    `yaw ${v3.yaw}° pitch ${v3.pitch}° · ${clipTxt}` +
    (v3.vectors ? ` · 矢量 ${drawn}（按 |U| 着色 ${magLo.toPrecision(3)}–${magHi.toPrecision(3)} m/s）` : '') +
    (v3.streamlines ? ` · 流线 ${streamCount}` : '') +
    (v3.iso.on ? ` · 等值面 ${isoCount} 片 @ ${Number(isoLevel.toPrecision(4))} ${unit}` : '') +
    (v3.contours.on ? ` · 等值线 ${contourLines} 段（${v3.contours.levels} 层）` : '') +
    isoNote;
  recordReports();
}

// ---------------------------------------------------------------------------
// report definitions (U8): face/volume integrals from the kernel's own VTU
// ---------------------------------------------------------------------------

function recordReports() {
  const vtu = state.vtu;
  if (!vtu || !vtu.cellFaces || vtu.cellFaces.length !== vtu.cells.length) return;
  if (state.reports.key !== vtu) {
    state.reports.key = vtu;
    state.reports.history = [];
  }
  if (state.reports.history.length && state.reports.history.at(-1).__vtu === vtu) return;
  const values = evaluateReports(vtu, state.reports.defs);
  const entry = { __vtu: vtu };
  for (const r of values) entry[r.name] = r.value;
  state.reports.history.push(entry);
  renderReports();
}

function renderReports() {
  const defs = state.reports.defs;
  const hist = state.reports.history;
  $('reportDefs').innerHTML = defs.map((d, i) =>
    `<div class="status">#${i + 1} <b>${d.name}</b> · ${d.field} · ${d.op}` +
    (d.normal ? ` · n=[${d.normal.join(',')}] @${(d.frac ?? 0.5).toFixed(2)}` : '') +
    ` <button data-report-remove="${i}" title="移除此定义">✕</button></div>`).join('')
    || '<p class="empty">尚无定义：选场与操作后点「＋ 定义」。</p>';
  const last = hist.at(-1);
  $('reportTable').innerHTML = defs.map((d) => {
    const v = last ? last[d.name] : NaN;
    return `<div class="status">${d.name} = ${Number.isFinite(v) ? v.toPrecision(8) : '—'}` +
      (v === undefined && last ? '（未定义）' : '') + '</div>';
  }).join('') || '<p class="empty">—</p>';

  const W = 520, H = 120;
  const palette = ['var(--accent)', '#e0a33e', '#4bb3c4', '#c96a9b'];
  if (hist.length < 2 || !defs.length) {
    $('reportChart').innerHTML = '<p class="empty">载入 VTU 场后记录历史（每次场更新一点）。</p>';
    return;
  }
  const ys = defs.map((d) => hist.map((h) => h[d.name]).filter(Number.isFinite));
  const all = ys.flat();
  if (!all.length) {
    $('reportChart').innerHTML = '<p class="empty">数据不足。</p>';
    return;
  }
  const lo = Math.min(...all), hi = Math.max(...all);
  const polys = ys.map((vals, i) => {
    const pts = vals.map((v, k) => {
      const x = 36 + (k / Math.max(1, hist.length - 1)) * (W - 56);
      const y = H - 18 - ((v - lo) / (hi - lo || 1)) * (H - 42);
      return `${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(' ');
    return `<polyline points="${pts}" fill="none" stroke="${palette[i % palette.length]}" stroke-width="1.8"/>`;
  }).join('');
  $('reportChart').innerHTML =
    `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">` +
    `<line x1="36" y1="${H - 18}" x2="${W - 10}" y2="${H - 18}" stroke="var(--line)"/>` +
    polys + `</svg><p class="hint">值域 ${lo.toPrecision(4)}..${hi.toPrecision(4)}</p>`;
}

$('reportAddBtn').addEventListener('click', () => {
  const field = $('reportField').value;
  const op = $('reportOp').value;
  const name = `${field}_${op}_${state.reports.defs.length + 1}`;
  const needsNormal = op === 'flux' || op === 'surfaceAverage';
  state.reports.defs.push({
    name, field, op,
    ...(needsNormal ? {
      normal: [Number($('reportNx').value) || 0, Number($('reportNy').value) || 0,
        Number($('reportNz').value) || 0],
      frac: Number($('reportFrac').value),
    } : {}),
  });
  state.reports.history = [];
  recordReports();
});
$('reportDefs').addEventListener('click', (e) => {
  const btn = e.target.closest('[data-report-remove]');
  if (!btn) return;
  state.reports.defs.splice(Number(btn.dataset.reportRemove), 1);
  state.reports.history = [];
  recordReports();
});
$('reportExportBtn').addEventListener('click', async () => {
  if (!state.reports.history.length) {
    setStatus('尚无 report 历史：先载入 VTU 场。', 'error');
    return;
  }
  const path = `/tmp/chflow_gui_reports_${++uidTimer}.csv`;
  try {
    await invoke('write_text_file', {
      path, content: reportsToCsv(state.reports.history, state.reports.defs),
    });
    setStatus(`report 历史已导出：${path}`, 'ok');
  } catch (e) {
    setStatus(`report CSV 导出失败：${e}`, 'error');
  }
});

$('loadVtuBtn').addEventListener('click', loadVtu);
$('fieldMode').addEventListener('change', (e) => {
  state.fieldMode = e.target.value;
  renderField();
});
$('viewYaw').addEventListener('input', (e) => {
  state.view3d.yaw = Number(e.target.value);
  renderField();
});
$('viewPitch').addEventListener('input', (e) => {
  state.view3d.pitch = Number(e.target.value);
  renderField();
});
$('clipAxis').addEventListener('change', (e) => {
  state.view3d.clipAxis = e.target.value;
  $('clipNormalRow').hidden = e.target.value !== 'custom';
  renderField();
});
// custom plane normal: read the three components, normalized in the pure
// helper (a zero vector keeps every cell, so a typo never blanks the view)
function readClipNormal() {
  state.view3d.clipNormal = [
    Number($('clipNx').value) || 0,
    Number($('clipNy').value) || 0,
    Number($('clipNz').value) || 0,
  ];
  renderField();
}
$('clipNx').addEventListener('input', readClipNormal);
$('clipNy').addEventListener('input', readClipNormal);
$('clipNz').addEventListener('input', readClipNormal);
$('clipFrac').addEventListener('input', (e) => {
  state.view3d.clipFrac = Number(e.target.value);
  renderField();
});
$('showVectors').addEventListener('change', (e) => {
  state.view3d.vectors = e.target.checked;
  renderField();
});
$('showStreamlines').addEventListener('change', (e) => {
  state.view3d.streamlines = e.target.checked;
  renderField();
});
$('isoOn').addEventListener('change', (e) => {
  state.view3d.iso.on = e.target.checked;
  renderField();
});
$('isoLevel').addEventListener('input', (e) => {
  state.view3d.iso.frac = Number(e.target.value);
  renderField();
});

// ---------------------------------------------------------------------------
// field export (U9): CFD-Post-style CSV of all cells or of a plane profile
// ---------------------------------------------------------------------------

function exportNormal() {
  const a = $('exportPlaneAxis').value;
  return a === 'x' ? [1, 0, 0]
    : a === 'y' ? [0, 1, 0]
      : a === 'z' ? [0, 0, 1] : state.view3d.clipNormal;
}

async function writeExport(csv, label) {
  const path = $('exportPath').value.trim() || '/tmp/chflow_gui_field_export.csv';
  try {
    await invoke('write_text_file', { path, content: csv });
    setStatus(`${label} -> ${path}`, 'ok');
  } catch (e) {
    setStatus(`导出失败: ${e}`, 'error');
  }
}

$('exportAllCsvBtn').addEventListener('click', () => {
  const vtu = state.vtu;
  if (!vtu) { setStatus('尚无场数据，请先载入 VTU。', 'error'); return; }
  writeExport(fieldsToCsv(vtu.points, vtu.cells, vtu.cellU, vtu.cellP),
    `已导出全部 ${vtu.cells.length} 个单元`);
});
$('exportPlaneCsvBtn').addEventListener('click', () => {
  const vtu = state.vtu;
  if (!vtu) { setStatus('尚无场数据，请先载入 VTU。', 'error'); return; }
  const csv = planeProfileCsv(vtu.points, vtu.cells, vtu.cellFaces, vtu.cellU, vtu.cellP,
    exportNormal(), Number($('exportPlaneFrac').value));
  if (csv === null) { setStatus('该 VTU 无 faces 数组，剖面无法定位。', 'error'); return; }
  writeExport(csv, `已导出剖面 ${csv.trim().split('\n').length - 1} 个单元`);
});

// ---------------------------------------------------------------------------
// slug animation (chflow-1d-twofluid terrain-slug timeline)
// ---------------------------------------------------------------------------

function stopSlugPlayback() {
  const s = state.slug;
  s.playing = false;
  if (s.timer) {
    clearInterval(s.timer);
    s.timer = null;
  }
  $('slugPlayBtn').textContent = '▶ 播放';
}

// U3: overlay palette for the selectable per-cell series (kernel cell keys).
const SLUG_SERIES_COLORS = {
  p: '#e0a33e', T: '#ff8f6b', alpha: '#8fd3a8', water: '#7fb3ff',
  vsl: '#c96a9b', vsg: '#9b8ad1', muEff: '#7db86a', wax: '#c9a06a',
  h: '#5fd3c4', inh: '#4bb3c4', sand: '#d9c26a', bed: '#b08d57',
};
const SLUG_SERIES_FALLBACK = ['#e0a33e', '#ff8f6b', '#c96a9b', '#9b8ad1', '#7db86a'];

function slugSeriesColor(key, i) {
  return SLUG_SERIES_COLORS[key] || SLUG_SERIES_FALLBACK[i % SLUG_SERIES_FALLBACK.length];
}

function renderSlugFrame() {
  const cv = $('slugCanvas');
  const s = state.slug;
  const fr = s.frames[s.idx];
  if (!fr) return;
  const dpr = window.devicePixelRatio || 1;
  cv.width = cv.clientWidth * dpr;
  cv.height = cv.clientHeight * dpr;
  const ctx = cv.getContext('2d');
  const dark = state.theme === 'dark';
  ctx.fillStyle = dark ? '#0c1220' : '#fbfcfe';
  ctx.fillRect(0, 0, cv.width, cv.height);

  const cells = fr.cells;
  const n = cells.length || 1;
  const bands = slugBands(fr);
  // F16 reuse: 1-D holdup contour bands over the same cell axis
  const cBands = state.slug.bands
    ? profileBands(cells.map((c) => (Number.isFinite(c.alpha) ? c.alpha : NaN)), [0.25, 0.5, 0.75])
    : [];
  const bandRow = 12 * dpr;
  const padL = 50 * dpr;
  const padR = 90 * dpr;
  const padT = 30 * dpr;
  const padB = (40 + (bands.length + cBands.length) * 14) * dpr;
  const bw = (cv.width - padL - padR) / n;
  const bh = cv.height - padT - padB;

  cells.forEach((c, i) => {
    const x = padL + i * bw;
    const alpha = Number.isFinite(c.alpha) ? Math.max(0, Math.min(1, c.alpha)) : 0;
    const water = Number.isFinite(c.water) ? Math.max(0, Math.min(1, c.water)) : 0;
    const slug = Number.isFinite(c.s) ? c.s > 0 : false;
    const hA = alpha * bh;
    const hW = Math.min(water, alpha) * bh;
    // gas headspace
    ctx.fillStyle = dark ? '#1a2536' : '#e6ecf3';
    ctx.fillRect(x, padT, bw - 1 * dpr, bh);
    // oil/first-phase holdup
    ctx.fillStyle = '#3fbf6f';
    ctx.fillRect(x, padT + bh - hA, bw - 1 * dpr, hA);
    // water column (bottom, inside the oil bar)
    ctx.fillStyle = '#3a7bd5';
    ctx.fillRect(x, padT + bh - hW, bw - 1 * dpr, hW);
    if (slug) {
      ctx.strokeStyle = '#ff6b6b';
      ctx.lineWidth = 2 * dpr;
      ctx.strokeRect(x, padT, bw - 1 * dpr, bh);
    }
  });

  // U3: selected quantities overlaid as per-frame-normalized polylines.
  // 'p' falls back to the kernel's node pressure profile when the cell rows
  // carry no p (the terrain-slug output shape).
  const chosen = s.series || [];
  const legend = [];
  chosen.forEach((key, i) => {
    const cellHas = cells.some((c) => Number.isFinite(c[key]));
    let pts;
    if (cellHas) {
      pts = cells.map((c, j) => ({ x: padL + (j + 0.5) * bw, v: c[key] }));
    } else if (key === 'p' && fr.nodeP && fr.nodeP.length) {
      pts = fr.nodeP.map((v, j) => ({ x: padL + j * bw, v }));
    } else return;
    if (!pts.some((q) => Number.isFinite(q.v))) return;
    let lo = Infinity;
    let hi = -Infinity;
    for (const q of pts) {
      if (Number.isFinite(q.v)) {
        lo = Math.min(lo, q.v);
        hi = Math.max(hi, q.v);
      }
    }
    const color = slugSeriesColor(key, i);
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.8 * dpr;
    ctx.beginPath();
    let started = false;
    for (const q of pts) {
      if (!Number.isFinite(q.v)) {
        started = false;
        continue;
      }
      const y = padT + bh - ((q.v - lo) / (hi - lo || 1)) * (bh - 6 * dpr) - 3 * dpr;
      if (!started) {
        ctx.moveTo(q.x, y);
        started = true;
      } else ctx.lineTo(q.x, y);
    }
    ctx.stroke();
    legend.push({ key, color, lo, hi });
  });
  ctx.font = `${10 * dpr}px sans-serif`;
  legend.forEach((lg, i) => {
    const y = padT + 12 * dpr + i * 13 * dpr;
    ctx.fillStyle = lg.color;
    ctx.fillText(`${lg.key} ${lg.lo.toPrecision(3)}..${lg.hi.toPrecision(3)}`, cv.width - padR + 8 * dpr, y);
  });

  // U3: flow-assurance band strips (kernel cell flags, contiguous runs)
  bands.forEach((b, r) => {
    const y = padT + bh + (6 + r * 14) * dpr;
    ctx.fillStyle = dark ? '#1a2536' : '#e6ecf3';
    ctx.fillRect(padL, y, (cv.width - padL - padR), bandRow);
    for (const run of b.runs) {
      ctx.fillStyle = b.color;
      ctx.fillRect(padL + run.from * bw, y, (run.to - run.from + 1) * bw, bandRow);
    }
    ctx.fillStyle = dark ? '#9aa7b8' : '#5b6b7d';
    ctx.fillText(b.label, padL - 44 * dpr, y + bandRow - 2 * dpr);
  });

  // F16: holdup contour bands (α >= level), one strip per level
  cBands.forEach((b, r) => {
    const y = padT + bh + (6 + (bands.length + r) * 14) * dpr;
    ctx.fillStyle = dark ? '#1a2536' : '#e6ecf3';
    ctx.fillRect(padL, y, (cv.width - padL - padR), bandRow);
    const col = ['#3fbf6f', '#4bb3c4', '#8fd3a8'][r % 3];
    for (const run of b.runs) {
      ctx.fillStyle = col;
      ctx.fillRect(padL + run.from * bw, y, (run.to - run.from + 1) * bw, bandRow);
    }
    ctx.fillStyle = dark ? '#9aa7b8' : '#5b6b7d';
    ctx.fillText(`α≥${b.level}`, padL - 44 * dpr, y + bandRow - 2 * dpr);
  });

  // axis labels every ~10 cells
  ctx.fillStyle = dark ? '#9aa7b8' : '#5b6b7d';
  ctx.font = `${10 * dpr}px sans-serif`;
  const step = Math.max(1, Math.round(n / 10));
  for (let i = 0; i < n; i += step) {
    const c = cells[i];
    if (!c) continue;
    ctx.fillText(String(c.index ?? i), padL + i * bw, padT + bh + (6 + (bands.length + cBands.length) * 14 + 12) * dpr);
  }
  $('slugReadout').textContent =
    `t=${fr.t} · vm=${fr.vm} · dp=${fr.dp} · 帧 ${s.idx + 1}/${s.frames.length}` +
    (legend.length ? ` · 叠画 ${legend.map((l) => l.key).join('/')}` : '');
  // E2: keep the flow-assurance event playhead in sync with the frame.
  renderFlowAssuranceTimeline();
}

/** U3: rebuild the series checkboxes from the keys the frames actually carry. */
function renderSlugSeries() {
  const box = $('slugSeries');
  const keys = slugSeriesKeys(state.slug.frames);
  const known = ['p', 'T', 'alpha', 'water', 'vsl', 'vsg', 'muEff', 'wax', 'h', 'inh', 'sand', 'bed'];
  const ordered = [...known.filter((k) => keys.includes(k)), ...keys.filter((k) => !known.includes(k))];
  const keep = new Set(state.slug.series.filter((k) => ordered.includes(k)));
  if (!keep.size) {
    for (const k of ['p', 'T', 'alpha']) if (ordered.includes(k)) keep.add(k);
    if (!keep.size && ordered.length) keep.add(ordered[0]);
  }
  state.slug.series = [...keep];
  box.innerHTML = ordered.map((k) =>
    `<label class="chk"><input type="checkbox" data-key="${k}" ` +
    `${keep.has(k) ? 'checked' : ''}/> ${k}</label>`).join('');
}

$('slugSeries').addEventListener('change', (e) => {
  const cb = e.target;
  if (cb.type !== 'checkbox') return;
  const set = new Set(state.slug.series);
  if (cb.checked) set.add(cb.dataset.key);
  else set.delete(cb.dataset.key);
  state.slug.series = [...set];
  renderSlugFrame();
});

function slugStep() {
  const s = state.slug;
  if (!s.frames.length) return;
  s.idx = (s.idx + 1) % s.frames.length;
  $('slugFrame').value = String(s.idx);
  renderSlugFrame();
}

function startSlugPlayback() {
  const s = state.slug;
  if (!s.frames.length) return;
  stopSlugPlayback();
  s.playing = true;
  $('slugPlayBtn').textContent = '‖ 暂停';
  const speed = Number($('slugSpeed').value) || 1;
  s.timer = setInterval(slugStep, 400 / speed);
}

$('slugPlayBtn').addEventListener('click', () => {
  if (state.slug.playing) stopSlugPlayback();
  else startSlugPlayback();
});
$('slugSpeed').addEventListener('change', () => {
  if (state.slug.playing) startSlugPlayback();
});
$('slugFrame').addEventListener('input', (e) => {
  stopSlugPlayback();
  state.slug.idx = Number(e.target.value);
  renderSlugFrame();
});

async function runSlug() {
  let path = $('slugCasePath').value.trim();
  if (!path) {
    // fall back to the kernel-verified twoFluidTerrainSlug template
    try {
      const tpls = await invoke('case_templates');
      const t = (tpls || []).find((x) => x.id === 'twoFluidTerrainSlug');
      if (t) {
        path = t.path;
        $('slugCasePath').value = path;
      }
    } catch (_) { /* keep empty -> error below */ }
  }
  if (!path) {
    setStatus('请填写 twofluid 段塞算例路径（或从模板库载入）。', 'error');
    return;
  }
  stopSlugPlayback();
  await startJob('slug', path, true);
}
$('runSlugBtn').addEventListener('click', runSlug);

// ---------------------------------------------------------------------------
// E2: structured flow-assurance declarations + event timeline
// ---------------------------------------------------------------------------

/** Build the five block editors once; each has an enable checkbox and per-key
 * inputs typed from FA_BLOCKS (num -> number input, word -> select). */
function buildFlowAssuranceGrid() {
  const grid = $('faGrid');
  if (!grid || grid.childElementCount) return;
  for (const name of FA_BLOCK_NAMES) {
    const spec = FA_BLOCKS[name];
    const box = document.createElement('div');
    box.className = 'fa-block';
    box.dataset.block = name;
    const head = document.createElement('div');
    head.className = 'fa-head';
    const chk = document.createElement('input');
    chk.type = 'checkbox';
    chk.id = `fa-${name}-on`;
    const label = document.createElement('label');
    label.htmlFor = chk.id;
    label.textContent = spec.label;
    head.append(chk, label);
    box.append(head);
    const fields = document.createElement('div');
    fields.className = 'fa-fields';
    for (const f of spec.fields) {
      const wrap = document.createElement('label');
      const key = document.createElement('span');
      key.className = 'fa-key';
      key.textContent = f.key;
      key.title = f.hint || '';
      let input;
      if (f.kind === 'word' && f.options) {
        input = document.createElement('select');
        const blank = document.createElement('option');
        blank.value = '';
        blank.textContent = '—';
        input.append(blank);
        for (const opt of f.options) {
          const o = document.createElement('option');
          o.value = opt;
          o.textContent = opt;
          input.append(o);
        }
      } else {
        input = document.createElement('input');
        input.type = 'text';
        input.spellcheck = false;
      }
      input.id = `fa-${name}-${f.key}`;
      input.placeholder = f.required ? '必填' : '默认';
      wrap.append(key, input);
      fields.append(wrap);
    }
    box.append(fields);
    chk.addEventListener('change', () => {
      box.classList.toggle('off', !chk.checked);
    });
    box.classList.add('off');
    grid.append(box);
  }
}

/** Read the form model from the grid inputs. */
function readFlowAssuranceForm() {
  const form = {};
  for (const name of FA_BLOCK_NAMES) {
    const enabled = !!$(`fa-${name}-on`)?.checked;
    const values = {};
    for (const f of FA_BLOCKS[name].fields) {
      const el = $(`fa-${name}-${f.key}`);
      if (el && el.value.trim()) values[f.key] = el.value.trim();
    }
    form[name] = { enabled, values };
  }
  return form;
}

/** Populate the grid inputs from a parsed form model. */
function renderFlowAssuranceForm(form) {
  for (const name of FA_BLOCK_NAMES) {
    const entry = form[name] || { enabled: false, values: {} };
    const chk = $(`fa-${name}-on`);
    if (!chk) continue;
    chk.checked = entry.enabled;
    const box = chk.closest('.fa-block');
    if (box) box.classList.toggle('off', !entry.enabled);
    for (const f of FA_BLOCKS[name].fields) {
      const el = $(`fa-${name}-${f.key}`);
      if (el) el.value = entry.values[f.key] ?? '';
    }
  }
}

/** Apply the current form to the case dictionary textarea. */
function applyFlowAssuranceToCase() {
  const form = readFlowAssuranceForm();
  const text = $('caseText').value;
  let next;
  try {
    next = applyFlowAssurance(text, form);
  } catch (e) {
    setStatus(`流动保障声明校验失败：${e.message}`, 'error');
    return null;
  }
  $('caseText').value = next;
  state.flowAssurance = form;
  return next;
}

/** Apply the declarations, write the resulting dictionary to a temp file and
 * run it through the two-fluid engine (the slug job takes a path). */
async function applyAndRunSlug() {
  const text = applyFlowAssuranceToCase();
  if (text === null) return;
  const path = `/tmp/chflow_gui_flow_assurance_${Date.now()}.dict`;
  try {
    await invoke('write_text_file', { path, content: text });
    $('slugCasePath').value = path;
    stopSlugPlayback();
    await startJob('slug', path, true);
    setStatus('已应用流动保障声明并提交段塞求解。', 'ok');
  } catch (e) {
    setStatus(`应用并运行段塞失败：${e}`, 'error');
  }
}

/** Draw the event timeline onto the small canvas: one lane for slug events,
 * one for purge phases; a playhead tracks the current frame. */
function renderFlowAssuranceTimeline() {
  const cv = $('faTimelineCanvas');
  if (!cv) return;
  const ctx = cv.getContext('2d');
  const w = cv.width = Math.max(320, cv.clientWidth || 640);
  const h = cv.height = 42;
  const light = document.documentElement.dataset.theme !== 'dark';
  ctx.clearRect(0, 0, w, h);
  ctx.fillStyle = light ? '#fafbfc' : '#1b1e23';
  ctx.fillRect(0, 0, w, h);
  const s = state.slug;
  const t0 = s.tMin, t1 = s.tMax;
  const xOf = (t) => 8 + (w - 16) * ((t - t0) / Math.max(1e-9, t1 - t0));
  // lanes: slug events (y 10), purge phases (y 28)
  const colors = light
    ? { GEN: '#2e9e5b', MERGE: '#d48a1e', OUT: '#c0392b', purge: '#3a6ea5' }
    : { GEN: '#5fd08a', MERGE: '#f0b45c', OUT: '#ef6f63', purge: '#6fa3dc' };
  for (const ev of s.events || []) {
    if (!Number.isFinite(ev.t)) continue;
    ctx.fillStyle = colors[ev.kind] || '#888';
    ctx.fillRect(xOf(ev.t) - 2, 6, 4, 10);
  }
  for (const p of s.purge || []) {
    if (!Number.isFinite(p.t)) continue;
    ctx.fillStyle = colors.purge;
    if (p.kind === 'event' && p.phase === 'COMPLETE') {
      // completed cycle: a taller bar so the period is legible on the lane
      ctx.fillRect(xOf(p.t) - 2, 24, 4, 12);
    } else {
      ctx.fillRect(xOf(p.t) - 1, 26, 2, 8);
    }
  }
  // time axis baseline + playhead
  ctx.strokeStyle = light ? '#d0d7de' : '#444b55';
  ctx.beginPath();
  ctx.moveTo(8, 21); ctx.lineTo(w - 8, 21);
  ctx.stroke();
  if (s.frames.length && Number.isFinite(s.frames[s.idx]?.t)) {
    ctx.strokeStyle = light ? '#333' : '#eee';
    ctx.beginPath();
    const px = xOf(s.frames[s.idx].t);
    ctx.moveTo(px, 2); ctx.lineTo(px, h - 2);
    ctx.stroke();
  }
}

buildFlowAssuranceGrid();
$('faLoadBtn')?.addEventListener('click', () => {
  renderFlowAssuranceForm(parseFlowAssurance($('caseText').value));
  setStatus('已从 case 字典载入流动保障声明。');
});
$('faApplyBtn')?.addEventListener('click', () => {
  if (applyFlowAssuranceToCase() !== null) setStatus('已应用流动保障声明到 case 字典。', 'ok');
});
$('faRunBtn')?.addEventListener('click', () => { applyAndRunSlug(); });
// form initialization happens at startup, after caseText gets its default

// ---------------------------------------------------------------------------
// E3: multi-segment schedule timeline editor (valveSchedule / sourceSchedule)
// ---------------------------------------------------------------------------

/** Form model for the schedule panel: valve tables keyed by pipe name (the
 * kernel's own {t, k} rows), plus the top-level source schedule ({t, rate}
 * rows targeting one node). Empty string cells are read as-is and rejected
 * loudly by applySchedules' kernel-aligned guards. */
let schedModel = { valves: {}, source: null };

const escAttr = (v) => String(v).replace(/[&<>"']/g, (c) => (
  { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function schedEventRowsHtml(events, key) {
  return events.map((ev, i) => `
      <div class="sched-evrow" data-i="${i}">
        <input class="sched-t" type="text" spellcheck="false" value="${escAttr(ev.t)}" />
        <input class="sched-v" type="text" spellcheck="false" value="${escAttr(ev[key])}" />
        <button class="mini sched-del-ev" title="删除事件">×</button>
      </div>`).join('');
}

/** Rebuild the panel body from schedModel (timelines are pure SVG strings). */
function renderSchedulePanel() {
  const body = $('schedBody');
  if (!body) return;
  const m = schedModel || { valves: {}, source: null };
  const parts = [];
  for (const [pipe, events] of Object.entries(m.valves)) {
    parts.push(`
      <div class="sched-valve" data-pipe="${escAttr(pipe)}">
        <div class="sched-head">
          <b>阀日程 · ${escAttr(pipe)}</b>
          ${timelineSvg(events, 'k', 400, 80)}
          <button class="mini sched-del-block" title="从模型移除，应用后从字典删除该管段的 valveSchedule 块">移除日程</button>
        </div>
        <div class="sched-evs">
          <div class="sched-evrow sched-evhead"><span>t (s)</span><span>k</span><span></span></div>
          ${schedEventRowsHtml(events, 'k')}
          <button class="mini sched-add-ev">＋ 事件</button>
        </div>
      </div>`);
  }
  const src = m.source || { enabled: false, node: '', events: [] };
  parts.push(`
      <div class="sched-source${src.enabled ? '' : ' off'}">
        <div class="sched-head">
          <label class="chk"><input type="checkbox" class="sched-src-on"${src.enabled ? ' checked' : ''} /> sourceSchedule（顶层注入率事件表）</label>
          <label>node <input class="sched-node" type="text" spellcheck="false" value="${escAttr(src.node)}" /></label>
          ${timelineSvg(src.events, 'rate', 400, 80)}
        </div>
        <div class="sched-evs">
          <div class="sched-evrow sched-evhead"><span>t (s)</span><span>rate</span><span></span></div>
          ${schedEventRowsHtml(src.events, 'rate')}
          <button class="mini sched-add-ev">＋ 事件</button>
        </div>
      </div>`);
  const declared = new Set(Object.keys(m.valves));
  const options = pipeNames($('caseText').value)
    .filter((n) => !declared.has(n))
    .map((n) => `<option value="${escAttr(n)}">${escAttr(n)}</option>`).join('');
  parts.push(`
      <div class="canvas-toolbar sched-addrow">
        <label>给管段加阀日程
          <select id="schedPipeSel">${options || '<option value="">（没有未声明日程的管段）</option>'}</select>
        </label>
        <button id="schedAddValveBtn" class="mini">＋ 添加阀日程</button>
        <span class="hint">留空/非法数在应用时被守卫拒收；相等相邻 t 为 STEP 关断；首尾值保持。</span>
      </div>`);
  body.innerHTML = parts.join('');
}

/** Read the panel rows back into the schedule model. */
function readScheduleForm() {
  const valves = {};
  for (const box of $('schedBody').querySelectorAll('.sched-valve')) {
    valves[box.dataset.pipe] = [...box.querySelectorAll('.sched-evrow[data-i]')]
      .map((row) => ({
        t: row.querySelector('.sched-t').value.trim(),
        k: row.querySelector('.sched-v').value.trim(),
      }));
  }
  const src = $('schedBody').querySelector('.sched-source');
  const source = src ? {
    enabled: src.querySelector('.sched-src-on').checked,
    node: src.querySelector('.sched-node').value.trim(),
    events: [...src.querySelectorAll('.sched-evrow[data-i]')].map((row) => ({
      t: row.querySelector('.sched-t').value.trim(),
      rate: row.querySelector('.sched-v').value.trim(),
    })),
  } : null;
  return { valves, source };
}

/** Apply the panel model to the case dictionary textarea. */
function applySchedulesToCase() {
  let next;
  try {
    next = applySchedules($('caseText').value, readScheduleForm());
  } catch (e) {
    setStatus(`日程事件表校验失败：${e.message}`, 'error');
    return null;
  }
  $('caseText').value = next;
  return next;
}

// structural edits read the current rows first, mutate the model, re-render
$('schedBody').addEventListener('click', (e) => {
  const btn = e.target.closest('button');
  if (!btn) return;
  if (btn.id === 'schedAddValveBtn') {
    const sel = $('schedPipeSel');
    if (!sel || !sel.value) return;
    schedModel = readScheduleForm();
    if (!schedModel.valves[sel.value]) schedModel.valves[sel.value] = [{ t: 0, k: 20 }];
    renderSchedulePanel();
    return;
  }
  const box = btn.closest('[data-pipe], .sched-source');
  schedModel = readScheduleForm();
  if (btn.classList.contains('sched-add-ev')) {
    const ev = { t: '' };
    ev[box.dataset.pipe !== undefined ? 'k' : 'rate'] = '';
    if (box.dataset.pipe !== undefined) schedModel.valves[box.dataset.pipe].push(ev);
    else schedModel.source.events.push(ev);
    renderSchedulePanel();
  } else if (btn.classList.contains('sched-del-ev')) {
    const i = Number(btn.closest('.sched-evrow').dataset.i);
    if (box.dataset.pipe !== undefined) schedModel.valves[box.dataset.pipe].splice(i, 1);
    else schedModel.source.events.splice(i, 1);
    renderSchedulePanel();
  } else if (btn.classList.contains('sched-del-block')) {
    delete schedModel.valves[box.dataset.pipe];
    renderSchedulePanel();
  }
});

$('schedBody').addEventListener('change', (e) => {
  if (e.target.classList.contains('sched-src-on')) {
    e.target.closest('.sched-source').classList.toggle('off', !e.target.checked);
  }
});

$('schedLoadBtn')?.addEventListener('click', () => {
  schedModel = parseSchedules($('caseText').value);
  renderSchedulePanel();
  setStatus('已从 case 字典载入多段日程事件表。');
});
$('schedApplyBtn')?.addEventListener('click', () => {
  if (applySchedulesToCase() !== null) setStatus('已应用多段日程事件表到 case 字典。', 'ok');
});
$('schedEsdBtn')?.addEventListener('click', () => {
  try {
    schedModel = esdPreset(readScheduleForm(),
      Number($('esdTime').value), Number($('esdClosedK').value));
    renderSchedulePanel();
    setStatus('已生成 ESD 预设（全部阀日程在给定时刻 STEP 关闭），点「应用到 case 字典」写入。', 'ok');
  } catch (e) {
    setStatus(`ESD 预设失败：${e.message}`, 'error');
  }
});

function onSlugDone(fullOutput) {
  const frames = parseTwoFluidFrames(fullOutput);
  if (!frames.length) {
    setStatus('段塞输出解析失败（日志无 t/vm/dp 帧行）。', 'error');
    return;
  }
  const s = state.slug;
  s.frames = frames;
  s.idx = 0;
  // E2: parse the kernel's flow-assurance report rows for the event timeline.
  s.events = parseSlugEvents(fullOutput);
  s.purge = parsePurgePhases(fullOutput);
  s.tMin = frames[0].t;
  s.tMax = frames.at(-1).t;
  $('faEventCount').textContent =
    `${s.events.length} 个段塞事件 / ${s.purge.length} 个吹扫相位`;
  renderFlowAssuranceTimeline();
  const slider = $('slugFrame');
  slider.max = String(frames.length - 1);
  slider.value = '0';
  slider.disabled = false;
  $('slugPlayBtn').disabled = false;
  renderSlugSeries();
  renderRiskPanel(frames);
  const last = frames.at(-1);
  const withSlug = last.cells.filter((c) => c.s > 0).length;
  renderSlugFrame();
  startSlugPlayback();
  setStatus(
    `段塞完成：${frames.length} 帧（t ${frames[0].t} → ${last.t} s），` +
    `末帧段塞单元 ${withSlug}/${last.cells.length}。`,
    'ok',
  );
}

// ---------------------------------------------------------------------------
// nodal analysis (U2): IPR-VLP dual curves + working point from the kernel
// ---------------------------------------------------------------------------

async function runNodal() {
  let path = $('nodalCasePath').value.trim();
  if (!path) {
    try {
      const tpls = await invoke('case_templates');
      const t = (tpls || []).find((x) => x.id === 'nodalVogel');
      if (t) {
        path = t.path;
        $('nodalCasePath').value = path;
      }
    } catch (_) { /* keep empty -> error below */ }
  }
  if (!path) {
    setStatus('请填写 nodal 算例路径（或从模板库载入 nodalVogel）。', 'error');
    return;
  }
  let text;
  try {
    text = await invoke('read_text_file', { path });
  } catch (e) {
    setStatus(`读取 nodal 算例失败: ${e}`, 'error');
    return;
  }
  let spec;
  try {
    spec = parseNodalIpr(text);
  } catch (e) {
    setStatus(`解析 ipr{} 失败: ${e}`, 'error');
    return;
  }
  const count = Math.max(5, Math.min(25, Number($('nodalSamples').value) || 13));
  const js = vlpSampleJs(spec, count);
  const texts = [text, ...js.map((j) => iprLinearVariant(text, spec.pRes, j))];
  $('runNodalBtn').disabled = true;
  setStatus(`节点分析运行中：${texts.length} 次内核求解（1 工作点 + ${js.length} 个 VLP 采样）…`, 'busy');
  let batch;
  try {
    batch = await invoke('run_nodal_texts', { texts, exePath: engines.nodal });
  } catch (e) {
    $('runNodalBtn').disabled = false;
    setStatus(`节点分析失败: ${e}`, 'error');
    return;
  }
  $('runNodalBtn').disabled = false;
  const wp = batch[0]?.ok ? batch[0].summary : null;
  if (!wp) {
    setStatus(`节点分析失败：${batch[0]?.error || '无输出'}`, 'error');
    return;
  }
  const vlp = [];
  for (let k = 0; k < js.length; k += 1) {
    const row = batch[k + 1];
    if (row?.ok && row.summary.rate > 0) {
      vlp.push({ j: js[k], q: row.summary.rate, p: row.summary.pwf });
    }
  }
  vlp.sort((a, b) => a.q - b.q);
  state.nodal = { spec, wp, vlp, path };
  renderNodal();
  setStatus(
    `节点分析完成：工作点 q=${wp.rate.toPrecision(4)} sm3/s、p_wf=${wp.pwf.toPrecision(4)} Pa` +
    `（${wp.ipr_form} IPR，${wp.flows ? '可自喷' : '不自喷'}），VLP 采样 ${vlp.length}/${js.length} 点。`,
    'ok',
  );
}

function renderNodal() {
  const d = state.nodal;
  const el = $('nodalChart');
  if (!d) {
    el.innerHTML = '<p class="empty">运行节点分析后绘制 IPR-VLP 双曲线。</p>';
    $('nodalTable').querySelector('tbody').innerHTML = '';
    return;
  }
  const u = state.units;
  const ipr = iprCurvePoints(d.spec, 48);
  const pts = [...ipr, ...d.vlp.map((v) => ({ q: v.q, p: v.p })), { q: d.wp.rate, p: d.wp.pwf }];
  const valid = pts.filter((s) => Number.isFinite(s.q) && Number.isFinite(s.p));
  if (valid.length < 3) {
    el.innerHTML = '<p class="empty">内核未给出可绘采样点（井筒不自喷？检查井口压力）。</p>';
    return;
  }
  const qHi = Math.max(...valid.map((s) => s.q), 1e-12) * 1.05;
  const pLo = Math.min(...valid.map((s) => s.p), d.wp.staticP || Infinity) * 0.9;
  const pHi = d.spec.pRes * 1.02;
  const W = 560;
  const H = 300;
  const X = (q) => 56 + (q / (qHi || 1)) * (W - 76);
  const Y = (p) => H - 30 - ((p - pLo) / (pHi - pLo || 1)) * (H - 54);
  const iprPts = ipr.map((s) => `${X(s.q).toFixed(1)},${Y(s.p).toFixed(1)}`).join(' ');
  const vlpPts = d.vlp.map((s) => `${X(s.q).toFixed(1)},${Y(s.p).toFixed(1)}`).join(' ');
  const wx = X(d.wp.rate);
  const wy = Y(d.wp.pwf);
  const vlpDots = d.vlp.map((s) =>
    `<circle cx="${X(s.q).toFixed(1)}" cy="${Y(s.p).toFixed(1)}" r="2.6" fill="#e0a33e"/>`).join('');
  const staticLine = Number.isFinite(d.wp.staticP)
    ? `<line x1="${X(0)}" y1="${Y(d.wp.staticP).toFixed(1)}" x2="${X(qHi)}" ` +
      `y2="${Y(d.wp.staticP).toFixed(1)}" stroke="var(--muted)" stroke-dasharray="4 3"/>` +
      `<text x="${X(qHi) - 4}" y="${(Y(d.wp.staticP) - 4).toFixed(1)}" text-anchor="end" ` +
      `font-size="10" fill="var(--muted)">静液柱 ${fmtDisplay(d.wp.staticP, 'pressure', u, 4)}</text>`
    : '';
  el.innerHTML =
    `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">` +
    `<line x1="56" y1="${H - 30}" x2="${W - 20}" y2="${H - 30}" stroke="var(--line)"/>` +
    `<line x1="56" y1="12" x2="56" y2="${H - 30}" stroke="var(--line)"/>` +
    staticLine +
    `<polyline points="${iprPts}" fill="none" stroke="#4bb3c4" stroke-width="2"/>` +
    (d.vlp.length > 1
      ? `<polyline points="${vlpPts}" fill="none" stroke="#e0a33e" stroke-width="1.6" stroke-dasharray="1 0"/>`
      : '') +
    vlpDots +
    `<line x1="${(wx - 7).toFixed(1)}" y1="${wy.toFixed(1)}" x2="${(wx + 7).toFixed(1)}" ` +
    `y2="${wy.toFixed(1)}" stroke="var(--bad)" stroke-width="2.4"/>` +
    `<line x1="${wx.toFixed(1)}" y1="${(wy - 7).toFixed(1)}" x2="${wx.toFixed(1)}" ` +
    `y2="${(wy + 7).toFixed(1)}" stroke="var(--bad)" stroke-width="2.4"/>` +
    `<text x="${(wx + 10).toFixed(1)}" y="${(wy - 6).toFixed(1)}" font-size="11" fill="var(--bad)">` +
    `工作点 q=${fmtDisplay(d.wp.rate, 'rate', u, 4)}</text>` +
    `<text x="${W - 20}" y="${H - 12}" text-anchor="end" font-size="10" fill="var(--muted)">` +
    `q (${unitLabel('rate', u)})</text>` +
    `<text x="12" y="20" font-size="10" fill="var(--muted)">p (${unitLabel('pressure', u)})</text>` +
    `<text x="${X(qHi * 0.55)}" y="${Y(ipr[Math.floor(ipr.length / 2)].p) - 6}" font-size="10" fill="#4bb3c4">IPR(${d.spec.form})</text>` +
    (d.vlp.length > 1
      ? `<text x="${X(d.vlp[Math.floor(d.vlp.length / 2)].q) + 8}" y="${Y(d.vlp[Math.floor(d.vlp.length / 2)].p) - 6}" font-size="10" fill="#e0a33e">VLP</text>`
      : '') +
    `</svg>`;
  $('nodalReadout').textContent =
    `工作点：q=${d.wp.rate.toPrecision(6)} sm3/s · p_wf=${d.wp.pwf.toPrecision(6)} Pa · ` +
    `井口=${d.wp.regime || '—'} 持液率=${Number.isFinite(d.wp.holdup) ? d.wp.holdup.toPrecision(4) : '—'} · ` +
    `迭代=${d.wp.iterations} flows=${d.wp.flows} ${d.wp.converged ? '收敛' : '未收敛'}`;
  $('nodalTable').querySelector('tbody').innerHTML = d.vlp.map((s) =>
    `<tr><td>${s.j.toExponential(2)}</td><td>${s.q.toPrecision(4)}</td><td>${s.p.toPrecision(5)}</td></tr>`).join('')
    || '<tr><td>—</td><td>—</td><td>—</td></tr>';
}

$('runNodalBtn').addEventListener('click', runNodal);

// ---------------------------------------------------------------------------
// scenario manager (U1): named boundary/equipment patches, batch run,
// comparison table + overlay curves (steady via run_case_texts, transient
// via the U10 run_transient_texts parallel engine)
// ---------------------------------------------------------------------------

const STEADY_SCEN_PARAMS = [
  { value: 'sourceRate', label: '源汇流量 (sm3/s)' },
  { value: 'pipeDiameter', label: '管径 (m)' },
  { value: 'pipeLength', label: '管长 (m)' },
];
let scenSeq = 0;

function scenarioTargetOptions(engine, param) {
  if (engine === 'steady') {
    if (param === 'sourceRate') return state.model.sources.map((s) => s.name);
    return state.model.pipes.map((p) => p.name);
  }
  const spec = TRANSIENT_WHAT_IF_PARAMS.find((p) => p.value === param);
  if (!spec || spec.section === 'time') return [];
  if (spec.section === 'sources') {
    const names = [];
    const re = /^\s*([A-Za-z_]\w*)\s*\{[\s\S]*?schedule\s*\{/gm;
    let m;
    while ((m = re.exec(caseText.value)) !== null) names.push(m[1]);
    return names;
  }
  return state.model.pipes.map((p) => p.name);
}

function applySteadyScenario(model, overrides) {
  const m = JSON.parse(JSON.stringify(model));
  for (const o of overrides) {
    if (o.param === 'sourceRate') {
      const s = m.sources.find((x) => x.name === o.target);
      if (!s) throw new Error(`稳态工况：源汇 ${o.target} 不存在`);
      s.rate = Number(o.value);
    } else if (o.param === 'pipeDiameter' || o.param === 'pipeLength') {
      const p = m.pipes.find((x) => x.name === o.target);
      if (!p) throw new Error(`稳态工况：管段 ${o.target} 不存在`);
      p[o.param === 'pipeDiameter' ? 'diameter' : 'length'] = Number(o.value);
    } else {
      throw new Error(`稳态工况不支持参数 ${o.param}（边界/阀门时序请切瞬态引擎）`);
    }
  }
  return m;
}

function renderScenarios() {
  const list = $('scenList');
  list.innerHTML = state.scenarios.map((s) =>
    `<div class="scen-item${s.id === state.scenarioSel ? ' sel' : ''}" data-id="${s.id}">` +
    `<span class="nm">${s.name}</span>` +
    `<span class="badge">${s.engine === 'transient' ? '瞬态' : '稳态'}</span>` +
    `<span class="badge">${s.overrides.length} 项</span>` +
    `<button class="mini" data-act="dup" title="复制">⧉</button>` +
    `<button class="mini" data-act="del" title="删除">✕</button></div>`).join('')
    || '<p class="empty">尚无工况，点「＋ 新工况」。</p>';
  renderScenarioEdit();
  renderTree();
}

function renderScenarioEdit() {
  const box = $('scenEdit');
  const sc = state.scenarios.find((s) => s.id === state.scenarioSel);
  if (!sc) {
    box.innerHTML = '<p class="empty">＋ 新工况，或点模型树/列表中的工况。</p>';
    return;
  }
  const params = sc.engine === 'steady' ? STEADY_SCEN_PARAMS : TRANSIENT_WHAT_IF_PARAMS;
  const rows = sc.overrides.map((o, i) => {
    const targets = scenarioTargetOptions(sc.engine, o.param);
    return `<tr><td><select data-i="${i}" data-f="param">` +
      params.map((p) => `<option value="${p.value}"${p.value === o.param ? ' selected' : ''}>${p.label}</option>`).join('') +
      `</select></td><td>${targets.length
        ? `<select data-i="${i}" data-f="target">` +
          targets.map((t) => `<option value="${t}"${t === o.target ? ' selected' : ''}>${t}</option>`).join('') +
          `</select>` : '<i>时间参数</i>'}</td>` +
      `<td><input data-i="${i}" data-f="value" type="number" step="any" value="${o.value}" style="width:110px"/></td>` +
      `<td><button class="mini" data-rm="${i}">✕</button></td></tr>`;
  }).join('');
  box.innerHTML =
    `<label>名称 <input id="scenName" type="text" value="${sc.name}" style="width:160px"/></label>` +
    ` <label>引擎 <select id="scenEngine">` +
    `<option value="steady"${sc.engine === 'steady' ? ' selected' : ''}>稳态</option>` +
    `<option value="transient"${sc.engine === 'transient' ? ' selected' : ''}>瞬态</option></select></label>` +
    `<table class="data-table"><thead><tr><th>参数</th><th>对象</th><th>取值</th><th></th></tr></thead>` +
    `<tbody>${rows}</tbody></table>` +
    `<button id="scenAddOv" class="mini">＋ 加一项</button>`;
}

$('scenAddBtn').addEventListener('click', () => {
  const id = `scen${++scenSeq}`;
  state.scenarios.push({
    id, name: `工况 ${state.scenarios.length + 1}`, engine: 'steady', overrides: [],
  });
  state.scenarioSel = id;
  renderScenarios();
  switchTab('scenarios');
});

$('scenList').addEventListener('click', (e) => {
  const item = e.target.closest('.scen-item');
  if (!item) return;
  const id = item.dataset.id;
  const act = e.target.dataset.act;
  const idx = state.scenarios.findIndex((s) => s.id === id);
  if (act === 'del') {
    state.scenarios.splice(idx, 1);
    if (state.scenarioSel === id) state.scenarioSel = null;
  } else if (act === 'dup') {
    const copy = JSON.parse(JSON.stringify(state.scenarios[idx]));
    copy.id = `scen${++scenSeq}`;
    copy.name = `${copy.name} 副本`;
    state.scenarios.splice(idx + 1, 0, copy);
    state.scenarioSel = copy.id;
  } else {
    state.scenarioSel = id;
  }
  renderScenarios();
});

$('scenEdit').addEventListener('click', (e) => {
  const sc = state.scenarios.find((s) => s.id === state.scenarioSel);
  if (!sc) return;
  if (e.target.id === 'scenAddOv') {
    const params = sc.engine === 'steady' ? STEADY_SCEN_PARAMS : TRANSIENT_WHAT_IF_PARAMS;
    const p0 = params[0];
    const targets = scenarioTargetOptions(sc.engine, p0.value);
    sc.overrides.push({ param: p0.value, target: targets[0] || '', value: 0 });
    renderScenarios();
  } else if (e.target.dataset.rm !== undefined) {
    sc.overrides.splice(Number(e.target.dataset.rm), 1);
    renderScenarios();
  }
});

$('scenEdit').addEventListener('change', (e) => {
  const sc = state.scenarios.find((s) => s.id === state.scenarioSel);
  if (!sc) return;
  const i = e.target.dataset.i;
  const f = e.target.dataset.f;
  if (i !== undefined && f !== undefined) {
    if (f === 'param') {
      const targets = scenarioTargetOptions(sc.engine, e.target.value);
      sc.overrides[Number(i)] = { param: e.target.value, target: targets[0] || '', value: sc.overrides[Number(i)].value };
      renderScenarioEdit();
    } else if (f === 'target') {
      sc.overrides[Number(i)].target = e.target.value;
    } else {
      sc.overrides[Number(i)].value = Number(e.target.value);
    }
    renderTree();
  } else if (e.target.id === 'scenName') {
    sc.name = e.target.value || sc.name;
    renderScenarios();
  } else if (e.target.id === 'scenEngine') {
    sc.engine = e.target.value;
    sc.overrides = [];
    renderScenarios();
  }
});

async function runScenarios() {
  const scs = state.scenarios;
  if (!scs.length) {
    $('scenStatus').textContent = '无工况可运行。';
    return;
  }
  const texts = { steady: [], transient: [] };
  const order = [];
  try {
    for (const s of scs) {
      if (s.engine === 'steady') {
        texts.steady.push(networkToCaseText(applySteadyScenario(state.model, s.overrides)));
      } else {
        if (!caseText.value.trim()) throw new Error(`瞬态工况「${s.name}」需要 case 字典（在「case 字典」页载入）`);
        texts.transient.push(applyScenarioOverrides(caseText.value, s.overrides));
      }
      order.push(s);
    }
  } catch (e) {
    $('scenStatus').textContent = `构建工况字典失败: ${e.message || e}`;
    return;
  }
  $('scenRunAllBtn').disabled = true;
  $('scenStatus').textContent = `批量运行 ${scs.length} 个工况…`;
  const byId = new Map();
  try {
    if (texts.steady.length) {
      const rows = scenarioRowsSteady(await invoke('run_case_texts', {
        texts: texts.steady, exePath: engines.steady,
      }));
      order.filter((s) => s.engine === 'steady').forEach((s, i) => byId.set(s.id, rows[i]));
    }
    if (texts.transient.length) {
      const rows = scenarioRowsTransient(await invoke('run_transient_texts', {
        texts: texts.transient, exePath: engines.transient,
      }));
      order.filter((s) => s.engine === 'transient').forEach((s, i) => byId.set(s.id, rows[i]));
    }
  } catch (e) {
    $('scenRunAllBtn').disabled = false;
    $('scenStatus').textContent = `批量运行失败: ${e}`;
    return;
  }
  $('scenRunAllBtn').disabled = false;
  state.scenarioResults = order.map((s) => ({ name: s.name, engine: s.engine, row: byId.get(s.id) }));
  const done = state.scenarioResults.filter((r) => r.row?.ok).length;
  $('scenStatus').textContent = `批量运行完成：${done}/${scs.length} 成功。`;
  populateScenSeries();
  renderScenarioResults();
}

function populateScenSeries() {
  const sel = $('scenSeriesSel');
  const tRow = (state.scenarioResults || []).find((r) => r.engine === 'transient' && r.row?.series);
  if (tRow) {
    const s = tRow.row.series;
    const opts = s.node_names.map((n) => `p:${n}`)
      .concat(s.pipe_names.map((n) => `Q:${n}`));
    sel.innerHTML = opts.map((o) => `<option value="${o}">${o}</option>`).join('');
    sel.disabled = false;
  } else {
    sel.innerHTML = '<option value="">（稳态：末管 Q 柱状）</option>';
    sel.disabled = true;
  }
}

function renderScenarioResults() {
  const box = $('scenResults');
  const res = state.scenarioResults;
  if (!res || !res.length) {
    box.innerHTML = '<p class="empty">批量运行后显示对比表与叠画曲线。</p>';
    return;
  }
  const u = state.units;
  const f = (v, fam, sig = 5) => (Number.isFinite(v) ? `${Number(toDisplay(v, fam, u)).toPrecision(sig)}` : '—');
  const rows = res.map((r) => {
    const w = r.row;
    if (!w) return `<tr><td>${r.name}</td><td>${r.engine}</td><td colspan="5">未运行</td></tr>`;
    if (!w.ok) return `<tr><td>${r.name}</td><td>${r.engine}</td><td colspan="5" class="bad">${w.error}</td></tr>`;
    if (r.engine === 'transient') {
      return `<tr><td>${r.name}</td><td>瞬态</td><td>${w.steps} 帧</td>` +
        `<td>${f(w.tEnd, 'time')}</td><td>${f(w.pEnd, 'pressure')}</td>` +
        `<td>${f(w.qEnd, 'rate')}</td><td>${f(w.pPeak, 'pressure')}</td></tr>`;
    }
    return `<tr><td>${r.name}</td><td>稳态</td><td>${w.converged ? '✓' : '✗'} ${w.iterations} 次</td>` +
      `<td>${Number.isFinite(w.residual) ? w.residual.toExponential(1) : '—'}</td>` +
      `<td>—</td><td>${f(w.qEnd, 'rate')}</td><td>${f(w.maxDp, 'pressure')}</td></tr>`;
  }).join('');
  // overlay chart
  const W = 600;
  const H = 150;
  let svg = `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">`;
  const legend = [];
  const seriesId = $('scenSeriesSel').value;
  const tRes = res.filter((r) => r.engine === 'transient' && r.row?.ok && r.row.series);
  if (tRes.length && seriesId) {
    tRes.forEach((r, i) => {
      const vals = monitorSeriesValues(r.row.series.samples, seriesId);
      if (!vals || !vals.some(Number.isFinite)) return;
      const fin = vals.filter(Number.isFinite);
      const lo = Math.min(...fin);
      const hi = Math.max(...fin);
      const color = MONITOR_PALETTE[i % MONITOR_PALETTE.length];
      const n = vals.length;
      const pts = vals.map((v, j) => Number.isFinite(v)
        ? `${(10 + (j / Math.max(1, n - 1)) * (W - 20)).toFixed(1)},${(H - 14 - ((v - lo) / (hi - lo || 1)) * (H - 28)).toFixed(1)}`
        : null).filter(Boolean).join(' ');
      svg += `<polyline points="${pts}" fill="none" stroke="${color}" stroke-width="1.5"/>`;
      legend.push({ label: `${r.name} ${lo.toPrecision(3)}..${hi.toPrecision(3)}`, color });
    });
  } else {
    const sRes = res.filter((r) => r.engine === 'steady' && r.row?.ok);
    const use = sRes.length ? sRes : tRes;
    const vals = use.map((r) => r.row.qEnd);
    const hi = Math.max(...vals.filter(Number.isFinite).map(Math.abs), 1e-12);
    const bw = Math.max(10, Math.min(50, (W - 40) / Math.max(1, use.length) - 10));
    use.forEach((r, i) => {
      const v = r.row.qEnd;
      if (!Number.isFinite(v)) return;
      const h = (Math.abs(v) / hi) * (H - 40);
      const x = 30 + i * ((W - 60) / Math.max(1, use.length));
      const color = MONITOR_PALETTE[i % MONITOR_PALETTE.length];
      svg += `<rect x="${x.toFixed(1)}" y="${(H - 20 - h).toFixed(1)}" width="${bw.toFixed(1)}" height="${h.toFixed(1)}" fill="${color}"/>`;
      svg += `<text x="${(x + bw / 2).toFixed(1)}" y="${H - 6}" font-size="9" text-anchor="middle" fill="var(--muted)">${r.name}</text>`;
      legend.push({ label: `${r.name} Q=${v.toPrecision(4)}`, color });
    });
  }
  svg += '</svg>';
  box.innerHTML =
    `<table class="data-table"><thead><tr><th>工况</th><th>引擎</th><th>收敛</th>` +
    `<th>末帧 t</th><th>末 p</th><th>末 Q</th><th>峰值/Δ</th></tr></thead><tbody>${rows}</tbody></table>` +
    svg +
    `<p class="hint">${legend.map((l) => `<span style="color:${l.color}">■</span> ${l.label}`).join(' · ') || '无可叠画序列'}</p>`;
}

$('scenRunAllBtn').addEventListener('click', runScenarios);
$('scenSeriesSel').addEventListener('change', renderScenarioResults);

// ---------------------------------------------------------------------------
// PVT fluid wizard (U4): composition -> chflow-pvt flash -> preview + curves
// ---------------------------------------------------------------------------

function renderPvtComps() {
  $('pvtCompNames').innerHTML = PVT_COMPONENTS.map((n) => `<option value="${n}">`).join('');
  const tb = $('pvtCompTable').querySelector('tbody');
  tb.innerHTML = state.pvt.comps.map((c, i) =>
    `<tr><td><input list="pvtCompNames" data-i="${i}" data-f="name" value="${c.name}" style="width:120px"/></td>` +
    `<td><input type="number" step="any" min="0" data-i="${i}" data-f="z" value="${c.z}" style="width:90px"/></td>` +
    `<td><button class="mini" data-rm="${i}">✕</button></td></tr>`).join('');
}

$('pvtAddComp').addEventListener('click', () => {
  state.pvt.comps.push({ name: '', z: 0 });
  renderPvtComps();
});

$('pvtCompTable').addEventListener('input', (e) => {
  const i = Number(e.target.dataset.i);
  const f = e.target.dataset.f;
  if (Number.isNaN(i) || !f || !state.pvt.comps[i]) return;
  state.pvt.comps[i][f] = f === 'z' ? Number(e.target.value) : e.target.value.trim();
});

$('pvtCompTable').addEventListener('click', (e) => {
  if (e.target.dataset.rm === undefined) return;
  state.pvt.comps.splice(Number(e.target.dataset.rm), 1);
  renderPvtComps();
});

async function runPvt() {
  const temps = $('pvtTemps').value.trim().split(/[\s,]+/).map(Number).filter(Number.isFinite);
  const pressures = $('pvtPressures').value.trim().split(/[\s,]+/).map(Number).filter(Number.isFinite);
  if (!temps.length || !pressures.length) {
    $('pvtReadout').textContent = '温度与压力网格不能为空。';
    return;
  }
  const base = {
    components: state.pvt.comps,
    pressures,
    bubble: $('pvtBubble').checked,
    viscosity: $('pvtVisc').checked,
    surfaceTension: $('pvtSigma').checked,
  };
  $('runPvtBtn').disabled = true;
  $('pvtReadout').textContent = `chflow-pvt 闪蒸中：${temps.length} 个温度 × ${pressures.length} 个压力…`;
  const merged = { nComponents: 0, components: [], molarMass: NaN, flashes: [], bubble: null, tables: [] };
  try {
    for (const T of temps) {
      const text = pvtFlashCaseText({ ...base, temperature: T });
      const out = await invoke('run_pvt_text', { caseText: text, exePath: engines.pvt });
      const one = parsePvtOutput(out);
      merged.nComponents = one.nComponents;
      merged.components = one.components;
      merged.molarMass = one.molarMass;
      merged.flashes.push(...one.flashes);
      merged.tables.push(...one.tables);
      if (one.bubble && !merged.bubble) merged.bubble = one.bubble;
    }
  } catch (e) {
    $('runPvtBtn').disabled = false;
    $('pvtReadout').textContent = `PVT 运行失败: ${e}`;
    return;
  }
  $('runPvtBtn').disabled = false;
  state.pvt.parsed = merged;
  renderPvt();
  $('pvtReadout').textContent =
    `${merged.nComponents} 组分 · M=${merged.molarMass.toPrecision(5)} kg/mol · ` +
    `${merged.flashes.length} 次闪蒸` +
    (merged.bubble ? ` · 泡点 p=${merged.bubble.p.toPrecision(5)} Pa @ T=${merged.bubble.T}` : '');
}

function pvtSvgFamily(title, groups, key, colorSuffix) {
  const W = 300;
  const H = 130;
  const all = groups.flatMap((g) => g.rows.map((r) => r[key])).filter((v) => Number.isFinite(v));
  const ps = groups.flatMap((g) => g.rows.map((r) => r.p)).filter(Number.isFinite);
  if (all.length < 2 || ps.length < 2) {
    return `<h4>${title}：无数据</h4>`;
  }
  const lo = Math.min(...all);
  const hi = Math.max(...all);
  const pLo = Math.min(...ps);
  const pHi = Math.max(...ps);
  const X = (p) => 34 + ((p - pLo) / (pHi - pLo || 1)) * (W - 44);
  const Y = (v) => H - 18 - ((v - lo) / (hi - lo || 1)) * (H - 30);
  let svg = `<h4 style="margin:6px 0 2px;font-size:12px">${title}（${lo.toPrecision(3)}..${hi.toPrecision(3)}）</h4>` +
    `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">` +
    `<line x1="34" y1="${H - 18}" x2="${W - 10}" y2="${H - 18}" stroke="var(--line)"/>` +
    `<line x1="34" y1="8" x2="34" y2="${H - 18}" stroke="var(--line)"/>`;
  groups.forEach((g, i) => {
    const color = MONITOR_PALETTE[i % MONITOR_PALETTE.length];
    const pts = g.rows
      .filter((r) => Number.isFinite(r[key]))
      .map((r) => `${X(r.p).toFixed(1)},${Y(r[key]).toFixed(1)}`)
      .join(' ');
    if (!pts) return;
    svg += `<polyline points="${pts}" fill="none" stroke="${color}" stroke-width="1.5"` +
      `${colorSuffix === 'dash' ? ' stroke-dasharray="4 3"' : ''}/>`;
    svg += `<text x="${W - 12}" y="${12 + i * 11}" font-size="9" text-anchor="end" fill="${color}">T=${g.T}</text>`;
  });
  svg += `<text x="${W - 12}" y="${H - 4}" font-size="9" text-anchor="end" fill="var(--muted)">p (Pa)</text></svg>`;
  return svg;
}

function renderPvt() {
  const d = state.pvt.parsed;
  if (!d) return;
  const tb = $('pvtFlashTable').querySelector('tbody');
  tb.innerHTML = d.flashes.map((f) =>
    `<tr><td>${f.T}</td><td>${f.p.toPrecision(4)}</td><td>${f.beta.toPrecision(3)}${f.single ? ' 单相' : ''}</td>` +
    `<td>${Number.isFinite(f.rhoL) ? f.rhoL.toPrecision(4) : '—'}</td>` +
    `<td>${Number.isFinite(f.rhoV) ? f.rhoV.toPrecision(4) : '—'}</td>` +
    `<td>${Number.isFinite(f.muL) ? f.muL.toExponential(2) : '—'}</td>` +
    `<td>${Number.isFinite(f.sigma) ? f.sigma.toExponential(2) : '—'}</td></tr>`).join('');
  const groups = pvtCurvesByT(d);
  $('pvtChart').innerHTML =
    pvtSvgFamily('ρ vs p · 液相', groups, 'rhoL', '') +
    pvtSvgFamily('ρ vs p · 气相', groups, 'rhoV', 'dash') +
    (groups.some((g) => g.rows.some((r) => Number.isFinite(r.muL)))
      ? pvtSvgFamily('μ vs p · 液相', groups, 'muL', '') +
        pvtSvgFamily('μ vs p · 气相', groups, 'muV', 'dash')
      : '') +
    (groups.some((g) => g.rows.some((r) => Number.isFinite(r.sigma)))
      ? pvtSvgFamily('σ vs p', groups, 'sigma', '')
      : '');
}

$('runPvtBtn').addEventListener('click', runPvt);
renderPvtComps();

// ---------------------------------------------------------------------------
// parameter sweep: steady model serialization or transient raw-dict what-if
// ---------------------------------------------------------------------------

function sweepTargetOptions() {
  const engine = $('sweepEngine').value;
  const param = $('sweepParam').value;
  if (engine === 'transient') {
    const spec = TRANSIENT_WHAT_IF_PARAMS.find((p) => p.value === param);
    if (!spec || spec.section === 'time') return [];
    if (spec.section === 'sources') {
      // Parse the raw dictionary so schedules kept only there stay addressable.
      const names = [];
      const re = /^\s*([A-Za-z_]\w*)\s*\{[\s\S]*?schedule\s*\{/gm;
      let m;
      while ((m = re.exec(caseText.value)) !== null) names.push(m[1]);
      return names;
    }
    return state.model.pipes.map((p) => p.name);
  }
  if (param === 'sourceRate') return state.model.sources.map((s) => s.name);
  return state.model.pipes.map((p) => p.name);
}

function renderSweepTargets(keep = true) {
  const sel = $('sweepTarget');
  const prev = keep ? sel.value : null;
  const opts = sweepTargetOptions();
  sel.innerHTML = opts.map((o) => `<option value="${o}">${o}</option>`).join('');
  sel.disabled = opts.length === 0;
  if (prev && opts.includes(prev)) sel.value = prev;
  const transient = $('sweepEngine').value === 'transient';
  $('sweepStatus').textContent = transient
    ? '瞬态 what-if：对每个取值精确替换当前 case 字典的一个声明（schedule / 阀日程 / time），其余字节不变；并行调真实瞬态引擎，多工况时间序列叠画。'
    : '稳态扫描：基于当前编辑器模型序列化的字典逐值求解（编辑器未建模的字典键不参与）。';
  $('sweepTable').querySelector('thead').innerHTML = transient
    ? '<tr><th>取值</th><th>样本</th><th>终态 p</th><th>终态 Q</th><th>峰值 p</th><th>末时刻</th></tr>'
    : '<tr><th>取值</th><th>收敛</th><th>迭代</th><th>残差</th><th>出口 Q</th><th>最大 dp</th></tr>';
}

$('sweepParam').addEventListener('change', () => renderSweepTargets(false));
$('sweepEngine').addEventListener('change', (e) => {
  state.sweep.engine = e.target.value;
  const sel = $('sweepParam');
  if (state.sweep.engine === 'transient') {
    sel.innerHTML = TRANSIENT_WHAT_IF_PARAMS.map(
      (p) => `<option value="${p.value}">${p.label}</option>`).join('');
  } else {
    sel.innerHTML = '<option value="sourceRate">源注入率</option>' +
      '<option value="pipeDiameter">管段直径</option>' +
      '<option value="pipeLength">管段长度</option>';
  }
  state.sweep.rows = [];
  renderSweepTargets(false);
});

async function runTransientSweep(values, param, target) {
  const baseText = caseText.value;
  let texts;
  try {
    texts = values.map((v) => applyWhatIfOverrides(baseText, param, target, v));
  } catch (e) {
    setStatus(`瞬态 what-if 参数化失败：${e.message || e}`, 'error');
    return;
  }
  $('sweepStatus').textContent = `瞬态 what-if：${values.length} 项已提交真实引擎…`;
  let batch;
  try {
    batch = await invoke('run_transient_texts', { texts, exePath: engines.transient });
  } catch (e) {
    state.sweep.running = false;
    setStatus(`瞬态批量扫描失败：${e}`, 'error');
    return;
  }
  state.sweep.rows = batch.map((r, i) => {
    const value = values[i];
    if (!r || !r.ok) return { value, error: String((r && r.error) || '未知错误'), series: null };
    const s = r.series;
    const samples = s.samples || [];
    const last = samples[samples.length - 1] || {};
    const firstNode = s.node_names[0];
    const firstPipe = s.pipe_names[0];
    return {
      value,
      error: null,
      series: s,
      steps: s.steps,
      tEnd: last.t ?? NaN,
      pEnd: firstNode ? last.pressure?.[firstNode] : NaN,
      qEnd: firstPipe ? last.flow?.[firstPipe] : NaN,
      pPeak: firstNode ? Math.max(...samples.map((x) => x.pressure?.[firstNode] ?? -Infinity)) : NaN,
    };
  });
  state.sweep.running = false;
  const failed = state.sweep.rows.filter((r) => r.error).length;
  $('sweepStatus').textContent =
    `瞬态 what-if 完成：${values.length} 项（失败 ${failed}），时间序列在下方叠画。`;
  renderTransientSweepResults();
}

async function runSweep() {
  if (state.sweep.running) {
    setStatus('扫描进行中，请等待完成。', 'error');
    return;
  }
  if (state.job.running) {
    setStatus('已有求解作业在运行，请先停止。', 'error');
    return;
  }
  const engine = $('sweepEngine').value;
  if (engine === 'steady' && isPartialFidelity(state.loadedTemplate)) {
    setStatus('当前模板为黑油/PVT 偏保真算例：编辑器模型无法完整复刻其字典，扫描已禁用。', 'error');
    return;
  }
  const param = $('sweepParam').value;
  const target = $('sweepTarget').value;
  const spec = engine === 'transient'
    ? TRANSIENT_WHAT_IF_PARAMS.find((p) => p.value === param) : null;
  const values = $('sweepValues').value
    .split(',')
    .map((s) => Number(s.trim()))
    .filter((v) => Number.isFinite(v));
  if (!target && !(spec && spec.section === 'time')) {
    setStatus('扫描作用对象为空：请先建好模型（源或管段），或选一个全局时间参数。', 'error');
    return;
  }
  if (!values.length) {
    setStatus('取值列表为空或非数值。', 'error');
    return;
  }

  state.sweep = { rows: [], running: true, out: state.sweep.out, engine };
  if (engine === 'transient') {
    await runTransientSweep(values, param, target);
    return;
  }
  const maxRef = Math.max(...values, 0);

  // Build every per-value case text up front (clone the model, apply the
  // parameter, serialize), then hand the whole batch to the parallel
  // run_case_texts command: scoped-thread worker pool, one engine process
  // per row, results mapped back by index.
  const texts = values.map((v) => {
    const m = JSON.parse(JSON.stringify(state.model));
    if (param === 'sourceRate') {
      for (const s of m.sources) if (s.name === target) s.rate = v;
    } else {
      for (const p of m.pipes) if (p.name === target) {
        if (param === 'pipeDiameter') p.diameter = v;
        else p.length = v;
      }
    }
    return networkToCaseText(m);
  });
  $('sweepStatus').textContent = `并行扫描：${values.length} 项已提交…`;

  let batch;
  try {
    batch = await invoke('run_case_texts', { texts, exePath: engines.steady });
  } catch (e) {
    state.sweep.running = false;
    setStatus(`并行扫描失败：${e}`, 'error');
    return;
  }

  const rows = batch.map((r, i) => {
    const v = values[i];
    if (!r || !r.ok) {
      return {
        value: v, converged: false, iterations: 0, residual: NaN,
        q: NaN, maxDp: NaN, pipes: [], nodes: [], error: String(r && r.error || '未知错误'),
      };
    }
    const sum = r.summary;
    const pipes = sum.pipes || [];
    const nodes = sum.nodes || [];
    const maxDp = pipes.length ? Math.max(...pipes.map((p) => Math.abs(p.dp || 0))) : 0;
    return {
      value: v,
      converged: !!sum.converged,
      iterations: sum.iterations || 0,
      residual: sum.residual || 0,
      q: pipes.length ? pipes[pipes.length - 1].q : NaN,
      maxDp,
      pipes,
      nodes,
      error: null,
    };
  });
  state.sweep.rows = rows;
  const failed = rows.filter((r) => r.error).length;
  const unconverged = rows.filter((r) => !r.error && !r.converged).length;
  $('sweepStatus').textContent =
    `并行扫描完成：${values.length} 项（失败 ${failed}、未收敛 ${unconverged}）`;
  state.sweep.running = false;
  renderSweepResults(maxRef);
}

function renderSweepResults(valueMax) {
  const u = state.units;
  const rows = state.sweep.rows;
  const tbody = $('sweepTable').querySelector('tbody');
  tbody.innerHTML = rows.map((r) =>
    `<tr>` +
    `<td>${Number(r.value.toPrecision(6))}</td>` +
    `<td>${r.error ? `<span style="color:var(--bad)">错误</span>` : r.converged ? `<span style="color:var(--good)">✓</span>` : '未收敛'}</td>` +
    `<td>${r.iterations}</td>` +
    `<td>${Number.isFinite(r.residual) ? r.residual.toExponential(2) : '—'}</td>` +
    `<td>${Number.isFinite(r.q) ? fmtDisplay(r.q, 'rate', u, 4) : '—'}</td>` +
    `<td>${Number.isFinite(r.maxDp) ? fmtDisplay(r.maxDp, 'pressure', u, 4) : '—'}</td>` +
    `</tr>`,
  ).join('') || '<tr><td>—</td></tr>';

  // output-series select: last-pipe Q / each pipe Q / each node p
  const first = rows.find((r) => !r.error && (r.pipes.length || r.nodes.length));
  const opts = [];
  if (first) {
    for (const p of first.pipes) opts.push(`Q:${p.name}`, `dp:${p.name}`);
    for (const n of first.nodes) opts.push(`p:${n.name}`);
  }
  const sel = $('sweepOut');
  const prev = state.sweep.out;
  sel.innerHTML = opts.map((o) => `<option value="${o}">${o}</option>`).join('');
  if (prev && opts.includes(prev)) { sel.value = prev; state.sweep.out = prev; }
  else if (opts.length) state.sweep.out = opts[0];

  const pick = (r, key) => {
    const [kind, name] = state.sweep.out.split(':');
    if (kind === 'Q') return r.pipes.find((p) => p.name === name)?.q;
    if (kind === 'dp') return r.pipes.find((p) => p.name === name)?.dp;
    return r.nodes.find((n) => n.name === name)?.pressure;
  };
  const series = rows.filter((r) => !r.error).map((r) => ({ t: r.value, v: pick(r) ?? NaN }));
  const valid = series.filter((s) => Number.isFinite(s.v));
  if (valid.length < 2 || !Number.isFinite(valueMax) || valueMax === 0) {
    $('sweepChart').innerHTML = '<p class="empty">数据不足，无法绘制响应曲线。</p>';
    return;
  }
  const vLo = Math.min(...rows.map((r) => r.value).filter(Number.isFinite), 0);
  const yLo = Math.min(...valid.map((s) => s.v));
  const yHi = Math.max(...valid.map((s) => s.v));
  const W = 560;
  const H = 140;
  const pts = valid.map((s) => {
    const x = 40 + ((s.t - vLo) / (valueMax - vLo || 1)) * (W - 60);
    const y = H - 20 - ((s.v - yLo) / (yHi - yLo || 1)) * (H - 50);
    return `${x.toFixed(1)},${y.toFixed(1)}`;
  }).join(' ');
  $('sweepChart').innerHTML =
    `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">` +
    `<line x1="40" y1="${H - 20}" x2="${W - 10}" y2="${H - 20}" stroke="var(--line)"/>` +
    `<polyline points="${pts}" fill="none" stroke="var(--accent)" stroke-width="2"/>` +
    `</svg>` +
    `<p class="hint">${state.sweep.out} 随取值变化 · 值域 ${yLo.toPrecision(4)}..${yHi.toPrecision(4)}</p>`;
}

/** Transient what-if results: one time-series polyline per parameter value. */
function renderTransientSweepResults() {
  const u = state.units;
  const rows = state.sweep.rows;
  const tbody = $('sweepTable').querySelector('tbody');
  tbody.innerHTML = rows.map((r) => `<tr>` +
    `<td>${Number(r.value.toPrecision(6))}</td>` +
    `<td>${r.error ? `<span style="color:var(--bad)">${r.error}</span>` : r.series.samples.length}</td>` +
    `<td>${Number.isFinite(r.pEnd) ? fmtDisplay(r.pEnd, 'pressure', u, 4) : '—'}</td>` +
    `<td>${Number.isFinite(r.qEnd) ? fmtDisplay(r.qEnd, 'rate', u, 4) : '—'}</td>` +
    `<td>${Number.isFinite(r.pPeak) ? fmtDisplay(r.pPeak, 'pressure', u, 4) : '—'}</td>` +
    `<td>${Number.isFinite(r.tEnd) ? r.tEnd.toPrecision(6) : '—'}</td></tr>`).join('')
    || '<tr><td>—</td></tr>';

  const first = rows.find((r) => !r.error && r.series);
  const opts = [];
  if (first) {
    for (const n of first.series.node_names) opts.push(`p:${n}`);
    for (const p of first.series.pipe_names) opts.push(`Q:${p}`);
  }
  const sel = $('sweepOut');
  const prev = state.sweep.out;
  sel.innerHTML = opts.map((o) => `<option value="${o}">${o}</option>`).join('');
  if (prev && opts.includes(prev)) sel.value = prev;
  if (opts.length) state.sweep.out = sel.value;

  const [kind, name] = (state.sweep.out || '').split(':');
  const pick = (s) => {
    const samples = s.samples || [];
    if (kind === 'Q') return samples.map((x) => ({ t: x.t, v: x.flow?.[name] ?? NaN }));
    return samples.map((x) => ({ t: x.t, v: x.pressure?.[name] ?? NaN }));
  };
  const series = rows.filter((r) => !r.error && r.series).map((r) => ({
    label: Number(r.value.toPrecision(5)), points: pick(r.series),
  }));
  const valid = series.flatMap((s) => s.points).filter((p) => Number.isFinite(p.v));
  if (valid.length < 2 || !opts.length) {
    $('sweepChart').innerHTML = '<p class="empty">数据不足，无法绘制工况叠画。</p>';
    return;
  }
  const tHi = Math.max(...valid.map((p) => p.t), 1e-12);
  const yLo = Math.min(...valid.map((p) => p.v));
  const yHi = Math.max(...valid.map((p) => p.v));
  const W = 560, H = 180;
  const palette = ['var(--accent)', '#e0a33e', '#4bb3c4', '#c96a9b', '#7db86a', '#9b8ad1'];
  const polys = series.map((s, i) => {
    const pts = s.points.filter((p) => Number.isFinite(p.v)).map((p) => {
      const x = 40 + (p.t / tHi) * (W - 60);
      const y = H - 24 - ((p.v - yLo) / (yHi - yLo || 1)) * (H - 54);
      return `${x.toFixed(1)},${y.toFixed(1)}`;
    }).join(' ');
    return `<polyline points="${pts}" fill="none" stroke="${palette[i % palette.length]}" stroke-width="1.8"/>` +
      `<text x="${W - 8}" y="${28 + i * 14}" text-anchor="end" font-size="10" ` +
      `fill="${palette[i % palette.length]}">${s.label}</text>`;
  }).join('');
  $('sweepChart').innerHTML =
    `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">` +
    `<line x1="40" y1="${H - 24}" x2="${W - 10}" y2="${H - 24}" stroke="var(--line)"/>` +
    polys + `</svg>` +
    `<p class="hint">${state.sweep.out} 随时间的多工况叠画（每色一个取值） · ` +
    `值域 ${yLo.toPrecision(4)}..${yHi.toPrecision(4)}</p>`;
}

$('sweepRunBtn').addEventListener('click', runSweep);
$('sweepOut').addEventListener('change', (e) => {
  state.sweep.out = e.target.value;
  if (state.sweep.engine === 'transient') {
    renderTransientSweepResults();
    return;
  }
  const maxRef = Math.max(...state.sweep.rows.map((r) => r.value), 0);
  renderSweepResults(maxRef);
});

// ---------------------------------------------------------------------------
// built-in template library (kernel-verified CTest cases)
// ---------------------------------------------------------------------------

async function loadTemplate(entry) {
  let text;
  try {
    text = await invoke('read_text_file', { path: entry.path });
  } catch (e) {
    setStatus(`读取模板失败: ${e}`, 'error');
    return;
  }
  const target = templateTarget(entry.id);
  state.loadedTemplate = entry.id;
  if (target.kind === 'slug') {
    $('slugCasePath').value = entry.path;
    switchTab('slug');
    setStatus(`已载入模板「${entry.name}」：点「▶ 运行段塞」开始动画。`, 'ok');
    return;
  }
  if (target.kind === '3d') {
    // E1 (§5.181): the template names its family member; the ▶ 3D run and
    // the dictionary-page path box follow it instead of hardcoded PISO.
    active3DEngine = target.engine || 'threeD';
    $('threeDEngineSel').value = active3DEngine;
    $('threeDExePath').value = engines[active3DEngine] || engines.threeD;
    $('caseText').value = text;
    const famLabel = (THREE_D_FAMILY.find((e) => e.key === active3DEngine) || {}).label || active3DEngine;
    switchTab('dict');
    setStatus(`已载入模板「${entry.name}」：点「▶ 3D」提交 3D 族引擎（${famLabel}）求解，完成后自动载入 3D 场视图。`, 'ok');
    return;
  }
  if (target.kind === 'coupled') {
    activeCoupledEngine = target.engine || 'coupled';
    $('coupledCasePath').value = entry.path;
    $('coupledExePath').value = engines[activeCoupledEngine] || engines.coupled;
    const famLabel = (COUPLED_FAMILY.find((e) => e.key === activeCoupledEngine) || {}).label || activeCoupledEngine;
    switchTab('coupled');
    setStatus(`已载入模板「${entry.name}」：点「▶ 运行耦合」提交耦合族引擎（${famLabel}）求解。`, 'ok');
    return;
  }
  if (target.kind === 'nodal') {
    $('nodalCasePath').value = entry.path;
    switchTab('nodal');
    setStatus(`已载入模板「${entry.name}」：点「▶ 节点分析」绘制 IPR-VLP 双曲线。`, 'ok');
    return;
  }
  if (target.kind === 'pvt') {
    const spec = parsePvtCaseDict(text);
    if (spec.components.length) state.pvt.comps = spec.components;
    if (Number.isFinite(spec.temperature)) $('pvtTemps').value = String(spec.temperature);
    if (spec.pressures.length) $('pvtPressures').value = spec.pressures.join(' ');
    $('pvtBubble').checked = spec.bubble;
    $('pvtVisc').checked = spec.viscosity;
    $('pvtSigma').checked = spec.surfaceTension;
    renderPvtComps();
    switchTab('pvt');
    setStatus(`已载入模板「${entry.name}」：组分与网格已回填向导，点「▶ 闪蒸预览」。`, 'ok');
    return;
  }
  freshModel(text);
  state.loadedTemplate = entry.id;
  state.rawDictText = text;
  syncDictText();
  switchTab('canvas');
  const fid = isPartialFidelity(entry.id)
    ? '（黑油/PVT 模板：运行与渲染完全可用；拖拽编辑与参数扫描对其禁用，避免静默丢物理。）'
    : '';
  setStatus(`已载入模板「${entry.name}」${fid}`, 'ok');
}

async function renderTemplateMenu() {
  const menu = $('templateMenu');
  menu.innerHTML = '';
  let tpls = [];
  try {
    tpls = await invoke('case_templates') || [];
  } catch (e) {
    const b = document.createElement('button');
    b.textContent = `模板加载失败: ${e}`;
    b.disabled = true;
    menu.appendChild(b);
    return;
  }
  for (const t of tpls) {
    const row = document.createElement('div');
    row.className = 'tpl-row';
    const b = document.createElement('button');
    b.textContent = (isPartialFidelity(t.id) ? '◇ ' : '◆ ') + t.name;
    b.title = t.path;
    b.addEventListener('click', () => loadTemplate(t));
    row.appendChild(b);
    if (isSteadyComparable(t.id)) {
      const c = document.createElement('button');
      c.className = 'mini';
      c.textContent = '对比';
      c.title = `运行并对比：模板「${t.name}」与当前模型并行求解、叠画曲线`;
      c.addEventListener('click', () => runTemplateCompare(t));
      row.appendChild(c);
    }
    menu.appendChild(row);
  }
}

// U13: run a steady template beside the current model and overlay the results
async function runTemplateCompare(t) {
  if (state.job.running) {
    setStatus('已有求解作业在运行，请先停止。', 'error');
    return;
  }
  let tplText;
  try {
    tplText = await invoke('read_text_file', { path: t.path });
  } catch (e) {
    setStatus(`读取模板失败: ${e}`, 'error');
    return;
  }
  const curText = networkToCaseText(state.model);
  $('compareStatus').textContent = '并行求解中…';
  document.querySelector('#bottomTabs [data-btab="curves"]').click();
  let batch;
  try {
    batch = await invoke('run_case_texts', {
      texts: [tplText, curText], exePath: engines.steady,
    });
  } catch (e) {
    $('compareStatus').textContent = `失败: ${e}`;
    return;
  }
  const [a, b] = batch;
  if (!a?.ok || !b?.ok) {
    $('compareStatus').textContent =
      `求解失败：${!a?.ok ? `模板 ${a?.error}` : ''} ${!b?.ok ? `当前模型 ${b?.error}` : ''}`;
    return;
  }
  state.compare = { name: t.name, a: a.summary, b: b.summary };
  renderTemplateCompare();
  $('compareStatus').textContent = `模板「${t.name}」 vs 当前模型：完成`;
}

function renderTemplateCompare() {
  const box = $('tmplCompare');
  const cmp = state.compare;
  if (!cmp) {
    box.innerHTML = '<p class="empty">从模板库点「对比」：模板与当前模型并行求解后在此叠画。</p>';
    return;
  }
  const u = state.units;
  const rows = compareSteadyRuns(cmp.a, cmp.b);
  const fmt = (v, fam) => (Number.isFinite(v) ? `${Number(toDisplay(v, fam, u)).toPrecision(5)}` : '—');
  const table = (list, fam) => list.map((r) =>
    `<tr><td>${r.name}</td><td>${fmt(r.a, fam)}</td><td>${fmt(r.b, fam)}</td>` +
    `<td>${Number.isFinite(r.d) ? r.d.toPrecision(3) : '—'}</td></tr>`).join('');
  // overlay bars: template (accent) vs current (pipe), per-node pressure
  const H = 120;
  const names = rows.nodes.map((r) => r.name);
  const maxP = Math.max(...rows.nodes.flatMap((r) => [r.a, r.b]).filter(Number.isFinite), 1);
  const bw = Math.max(6, Math.min(24, 500 / Math.max(1, names.length * 2)));
  let svg = `<svg class="chart" viewBox="0 0 600 ${H + 20}" preserveAspectRatio="xMidYMid meet">`;
  rows.nodes.forEach((r, i) => {
    const x = 20 + i * (560 / Math.max(1, names.length));
    for (const [v, color, off] of [[r.a, 'var(--accent)', 0], [r.b, 'var(--pipe)', bw + 2]]) {
      if (!Number.isFinite(v)) continue;
      const h = (v / maxP) * (H - 20);
      svg += `<rect x="${(x + off).toFixed(1)}" y="${(H - h + 5).toFixed(1)}" width="${bw.toFixed(1)}" height="${h.toFixed(1)}" fill="${color}"/>`;
    }
    svg += `<text x="${(x + bw).toFixed(1)}" y="${H + 16}" font-size="9" text-anchor="middle" fill="var(--muted)">${r.name}</text>`;
  });
  svg += '</svg>';
  box.innerHTML =
    `<h3>模板并排对比（${cmp.name}） <span class="hint">` +
    `<span style="color:var(--accent)">■</span> 模板 · <span style="color:var(--pipe)">■</span> 当前模型</span></h3>` +
    svg +
    `<div class="cmp-grid">` +
    `<table class="data-table"><thead><tr><th>节点</th><th>p 模板</th><th>p 当前</th><th>Δ</th></tr></thead><tbody>${table(rows.nodes, 'pressure')}</tbody></table>` +
    `<table class="data-table"><thead><tr><th>管段</th><th>Q 模板</th><th>Q 当前</th><th>Δ</th></tr></thead><tbody>${table(rows.pipes, 'rate')}</tbody></table>` +
    `</div>`;
}

// ---------------------------------------------------------------------------
// charts (bottom curves tab)
// ---------------------------------------------------------------------------

function renderBarChart(el, labels, values, unitText, color) {
  if (!values.length) {
    el.innerHTML = '<p class="empty">暂无数据。</p>';
    return;
  }
  const H = 120;
  const bars = toBarHeights(values, H - 30);
  const bw = Math.max(8, Math.min(40, 560 / Math.max(1, labels.length)));
  let svgStr = `<svg class="chart" viewBox="0 0 600 ${H + 20}" preserveAspectRatio="xMidYMid meet">`;
  labels.forEach((lb, i) => {
    const x = 20 + i * (560 / labels.length);
    const h = bars[i];
    svgStr += `<rect x="${x.toFixed(1)}" y="${(H - h + 5).toFixed(1)}" width="${bw.toFixed(1)}" height="${h.toFixed(1)}" fill="${color}"/>`;
    svgStr += `<text x="${(x + bw / 2).toFixed(1)}" y="${H + 16}" font-size="9" text-anchor="middle" fill="var(--muted)">${lb}</text>`;
    svgStr += `<text x="${(x + bw / 2).toFixed(1)}" y="${(H - h).toFixed(1)}" font-size="8" text-anchor="middle" fill="var(--text)">${Number(toDisplay(values[i], unitText.family, state.units).toPrecision(4))}</text>`;
  });
  svgStr += '</svg>';
  el.innerHTML = svgStr;
}

function renderCharts() {
  const u = state.units;
  if (state.steady) {
    renderBarChart(
      $('nodeChart'),
      state.steady.nodes.map((n) => n.name),
      state.steady.nodes.map((n) => n.p),
      { family: 'pressure' },
      'var(--accent)',
    );
    renderBarChart(
      $('pipeChart'),
      state.steady.pipes.map((p) => p.name),
      state.steady.pipes.map((p) => p.q),
      { family: 'rate' },
      'var(--pipe)',
    );
  } else {
    $('nodeChart').innerHTML = '<p class="empty">运行稳态求解后显示。</p>';
    $('pipeChart').innerHTML = '<p class="empty">运行稳态求解后显示。</p>';
  }

  // transient multi-series overlay driven by the monitor tree (U5)
  const t = state.transient;
  if (t && t.samples.length) {
    $('timeSlider').max = String(t.samples.length - 1);
    const idx = Math.min(Number($('timeSlider').value) || 0, t.samples.length - 1);
    const frame = t.samples[idx];
    $('timeReadout').textContent = `t=${frame.t}`;
    renderMonitorChart(t, idx);
  } else {
    $('transientChart').innerHTML = '<p class="empty">运行瞬态求解后显示时序曲线。</p>';
    $('timeReadout').textContent = '—';
  }
}

// U5: monitor tree -> layered checkbox selection, curves generated per the
// quantities the kernel actually declared.
const MONITOR_PALETTE = ['#4bb3c4', '#e0a33e', '#8fd3a8', '#ff8f6b', '#9b8ad1', '#c96a9b', '#7db86a', '#d9c26a'];

/** Series ids the current transient result can produce, pruned + defaulted. */
function syncMonitors(t) {
  const groups = buildMonitorTree(t);
  const avail = new Set();
  for (const g of groups) {
    for (const it of g.items) for (const s of it.series) avail.add(s);
  }
  state.monitors = state.monitors.filter((s) => avail.has(s));
  if (!state.monitors.length) {
    const firstLine = [...avail].find((s) => s !== 'marker' && s !== 'env:p');
    state.monitors = firstLine ? [firstLine, 'marker'] : [...avail];
  }
  return groups;
}

function renderMonitorTree() {
  const box = $('monitorTree');
  const t = state.transient;
  if (!t || !t.samples.length) {
    box.innerHTML = '<p class="empty">运行瞬态求解后生成变量选择树。</p>';
    return;
  }
  const groups = syncMonitors(t);
  box.innerHTML = groups.map((g) => {
    const items = g.items.map((it) => it.series.map((s) => {
      // E10: show the engine-specific hydrate semantics inline so the transient
      // dh/dt=k(1-h) fraction is never read as the two-fluid ln(p) h.
      const note = hydrateSeriesNote(s);
      const tag = note ? ` <span class="model-note" title="${note}">常速率 h ⓘ</span>` : '';
      return `<label class="chk"><input type="checkbox" data-series="${s}" ` +
        `${state.monitors.includes(s) ? 'checked' : ''}/> ${it.series.length > 1 ? s : it.label}${tag}</label>`;
    }).join(''));
    return `<div class="mgroup"><b>${g.label}</b>${items.join('')}</div>`;
  }).join('');
}

function renderMonitorChart(t, idx) {
  syncMonitors(t);
  const W = 600;
  const H = 120;
  const sel = state.monitors.filter((s) => s !== 'marker');
  const X = (i) => (i / Math.max(1, t.samples.length - 1)) * (W - 20) + 10;
  const parts = [`<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">`];
  const legend = [];
  // env:p first, so lines draw on top of the band
  if (sel.includes('env:p')) {
    const band = monitorEnvBand(t.samples, 'p');
    const finite = band.filter((b) => Number.isFinite(b.lo));
    if (finite.length) {
      const lo = Math.min(...finite.map((b) => b.lo));
      const hi = Math.max(...finite.map((b) => b.hi));
      const Y = (v) => H - 10 - ((v - lo) / (hi - lo || 1)) * (H - 20);
      const up = band.map((b, i) => `${X(i).toFixed(1)},${Y(b.hi).toFixed(1)}`).join(' ');
      const dn = band.map((b, i) => `${X(i).toFixed(1)},${Y(b.lo).toFixed(1)}`).reverse().join(' ');
      parts.push(`<polygon points="${up} ${dn}" fill="#4bb3c4" opacity="0.18"/>`);
      legend.push({ label: `压力包络 ${lo.toPrecision(3)}..${hi.toPrecision(3)}`, color: '#4bb3c4' });
    }
  }
  sel.filter((s) => s !== 'env:p').forEach((series, i) => {
    const vals = monitorSeriesValues(t.samples, series);
    if (!vals || !vals.some(Number.isFinite)) return;
    const finite = vals.filter(Number.isFinite);
    const lo = Math.min(...finite);
    const hi = Math.max(...finite);
    const color = MONITOR_PALETTE[i % MONITOR_PALETTE.length];
    const pts = vals.map((v, j) => Number.isFinite(v)
      ? `${X(j).toFixed(1)},${(H - 10 - ((v - lo) / (hi - lo || 1)) * (H - 20)).toFixed(1)}`
      : null).filter(Boolean).join(' ');
    parts.push(`<polyline points="${pts}" fill="none" stroke="${color}" stroke-width="1.5"/>`);
    legend.push({ label: `${series} ${lo.toPrecision(3)}..${hi.toPrecision(3)}`, color });
    if (state.monitors.includes('marker')) {
      const v = vals[idx];
      if (Number.isFinite(v)) {
        parts.push(`<circle cx="${X(idx).toFixed(1)}" cy="${(H - 10 - ((v - lo) / (hi - lo || 1)) * (H - 20)).toFixed(1)}" r="3" fill="var(--coupled)"/>`);
      }
    }
  });
  if (state.monitors.includes('marker')) {
    parts.push(`<line x1="${X(idx).toFixed(1)}" y1="6" x2="${X(idx).toFixed(1)}" y2="${H - 6}" stroke="var(--coupled)" stroke-dasharray="3 3"/>`);
  }
  parts.push('</svg>');
  const u = state.units;
  const kindUnit = (s) => (s.startsWith('p:') || s === 'env:p' ? unitLabel('pressure', u)
    : s.startsWith('Q:') ? unitLabel('rate', u) : '');
  parts.push(`<p class="hint">` + (legend.length
    ? legend.map((l) => `<span style="color:${l.color}">■</span> ${l.label} ${kindUnit(l.label.split(' ')[0])}`).join(' · ')
    : '未勾选变量') + ` · t ${t.samples[0].t}..${t.samples.at(-1).t}</p>`);
  $('transientChart').innerHTML = parts.join('');
}

$('monitorTreeBtn').addEventListener('click', () => {
  const box = $('monitorTree');
  box.classList.toggle('hidden');
  if (!box.classList.contains('hidden')) renderMonitorTree();
});

$('monitorTree').addEventListener('change', (e) => {
  const cb = e.target;
  if (cb.type !== 'checkbox') return;
  const set = new Set(state.monitors);
  if (cb.checked) set.add(cb.dataset.series);
  else set.delete(cb.dataset.series);
  state.monitors = [...set];
  renderCharts();
});

$('timeSlider').addEventListener('input', renderCharts);

// ---------------------------------------------------------------------------
// convergence panel
// ---------------------------------------------------------------------------

function renderConvergence() {
  const parts = [];
  if (state.steady) {
    const s = state.steady;
    parts.push(
      `<p>稳态：iterations=<b>${s.iterations}</b> residual=<b>${s.residual.toExponential(3)}</b> ` +
      `<b style="color:${s.converged ? 'var(--good)' : 'var(--bad)'}">${s.converged ? '收敛' : '未收敛'}</b></p>`,
    );
  }
  if (state.transient && state.transient.samples.length) {
    const t = state.transient;
    parts.push(`<p>瞬态：steps=<b>${t.steps}</b> dt=<b>${t.dt}</b> endTime=<b>${t.endTime}</b> 帧=<b>${t.samples.length}</b></p>`);
  }
  if (state.coupled && state.coupled.solved) {
    const ex = state.coupled.solved;
    parts.push(
      `<p>耦合：edge=<b>${ex.edge}</b> iterations=<b>${ex.iterations}</b> ` +
      `residual=<b>${ex.residual.toExponential(3)}</b> R=<b>${ex.r.toPrecision(7)}</b></p>`,
    );
  }
  convEl.innerHTML = parts.length ? parts.join('') : '<p>等待运行…</p>';
}

// ---------------------------------------------------------------------------
// units / theme / docks
// ---------------------------------------------------------------------------

function applyUnits(system) {
  state.units = system;
  $('unitChip').textContent = system === UNIT_SYSTEMS.SI ? '单位 SI' : '单位 油田';
  $('tbUnits').textContent = system === UNIT_SYSTEMS.SI ? '单位：SI' : '单位：油田';
  renderTree();
  renderProps();
  renderCanvas();
  renderCharts();
  renderConvergence();
  renderCoupled();
}

function applyTheme(theme) {
  state.theme = theme;
  document.documentElement.dataset.theme = theme;
  $('tbTheme').textContent = theme === 'dark' ? '主题：深' : '主题：浅';
  if (state.slug.frames.length) renderSlugFrame();
}

$('tbUnits').addEventListener('click', () => {
  applyUnits(state.units === UNIT_SYSTEMS.SI ? UNIT_SYSTEMS.OILFIELD : UNIT_SYSTEMS.SI);
});
$('tbTheme').addEventListener('click', () => {
  applyTheme(state.theme === 'dark' ? 'light' : 'dark');
});
$('menuUnitsSI').addEventListener('click', () => applyUnits(UNIT_SYSTEMS.SI));
$('menuUnitsOil').addEventListener('click', () => applyUnits(UNIT_SYSTEMS.OILFIELD));
$('menuThemeDark').addEventListener('click', () => applyTheme('dark'));
$('menuThemeLight').addEventListener('click', () => applyTheme('light'));

$('collapseLeft').addEventListener('click', () => {
  state.collapsedLeft = true;
  document.querySelector('.workspace').classList.add('no-left');
});
$('expandLeft').addEventListener('click', () => {
  state.collapsedLeft = false;
  document.querySelector('.workspace').classList.remove('no-left');
});
$('collapseRight').addEventListener('click', () => {
  state.collapsedRight = true;
  document.querySelector('.workspace').classList.add('no-right');
});
$('expandRight').addEventListener('click', () => {
  state.collapsedRight = false;
  document.querySelector('.workspace').classList.remove('no-right');
});
$('menuToggleLeft').addEventListener('click', () => {
  const ws = document.querySelector('.workspace');
  ws.classList.toggle('no-left');
  state.collapsedLeft = ws.classList.contains('no-left');
});
$('menuToggleRight').addEventListener('click', () => {
  const ws = document.querySelector('.workspace');
  ws.classList.toggle('no-right');
  state.collapsedRight = ws.classList.contains('no-right');
});

// ---------------------------------------------------------------------------
// P4 流型图（Beggs-Brill 图，内核 beggsBrillPattern 真值源）
// ---------------------------------------------------------------------------

function renderRegimeMap() {
  const cv = $('regimeCanvas');
  const dpr = window.devicePixelRatio || 1;
  cv.width = cv.clientWidth * dpr;
  cv.height = cv.clientHeight * dpr;
  const ctx = cv.getContext('2d');
  const dark = state.theme === 'dark';
  const D = Number($('regimeD').value) || 0.05;
  const pts = [];
  const op = mapPoint(Number($('regimeVsl').value), Number($('regimeVsg').value), D);
  pts.push({ lambda: op.lambda, nfr: op.nfr, label: '工况点' });
  if ($('regimeSlug').checked) {
    const fr = state.slug.frames[state.slug.idx];
    if (fr && Array.isArray(fr.cells)) {
      for (const c of fr.cells) {
        if (!Number.isFinite(c.vsl) || !Number.isFinite(c.vsg)) continue;
        const m = mapPoint(c.vsl, c.vsg, D);
        pts.push({ lambda: m.lambda, nfr: m.nfr, label: `#${c.index ?? ''}` });
      }
    }
  }
  const { commands } = regimeMapDrawCommands({
    width: cv.width, height: cv.height, points: pts,
  });
  ctx.fillStyle = dark ? '#0c1220' : '#fbfcfe';
  ctx.fillRect(0, 0, cv.width, cv.height);
  ctx.font = `${10 * dpr}px sans-serif`;
  for (const c of commands) {
    if (c.type === 'rect') {
      ctx.globalAlpha = 0.42;
      ctx.fillStyle = c.color;
      ctx.fillRect(c.x, c.y, c.w + 1, c.h + 1);
      ctx.globalAlpha = 1;
    } else if (c.type === 'polyline') {
      ctx.strokeStyle = dark ? '#e8edf4' : '#22303f';
      ctx.lineWidth = 1.4 * dpr;
      ctx.beginPath();
      c.points.forEach((p, i) => (i === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
      ctx.stroke();
    } else if (c.type === 'label') {
      ctx.fillStyle = dark ? '#e8edf4' : '#22303f';
      ctx.fillText(c.text, c.x + 4 * dpr, c.y + 3 * dpr);
    } else if (c.type === 'marker') {
      ctx.strokeStyle = dark ? '#ffffff' : '#111111';
      ctx.lineWidth = 1.6 * dpr;
      ctx.beginPath();
      ctx.arc(c.x, c.y, 4 * dpr, 0, 7);
      ctx.stroke();
      ctx.fillStyle = dark ? '#ffffff' : '#111111';
      ctx.fillText(c.label, c.x + 6 * dpr, c.y - 4 * dpr);
    }
  }
  const slugN = pts.length - 1;
  $('regimeStatus').textContent =
    `工况点：λ=${op.lambda.toPrecision(3)} Nfr=${op.nfr.toPrecision(3)} → ` +
    `${REGIME_NAMES[regimeHit(op.lambda, op.nfr)]}` +
    (slugN ? ` · 叠加段塞逐段点 ${slugN} 个` : ' · 无段塞计算点（先运行段塞）');
  $('regimeLegend').innerHTML = REGIME_NAMES
    .map((n, i) => `<span style="color:${REGIME_COLORS[i]}">■</span> ${n}`).join(' · ')
    + ' · 横轴 λ=v_sl/v_m（无滑移持液率），纵轴 Nfr=v_m²/(gD)（对数）';
}

$('regimeDrawBtn').addEventListener('click', renderRegimeMap);

// ---------------------------------------------------------------------------
// P20 敏感性龙卷风图（复用稳态批量扫描调度）
// ---------------------------------------------------------------------------

function tornadoCurrentValue(param, target) {
  if (param === 'sourceRate') {
    const s = state.model.sources.find((x) => x.name === target);
    return s ? Number(s.rate) : NaN;
  }
  const p = state.model.pipes.find((x) => x.name === target);
  if (!p) return NaN;
  return Number(param === 'pipeDiameter' ? p.diameter : p.length);
}

function steadyMetric(r, metric) {
  if (!r || !r.ok || !r.summary) return NaN;
  const pipes = r.summary.pipes || [];
  if (!pipes.length) return NaN;
  if (metric === 'qEnd') return pipes[pipes.length - 1].q;
  return Math.max(...pipes.map((p) => Math.abs(p.dp || 0)));
}

async function runTornado() {
  if (state.tornado.running) return;
  if ($('sweepEngine').value !== 'steady') {
    $('tornadoStatus').textContent = '龙卷风分析当前基于稳态批量扫描：请把上方引擎切到「稳态管网」。';
    return;
  }
  const rows = state.tornado.rows.length
    ? state.tornado.rows
    : [{ param: $('sweepParam').value, target: $('sweepTarget').value }];
  const lowFrac = (Number($('tornadoLow').value) || 0) / 100;
  const highFrac = (Number($('tornadoHigh').value) || 0) / 100;
  const metric = $('tornadoMetric').value;
  const specs = [];
  for (const row of rows) {
    const cur = tornadoCurrentValue(row.param, row.target);
    if (!Number.isFinite(cur)) {
      $('tornadoStatus').textContent = `参数 ${row.param}@${row.target} 在当前模型中不存在，已跳过。`;
      continue;
    }
    specs.push({ ...row, cur });
  }
  if (!specs.length) {
    $('tornadoStatus').textContent = '没有可分析的参数（当前模型缺少对应对象）。';
    return;
  }
  const texts = [networkToCaseText(state.model)];
  for (const s of specs) {
    for (const frac of [lowFrac, highFrac]) {
      texts.push(networkToCaseText(applySteadyScenario(state.model, [{
        param: s.param, target: s.target, value: s.cur * (1 + frac),
      }])));
    }
  }
  state.tornado.running = true;
  $('tornadoStatus').textContent = `龙卷风分析：基准 + ${specs.length} 参数 × 2 水平已提交并行求解…`;
  let batch;
  try {
    batch = await invoke('run_case_texts', { texts, exePath: engines.steady });
  } catch (e) {
    state.tornado.running = false;
    $('tornadoStatus').textContent = `龙卷风分析失败：${e}`;
    return;
  }
  const cases = specs.map((s, i) => ({
    label: `${s.param}@${s.target}`,
    low: steadyMetric(batch[1 + 2 * i], metric),
    high: steadyMetric(batch[2 + 2 * i], metric),
  }));
  state.tornado.data = tornadoData(steadyMetric(batch[0], metric), cases);
  state.tornado.running = false;
  $('tornadoStatus').textContent =
    `龙卷风分析完成：${cases.length} 个参数（指标 ${metric === 'qEnd' ? '出口 Q' : '最大 dp'}）。`;
  renderTornadoChart();
}

function renderTornadoChart() {
  const el = $('tornadoChart');
  const rows = state.tornado.data;
  if (!rows || !rows.length) {
    el.innerHTML = '<p class="empty">点「▶ 龙卷风分析」：基准值居中，按 |高−低| 降序的双向条形。</p>';
    return;
  }
  const W = 560;
  const H = 40 + rows.length * 34;
  const { commands } = tornadoCommands(rows, { width: W, height: H, marginL: 170 });
  let svg = `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">`;
  for (const c of commands) {
    if (c.type === 'line') {
      svg += `<line x1="${c.x1.toFixed(1)}" y1="${c.y1.toFixed(1)}" x2="${c.x2.toFixed(1)}" `
        + `y2="${c.y2.toFixed(1)}" stroke="var(--text)" stroke-dasharray="4 3"/>`;
    } else if (c.type === 'rect' && c.w > 0) {
      svg += `<rect x="${c.x.toFixed(1)}" y="${c.y.toFixed(1)}" width="${c.w.toFixed(1)}" `
        + `height="${c.h.toFixed(1)}" fill="${c.color}" opacity="0.8"/>`;
    } else if (c.type === 'text') {
      svg += `<text x="${c.x.toFixed(1)}" y="${(c.y + 3).toFixed(1)}" font-size="10" `
        + `text-anchor="${c.anchor}" fill="var(--text)">${c.text}</text>`;
    }
  }
  svg += '</svg>';
  el.innerHTML = svg;
}

$('tornadoAddBtn').addEventListener('click', () => {
  const row = { param: $('sweepParam').value, target: $('sweepTarget').value };
  if (!row.target) {
    $('tornadoStatus').textContent = '当前扫描参数需要作用对象但对象为空：先建模型。';
    return;
  }
  if (state.tornado.rows.some((r) => r.param === row.param && r.target === row.target)) {
    $('tornadoStatus').textContent = `${row.param}@${row.target} 已在敏感性列表中。`;
    return;
  }
  state.tornado.rows.push(row);
  $('tornadoStatus').textContent =
    `敏感性参数：${state.tornado.rows.map((r) => `${r.param}@${r.target}`).join('、')}`;
});
$('tornadoClearBtn').addEventListener('click', () => {
  state.tornado.rows = [];
  state.tornado.data = null;
  renderTornadoChart();
  $('tornadoStatus').textContent = '已清空敏感性参数列表。';
});
$('tornadoRunBtn').addEventListener('click', runTornado);

// ---------------------------------------------------------------------------
// F16 等值线 / 等值带控件
// ---------------------------------------------------------------------------

$('contourOn').addEventListener('change', (e) => {
  state.view3d.contours.on = e.target.checked;
  renderField();
});
$('contourLevels').addEventListener('change', (e) => {
  state.view3d.contours.levels = Number(e.target.value) || 5;
  renderField();
});
$('slugBandsOn').addEventListener('change', (e) => {
  state.slug.bands = e.target.checked;
  if (state.slug.frames.length) renderSlugFrame();
});

// ---------------------------------------------------------------------------
// F18 参数集（Workbench 级 parameter set）
// ---------------------------------------------------------------------------

function renderParamSets() {
  const list = $('psList');
  list.innerHTML = state.paramSets.map((s) =>
    `<div class="scen-item${s.id === state.activeParamSet ? ' sel' : ''}" data-id="${s.id}">` +
    `<span class="nm">${s.name}</span>` +
    `<span class="badge">${s.engine === 'transient' ? '瞬态' : '稳态'}</span>` +
    `<span class="badge">${(s.overrides || []).length} 项</span>` +
    `<button class="mini" data-act="apply" title="激活：精确应用到当前模型/字典">✔</button>` +
    `<button class="mini" data-act="del" title="删除">✕</button></div>`).join('')
    || '<p class="empty">尚无参数集，点「＋ 新参数集」。</p>';
  renderParamSetEdit();
  const conflicts = findConflicts(state.paramSets);
  if (conflicts.length) {
    $('psStatus').textContent = `⚠ 参数集冲突 ${conflicts.length} 处：` +
      conflicts.map((c) => `${c.key}（${c.valueA} vs ${c.valueB}）`).join('、');
  }
}

function renderParamSetEdit() {
  const box = $('psEdit');
  const ps = state.paramSets.find((s) => s.id === state.activeParamSet);
  if (!ps) {
    box.innerHTML = '<p class="empty">＋ 新参数集，或点列表中的参数集。</p>';
    return;
  }
  const params = ps.engine === 'steady' ? STEADY_SCEN_PARAMS : TRANSIENT_WHAT_IF_PARAMS;
  const rows = (ps.overrides || []).map((o, i) => {
    const targets = scenarioTargetOptions(ps.engine, o.param);
    return `<tr><td><select data-i="${i}" data-f="param">` +
      params.map((p) => `<option value="${p.value}"${p.value === o.param ? ' selected' : ''}>${p.label}</option>`).join('') +
      `</select></td><td><select data-i="${i}" data-f="target">` +
      (targets.length
        ? targets.map((t) => `<option${t === o.target ? ' selected' : ''}>${t}</option>`).join('')
        : '<option value="">（全局）</option>') +
      `</select></td><td><input data-i="${i}" data-f="value" type="number" step="any" value="${o.value}"/></td>` +
      `<td><button class="mini" data-act="rm" data-i="${i}">✕</button></td></tr>`;
  }).join('');
  box.innerHTML =
    `<label class="field"><span>名称</span><input id="psName" value="${ps.name}"/></label>` +
    `<label class="field"><span>描述</span><input id="psDesc" value="${ps.description || ''}"/></label>` +
    `<label class="field"><span>引擎</span><select id="psEngine">` +
    `<option value="steady"${ps.engine === 'steady' ? ' selected' : ''}>稳态管网</option>` +
    `<option value="transient"${ps.engine === 'transient' ? ' selected' : ''}>瞬态 what-if</option></select></label>` +
    `<table class="data-table"><thead><tr><th>参数</th><th>对象</th><th>覆盖值</th><th></th></tr></thead>` +
    `<tbody>${rows}</tbody></table>` +
    `<div class="actions"><button id="psAddOverride" class="mini">＋ 覆盖项</button>` +
    `<button id="psMergeBtn" class="mini" title="把上方当前扫描选择并入本参数集">＋ 并入扫描参数</button>` +
    `<button id="psApplyBtn" class="primary">激活（应用到当前模型）</button></div>`;
}

function applyParamSet(ps) {
  try {
    if (ps.engine === 'transient') {
      caseText.value = activateSet(ps, state.model, caseText.value);
      applyDict(false);
    } else {
      state.model = applySetToModel(state.model, ps);
      afterModelChange();
    }
    state.activeParamSet = ps.id;
    setStatus(`参数集「${ps.name}」已激活（${(ps.overrides || []).length} 项覆盖精确应用）。`, 'ok');
    renderParamSets();
  } catch (e) {
    setStatus(`激活参数集失败：${e.message || e}`, 'error');
  }
}

async function runAllParamSets() {
  if (!state.paramSets.length) {
    $('psStatus').textContent = '尚无参数集可求解。';
    return;
  }
  const steady = state.paramSets.filter((s) => s.engine !== 'transient');
  const transient = state.paramSets.filter((s) => s.engine === 'transient');
  const parts = [];
  try {
    if (steady.length) {
      const texts = steady.map((s) => networkToCaseText(applySetToModel(state.model, s)));
      const batch = await invoke('run_case_texts', { texts, exePath: engines.steady });
      parts.push({ engine: 'steady', sets: steady, batch });
    }
    if (transient.length) {
      const texts = transient.map((s) => activateSet(s, state.model, caseText.value));
      const batch = await invoke('run_transient_texts', { texts, exePath: engines.transient });
      parts.push({ engine: 'transient', sets: transient, batch });
    }
  } catch (e) {
    $('psStatus').textContent = `批量求解失败：${e}`;
    return;
  }
  const rows = [];
  for (const part of parts) {
    const norm = part.engine === 'steady' ? scenarioRowsSteady(part.batch) : scenarioRowsTransient(part.batch);
    part.sets.forEach((s, i) => {
      const r = norm[i];
      if (part.engine === 'steady') {
        rows.push({
          name: s.name, engine: '稳态',
          ok: r.ok ? (r.converged ? '✓ 收敛' : '未收敛') : `错误：${r.error}`,
          iters: r.ok ? r.iterations : '—',
          resid: r.ok ? Number(r.residual).toExponential(2) : '—',
          a: r.ok ? fmtDisplay(r.maxDp, 'pressure', state.units, 4) : '—',
          b: r.ok ? fmtDisplay(r.qEnd, 'rate', state.units, 4) : '—',
        });
      } else {
        rows.push({
          name: s.name, engine: '瞬态',
          ok: r.ok ? `${r.steps} 步` : `错误：${r.error}`,
          iters: r.ok ? Number(r.tEnd).toPrecision(4) : '—',
          resid: '—',
          a: r.ok ? fmtDisplay(r.pEnd, 'pressure', state.units, 4) : '—',
          b: r.ok ? fmtDisplay(r.qEnd, 'rate', state.units, 4) : '—',
        });
      }
    });
  }
  $('psResults').innerHTML =
    `<table class="data-table"><thead><tr><th>参数集</th><th>引擎</th><th>状态</th>` +
    `<th>迭代/时长</th><th>残差</th><th>最大 dp / 终态 p</th><th>出口 Q</th></tr></thead><tbody>` +
    rows.map((r) => `<tr><td>${r.name}</td><td>${r.engine}</td><td>${r.ok}</td>` +
      `<td>${r.iters}</td><td>${r.resid}</td><td>${r.a}</td><td>${r.b}</td></tr>`).join('') +
    '</tbody></table>';
  $('psStatus').textContent = `批量求解完成：${rows.length} 个参数集。`;
}

$('psAddBtn').addEventListener('click', () => {
  const ps = {
    id: nextParamSetId(),
    name: `参数集 ${state.paramSets.length + 1}`,
    description: '',
    engine: $('sweepEngine') ? $('sweepEngine').value : 'steady',
    overrides: [],
  };
  state.paramSets.push(ps);
  state.activeParamSet = ps.id;
  renderParamSets();
});
$('psRunAllBtn').addEventListener('click', runAllParamSets);
$('psExportBtn').addEventListener('click', async () => {
  try {
    await invoke('write_text_file', {
      path: '/tmp/chflow_paramsets.json',
      content: paramSetsToJson(state.paramSets, state.activeParamSet),
    });
    $('psStatus').textContent = '参数集已导出 /tmp/chflow_paramsets.json。';
  } catch (e) {
    $('psStatus').textContent = `导出失败：${e}`;
  }
});
$('psImportBtn').addEventListener('click', async () => {
  try {
    const text = await invoke('read_text_file', { path: '/tmp/chflow_paramsets.json' });
    const data = parseParamSets(text);
    state.paramSets = data.sets;
    state.activeParamSet = data.activeId;
    renderParamSets();
    $('psStatus').textContent = `已导入 ${data.sets.length} 个参数集。`;
  } catch (e) {
    $('psStatus').textContent = `导入失败：${e}`;
  }
});
$('psList').addEventListener('click', (e) => {
  const item = e.target.closest('.scen-item');
  if (!item) return;
  const ps = state.paramSets.find((s) => s.id === item.dataset.id);
  if (!ps) return;
  const act = e.target.dataset.act;
  if (act === 'del') {
    state.paramSets = state.paramSets.filter((s) => s.id !== ps.id);
    if (state.activeParamSet === ps.id) state.activeParamSet = null;
    renderParamSets();
  } else if (act === 'apply') {
    applyParamSet(ps);
  } else {
    state.activeParamSet = ps.id;
    renderParamSets();
  }
});
$('psEdit').addEventListener('change', (e) => {
  const ps = state.paramSets.find((s) => s.id === state.activeParamSet);
  if (!ps) return;
  const el = e.target;
  if (el.id === 'psName') ps.name = el.value;
  else if (el.id === 'psDesc') ps.description = el.value;
  else if (el.id === 'psEngine') { ps.engine = el.value; renderParamSetEdit(); return; }
  else if (el.dataset.i !== undefined) {
    const o = ps.overrides[Number(el.dataset.i)];
    if (!o) return;
    if (el.dataset.f === 'value') o.value = Number(el.value);
    else if (el.dataset.f === 'param') { o.param = el.value; o.target = ''; renderParamSetEdit(); return; }
    else if (el.dataset.f === 'target') o.target = el.value;
  }
  renderParamSets();
});
$('psEdit').addEventListener('click', (e) => {
  const ps = state.paramSets.find((s) => s.id === state.activeParamSet);
  if (!ps) return;
  if (e.target.id === 'psAddOverride') {
    const params = ps.engine === 'steady' ? STEADY_SCEN_PARAMS : TRANSIENT_WHAT_IF_PARAMS;
    ps.overrides.push({ param: params[0].value, target: '', value: 0 });
    renderParamSetEdit();
  } else if (e.target.dataset.act === 'rm') {
    ps.overrides.splice(Number(e.target.dataset.i), 1);
    renderParamSetEdit();
  } else if (e.target.id === 'psApplyBtn') {
    applyParamSet(ps);
  } else if (e.target.id === 'psMergeBtn') {
    const row = {
      param: $('sweepParam').value,
      target: $('sweepTarget').value,
      value: tornadoCurrentValue($('sweepParam').value, $('sweepTarget').value),
    };
    if (!Number.isFinite(row.value)) {
      $('psStatus').textContent = '当前扫描参数在模型中无对应对象，无法并入。';
      return;
    }
    ps.overrides = mergeOverrides(ps.overrides, [row]);
    renderParamSetEdit();
    renderParamSets();
  }
});

// ---------------------------------------------------------------------------
// F17 在线监控面板（ITER / MONITOR 增量流）
// ---------------------------------------------------------------------------

const MONITOR_PALETTE2 = ['#4bb3c4', '#e0a33e', '#8fd3a8', '#ff8f6b', '#9b8ad1', '#c96a9b', '#7db86a'];

function renderMonitorPanel() {
  const st = state.monitorStream;
  const el = $('monitorChart');
  const nMon = Object.keys(st.monitors).length;
  $('monitorStatus').textContent =
    `iter=${st.lastIter} · 残差采样=${st.iters.length} · monitor 曲线=${nMon}` +
    ` · 乱序拒绝=${st.rejected} · 坏行=${st.malformed}`;
  if (!st.iters.length && !nMon) {
    el.innerHTML = '<p class="empty">暂无 ITER/MONITOR 行：运行任意求解后此面板按流增量解析逐迭代监控。</p>';
    $('monitorLegend').textContent = '';
    return;
  }
  const W = 560;
  const H = 160;
  const { commands, seriesNames } = monitorDrawCommands(st, W, H);
  let svg = `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="none">`;
  let idx = 0;
  for (const c of commands) {
    if (c.role === 'tick') {
      svg += `<line x1="${c.x1}" y1="${c.y1.toFixed(1)}" x2="${c.x2}" y2="${c.y2.toFixed(1)}" `
        + 'stroke="var(--muted)" stroke-dasharray="2 3" opacity="0.5"/>';
      svg += `<text x="2" y="${(c.y1 - 2).toFixed(1)}" font-size="8" fill="var(--muted)">`
        + `${c.value.toExponential(0)}</text>`;
    } else if (c.type === 'polyline') {
      const color = MONITOR_PALETTE2[idx % MONITOR_PALETTE2.length];
      idx += 1;
      svg += `<polyline points="${c.points}" fill="none" stroke="${color}" stroke-width="1.4"/>`;
    }
  }
  svg += '</svg>';
  el.innerHTML = svg;
  let j = 0;
  $('monitorLegend').innerHTML = seriesNames.map((n) => {
    const color = MONITOR_PALETTE2[j++ % MONITOR_PALETTE2.length];
    return `<span style="color:${color}">■</span> ${n}`;
  }).join(' · ');
}

$('monitorClearBtn').addEventListener('click', () => {
  state.monitorStream = createMonitorState(1500);
  renderMonitorPanel();
});

// ---------------------------------------------------------------------------
// project save / load / new
// ---------------------------------------------------------------------------

function projectJson() {
  return JSON.stringify({
    version: 1,
    kind: 'chuanflow-project',
    model: state.model,
    units: state.units,
    theme: state.theme,
    paramSets: state.paramSets,
    activeParamSet: state.activeParamSet,
  }, null, 2);
}

async function saveProject() {
  const path = $('projPath').value.trim();
  if (!path) {
    setStatus('请先填写工程文件路径。', 'error');
    return;
  }
  try {
    await invoke('write_text_file', { path, content: projectJson() });
    setStatus(`工程已保存：${path}`, 'ok');
    appendLog(`[project] 已保存 ${path}`);
  } catch (e) {
    setStatus(`保存工程失败: ${e}`, 'error');
  }
}

async function openProject() {
  const path = $('projPath').value.trim();
  if (!path) {
    setStatus('请先填写工程文件路径。', 'error');
    return;
  }
  try {
    const text = await invoke('read_text_file', { path });
    const data = JSON.parse(text);
    if (data.kind !== 'chuanflow-project') throw new Error('不是 ChuanFlow 工程文件');
    state.model = { ...emptyModel(), ...data.model };
    state.sel = null;
    state.steady = null;
    state.transient = null;
    state.coupled = null;
    state.coupledEdge = null;
    if (data.units) applyUnits(data.units);
    if (data.theme) applyTheme(data.theme);
    state.paramSets = Array.isArray(data.paramSets) ? data.paramSets : [];
    state.activeParamSet = data.activeParamSet ?? null;
    afterModelChange();
    setStatus(`工程已载入：${path}`, 'ok');
    appendLog(`[project] 已载入 ${path}`);
  } catch (e) {
    setStatus(`打开工程失败: ${e}`, 'error');
  }
}

function newProject() {
  state.paramSets = [];
  state.activeParamSet = null;
  freshModel(DEFAULT_CASE);
  setStatus('已新建工程（示例算例）。', 'ok');
}

$('tbNew').addEventListener('click', newProject);
$('tbSave').addEventListener('click', saveProject);
$('menuNewProject').addEventListener('click', newProject);
$('menuSaveProject').addEventListener('click', saveProject);
$('menuOpenProject').addEventListener('click', openProject);
$('menuExportDict').addEventListener('click', () => switchTab('dict'));
$('menuRunSteady').addEventListener('click', () => $('tbRunSteady').click());
$('menuRunTransient').addEventListener('click', () => $('tbRunTransient').click());
$('menuRunCoupled').addEventListener('click', () => $('runCoupledBtn').click());
$('menuStop').addEventListener('click', () => $('tbStop').click());

// ---------------------------------------------------------------------------
// toolbar run buttons
// ---------------------------------------------------------------------------

$('tbRunSteady').addEventListener('click', () => {
  const text = state.rawDictText || networkToCaseText(state.model);
  // Audit fix (P1-2): the steady network app silently ignores the E4 device
  // blocks — point the user at ▶点元 instead of dropping the physics.
  if (POINT_ELEMENT_BLOCK_RE.test(text)) {
    appendLog('[steady] 检测到点元设备块（jetPumps/pcpPumps/plungerGasLifts/separators）：稳态求解不消费其物理，请用 ▶点元 运行。');
    setStatus('检测到点元设备块：稳态求解不含设备效应，请使用 ▶点元。');
  }
  startJob('steady', text);
});

// Audit fix (P1-2): point-element run — same serialized dictionary, device
// blocks consumed by chflow-1d-point-elements.
$('tbRunPointElements')?.addEventListener('click', () => {
  const text = state.rawDictText || networkToCaseText(state.model);
  if (!POINT_ELEMENT_BLOCK_RE.test(text)) {
    setStatus('当前画布没有点元设备（射流泵/螺杆泵/柱塞气举/分离器），无需运行点元引擎。', 'error');
    return;
  }
  startJob('pointElements', text);
});
$('menuRunPointElements')?.addEventListener('click', () => $('tbRunPointElements').click());
$('tbRunTransient').addEventListener('click', () => {
  const text = state.rawDictText || networkToCaseText(state.model);
  startJob('transient', text);
});
$('tbRunCoupled').addEventListener('click', () => $('runCoupledBtn').click());
$('tbRun3d').addEventListener('click', () => startJob('3d', $('caseText').value));
$('menuRun3d').addEventListener('click', () => startJob('3d', $('caseText').value));

// ---------------------------------------------------------------------------
// model change plumbing
// ---------------------------------------------------------------------------

function afterModelChange(syncDict = true) {
  // a manual model edit supersedes the raw dict: runs serialize from the model
  state.rawDictText = null;
  if (syncDict) syncDictText();
  renderTree();
  renderProps();
  renderCanvas();
  renderCharts();
  renderSweepTargets();
}

// ---------------------------------------------------------------------------
// startup
// ---------------------------------------------------------------------------

caseText.value = DEFAULT_CASE;

// E2: initialize the flow-assurance form from the case dictionary now that
// it carries its default text (the grid was built earlier at module load).
renderFlowAssuranceForm(parseFlowAssurance(caseText.value));

// Audit fix (P2-6): the schedule panel must initialize from the case dictionary
// too — an empty form applied by click would removeDictBlock the dictionary's
// existing sourceSchedule{} silently. Mirrors the E2 flow-assurance startup load.
schedModel = parseSchedules(caseText.value);
renderSchedulePanel();

reflectEngines();

(async () => {
  try {
    const paths = await invoke('default_engines');
    if (Array.isArray(paths) && paths.length === 8 && paths[0] && paths[1]) {
      engines.steady = paths[0];
      engines.transient = paths[1];
      engines.coupled = paths[2];
      engines.coupledCase = paths[3];
      engines.twofluid = paths[4];
      engines.nodal = paths[5];
      engines.pvt = paths[6];
      if (paths[7]) engines.threeD = paths[7];
      persistEngines();
      reflectEngines();
    } else {
      console.warn('default_engines returned an unexpected payload', paths);
      setStatus('引擎路径解析返回异常（保留设置或开发树默认值），可在 case 字典页手动核对求解器路径。');
    }
  } catch (err) {
    // E9: never silently swallow this — the user must know the packaged
    // layout did not resolve and the developer-tree fallback is in use.
    console.warn('default_engines failed:', err);
    setStatus('未解析到随包引擎路径，已回退到设置/开发树默认值；如求解失败请在 case 字典页核对路径。');
  }
  // E1 (§5.181): resolve the 3D engine family (binaries + case paths). The
  // Rust resolve beats the hardcoded developer-tree defaults, same contract
  // as the 8 pre-existing keys; an older backend without the command keeps
  // the defaults quietly.
  try {
    const fam = await invoke('three_d_engines') || [];
    for (const e of fam) {
      if (THREE_D_FAMILY_KEYS.includes(e.key) && e.exe) engines[e.key] = e.exe;
    }
  } catch (err) {
    console.warn('three_d_engines failed:', err);
  }
  // Coupling family binaries + case paths from the Rust `coupled_engines`
  // resolve; an older backend without the command keeps the dev-tree defaults.
  try {
    const fam = await invoke('coupled_engines') || [];
    for (const e of fam) {
      if (COUPLED_FAMILY_KEYS.includes(e.key) && e.exe) engines[e.key] = e.exe;
    }
  } catch (err) {
    console.warn('coupled_engines failed:', err);
  }
  populate3DEngineSelect();
  $('threeDExePath').value = engines[active3DEngine] || engines.threeD;
  persistEngines();
  reflectEngines();
  renderTemplateMenu();
  reconcileRegimeBoundaries();
})();

// E1 (§5.181): fill the dictionary-page 3D family selector from the static
// registry (labels), so the family is usable before the backend resolves.
function populate3DEngineSelect() {
  const sel = $('threeDEngineSel');
  if (!sel || sel.options.length) return;
  for (const e of THREE_D_FAMILY) {
    const opt = document.createElement('option');
    opt.value = e.key;
    opt.textContent = `${e.label} (${e.bin})`;
    sel.appendChild(opt);
  }
  sel.value = active3DEngine;
}

$('threeDEngineSel').addEventListener('change', (e) => {
  active3DEngine = e.target.value;
  $('threeDExePath').value = engines[active3DEngine] || engines.threeD;
  const famLabel = (THREE_D_FAMILY.find((en) => en.key === active3DEngine) || {}).label || active3DEngine;
  setStatus(`3D 求解器已切换到族成员：${famLabel}。载入对应模板后点「▶ 3D」求解。`, 'ok');
});

/**
 * E5: reconcile the frontend's hardcoded Beggs-Brill boundary copy against
 * the Rust port of the kernel formula (which itself is reconciled against
 * the C++ engine report by a unit test). A mismatch means someone changed
 * one side without the other — surface it instead of failing silently.
 */
async function reconcileRegimeBoundaries() {
  const probes = [0.05, 0.2, 0.386561092748, 0.6, 0.9];
  try {
    const rows = await Promise.all(
      probes.map((lam) => invoke('beggs_brill_boundaries', { lambda: lam })),
    );
    let bad = 0;
    rows.forEach((r, i) => {
      const lam = probes[i];
      const pat = beggsBrillPattern(lam, 0);
      const want = [pat.l1, pat.l2, pat.l3, pat.l4];
      const got = Array.isArray(r) ? r : [r[0], r[1], r[2], r[3]];
      got.forEach((v, k) => {
        const scale = Math.abs(want[k]) || 1;
        if (Math.abs(v - want[k]) > 1e-9 * scale) bad++;
      });
    });
    if (bad) {
      console.error(`regime boundary reconciliation failed on ${bad} values`);
      setStatus('流型图边界与内核公式失配（E5 对账失败），请在 regimemap.js 与 beggsBrill.cpp 之间核对。', 'error');
    }
  } catch (_) {
    /* backend without the command (older binary) -> skip quietly */
  }
}

['exePath', 'transientExePath', 'coupledExePath', 'twofluidExePath', 'nodalExePath', 'pvtExePath', 'threeDExePath', 'coupledCasePath'].forEach((id) => {
  $(id).addEventListener('change', (e) => {
    const v = e.target.value.trim();
    if (id === 'exePath') engines.steady = v;
    if (id === 'transientExePath') engines.transient = v;
    if (id === 'coupledExePath') {
      // edits the active coupling family member (mirrors the 3D exe box)
      engines[activeCoupledEngine] = v;
    }
    if (id === 'twofluidExePath') engines.twofluid = v;
    if (id === 'nodalExePath') engines.nodal = v;
    if (id === 'pvtExePath') engines.pvt = v;
    if (id === 'threeDExePath') {
      // E1 (§5.181): the box edits the *active* family member's binary
      engines[active3DEngine] = v;
    }
    if (id === 'coupledCasePath') engines.coupledCase = v;
    persistEngines();
  });
});

freshModel(DEFAULT_CASE);
applyUnits(UNIT_SYSTEMS.SI);
renderSweepTargets(false);
renderCharts();
renderConvergence();
setStatus('就绪。左侧项目树选择元素，工具栏运行稳态/瞬态/耦合求解；文件菜单载入算例模板。');
