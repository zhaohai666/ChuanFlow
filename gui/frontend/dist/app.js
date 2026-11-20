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
  clipKeepsByPlane, adaptiveVectorIndices,
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
  extractCaseDetails, parseTwoFluidFrames,
} from './workspace.js';
import {
  templateTarget, isPartialFidelity,
} from './templates.js';

// ---------------------------------------------------------------------------
// state
// ---------------------------------------------------------------------------

const DEFAULT_EXE = '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-steady';
const DEFAULT_TRANSIENT_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-transient';
const DEFAULT_COUPLED_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/coupling/apps/chflow-coupled-1d3d';
const DEFAULT_TWOFLUID_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-twofluid';
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
  transientSeries: 'p:A',
  coupled: null,
  coupledEdge: null,
  vtu: null,
  colorMode: 'none',
  fieldMode: 'p',
  view3d: { yaw: 0, pitch: 0, clipAxis: 'none', clipFrac: 0.5, clipNormal: [1, 0, 0], vectors: false },
  slug: { frames: [], idx: 0, playing: false, timer: null },
  sweep: { rows: [], running: false, out: '' },
  loadedTemplate: null,
  job: { id: null, label: null, running: false, offset: 0, full: '' },
  collapsedLeft: false,
  collapsedRight: false,
};

const engines = {
  steady: DEFAULT_EXE,
  transient: DEFAULT_TRANSIENT_EXE,
  coupled: DEFAULT_COUPLED_EXE,
  coupledCase: DEFAULT_COUPLED_CASE,
  twofluid: DEFAULT_TWOFLUID_EXE,
};

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
  const NS = 'http://www.w3.org/2000/svg';
  const g = document.createElementNS(NS, 'g');
  g.dataset.pipeId = p.id;
  const title = document.createElementNS(NS, 'title');
  const sel = () => {
    state.sel = { kind: 'pipe', ref: p.id };
    state.model.selection = { kind: 'pipe', id: p.id };
    renderCanvas();
    renderProps();
    renderTree();
  };

  if (p.liftDomePressure) {
    // gas-lift valve (§5.89): dome charge + port, drawn as a bell with stem
    g.setAttribute('class', 'dev-gaslift');
    const bell = document.createElementNS(NS, 'path');
    bell.setAttribute(
      'd',
      `M ${mx - 8},${my + 6} Q ${mx - 8},${my - 9} ${mx},${my - 9} ` +
      `Q ${mx + 8},${my - 9} ${mx + 8},${my + 6} Z`,
    );
    const port = document.createElementNS(NS, 'circle');
    port.setAttribute('cx', mx);
    port.setAttribute('cy', my + 8);
    port.setAttribute('r', 2.5);
    g.append(bell, port);
    const bits = [`气举阀 ${p.name}: dome=${p.liftDomePressure}`];
    if (p.liftAreaRatio) bits.push(`Ar=${p.liftAreaRatio}`);
    if (p.liftSpringPressure) bits.push(`spring=${p.liftSpringPressure}`);
    if (p.liftDischargeCoeff) bits.push(`Cd=${p.liftDischargeCoeff}`);
    if (p.liftKappa) bits.push(`kappa=${p.liftKappa}`);
    title.textContent = bits.join(' ');
  } else if (p.rodPlungerDiameter) {
    // rod pump (§5.90): barrel + plunger, drawn as a piston jack
    g.setAttribute('class', 'dev-rodpump');
    const barrel = document.createElementNS(NS, 'rect');
    barrel.setAttribute('x', mx - 10);
    barrel.setAttribute('y', my - 5);
    barrel.setAttribute('width', 20);
    barrel.setAttribute('height', 10);
    const horse = document.createElementNS(NS, 'path');
    horse.setAttribute(
      'd',
      `M ${mx},${my - 5} L ${mx},${my - 12} L ${mx + 7},${my - 12}`,
    );
    g.append(barrel, horse);
    const bits = [`有杆泵 ${p.name}: bore=${p.rodPlungerDiameter}`];
    if (p.rodStrokeLength) bits.push(`stroke=${p.rodStrokeLength}`);
    if (p.rodStrokesPerMinute) bits.push(`spm=${p.rodStrokesPerMinute}`);
    if (p.rodDiameter) bits.push(`rodD=${p.rodDiameter}`);
    if (p.rodLength) bits.push(`rodL=${p.rodLength}`);
    title.textContent = bits.join(' ');
  } else if (p.valve) {
    g.setAttribute('class', 'dev-valve');
    const t1 = document.createElementNS(NS, 'polygon');
    t1.setAttribute('points', `${mx - 7},${my - 7} ${mx - 7},${my + 7} ${mx},${my}`);
    const t2 = document.createElementNS(NS, 'polygon');
    t2.setAttribute('points', `${mx + 7},${my - 7} ${mx + 7},${my + 7} ${mx},${my}`);
    const stem = document.createElementNS(NS, 'line');
    stem.setAttribute('x1', mx);
    stem.setAttribute('y1', my - 7);
    stem.setAttribute('x2', mx);
    stem.setAttribute('y2', my - 12);
    g.append(t1, t2, stem);
    title.textContent = `节流阀 ${p.name}: valve K=${p.valve}`;
  } else if (p.pumpHead || p.pumpRate) {
    g.setAttribute('class', 'dev-pump');
    const c = document.createElementNS(NS, 'circle');
    c.setAttribute('cx', mx);
    c.setAttribute('cy', my);
    c.setAttribute('r', 9);
    const dx = bx - ax;
    const dy = by - ay;
    const len = Math.hypot(dx, dy) || 1;
    const ux = dx / len;
    const uy = dy / len;
    const tri = document.createElementNS(NS, 'polygon');
    tri.setAttribute(
      'points',
      `${mx + ux * 6},${my + uy * 6} ${mx - ux * 4 - uy * 4},${my - uy * 4 + ux * 4} ${mx - ux * 4 + uy * 4},${my - uy * 4 - ux * 4}`,
    );
    g.append(c, tri);
    const bits = [`泵 ${p.name}:`];
    if (p.pumpHead) bits.push(`pumpHead=${p.pumpHead}`);
    if (p.pumpRate) bits.push(`pumpRate=${p.pumpRate}`);
    if (p.pumpShutoffHead) bits.push(`shutoff=${p.pumpShutoffHead}`);
    if (p.pumpStages) bits.push(`stages=${p.pumpStages}`);
    if (p.pumpSpeed) bits.push(`speed=${p.pumpSpeed}`);
    if (p.pumpGasFraction) bits.push(`gasFrac=${p.pumpGasFraction}`);
    if (p.pumpIntakePressure) bits.push(`intake=${p.pumpIntakePressure}`);
    title.textContent = bits.join(' ');
  }
  g.appendChild(title);
  g.addEventListener('click', sel);
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

    if (p.valve || p.pumpHead || p.pumpRate) {
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
$('removeDeviceBtn').addEventListener('click', () => {
  withSelectedPipe((p) => {
    state.model = setPipeFields(state.model, p.id, {
      valve: 0, pumpHead: 0, pumpRate: 0,
      pumpShutoffHead: 0, pumpStages: 0, pumpSpeed: 0,
      pumpGasFraction: 0, pumpIntakePressure: 0,
      liftDomePressure: 0, liftAreaRatio: 0, liftSpringPressure: 0,
      liftDischargeCoeff: 0, liftKappa: 0,
      rodPlungerDiameter: 0, rodStrokeLength: 0, rodStrokesPerMinute: 0,
      rodDiameter: 0, rodLength: 0, rodYoungModulus: 0, rodDensity: 0,
    });
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
        valve: geo.valve ?? 0,
        pumpHead: geo.pumpHead ?? 0,
        pumpRate: geo.pumpRate ?? 0,
        pumpShutoffHead: geo.pumpShutoffHead ?? 0,
        pumpStages: geo.pumpStages ?? 0,
        pumpSpeed: geo.pumpSpeed ?? 0,
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
    : label === 'transient' ? engines.transient
    : label === 'slug' ? engines.twofluid
    : engines.coupled;
  let casePath;
  if (label === 'coupled') {
    casePath = await prepareCoupledCase();
    if (!casePath) return;
  } else if (label === 'slug') {
    casePath = isPath ? caseText_ : await writeCaseFile(caseText_, label);
  } else {
    casePath = await writeCaseFile(caseText_, label);
  }
  try {
    const id = await invoke('start_case', { label, exePath: exe, casePath });
    state.job = { id, label, running: true, offset: 0, full: '' };
    setJobChip(true, label);
    appendLog(`[job ${id}] 启动 ${label}: ${exe} ${casePath}`);
    setStatus(`作业 ${label} 已启动。`);
    switchTab(label === 'coupled' ? 'coupled' : label === 'slug' ? 'slug' : 'canvas');
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
  }
  if (st.finished) {
    const label = state.job.label;
    const full = state.job.full;
    setJobChip(false);
    state.job = { id: null, label: null, running: false, offset: 0, full: '' };
    if (st.exit === 0) {
      appendLog(`[job] ${label} 完成`);
      if (label === 'steady') onSteadyDone(full);
      else if (label === 'transient') onTransientDone(full);
      else if (label === 'coupled') onCoupledDone(full);
      else if (label === 'slug') onSlugDone(full);
    } else {
      appendLog(`[job] ${label} 退出码 ${st.exit}`);
      setStatus(`${label} 求解失败（退出码 ${st.exit}），详见求解日志。`, 'error');
    }
  }
}, 300);

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

function onTransientDone(fullOutput) {
  const t = parseTransientOutput(fullOutput);
  if (!t.samples.length) {
    setStatus('瞬态输出解析失败（日志无时间帧）。', 'error');
    return;
  }
  state.transient = t;
  renderCharts();
  renderConvergence();
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
  const { points, cells, cellU, cellP } = vtu;
  const values = state.fieldMode === 'p'
    ? cellP
    : cellU.map((u) => Math.hypot(u[0], u[1], u[2]));
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
    const vecColor = state.theme === 'dark' ? '#e6edf6' : '#1c2733';
    ctx.strokeStyle = vecColor;
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

  const unit = state.fieldMode === 'p' ? unitLabel('pressure', state.units) : 'm/s';
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
    `渲染场 ${state.fieldMode === 'p' ? '压力 p' : '流速 |U|'} · ` +
    `yaw ${v3.yaw}° pitch ${v3.pitch}° · ${clipTxt}` +
    (v3.vectors ? ` · 矢量 ${drawn}` : '');
}

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
  const padL = 50 * dpr;
  const padR = 30 * dpr;
  const padT = 30 * dpr;
  const padB = 40 * dpr;
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

  // axis labels every ~10 cells
  ctx.fillStyle = dark ? '#9aa7b8' : '#5b6b7d';
  ctx.font = `${10 * dpr}px sans-serif`;
  const step = Math.max(1, Math.round(n / 10));
  for (let i = 0; i < n; i += step) {
    const c = cells[i];
    if (!c) continue;
    ctx.fillText(String(c.index ?? i), padL + i * bw, cv.height - padB + 14 * dpr);
  }
  $('slugReadout').textContent =
    `t=${fr.t} · vm=${fr.vm} · dp=${fr.dp} · 帧 ${s.idx + 1}/${s.frames.length}`;
}

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

function onSlugDone(fullOutput) {
  const frames = parseTwoFluidFrames(fullOutput);
  if (!frames.length) {
    setStatus('段塞输出解析失败（日志无 t/vm/dp 帧行）。', 'error');
    return;
  }
  const s = state.slug;
  s.frames = frames;
  s.idx = 0;
  const slider = $('slugFrame');
  slider.max = String(frames.length - 1);
  slider.value = '0';
  slider.disabled = false;
  $('slugPlayBtn').disabled = false;
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
// parameter sweep (steady, one run_case_text per value)
// ---------------------------------------------------------------------------

function sweepTargetOptions() {
  const param = $('sweepParam').value;
  if (param === 'sourceRate') {
    return state.model.sources.map((s) => s.name);
  }
  return state.model.pipes.map((p) => p.name);
}

function renderSweepTargets(keep = true) {
  const sel = $('sweepTarget');
  const prev = keep ? sel.value : null;
  const opts = sweepTargetOptions();
  sel.innerHTML = opts.map((o) => `<option value="${o}">${o}</option>`).join('');
  if (prev && opts.includes(prev)) sel.value = prev;
}

$('sweepParam').addEventListener('change', () => renderSweepTargets(false));

async function runSweep() {
  if (state.sweep.running) {
    setStatus('扫描进行中，请等待完成。', 'error');
    return;
  }
  if (state.job.running) {
    setStatus('已有求解作业在运行，请先停止。', 'error');
    return;
  }
  if (isPartialFidelity(state.loadedTemplate)) {
    setStatus('当前模板为黑油/PVT 偏保真算例：编辑器模型无法完整复刻其字典，扫描已禁用。', 'error');
    return;
  }
  const param = $('sweepParam').value;
  const target = $('sweepTarget').value;
  const values = $('sweepValues').value
    .split(',')
    .map((s) => Number(s.trim()))
    .filter((v) => Number.isFinite(v));
  if (!target) {
    setStatus('扫描作用对象为空：请先建好模型（源或管段）。', 'error');
    return;
  }
  if (!values.length) {
    setStatus('取值列表为空或非数值。', 'error');
    return;
  }

  state.sweep = { rows: [], running: true, out: state.sweep.out };
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
  const vLo = Math.min(...values.filter(Number.isFinite), 0);
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

$('sweepRunBtn').addEventListener('click', runSweep);
$('sweepOut').addEventListener('change', (e) => {
  state.sweep.out = e.target.value;
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
  if (target.kind === 'coupled') {
    $('coupledCasePath').value = entry.path;
    switchTab('coupled');
    setStatus(`已载入模板「${entry.name}」：点「▶ 运行耦合」求解。`, 'ok');
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
    const b = document.createElement('button');
    b.textContent = (isPartialFidelity(t.id) ? '◇ ' : '◆ ') + t.name;
    b.title = t.path;
    b.addEventListener('click', () => loadTemplate(t));
    menu.appendChild(b);
  }
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

  // transient series select
  const t = state.transient;
  if (t && t.samples.length) {
    const opts = [];
    for (const name of t.nodeNames) {
      opts.push(`p:${name}`, `alpha:${name}`, `hydrate:${name}`);
    }
    for (const name of t.pipeNames) opts.push(`Q:${name}`);
    const sel = $('seriesSelect');
    if (sel.dataset.opts !== opts.join(',')) {
      sel.innerHTML = opts.map((o) => `<option value="${o}">${o}</option>`).join('');
      sel.dataset.opts = opts.join(',');
      state.transientSeries = opts[0];
      sel.value = opts[0];
    }
    if (!opts.includes(state.transientSeries)) state.transientSeries = opts[0];

    $('timeSlider').max = String(t.samples.length - 1);
    const idx = Math.min(Number($('timeSlider').value) || 0, t.samples.length - 1);
    const frame = t.samples[idx];
    $('timeReadout').textContent = `t=${frame.t}`;

    const [kind, name] = state.transientSeries.split(':');
    const getter =
      kind === 'p' ? (s) => s.pressure[name]
      : kind === 'alpha' ? (s) => s.alpha[name]
      : kind === 'hydrate' ? (s) => s.hydrate[name]
      : (s) => s.flow[name];
    const pts = toTimeSeriesPoints(t.samples, getter, 560, 110);
    $('transientChart').innerHTML =
      `<svg class="chart" viewBox="0 0 600 120" preserveAspectRatio="none">` +
      `<polyline points="${pts.points}" fill="none" stroke="var(--accent)" stroke-width="1.5"/>` +
      `<circle cx="${(pts.points.split(' ')[idx] || '0,60').split(',')[0]}" cy="${(pts.points.split(' ')[idx] || '0,60').split(',')[1]}" r="3" fill="var(--coupled)"/>` +
      `</svg>` +
      `<p class="hint">${state.transientSeries} · t ${pts.tMin}..${pts.tMax} · 值域 ${pts.min.toPrecision(4)}..${pts.max.toPrecision(4)} ${kind === 'p' ? unitLabel('pressure', u) : kind === 'Q' ? unitLabel('rate', u) : ''}</p>`;
  } else {
    $('transientChart').innerHTML = '<p class="empty">运行瞬态求解后显示时序曲线。</p>';
    $('timeReadout').textContent = '—';
  }
}

$('seriesSelect').addEventListener('change', (e) => {
  state.transientSeries = e.target.value;
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
// project save / load / new
// ---------------------------------------------------------------------------

function projectJson() {
  return JSON.stringify({
    version: 1,
    kind: 'chuanflow-project',
    model: state.model,
    units: state.units,
    theme: state.theme,
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
    afterModelChange();
    setStatus(`工程已载入：${path}`, 'ok');
    appendLog(`[project] 已载入 ${path}`);
  } catch (e) {
    setStatus(`打开工程失败: ${e}`, 'error');
  }
}

function newProject() {
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
  startJob('steady', text);
});
$('tbRunTransient').addEventListener('click', () => {
  const text = state.rawDictText || networkToCaseText(state.model);
  startJob('transient', text);
});
$('tbRunCoupled').addEventListener('click', () => $('runCoupledBtn').click());

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
$('exePath').value = engines.steady;
$('transientExePath').value = engines.transient;
$('coupledExePath').value = engines.coupled;
$('twofluidExePath').value = engines.twofluid;
$('coupledCasePath').value = engines.coupledCase;

(async () => {
  try {
    const paths = await invoke('default_engines');
    if (Array.isArray(paths) && paths.length === 5 && paths[0] && paths[1]) {
      engines.steady = paths[0];
      engines.transient = paths[1];
      engines.coupled = paths[2];
      engines.coupledCase = paths[3];
      engines.twofluid = paths[4];
      $('exePath').value = paths[0];
      $('transientExePath').value = paths[1];
      $('coupledExePath').value = paths[2];
      $('coupledCasePath').value = paths[3];
      $('twofluidExePath').value = paths[4];
    }
  } catch (_) {
    /* keep developer defaults */
  }
  renderTemplateMenu();
})();

['exePath', 'transientExePath', 'coupledExePath', 'twofluidExePath', 'coupledCasePath'].forEach((id) => {
  $(id).addEventListener('change', (e) => {
    const v = e.target.value.trim();
    if (id === 'exePath') engines.steady = v;
    if (id === 'transientExePath') engines.transient = v;
    if (id === 'coupledExePath') engines.coupled = v;
    if (id === 'twofluidExePath') engines.twofluid = v;
    if (id === 'coupledCasePath') engines.coupledCase = v;
  });
});

freshModel(DEFAULT_CASE);
applyUnits(UNIT_SYSTEMS.SI);
renderSweepTargets(false);
renderCharts();
renderConvergence();
setStatus('就绪。左侧项目树选择元素，工具栏运行稳态/瞬态/耦合求解；文件菜单载入算例模板。');
