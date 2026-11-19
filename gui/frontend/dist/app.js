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
  extractCaseDetails,
} from './workspace.js';

// ---------------------------------------------------------------------------
// state
// ---------------------------------------------------------------------------

const DEFAULT_EXE = '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-steady';
const DEFAULT_TRANSIENT_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-transient';
const DEFAULT_COUPLED_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/coupling/apps/chflow-coupled-1d3d';
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
  job: { id: null, label: null, running: false, offset: 0, full: '' },
  collapsedLeft: false,
  collapsedRight: false,
};

const engines = {
  steady: DEFAULT_EXE,
  transient: DEFAULT_TRANSIENT_EXE,
  coupled: DEFAULT_COUPLED_EXE,
  coupledCase: DEFAULT_COUPLED_CASE,
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
      : `${p.name}`;
    svg.appendChild(label);
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

async function startJob(label, caseText_) {
  if (state.job.running) {
    setStatus('已有作业在运行，请先停止。', 'error');
    return;
  }
  const exe =
    label === 'steady' ? engines.steady
    : label === 'transient' ? engines.transient
    : engines.coupled;
  let casePath;
  if (label === 'coupled') {
    casePath = await prepareCoupledCase();
    if (!casePath) return;
  } else {
    casePath = await writeCaseFile(caseText_, label);
  }
  try {
    const id = await invoke('start_case', { label, exePath: exe, casePath });
    state.job = { id, label, running: true, offset: 0, full: '' };
    setJobChip(true, label);
    appendLog(`[job ${id}] 启动 ${label}: ${exe} ${casePath}`);
    setStatus(`作业 ${label} 已启动。`);
    switchTab(label === 'coupled' ? 'coupled' : 'canvas');
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

  const dpr = window.devicePixelRatio || 1;
  cv.width = cv.clientWidth * dpr;
  cv.height = cv.clientHeight * dpr;
  const ctx = cv.getContext('2d');
  ctx.clearRect(0, 0, cv.width, cv.height);

  // x-y projection bounds over used vertices
  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  for (const c of cells) {
    for (const v of c) {
      xmin = Math.min(xmin, points[v].x);
      xmax = Math.max(xmax, points[v].x);
      ymin = Math.min(ymin, points[v].y);
      ymax = Math.max(ymax, points[v].y);
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
    const vids = cells[k];
    if (vids.length < 3) return;
    ctx.beginPath();
    ctx.moveTo(px(points[vids[0]].x), py(points[vids[0]].y));
    for (let i = 1; i < vids.length; i += 1) {
      ctx.lineTo(px(points[vids[i]].x), py(points[vids[i]].y));
    }
    ctx.closePath();
    ctx.fillStyle = rainbowColor(val, lo, hi);
    ctx.fill();
    ctx.strokeStyle = 'rgba(0,0,0,0.15)';
    ctx.lineWidth = Math.max(0.5, 0.5 * dpr);
    ctx.stroke();
  });

  const unit = state.fieldMode === 'p' ? unitLabel('pressure', state.units) : 'm/s';
  const conv = state.fieldMode === 'p'
    ? (v) => toDisplay(v, 'pressure', state.units)
    : (v) => v;
  $('legendLo').textContent = `${Number(conv(lo).toPrecision(5))} ${unit}`;
  $('legendHi').textContent = `${Number(conv(hi).toPrecision(5))} ${unit}`;
  $('fieldInfo').textContent =
    `块网格 ${cells.length} 单元 / ${points.length} 顶点 · ` +
    `渲染场 ${state.fieldMode === 'p' ? '压力 p' : '流速 |U|'} · x-y 投影`;
}

$('loadVtuBtn').addEventListener('click', loadVtu);
$('fieldMode').addEventListener('change', (e) => {
  state.fieldMode = e.target.value;
  renderField();
});

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
}

// ---------------------------------------------------------------------------
// startup
// ---------------------------------------------------------------------------

caseText.value = DEFAULT_CASE;
$('exePath').value = engines.steady;
$('transientExePath').value = engines.transient;
$('coupledExePath').value = engines.coupled;
$('coupledCasePath').value = engines.coupledCase;

(async () => {
  try {
    const paths = await invoke('default_engines');
    if (Array.isArray(paths) && paths.length === 4 && paths[0] && paths[1]) {
      engines.steady = paths[0];
      engines.transient = paths[1];
      engines.coupled = paths[2];
      engines.coupledCase = paths[3];
      $('exePath').value = paths[0];
      $('transientExePath').value = paths[1];
      $('coupledExePath').value = paths[2];
      $('coupledCasePath').value = paths[3];
    }
  } catch (_) {
    /* keep developer defaults */
  }
})();

['exePath', 'transientExePath', 'coupledExePath', 'coupledCasePath'].forEach((id) => {
  $(id).addEventListener('change', (e) => {
    const v = e.target.value.trim();
    if (id === 'exePath') engines.steady = v;
    if (id === 'transientExePath') engines.transient = v;
    if (id === 'coupledExePath') engines.coupled = v;
    if (id === 'coupledCasePath') engines.coupledCase = v;
  });
});

freshModel(DEFAULT_CASE);
applyUnits(UNIT_SYSTEMS.SI);
renderCharts();
renderConvergence();
setStatus('就绪。左侧项目树选择元素，工具栏运行稳态/瞬态/耦合求解。');
