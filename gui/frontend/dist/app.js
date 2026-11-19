import { invoke } from './tauri-api.js';
import {
  validateCaseText,
  toBarHeights,
  toLinePoints,
  toTimeSeriesPoints,
} from './chart.js';
import {
  emptyModel,
  addNode,
  moveNode,
  removeNode,
  setNodePressure,
  addPipe,
  removePipe,
  addSource,
  removeSource,
  setTime,
  hitTestNode,
  networkToCaseText,
} from './editor.js';

const DEFAULT_EXE = '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-steady';
const DEFAULT_TRANSIENT_EXE =
  '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-transient';

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

const caseText = document.getElementById('caseText');
const exePath = document.getElementById('exePath');
const transientExePath = document.getElementById('transientExePath');
const statusEl = document.getElementById('status');

const W = 480;
const H = 180;
const PAD_L = 34;
const PAD_R = 10;
const PAD_T = 10;
const PAD_B = 30;
const PLOT_W = W - PAD_L - PAD_R;
const PLOT_H = H - PAD_T - PAD_B;

const CW = 640;
const CH = 420;
const canvas = document.getElementById('networkCanvas');
const nodePressureInput = document.getElementById('nodePressure');
const timeDtInput = document.getElementById('timeDt');
const timeEndInput = document.getElementById('timeEnd');
const seriesSelect = document.getElementById('seriesSelect');
const timeSlider = document.getElementById('timeSlider');
const timeReadout = document.getElementById('timeReadout');

let model = emptyModel();
let dragId = null;
let transientSeries = null; // last parsed TransientSeries

function setStatus(msg, kind = '') {
  statusEl.textContent = msg;
  statusEl.className = 'status' + (kind ? ' ' + kind : '');
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;',
  }[c]));
}

function axis() {
  return `<line x1="${PAD_L}" y1="${PAD_T + PLOT_H}" x2="${W - PAD_R}" y2="${PAD_T + PLOT_H}" stroke="#263349"/>`;
}

function renderBars(elId, items, valueKey) {
  const el = document.getElementById(elId);
  if (!items.length) {
    el.innerHTML = '<p class="empty">暂无数据</p>';
    return;
  }
  const heights = toBarHeights(items.map((it) => Number(it[valueKey])), PLOT_H);
  const slot = PLOT_W / items.length;
  const barW = Math.min(slot * 0.6, 42);
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">${axis()}`;
  items.forEach((it, i) => {
    const h = Math.max(0, heights[i]);
    const x = PAD_L + i * slot + (slot - barW) / 2;
    const y = PAD_T + PLOT_H - h;
    const label = escapeHtml(it.name);
    const val = Number(it[valueKey]);
    s += `<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${barW.toFixed(1)}" height="${h.toFixed(1)}" fill="#4ea1ff"><title>${label}: ${val}</title></rect>`;
    s += `<text x="${(x + barW / 2).toFixed(1)}" y="${(PAD_T + PLOT_H + 14).toFixed(1)}" text-anchor="middle" font-size="10" fill="#9aa7b8">${label}</text>`;
    s += `<text x="${(x + barW / 2).toFixed(1)}" y="${(y - 4).toFixed(1)}" text-anchor="middle" font-size="9" fill="#e6edf6">${val.toPrecision(4)}</text>`;
  });
  el.innerHTML = s + '</svg>';
}

function renderLine(elId, items, valueKey) {
  const el = document.getElementById(elId);
  if (!items.length) {
    el.innerHTML = '<p class="empty">暂无数据</p>';
    return;
  }
  const values = items.map((it) => Number(it[valueKey]));
  const points = toLinePoints(values, PLOT_W, PLOT_H)
    .split(' ')
    .map((p) => {
      const [x, y] = p.split(',').map(Number);
      return `${(x + PAD_L).toFixed(2)},${(y + PAD_T).toFixed(2)}`;
    })
    .join(' ');
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">${axis()}`;
  s += `<polyline points="${points}" fill="none" stroke="#4ad991" stroke-width="2"/>`;
  const dots = points.split(' ');
  items.forEach((it, i) => {
    const [cx, cy] = dots[i].split(',').map(Number);
    s += `<circle cx="${cx}" cy="${cy}" r="3" fill="#4ad991"><title>${escapeHtml(it.name)}: ${values[i]}</title></circle>`;
    s += `<text x="${cx.toFixed(1)}" y="${(PAD_T + PLOT_H + 14).toFixed(1)}" text-anchor="middle" font-size="10" fill="#9aa7b8">${escapeHtml(it.name)}</text>`;
  });
  el.innerHTML = s + '</svg>';
}

function renderConvergence(summary) {
  const ok = summary.converged;
  document.getElementById('convergence').innerHTML =
    '<h2>求解状态</h2>' +
    `<p>收敛: <b style="color:${ok ? 'var(--good)' : 'var(--bad)'}">${ok ? '是' : '否'}</b></p>` +
    `<p>迭代次数: <b>${summary.iterations}</b></p>` +
    `<p>相对残差: <b>${Number(summary.residual).toExponential(3)}</b></p>` +
    `<p>节点数: <b>${summary.nodes.length}</b> · 管道数: <b>${summary.pipes.length}</b></p>`;
}

// --- Drag editor -----------------------------------------------------------

function syncTextFromModel() {
  caseText.value = networkToCaseText(model);
}

function renderCanvas() {
  const byName = new Map(model.nodes.map((n) => [n.name, n]));
  let s = '';
  for (const p of model.pipes) {
    const a = byName.get(p.from);
    const b = byName.get(p.to);
    if (!a || !b) continue;
    s += `<line class="edge" x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}"/>`;
  }
  for (const n of model.nodes) {
    const sel = model.selection && model.selection.kind === 'node' && model.selection.id === n.id;
    s += `<g class="node${sel ? ' selected' : ''}" data-id="${n.id}">` +
      `<circle cx="${n.x}" cy="${n.y}" r="12" fill="#1d2940" stroke="#4ea1ff"/>` +
      `<text x="${n.x + 16}" y="${n.y + 4}" fill="#e6edf6">${escapeHtml(n.name)}${n.pressure !== null ? ' *' : ''}</text>` +
      `</g>`;
  }
  canvas.innerHTML = s;
  nodePressureInput.value =
    model.selection && model.selection.kind === 'node'
      ? (byName.get(model.nodes.find((n) => n.id === model.selection.id)?.name)?.pressure ?? '')
      : '';
}

function canvasPoint(evt) {
  const rect = canvas.getBoundingClientRect();
  return {
    x: ((evt.clientX - rect.left) / rect.width) * CW,
    y: ((evt.clientY - rect.top) / rect.height) * CH,
  };
}

canvas.addEventListener('pointerdown', (evt) => {
  const { x, y } = canvasPoint(evt);
  const hit = hitTestNode(model, x, y, 16);
  if (hit) {
    if (model.pendingLink && model.pendingLink !== hit.id) {
      model = addPipe(model, model.pendingLink, hit.id);
      model = { ...model, pendingLink: null };
    } else {
      dragId = hit.id;
      model = { ...model, selection: { kind: 'node', id: hit.id } };
    }
  } else {
    const added = addNode(model, x, y);
    model = { ...added.model, selection: { kind: 'node', id: added.id } };
    dragId = added.id;
  }
  renderCanvas();
  syncTextFromModel();
});

canvas.addEventListener('pointermove', (evt) => {
  if (!dragId) return;
  const { x, y } = canvasPoint(evt);
  model = moveNode(model, dragId, x, y, CW, CH);
  renderCanvas();
});

window.addEventListener('pointerup', () => {
  if (dragId) {
    dragId = null;
    syncTextFromModel();
  }
});

nodePressureInput.addEventListener('change', () => {
  if (!model.selection || model.selection.kind !== 'node') return;
  const v = nodePressureInput.value.trim();
  model = setNodePressure(model, model.selection.id, v === '' ? null : Number(v));
  renderCanvas();
  syncTextFromModel();
});

document.getElementById('addSourceBtn').addEventListener('click', () => {
  if (!model.selection || model.selection.kind !== 'node') return;
  const n = model.nodes.find((x) => x.id === model.selection.id);
  if (n) {
    model = addSource(model, n.name, 1e-6);
    syncTextFromModel();
    setStatus(`已给节点 ${n.name} 加源。`, 'ok');
  }
});

document.getElementById('deleteBtn').addEventListener('click', () => {
  if (!model.selection) return;
  model = model.selection.kind === 'node'
    ? removeNode(model, model.selection.id)
    : removePipe(model, model.selection.id);
  renderCanvas();
  syncTextFromModel();
});

window.addEventListener('keydown', (evt) => {
  if ((evt.key === 'Delete' || evt.key === 'Backspace') && document.activeElement.tagName !== 'TEXTAREA' && document.activeElement.tagName !== 'INPUT') {
    document.getElementById('deleteBtn').click();
  }
});

timeDtInput.addEventListener('change', () => {
  model = setTime(model, { dt: Number(timeDtInput.value) || 0.01 });
  syncTextFromModel();
});
timeEndInput.addEventListener('change', () => {
  model = setTime(model, { endTime: Number(timeEndInput.value) || 1 });
  syncTextFromModel();
});

// --- Run steady / transient ------------------------------------------------

async function run() {
  const text = caseText.value;
  const check = validateCaseText(text);
  if (!check.valid) {
    setStatus('算例校验失败:\n' + check.errors.join('\n'), 'error');
    return;
  }
  setStatus('正在运行求解器…');
  try {
    const summary = await invoke('run_case_text', {
      caseText: text,
      exePath: exePath.value.trim(),
    });
    renderBars('nodeChart', summary.nodes, 'pressure');
    renderLine('pipeChart', summary.pipes, 'q');
    renderConvergence(summary);
    setStatus(
      `运行完成：${summary.nodes.length} 节点 / ${summary.pipes.length} 管道，迭代 ${summary.iterations}。`,
      'ok',
    );
  } catch (e) {
    setStatus('运行失败: ' + ((e && e.message) || String(e)), 'error');
  }
}

function seriesValue(sample, key) {
  const sep = key.indexOf('|');
  const kind = key.slice(0, sep);
  const name = key.slice(sep + 1);
  const table = kind === 'pressure' ? sample.pressure
    : kind === 'alpha' ? sample.alpha
    : kind === 'hydrate' ? sample.hydrate
    : sample.flow;
  return table[name];
}

function renderTransient() {
  if (!transientSeries || !transientSeries.samples.length) {
    document.getElementById('transientChart').innerHTML = '<p class="empty">暂无瞬态数据</p>';
    return;
  }
  const key = seriesSelect.value;
  const result = toTimeSeriesPoints(
    transientSeries.samples,
    (s) => seriesValue(s, key),
    PLOT_W,
    PLOT_H,
  );
  if (!result.points) {
    document.getElementById('transientChart').innerHTML = '<p class="empty">该序列无数据</p>';
    return;
  }
  const pts = result.points.split(' ').map((p) => {
    const [x, y] = p.split(',').map(Number);
    return `${(x + PAD_L).toFixed(2)},${(y + PAD_T).toFixed(2)}`;
  }).join(' ');
  let s = `<svg class="chart" viewBox="0 0 ${W} ${H}" preserveAspectRatio="xMidYMid meet">${axis()}`;
  s += `<polyline points="${pts}" fill="none" stroke="#ffb454" stroke-width="2"/>`;
  // 时间滑条当前帧标记
  const idx = Number(timeSlider.value);
  const sample = transientSeries.samples[idx];
  if (sample) {
    const x = ((sample.t - result.tMin) / (result.tMax - result.tMin || 1)) * PLOT_W + PAD_L;
    const v = seriesValue(sample, key);
    const y = PLOT_H - ((v - result.min) / (result.max - result.min || 1)) * PLOT_H + PAD_T;
    s += `<line x1="${x}" y1="${PAD_T}" x2="${x}" y2="${PAD_T + PLOT_H}" stroke="#263349" stroke-dasharray="3 3"/>`;
    s += `<circle cx="${x}" cy="${y.toFixed(2)}" r="4" fill="#ff6b6b"/>`;
  }
  s += '</svg>';
  document.getElementById('transientChart').innerHTML = s;
  timeReadout.textContent =
    sample ? `t = ${Number(sample.t).toPrecision(6)} s · 值 = ${Number(seriesValue(sample, key)).toPrecision(6)}` : '—';
  document.getElementById('seriesLabel').textContent = key.replace('|', ' · ');
}

function populateSeries(series) {
  const opts = [];
  for (const n of series.node_names) {
    opts.push(['pressure', n], ['alpha', n], ['hydrate', n]);
  }
  for (const p of series.pipe_names) opts.push(['flow', p]);
  seriesSelect.innerHTML = opts
    .map(([kind, name]) => `<option value="${kind}|${escapeHtml(name)}">${kind} · ${escapeHtml(name)}</option>`)
    .join('');
}

async function runTransient() {
  if (!caseText.value.includes('time')) {
    setStatus('瞬态运行需要 case 含 time{} 块（dt / endTime）。', 'error');
    return;
  }
  setStatus('正在运行瞬态求解器…');
  try {
    const series = await invoke('run_transient_text', {
      caseText: caseText.value,
      exePath: transientExePath.value.trim(),
    });
    transientSeries = series;
    populateSeries(series);
    timeSlider.max = String(Math.max(0, series.samples.length - 1));
    timeSlider.value = '0';
    renderTransient();
    setStatus(
      `瞬态完成：${series.samples.length} 帧，t=0..${series.samples[series.samples.length - 1].t}，steps=${series.steps}。`,
      'ok',
    );
  } catch (e) {
    setStatus('瞬态运行失败: ' + ((e && e.message) || String(e)), 'error');
  }
}

seriesSelect.addEventListener('change', renderTransient);
timeSlider.addEventListener('input', renderTransient);

async function validateOnly() {
  const text = caseText.value;
  const local = validateCaseText(text);
  let msg = local.valid
    ? `本地校验通过：${local.nodes.length} 节点 / ${local.pipes.length} 管道 / ${local.sources.length} 源`
    : '本地校验失败:\n' + local.errors.join('\n');
  try {
    const data = await invoke('validate_case', { caseText: text });
    msg += `\nRust 校验通过：${data.nodes.length} 节点 / ${data.pipes.length} 管道`;
    setStatus(msg, 'ok');
  } catch (e) {
    setStatus(msg + '\nRust 校验失败: ' + ((e && e.message) || String(e)), 'error');
  }
}

caseText.value = DEFAULT_CASE;
exePath.value = DEFAULT_EXE;
transientExePath.value = DEFAULT_TRANSIENT_EXE;
document.getElementById('runBtn').addEventListener('click', run);
document.getElementById('runTransientBtn').addEventListener('click', runTransient);
document.getElementById('validateBtn').addEventListener('click', validateOnly);
document.getElementById('resetBtn').addEventListener('click', () => {
  caseText.value = DEFAULT_CASE;
  transientSeries = null;
  document.getElementById('transientChart').innerHTML = '<p class="empty">暂无瞬态数据</p>';
  setStatus('已恢复示例算例。');
});
renderCanvas();
