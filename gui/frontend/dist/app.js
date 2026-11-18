import { invoke } from './tauri-api.js';
import { validateCaseText, toBarHeights, toLinePoints } from './chart.js';

// The C++ engine built by the top-level CMake project; editable in the UI.
const DEFAULT_EXE = '/Users/haizai/project/ChuanFlow/build/src/1d/apps/chflow-1d-steady';

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
`;

const caseText = document.getElementById('caseText');
const exePath = document.getElementById('exePath');
const statusEl = document.getElementById('status');

const W = 480;
const H = 180;
const PAD_L = 34;
const PAD_R = 10;
const PAD_T = 10;
const PAD_B = 30;
const PLOT_W = W - PAD_L - PAD_R;
const PLOT_H = H - PAD_T - PAD_B;

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
document.getElementById('runBtn').addEventListener('click', run);
document.getElementById('validateBtn').addEventListener('click', validateOnly);
document.getElementById('resetBtn').addEventListener('click', () => {
  caseText.value = DEFAULT_CASE;
  setStatus('已恢复示例算例。');
});
