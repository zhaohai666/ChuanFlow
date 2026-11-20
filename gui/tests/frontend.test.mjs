import test from 'node:test';
import assert from 'node:assert/strict';
import {
  validateCaseText,
  toBarHeights,
  toLinePoints,
  toTimeSeriesPoints,
} from '../frontend/dist/chart.js';
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
  setFluid,
  setPipeFields,
  updateSourceRate,
  hitTestNode,
  networkToCaseText,
} from '../frontend/dist/editor.js';
import {
  UNIT_SYSTEMS,
  toDisplay,
  fromDisplay,
  unitLabel,
  fmtDisplay,
  parseInput,
} from '../frontend/dist/units.js';
import {
  buildProjectTree,
  parseCoupledOutput,
  parseVtu,
  cellCentroid,
  fieldRange,
  rainbowColor,
  classifyLogLine,
  parseSteadyOutput,
  parseTransientOutput,
  parseTwoFluidFrames,
  extractCaseDetails,
} from '../frontend/dist/workspace.js';
import {
  TEMPLATE_IDS,
  templateTarget,
  isPartialFidelity,
} from '../frontend/dist/templates.js';

const SAMPLE = `// 三节点串联管网
solver { type steadyNetwork; }
fluid
{
    rho 850;
    mu 0.1;
}

nodes
{
    A { pressure 3.0e4; }
    B { pressure 2.0e4; }
    E {};
}

pipes
{
    P1 { from A; to E; length 20; diameter 0.02; }
    P2 { from E; to B; length 20; diameter 0.02; }
}

sources
{
    E { rate 1e-6; }
}
`;

test('validateCaseText: 正常算例（分行块风格）', () => {
  const r = validateCaseText(SAMPLE);
  assert.equal(r.valid, true, r.errors.join('; '));
  assert.equal(r.nodes.length, 3);
  assert.equal(r.nodes[0].pressure, 30000);
  assert.equal(r.nodes[2].pressure, null);
  assert.equal(r.pipes.length, 2);
  assert.equal(r.sources.length, 1);
  assert.equal(r.sources[0].rate, 1e-6);
});

test('validateCaseText: 单行块风格', () => {
  const r = validateCaseText('nodes { A { pressure 1; } }\npipes { P { from A; to A; } }');
  assert.equal(r.valid, true, r.errors.join('; '));
  assert.equal(r.nodes.length, 1);
  assert.equal(r.pipes.length, 1);
});

test('validateCaseText: 空文本', () => {
  const r = validateCaseText('');
  assert.equal(r.valid, false);
  assert.ok(r.errors.some((e) => e.includes('为空')));
});

test('validateCaseText: 缺少 pipes', () => {
  const r = validateCaseText('nodes\n{\n A { pressure 1; }\n}\n');
  assert.equal(r.valid, false);
  assert.ok(r.errors.some((e) => e.includes('pipes')));
});

test('validateCaseText: 未知节点引用', () => {
  const bad = SAMPLE.replace('P2 { from E; to B;', 'P2 { from E; to MISSING;');
  const r = validateCaseText(bad);
  assert.equal(r.valid, false);
  assert.ok(r.errors.some((e) => e.includes('MISSING')));
});

test('toBarHeights: 按最大值归一化', () => {
  const h = toBarHeights([10, 20, 5], 100);
  assert.deepEqual(h.map((x) => Math.round(x)), [50, 100, 25]);
});

test('toBarHeights: 全零与空输入', () => {
  assert.deepEqual(toBarHeights([0, 0], 100), [0, 0]);
  assert.deepEqual(toBarHeights([], 100), []);
});

test('toLinePoints: 等距映射到 SVG 坐标', () => {
  const pts = toLinePoints([0, 10, 20], 100, 50).split(' ');
  assert.equal(pts.length, 3);
  assert.deepEqual(pts[0].split(',').map(Number), [0, 50]);
  assert.deepEqual(pts[2].split(',').map(Number), [100, 0]);
  const mid = pts[1].split(',').map(Number);
  assert.ok(Math.abs(mid[0] - 50) < 0.01);
  assert.ok(Math.abs(mid[1] - 25) < 0.01);
});

test('toLinePoints: 负值基线在 0', () => {
  const pts = toLinePoints([-5, 0, 5], 100, 100, 5).split(' ');
  assert.deepEqual(pts[0].split(',').map(Number), [0, 100]);
  assert.deepEqual(pts[2].split(',').map(Number), [100, 0]);
});

test('toLinePoints: 单点与空输入', () => {
  assert.deepEqual(toLinePoints([3], 100, 50).split(',').map(Number), [0, 0]);
  assert.equal(toLinePoints([], 100, 50), '');
});

// ---------------------------------------------------------------------------
// D12 增量2：拖拽建模纯函数 + 瞬态时间序列
// ---------------------------------------------------------------------------

function buildModel() {
  let m = emptyModel();
  const a = addNode(m, 100, 80, 'A'); m = a.model;
  const b = addNode(m, 400, 80, 'B'); m = b.model;
  m = setNodePressure(m, a.id, 5e5);
  m = addPipe(m, a.id, b.id);
  m = addSource(m, 'A', 0.002);
  return m;
}

test('editor: 拖拽布局 moveNode 限制在画布内', () => {
  let m = buildModel();
  const id = m.nodes[0].id;
  m = moveNode(m, id, 9999, -5, 640, 420);
  assert.deepEqual([m.nodes[0].x, m.nodes[0].y], [628, 12]);
});

test('editor: 命中测试返回最近节点，边界外为 null', () => {
  const m = buildModel();
  assert.equal(hitTestNode(m, 100, 80).name, 'A');
  assert.equal(hitTestNode(m, 400, 80).name, 'B');
  assert.equal(hitTestNode(m, 250, 250), null);
});

test('editor: removeNode 级联删除关联的管与源', () => {
  let m = buildModel();
  m = removeNode(m, m.nodes[0].id);
  assert.equal(m.nodes.length, 1);
  assert.equal(m.pipes.length, 0);
  assert.equal(m.sources.length, 0);
});

test('editor: 拒绝自环与重复管', () => {
  let m = emptyModel();
  const a = addNode(m, 10, 10, 'A'); m = a.model;
  const b = addNode(m, 50, 10, 'B'); m = b.model;
  assert.equal(addPipe(m, a.id, a.id).pipes.length, 0);   // 自环拒绝
  m = addPipe(m, a.id, b.id);
  assert.equal(addPipe(m, a.id, b.id).pipes.length, 1);   // 重复同向拒绝
  assert.equal(addPipe(m, b.id, a.id).pipes.length, 2);   // 反向是不同管，允许
});

test('editor: networkToCaseText 可被 validateCaseText 接受', () => {
  const text = networkToCaseText(buildModel());
  const r = validateCaseText(text);
  assert.equal(r.valid, true, r.errors.join('; '));
  assert.equal(r.nodes.length, 2);
  assert.equal(r.pipes.length, 1);
  assert.equal(r.sources[0].node, 'A');
  assert.ok(text.includes('node A;'), '必须输出显式 node= 键');
  assert.ok(text.includes('length 1000'), '必须输出管几何');
  assert.ok(text.includes('time {'), '必须输出 time 块');
});

test('editor: 源/管/时间编辑反映到序列化文本', () => {
  let m = buildModel();
  m = removeSource(m, m.sources[0].name);
  assert.equal(m.sources.length, 0);
  m = addSource(m, 'B', -1e-6);
  m = setTime(m, { dt: 0.02, endTime: 3 });
  const text = networkToCaseText(m);
  assert.ok(text.includes('node B; rate -0.000001'));
  assert.ok(text.includes('dt 0.02; endTime 3'));
});

test('toTimeSeriesPoints: 按真实 t 作横轴（非等距）', () => {
  const samples = [{ t: 0, p: 0 }, { t: 0.1, p: 10 }, { t: 1.0, p: 5 }];
  const r = toTimeSeriesPoints(samples, (s) => s.p, 100, 50);
  const pts = r.points.split(' ').map((p) => p.split(',').map(Number));
  assert.equal(pts.length, 3);
  assert.deepEqual(pts[0], [0, 50]);       // p=0 在底
  assert.deepEqual(pts[1], [10, 0]);       // t=0.1/1.0 -> x=10
  assert.ok(Math.abs(pts[2][0] - 100) < 0.01); // t=1.0 在最右
  assert.equal(r.min, 0);
  assert.equal(r.max, 10);
  assert.equal(r.tMin, 0);
  assert.equal(r.tMax, 1);
});

test('toTimeSeriesPoints: 空/单点/非有限值', () => {
  assert.equal(toTimeSeriesPoints([], (s) => s.v, 100, 50).points, '');
  const one = toTimeSeriesPoints([{ t: 2, v: 7 }], (s) => s.v, 100, 50);
  assert.deepEqual(one.points.split(',').map(Number), [0, 0]); // 单点在顶部
  const bad = toTimeSeriesPoints([{ t: 0, v: NaN }], (s) => s.v, 100, 50);
  assert.equal(bad.points, '');
});

test('editor: 序列化往返（networkToCaseText -> validate 结构一致）', () => {
  const m = buildModel();
  const r = validateCaseText(networkToCaseText(m));
  assert.deepEqual(r.nodes.map((n) => n.name), m.nodes.map((n) => n.name));
  assert.deepEqual(r.pipes.map((p) => [p.from, p.to]), m.pipes.map((p) => [p.from, p.to]));
  assert.equal(r.nodes[0].pressure, 5e5);
});


// ---------------------------------------------------------------------------
// M7：单位系统（显示层换算，字典恒 SI）
// ---------------------------------------------------------------------------

test('units: 长度/压强/流量 SI↔油田往返逐位回到 SI', () => {
  for (const family of ['length', 'pressure', 'rate', 'viscosity', 'density']) {
    const si = 123456.789;
    const oil = toDisplay(si, family, UNIT_SYSTEMS.OILFIELD);
    const back = fromDisplay(oil, family, UNIT_SYSTEMS.OILFIELD);
    assert.ok(Math.abs(back - si) <= 1e-9 * Math.abs(si), family + ' 往返失败');
  }
});

test('units: 温度 K↔°F 锚点换算', () => {
  assert.ok(Math.abs(toDisplay(273.15, 'temperature', UNIT_SYSTEMS.OILFIELD) - 32) < 1e-9);
  assert.ok(Math.abs(toDisplay(373.15, 'temperature', UNIT_SYSTEMS.OILFIELD) - 212) < 1e-9);
  assert.ok(Math.abs(fromDisplay(212, 'temperature', UNIT_SYSTEMS.OILFIELD) - 373.15) < 1e-9);
});

test('units: 已知换算系数', () => {
  assert.ok(Math.abs(toDisplay(1, 'length', UNIT_SYSTEMS.OILFIELD) - 3.28083989501312) < 1e-12);
  assert.ok(Math.abs(toDisplay(6894.757293168361, 'pressure', UNIT_SYSTEMS.OILFIELD) - 1) < 1e-12);
  assert.ok(Math.abs(toDisplay(1, 'viscosity', UNIT_SYSTEMS.OILFIELD) - 1000) < 1e-12);
  assert.equal(unitLabel('rate', UNIT_SYSTEMS.OILFIELD), 'bbl/d');
  assert.equal(unitLabel('pressure', UNIT_SYSTEMS.SI), 'Pa');
});

test('units: fmtDisplay 与 parseInput 互逆', () => {
  const text = fmtDisplay(2000, 'pressure', UNIT_SYSTEMS.OILFIELD, 6);
  assert.ok(text.endsWith(' psi'), text);
  const value = parseInput(text.split(' ')[0], 'pressure', UNIT_SYSTEMS.OILFIELD);
  assert.ok(Math.abs(value - 2000) < 0.01);
  assert.ok(Number.isNaN(parseInput('abc', 'pressure', UNIT_SYSTEMS.SI)));
});

// ---------------------------------------------------------------------------
// M7：编辑器新函数（管段几何 / 源流量 / 流体）
// ---------------------------------------------------------------------------

test('editor: setPipeFields / updateSourceRate / setFluid 反映到序列化', () => {
  let m = buildModel();
  const pid = m.pipes[0].id;
  m = setPipeFields(m, pid, { length: 2400, diameter: 0.073, elevation: -12.5 });
  m = updateSourceRate(m, m.sources[0].name, 3.5e-4);
  m = setFluid(m, { rho: 930, mu: 0.042 });
  const text = networkToCaseText(m);
  assert.ok(text.includes('length 2400'), text);
  assert.ok(text.includes('diameter 0.073'));
  assert.ok(text.includes('elevation -12.5'));
  assert.ok(text.includes('rate 0.00035'));
  assert.ok(text.includes('fluid { rho 930; mu 0.042; }'), '流体块必须来自模型');
  const r = validateCaseText(text);
  assert.equal(r.valid, true, r.errors.join('; '));
});

test('editor: setFluid 缺省模型也能回落到 network.cpp 默认值', () => {
  const m = { ...emptyModel(), fluid: undefined };
  const m2 = setFluid(m, { mu: 0.3 });
  assert.deepEqual(m2.fluid, { rho: 850, mu: 0.3 });
});

// ---------------------------------------------------------------------------
// M7：项目树
// ---------------------------------------------------------------------------

test('workspace: buildProjectTree 分层覆盖流体/1D/监测/结果', () => {
  const m = buildModel();
  const tree = buildProjectTree(m, { hasSteady: true });
  const top = tree.map((n) => n.id);
  assert.deepEqual(top, ['fluid', 'net1d', 'monitor', 'results']);
  const net = tree.find((n) => n.id === 'net1d');
  assert.equal(net.children.length, 4);
  const nodesGroup = net.children.find((c) => c.id === 'group.nodes');
  assert.equal(nodesGroup.children.length, m.nodes.length);
  const results = tree.find((n) => n.id === 'results');
  const steady = results.children.find((c) => c.id === 'res.steady');
  assert.equal(steady.ok, true);
  const transient = results.children.find((c) => c.id === 'res.transient');
  assert.equal(transient.ok, false);
});

test('workspace: buildProjectTree 有耦合边时出现耦合接口组', () => {
  const m = buildModel();
  const tree = buildProjectTree(m, { coupledEdge: 'P1' });
  assert.ok(tree.some((n) => n.id === 'coupled'));
  const none = buildProjectTree(m);
  assert.ok(!none.some((n) => n.id === 'coupled'));
});

// ---------------------------------------------------------------------------
// M7：耦合求解输出解析
// ---------------------------------------------------------------------------

const COUPLED_SAMPLE = `[chflow] node final p 1000669.324963046
[chflow] node final p 1000048.8923985178
[chflow] node final p 1000000
[chflow] pipe final A Q 0.00023999999999990885 v 0.030557749073632297 Re 24.446199258905839 dp 48.892398517811671
[chflow] pipe final X Q 0.00015777038680818138 v 0.50219873868084841 Re 80.351798188935746 dp 620.43256452819332
[chflow] coupled solved: edge X q 0.00015777038680818138 dp 620.4459384124541 R 3932503.2858195421 iterations 4 residual 2.1555277313915015e-05
[chflow] coupled seed 1: iterations 6 (converged)
[chflow] interface seed 1 iter 1 q 0.0002 dp 800 R 2e7
[chflow] interface seed 1 iter 2 q 0.00016 dp 640 R 4e6
`;

test('workspace: parseCoupledOutput 提取节点/管段/耦合边/交界面迭代', () => {
  const c = parseCoupledOutput(COUPLED_SAMPLE);
  assert.equal(c.nodes.length, 3);
  assert.equal(c.nodes[0].p, 1000669.324963046);
  assert.equal(c.pipes.length, 2);
  assert.equal(c.pipes[0].name, 'A');
  assert.equal(c.pipes[0].dp, 48.892398517811671);
  assert.equal(c.solved.edge, 'X');
  assert.equal(c.solved.iterations, 4);
  assert.equal(c.seeds[0].converged, true);
  assert.equal(c.interfaces.length, 2);
  assert.equal(c.interfaces[1].dp, 640);
  assert.equal(c.interfaces[1].resistance, 4e6);
});

// ---------------------------------------------------------------------------
// M7：稳态/瞬态 stdout 解析（JS 侧，作业式运行）
// ---------------------------------------------------------------------------

test('workspace: parseSteadyOutput 提取收敛与场', () => {
  const out = [
    '[chflow] node A p 30000',
    '[chflow] node B p 20000',
    '[chflow] pipe P1 Q 1e-6 v 0.003 Re 50 f 0.5 dp 5000',
    '[chflow] steady solved: iterations=7 residual=1e-9 (converged)',
  ].join('\n');
  const s = parseSteadyOutput(out);
  assert.equal(s.nodes.length, 2);
  assert.equal(s.pipes[0].dp, 5000);
  assert.equal(s.iterations, 7);
  assert.equal(s.converged, true);
  const bad = parseSteadyOutput('nothing useful');
  assert.equal(bad.nodes.length, 0);
  assert.equal(bad.converged, false);
});

test('workspace: parseTransientOutput 时间帧与序列', () => {
  const out = [
    '[chflow] transient: dt=0.01 endTime=0.3',
    '[chflow] t 0',
    '[chflow] node A p 200000 alpha 0.5 hydrate 0',
    '[chflow] pipe P Q 0',
    '[chflow] t 0.01',
    '[chflow] node A p 199000 alpha 0.49 hydrate 0.01',
    '[chflow] pipe P Q 1e-6',
    '[chflow] transient solved: steps=30',
  ].join('\n');
  const t = parseTransientOutput(out);
  assert.equal(t.dt, 0.01);
  assert.equal(t.steps, 30);
  assert.equal(t.samples.length, 2);
  assert.deepEqual(t.nodeNames, ['A']);
  assert.equal(t.samples[1].alpha.A, 0.49);
  assert.equal(t.samples[1].flow.P, 1e-6);
});

// ---------------------------------------------------------------------------
// M7：ASCII VTU 解析 + 色图
// ---------------------------------------------------------------------------

const VTU_SAMPLE = `<VTKFile type="UnstructuredGrid" version="0.1">
<Piece NumberOfPoints="4" NumberOfCells="1">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0 0 0  1 0 0  1 1 0  0 1 0
</DataArray>
</Points>
<Cells>
<DataArray type="Int64" Name="connectivity" format="ascii">
0 1 2 3
</DataArray>
<DataArray type="Int64" Name="offsets" format="ascii">
4
</DataArray>
</Cells>
<CellData>
<DataArray type="Float64" Name="U" NumberOfComponents="3" format="ascii">
0.5 0 0
</DataArray>
<DataArray type="Float64" Name="p" format="ascii">
101325
</DataArray>
</CellData>
</Piece>
</VTKFile>`;

test('workspace: parseVtu 读取 ASCII 场并算质心', () => {
  const v = parseVtu(VTU_SAMPLE);
  assert.ok(v, '必须解析成功');
  assert.equal(v.points.length, 4);
  assert.equal(v.cells.length, 1);
  assert.deepEqual(v.cells[0], [0, 1, 2, 3]);
  assert.equal(v.cellP[0], 101325);
  const c = cellCentroid(v.points, v.cells[0]);
  assert.equal(c.x, 0.5);
  assert.equal(c.y, 0.5);
  assert.equal(parseVtu('not a vtu'), null);
});

test('workspace: fieldRange / rainbowColor 端点与钳制', () => {
  assert.deepEqual(fieldRange([3, 1, 2]), [1, 3]);
  assert.deepEqual(fieldRange([]), [0, 1]);
  assert.equal(rainbowColor(3, 0, 3), rainbowColor(9, 0, 3), '超界值钳制到同一端色');
  assert.ok(rainbowColor(0, 0, 3).startsWith('rgb('));
});

// ---------------------------------------------------------------------------
// M7：日志分级
// ---------------------------------------------------------------------------

test('workspace: classifyLogLine 红黄绿分级', () => {
  assert.equal(classifyLogLine('cfFatal: bad input'), 'error');
  assert.equal(classifyLogLine('WARNING: something odd'), 'warn');
  assert.equal(classifyLogLine('steady solved: iterations=7 (converged)'), 'ok');
  assert.equal(classifyLogLine('[chflow] t 0.5 vm 1.2'), 'info');
});

// ---------------------------------------------------------------------------
// M7：字典几何提取（dict -> 模型忠实重建）
// ---------------------------------------------------------------------------

test('workspace: extractCaseDetails 提取流体/时间/管段几何', () => {
  const d = extractCaseDetails(SAMPLE);
  assert.deepEqual(d.fluid, { rho: 850, mu: 0.1 });
  assert.equal(d.time, null);
  assert.equal(d.pipes.P1.length, 20);
  assert.equal(d.pipes.P1.diameter, 0.02);
  assert.equal(d.pipes.P2.length, 20);
});

test('workspace: extractCaseDetails 时间块与单行管段', () => {
  const text = 'fluid { rho 930; mu 0.042; }\npipes { P { from A; to B; length 42; diameter 0.073; } }\ntime { dt 0.02; endTime 3; reportEvery 5; }';
  const d = extractCaseDetails(text);
  assert.deepEqual(d.fluid, { rho: 930, mu: 0.042 });
  assert.deepEqual(d.time, { dt: 0.02, endTime: 3, reportEvery: 5 });
  assert.equal(d.pipes.P.length, 42);
});


// ---------------------------------------------------------------------------
// M8：设备组件（阀/泵）字段进出字典
// ---------------------------------------------------------------------------

test('editor: 阀/泵键写入序列化并忠实提取', () => {
  let m = buildModel();
  const pid = m.pipes[0].id;
  m = setPipeFields(m, pid, {
    pumpHead: 25, pumpRate: 3.0e-4, pumpShutoffHead: 40,
  });
  const pid2 = m.pipes[1] ? m.pipes[1].id : null;
  if (!pid2) {
    const r = validateCaseText(networkToCaseText(m));
    assert.equal(r.valid, true);
  }
  const text = networkToCaseText(m);
  assert.ok(text.includes('pumpHead 25'), text);
  assert.ok(text.includes('pumpRate 0.0003'));
  assert.ok(text.includes('pumpShutoffHead 40'));
  assert.ok(!text.includes('valve'), '未设置时不得输出 valve 键');
  const d = extractCaseDetails(text);
  assert.equal(d.pipes.P1.pumpHead, 25);
  assert.equal(d.pipes.P1.valve, 0);
});

test('editor: 阀键序列化与提取', () => {
  let m = buildModel();
  m = setPipeFields(m, m.pipes[0].id, { valve: 0.012 });
  const text = networkToCaseText(m);
  assert.ok(text.includes('valve 0.012'));
  assert.ok(!text.includes('pumpHead'), '未设置时不得输出 pumpHead 键');
  const d = extractCaseDetails(text);
  assert.equal(d.pipes.P1.valve, 0.012);
  assert.equal(d.pipes.P1.pumpHead, 0);
});

// ---------------------------------------------------------------------------
// M8：twoFluid 帧解析（段塞动画数据源）
// ---------------------------------------------------------------------------

const SLUG_SAMPLE = `[chflow] twoFluid: 50 sections, dt=0.01 endTime=120
[chflow] t 0 vm 0 dp -23779.94
[chflow] cell 0 alpha 0.55 water 0 vsl 0 vsg 0 s 0
[chflow] cell 1 alpha 0.55 water 0 vsl 0 vsg 0 s 0
[chflow] t 1 vm 0 dp -23779.94
[chflow] cell 0 alpha 0.61 water 0 vsl 0 vsg 0 s 0
[chflow] cell 1 alpha 0.55 water 0 vsl 0 vsg 0 s 1
[chflow] t 2 vm 0.3 dp -23800
[chflow] cell 0 alpha 0.7 water 0 vsl 0.2 vsg 0.1 T 300 p 1e6 s 0
[chflow] cell 1 alpha 0.5 water 0.1 vsl 0.15 vsg 0.15 T 301 p 1.1e6 s 2
`;

test('workspace: parseTwoFluidFrames 帧与逐格键值', () => {
  const fr = parseTwoFluidFrames(SLUG_SAMPLE);
  assert.equal(fr.length, 3);
  assert.equal(fr[0].t, 0);
  assert.equal(fr[0].vm, 0);
  assert.equal(fr[1].cells.length, 2);
  assert.equal(fr[1].cells[1].alpha, 0.55);
  assert.equal(fr[1].cells[1].s, 1);
  // tail keys vary by case features (T/p appear only when thermal/EOS on)
  assert.equal(fr[2].cells[0].T, 300);
  assert.equal(fr[2].cells[1].p, 1.1e6);
  assert.equal(fr[2].cells[1].s, 2);
  assert.equal(parseTwoFluidFrames('no frames here').length, 0);
});

// ---------------------------------------------------------------------------
// M8：模板库元数据
// ---------------------------------------------------------------------------

test('templates: 模板目标页与保真度标注', () => {
  assert.equal(TEMPLATE_IDS.length, 7);
  assert.deepEqual(templateTarget('twoFluidTerrainSlug'), { tab: 'slug', kind: 'slug' });
  assert.deepEqual(templateTarget('coupled1d3d'), { tab: 'coupled', kind: 'coupled' });
  assert.deepEqual(templateTarget('steadyEquipmentPump'), { tab: 'canvas', kind: 'dict' });
  assert.equal(isPartialFidelity('steadyNetworkMultiphase'), true);
  assert.equal(isPartialFidelity('steadyEquipmentValve'), false);
});
