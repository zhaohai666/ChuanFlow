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
  clipKeepsByPlane,
  adaptiveVectorIndices,
  applyWhatIfOverrides,
  TRANSIENT_WHAT_IF_PARAMS,
  iprLinearVariant,
  applyScenarioOverrides,
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
  traceStreamlines,
  cellVolumes,
  evaluateReports,
  reportsToCsv,
  classifyLogLine,
  parseSteadyOutput,
  parseTransientOutput,
  parseTwoFluidFrames,
  extractCaseDetails,
  parseNodalIpr,
  iprCurvePoints,
  vlpSampleJs,
  slugSeriesKeys,
  slugBands,
  buildMonitorTree,
  monitorSeriesValues,
  monitorEnvBand,
  compareSteadyRuns,
  scenarioRowsTransient,
  scenarioRowsSteady,
  pvtFlashCaseText,
  parsePvtOutput,
  pvtCurvesByT,
  parsePvtCaseDict,
  vorticityAndQ,
  isoSurfacePatches,
  faceAdjacency,
  fieldsToCsv,
  planeProfileCsv,
} from '../frontend/dist/workspace.js';
import {
  TEMPLATE_IDS,
  templateTarget,
  isPartialFidelity,
  isSteadyComparable,
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
  assert.equal(TEMPLATE_IDS.length, 9);
  assert.deepEqual(templateTarget('twoFluidTerrainSlug'), { tab: 'slug', kind: 'slug' });
  assert.deepEqual(templateTarget('coupled1d3d'), { tab: 'coupled', kind: 'coupled' });
  assert.deepEqual(templateTarget('nodalVogel'), { tab: 'nodal', kind: 'nodal' });
  assert.deepEqual(templateTarget('pvtFlash'), { tab: 'pvt', kind: 'pvt' });
  assert.deepEqual(templateTarget('steadyEquipmentPump'), { tab: 'canvas', kind: 'dict' });
  assert.equal(isPartialFidelity('steadyNetworkMultiphase'), true);
  assert.equal(isPartialFidelity('steadyEquipmentValve'), false);
});

test('clipKeepsByPlane: 轴对齐等价旧语义 + 任意法向 + 退化', () => {
  const points = [
    { x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 1, y: 1, z: 0 },
    { x: 0, y: 1, z: 0 }, { x: 0, y: 0, z: 2 }, { x: 1, y: 0, z: 2 },
    { x: 1, y: 1, z: 2 }, { x: 0, y: 1, z: 2 },
  ];
  const cells = [
    [0, 1, 2, 3], // z 层 0（质心 z=0）
    [4, 5, 6, 7], // z 层 2（质心 z=2）
  ];
  // frac=0.5 在两格质心中点处剖切：x 法向两格质心 x 相同 → 保留语义不炸
  assert.deepEqual(clipKeepsByPlane(cells, points, [1, 0, 0], 0.5), [true, true]);
  // z 法向、frac=0.5：只保留低侧（z=0 层）
  assert.deepEqual(clipKeepsByPlane(cells, points, [0, 0, 1], 0.5), [true, false]);
  // 非单位长度法向：归一化后同一结果
  assert.deepEqual(clipKeepsByPlane(cells, points, [0, 0, 7], 0.5), [true, false]);
  // 反向法向翻转保留侧
  assert.deepEqual(clipKeepsByPlane(cells, points, [0, 0, -1], 0.5), [false, true]);
  // 退化（零向量）全保留
  assert.deepEqual(clipKeepsByPlane(cells, points, [0, 0, 0], 0.5), [true, true]);
});

test('adaptiveVectorIndices: 上限、最强必保、确定性与覆盖', () => {
  const mags = Array.from({ length: 1000 }, (_, i) => Math.abs(Math.sin(i)));
  mags[777] = 1e9; // 最强格
  const idx = adaptiveVectorIndices(mags, 400);
  assert.equal(idx.length <= 400, true);
  assert.equal(idx.includes(777), true, '最强格必须保留');
  assert.deepEqual(idx, [...idx].sort((a, b) => a - b), '输出升序');
  assert.deepEqual(adaptiveVectorIndices(mags, 400), idx, '确定性');
  // 小于上限时全保留
  assert.deepEqual(adaptiveVectorIndices([1, 2, 3], 10), [0, 1, 2]);
  // 非法上限
  assert.deepEqual(adaptiveVectorIndices([1, 2], 0), []);
  // 覆盖：非热点也有被抽中的（rest 半额均匀采样）
  const hot = adaptiveVectorIndices(Array.from({ length: 1000 }, () => 0.5), 400);
  assert.equal(hot.length <= 400, true);
  assert.ok(hot.length > 200, '覆盖半额应填满采样位');
});

test('设备族：序列化往返 + 平管不泄漏 + details 解析', () => {
  let m = emptyModel();
  const { model: m1, id: a } = addNode(m, 0, 0, 'A');
  const { model: m2, id: b } = addNode(m1, 100, 100, 'B');
  const m3 = addPipe(m2, a, b); // addPipe returns the model directly
  m = setPipeFields(m3, m3.pipes[0].id, {
    liftDomePressure: 1.6e7, liftAreaRatio: 0.31,
    rodPlungerDiameter: 0.038, rodStrokesPerMinute: 6,
    pumpGasFraction: 0.25, pumpIntakePressure: 3.5e6,
  });
  const text = networkToCaseText(m);
  for (const key of ['liftDomePressure', 'liftAreaRatio', 'rodPlungerDiameter',
    'rodStrokesPerMinute', 'pumpGasFraction', 'pumpIntakePressure']) {
    assert.ok(text.includes(key), `missing ${key}`);
  }
  // 平管字节恒等：无设备键时不泄漏（M5 语义保持）
  const plain = networkToCaseText(emptyModel());
  assert.equal(plain.includes('liftDomePressure'), false);
  // details 解析读回全部设备键
  const details = extractCaseDetails(text);
  const p = details.pipes[Object.keys(details.pipes)[0]];
  assert.equal(p.liftDomePressure, 1.6e7);
  assert.equal(p.rodStrokesPerMinute, 6);
  assert.equal(p.pumpIntakePressure, 3.5e6);
});

// ---------------------------------------------------------------------------
// U10 / §5.137: transient what-if overrides on the raw dictionary
// ---------------------------------------------------------------------------

const WHAT_IF_CASE = `// real transientValve shape (volume/bulk/schedules preserved)
fluid { rho 850; mu 0.1; }
nodes
{
    TANK { volume 5; bulkModulus 4.0e7; pressure 5.0e5; liquidFraction 1; }
    OUT  { fixedPressure 0; }
}
pipes
{
    LINE
    {
        from TANK; to OUT; length 100; diameter 0.05;
        valve 20;
        valveClose { start 2; end 6; closedK 1.0e12; }
        valveOpen  { start 10; end 14; }
    }
}
sources
{
    inj { node TANK; rate 0;
      schedule { openRate 2.0e-6;
                 close { start 1.2; end 1.4; }
                 open  { start 2.0; end 2.2; } } }
}
time { dt 0.02; endTime 20; reportEvery 1; }
`;

test('what-if：关阀/开阀/关闭 K 精确替换嵌套块且不碰其他字节', () => {
  const a = applyWhatIfOverrides(WHAT_IF_CASE, 'valveCloseStart', 'LINE', 3.5);
  assert.ok(a.includes('valveClose { start 3.5; end 6; closedK 1.0e12; }'));
  assert.ok(a.includes('valveOpen  { start 10; end 14; }'));
  const b = applyWhatIfOverrides(WHAT_IF_CASE, 'valveOpenEnd', 'LINE', 15);
  assert.ok(b.includes('valveOpen  { start 10; end 15; }'));
  assert.ok(b.includes('valveClose { start 2; end 6;'));
  const c = applyWhatIfOverrides(WHAT_IF_CASE, 'closedK', 'LINE', 5e11);
  assert.ok(c.includes('closedK 500000000000;'));
  // 只有目标数值变化：把改后的值换回原值必须逐字节回到输入
  assert.equal(applyWhatIfOverrides(a, 'valveCloseStart', 'LINE', 2), WHAT_IF_CASE);
});

test('what-if：源日程与全局 time 参数', () => {
  assert.ok(applyWhatIfOverrides(WHAT_IF_CASE, 'scheduleOpenRate', 'inj', 3e-6)
    .includes('openRate 0.000003;'));
  const closed = applyWhatIfOverrides(WHAT_IF_CASE, 'scheduleCloseStart', 'inj', 0.8);
  assert.ok(closed.includes('close { start 0.8; end 1.4; }'));
  assert.ok(closed.includes('open  { start 2.0; end 2.2; }'));
  const opened = applyWhatIfOverrides(WHAT_IF_CASE, 'scheduleOpenStart', 'inj', 2.5);
  assert.ok(opened.includes('open  { start 2.5; end 2.2; }')); // 引擎负责 start>end 拒收
  assert.ok(applyWhatIfOverrides(WHAT_IF_CASE, 'endTime', null, 30)
    .includes('time { dt 0.02; endTime 30; reportEvery 1; }'));
  assert.ok(applyWhatIfOverrides(WHAT_IF_CASE, 'dt', null, 0.01)
    .includes('time { dt 0.01; endTime 20;'));
  assert.ok(applyWhatIfOverrides(WHAT_IF_CASE, 'reportEvery', null, 5)
    .includes('reportEvery 5;'));
});

test('what-if：参数表完备 + 非法值/缺目标/缺键大声拒收', () => {
  assert.equal(TRANSIENT_WHAT_IF_PARAMS.length, 13);
  assert.throws(() => applyWhatIfOverrides(WHAT_IF_CASE, 'noSuchParam', 'LINE', 1), /unsupported/);
  assert.throws(() => applyWhatIfOverrides(WHAT_IF_CASE, 'valveCloseStart', 'LINE', NaN), /not finite/);
  assert.throws(() => applyWhatIfOverrides(WHAT_IF_CASE, 'valveCloseStart', null, 3), /needs a target/);
  assert.throws(() => applyWhatIfOverrides(WHAT_IF_CASE, 'valveCloseStart', 'NOPE', 3), /no pipes entry/);
  assert.throws(() => applyWhatIfOverrides(WHAT_IF_CASE, 'scheduleOpenRate', 'LINE', 3), /no sources entry/);
  assert.throws(() => applyWhatIfOverrides('time { dt 1; }', 'endTime', null, 3), /no endTime key/);
  assert.throws(() => applyWhatIfOverrides('nodes {}', 'endTime', null, 3), /no time\{\} block/);
  // 注释里的同名块不得被误认（braceEnd 跳过 // 行）
  assert.throws(() => applyWhatIfOverrides('// valveClose { start 1; }\nnodes {}',
    'valveCloseStart', 'LINE', 3), /no pipes entry/);
});

// ---------------------------------------------------------------------------
// U7 / §5.138: streamline integration through the cell velocity field
// ---------------------------------------------------------------------------

test('流线：均匀场严格直线、旋转场切向、零速/出界终止', () => {
  const grid = (nx, ny, nz) => {
    const points = [], cells = [];
    // hex vertices, one shared point list per lattice node; cell = 8 corners
    for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++)
      points.push({ x: i, y: j, z: k });
    const id = (i, j, k) => i + (nx + 1) * (j + (ny + 1) * k);
    for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++)
      cells.push([id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
        id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)]);
    return { points, cells };
  };
  const { points, cells } = grid(6, 4, 2);
  const uniform = cells.map(() => [2, 0, 0]);
  const lines = traceStreamlines(points, cells, uniform, [[0.5, 1.5, 0.5]], { maxSteps: 20 });
  assert.equal(lines.length, 1);
  const line = lines[0];
  assert.ok(line.length > 10);
  assert.ok(line.every((p) => Math.abs(p[1] - 1.5) < 1e-12 && Math.abs(p[2] - 0.5) < 1e-12));
  assert.ok(line[1][0] > line[0][0]);
  assert.ok(line[line.length - 1][0] <= 6.2, 'stops at the domain boundary');
  assert.equal(traceStreamlines(points, cells, cells.map(() => [0, 0, 0]),
    [[0.5, 0.5, 0.5]], { maxSteps: 20 }).length, 0);

  const rotCells = grid(9, 9, 1);
  const rotU = rotCells.cells.map((vids) => {
    const c = cellCentroid(rotCells.points, vids);
    return [-c.y + 4, c.x - 4, 0];
  });
  const circle = traceStreamlines(rotCells.points, rotCells.cells, rotU,
    [[5, 4, 0.5]], { maxSteps: 60, stepScale: 0.05 });
  assert.equal(circle.length, 1);
  assert.ok(circle[0].length > 20);
  assert.ok(circle[0][1][1] > 4, 'counter-clockwise tangent at seed');
  const radii = circle[0].map((p) => Math.hypot(p[0] - 4, p[1] - 4));
  assert.ok(Math.max(...radii) < 1.5 && Math.min(...radii) > 0.6,
    `Euler drift stays bounded: ${Math.min(...radii)}..${Math.max(...radii)}`);
});

test('流线：大网格 + 多种子的性能预算', () => {
  const nx = 20, ny = 20, nz = 8;
  const points = [], cells = [];
  for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++)
    points.push({ x: i, y: j, z: k });
  const id = (i, j, k) => i + (nx + 1) * (j + (ny + 1) * k);
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++)
    cells.push([id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k),
      id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)]);
  const u = cells.map(() => [1, 0.2, 0]);
  const seeds = Array.from({ length: 24 }, (_, i) => [i % 20 + 0.5, Math.floor(i / 20) * 9 + 0.5, 4.5]);
  const t0 = Date.now();
  const lines = traceStreamlines(points, cells, u, seeds, { maxSteps: 60 });
  assert.ok(lines.length >= 20);
  assert.ok(Date.now() - t0 < 3000, 'streamline tracing stays under the render budget');
});

// ---------------------------------------------------------------------------
// U8 / §5.139: report definitions from real VTU face geometry
// ---------------------------------------------------------------------------

test('report definitions：六面体闭式体积/体积分/面通量/面平均 + CSV', () => {
  const points = [
    { x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 1, y: 1, z: 0 }, { x: 0, y: 1, z: 0 },
    { x: 0, y: 0, z: 0.5 }, { x: 1, y: 0, z: 0.5 }, { x: 1, y: 1, z: 0.5 }, { x: 0, y: 1, z: 0.5 },
  ];
  const vtu = {
    points,
    cells: [[0, 1, 2, 3, 4, 5, 6, 7]],
    cellFaces: [[[0, 3, 7, 4], [1, 5, 6, 2], [0, 4, 5, 1], [3, 2, 6, 7], [0, 1, 2, 3], [4, 7, 6, 5]]],
    cellU: [[2, 0, 0]],
    cellP: [10],
  };
  const vols = cellVolumes(points, vtu.cells, vtu.cellFaces);
  assert.ok(Math.abs(vols[0] - 0.5) < 1e-14, `volume ${vols[0]}`);
  const defs = [
    { name: 'pInt', field: 'p', op: 'volumeIntegral' },
    { name: 'pAvg', field: 'p', op: 'volumeAverage' },
    { name: 'pMax', field: 'p', op: 'max' },
    { name: 'pMin', field: 'p', op: 'min' },
    { name: 'ux', field: 'Ux', op: 'volumeAverage' },
    { name: 'flux', field: '|U|', op: 'flux', normal: [1, 0, 0], frac: 0 },
    { name: 'surf', field: 'p', op: 'surfaceAverage', normal: [1, 0, 0], frac: 0 },
  ];
  const r = evaluateReports(vtu, defs).reduce((acc, x) => ({ ...acc, [x.name]: x.value }), {});
  assert.ok(Math.abs(r.pInt - 5) < 1e-14);
  assert.ok(Math.abs(r.pAvg - 10) < 1e-14);
  assert.equal(r.pMax, 10);
  assert.equal(r.pMin, 10);
  assert.equal(r.ux, 2);
  assert.ok(Math.abs(r.flux - 2 * 0.5) < 1e-14, `flux ${r.flux}`);
  assert.ok(Math.abs(r.surf - 10) < 1e-14, `surface average ${r.surf}`);
  const csv = reportsToCsv([{ pInt: 5, pAvg: 10 }, { pInt: 6, pAvg: 12 }], defs.slice(0, 2));
  assert.equal(csv, 'step,pInt,pAvg\n1,5,10\n2,6,12\n');
  assert.ok(Number.isNaN(evaluateReports({ ...vtu, cellFaces: null }, defs)[0].value));
});

test('report definitions：两单元体积加权 ≠ 算术平均', () => {
  const points = [
    { x: 0, y: 0, z: 0 }, { x: 1, y: 0, z: 0 }, { x: 1, y: 1, z: 0 }, { x: 0, y: 1, z: 0 },
    { x: 2, y: 0, z: 0 }, { x: 2, y: 1, z: 0 },
    { x: 0, y: 0, z: 1 }, { x: 1, y: 0, z: 1 }, { x: 1, y: 1, z: 1 }, { x: 0, y: 1, z: 1 },
    { x: 2, y: 0, z: 1 }, { x: 2, y: 1, z: 1 },
  ];
  const facesA = [[0, 3, 9, 6], [1, 2, 8, 7], [0, 6, 7, 1], [3, 2, 8, 9],
    [0, 1, 2, 3], [6, 9, 8, 7]];
  const facesB = [[1, 7, 8, 2], [4, 10, 11, 5], [1, 4, 10, 7], [2, 8, 11, 5],
    [1, 2, 5, 4], [7, 10, 11, 8]];
  const vtu = {
    points,
    cells: [[0, 1, 2, 3, 6, 7, 8, 9], [1, 4, 5, 2, 7, 10, 11, 8]],
    cellFaces: [facesA, facesB],
    cellU: [[0, 0, 0], [0, 0, 0]],
    cellP: [100, 1000],
  };
  const vols = cellVolumes(points, vtu.cells, vtu.cellFaces);
  assert.ok(Math.abs(vols[0] - 1) < 1e-14 && Math.abs(vols[1] - 1) < 1e-14);
  const [vInt, vAvg] = evaluateReports(vtu, [
    { name: 'i', field: 'p', op: 'volumeIntegral' },
    { name: 'a', field: 'p', op: 'volumeAverage' },
  ]);
  assert.ok(Math.abs(vInt.value - 1100) < 1e-12);
  assert.ok(Math.abs(vAvg.value - 550) < 1e-12);
});

// ---------------------------------------------------------------------------
// U2：节点分析（IPR 声明解析 + 闭式曲线 + VLP 采样字典）
// ---------------------------------------------------------------------------

const NODAL_DICT = `fluid { blackoil }
ipr
{
    form              vogel;
    reservoirPressure 2.2e7;  // Pa
    maxRate           0.006;  // sm3/s
}
well { length 2400; diameter 0.062; }
wellheadPressure 2.5e6;
`;

test('parseNodalIpr: vogel 声明与等效 J', () => {
  const spec = parseNodalIpr(NODAL_DICT);
  assert.equal(spec.form, 'vogel');
  assert.equal(spec.pRes, 2.2e7);
  assert.equal(spec.qmax, 0.006);
  assert.equal(spec.qAof, 0.006);
  assert.ok(Math.abs(spec.jRef - 0.006 / 2.2e7) < 1e-16);
  assert.throws(() => parseNodalIpr('well { length 1; }'), /ipr/);
});

test('parseNodalIpr: linear 与 jones', () => {
  const lin = parseNodalIpr('ipr { form linear; reservoirPressure 2e7; productivityIndex 1e-10; }');
  assert.equal(lin.qAof, 2e7 * 1e-10);
  const jon = parseNodalIpr('ipr { form jones; reservoirPressure 2e7; deliverability 1e-12; exponent 0.5; }');
  assert.ok(Math.abs(jon.qAof - 1e-12 * 2e7) < 1e-12);
});

test('iprCurvePoints: 三种闭式形态的端点与单调性', () => {
  const lin = iprCurvePoints({ form: 'linear', pRes: 2e7, J: 1e-10, qAof: 2e-3 }, 10);
  assert.equal(lin.length, 11);
  assert.equal(lin[0].p, 2e7);
  assert.ok(Math.abs(lin.at(-1).p) < 1e-6);
  const vog = iprCurvePoints({ form: 'vogel', pRes: 2.2e7, qmax: 0.006, qAof: 0.006 }, 10);
  assert.ok(Math.abs(vog[0].q) < 1e-12); // p=pRes -> q=0
  assert.ok(Math.abs(vog.at(-1).q - 0.006) < 1e-12); // p=0 -> qmax
  for (let i = 1; i < vog.length; i += 1) {
    assert.ok(vog[i].q > vog[i - 1].q && vog[i].p < vog[i - 1].p, 'q 升 p 降');
  }
  const jon = iprCurvePoints({ form: 'jones', pRes: 2e7, C: 1e-10, n: 0.5, qAof: 2e-3 }, 10);
  assert.ok(Math.abs(jon[0].q) < 1e-12 && Math.abs(jon.at(-1).q - 2e-3) < 1e-9);
});

test('vlpSampleJs + iprLinearVariant: 采样网格与字典替换', () => {
  const spec = parseNodalIpr(NODAL_DICT);
  const js = vlpSampleJs(spec, 7);
  assert.equal(js.length, 7);
  assert.ok(js[0] < js.at(-1) && js.every(Number.isFinite));
  const variant = iprLinearVariant(NODAL_DICT, spec.pRes, 1.23e-10);
  assert.match(variant, /form\s+linear;/);
  assert.match(variant, /productivityIndex 1\.23e-10;/);
  assert.doesNotMatch(variant, /vogel|maxRate/);
  assert.ok(variant.includes('wellheadPressure 2.5e6;'), '其余字节不动');
  assert.ok(variant.includes('well { length 2400; diameter 0.062; }'));
  const round = parseNodalIpr(variant);
  assert.equal(round.form, 'linear');
  assert.equal(round.J, 1.23e-10);
  assert.throws(() => iprLinearVariant('no ipr here', 1, 1), /ipr/);
});

// ---------------------------------------------------------------------------
// U3：段塞多量叠加（可用键 + 流型/流动保障色带）
// ---------------------------------------------------------------------------

test('slugSeriesKeys: 收集内核 cell 行的全部数值键', () => {
  const out = `t 1.0 vm 0.5 dp 100
cell 0 alpha 0.6 water 0.2 p 5e5 T 320 s 0.1
cell 1 alpha 0.4 water 0.1 p 4.9e5 T 321 s 0
t 2.0 vm 0.6 dp 120
cell 0 alpha 0.7 water 0.3 p 5.1e5 muEff 0.02 wax 0.01
`;
  const frames = parseTwoFluidFrames(out);
  assert.equal(frames.length, 2);
  const keys = slugSeriesKeys(frames);
  for (const k of ['alpha', 'water', 'p', 'T', 's', 'muEff', 'wax']) {
    assert.ok(keys.includes(k), `${k} 应可用`);
  }
  assert.ok(!keys.includes('index'));
});

test('slugBands: 连续段合并 + 仅输出帧内存在的键', () => {
  const frame = {
    cells: [
      { index: 0, s: 0.2, h: 0 },
      { index: 1, s: 0.3, h: 0 },
      { index: 2, s: 0, h: 0 },
      { index: 3, s: 0.1, h: 0.5 },
      { index: 4, s: 0, h: 0.4 },
    ],
  };
  const bands = slugBands(frame);
  assert.equal(bands.length, 2); // s 与 h（wax/sand/inh 帧内不存在）
  const s = bands.find((b) => b.key === 's');
  assert.deepEqual(s.runs, [{ from: 0, to: 1 }, { from: 3, to: 3 }]);
  const h = bands.find((b) => b.key === 'h');
  assert.deepEqual(h.runs, [{ from: 3, to: 4 }]);
  assert.equal(s.label, '段塞');
});

// ---------------------------------------------------------------------------
// U5：变量选择树（按声明量分层的勾选序列）
// ---------------------------------------------------------------------------

const MON_T = {
  nodeNames: ['A', 'B'],
  pipeNames: ['P1'],
  samples: [
    { t: 0, pressure: { A: 1e6, B: 2e6 }, alpha: { A: 0.4, B: 0.6 }, hydrate: {}, flow: { P1: 0.01 } },
    { t: 1, pressure: { A: 1.2e6, B: 1.8e6 }, alpha: { A: 0.5, B: 0.7 }, hydrate: {}, flow: { P1: 0.02 } },
  ],
};

test('buildMonitorTree: 只生成声明量，分节点/管段/瞬态帧三层', () => {
  const groups = buildMonitorTree(MON_T);
  assert.deepEqual(groups.map((g) => g.id), ['nodes', 'pipes', 'frames']);
  const nodes = groups[0];
  assert.equal(nodes.items.length, 2);
  assert.deepEqual(nodes.items[0].series, ['p:A', 'alpha:A']); // hydrate 未声明 -> 不出现
  assert.deepEqual(groups[1].items[0].series, ['Q:P1']);
  const frames = groups[2];
  const fSeries = frames.items.flatMap((it) => it.series);
  assert.ok(fSeries.includes('env:p') && fSeries.includes('avg:alpha') && fSeries.includes('marker'));
});

test('monitorSeriesValues / monitorEnvBand: 逐帧取值与包络', () => {
  assert.deepEqual(monitorSeriesValues(MON_T.samples, 'p:A'), [1e6, 1.2e6]);
  assert.deepEqual(monitorSeriesValues(MON_T.samples, 'Q:P1'), [0.01, 0.02]);
  assert.deepEqual(monitorSeriesValues(MON_T.samples, 'hydrate:A'), [undefined, undefined]);
  assert.deepEqual(monitorSeriesValues(MON_T.samples, 'bogus:X'), null);
  assert.deepEqual(monitorSeriesValues(MON_T.samples, 'avg:alpha'), [0.5, 0.6]);
  const band = monitorEnvBand(MON_T.samples, 'p');
  assert.deepEqual(band, [{ lo: 1e6, hi: 2e6 }, { lo: 1.2e6, hi: 1.8e6 }]);
});

// ---------------------------------------------------------------------------
// U13：模板并排对比（按名称对齐 + Δ）
// ---------------------------------------------------------------------------

test('compareSteadyRuns: 名称并集对齐、缺失侧为 NaN', () => {
  const a = {
    nodes: [{ name: 'A', pressure: 1e6 }, { name: 'B', pressure: 2e6 }],
    pipes: [{ name: 'P1', q: 0.01 }],
  };
  const b = {
    nodes: [{ name: 'A', pressure: 1.1e6 }, { name: 'C', pressure: 3e6 }],
    pipes: [{ name: 'P1', q: 0.02 }],
  };
  const r = compareSteadyRuns(a, b);
  assert.equal(r.nodes.length, 3);
  const ra = r.nodes.find((x) => x.name === 'A');
  assert.equal(ra.a, 1e6);
  assert.equal(ra.b, 1.1e6);
  assert.ok(Math.abs(ra.d + 1e5) < 1e-3);
  const rc = r.nodes.find((x) => x.name === 'C');
  assert.ok(Number.isNaN(rc.a) && rc.b === 3e6 && Number.isNaN(rc.d));
  assert.equal(r.pipes[0].name, 'P1');
  assert.deepEqual(compareSteadyRuns(null, null), { nodes: [], pipes: [] });
});

test('isSteadyComparable: 仅 steady 字典模板可并排运行', () => {
  assert.ok(isSteadyComparable('steadyFlowSimple'));
  assert.ok(isSteadyComparable('steadyNetworkMultiphase'));
  assert.ok(!isSteadyComparable('transientBlowdown'));
  assert.ok(!isSteadyComparable('twoFluidTerrainSlug'));
  assert.ok(!isSteadyComparable('nodalVogel'));
});

// ---------------------------------------------------------------------------
// U1：工况管理器（多参数 patch 组合 + 批量行归一 + 模型树挂载）
// ---------------------------------------------------------------------------

const SCEN_DICT = `time { dt 1; endTime 10; reportEvery 1; }
pipes {
  P1 { from A; to B; valveClose { start 2; end 4; closedK 1e8; } }
}
sources {
  S1 { schedule { openRate 0.5; close { start 1; end 2; }; open { start 3; end 4; }; } }
}
`;

test('applyScenarioOverrides: 多项逐条组合、错误带序号', () => {
  const out = applyScenarioOverrides(SCEN_DICT, [
    { param: 'closedK', target: 'P1', value: 2e8 },
    { param: 'scheduleOpenRate', target: 'S1', value: 0.25 },
    { param: 'endTime', target: '', value: 20 },
  ]);
  assert.match(out, /closedK 200000000/);
  assert.match(out, /openRate 0\.25/);
  assert.match(out, /endTime 20/);
  assert.doesNotMatch(out, /closedK 1e8/);
  assert.throws(
    () => applyScenarioOverrides(SCEN_DICT, [
      { param: 'closedK', target: 'P1', value: 1 },
      { param: 'bogusParam', target: '', value: 1 },
    ]),
    /第 2 项/,
  );
});

test('scenarioRowsTransient / scenarioRowsSteady: 批量行归一', () => {
  const t = scenarioRowsTransient([
    {
      ok: true,
      series: {
        steps: 3, node_names: ['A'], pipe_names: ['P1'],
        samples: [
          { t: 0, pressure: { A: 1e6 }, flow: { P1: 0.1 } },
          { t: 1, pressure: { A: 2e6 }, flow: { P1: 0.2 } },
          { t: 2, pressure: { A: 1.5e6 }, flow: { P1: 0.3 } },
        ],
      },
    },
    { ok: false, error: 'boom' },
  ]);
  assert.equal(t[0].steps, 3);
  assert.equal(t[0].tEnd, 2);
  assert.equal(t[0].pEnd, 1.5e6);
  assert.equal(t[0].qEnd, 0.3);
  assert.equal(t[0].pPeak, 2e6);
  assert.ok(t[0].series);
  assert.equal(t[1].ok, false);
  assert.match(t[1].error, /boom/);

  const s = scenarioRowsSteady([
    {
      ok: true,
      summary: {
        converged: true, iterations: 4, residual: 1e-8,
        pipes: [{ name: 'P1', q: 0.5, dp: -100 }],
      },
    },
    { ok: false, error: 'bad' },
  ]);
  assert.ok(s[0].converged);
  assert.equal(s[0].iterations, 4);
  assert.equal(s[0].qEnd, 0.5);
  assert.equal(s[0].maxDp, 100);
  assert.equal(s[1].ok, false);
});

test('buildProjectTree: 工况组挂模型树', () => {
  const tree = buildProjectTree(emptyModel(), {
    scenarios: [{ id: 'x1', name: '高产量', engine: 'transient', nOverrides: 2 }],
  });
  const g = tree.find((n) => n.id === 'group.scenarios');
  assert.ok(g);
  assert.equal(g.children.length, 1);
  assert.equal(g.children[0].id, 'scen:x1');
  assert.equal(g.children[0].kind, 'scenario');
  assert.match(g.children[0].label, /高产量 · 瞬态 · 2 项/);
  const none = buildProjectTree(emptyModel(), {});
  assert.ok(!none.find((n) => n.id === 'group.scenarios'));
});

// ---------------------------------------------------------------------------
// U4：PVT 流体配置向导（字典序列化 + chflow-pvt 报告解析，样本为真机输出）
// ---------------------------------------------------------------------------

const PVT_OUT = `[chflow] pvt: 4 components, feed methane=0.340000 ethane=0.180000 propane=0.260000 nPentane=0.220000  molarMass 0.03820544
[chflow] pvt flash T 320 p 1000000 beta 0.859839411564 single 0 iterations 10 converged 1 residual 3.93018950717e-14
[chflow] pvt liquid rho 588.35051233 M 0.0637106684302
[chflow] pvt vapour rho 13.968565418 M 0.0340478871164
[chflow] pvt x methane 0.0217948200917
[chflow] pvt y methane 0.391869947643
[chflow] pvt viscosity liquid 1.23456789 vapour 9.87e-6
[chflow] pvt surfaceTension parachors fromTable (Weinaug-Katz 1943)
[chflow] pvt surfaceTension 0.00456789012
[chflow] pvt flash T 320 p 12000000 beta 0 single 1 iterations 36 converged 1 residual 4.84945417156e-13
[chflow] pvt liquid rho 436.364143698 M 0.03820544
[chflow] pvt bubble T 320 p 8250000
[chflow] pvt table liquidDensity T 310 p 2000000 value 563.753863338
`;

test('pvtFlashCaseText: 向导规格 -> flash 字典', () => {
  const text = pvtFlashCaseText({
    temperature: 320,
    pressures: [1e6, 4e6],
    components: [{ name: 'methane', z: 0.34 }, { name: '', z: 0 }, { name: 'nPentane', z: 0.66 }],
    bubble: true, viscosity: true, surfaceTension: false,
  });
  assert.match(text, /temperature 320;/);
  assert.match(text, /pressures \(1000000 4000000\);/);
  assert.match(text, /methane 0\.34;/);
  assert.match(text, /nPentane 0\.66;/);
  assert.ok(!/^\s+;/m.test(text), '空组分不写入');
  assert.match(text, /bubblePoint yes;/);
  assert.match(text, /model lbc;/);
  assert.ok(!text.includes('surfaceTension'));
  assert.throws(() => pvtFlashCaseText({ temperature: 320, pressures: [], components: [{ name: 'a', z: 1 }] }), /压力网格/);
  assert.throws(() => pvtFlashCaseText({ temperature: 320, pressures: [1], components: [] }), /组分/);
});

test('parsePvtOutput: 真机报告全量解析', () => {
  const r = parsePvtOutput(PVT_OUT);
  assert.equal(r.nComponents, 4);
  assert.equal(r.components[0].name, 'methane');
  assert.ok(Math.abs(r.molarMass - 0.03820544) < 1e-12);
  assert.equal(r.flashes.length, 2);
  const f0 = r.flashes[0];
  assert.equal(f0.T, 320);
  assert.equal(f0.p, 1e6);
  assert.ok(Math.abs(f0.beta - 0.859839411564) < 1e-15);
  assert.equal(f0.single, false);
  assert.equal(f0.rhoL, 588.35051233);
  assert.equal(f0.rhoV, 13.968565418);
  assert.equal(f0.x.methane, 0.0217948200917);
  assert.equal(f0.y.methane, 0.391869947643);
  assert.equal(f0.muL, 1.23456789);
  assert.equal(f0.muV, 9.87e-6);
  assert.equal(f0.sigma, 0.00456789012); // 声明行不得被当成数值
  const f1 = r.flashes[1];
  assert.equal(f1.single, true);
  assert.ok(!Number.isFinite(f1.rhoV));
  assert.deepEqual(r.bubble, { T: 320, p: 8250000 });
  assert.deepEqual(r.tables[0], { prop: 'liquidDensity', T: 310, p: 2e6, value: 563.753863338 });
  assert.throws(() => parsePvtOutput('steady solved: converged'), /pvt flash/);
});

test('pvtCurvesByT + parsePvtCaseDict: 成族与模板回填', () => {
  const r = parsePvtOutput(PVT_OUT);
  const groups = pvtCurvesByT(r);
  assert.equal(groups.length, 1);
  assert.equal(groups[0].T, 320);
  assert.deepEqual(groups[0].rows.map((x) => x.p), [1e6, 1.2e7]);
  const spec = parsePvtCaseDict(
    'temperature 320.0;\npressures (1.0e6 4.0e6 1.2e7);\n' +
    'composition { methane 0.34; nPentane 0.66; }\n' +
    'solve { bubblePoint yes; }\nviscosity { model lbc; }\n');
  assert.equal(spec.temperature, 320);
  assert.deepEqual(spec.pressures, [1e6, 4e6, 1.2e7]);
  assert.deepEqual(spec.components, [{ name: 'methane', z: 0.34 }, { name: 'nPentane', z: 0.66 }]);
  assert.ok(spec.bubble && spec.viscosity && !spec.surfaceTension);
});

// ---------------------------------------------------------------------------
// U6（§5.152）：Green-Gauss 涡量/Q 判据 + 等值面片提取
// ---------------------------------------------------------------------------

/** nx x ny x nz unit-cube hex mesh in the same faces layout as VtkWriter. */
function hexMesh(nx, ny, nz) {
  const P = (i, j, k) => i + (nx + 1) * (j + (ny + 1) * k);
  const points = [];
  for (let k = 0; k <= nz; k += 1)
    for (let j = 0; j <= ny; j += 1)
      for (let i = 0; i <= nx; i += 1)
        points.push({ x: i, y: j, z: k });
  const cells = [];
  const cellFaces = [];
  for (let k = 0; k < nz; k += 1)
    for (let j = 0; j < ny; j += 1)
      for (let i = 0; i < nx; i += 1) {
        const c = [P(i, j, k), P(i + 1, j, k), P(i + 1, j + 1, k), P(i, j + 1, k),
          P(i, j, k + 1), P(i + 1, j, k + 1), P(i + 1, j + 1, k + 1), P(i, j + 1, k + 1)];
        cells.push(c);
        cellFaces.push([
          [c[0], c[1], c[5], c[4]], [c[3], c[2], c[6], c[7]],
          [c[0], c[1], c[2], c[3]], [c[4], c[5], c[6], c[7]],
          [c[0], c[3], c[7], c[4]], [c[1], c[2], c[6], c[5]],
        ]);
      }
  return { points, cells, cellFaces };
}

test('vorticityAndQ: 刚性旋转场的中心单元给出精确 ω_z=2、Q=1', () => {
  const { points, cells, cellFaces } = hexMesh(3, 3, 1);
  const cellU = cells.map((c) => {
    const cy = points[c[0]].y + 0.5;
    const cx = points[c[0]].x + 0.5;
    return [-cy, cx, 0];
  });
  const r = vorticityAndQ(points, cells, cellFaces, cellU);
  assert.ok(r);
  const mid = 4; // k=0, j=1, i=1 -> all four side faces internal
  assert.ok(Math.abs(r.omega[mid][2] - 2) < 1e-9, `wz=${r.omega[mid][2]}`);
  assert.ok(Math.abs(r.omega[mid][0]) < 1e-9 && Math.abs(r.omega[mid][1]) < 1e-9);
  assert.ok(Math.abs(r.mag[mid] - 2) < 1e-9);
  assert.ok(Math.abs(r.q[mid] - 1) < 1e-9, `q=${r.q[mid]}`);
  // 纯应变场 U=(x,-y,0)：无涡但 Q<0（应变率主导）
  const strain = cells.map((c) => {
    const cy = points[c[0]].y + 0.5;
    const cx = points[c[0]].x + 0.5;
    return [cx, -cy, 0];
  });
  const s = vorticityAndQ(points, cells, cellFaces, strain);
  assert.ok(Math.abs(s.omega[mid][2]) < 1e-9);
  assert.ok(Math.abs(s.mag[mid]) < 1e-9);
  assert.ok(Math.abs(s.q[mid] + 1) < 1e-9, `q=${s.q[mid]}`);
});

test('faceAdjacency: 中心单元 4 个内面 + 2 个边界（顶底）', () => {
  const { cells, cellFaces } = hexMesh(3, 3, 1);
  const adj = faceAdjacency(cells, cellFaces);
  assert.ok(adj);
  const internal = adj.nbrs[4].filter((n) => n >= 0).length;
  assert.equal(internal, 4);
  assert.equal(adj.nbrs[0].filter((n) => n >= 0).length, 2); // corner cell
});

test('vorticityAndQ: 无 faces 数组返回 null（UI 回退 |U|）', () => {
  const { points, cells } = hexMesh(2, 1, 1);
  assert.equal(vorticityAndQ(points, cells, null, cells.map(() => [1, 0, 0])), null);
});

test('isoSurfacePatches: 跨越阈值的内面成壳、位置按线性穿越分数偏移', () => {
  const { points, cells, cellFaces } = hexMesh(3, 3, 1);
  const values = cells.map((c) => points[c[0]].x + 0.5); // p = 单元中心 x
  const r = isoSurfacePatches(points, cells, cellFaces, values, 1.25);
  assert.equal(r.patches.length, 3); // 三行 j 各一片，位于 i=0/i=1 之间
  assert.deepEqual(r.f, [0.75, 0.75, 0.75]);
  for (const ring of r.patches) {
    const cx = ring.reduce((s, p) => s + p[0], 0) / ring.length;
    assert.ok(Math.abs(cx - 1.25) < 1e-12, `centroid x=${cx}`);
  }
  // 阈值在全部值之外 -> 无片；恰等于某 cell 值 -> 不算跨越
  assert.equal(isoSurfacePatches(points, cells, cellFaces, values, 0.25).patches.length, 0);
  assert.equal(isoSurfacePatches(points, cells, cellFaces, values, 1.5).patches.length, 0);
});

// ---------------------------------------------------------------------------
// U9（§5.153）：场 CSV 导出 + 剖面 CSV
// ---------------------------------------------------------------------------

test('fieldsToCsv: 每 cell 一行（质心 + p + U + |U|）', () => {
  const { points, cells, cellFaces } = hexMesh(3, 3, 1);
  const cellU = cells.map((c) => [1, 2, 2]);
  const cellP = cells.map((c) => 100 + points[c[0]].x);
  const lines = fieldsToCsv(points, cells, cellU, cellP).trim().split('\n');
  assert.equal(lines[0], 'cell,cx,cy,cz,p,Ux,Uy,Uz,magU');
  assert.equal(lines.length, 10); // 1 表头 + 9 单元
  const first = lines[1].split(',').map(Number);
  assert.deepEqual(first.slice(1, 4), [0.5, 0.5, 0.5]);
  assert.equal(first[4], 100);
  assert.equal(first[8], 3);
});

test('planeProfileCsv: 剖面取落在平面上的面所属单元；无 faces 返回 null', () => {
  const { points, cells, cellFaces } = hexMesh(3, 3, 1);
  const cellU = cells.map(() => [0, 0, 1]);
  const cellP = cells.map(() => 1);
  // 法向 x、frac=1/3 -> 平面 x=1：i=0 的右面与 i=1 的左面都贴住 -> 6 个单元
  const csv = planeProfileCsv(points, cells, cellFaces, cellU, cellP, [1, 0, 0], 1 / 3);
  const rows = csv.trim().split('\n').slice(1).map((l) => l.split(',').map(Number));
  assert.equal(rows.length, 6);
  for (const r of rows) assert.ok(r[0] === 0.5 || r[0] === 1.5, `cx=${r[0]}`);
  // 顶面 frac=1：单层网格所有 cell 都有顶面
  const top = planeProfileCsv(points, cells, cellFaces, cellU, cellP, [0, 0, 1], 1);
  assert.equal(top.trim().split('\n').length - 1, 9);
  assert.equal(planeProfileCsv(points, cells, null, cellU, cellP, [1, 0, 0], 0.5), null);
});
