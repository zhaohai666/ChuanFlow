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
  hitTestNode,
  networkToCaseText,
} from '../frontend/dist/editor.js';

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
