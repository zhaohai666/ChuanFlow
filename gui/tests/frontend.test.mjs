import test from 'node:test';
import assert from 'node:assert/strict';
import {
  validateCaseText,
  toBarHeights,
  toLinePoints,
} from '../frontend/dist/chart.js';

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
