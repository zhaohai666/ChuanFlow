/**
 * E2b tests for the seven structured dictionary blocks:
 *   wax / pigging / sand / acoustics (top-level two-fluid blocks, flowAssurance.js)
 *   correlations / rheology         (per-pipe steady keys, dictBlocks.js)
 *   blackoil                         (nested fluid.blackoil, dictBlocks.js)
 *
 * Three kinds of case per block:
 *   1. schema alignment — the spec's keys are the kernel's own spellings,
 *      hard-coded below; a renamed kernel key turns the suite red.
 *   2. round-trip — render the fixture to dictionary text, parse it back with
 *      the GUI's own dict parser, every declared field survives with the
 *      correct type.
 *   3. rejection — a missing required key or a kernel mutex/range conflict is
 *      refused up front with an explicit message (asserted as a throw, not as
 *      the existence of a helper).
 */

import test from 'node:test';
import assert from 'node:assert/strict';
import {
  FA_BLOCKS,
  flowAssuranceBlockLines,
  parseFlowAssurance,
  applyFlowAssurance,
} from '../frontend/dist/flowAssurance.js';
import {
  PIPE_BLOCKS,
  BLACKOIL_SPEC,
  applyPipeBlock,
  parsePipeEntry,
  applyBlackoil,
  parseBlackoil,
} from '../frontend/dist/dictBlocks.js';

// Kernel keys copied verbatim from the source reads:
//   wax         twoFluid.cpp:3481-3498
//   pigging     twoFluid.cpp:3866-3902
//   sand        twoFluid.cpp:3589-3605 (guards 454-464)
//   acoustics   twoFluid.cpp:3904-4004
//   correlations network.cpp:3431-3469 (scope guards 1228-1269)
//   rheology    network.cpp:3813-3836
//   blackoil     network.cpp:3048-3053 (BlackOil.h:82-89)
const KERNEL_KEYS = {
  wax: ['meltingTemperature', 'fusionEnthalpy', 'waxFraction', 'growthRate',
    'dissolutionRate', 'latentHeat', 'initial', 'feedback', 'depositRoughness'],
  pigging: ['enable', 'startX', 'endX', 'pigVelocity', 'blockage',
    'pushEfficiency', 'frontCapacity'],
  sand: ['terminalVelocity', 'criticalVelocity', 'erodeRate', 'initial',
    'bedInitial', 'feedback'],
  acoustics: ['enable', 'initialGauge', 'soundSpeed', 'aModel',
    'liquidSoundSpeed', 'cfl', 'openEnd'],
  correlations: ['correlation', 'flowPattern', 'hbGriffith'],
  rheology: ['flowIndex', 'consistency', 'yieldStress'],
  blackoil: ['gasGravity', 'oilDensity', 'gor', 'waterCut', 'temperature',
    'surfaceTension'],
};

// Valid boundary/default fixtures, one per block.
const FIXTURES = {
  wax: {
    meltingTemperature: 323.75, fusionEnthalpy: 54890, waxFraction: 0.5,
    growthRate: 0.01, dissolutionRate: 0.02, latentHeat: 0, initial: 0,
    feedback: 'false', depositRoughness: 0,
  },
  pigging: {
    enable: 'true', startX: 0, endX: 20, pigVelocity: 2.0, blockage: 0.5,
    pushEfficiency: 0.6, frontCapacity: 8,
  },
  sand: {
    terminalVelocity: 0.3, criticalVelocity: 1.0, erodeRate: 0, initial: 0.2,
    bedInitial: 0, feedback: 'true',
  },
  acoustics: {
    enable: 'true', initialGauge: 1e6, soundSpeed: 0, aModel: 'twoPhase',
    cfl: 0.5, openEnd: 'outlet',
  },
  correlations: { correlation: 'hagedornBrown', hbGriffith: 'true' },
  rheology: { flowIndex: 0.5, consistency: 0.0447, yieldStress: 0.05 },
  blackoil: {
    gasGravity: 0.75, oilDensity: 850, gor: 0, waterCut: 0, temperature: 340,
    surfaceTension: 0.025,
  },
};

const FA_BLOCK_NAMES = ['wax', 'pigging', 'sand', 'acoustics'];
const PIPE_BLOCK_NAMES = ['correlations', 'rheology'];

const BASE_TWOFLUID = [
  'fluid { rhoL 1000; rhoG 20; muL 0.01; muG 1e-5; }',
  'pipe { length 10; diameter 0.1; sections 10; }',
  'closure { model wallis; }',
  '',
].join('\n');

const BASE_STEADY = [
  'fluid { rho 850; mu 0.1; }',
  'nodes { A { pressure 1.005e6; } B { pressure 1.0e6; } }',
  'pipes {',
  '    P1 { from A; to B; length 100; diameter 0.05; }',
  '}',
  '',
].join('\n');

// ---------------------------------------------------------------------------
// 1. Schema alignment: spec field keys must equal the kernel keys, verbatim.
// ---------------------------------------------------------------------------

for (const name of FA_BLOCK_NAMES) {
  test(`schema: ${name} 键名与内核逐字对齐`, () => {
    const spec = FA_BLOCKS[name];
    assert.ok(spec, `缺少块 ${name} 的 schema`);
    assert.deepEqual(spec.fields.map((f) => f.key), KERNEL_KEYS[name]);
  });
}

for (const name of PIPE_BLOCK_NAMES) {
  test(`schema: ${name} 键名与内核逐字对齐`, () => {
    assert.deepEqual(PIPE_BLOCKS[name].fields.map((f) => f.key), KERNEL_KEYS[name]);
  });
}

test('schema: blackoil 键名与内核逐字对齐', () => {
  assert.deepEqual(BLACKOIL_SPEC.fields.map((f) => f.key), KERNEL_KEYS.blackoil);
});

// ---------------------------------------------------------------------------
// 2. Round-trip: fixture -> dict text -> GUI parser -> values survive, typed.
// ---------------------------------------------------------------------------

for (const name of FA_BLOCK_NAMES) {
  test(`round-trip: ${name} 生成文本经解析往返不丢字段`, () => {
    const fixture = FIXTURES[name];
    const form = { [name]: { enabled: true, values: { ...fixture } } };
    const out = applyFlowAssurance(BASE_TWOFLUID, form);
    // every declared key appears literally in the emitted dictionary text
    for (const k of KERNEL_KEYS[name]) {
      if (fixture[k] !== undefined) assert.match(out, new RegExp(`\\b${k}\\b`));
    }
    const parsed = parseFlowAssurance(out)[name].values;
    for (const [k, raw] of Object.entries(fixture)) {
      assert.ok(parsed[k] !== undefined, `${name}.${k} 往返后丢失`);
      if (FA_BLOCKS[name].fields.find((f) => f.key === k).kind === 'num') {
        assert.ok(Number.isFinite(Number(parsed[k])), `${name}.${k} 不是有限数值`);
      } else {
        assert.equal(parsed[k], String(raw));
      }
    }
  });
}

for (const name of PIPE_BLOCK_NAMES) {
  test(`round-trip: ${name} 管段键往返不丢字段`, () => {
    const fixture = FIXTURES[name];
    const out = applyPipeBlock(BASE_STEADY, 'P1', name, { ...fixture });
    for (const k of KERNEL_KEYS[name]) {
      if (fixture[k] !== undefined) assert.match(out, new RegExp(`\\b${k}\\b`));
    }
    const parsed = parsePipeEntry(out, 'P1');
    for (const [k, raw] of Object.entries(fixture)) {
      assert.ok(parsed[k] !== undefined, `${name}.${k} 往返后丢失`);
      if (PIPE_BLOCKS[name].fields.find((f) => f.key === k).kind === 'num') {
        assert.ok(Number.isFinite(Number(parsed[k])), `${name}.${k} 不是有限数值`);
      } else {
        assert.equal(parsed[k], String(raw));
      }
    }
    // unrelated pipe geometry survives verbatim
    assert.equal(parsed.from, 'A');
    assert.equal(parsed.to, 'B');
    assert.equal(parsed.length, '100');
  });

  test(`removal: ${name} 置空后管段键被清除、几何保留`, () => {
    const seeded = applyPipeBlock(BASE_STEADY, 'P1', name, { ...FIXTURES[name] });
    const out = applyPipeBlock(seeded, 'P1', name, null);
    const parsed = parsePipeEntry(out, 'P1');
    for (const k of KERNEL_KEYS[name]) assert.equal(parsed[k], undefined);
    assert.equal(parsed.from, 'A');
    assert.equal(parsed.length, '100');
  });
}

test('round-trip: blackoil 嵌套块往返不丢字段', () => {
  const fixture = FIXTURES.blackoil;
  const out = applyBlackoil('fluid { rho 850; }\n', { ...fixture });
  for (const k of KERNEL_KEYS.blackoil) assert.match(out, new RegExp(`\\b${k}\\b`));
  const parsed = parseBlackoil(out);
  for (const [k, raw] of Object.entries(fixture)) {
    assert.ok(parsed[k] !== undefined, `blackoil.${k} 往返后丢失`);
    if (BLACKOIL_SPEC.fields.find((f) => f.key === k).kind === 'num') {
      assert.ok(Number.isFinite(Number(parsed[k])), `blackoil.${k} 不是有限数值`);
    } else {
      assert.equal(parsed[k], String(raw));
    }
  }
});

test('removal: blackoil 置空后子块清除，fluid 其他键保留', () => {
  const seeded = applyBlackoil('fluid { rho 850; mu 0.1; }\n', { ...FIXTURES.blackoil });
  const out = applyBlackoil(seeded, null);
  assert.equal(parseBlackoil(out), null);
  assert.match(out, /\brho\s+850/);
});

// ---------------------------------------------------------------------------
// 3. Rejection: missing required / range / mutex conflicts are refused loud.
// ---------------------------------------------------------------------------

function expectReject(fn, fragment) {
  assert.throws(fn, (err) => {
    assert.ok(err.message.includes(fragment), `提示应包含「${fragment}」，got: ${err.message}`);
    return true;
  });
}

test('拒绝: wax 越界分数 / 负生长速率', () => {
  expectReject(() => flowAssuranceBlockLines('wax', { ...FIXTURES.wax, waxFraction: 1.5 }),
    'waxFraction');
  expectReject(() => flowAssuranceBlockLines('wax', { ...FIXTURES.wax, growthRate: -1 }),
    'growthRate');
});

test('拒绝: pigging 缺 frontCapacity / 堵塞比越界 / 启用零速度 / 端点倒置', () => {
  expectReject(() => flowAssuranceBlockLines('pigging', {
    ...FIXTURES.pigging, frontCapacity: undefined,
  }), 'frontCapacity');
  expectReject(() => flowAssuranceBlockLines('pigging', { ...FIXTURES.pigging, blockage: 1 }),
    'blockage');
  expectReject(() => flowAssuranceBlockLines('pigging', {
    ...FIXTURES.pigging, enable: 'true', pigVelocity: 0,
  }), 'pigVelocity');
  expectReject(() => flowAssuranceBlockLines('pigging', {
    ...FIXTURES.pigging, startX: 20, endX: 10,
  }), 'endX');
});

test('拒绝: sand 缺必填速度 / 临界速度非正 / 床层分数越界', () => {
  expectReject(() => flowAssuranceBlockLines('sand', {
    ...FIXTURES.sand, terminalVelocity: undefined,
  }), 'terminalVelocity');
  expectReject(() => flowAssuranceBlockLines('sand', {
    ...FIXTURES.sand, criticalVelocity: 0,
  }), 'criticalVelocity');
  expectReject(() => flowAssuranceBlockLines('sand', { ...FIXTURES.sand, bedInitial: 1 }),
    'bedInitial');
});

test('拒绝: acoustics wood 与 soundSpeed 互斥 / wood 缺 aL / CFL 越界 / 未知 openEnd', () => {
  expectReject(() => flowAssuranceBlockLines('acoustics', {
    ...FIXTURES.acoustics, aModel: 'wood', soundSpeed: 300, liquidSoundSpeed: 1500,
  }), '互斥');
  expectReject(() => flowAssuranceBlockLines('acoustics', {
    ...FIXTURES.acoustics, aModel: 'wood', soundSpeed: 0,
  }), 'liquidSoundSpeed');
  expectReject(() => flowAssuranceBlockLines('acoustics', { ...FIXTURES.acoustics, cfl: 0.6 }),
    'cfl');
  expectReject(() => flowAssuranceBlockLines('acoustics', {
    ...FIXTURES.acoustics, openEnd: 'side',
  }), 'openEnd');
});

test('拒绝: correlations 定律与流型互斥 / hbGriffith 配错定律 / 非法流型索引', () => {
  expectReject(() => applyPipeBlock(BASE_STEADY, 'P1', 'correlations', {
    correlation: 'beggsBrill', flowPattern: 2,
  }), '互斥');
  expectReject(() => applyPipeBlock(BASE_STEADY, 'P1', 'correlations', {
    correlation: 'beggsBrill', hbGriffith: 'true',
  }), 'hagedornBrown');
  expectReject(() => applyPipeBlock(BASE_STEADY, 'P1', 'correlations', {
    flowPattern: 1,
  }), 'flowPattern');
});

test('拒绝: rheology n>0 缺 K / n=0 带 K / 负屈服应力', () => {
  expectReject(() => applyPipeBlock(BASE_STEADY, 'P1', 'rheology', {
    flowIndex: 0.5, consistency: 0,
  }), 'consistency');
  expectReject(() => applyPipeBlock(BASE_STEADY, 'P1', 'rheology', {
    flowIndex: 0, consistency: 0.1,
  }), 'flowIndex');
  expectReject(() => applyPipeBlock(BASE_STEADY, 'P1', 'rheology', {
    yieldStress: -0.05,
  }), 'yieldStress');
});

test('拒绝: blackoil 油密度非正 / 含水率越界', () => {
  expectReject(() => applyBlackoil('fluid { rho 850; }\n', {
    ...FIXTURES.blackoil, oilDensity: 0,
  }), 'oilDensity');
  expectReject(() => applyBlackoil('fluid { rho 850; }\n', {
    ...FIXTURES.blackoil, waterCut: 1,
  }), 'waterCut');
});

test('拒绝: 未知块 id 被拒绝（不静默写坏字典）', () => {
  assert.throws(() => applyPipeBlock(BASE_STEADY, 'P1', 'nope', {}), /unknown/);
  assert.throws(() => flowAssuranceBlockLines('nope', {}), /unknown/);
});
