/**
 * Built-in case template metadata for the 模板库 menu. The actual dictionary
 * text always comes from kernel-verified case files (CTest regression cases)
 * resolved by the Rust `case_templates` command — a template is by
 * construction a case the engine accepts, never a hand-typed guess.
 *
 * `tab` is the workspace tab a template opens into; `slug`/`coupled` mark the
 * templates that drive the animation / coupling views instead of the editor.
 * 3D templates additionally carry an `engine` key naming the member of the
 * 3D engine family (Rust `three_d_engines` registry, §5.181) that must run
 * them — a template whose dict has no piso{} block must never be submitted
 * to chflow-piso-foam (and vice versa).
 * Pure metadata + helper mapping, covered by node --test.
 */

/**
 * E1 (§5.181): the 3D engine family, mirrored from the Rust
 * THREE_D_FAMILY registry. Each row pairs a persistent engine key
 * (chflow.engines.v1 style) with the binary name and the template that
 * drives it; the Rust unit test `three_d_family_cases_are_all_registered_
 * templates` keeps both sides from drifting.
 */
export const THREE_D_FAMILY = [
  { key: 'threeD', bin: 'chflow-piso-foam', template: 'startedChannel', label: 'PISO 瞬态通道启动' },
  { key: 'threeDSimple', bin: 'chflow-simple-foam', template: 'poiseuilleChannel', label: 'SIMPLE 稳态泊肃叶通道' },
  { key: 'threeDVof', bin: 'chflow-vof', template: 'vofDroplet', label: 'VOF 液滴界面捕捉' },
  { key: 'threeDCompressible', bin: 'chflow-3d-compressible', template: 'compressibleNsSod', label: '可压缩 NS（Sod 激波管）' },
  { key: 'threeDScalar', bin: 'chflow-scalar-transport', template: 'advectionDiffusion', label: '标量对流-扩散' },
  { key: 'threeDDpm', bin: 'chflow-3d-dpm', template: 'dpmChannel', label: 'DPM 颗粒通道（单向耦合）' },
  { key: 'threeDP1', bin: 'chflow-3d-p1', template: 'radiationP1Slab', label: 'P1 辐射灰色介质板' },
  { key: 'threeDCombustion', bin: 'chflow-3d-combustion', template: 'combustionFixedOxidiser', label: '一步化学燃烧' },
  { key: 'threeDLes', bin: 'chflow-3d-les', template: 'lesSgs', label: 'LES 亚格子模型' },
  { key: 'threeDMelting', bin: 'chflow-3d-melting', template: 'meltingFront', label: '熔化前沿（相变）' },
  { key: 'threeDCavitation', bin: 'chflow-3d-cavitation', template: 'cavitationVenturi', label: 'Venturi 空化' },
];

const THREE_D_BY_TEMPLATE = new Map(
  THREE_D_FAMILY.map((e) => [e.template, e.key]),
);

export const TEMPLATE_IDS = [
  'steadyNetwork',
  'steadyEquipmentValve',
  'steadyEquipmentPump',
  'steadyNetworkMultiphase',
  'transientBlowdown',
  'twoFluidTerrainSlug',
  'coupled1d3d',
  'nodalVogel',
  'pvtFlash',
  'startedChannel',
  // E1 (§5.181): 3D family templates (one per wired 3D app)
  ...THREE_D_FAMILY.filter((e) => e.key !== 'threeD').map((e) => e.template),
];

/** Which tab a template opens into and how it feeds the workspace. */
export function templateTarget(id) {
  if (id === 'twoFluidTerrainSlug') return { tab: 'slug', kind: 'slug' };
  if (id === 'coupled1d3d') return { tab: 'coupled', kind: 'coupled' };
  if (id === 'nodalVogel') return { tab: 'nodal', kind: 'nodal' };
  if (id === 'pvtFlash') return { tab: 'pvt', kind: 'pvt' };
  if (THREE_D_BY_TEMPLATE.has(id)) return { tab: 'dict', kind: '3d', engine: THREE_D_BY_TEMPLATE.get(id) };
  if (id === 'transientBlowdown') return { tab: 'canvas', kind: 'dict' };
  return { tab: 'canvas', kind: 'dict' };
}

/** Steady dictionary templates can be batch-run against the current model. */
export function isSteadyComparable(id) {
  return /^steady/.test(id);
}

/**
 * Templates whose dictionaries carry blocks the drag-editor model cannot
 * round-trip (blackoil / PVT compositional). Loading one still runs and
 * renders exactly (the raw dict is the run authority), but model-fidelity
 * features — the drag editor and the parameter sweep — are disabled for it
 * and the UI says so, instead of silently dropping physics.
 */
export const PARTIAL_FIDELITY_IDS = ['steadyNetworkMultiphase'];

export function isPartialFidelity(id) {
  return PARTIAL_FIDELITY_IDS.includes(id);
}
