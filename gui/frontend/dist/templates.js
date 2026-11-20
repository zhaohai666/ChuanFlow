/**
 * Built-in case template metadata for the 模板库 menu. The actual dictionary
 * text always comes from kernel-verified case files (CTest regression cases)
 * resolved by the Rust `case_templates` command — a template is by
 * construction a case the engine accepts, never a hand-typed guess.
 *
 * `tab` is the workspace tab a template opens into; `slug`/`coupled` mark the
 * templates that drive the animation / coupling views instead of the editor.
 * Pure metadata + helper mapping, covered by node --test.
 */

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
];

/** Which tab a template opens into and how it feeds the workspace. */
export function templateTarget(id) {
  if (id === 'twoFluidTerrainSlug') return { tab: 'slug', kind: 'slug' };
  if (id === 'coupled1d3d') return { tab: 'coupled', kind: 'coupled' };
  if (id === 'nodalVogel') return { tab: 'nodal', kind: 'nodal' };
  if (id === 'pvtFlash') return { tab: 'pvt', kind: 'pvt' };
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
