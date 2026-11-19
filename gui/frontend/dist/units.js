/**
 * Unit systems for the engineering workspace: SI (国际单位) and oilfield
 * (油田英制). The case dictionary is ALWAYS SI — the unit system is a
 * display layer only, exactly like the commercial tools: every input box,
 * label and chart axis converts on the fly and converts back on apply.
 *
 * Pure functions, no DOM and no Tauri imports, so `node --test` covers the
 * conversion round trips bit-for-bit against the kernel-facing values.
 */

export const UNIT_SYSTEMS = {
  SI: 'SI',
  OILFIELD: 'oilfield',
};

/** Conversion table per quantity family. `to` converts SI -> display. */
const FAMILIES = {
  length: {
    label: { SI: 'm', oilfield: 'ft' },
    to: {
      SI: (v) => v,
      oilfield: (v) => v * 3.28083989501312,
    },
    from: {
      SI: (v) => v,
      oilfield: (v) => v / 3.28083989501312,
    },
  },
  pressure: {
    label: { SI: 'Pa', oilfield: 'psi' },
    to: {
      SI: (v) => v,
      oilfield: (v) => v / 6894.757293168361,
    },
    from: {
      SI: (v) => v,
      oilfield: (v) => v * 6894.757293168361,
    },
  },
  rate: {
    label: { SI: 'm³/s', oilfield: 'bbl/d' },
    to: {
      SI: (v) => v,
      oilfield: (v) => v * 543439.668672, // 6.28981098 bbl/m³ x 86400 s/d
    },
    from: {
      SI: (v) => v,
      oilfield: (v) => v / 543439.668672,
    },
  },
  temperature: {
    label: { SI: 'K', oilfield: '°F' },
    to: {
      SI: (v) => v,
      oilfield: (v) => ((v - 273.15) * 9) / 5 + 32,
    },
    from: {
      SI: (v) => v,
      oilfield: (v) => ((v - 32) * 5) / 9 + 273.15,
    },
  },
  viscosity: {
    label: { SI: 'Pa·s', oilfield: 'cP' },
    to: {
      SI: (v) => v,
      oilfield: (v) => v * 1000,
    },
    from: {
      SI: (v) => v,
      oilfield: (v) => v / 1000,
    },
  },
  density: {
    label: { SI: 'kg/m³', oilfield: 'lb/ft³' },
    to: {
      SI: (v) => v,
      oilfield: (v) => v * 0.062427961,
    },
    from: {
      SI: (v) => v,
      oilfield: (v) => v / 0.062427961,
    },
  },
};

export const QUANTITY_FAMILIES = Object.keys(FAMILIES);

/** Convert a kernel-side (SI) value into the display system. */
export function toDisplay(value, family, system) {
  const f = FAMILIES[family];
  if (!f) throw new Error(`unknown quantity family: ${family}`);
  if (!Number.isFinite(Number(value))) return NaN;
  return f.to[system](Number(value));
}

/** Convert a display-system value back into kernel-side (SI). */
export function fromDisplay(value, family, system) {
  const f = FAMILIES[family];
  if (!f) throw new Error(`unknown quantity family: ${family}`);
  if (!Number.isFinite(Number(value))) return NaN;
  return f.from[system](Number(value));
}

/** Unit label for a family in the given system ('m', 'psi', 'bbl/d', ...). */
export function unitLabel(family, system) {
  const f = FAMILIES[family];
  if (!f) throw new Error(`unknown quantity family: ${family}`);
  return f.label[system];
}

/**
 * Format a kernel-side (SI) value for display: converted, compact
 * significant digits, with the unit appended.
 */
export function fmtDisplay(value, family, system, sigDigits = 6) {
  const v = toDisplay(value, family, system);
  if (!Number.isFinite(v)) return '—';
  return `${Number(v.toPrecision(sigDigits))} ${unitLabel(family, system)}`;
}

/**
 * Parse a display-system text input into a kernel-side (SI) number.
 * Returns NaN for unparsable input so callers can reject loudly.
 */
export function parseInput(text, family, system) {
  const v = Number(String(text).trim());
  if (!Number.isFinite(v)) return NaN;
  return fromDisplay(v, family, system);
}
