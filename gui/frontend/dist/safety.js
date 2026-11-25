/**
 * Flow-assurance / well-control readouts for the workbench.
 *
 * Two kernel report families are folded here, both of them *read-only* on the
 * solver side:
 *
 *   §5.185  `chflow-1d-transient` prints one `reservoir <name> node <n> q …
 *           gasq … vol … gasvol … ub … dp … killmud … sidpp` line per frame
 *           per declared reservoir. The numbers are the PI inflow, the gas
 *           share it carries, the volumes it has pushed into the node, the
 *           underbalance flag and the two static kill readings.
 *   §5.187  `chflow-1d-twofluid` prints `risk erosion …`, `risk hydrate …`,
 *           `risk corrosion …` per frame: API RP 14E velocity ratio, the
 *           Clausius-Clapeyron subcooling grade and the CO2/H2S partials.
 *
 * Nothing here recomputes physics; every displayed number is copied from the
 * kernel line, and the labels are the kernel's own band indices.
 */

const H = {
  num: (s) => {
    const v = Number(s);
    return Number.isFinite(v) ? v : NaN;
  },
  fmt: (v, digits = 4) => (Number.isFinite(v) ? v.toPrecision(digits) : '—'),
  si: (v, unit) => (Number.isFinite(v) ? `${v.toPrecision(4)} ${unit}` : '—'),
};

/** Kernel band index → Chinese label (0/1/2/3 as twoFluid.cpp grades them). */
export const HYDRATE_GRADE = ['无', '轻', '中', '重'];
export const CO2_BAND = ['低', '中', '高'];

/**
 * Fold a `reservoir …` line into its fields. Exported so tests can pin the
 * grammar without going through parseTransientOutput.
 */
export function parseReservoirLine(line) {
  const m = /^reservoir\s+(\S+)\s+node\s+(\S+)\s+q\s+(\S+)\s+gasq\s+(\S+)\s+vol\s+(\S+)\s+gasvol\s+(\S+)\s+ub\s+(\d+)\s+dp\s+(\S+)\s+killmud\s+(\S+)\s+sidpp\s+(\S+)$/
    .exec(String(line).trim());
  if (!m) return null;
  return {
    name: m[1],
    node: m[2],
    q: H.num(m[3]),
    gasq: H.num(m[4]),
    vol: H.num(m[5]),
    gasvol: H.num(m[6]),
    underbalanced: Number(m[7]) !== 0,
    dp: H.num(m[8]),
    killMudWeight: H.num(m[9]),
    sidpp: H.num(m[10]),
  };
}

/**
 * The well-control panel of one transient run: the LAST reported frame, plus
 * the extremum of the gas share over the whole history (the kick grows even
 * when the last frame happens to be after shut-in).
 *
 * `series` is what parseTransientOutput returns; a series without reservoirs
 * (the common case) yields `null` so the panel can stay hidden.
 */
export function wellControlPanel(series) {
  if (!series || !series.samples || !series.samples.length) return null;
  const names = series.reservoirNames || [];
  if (!names.length) return null;
  const last = series.samples.at(-1);
  const entries = [];
  let peakGasShare = 0;
  let peakGasShareAt = NaN;
  let anyUnderbalanced = false;
  for (const name of names) {
    const r = last.reservoir?.[name];
    if (!r) continue;
    if (r.underbalanced) anyUnderbalanced = true;
    for (const s of series.samples) {
      const e = s.reservoir?.[name];
      if (!e || !(e.q > 0)) continue;
      const share = e.gasq / e.q;
      if (share > peakGasShare) {
        peakGasShare = share;
        peakGasShareAt = s.t;
      }
    }
    entries.push({ ...r });
  }
  if (!entries.length) return null;

  const rows = [];
  for (const e of entries) {
    const share = e.q > 0 ? e.gasq / e.q : 0;
    rows.push(
      { label: `${e.name}·节点${e.node} 流量 q`, value: H.si(e.q, 'm³/s'), raw: e.q },
      { label: `${e.name}·气侵 qg`, value: H.si(e.gasq, 'm³/s'), raw: e.gasq,
        note: `气相份额 ${(share * 100).toFixed(2)}%` },
      { label: `${e.name}·累计侵入`, value: `${H.fmt(e.vol)} m³`,
        note: `其中气 ${H.fmt(e.gasvol)} m³`, raw: e.vol },
      { label: `${e.name}·欠平衡`, value: e.underbalanced ? '是（储层向井筒流入）' : '否（过平衡，q≡0）',
        flag: e.underbalanced ? 'warn' : 'ok' },
      { label: `${e.name}·驱动压差`, value: `${H.fmt(e.dp)} Pa`, raw: e.dp },
      { label: `${e.name}·压井密度`, value: `${H.fmt(e.killMudWeight, 6)} kg/m³`, raw: e.killMudWeight,
        note: '静态读数 p_res/(g·depth)' },
      { label: `${e.name}·SIDPP`, value: `${H.fmt(e.sidpp, 6)} Pa`, raw: e.sidpp,
        note: '静态代数 p_res − ρ_mud·g·depth' },
    );
  }
  return {
    t: last.t,
    entries,
    rows,
    underbalanced: anyUnderbalanced,
    peakGasShare,
    peakGasShareAt,
    summary:
      `${entries.length} 个储层，t=${last.t} s：` +
      entries.map((e) => `${e.name} q=${H.fmt(e.q, 3)} m³/s`).join('；') +
      (peakGasShare > 0
        ? `；峰值气相份额 ${(peakGasShare * 100).toFixed(2)}% @ t=${H.fmt(peakGasShareAt, 3)} s`
        : ''),
  };
}

/** One `risk …` line → its structured entry. */
export function parseRiskLine(line) {
  const l = String(line).trim();
  let m = /^risk erosion seg (\d+) vm (\S+) rhom (\S+) velim (\S+) ratio (\S+) over (\d+)(?:\s+(\S+))?/
    .exec(l);
  if (m) {
    return {
      kind: 'erosion',
      seg: Number(m[1]),
      vm: H.num(m[2]),
      rhoM: H.num(m[3]),
      velim: H.num(m[4]),
      ratio: H.num(m[5]),
      over: Number(m[6]),
      overFlag: (m[7] || '').length > 0,
    };
  }
  m = /^risk hydrate seg (\d+) p (\S+) T (\S+) Teq (\S+) dTsub (\S+) grade (\d+)/.exec(l);
  if (m) {
    return {
      kind: 'hydrate',
      seg: Number(m[1]),
      p: H.num(m[2]),
      T: H.num(m[3]),
      teq: H.num(m[4]),
      subcool: H.num(m[5]),
      grade: Number(m[6]),
    };
  }
  m = /^risk corrosion co2 (\S+) h2s (\S+) band (\d+) sour (\d+)/.exec(l);
  if (m) {
    return {
      kind: 'corrosion',
      co2Partial: H.num(m[1]),
      h2sPartial: H.num(m[2]),
      band: Number(m[3]),
      sour: Number(m[4]) !== 0,
    };
  }
  return null;
}

/**
 * The flow-assurance panel of one two-fluid run, read off the LAST frame that
 * carries a `risk …` line (a case without `riskReport{}` prints none, so the
 * panel stays hidden).
 */
export function riskPanel(frames) {
  if (!Array.isArray(frames) || !frames.length) return null;
  const withRisk = [...frames].reverse().find(
    (f) => f && f.risk && (f.risk.erosion || f.risk.hydrate || f.risk.corrosion),
  );
  if (!withRisk) return null;
  const { erosion, hydrate, corrosion } = withRisk.risk;
  if (!erosion || !hydrate || !corrosion) return null;

  const lines = [
    {
      label: '冲蚀（API RP 14E）',
      value: `Ve ${H.fmt(erosion.velim, 4)} m/s · vm ${H.fmt(erosion.vm, 4)} m/s`,
      detail:
        `最坏段 #${erosion.seg}，ρ_m ${H.fmt(erosion.rhoM, 4)} kg/m³，` +
        `vm/Ve = ${H.fmt(erosion.ratio, 4)}，超限 ${erosion.over} 段`,
      flag: erosion.over > 0 ? 'warn' : 'ok',
    },
    {
      label: '水合物（无阻塞读数）',
      value: `亚冷却 ${H.fmt(hydrate.subcool, 4)} K · ${HYDRATE_GRADE[hydrate.grade] || '—'}`,
      detail:
        `段 #${hydrate.seg}：p ${H.fmt(hydrate.p, 4)} Pa，T ${H.fmt(hydrate.T, 4)} K，` +
        `Teq ${H.fmt(hydrate.teq, 4)} K。h 不进 d_eff、不增黏（§5.186）`,
      flag: hydrate.grade >= 2 ? 'warn' : hydrate.grade === 1 ? 'note' : 'ok',
    },
    {
      label: '腐蚀（分压读数）',
      value: `CO₂ ${H.fmt(corrosion.co2Partial, 4)} Pa · 带 ${CO2_BAND[corrosion.band] || '—'}` +
        (corrosion.sour ? ' · 含 H₂S' : ''),
      detail: `H₂S ${H.fmt(corrosion.h2sPartial, 4)} Pa（0.5 psia 酸线）`,
      flag: corrosion.sour || corrosion.band === 2 ? 'warn' : 'ok',
    },
  ];
  return {
    t: withRisk.t,
    erosion,
    hydrate,
    corrosion,
    lines,
    summary:
      `t=${withRisk.t} s：冲蚀比 ${H.fmt(erosion.ratio, 3)}（超限 ${erosion.over} 段）、` +
      `水合物 ${HYDRATE_GRADE[hydrate.grade] || '—'}、` +
      `CO₂ 带 ${CO2_BAND[corrosion.band] || '—'}${corrosion.sour ? ' + H₂S' : ''}`,
  };
}
