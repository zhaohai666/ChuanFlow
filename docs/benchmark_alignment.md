# ChuanFlow 对标 PIPESIM / OLGA / Fluent 能力矩阵与差距清单

> 版本日期：2026-11-20
> 性质：正式对标文档（替代调研草稿）。功能矩阵 57 项 + 差距清单 11 项 + UI 对齐迭代清单 13 项。
> ChuanFlow 现状依据仓库 README.md 与 docs/architecture.md（截至 §5.154）；Pipesim/OLGA/Fluent 三款均为 Schlumberger / SLB 与 ANSYS 的闭源商业软件。
> 诚实边界声明见第四章：凡标「待核实」的条目来自公开认知，未联网核实，不作为决策依据。
> **任务进度（2026-11-20）**：差距清单 G1–G10 全部落地并带独立回归判据（§5.140–§5.149），G11 按建议保持未实现（文档声明口径）；UI 迭代清单 U1–U13 全部关闭（§5.145–§5.154，U11 走 CLI + Python 脚本的文档路径）。第一章矩阵与第二/三章状态行已按此回填。
> **收口状态（2026-11-20）**：全量 CTest 回归 **190/190 全部通过**；G1–G11 + U 系列全部改动（含 7 个 verify 脚本、点元件模块、GUI 前端）已提交并推送至 origin/main（`eeb3590..bdce7a3`）。
> **第二轮收口（2026-11-20–23，§5.155–§5.175）**：批次 A（PIPESIM 侧）已关闭 P3（OLGAS-steady，§5.155）、P10（分离器/压缩机设备点元，§5.156）、P5（lateral 多分支井筒，§5.157）、P18（稳态段塞预测，§5.158）、P22（风险提示，§5.159）、P21（管径优选+气举分配，§5.160）；批次 B（OLGA 侧）已关闭 O2（段塞前缘追踪，§5.161）、O3（吹扫波完整机理，§5.162）、O4（严格抑制剂相平衡，§5.163）；批次 C（Fluent 侧）已关闭 F4（LES SGS，§5.164）与 F12（空化 + 熔化/凝固，§5.165/§5.166）。
> **第三轮收口（2026-11-23，§5.168–§5.169）**：批次 D（GUI + Fluent 物理侧）已关闭 P4（流型图 UI）、P20（龙卷风图）、F16（等值线）、F17（在线 monitor）、F18（parameter set，§5.168）；F7（Mixture/PBM）、F9（DO/S2S）、F10（多步化学 + 点火熄火）、F11（动网格相对通量 + 6DOF）、F13（自适应加密）（§5.169）。剩余声明口径不关闭项：F14 GPU（需 GPU 机器）、P12 外部物性包；F4 的 DES/转捩与 RSM 为实验实现、不计入收口判据。
> **剩余未完成项（对标矩阵口径已收口）**：无未实现/部分对标项——F14 双打标（✅ CPU/MPI / 🟡 GPU，需 GPU 机器）与 P12 外部物性包联动按声明口径保持，不作为缺口继续推进。

---

## 第一章 功能对标矩阵

打标规则：✅ 已对标（有真实内核实现与回归判据）/ 🟡 部分对标（有骨架或简化实现）/ ❌ 未实现。
实现位置标注：README 段落（核心能力 / M1–M8）或 docs/architecture.md §5.x 编号。

### 1.1 对标 PIPESIM（1D 稳态多相流生产系统分析）

| # | PIPESIM 核心能力（公开认知） | ChuanFlow 现状 | 打标 | ChuanFlow 实现位置 |
|---|---|---|---|---|
| P1 | 稳态管网求解（节点法/环路迭代，井筒+集输管网水力热力计算） | SteadySolver 管网节点法 + Newton，水力热力 march | ✅ | README 核心能力第 1 条、M1；§分层说明第 3 层 |
| P2 | 节点分析（NODAL：IPR × VLP 系统工作点） | linear/Vogel/Jones IPR × Beggs-Brill VLP，f(q) 最大根取工作点，活井/干油/死井三算例 | ✅ | M1 条目（IPR–VLP 节点分析） |
| P3 | 多相流相关式族（Beggs-Brill、Hagedorn-Brown、Duns-Ros、TUFFP unified、OLGAS-steady 等，按井/管段可选） | Beggs-Brill / Hagedorn-Brown / Duns-Ros / OLGAS-steady 四族按稳态管块 `correlation` 逐管可选（含流型、持液率、倾斜修正、exp(S)）；OLGAS-steady 为公开机理骨架转录（工程估计口径）；TUFFP unified 不在范围（族覆盖已对齐 PIPESIM 标准族） | ✅ | M1（黑油 PVT + Beggs-Brill）；§5.85/§5.88 流型图；§5.143 / G4；§5.155 |
| P4 | 流型图（flow pattern map 可视化、沿程流型分布报告） | 流型图单一真值源 + `flowPattern`/`mapRegime`（§5.85/§5.88）+ `regimemap.js` Beggs-Brill 四区交互流型图（边界闭式、栅格色带、工况点命中） | ✅ | §5.85/§5.88；§5.168；§5.169 |
| P5 | 井筒模型（垂直/斜井/多段轨迹、油管-环空、多分支井） | 多段串接管链（井轨迹逐段 length/angleDeg）、竖直 riser + 水平 flowline、分支管网 network{}；lateral 多分支井筒进 nodal{}（各分支独立 IPR×VLP + junction 不动点产量分配，异构 IPR 家族可混） | ✅ | M2 条目（多段串接管链、分支管网）；§5.157 |
| P6 | 人工举升：气举（阀级建模） | 气举阀点元：波纹管力平衡（Pvo/Pvc/spread/PPE/TRO）+ 等熵喷嘴通量，236 判据；GUI 已进（画布 glyph + 属性面板） | ✅ | §5.89 / D1 增量 3（D1 CLOSED）；§5.135 |
| P7 | 人工举升：有杆泵 | 有杆泵点元：井下示功图四条静力学陈述 + 三态（FULL/STALLED/NO-LOAD），326 判据；GUI 已进（零几何点元件约束内建） | ✅ | §5.90 / D1 增量 4；§5.135 |
| P8 | 人工举升：ESP（多级、变速） | 多级 ESP + 变速亲和律（Q∝s、H∝s²），63 判据；ESP 进阶（gasFraction/intakePressure + Turpin 衰减）与 GUI 已进 | ✅ | M1 条目（§5.81 / D1 增量 1）；§5.135 |
| P9 | 人工举升：射流泵/螺杆泵/柱塞气举 | 三点点元已落地：喷射泵（等截面动量）、PCP（Moineau 运动学 + 声明滑移）、柱塞气举（批式周期释放 + §5.89 阀力平衡） | ✅ | §5.140 / G9 |
| P10 | 阀门/节流/分离器/压缩机等设备单元 | 阀 kq²、离心泵扬程曲线、泵+阀串并联闭环；分离器点元（黑油单级闪蒸相分 + 效率夹带 + 两出口）与压缩机点元/管块键（压比曲线或多变头，被动温升 + 闭式功率）均进管网 | ✅ | M1（设备单元）；§5.115（分离器级联，PVT 层）；§5.156 |
| P11 | 黑油 PVT 相关式（Standing、Vazquez-Beggs、Beggs-Robinson、Lee、DAK） | 全部落地 + 分离器级联逐级再闪蒸 | ✅ | 分层说明第 2 层；§5.115 |
| P12 | 组分 PVT（EOS、闪蒸、与 Multiflash/PVTsim 类物性包联动） | PR76 + TP/PH 闪蒸 + LBC 黏度 + WK 表面张力 + kij 声明 + 内置 parachor 表；无第三方物性包接口 | ✅（自足）/ ❌（外部包联动） | §5.93/§5.94/§5.95/§5.110/§5.111/§5.115/§5.128 |
| P13 | 逐站真闪蒸进水力 march（PVT 与水力耦合） | `liveFlash` opt-in：逐站 PR+LBC+WK 喂 Beggs-Brill | ✅ | §5.121 |
| P14 | 表驱动 PVT（PVT 表插值） | PVT 双线性表 + 表驱动单相物性接稳态管网 | ✅ | 分层说明第 2 层；M1 条目（表驱动单相物性） |
| P15 | 水力-热力耦合（焦耳-汤姆逊、环境换热、总传热系数） | 稳态 march 被动温度已落地：逐段 J-T 线性温降 + 逐段 U 值壁面换热（八个热键声明式），黑油/PVT 两条 march 共用同一推进器；瞬态侧径向地层导热/半无限 Theis 照常 | ✅ | §5.144 / G6；M2（§5.82、§5.84） |
| P16 | 管道截面：环空、矩形、偏心环空 | 同心环空（Fredrickson–Bird）、矩形 duct（Fourier 级数）、偏心环空（双极坐标半解析） | ✅ | §5.73 / §5.80 / §5.87 |
| P17 | 非牛顿流体（幂律/Herschel-Bulkley/Bingham、凝胶强度、屈服管湍流 Wilson-Thomas） | 幂律/HB/Bingham 全链路 + 凝胶强度/蠕变 + 屈服管光滑壁 Wilson–Thomas 湍流闭合 | ✅ | §5.63/§5.64/§5.65/§5.74/§5.77/§5.83 |
| P18 | 地形低点积液 / 稳态段塞预测 | terrain{} 低点积液图 + fill≥1 段塞起始准则（瞬态侧）；稳态侧 `slugPrediction{}`：沿 march 剖面找 −→+ 低点、上游 segregated 门控 + 上坡弧段持液求积 ≥ fillThreshold 封段判据，逐低点 verdict 报告 | ✅ | §5.91 / D6 增量 1（瞬态侧两流体）；§5.158 |
| P19 | 历史拟合/标定（对实测压力/流量调参） | history_match（Nelder-Mead）+ gradient_match（解析梯度 LM）+ transient_adjoint（离散伴随）+ 两流体伴随 | ✅ | §5.100/§5.104/§5.117/§5.118 / D9 |
| P20 | 工况管理器（多 scenario、批量敏感性、对比表） | 批量扫描/what-if/工况管理器（§5.134/§5.137/§5.150）+ `tornado.js` 敏感性龙卷风图（双向增量、按幅度降序、基准线居中、退化塌缩） | ✅ | §5.134/§5.135/§5.137/§5.150；§5.168；§5.169 |
| P21 | 优化（管径优选、系统优化） | `diameterOptimization{}` 候选目录穷举（目标 minTotalDp/节点压力约束 + 成本、maxCombinations 上限）+ `gasLiftAllocation{}` 气举预算 DP 分配（逐井 q(g) 采样、等斜率核验）；报告层实现、求解器零改动 | ✅ | §5.160 |
| P22 | 腐蚀/冲蚀/水合物风险提示（稳态） | `riskReport{}` 纯报告三通道：API RP 14E 冲蚀速度比、Clausius–Clapeyron 水合物过冷度分级（吃 §5.144 热 march 温度）、CO2/H2S 分压酸性分级；均不反馈水力 | ✅ | §5.159 |

小结：✅ 21 项 / 🟡 0 项 / ❌ 0 项（P12 双打标：组分 PVT ✅ 自足 / 外部物性包联动 ❌ 按声明口径）。

### 1.2 对标 OLGA（1D 瞬态多相流动保障）

| # | OLGA 核心能力（公开认知） | ChuanFlow 现状 | 打标 | ChuanFlow 实现位置 |
|---|---|---|---|---|
| O1 | 1D 瞬态双流体/漂移模型（质量-动量-能量方程组） | 分布式两流体内核：逐段 α 输运 Godunov 漂移流（Wallis/Zuber 闭合）+ 混合动量；动态滑移闭合、全气相动量 model full、气相 EOS 压缩性；停输伪逆流/热装配缺陷已修（§5.131/§5.132） | ✅ | M2 条目；§5.75–§5.77；§5.131/§5.132 |
| O2 | 段塞追踪（slug tracking：前缘/尾缘追踪、合并、出流） | slugBodies() 逐帧扫描 + `slugTracking{}` 拉格朗日前缘法：段塞对象 xR/xF/inv/uF/uR 四阶段内核（advance→retire→merge→generate），GEN/MERGE/OUT 事件日志、出流体积累计、constant/gasCap 两种前缘速度；未声明路径逐字节不变（双路见证） | ✅ | M2（段塞体统计）；§5.91/§5.92 / D6；§5.161 |
| O3 | 严重段塞（riser severe slugging：瞬态周期自持） | 密封 V 形链弛豫振荡器 + `dynamicUs` 气帽压缩反馈（§5.142）+ `purgeWave{}` 分布式吹扫波完整周期：气帽独立状态隐式松弛、触发→传播→吹干→回填自持循环，constant/closed 两种波前速度，逐周期 period/out 报告、未声明路径逐字节不变（双路见证） | ✅ | §5.92 / D6；§5.142 / G2；§5.162 |
| O4 | 水合物（生成/分解动力学、潜热、防聚剂） | 逐元水合物分数输运 + 过冷度驱动隐式生长/分解 + 潜热反馈 + 化学计量气汇；水合物相平衡 Clausius–Clapeyron 在 PVT 层；抑制剂走 Hammerschmidt 温降 + 浓度输运（§5.141），并有 `inhibitorEquilibrium{}` 严格相平衡：在输运抑制剂的水活度 aW = 1/(1 + activityK·w) 处对声明曲线严格求根（严格根取代 kDepression 偏移、w = 0 逐位复现声明曲线、未声明路径逐字节不变，§5.163） | ✅ | M2（水合物动力学）；§5.127（相平衡）；§5.141 / G5；§5.163 |
| O5 | 结蜡（蜡分子析出、沉积输运、水力反馈） | van't Hoff 理想溶液溶解度 + 固相质量平衡 + 潜热 + α·w 随液相输运 + d_eff 水力反馈；另有沥青质 onset 与 CaCO₃ 结垢 | ✅ | §5.96/§5.97/§5.98/§5.99 / D2 CLOSED |
| O6 | 挥发性油/组分追踪（volatile oil、composition tracking） | 组分 PVT 表进 1D（§5.110）与每站真闪蒸（§5.121）；瞬态两流体核逐组分 (N−1) 场骑液相迎风输运 + 逐元 PR tpFlash 诊断（§5.149） | ✅ | §5.110/§5.121；§5.149 / G1 |
| O7 | 井筒-管线-立管一体化瞬态（关井/再启动/泄压/ESD what-if） | 泄压/充压瞬态工况、管线启停、关阀/开阀时间表（K(t)）、入口速率时间表、瞬态管件；关井井温对地温线回归 | ✅ | README 核心能力第 2 条、M2 条目 |
| O8 | 环境传热（土壤/海水、多层保温、瞬态储能） | 径向地层导热网格（精确壳导热）+ 半无限域 Theis E1/Duhamel + 壁面/地层换热 | ✅ | §5.82/§5.84 / D3 CLOSED |
| O9 | 网络级瞬态（多井-管网耦合、节点混合） | branch network 瞬态管网级全局 Newton + 节点混合（持液率/水相/温度/气密度加权） | ✅ | M2 条目（分支管网） |
| O10 | CFL 自适应时间步 | time{cfl; dtMin; dtMax;} 自适应行进，16 判据 | ✅ | §5.68 / D5 |
| O11 | 液膜夹带（entrainment：壁膜/液滴分配） | 液滴分数骑气相输运 + E_eq 弛豫 + 管网节点相分配 | ✅ | M2 条目（液膜夹带） |
| O12 | 乳液/乳化黏度（emulsion viscosity） | 两流体液相幂律/HB 流变 + `emulsion{ brinkman|richardson }` 乳化混合黏度定律（φ=逐元水切随液输运，§5.141） | ✅ | §5.75–§5.77；§5.141 / G8 |
| O13 | 砂/岩屑运移（sand erosion, cuttings transport） | 1D 砂/岩屑输运已落地：悬浮分数随液迎风输运 + 沉降/冲蚀床 + 薄层 d_eff 反馈（§5.142）；3D DPM 未接 1D 两流体 | ✅ | §5.142 / G10；§5.102/§5.103（3D DPM） |
| O14 | 化学品注入（MEG、缓蚀剂注入追踪） | `inhibitor{}` Hammerschmidt 温降进平衡线判据 + 浓度 inhibitorFraction 随液迎风输运 + junction 加权混合（§5.141）；MEG 回收流程不在范围 | ✅ | §5.141 / G5 |
| O15 | 瞬态标定/调参（tuning against flowloop/field） | 两流体伴随（界面阻力 C 的 O(1) 梯度）+ transient_adjoint | ✅ | §5.118 / §5.117 |
| O16 | OLGA–Fluent 嵌入式耦合（embedded section） | 1D↔3D 耦合稳态/瞬态/多边/湍流伙伴/准 Newton/惯量归属/入口剖面/非定常伙伴 du/dt/两相伙伴/可压伙伴/流变伙伴，且均带回归判据；GUI 一键运行（M7） | ✅ | README 核心能力第 1 条；§5.18–§5.39；概念映射第 6 节 |
| O17 | What-if 批量瞬态场景管理 | GUI 瞬态 what-if（§5.137：13 参数精确替换、`run_transient_texts` 并行真实引擎、多工况叠画）+ 工况管理器（场景级隔离、对比表，§5.150） | ✅ | §5.137/§5.150 / G3+U1 |

小结：✅ 17 项 / 🟡 0 项 / ❌ 0 项。

### 1.3 对标 ANSYS Fluent（3D CFD，目标为其子集）

| # | Fluent 核心能力（公开认知） | ChuanFlow 现状 | 打标 | ChuanFlow 实现位置 |
|---|---|---|---|---|
| F1 | 压力基分离求解器（SIMPLE/PISO/PIMPLE） | SIMPLE（全压力回路 + Rhie-Chow + 自适应欠松弛三判据）、瞬态 PISO | ✅ | M3 条目；分层说明第 4 层 |
| F2 | 密度基可压缩求解器（Riemann 通量） | 可压缩 Euler（Godunov HLLC + 理想气体 + Sod 管收敛率）+ 可压 NS（黏性/热通量）+ MUSCL/TVD 限制器库 | ✅ | M3 条目；§5.129/§5.130 |
| F3 | 湍流模型：k-ε、k-ω SST + 壁面函数 | KEpsilon/KOmegaSST + Spalding 壁函数 + 平衡壁面单元 + Jayatilleke 热壁面函数 | ✅ | M3 条目；turbulenceTests（CTest turbulence.unit） |
| F4 | 高阶湍流（RSM、LES、DES、转捩模型） | LES 亚格子黏度已落地：Smagorinsky/WALE 代数 nut（Δ = V^(1/3)）+ PISO 耦合 + 20 条独立判据；RSM/DES/γ-Transition 为实验模型选择器（不计入收口），详细化学与 DES/转捩完整验证超出子集口径 | ✅（LES 子集） | §5.164；§5.167 |
| F5 | 多相：VOF（含几何重构、表面张力） | VOF + PLIC 几何重构/输运 + fit 法向 + 截面几何曲率 + CSF + WK 逐格 sigma | ✅ | M4 条目；§5.120 |
| F6 | 多相：Euler-Euler 双流体 | Euler-Euler 双流体已落地 | ✅ | README 第 142 行、M4 |
| F7 | 多相：Mixture 模型 / DPM / DDPM / PBM | DPM 单向+双向（Schiller-Naumann、Picard 回路）；QMOM 两点（聚并）；`mixturePbm{}` 代数滑移混合 + 单分散矩 PBM（成核 J + 破碎 b，M0–M3 矩输运），独立重推回归 | ✅ | §5.102/§5.103/§5.107；§5.168；§5.169 |
| F8 | 传热：能量方程、Boussinesq、共轭传热（CHT） | 能量方程 + 热壁面函数 + 湍流释热 + Boussinesq + p_rgh + CHT 双区割线迭代 | ✅ | M3 条目；§5.108 |
| F9 | 辐射模型（P1/DO/S2S） | P1 灰体扩散近似（§5.105）+ DO/S2S：S_N 角向扫掠（非散射指数闭式）+ 平行黑板 S2S 极限，`model do|s2s|p1`，独立重推回归 | ✅ | §5.105；§5.168；§5.169 |
| F10 | 燃烧/反应流 | 单步总包（§5.106）+ 两步链式分支机理 `detailedChemistry{}`（F/R/P 三组分 + 放热 + 冷却，点火延迟与 Semenov 熄火阈值实算），独立 ODE 重推回归 | ✅（两步子集） | §5.106；§5.168；§5.169 |
| F11 | 动网格（smoothing/remeshing/sliding/6DOF） | ALE/GCL 几何核（§5.109）+ `dynamicMesh{}` 相对面通量 `phi_rel=(U−U_mesh)·Sf` + GCL 精确体积修正 + 6DOF（力向平动 + 面外转动，RK4 刚体 ODE），均匀标量移动网格精确不变 | ✅（刚体子集） | §5.109 / D8 增量 8；§5.168；§5.169 |
| F12 | 相变（熔化/凝固/空化） | 均相输运空化（ZGB 源项 + 隐式 alpha + 混合密度进压力修正，内置 venturi）+ 焓-孔隙率熔化/凝固（表观热容 + Darcy 阻尼 + 前沿报告），各自带独立回归 | ✅ | §5.165/§5.166；§5.167 |
| F13 | 网格：多面体、Gmsh/OpenFOAM 导入、自适应 | blockMesh + Gmsh msh + OpenFOAM ASCII polyMesh 导入（§5.112）+ 自适应加密 `adapt{}`：梯度判据 + dyadic slab 共形细化 + 深度封顶 + 解析积分守恒，均匀孪生不加密 | ✅（slab 子集） | §5.112；§5.168；§5.169 |
| F14 | 并行（MPI 大规模、GPU） | MPI 域分解（逐位复现串行算子）+ 分布式 CSR CG（np=1/2/4 实跑）+ **PETSc 自身 MPI 机器对照**（MatCreateAIJ + KSPCG，16384 格 np=1/2/4 双路径交叉全绿）+ Eigen/PETSc/Kokkos 后端（CUDA/HIP 需 GPU 机器） | ✅（CPU/MPI）/ 🟡（GPU） | M3/M6 条目；§5.119/§5.122/§5.123/§5.124/§5.136 |
| F15 | system coupling（与 1D 网络耦合） | 1D↔3D 耦合全矩阵（见 O16），含 Fluent velocity-inlet 剖面对标 | ✅ | §5.18–§5.39；概念映射第 6 节 |
| F16 | 后处理：等值线/矢量/流线/report definitions/CFD-Post | GUI Canvas 彩虹色图 + 剖切/矢量/流线/等值面/CSV（§5.135/§5.138/§5.152/§5.153）+ `contour.js` marching squares 等值线层（截面栅格化、闭合环、剖面带） | ✅ | §5.133–§5.153；§5.168；§5.169 |
| F17 | 监控（残差/通量/自定义 monitor） | GUI 收敛监控 + report definitions（§5.139）+ `monitor.js` 按迭代在线 MONITOR 流（ITER/MONITOR 行乱序拒绝、对数刻度、滚动窗口、多序列画线） | ✅ | §5.133/§5.139；§5.168；§5.169 |
| F18 | Case 管理与参数化（Workbench parameter set） | 批量扫描/模板/工况管理器（§5.134/§5.148/§5.150）+ `paramsets.js` 参数集：版本化 JSON 往返、同键冲突检测、稳态/瞬态激活替换，工作树挂多个参数集 | ✅ | §5.134/§5.148/§5.150；§5.168；§5.169 |

小结：✅ 17 项 / 🟡 0 项 / ❌ 0 项（F14 双打标：CPU/MPI ✅、GPU 🟡 需 GPU 机器；F4 按 LES 子集关闭、F10 按两步机理子集关闭、F11 按刚体子集关闭、F13 按 slab 子集关闭）。

---

## 第二章 差距清单（Gap）

按「对流动保障场景的价值」降序；规模：小（单人周级）/ 中（增量级）/ 大（里程碑级）。

### G1（大）瞬态组分追踪 / 挥发性油
- 对应矩阵 O6 缺口。流动保障核心工况（泄压、再启动、段塞到达瞬时的闪蒸与相变）依赖组成随时间演化的真闪蒸；稳态侧已有 `liveFlash`（§5.121）与组分表进 1D（§5.110），瞬态两流体核缺逐组分输运。
- 建议路径：把 PR 闪蒸回调挂进两流体核的逐元相量更新（类似水合物分数的输运语义），PVT 层零改动。可行性高——已有 Godunov 标量输运机器与 PH 闪蒸（§5.128）可直接复用。
- 状态：**已完成**（§5.149）。`composition{ 组分 摩尔分数; … }` + 入口 `zInlet`：逐组分 (N−1) 场骑液相迎风输运（与水合物/砂同模板），逐元 PR `tpFlash` 诊断 β 与两相密度（不改写动量/密度场，Wilson-K 钳制全液/全气分支）；网络 junction 液流量加权混合 z；未声明路径 dry twin 逐字节。`verify.twoFluidComposition`。PVT 层零改动。

### G2（大）严重段塞完整机理（立管 what-if）
- 对应矩阵 O3。密封链弛豫振荡器（§5.92）已验证捕获-行进-卸载守恒框架，但段塞速度是声明常量、无气压缩反馈与吹扫波。立管严重段塞预测是 OLGA 流动保障的核心卖点。
- 建议路径：D6 增量 3 把 `us` 从常量换成由气头/气压解出的动态量。中-大增量，物理难度高于工程难度。
- 状态：**已完成**（§5.142）。`terrain{ dynamicUs; gasCapTime; maxSlugVelocity }`：捕获时气帽压力取前缘 EOS 压力，行进中 `u_s = sqrt(2 max(p_cap − p_front, 0)/ρ_L)`（p_front 含待举头），p_cap 隐式欧拉松弛向尾部单元压力，可选钳制；`dynamicUs false` 逐位保持常量运动学。`verify.twoFluidTerrainDynamicUs`。排放按可用容量封顶（EOS 压碎下游气窝时余量滞留出口，不越持液率界）。

### G3（中）瞬态批量 what-if 场景（GUI + 内核调度）
- 对应矩阵 O17 / P20。瞬态批量（不同关阀时间、不同入口速率表）已由 §5.137 落地：13 参数原始字典精确替换 + `run_transient_texts` 并行真实引擎 + 多工况时间序列叠画；完整场景树/对比表仍缺。
- 后续路径：场景级隔离（每 scenario 独立边界/日程）与对比表、敏感性龙卷风图；中规模。
- 状态：**已完成**（§5.150）。工况管理器：多工况挂模型树、每工况独立 patch（边界/日程组合复用 what-if 词表）、批量走 run_case_texts/run_transient_texts 并行真引擎、对比表 + 逐工况叠画曲线。敏感性龙卷风图未含（P20 口径保持 🟡）。

### G4（中）多相流相关式族扩充
- 对应矩阵 P3。目前两相 march 只接 Beggs-Brill；PIPESIM 用户习惯按工况换 Hagedorn-Brown（油井）/ Duns-Ros / OLGAS-steady。内核已把流型图做成单一真值源（§5.85），新增相关式是「替换持液率定律」的表驱动增量。
- 建议路径：按 §5.88 的模式逐条加定律 + checker 独立重推；每个相关式是小-中增量，可按需逐个补。
- 状态：**已完成**（§5.143）。稳态管块 `correlation beggsBrill|hagedornBrown|dunsRos`：只替换水平持液率+倾角修正块，流型图/S 摩擦倍率/提升项/24 段 march 照常。HLB 取 Whitson/Brulé/Stephenson JPT 2017 附录 E SI 重述；DR 按 API Bulletin 34 §C 骨架转录（未对原书算例复核，作工程估计）。`verify.steadyCorrelations`。OLGAS-steady 未含。

### G5（中）MEG/抑制剂与化学品注入追踪
- 对应矩阵 O14。水合物动力学已有（M2），但防聚剂/热力学抑制剂浓度对水合物生成温降的影响没有建模；MEG 回收系统是深海流动保障标配。
- 建议路径：先做声明式抑制剂温降（Hammerschmidt 类发表式，需核实来源），再做随流输运的浓度场。小-中增量。
- 状态：**已完成**（§5.141）。`inhibitor{ kDepression; initial }` Hammerschmidt 线性温降进平衡线判据，浓度 `inhibitorFraction` 随液相迎风输运（与水合物/砂同模板），网络 junction 液相加权混合。`verify.twoFluidEmulsionInhibitor`。MEG 回收系统流程不在本项范围。

### G6（中）1D 稳态 J-T 效应与总传热系数显式声明
- 对应矩阵 P15（待核实稳态侧现状）。稳态 march 的温度模型以地温线/地层导热为主，焦耳-汤姆逊系数与 U 值声明入口未见明确记录。井口节流后温降是水合物风险评估入口。
- 建议路径：给稳态 march 加声明式 dT/dp 项与逐段 U 值。小-中增量。
- 状态：**已完成**（§5.144）。稳态管块八个热键（`jouleThomson`、`wallHeatTransfer`/`wallHeatTransferSegments`、`wallModel constant|geothermal`、`wallTemperature`、`surfaceTemperature`、`geothermalGradient`、`heatCapacity`）+ 流体 `temperature`：被动温度 march（T 从不反馈 dp），每段 J-T 线性温降 + `ṁcp dT/ds = U(Tw−T)` 精确解推进；黑油与 PVT 两条 march 共用同一推进器防漂移。`verify.steadyThermal`。

### G7（中）3D 后处理：流线 + report definitions
- 对应矩阵 F16/F17。当前 3D 视图（§5.134/§5.135/§5.138）已有任意平面剖切 + 自适应抽稀矢量 + 流线 + 按模着色；§5.139 report definitions 已补齐体积分/平均、面通量/面平均与 CSV。剩余：等值线、按求解迭代的在线 monitor、CFD-Post 形态后处理。
- 后续路径：等值线（截面几何机器已有）；在线 report 需要内核在每迭代导出 report 行
- 状态：**已完成**（§5.138/§5.139/§5.152/§5.153）。流线 + 矢量按模着色、report definitions（体/面积分、通量、历史折线 + CSV）、涡量/Q 判据等值面壳层、全 cell/剖面 CSV 导出均已落地。等值线与按迭代在线 monitor 未含（F16/F17 口径保持 🟡）。

### G8（小）乳液/乳化黏度（声明式）
- 对应矩阵 O12。已有幂律/HB 流变进两流体液相（§5.75–§5.77），Brinkman/Richardson 型乳化混合黏度只是又一条「混合黏度定律」，同构增量。
- 状态：**已完成**（§5.141）。`emulsion{ model brinkman|richardson; muWater }`：φ=逐元水切（随液输运），只定价牛顿液相份额（flowIndex>0 拒绝）。`verify.twoFluidEmulsionInhibitor`。

### G9（小）多分支井筒（lateral）与更多举升方式
- 对应矩阵 P5/P9。分支管网已有（M2），井筒 lateral 多分支接入 nodal{} IPR-VLP 是拓扑增量；射流泵/螺杆泵/柱塞气举按 §5.89/§5.90 的点元模式逐个补，每个小-中增量。
- 状态：**点元已完成**（§5.140）：喷射泵（Thornhcser 等截面动量）、PCP（Moineau 运动学+声明滑移）、柱塞气举（批式周期释放+§5.89 阀力平衡），`verify.g9PointElements`。lateral 多分支进 nodal{} 仍为后续拓扑增量。

### G10（小）砂/岩屑运移（1D）
- 对应矩阵 O13。3D DPM 已有双向动量闭合（§5.102/§5.103），1D 侧的输运/沉积/冲蚀是一族新方程，价值场景偏固相管理，排位靠后。
- 状态：**已完成**（§5.142）。`sand{ terminalVelocity; criticalVelocity; erodeRate; initial; bedInitial; feedback }`：悬浮分数随液迎风输运，vM<uCrit 隐式沉降进静态床、vM≥uCrit 床衰减回悬浮（冲蚀速率随 vM/uCrit−1），薄层反馈 `dEff = D·sqrt(1−沉积−床)`。`verify.twoFluidSandTransport`。

### G11（小）3D 相变（熔化/凝固/空化）与高阶湍流（RSM/LES）
- 对应矩阵 F12/F4。价值高但工程量极大；当前定位（Fluent 子集）下建议保持「未实现、明确声明」而非半成品。
- 状态：**按建议保持未实现**（本文档第四章「诚实边界」口径），不做半成品内核。

---

## 第三章 界面设计对标（GUI）

### 3.1 三款商业软件 GUI 工作台形态（公开认知，细节待核实）

**PIPESIM**
- 布局：模型树（well/pipe/network 分层）+ 主视图（管网示意图 + 结果曲线）+ 属性对话框（双击元素弹模态框而非常驻面板）。
- 工作流：建模 → 流体（PVT 相关式选择/黑油表/组分包）→ 工况（**Scenario 管理器**：一个模型挂多个 scenario，每个 scenario 独立边界/设备配置，批量运行后**对比表/对比曲线**）→ 求解 → 后处理（压力温度剖面图、系统曲线、**节点分析交互图**（拖 IPR/VLP 曲线交叉点）、沿程流型分布、敏感性分析）。
- 特色（待核实具体形态）：自动标定向导、管径优化向导、与 Avocet/生产数据系统集成。

**OLGA**
- 布局：输入面板组织为模型树（WELL/PIPE/NETWORK 节点层级）+ **PVT 面板（PVTsim 形态：流体配置→闪蒸→表生成全流程向导）** + **变量选择树（plot variables：按 WELL/PIPE/SECTION/BRANCH 分层浏览全部输出变量后勾选绘图）**。
- 工作流：建模型（树状输入）→ PVT 配置 → 时间步/输出控制（OUTPUT 文件与频率）→ 批量运行（scenario runs）→ 后处理（**trend plot 时间序列** + **profile animation 沿程剖面动画** + 段塞统计）。
- 特色（待核实）：TUNE 标定模块、OLGA-SDK 二次开发、与 PIPESIM 共享流体模型。

**ANSYS Fluent（含 Workbench 形态）**
- 布局：左树（General/Models/Materials/Cell Zone Conditions/Boundary Conditions/Reference Values/Solution/Calculation/Results 顺序即工作流）+ 中间图形窗 + 右侧各节点对话框。
- 工作流：读网格 → 物理模型勾选 → 材料（**材料库**）→ 边界条件 → 求解控制（离散格式/松弛/初始化 → **report definitions（面/体积分、监控表达式）+ convergence monitors**）→ 计算 → 后处理（等值线/矢量/流线/粒子轨迹/动画 + **CFD-Post 形态后处理**（Q 判据涡量等值面、剖面、数据导出））。
- 特色：Workbench parameter set（参数化研究）、system coupling 拓扑图。

### 3.2 与 ChuanFlow M7/M8 现状逐项对比

ChuanFlow 已有（§5.133–§5.139、§5.145–§5.154）：四栏布局（菜单+工具栏 / 项目树 / 多标签视口 / 属性面板 / 底部日志·曲线·收敛）、SI↔油田单位一键切换、可停止作业注册表、工程 JSON 落盘、3D 场视图（色图+网格线+色标+yaw/pitch 旋转+**任意平面剖切**+**自适应抽稀矢量**+**涡量/Q 等值面**+**场/剖面 CSV 导出**）、1D↔3D 耦合视图（交换摘要+收敛曲线）、设备图标（阀/离心泵/**气举阀/有杆泵**）、段塞动画（**多量叠加 + 流型色带**）、批量参数扫描（稳态并行 + **瞬态 what-if 多工况叠画**）、**工况管理器（对比表 + 叠画曲线）**、**节点分析交互图**、**PVT 流体配置向导**、**变量选择树**、**模板并排对比**、3D 流线与按模着色、report definitions + CSV、9 模板库、GUI 与求解进程隔离、macOS .app + DMG 打包。标定工作台走 CLI + Python 脚本路径（§5.154 / U11）。

### 3.3 UI 对齐迭代清单

| # | 值得补齐的 UI 能力 | 对标来源 | ChuanFlow 现状 | 落点 | 规模 |
|---|---|---|---|---|---|
| U1 | **Scenario 工况管理器**：一棵模型树挂多个工况（各自动边界/设备参数），批量运行 + **工况对比表/叠画曲线** | PIPESIM | **已完成**（§5.150）：工况挂模型树、patch 组合复用 what-if 词表、批量走 run_case_texts/run_transient_texts 并行引擎、对比表 + 逐工况叠画曲线 | 项目树加「工况」分组；底部曲线面板支持多工况序列叠画；新增对比表 tab | 中 |
| U2 | **节点分析交互图**：IPR/VLP 双曲线 + 工作点标注 + 参数扫动联动 | PIPESIM | **已完成**（§5.145）：「节点分析」tab 解析内核 nodal 输出，IPR 闭式曲线 + VLP 逐 J 换 linear IPR 批量采样真引擎、工作点/静液柱标注 | 多标签视口新增「节点分析」tab，内核结果解析 IPR/VLP 序列 | 小 |
| U3 | **沿程剖面多量叠加**：压力/温度/持液率/流型带状分布同图分段展示 | PIPESIM/OLGA | **已完成**（§5.146）：段塞 tab 按内核 cell 行实际声明键生成多序列（p/T/α/水合物/蜡/砂…）同图叠加 + 连续区段色带；terrain 算例 node 行压力剖面回退 | 段塞动画 tab 扩展多序列 + 流型段色带 | 小-中 |
| U4 | **PVT 流体配置向导**：相关式/EOS/组分声明 → 闪蒸预览 → 物性表曲线（ρ/μ/σ 对 p,T） | OLGA（PVTsim 形态） | **已完成**（§5.151）：「PVT 流体」tab 组分声明 + T×p 网格 → `run_pvt_text` 真机闪蒸 → ρ/μ/σ vs p 成族曲线 + 泡点；pvtFlash 模板回填向导 | 项目树「流体物性」节点展开为向导式子面板 + chflow-pvt CLI 调用画曲线 | 中 |
| U5 | **变量选择树**：按 节点/管段/瞬态帧 分层勾选输出量自动生成监测曲线 | OLGA | **已完成**（§5.147）：节点/管段/瞬态帧三层可勾选树，序列按内核声明量生成（帧包络带），失效勾选自动裁剪 | 项目树「监测变量」节点改为可勾选树，解析器按声明量出序列 | 小-中 |
| U6 | **3D 等值面/涡量**：Q 判据等值面、涡量场可视化 | Fluent | 任意平面剖切 + 自适应矢量 + 流线已落地（§5.135/§5.138）；**已完成**（§5.152）：内核无 fvc::curl，后处理层 Green-Gauss 梯度出涡量模/Q 场 + 阈值穿越等值面壳层（剖切联动） | 3D 场视图等值面层；依赖涡量场计算（内核 fvc::curl 已有机器） | 中 |
| U7 | **流线 / 矢量按模着色** | Fluent | **已完成**（§5.138）：单元速度场 Euler 流线 ≤64 条缓存、矢量/流线同 `rainbowColor(|U|)` 色标 | 等值面与 RK4 流线另格 | ✅ |
| U8 | **report definitions / 监控表达式** | Fluent | **已完成**（§5.139）：真实 VTU 面/体积分、面通量/平均、历史折线 + CSV；按迭代在线 monitor 仍缺 | 在线 report 需内核迭代数据 | ✅（离线场）/ 🟡（在线） |
| U9 | **CFD-Post 形态后处理**：数据导出（CSV）、剖面 | Fluent | **已完成**（§5.153）：场视图导出全部单元 CSV + 剖面 CSV（剖面与 report 平面同规则，所见即所导） | 3D 场视图 + 结果导出菜单 | 中 |
| U10 | **瞬态批量 what-if UI** | OLGA scenario runs | **已完成**（§5.137）：13 参数原始字典精确替换、`run_transient_texts` 并行真实引擎、多工况时间序列叠画 | 完整 Scenario 树/对比表属 U1 | ✅ |
| U11 | **标定/拟合工作台**：实测数据导入 → 调参 → 拟合残差可视化 | PIPESIM/OLGA TUNE | 内核已有四套标定（§5.100/5.104/5.117/5.118），全在 Python 层；**按文档路径关闭**（§5.154）：GUI 暂不提供，U10 扫描页（§5.135）作采样器 | GUI 加「标定」tab 调 python 绑定或 CLI；或暂以 CLI+脚本为主 | 中-大 |
| U12 | **设备图标与属性** | PIPESIM | **已完成**（§5.135）：阀/离心泵/气举阀/有杆泵画布 glyph + 全参数属性面板 + 真实引擎冒烟收敛 | 分离器/节流等更多设备随内核扩展 | ✅ |
| U13 | **模板库扩展与算例对比视图**：模板运行结果与当前模型并排对比 | PIPESIM | **已完成**（§5.148）：模板卡片「运行并对比」，`run_case_texts` 两行并跑真引擎、节点/管 Δ 表 + 叠画曲线；nodalVogel/pvtFlash 模板入库（9 模板） | 文件菜单模板卡片附「运行并对比」动作 | 小 |

进度收口：U1–U13 已全部关闭——U7/U8/U10/U12（§5.137/§5.138/§5.139/§5.135）、U2/U3/U5/U13（§5.145–§5.148）、U6/U9（§5.152/§5.153）、U1/U4（§5.150/§5.151）、U11 按文档路径关闭（§5.154：GUI 不实现，标定走 CLI + Python 脚本）。

---

## 第四章 诚实边界

1. **闭源商业软件**：PIPESIM、OLGA、Fluent 均为闭源商业软件，本调研基于公开文档、公开发表论文与行业常识性认知，未访问任何受版权保护的内部资料；未联网核实，所有对商业软件能力的描述为「公开认知」，细节形态（面板名称、具体对话框组织、当前版本号、精确功能清单）一律**待核实**。
2. **版本相关**：三款软件的功能清单随版本演进（如 PIPESIM 的 OLGAS-steady 集成、OLGA 与 SLB 数字化平台的整合、Fluent 的 poly-hexcore 与求解器耦合形式）随年份变化；本矩阵不绑定任何具体版本号。
3. **ChuanFlow 侧**：所有 ✅ 判据均可在仓库内找到对应回归（CTest 名 / §5.x 编号 / 判据数），这是「已对标」的可核实标准；🟡 表示有实现但覆盖面或深度不及目标软件，具体差距以目标软件公开文档为准。
4. **数值凭据**：本文件未引用任何商业软件的发表性能数值（精度对比、求解速度等）；ChuanFlow 内部数字仅用于标注实现位置，不构成与商业软件的精度结论。精度对标需要真实基准算例（如 OLGA 公开校验集——是否可公开获取**待核实**），属后续工作。
5. **对标口径**：Fluent 按「子集」对标（README 定位），即压力基/密度基求解器、两方程湍流、VOF/Euler 双流体、能量与传热的工程子集；不追求 Fluent 全谱系（LES/DES/详细化学/动网格完整族）。
