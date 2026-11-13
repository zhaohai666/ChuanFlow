# ChuanFlow 对标 Pipesim / OLGA / Fluent 能力差异清单

> 本文档记录 ChuanFlow 相对 **Pipesim / OLGA / Fluent** 三款商业软件的现存能力差异，作为"继续排查还有哪些能力没有实现"工作线的单一可信来源。
>
> - 每完成一个增量，本文档同步更新：状态列从 `OPEN` → `IN_PROGRESS` → `CLOSED`，并在末尾的"实施日志"里追加一行实测结果。
> - 与 README 的关系：README 声明"已落地"的能力**不**列在本文档里；本文档只列**尚未落地**或**落地不完整**的能力。
> - 与 `docs/architecture.md` §5.x 的关系：架构文档里的"边界①/②/…"是已命名但未关闭的格子；本文档把它们按**实施成本**重新归组，不再按增量号排列。
> - 与根目录 memory 文件 `project-chuanflow-benchmark-gaps.md` 的关系：memory 文件记录"为什么某条路径走不通"等不可提交的中间知识；本文档只记录"还差什么"与"做到哪一步"。

## 0. 阅读方式

每个差异项包含以下字段：

| 字段 | 含义 |
|------|------|
| **ID** | 唯一编号；前缀字母表示所属类别（A / B / C / D） |
| **能力** | 商业软件里的名字 |
| **参考工具** | 在 Pipesim / OLGA / Fluent 中的对应功能 |
| **当前状态** | `OPEN` / `IN_PROGRESS` / `CLOSED` / `BLOCKED` |
| **阻塞项** | 状态为 `BLOCKED` 时给出阻塞原因 |
| **内核改动** | 需要修改的内核文件 / 函数（如已知） |
| **成本估计** | 以"增量"为单位（一个 §5.x 小节 = 1 增量） |
| **落点** | 完成后写入的 architecture.md 章节号 |

状态机：

```
        ┌──── OPEN ────┐
        │              │
        ▼              │
   IN_PROGRESS ───► CLOSED
        │
        ▼
     BLOCKED  ──（阻塞解除）──► IN_PROGRESS
```

> **落点编号说明（2026-11-11 审计）**：带「（预估）」的 §5.x 落点只是**占位编号**，实际按增量关闭顺序分配，不预留。架构文档当前已用到 **§5.82**（§5.63 = C1/A4，§5.64 = B6，§5.65 = C2，§5.66 = C4，§5.67 = B5，§5.68 = D5，§5.69 = A7，§5.70 = A3，§5.71 = A9，§5.72 = A10，§5.73 = A6，§5.74 = A2，§5.75 = A5 增量 1，§5.76 = A5 增量 2，§5.77 = A5 增量 3，§5.78 = A8 增量 1，§5.79 = A8 增量 2，§5.80 = C3 增量 1，§5.81 = D1 增量 1，§5.82 = D3 增量 1）；因此 A2 原写 §5.63 已被占用，C3 原写 §5.66 与 A6 的 §5.67 不一致。这几条在本轮审计里改为「待分配 / 顺延」，其余预估编号维持不动，等该项真正启动时再定。

## 1. A 类 —— 已命名且已记录的 open cell

> 这些差异在之前的 §5.x 增量里**已经被发现**、**已经被量化**、**已经被写进 architecture.md 的"边界①/②/…"小节**；它们不是新发现，只是尚未关闭。除非阻塞项解除，不要重新从零推导。

### A1. 屈服 + 湍流摩擦闭合（Darby–Metzner / Wilson–Thomas 类）

- **参考工具**：Pipesim / OLGA 非牛顿湍流管段默认行为
- **当前状态**：`BLOCKED`
- **阻塞项**：Wilson & Thomas 2006 (*Can. J. Chem. Eng.* 84:520–525) 与 Thomas & Wilson 2007 (Hydrotransport 17) 两份文献对 `β = τ_w/τ_B` 的区间收缩写法与 `11.6/χ` 钝化项的陈述不一致； adjudicator 已推得 `D U* ρ/μ_B = He^0.5·β^0.5`（对数项必须在 `τ_B → 0` 时塌回 `Re*`），但 β 自身的解方程形式仍待第三方核对。Task #84 保留已失败路径与两份已提取文本（`/tmp/wilson_thomas.pdf`、`/tmp/wilson_rough.pdf.txt`，本箱无 poppler，可用 stdlib zlib 重抽）。
- **内核改动**：`src/1d/src/network.cpp` 增加 `yieldedTurbulentFriction(Re, He, n)`；`src/1d/src/transient.cpp::pipeFrictionPa` 同路径；§5.51 的 `yieldedTransitionOf` 已报告 locus，此处只需写摩擦侧。
- **成本估计**：1 增量（设计已 settled，只差转录核对）
- **落点**：architecture.md §5.73（预估）
- **判据**：§5.51 的 `f·Re_B·g(xi) = 64` 在 4× 与 151× 超过 locus 的两格**必须由新公式复现**（替换而非删除）；"rate is the laminar root" 家族（最差 3.2e-12）同样替换。

### A2. 亚阈值蠕变 / 凝胶强度 / 重启压力

- **参考工具**：OLGA restart-pressure / Pipesim gel-strength
- **当前状态**：`CLOSED`（§5.74）
- **阻塞项**：无（需要一条本构律，不是文献问题）
- **内核改动**：~~`holdsStatic` 改成 `qPrev <= q_creep(t_rest, …)` 的时变判据，并在节点行加一笔凝胶累积项~~ **实施时收窄**：held 判据的形状不变（仍是闩锁处的 `|drive| <= threshold` 比较），变的是**阈值本身**（`4 (tau0 + tau_g(1 - e^{-t_rest/t_g})) L/d`，钟在 `commitStep` 走）；蠕变不进节点行也不进活定律——held 支新增一支**强迫值** `q = pi d^4 drive/(128 mu_eff L)`（与 `Q = 0` 同一行形状），节点行通过连续性自然看到它，无需额外的凝胶累积项。
- **成本估计**：1.5 增量（本构 + 一个启停算例）→ 实测 **1 增量**
- **落点**：architecture.md **§5.74**
- **备注**：§5.49 分支 `Q` 严格为 0、无蠕变、无触变性的缺席已由 §5.64 的注入臂量化；§5.74 关闭它——凝胶臂的开路层比无凝胶孪生臂**晚 22 层（11 s）**（「静止 N 小时重启压差」作为日历可观测），蠕变臂逐层闭式 6.0e-12 且线性 1.25 倍，不开的墙仍是 5.49 的 `Q = 0`，开路的墙停在 5.46 的 Buckingham 律；零声明逐位不变（全量回归为证）。

### A3. β（对流动量项）与静止支同行

- **参考工具**：Fluent 截面动量通量因子
- **当前状态**：`CLOSED`（2026-11-11，§5.70）
- **阻塞项**：无
- **内核改动**：`src/1d/include/cf/1d/Transient.h`（`TransientPipeConfig` 加 `betaIn`/`betaOut`，缺省 1.0；`PipeDyn` 存声明值）+ `transient.cpp`（loader 两键、两条拒收 β < 1 与耦合管声明 β、装填 `d.convective = ρ(βOut−βIn)/A²` 复用既有行路径）+ `transientApp.cpp`（`betaIn ≠ betaOut` 的管条件打印 `dpConv`）+ 绑定三个访问器。**零声明零行为变化**；`pipeFrictionPa` 未动（稳态侧无动量方程，β 不属于那里）。
- **成本估计**：1 增量；实测 ~0.7
- **落点**：architecture.md §5.70
- **被命名 5 次**：§5.34 边界① / §5.35 边界② / §5.36 边界② / §5.37 边界② / §5.38 边界③——全部由 §5.70 关闭（管行认识 β + 与静止支同台）。
- **已测量**：三根几何全同的管只差声明就分开（记账 −6.4e-4 / 退账 +6.4e-4 相对，方向即物理）；全部 501 态被行的隐式二次方程闭式复现（记账正根 / 退账小根，最差 1.32e-9）；账随 q|q| 反号（反驱 |q| 逐位相同）；held pipe 无账可报（Q=0 消掉 q|q|，§5.50 的陈述）。

### A4. 耦合边拒绝流变与屈服

- **参考工具**：Fluent system coupling + 非牛顿伙伴块
- **当前状态**：`CLOSED`（由 §5.63 / task #96 关闭）
- **阻塞项**：无
- **内核改动**：稳态 `network.cpp` / 瞬态 `transient.cpp` 的四处 `coupledResistance > 0 && rheology` 拒收已移除；声明改存 `NetworkPipe::partnerRheology`，1D 壁面律字段清零，四个 coupling apps 转发给 `DetailBlock`。
- **成本估计**：1 增量（与 C1 合并落地，实测 1 增量而非 2.5）
- **落点**：architecture.md §5.63
- **备注**：与 C1 本质相同；A 类与 C 类的区别是 A 在架构文档里有命名边界、C 没有。实际落法不是"把 `pipeFrictionPa` 用完整 key 调用"——耦合边的损失**就是**伙伴的陈述（§5.43），再算一次 1D 摩阻是重复计费；流变只在 3D 伙伴侧生效。
- **既有判据迁移**：5.45 / 5.46 / 5.48 / 5.49 四个 checker 的 `case_rejectCoupled*.dict` 原本断言"非零退出"，首次全量回归因此 4 条失败；现改为**接受 + 与删去三个 key 的孪生档 stdout 逐位相同**，四个 case 头注释同步改写。文件名保留 `reject` 只为兼容旧引用；`ctest -R` 复跑 4/4 通过。

### A5. PVT / 黑油 / 两流体分布式拒绝流变与屈服

- **参考工具**：OLGA + Pipesim 全管线非牛顿
- **当前状态**：**CLOSED**（三个增量全部落地：§5.75 两流体液相幂律、§5.76 黑油/PVT 混合规则决策、§5.77 两流体分布式屈服闩锁）
- **阻塞项**：无（表面比 A4 大）
- **内核改动**：~~`src/1d/src/twoFluid.cpp` + `src/pvt/*` 调度面~~ **增量 1 已落地**（`twoFluid.cpp`：`fluid{}` 键 `flowIndex`/`consistency` + 混合物壁律 `tau_w = alpha·K'γ̇^n + (1−alpha)·μ_G·γ̇`、表观混合黏度进 Churchill、四拒收、`alpha = 1` = 精确单相 MR 律、零声明逐字）；**增量 2 已落地**（共享 `wallRateApparentViscosityPaS`：教材 K 折成 Rabinowitsch `K'`、按段壁率 `γ̇ = 8|v|/D` 重定价液相份额；BB 保留已发表 `lambda` 无滑移加权，持液率/流型/倾斜/`exp(S)`/`Ek` 一字不动；PVT 表仍供给 `rho(p)`/`mu(p)` 压力依赖；loader 撤掉两条拒收、新增一条 `consistency` 无 `flowIndex` 拒收）；**增量 3 已落地**（`twoFluid.cpp`：分布式阈值 `Σ(4dx_i/d)α_i τ₀`、步首 latch 的 held 行 `u = 0`、rest 残差上的 crossing 见证、serial/network 两处 Newton 的单边差分与 Re 钳位邻域层流恒等式、network 单向释放棘轮；`TwoFluid.h` 加 `yieldStress` / `initial.vm` / `staticHold()` 等接口；`twoFluidApp.cpp` 报 `static hold`）。
- **成本估计**：2 增量 → 增量 1 实测 **~0.5**、增量 2 实测 **~1**、增量 3 实测 **~1.5**（原「黑油/PVT 本构律进不去」的判断收窄为「混合规则决策」：只重定价液相黏度数值，不换相关式；增量 3 多出的成本全在 `u = 0` 的不连续上——两处 Newton 的差分、Re 钳位与 network 启动棘轮，都是实测踩出来后补的）
- **落点**：architecture.md **§5.75**（增量 1）、**§5.76**（增量 2）、**§5.77**（增量 3）；§5.66（预估）作废
- **备注**：三行里只有两流体的摩擦是我们自己的律（单站点 `cellFriction`），所以先落地；黑油/PVT 的摩擦是经验相关式，但**黏度只有一个入口**（BB 的 `muM`、PVT 的 `mu(p)`），因此在入口处重定价即可——持液率、流型、倾斜相关式都不读黏度，天然不受影响。`tau_0` 在这两条相关式链里只是 Papanastasiou 有限正则化份额，**没有死支**；真静止阈值由 §5.77 以分布式版落地（`Σ(4dx_i/d)α_i τ₀`，100 判据）。零声明逐位不变由既有 9 项受影响测试 + 既有 twoFluid 全系列 + 全量回归见证。

### A6. 环空 / 矩形 duct（Dodge–Metzner）

- **当前状态**：`CLOSED`（2026-11-12，§5.73）
- **阻塞项**：无
- **内核改动**：已落地——`network.cpp` 截面几何 4 函数 + FB 幂律积分 + Churchill 环空变体 + addPipe 校验/解析 + 两相 march 入口拒收 + loader 键 `dInner`；transient/twoFluid loader 对 `dInner` 拒收。
- **成本估计**：1 增量；实测 ~1
- **落点**：architecture.md §5.73
- **已测量**：f·Re = C(κ) 与 Q = FB 闭式双 1e-9（κ = 2/3、0.1、1e-4）；圆管 f·Re = 64 + Poiseuille 1e-9；C 阶梯 95.74 > 89.37 > 71.78 > 64；幂律 Q 对独立 2048 带重解 5e-6、f·Re_GP = C_pl 1e-9、C_pl 阶梯 → 71.55；Dh 接线 1e-9 全档；四条拒收；18 判据全绿。矩形 duct 已由 C3 增量 1（§5.80）落地；偏心环空仍 OPEN（Dh 抽象已就位）。

### A7. 求解器退出门槛已测未立

- **参考工具**：所有商业软件都有用户可调收敛门槛
- **当前状态**：`CLOSED`（2026-11-11，§5.69）
- **阻塞项**：无
- **内核改动**：`src/1d/include/cf/1d/SolverGates.h`（新建，`cf::gates` 五个 `inline constexpr`：steadyPipeResidualGate 1e-13 / steadyPipeBracketGate 1e-9 / steadyPipeMaxIter 100 / transientNewtonGate 1e-11 / transientNewtonMaxIter 60）；三处站点改用常量名，零行为变化（`network.cpp` 残差门 `:343`、括号门 `:360`、`binghamWallStress` 帽 `:145` 也并入同一管解帽政策；`transient.cpp` 门 `:465`、帽 `:365`）。两相族（`solvePipeFlowTwoPhase`，残差门 1e-7/1e-12·pUp）按声明范围留在 header 之外。
- **成本估计**：0.5 增量（文档 + 一个 policy assertion）；实测 ~0.5
- **落点**：architecture.md §5.69
- **已测量**：§5.47 稳态括号门先于残差门触发，§5.49 瞬态牛顿出口不收敛到不动点而累积缓慢松弛；checker 实测六管闭式足迹 2.4e-13~4.9e-12（括号门名义 1e-9 的实际余量）、瞬态 march 漂移非零 < 1e-9（worst 1.46e-12）、tol 1e-16 逐位不变、maxIter 0 大声失败。

### A8. 1D↔3D 界面多相 / 可压伙伴块

- **当前状态**：`CLOSED`（2026-11-12，增量 1 §5.78 多相冻结相场 + 增量 2 §5.79 等温可压 duct，两增量均已全量回归）
- **阻塞项**：无（「密度随 1D 网络节点压强变化」是关闭时留下的**命名边界**，不是本条验收项）
- **内核改动**：**接口层 `ICoupledElement::lossAt(scalar qMag)` 与 `src/3d` 全程零改动**。增量 1 —— `CoupledNetwork.h` 新增 `PhaseConfig` / `class PhaseField : public IEnergyCoupling, public IEddyViscosity` / `DetailBlock` 的 `phaseCfg_`、`phase_`、`phases()`；`coupledNetwork.cpp` 新增 `readPhases()`、`PhaseField::PhaseField()`、ctor 接线（逐面黏度 hook 互斥）与 `lossAt()` 双钩子转发。增量 2 —— 新增 `CompressibleConfig` / `CompressibleMeasure` / `readCompressible()` / `DetailBlock::solveCompressibleDuct()`，声明时 `lossAt` 整条短路，沿 mesh 自己的 `sectionX_` 逐站以 `ρ = p/RsT`、`ṁ` 守恒推进；四条互斥拒收（phases / 非牛顿 / RANS / 非稳态）。四个 coupling app 转发两种声明，multi app 逐边 `BlockSpec`。**为什么没用 `chflow-3d-compressible`**：它是闭盒 Sod 激波管（无黏 Euler + HLLC，全在匿名命名空间，只有 `main`，无 inlet/outlet），没有库级入口；`src/3d` 的 `phi` 是体积通量、`p` 是 kinematic、`fvm::div`/`fvc::divPhi` 逐字假定常密度，也没有 EOS 类 / 质量通量。
- **成本估计**：2 增量；实测 2（各 ~1）
- **落点**：architecture.md §5.78（增量 1）/ §5.79（增量 2）
- **已测量（增量 1）**：单路径稳态网络固定 `q = 2.4e-4`，使 `dp − dp_plain` 只含相场贡献；dict-only 闭式 `alphaMean 0.5 / rhoMean 700 / muMean 0.1 / bodyDrop 686 Pa` 对四个打印量 1e-9，实测压降差 689.21 Pa（+0.47% 离散足迹）；重力线性 0.500005 / −0.249959、端点密度 200 < 700 < 1200 且 700 ≠ 零声明的 800（防止「忽略相场直接用 rho_ref」假通过）；冻结界面使 `probeSpread` 0.905 → 101.672 Pa；`viscosityMix` 均匀臂梯度 ×1.9955、分层臂 `muMean 0.125` 且 `dp` 落在 978.44 < 1207.44 < 1855.56；零声明孪生与声明-空块对 `dp/q/sweeps/probeSlope/dpPa` 逐位相同；**70 判据 + 7 条大声拒收**、CTest 132 → **133**（`verify.coupledMultiphase`，Test #53，单项 17.68 s）；全量 **133/133 绿，974.68 s**。
- **已测量（增量 2）**：等温可压平面 duct 闭式 `p_in² − p_out² = 24 μ L p_in q_in/(W H³)`（`L` 为 cell-centre 跨度 0.49 m）；主臂 `pRef 10000 / RsT 300000 / q_in 2.4e-4` ⇒ `pOut 9075.24104363 / dp 924.758956369 Pa / ρ 0.0333 → 0.03025 / q_out 2.6446e-4 / ṁ 8e-6 / dp²/dx −3.6e7`，16 个打印量全部 1e-9；5 个站各自的 `p²` 线性、`ρ = p/RsT`、`ρq = ṁ`、`u = q/A` 各 1e-9；七臂各 6 量对闭式；`p(x)` 对直线插值偏离 11.21 Pa（零声明 0.91 Pa）且随速率单调增大；`RsT` 加倍时 `dp` **逐位相同**而 `ρ_in`/`ṁ` 各减半；三档 `pRef` 给 977.56 > 924.76 > 902.36 Pa；三档 `rate` 给 `dp/q` = 3.7598e6 < 3.8532e6 < 4.0732e6（非线性）；双 seed 同一固定点；零声明 twin（同名块放在顶层、不属于 `detail{}`）十一量 **0 容差**；**126 判据 + 8 条大声拒收**、CTest 133 → **134**（`verify.coupledCompressible`，Test #54，单项 31.72 s）；全量 **134/134 绿，1077.29 s**。
- **关闭时留下的边界**：① `pRef` 是 case 声明的绝对入口压强，不随网络节点压强变——`lossAt(q)` 不收压强，要真随网络变必须扩接口并重做准 Newton 论证；② `qMag` 定为入口体积流量，出口是结果，1D 侧仍按常密度 Q 记账；③ 等温 + 局部充分发展 Poiseuille + 忽略轴向惯性，不是可压 Navier–Stokes（无能量方程 / 激波 / 壅塞后物理）；④ 与 `phases{}` / 非牛顿 / RANS / 非稳态互斥，「既分层又可压」未建；⑤ `chflow-3d-compressible` 仍是闭盒 Sod 管，未接成库级求解器。

### A9. 两条耦合边同处串联位置

- **当前状态**：`CLOSED`（2026-11-11，§5.71）
- **阻塞项**：无
- **成本估计**：1 增量；实测 ~0.8（0 内核改动——§5.25/§5.32 的多边内核收任意拓扑，case + checker 为主）
- **落点**：architecture.md §5.71
- **已测量**：串联绑定 q_X1 == q_X2 逐位（72 对）；行可加（rowDrop 逐位闭合、合并行 = 裸和 ≤ 8.5e-5 Pa）；恒定惯量总额恰按 M1/(M1+M2) 分配（1.1e-16）且移动中间节点不动速率；四表四斜率反比于 M 和（0.05%）、严格排序；QN 同 rung；worst 换手两臂都发生；122 判据全绿（Test #71，10.80 s）。

### A10. 行 `M·dq/dt` 双重计费已测未消

- **当前状态**：`CLOSED`（2026-11-12，§5.72）
- **阻塞项**：无（准则已发布：X 是两只括号的函数，不是层数的函数——首层按档阶梯界、塌缩日程 = 账单活着的层 + 一个沉降层 + 容差地板、M_dyn 跨档一致且排他）
- **成本估计**：1 增量；实测 ~0.8（0 内核改动、0 新 case——checker 复用 §5.44 算例把输出重读成政策）
- **落点**：architecture.md §5.72
- **已测量**：首层阶梯 0.6635/0.5602/0.0046 严格排序（对已发表值 1e-9）；M_dyn = 2.143538e6 三档同数（4.7e-10）、0 < M_dyn < M_partner（region 的 28.6%）；rung 在 dq=0 处逐位不可见（层 2..7 spread 0.0）；层 2 = 唯一超地板沉降层（1.07e-2 ≤ 0.012）、层 3 起全部 ≤ 1e-3（最差 8.07e-4）；slab 残账 13.19 Pa（0.46%）→ 推荐选项成文；18 判据全绿（Test #69，192.68 s）。

### A11. 真 MPI 传输 / 域耦合矩阵求解

- **当前状态**：`BLOCKED`
- **阻塞项**：用户 OK 之前不得安装 MPI 工具包（ standing 约束）
- **成本估计**：≥ 3 增量
- **落点**：architecture.md §5.72（预估）

## 2. B 类 —— 零内核改动增量（case + checker only）

> 与 §5.50 / §5.52 / §5.60 同一形态：内核已经能回答这个问题，只是没有一个新的拓扑或新的参数交叉去把它问出来。**首选从这一类开始**，因为风险最低、回归零成本。

### B1. 源 / 汇日程唤醒停泊臂

- **参考工具**：Pipesim scheduled injection / OLGA `MASSFLOW` 边界
- **当前状态**：`CLOSED`（由 §5.61 / task #95 关闭）
- **内核改动**：`src/1d/include/cf/1d/Transient.h::TransientSource` 加 ramp 字段；`src/1d/src/transient.cpp` 加 `sourceRateAt()` helper，`step()` 与 `commitStep()` 把 `nd.source` 改成按 `tNew` 求值；parser 加 `schedule{}` 子块。
- **成本估计**：1 增量（实质是 C 类小钩子，与阀日程同形态）
- **落点**：architecture.md §5.61
- **来源**：§5.56 边界④。
- **拓扑**：一根屈服壁 + 一根带 `valveClose/Open` 的泄放管 + 一个 `source{ schedule … }` 在同一自由节点；日程在 sealed 窗口内把节点抬过该壁的放开压力，让"本来一直 held"的壁被叫醒一次。

### B2. 停泊带下边缘：pB 下移使 NEVER 臂醒来

- **参考工具**：OLGA `INITIAL` 扫描 / Pipesim 边界压力灵敏度
- **当前状态**：`CLOSED`（由 §5.62 关闭）
- **内核改动**：无
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.62
- **来源**：§5.59 边界③ / §5.60 边界④。
- **做法**：把 `B { fixedPressure 2.0e4 }` 改成 `1.95e4`（pB 下移 500 Pa，节点地板被移走，n=0.5 臂 d 在第 1689 步第一次释放、313 个流动状态；臂 c 的完整重锁周期也回来了）。

### B3. HB × 阀日程 × 双向闩锁 三路同台

- **参考工具**：OLGA 带 Herschel-Bulkley 泥浆的阀行程
- **当前状态**：`CLOSED`（由 §5.59 / task #93 关闭）
- **内核改动**：无（§5.52 是 Bingham，§5.54 也是 Bingham；把 `flowIndex` 声明进同一份 dict 即成新算例）
- **成本估计**：1 增量
- **落点**：architecture.md §5.59
- **判据**：§5.52 的"阶梯只动翻步不动是否"必须仍在 HB 下成立；§5.54 的"谁能重睡由邻居阈值次序决定"必须仍在 HB 下成立。

### B4. 分支节点 + 双向闩锁（一节点 ≥ 三根管）

- **参考工具**：OLGA Tee + check-valve 重关
- **当前状态**：`CLOSED`（由 §5.54 关闭）
- **内核改动**：无
- **成本估计**：1 增量
- **落点**：architecture.md §5.54
- **审计**：§5.54 每臂 2 根屈服壁 + 1 根泄放管 = 3 根管 / 节点，11 根壁放开、14 次放开事件、3 次重锁、全文件 17 次翻；重闩由阈值次序决定（弱邻居还在流把节点推过天花板），不由日历决定。文档初版误记 §5.54 为"只一进一出"。

### B5. 三张独立的逐边表 × 逐边 `profile` 同台

- **参考工具**：Fluent system coupling 全特性
- **当前状态**：`CLOSED`（2026-11-11，§5.67）
- **口径**：`profile` 是 `starts.<arm>` 里的**字段**，不是独立表块。真正的第三张表是 §5.43 的 `convectiveToRow{ <臂>{ <边> charge|report } }`。本格要的是 `ddtFrom` + `inertiaFrom` + `convectiveToRow` **三张独立逐边表**连同逐边 `profile` 进**同一份 case**。
- **已落地但仍不够**：`tests/cases/transientCoupledMultiTurbulentProfileInertia/case.dict`（§5.41）已有 `ddtFrom` + `inertiaFrom` + 逐边 `profile`；`tests/cases/transientCoupledMultiTurbulentShapeAndInertiaSplit/case.dict`（§5.44）已有 `inertiaFrom` + `convectiveToRow` + `profile`。审计时全仓无一份 case 同时含三张表。
- **内核改动**：无（三张表两两已交叉：§5.40 / §5.41 / §5.44）
- **成本估计**：1 增量（实际 ~0.5：1 case 文件 + 1 checker，拓扑/网格/dt/seeds 全部复用 §5.41）
- **落点**：architecture.md §5.67
- **实测（2026-11-11）**：`tests/cases/transientCoupledMultiTurbulentTripleTable`，5 臂（base / mix / mix_charge / pin / both_chr）——4 臂是 §5.41/§5.44 的已发布格子，唯一新格 `mix_charge` 是**quasi + slab + charge 三表叠一边**（quasi 边的行用 slab 的 M 付 M·dq/dt、同时让出 dpConv、grade 恒 0，运行验证可行）；base 臂 1e-12 复现 §5.41 p_net 的 level-1/level-8 状态与 pass 序列；80 状态 row 闭合 slack 恰 0；charge 只动一行的一项（mix 的 X1 与 mix_charge 的 X1 17 字段 × 8 层逐位相同）；chord 保持在容差内（最差 8.9e-4）；21 判据 0 失败，`verify.transientCoupledMultiTurbulentTripleTable`。

### B6. 全网络每条边都低于自家阈值（稳态）

- **参考工具**：Pipesim gelled pipeline 启动研究
- **当前状态**：`CLOSED`（由 §5.64 关闭）
- **内核改动**：**无**
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.64
- **实测**：2 份入库 dict + 4 个临时派生臂、**45 条判据 / 0 失败**（CTest `verify.steadyAllBelowThreshold` = Test #45，0.44 s）；全死档 `iterations=0 residual=0 converged`、自由节点 E 随定压均值 25000→25500、阈值相等刀刃 `Q` 仍恰为 0、单管放开时 `iterations` 仍是 0 而该管按 Buckingham 闭式流出 1.37182879207e-6、注入 1e-6 m³/s 时 Newton 跑 19 步且 E = 33978.9475076 Pa 与闭式反解一致。`verify.steadyAllBelowThreshold`。
- **审计修正**：工作区曾有一份未跟踪草稿 `tests/cases/steadyAllBelowThreshold/`，头注释自称"已由 checker 验证 / 已注册 ctest"，但当时两者都不存在，且把 Jacobian 自由节点行写成"全零"（与 `network.cpp:519` 的 HP 电导帽相反）。本案把草稿补成正式算例并改正注释；全死时自由节点压力**不可辨识**（`0 = 0`），报出的 E 只是初始化值。

## 3. C 类 —— 小内核钩子 + 算例（与 §5.48 / §5.51 同形态）

> 与 B 类的差别：需要在某处加一条调度线或新增一个 key；与 D 类的差别：不动方程、不新增模块。

### C1. 流变参数（`flowIndex` / `consistency` / `yieldStress`）进入耦合边

- **当前状态**：`CLOSED`（由 §5.63 关闭，与 A4 同体）
- **内核改动**：4 处 fatal guard 删除；`ICoupledElement::rheology()`；`DetailBlock` 构造接受 `(n, K, τ_0)` 与数值控制；新增 `src/3d/src/solvers/generalizedNewtonian.cpp` 的逐面表观黏度（power-law / HB，Papanastasiou 正则化 + 壁面剪切率覆盖）；4 个 coupling apps 转发
- **成本估计**：2.5 增量 → **实测 1 增量**。原因：3D 侧**没有**改 SIMPLE 的方程装配——`IEddyViscosity::nuEffFace()` 早在 RANS 伙伴块里就是现成的逐面黏度钩子，`GeneralizedNewtonianViscosity` 只是它的第二个实现，`fvm::laplacianFace` 一字未动。原估计把"新增逐面黏度回路"算成了改求解器。
- **落点**：architecture.md §5.63（与 A4 合并）
- **实测**：306 条判据、5 个稳态臂 + 1 个瞬态臂 + 4 个拒收边界；声明但牛顿档与未声明档的求解日志**逐位相同**；1D 侧 `Re` 仍逐位是 `ρvd/μ`；`dp/dp_closed` 在 1.0203 ~ 1.0752
- **新边界**（§5.63 关闭当时的记录）：`turbulence{}` 与非牛顿流体互斥（两者争同一个逐面黏度钩子）；多边同时带流变已转发但未端到端验证；PVT/黑油/两流体仍拒收（A5）——**后续状态**：PVT/黑油两行的拒收已由 §5.76 撤除，两流体幂律由 §5.75 落地，两流体屈服闩锁由 §5.77 落地（A5 三增量全部完成，条目已 CLOSED），见 A5 条目

### C2. HB × 阀 `K`（无日程）

- **当前状态**：`CLOSED`（由 §5.65 关闭）
- **内核改动**：**无**。parser 自 §5.47 / §5.49 / §5.50 起就同时收一根管上的 `flowIndex` / `consistency` / `yieldStress` / `valve`，没有 guard 拒收这个组合——只是从未被问过。成本估计 0.5 → 实际 ~0.3 增量（2 case 文件 + 1 checker + 1 `add_test`）
- **落点**：architecture.md §5.65
- **实测**：56 判据 / 0 失败，16 管 × 2001 状态全部被「HB 剖面积分 + `K ρ|Q|Q/2A²`」递推复现（最差 4.347e-10）；持停管五个数量级 K 逐 token 相同；同一 `K = 500` 在 `n = 0.5` 管上拿走 2.33 Pa、在 `n = 1.5` 管上拿走 308.76 Pa（133×）——**阀的分成由壁律指数决定**；三条停止迹（第 9 / 9 / 5 步）逐状态复现；稳态核对与行进行在同一个两损失根上会合。`verify.transientYieldIndexValveK`
- **边界**：阀 K 仍按几何项计，不随流变变稠变稀（阀有自己的 K，流体有自己的 K_cons）；带日程的组合由 §5.58–5.62 覆盖

### C3. 矩形 duct / 偏心环空（稳态截面）

- **当前状态**：`IN_PROGRESS`（增量 1 矩形 duct 已落地 §5.80；增量 2 偏心环空待启动）
- **阻塞项**：无（矩形）；偏心环空需 bipolar 坐标下的数值 Poisson 锚（避免转录未确认的发表闭式）
- **内核改动**：已落地（增量 1）——`network.cpp` 矩形截面分支（`pipeIsRect` + `pipeArea`/`pipeHydraulicDiameter` 矩形成员 + Fourier 级数 `rectLaminarBracket`/`rectLaminarConductance`/`rectLaminarConstant` + Churchill 复用 + addPipe 校验/解析 + 两相/PVT march 拒收 + loader 键 `ductHeight`）；transient/twoFluid loader 对 `ductHeight` 拒收。
- **成本估计**：2 增量（矩形 ~1，实测吻合；偏心环空 ~1，待启动）
- **落点**：architecture.md §5.80（增量 1 矩形）；增量 2 偏心环空待分配（原写 §5.67 与 A6 撞号、§5.66 与 C4 撞号，均作废）
- **已测量（增量 1）**：五个 aspect（0.01/0.1/0.25/0.5/1）f·Re = C_rect 与 Q = conductance·dp 双 1e-9；C 阶梯 94.705 > 84.676 > 72.931 > 62.192 > 56.908 单调落向正方 56.91（aspect→0 落向平行平板 96）；旋转孪生 0.02×0.04 == 0.04×0.02 1e-9；圆管 f·Re = 64 + Poiseuille 1e-9 逐位；**级数无关第二锚**（纯 Python 二维有限差分 Poisson + Richardson）确认正方 bracket 5.2e-7；Dh 接线 1e-9 全档；湍流档两几何 Re > 4000 且 Q 分离；五条拒收；26 判据全绿。偏心环空仍 OPEN（Dh 抽象已就位）。

### C4. `reportEvery` 与 `endTime` 解耦

- **当前状态**：`CLOSED`（2026-11-11，§5.66）
- **内核改动**：`transient.cpp` + `transientApp.cpp` 步进循环（ceil + 末步裁剪）；parser 零改动（`time { dt; endTime; reportEvery; }` 语法不变）
- **成本估计**：0.25 增量（实际吻合）
- **落点**：architecture.md §5.66
- **备注**：§5.49 以来 1D 瞬态核心是固定 `dt`（`src/1d/src/transient.cpp` 由 `endTime/dt` 定步数、每 `reportEvery` 步采样），隐藏了"求解器自己的 `dt` 能不能跑在打印步之下"的问题。**落地（2026-11-11）**：步数改为 `ceil(endTime/dt − 1e-9)`（近整除容差，如 0.12/0.03 = 4.000000000000001 不误计），末步 `h = min(dt, endTime − (i−1)dt)` 裁剪，`tPrev` 用步索引不累加（末步 t 精确落在 endTime 上），末态恒打印（沿用 `transientApp.cpp` 已有的 `|| i == nSteps`）；`dt`/`endTime` 必须为正，否则拒收。endTime 为 dt 整倍数时逐位向后兼容（`transientYieldValve` 仍 2001 状态、steps=2000）。**审计补充（2026-11-11）**：1D app `transientApp.cpp` 的 `|| i == nSteps` 末步打印与 3D PISO `pisoFoamApp.cpp:148` 的末步裁剪当时已有，本项把二者补进 1D 核心。新增 case `transientTimeDecoupled`（dt 0.03 / endTime 0.1，三管 Newtonian×2 + Bingham）+ checker `check_transient_time_decoupled.py`（六臂 + 拒收，40 判据：全序列隐式欧拉递推复现含裁剪步、reportEvery 3/5 在 4 步运行上仍落在末态、近整除陷阱、compat 孪生、endTime 0 拒收），`verify.transientTimeDecoupled`（Test #83）。

## 4. D 类 —— 大型能力块（新模块或新方程组）

> 这些是 README M1/M2/M3 里命名但树里还没有的能力。标记并排期，不要当作一个增量启动。

### D1. 人工举升设备（ESP / 气举 / 有杆泵）

- **参考工具**：Pipesim / OLGA AL 模块
- **当前状态**：`IN_PROGRESS`（增量 1 多级 ESP + 亲和律已落地 §5.81；气体干涉待已发表相关式确认，气举阀/有杆泵后续增量）
- **成本估计**：3 增量（ESP 多级曲线 + 气体干涉；气举阀动力学；有杆泵井下示功图）
- **落点**：architecture.md §5.81（增量 1 多级 ESP + 亲和律）；增量 2 气体干涉待分配
- **README 现状**：离心泵 + 多级 ESP（`pumpStages`）+ 变速（`pumpSpeed`）。
- **已测量**：63 判据全绿——速率注入逐点曲线（级数阶梯 1:2:3、变速二次项相消）、闭合回路工作点（LoopBase 复现 §5.33 锚 3.085629e-4）、常扬程 N·s²·h0 = 97.2、三条拒收；零声明逐位不变（`steadyEquipmentPump` + 全量回归见证）。

### D2. 蜡 / 沥青质 / 结垢沉积

- **参考工具**：OLGA / Pipesim wax & scale 模块
- **当前状态**：`OPEN`
- **成本估计**：3 增量
- **落点**：architecture.md §6.x
- **依赖**：需要在 `src/pvt` 里增加固相质量平衡，并与热力场耦合。

### D3. 瞬态地层热耦合（径向地层网格）

- **参考工具**：OLGA `HEATFORM`
- **当前状态**：`IN_PROGRESS`（2026-11-13，增量 1 §5.82）
- **成本估计**：2 增量
- **落点**：architecture.md §5.82（增量 1：径向地层导热网格耦合两流体温度方程）
- **README 现状**：只有"恒壁温或地温线"，无径向地层离散。
- **已测量（增量 1）**：`formation{}` 块（conductivity / outerRadius / rhoC / shells）把固定壁换成有限热容径向网格；柱坐标 Fourier 面导纳 `g_j = 2πk·dx/ln(ρ_{j+1}/ρ_j)`，稳态节点温度精确落在对数剖面（几何节点上 = j 的线性），膜阻与地层阻 `R_form = ln(r_e/r_w)/(2πk·dx)` 串联。132 判据（稳态流体递推 + 壁分压 + 对数剖面 + 逐帧地层能量平衡 intake = storage + far + 五拒收）。`verify.twoFluidFormationRadial`（Test #132）。**边界**：逐段独立无轴向岩石导热；算子分裂 O(dt) 瞬态误差（稳态精确）；远场 Dirichlet，线源 E1 瞬态留待增量 2。D3 保持 IN_PROGRESS。

### D4. 流型图作为调度面

- **参考工具**：OLGA 流型图 / Pipesim Beggs-Brill pattern dispatch
- **当前状态**：`OPEN`
- **成本估计**：2 增量
- **落点**：architecture.md §6.x
- **README 现状**：`beggsBrill.cpp` 是单一闭合；流型图本身不是调度面。

### D5. 自适应时间步（CFL 驱动的 `dt`）

- **参考工具**：OLGA adaptive time stepping
- **当前状态**：`CLOSED`（2026-11-11，§5.68）
- **成本估计**：1.5 增量（实际 ~0.5：内核 ~40 行 + 1 case 文件 + 1 checker）
- **落点**：architecture.md §5.68
- **备注**：每个瞬态用例都用固定 `dt`；在最慢管上放一个 CFL 探针能直接消掉 §5.38 记录的 `dt 0.05` 发散墙。**落地（2026-11-11）**：`time { cfl <target>; dtMin; dtMax; }` 三个可选键——`cfl` 在场即切换自适应行进，`dt_next = clamp(cfl·min_k(L_k A_k/|q_k|), dtMin, dtMax)` 从刚提交的流量读出；首步是 case 自己的 `dt`；静止网络取 `dtMax`（静止管不强加约束）；`t` 累加（变步长没有索引可定位），末步裁到 `endTime − t` 且裁剪触发时 `t` 精确吸附到 endTime；`cfl` 缺席 = 5.66 的固定步长逐位不变。三条拒收：`cfl` 声明但 ≤ 0、`dtMin ≤ 0`、`dtMax < dtMin`。Python 绑定暴露三个新字段与 controls 版 `run_transient` 重载。**边界**：这是 1D 管网 Courant 规则，不触碰 §5.38 记录的耦合界面发散墙（那是 3D 块侧的预算，见 §5.38 边界）；case 设计上让全管保持层流（Re_inf ≈ 1345 / 1272），使 checker 的整条调度闭式可复现——湍流档的 CFL 判据未单列。新增 case `transientAdaptiveCfl`（调度访问全部四个 regime：初始 0.004 → 天花板 0.1 ×2 → 规则自定 ×3 → 地板 0.08 ×20 → 余量钳）+ checker `check_transient_adaptive_cfl.py` 16 判据，`verify.transientAdaptiveCfl`（Test #85，0.81 s）。

### D6. 地形诱发段塞（起伏管线低点积液）

- **参考工具**：OLGA terrain-induced slugging
- **当前状态**：`OPEN`
- **成本估计**：2 增量
- **落点**：architecture.md §6.x
- **README 现状**：`pipes{}` 各 section 已经可以倾角，但没有低点积液图与段塞起始准则。

### D7. 组分 PVT 扩展（LBC 黏度、表面张力、含水、CO₂ / H₂S）

- **参考工具**：Pipesim / OLGA 组分 PVT
- **当前状态**：`OPEN`
- **成本估计**：3 增量
- **落点**：architecture.md §6.x
- **README 现状**：明列"黑油分离器、黏度/导热/表面张力关联式、水合物相平衡与生成条件仍待实现"。

### D8. 3D DPM / 拉格朗日粒子 / 群平衡 / 燃烧 / 辐射 / CHT / 动网格

- **参考工具**：Fluent 全套
- **当前状态**：`OPEN`
- **成本估计**：每个子项 ≥ 2 增量
- **落点**：architecture.md §6.x
- **README 现状**：3D 清单止于 VOF + Euler-Euler + RANS。

### D9. 历史拟合 / 参数反演

- **参考工具**：Pipesim Prosper / GAP
- **当前状态**：`OPEN`
- **成本估计**：3 增量
- **落点**：architecture.md §6.x
- **依赖**：需要 `chflow-adjoint` 或 Python 层优化器包前向内核。

## 5. 优先级排序（推荐实施顺序）

1. ~~**B3**~~ — CLOSED（§5.59）
2. ~~**B2**~~ — CLOSED（§5.62）
3. ~~**B4**~~ — CLOSED（§5.54 审计）
4. ~~**C1 / A4**~~ — **CLOSED**（§5.63）：流变进入耦合边，稳态 + 瞬态 + 3D 逐面黏度同台；实测 1 增量而非 2.5。
5. **A1** —— 任务 #84，β 项被第三方文献确认后启动（仍 `BLOCKED`）。
6. ~~**D5**~~ — **CLOSED**（§5.68，`time{ cfl; dtMin; dtMax; }` 三个可选键 + 自适应行进，`cfl` 缺席 = 5.66 固定步长逐位不变；实测 ~0.5 增量而非 1.5）。
7. ~~**B6**~~ — **CLOSED**（§5.64）；~~**C2**~~ — **CLOSED**（§5.65，HB × 无日程阀 K，阀的分成由壁律指数决定）；~~**C4**~~ — **CLOSED**（§5.66，ceil + 末步裁剪 + 末态恒打印，endTime 为 dt 整倍数时逐位兼容）；~~**B5**~~ — **CLOSED**（§5.67，三张逐边表 × 逐边 profile 同台，quasi+slab+charge 三表叠一边）；~~**A7**~~ — **CLOSED**（§5.69，`SolverGates.h` 五常量单一来源 + policy assertion，零行为变化）；**C3** — 增量 1 矩形 duct 已 CLOSED（§5.80）；增量 2 偏心环空按需要穿插。
8. ~~**A9**~~ — **CLOSED**（§5.71，两条耦合边串联：一个速率逐位绑定、行可加、账按 M1/(M1+M2) 分配、斜率阶梯反比于 M 和；实测 ~0.8 增量而非 1）。~~**A10**~~ — **CLOSED**（§5.72，双重计费准则成文：首层按档阶梯界 + 塌缩日程 + M_dyn 排他界，`slab` 成为推荐选项；实测 ~0.8 增量而非 1）。~~**A6**~~ — **CLOSED**（§5.73，环空截面：FB 1958 层流锚 + Churchill 换层流成员 + FB 幂律环空，圆管逐位不变；实测 1 增量如估）。~~**A2**~~ — **CLOSED**（§5.74，凝胶强度 + 亚阈值蠕变：静止阈值增广成只在静止期走的钟、held 支可声明蠕变强迫值、开路层延迟 = 重启压差的日历可观测；实测 1 增量而非 1.5）。~~**A5**~~ — **CLOSED**（增量 1 §5.75 两流体液相幂律 ~0.5、增量 2 §5.76 黑油/PVT 混合规则决策 ~1、增量 3 §5.77 两流体分布式屈服闩锁 ~1.5，三增量全部落地）；~~**A8**~~ — **CLOSED**（增量 1 §5.78 多相冻结相场 + 双钩子 70 判据，增量 2 §5.79 等温可压 duct 126 判据；接口层与 `src/3d` 全程零改动；关闭时留下「pRef 不随网络节点压强变」等 5 条命名边界）；**A11** — 装 MPI 工具包需用户确认，仍 `BLOCKED`。
9. **D1 / D2 / D3 / D4 / D6 / D7 / D8 / D9** — 大型能力块，按产品规划排期。

> **成本估计的修正（§5.63 实测）**：本清单里"3D 侧需要完整非牛顿黏度回路"被估成 2.5 增量，实际 1 增量。教训是先看 `src/3d/include/cf/fvm/EddyViscosity.h`：RANS 伙伴块（§5.23）早已把逐面黏度做成了求解器钩子，`fvm::laplacianFace` 每趟读一次 `nuEffFace()`，因此新增一种流体只需实现该接口，不动方程装配。后续 C 类成本估计应先核对现成钩子再报数。

## 6. 实施日志

> 每完成一个增量，按 `| 日期 | ID | 增量 | commit | 实测 ctest | 备注 |` 追加一行。

| 日期 | ID | 增量 | commit | 实测 ctest | 备注 |
|------|----|------|--------|------------|------|
| 2026-11-11 | — | 文档初版 | 93b2c69 | — | A 类 11 条、B 类 6 条、C 类 4 条、D 类 9 条 |
| 2026-11-11 | B3 | — | — | — | **审计为 CLOSED**：§5.59（task #93, commit 已在闩锁线）= Herschel-Bulkley × 阀 `valveClose`/`valveOpen` 日程 × 双向闩锁三路同台，三档 (n,K)=(0.5,1.0)/(1.0,0.1)/(2.0,0.02)；本节 5.59 起即关闭，文档初版误记为 OPEN |
| 2026-11-11 | B1 | §5.61 | 6581c96 | 1311.16 s / 116 | 源日程唤醒停泊臂：`sources{ schedule{ openRate; close; open } }` 与阀日程同一 ramp helper；1 HB 壁 + 1 阀泄放 + 1 注入源同台；6423 判据、2 case 文件、Test #90 (2.23 s) |
| 2026-11-11 | B2 | §5.62 | ee21027 | 385.66 s / 117 | 停泊带下边缘：pB 20000→19500 使 n=0.5 臂 d 第一次醒来（313 流动状态）、臂 c 重锁周期回来（388/929/1390）；pB 扫描定位转捩 ≈19780；84 判据、3 case 文件、Test #91、0 内核行 |
| 2026-11-11 | B4 | — | — | — | **审计为 CLOSED**：§5.54 每臂 2 壁 + 1 泄放 = 3 管/节点，17 次翻（含 3 次重锁）；文档初版误记 §5.54 为"只一进一出" |
| 2026-11-11 | C1 / A4 | §5.63 | e6b757f | 375.04 s / 117 | 耦合边接受 `flowIndex`/`consistency`/`yieldStress`；新增 `GeneralizedNewtonianViscosity`（幂律 + HB 逐面表观黏度）；4 个 coupling apps 转发；声明但牛顿档求解日志逐位相同；新 checker 5 稳态臂 + 1 瞬态臂 + 4 拒收边界共 306 判据（Test #49）；4 个旧负例 checker 迁移为"接受 + 与无声明孪生档逐位相同"，全量 **117/117 绿**；成本 2.5→1 增量 |
| 2026-11-11 | 文档审计 | — | a7a3a0d | — | **按文档排查代码完成度**：架构 §5.1–5.64（`###` 小节 63 个，另 5.17 是 `##` 标题）、审计时 CTest 117 项（本轮新增 B6 后 118 项）；30 条能力差异中审计当时 **6 条 CLOSED**（A4 / B1 / B2 / B3 / B4 / C1）、**2 条 BLOCKED**（A1 / A11）、其余 OPEN —— 本轮关闭 B6 后为 7 条 CLOSED、21 条 OPEN。7 处文档漂移：B6 孤立草稿虚假声明 + Jacobian 注释反了、B5 标题口径与 §5.41 冲突、A7 三处行号漂移 5/5/147 行、A2 落点 §5.63 与 C1/A4 撞号、A3 落点 §5.64 预留、C3 与 A6 落点不一致、C4 漏记 1D/3D 已有的末步解耦；优先级第 6 项「D5 第一项」与 §2「B 类优先」互斥。B5 收窄为「三张独立逐边表 + 逐边 profile」（全仓无此 case），仍 OPEN |
| 2026-11-11 | B6 | §5.64 | a7a3a0d | 714.41 s / 118 | 全网络全死支的 0 迭代不动点：0 内核改动，2 份入库 dict + 4 个临时派生臂、45 判据；E 随定压均值 25000→25500、阈值刀刃 Q 仍为 0、单管放开时 iterations 仍为 0、注入臂 Newton 19 步且 E = 33978.9475076 Pa 与 Buckingham 闭式反解一致；`verify.steadyAllBelowThreshold`（Test #45，全量里 0.88 s）、全量 **118/118 绿** |
| 2026-11-11 | C2 | §5.65 | fb09958 | 885.06 s / 120 | HB 壁 × 无日程阀 K（`flowIndex`/`consistency`/`yieldStress`/`valve` 同一根管），**0 内核改动**：2 case 文件 + 1 checker、56 判据；16 管 × 2001 状态被「HB 剖面积分 + K ρ\|Q\|Q/2A²」递推复现（最差 4.347e-10）；同一 K=500 在 n=0.5 管拿 2.33 Pa、n=1.5 管拿 308.76 Pa（133×，阀分成由壁律指数决定）；三停止迹第 9/9/5 步逐状态复现；稳态核对同根会合；`verify.transientYieldIndexValveK` |
| 2026-11-11 | C4 | §5.66 | eb794d9 | 885.06 s / 120 | `reportEvery` 与 `endTime` 解耦：步数 `ceil(endTime/dt − 1e-9)`、末步 `h = min(dt, endTime − (i−1)dt)`、`tPrev` 用步索引、末态恒打印、dt/endTime 拒收 0；1 case 文件 + 1 checker、40 判据（六臂 + 拒收）；0.1/0.03 → 0.03/0.03/0.03/0.01 五态落在 t=0.1 精确、0.12/0.03（浮点 4.000000000000001）不误计、endTime 整倍数时 `transientYieldValve` 仍 2001 状态逐位兼容；全序列隐式欧拉递推复现最差 3.27e-11（Bingham 臂）；`verify.transientTimeDecoupled`（Test #83）。**checker 首跑 21 失败根因**：头注释里的示例字样（`endTime 0.1`、`reportEvery 3 and 5`）被 END/REPORT 正则先命中——键提取必须先剥 `//` 注释再读活字典文本 |
| 2026-11-11 | B5 | §5.67 | 95aefb6 | 563.25 s / 121 | 三张独立逐边表（`ddtFrom` + `inertiaFrom` + `convectiveToRow`）× 逐边 `profile` 同台，**0 内核改动**：1 case 文件 + 1 checker、21 判据；5 臂（base/mix/mix_charge/pin/both_chr，拓扑/网格/dt/seeds 复用 §5.41），4 臂是 §5.41/§5.44 的已发布格子，唯一新格 `mix_charge` = quasi + slab + charge 三表叠一边（quasi 边的行用 slab 的 M 付 M·dq/dt、同时让出 dpConv、grade 恒 0——运行验证可行）；base 臂 1e-12 复现 §5.41 p_net 的 level-1/8 状态与 pass 序列；80 状态 row 闭合 slack 恰 0；charge 只动一行的一项（mix 与 mix_charge 的 X1 逐位相同）；chord 最差 8.9e-4 在容差内；`verify.transientCoupledMultiTurbulentTripleTable` |
| 2026-11-11 | D5 | §5.68 | 74066e3 | 743.73 s / 122 | CFL 自适应时间步：`time{ cfl; dtMin; dtMax; }` 三个可选键，`cfl > 0` 切换自适应行进 `dt_next = clamp(cfl·min_k(L_k A_k/\|q_k\|), dtMin, dtMax)`（从刚提交的流量读出）、首步 = case 的 dt、静止网络取 dtMax、t 累加、末步裁到余量且 t 精确吸附 endTime、`cfl` 缺席 = 5.66 固定步长逐位不变；内核 `TransientNetwork::cflStep` + core/app 各一条自适应分支 + 三条拒收 + Python 绑定新字段与 controls 重载；1 case 文件 + 1 checker、16 判据——调度逐步被「上一态流量的 Courant 规则」复现（调度是物理的性质不是记住的表）、全部 27 态被逐步 h 的隐式欧拉行复现（最差 8.53e-10，全管层流 Re_inf ≈ 1345/1272）、固定行进 501 态闭式复现且与自适应末态不同、调度访问全部四个 regime、三条拒收；`verify.transientAdaptiveCfl`（Test #85，0.81 s）。**case 设计教训**：rule_inf = cfl ρL²/(drive·τ) 与 Re_inf = ρd³·drive/(32 μ²L) 联立定 drive——floor 绑定与层流可以同时成立，但 drive 要按这两条 identity 选（第一版 μ=0.01 时 Re 14280 进湍流区，checker 层流闭式全面失配 0.83 才暴露） |
| 2026-11-11 | A7 | §5.69 | 825b200 | 3471.66 s / 123 | 求解器退出门门槛策略：`src/1d/include/cf/1d/SolverGates.h` 五个 `inline constexpr` 单一来源（steadyPipeResidualGate 1e-13 / steadyPipeBracketGate 1e-9 / steadyPipeMaxIter 100 / transientNewtonGate 1e-11 / transientNewtonMaxIter 60），`network.cpp` 三处 + `transient.cpp` 两处站点改用常量名、`binghamWallStress` 的保险帽并入同一管解帽政策，**零行为变化**；两相族门按声明范围留在 header 外。checker 26 判据：header ↔ 文档值互证（名字与值同处一行）、站点引用 + 旧裸字面量消失、六管闭式足迹 2.4e-13~4.9e-12（括号门名义 1e-9 的实际余量）、tol 1e-16 逐位不变（内门是地板）、maxIter 0 大声失败、瞬态 march 漂移非零 < 1e-9（worst 1.46e-12，§5.49 声明的可测形式）；`verify.solverGates`（Test #86，0.29 s） |
| 2026-11-11 | A3 | §5.70 | dc7dd7f | 2579.04 s / 124 | β（对流动量项）可声明化：管块可选键 `betaIn`/`betaOut`（缺省 1.0），瞬态行装 `dpConv = ρ(βOut−βIn)/A²·q\|q\|`——与 §5.43 接口装填的 `convective` 同一笔账，推广到 case 自己拥有；两条拒收（β < 1 "no profile is fuller than flat"、耦合管声明 β——接口独占那笔账）；app 条件打印 + 绑定 3 访问器，**零声明零行为变化**、`pipeFrictionPa` 未动。checker 16 判据：恒等哨兵（βIn=βOut=1.0 与不声明逐位相同）、三根几何全同管只差声明就分开（记账 −6.4e-4 / 退账 +6.4e-4，方向即物理）、501 态全部落在行的隐式二次方程闭式根上（记账正根 / 退账小根 + 判别式断言，最差 1.32e-9）、账随 q\|q\| 反号（反驱 \|q\| 逐位相同）、屈服 + β 一行（tw 空间二分）、held pipe 无账可报（Q=0 消掉 q\|q\|）、四条拒收；`verify.transientBetaMomentum`（Test #87，0.70 s）。关闭 §5.34①/§5.35②/§5.36②/§5.37②/§5.38③ 五处命名 |
| 2026-11-11 | A9 | §5.71 | b835dcb | 580.48 s / 125 | 两条耦合边串联（拓扑 A → X1 → J1 → X2 → Snk，6 臂 = 四表 × GS + bore/mixX2 × QN），**0 内核改动**：1 case 文件 + 1 checker、122 判据。串联绑定 q_X1 == q_X2 逐位（72 对）、行可加（rowDrop == R·q + M·dq/dt 逐位 0，合并行 = 裸和 ≤ 8.5e-5 Pa）、惯量总额恰按 M1/(M1+M2) 分配（1.1e-16）且移动中间节点（mixX2 J1 1737.9 < mixX1 2529.7）、四表四斜率反比于 MΣ（0.05%，τ = MΣ/RΣ ≈ 0.52 s 线性段，行程 0.12%）、QN 同 rung（≤ 1%）、worst 换手两臂都发生；`verify.transientCoupledMultiBothSeries`（Test #71，10.80 s）。**case 设计教训**：drive 4.8e5 把 Q_inf 推出块 SIMPLE 60 扫稳定窗、爬升穿过渡 Re 让 kOmegaSST 块炸（continuity 1e55）→ 降 6.0e3；τ = MΣ/RΣ ≈ 0.52 s（真实湍流弦 RΣ = 3.2e7），dt 5e-5 × 12 层只走 0.12% 行程 → 阶梯读线性斜率而非指数 τ；uRef 对 initialFlow 高度敏感（4.5e-5 → 0.058 vs 6.5e-5 → 3.52 第一趟即炸） |
| 2026-11-12 | A10 | §5.72 | 0f8c0de | 811.17 s / 126 | 双重计费准则成文：「多少层才允许 rowError 低于 X」的解析答案——X 是两只括号的函数不是层数的函数。**0 内核改动、0 新 case**（checker `check_transient_row_error_criterion.py` 复用 §5.44 算例，CTest `verify.transientRowErrorCriterion`）、18 判据：账单 = (M_row − M_dyn)·dq/dt 且 M_dyn 是块的属性（三档同读 2.143538e6 kg/m⁴，散布 4.7e-10；0 < M_dyn < M_partner，region 的 28.6%）、首层按档阶梯界（0.70/0.60/0.006，实测 0.6635/0.5602/0.0046 严格排序、对已发表值 1e-9）、rung 在 dq=0 处逐位不可见（层 2..7 spread 0.0、层 8 移 1.1e-4）、塌缩日程（层 2 = 唯一超地板沉降层 1.07e-2 ≤ 0.012、层 3 起全部 ≤ 1e-3 最差 8.07e-4、沉降残差 64× 衰减）、charge 不串扰（层 1 ≤ 1.4e-7）、自由边恒等 2.3e-13 + 账单恒正 + 价格帽（余量 ≥ 36%）、弦差守容差（96 状态最差 0.0235 Pa）；**slab 推荐选项成文**（残账 13.19 Pa = 0.46%）。checker 全程 192.68 s（Test #69） |
| 2026-11-12 | A6 | §5.73 | 33fd1b2 | 452.09 s / 127 | 环空截面进稳态管网：可选键 `dInner`（`diameter` 保持外径），水力直径 Dh = d − dInner 进 Re/粗糙度/Darcy 换算，层流锚 = FB 1958 精确闭式（牛顿）与 FB 幂律环空（λ 二分 + 剖面积分换序，addPipe 解析一次），湍流 = Churchill 换层流成员（C(κ)/8 替代 8，C = 64 时 IEEE 逐位还原）；幂律 Re_GP 用 DM 原形（Dh 上、无 Rabinowitsch），C_pl 由"层流渐近线 = FB 剖面"反解、K 相消、κ→0 成员 71.55。内核 ~130 行 + 1 case + 1 checker、18 判据：f·Re = C(κ) 与 Q = FB 闭式双 1e-9（κ = 2/3、0.1、1e-4）、圆管 f·Re = 64 + Poiseuille 1e-9、C 阶梯 95.74 > 89.37 > 71.78 > 64、幂律 Q 对独立 2048 带重解 5e-6、C_pl 阶梯 90.28 > 84.77 > 71.55、Dh 接线 1e-9 全档、四条拒收大声点名；圆管分支逐位不变（全量回归为证）。拒收范围：屈服、耦合边、黑油/PVT march、transient/twoFluid loader。**case 设计教训**：幂律臂不能放 dp = 1 Pa 精确档——n = 0.5 的 Re_GP 落进内核 1e-12 的 Re 地板被扭曲，幂律精确锚在 dp = 1e4 档（湍流项 < 1e-100 相对） |
| 2026-11-12 | A2 | §5.74 | 3c6869c | 464.83 s / 128 | 凝胶强度 + 亚阈值蠕变进瞬态行：可选键 `gelStrength`/`gelTime`/`creepRate`（只对已声明 yieldStress 的管有意义），静止阈值增广成只在静止期走的钟 `4(τ0 + τ_g(1−e^(−t_rest/t_g)))L/d`（commitStep 累积 tRest、一步流动清零），held 支新增蠕变强迫值 `q = πd⁴·drive/(128·μ_eff·L)`（μ_eff = μ + τ0/γ̇_c，闩锁 drive 上、与 Q=0 同一行形状），crossing 管停在 Q=0 从下一层闩锁进蠕变支。内核 ~70 行 + 1 case + 1 checker、19 判据：Plain 第 60 层 / Gel 第 82 层开路（checker 用同一闩锁递推预测整数层，前一步 8200 ≤ thr(40) = 8202.7 闭式锁定）、Gel 晚 22 层（11 s）= 重启压差的日历可观测、蠕变两臂 600 步逐层闭式 6.0e-12 且线性 1.25 倍（9.1e-12）、不开的墙（80000 vs 可达 61500）全程 Q=0、开路的墙停在 Buckingham 根 dp* = 6142.6706 Pa（tw* = 3.0713353，非牛顿的 6.5 Pa）、五条拒收大声点名；零声明 = 5.49 逐位不变（全量回归 128/128 为证）。`verify.transientGelCreep`（Test #90，1.15 s）。**case 设计教训**：固定节点压下凝胶只影响闩锁决策——可观测面是"源注入抬压 → 开路层延迟"而非稳态压差，五个独立源腔（held 期恰 102.5 Pa/步）让开路层成为整数；注入速率取非整值（2.05e-7）让台阶踩不到阈值的逐位边界（Plain 在 6000 整处的逐位相等因此被绕开） |
| 2026-11-12 | A5（增量 1） | §5.75 | 046cf79 | 459.68 s / 129 | 两流体分布式行认识液相幂律流变：`fluid{}` 可选键 `flowIndex`/`consistency`（液相的），混合物壁应力 `τ_w = α·K'γ̇^n + (1−α)·μ_G·γ̇`（γ̇ = 8v_m/d，K' = K((3n+1)/4n)^n Rabinowitsch 修正），表观混合黏度 τ_w/γ̇ 进同一支 Churchill——层流成员 f = 64/Re 恰好读回 (4dx/d)·τ_w（与未声明分支同一条恒等式）；α = 1 = 精确单相 MR 律、零声明走旧路径逐字（17 个既有 twoFluid case 为证）。内核 ~45 行 + 1 case 目录（4 文件） + checker 新模式、17 判据：dp = 240.16 Pa 对全 blend 律与纯层流闭式双 1e-9、K' 不是教材 K（差 = 液体份额修正 25.3375 Pa 1e-9）、flowIndex 1/K = μ_L 与无声明孪生全场一致、纯液臂 dp = 300 精确、四条拒收（yieldStress"does not take a threshold yet"——分布式阈值要 5.49 式闩锁自成增量）。`verify.twoFluidRheology`（Test #121，0.58 s）。**case 设计教训**：① wallis us = 0 让 α* = vsl/vm 成为通量恒等式，律的全部输入都是 case 常数；② 初始 holdup 冲刷是 vm 上的 sharp front + 数值扩散尾——endTime 60 留 5e-4 缓弛尾（闭式判据假红），endTime 600 才精确落位；③ 黑油（Beggs-Brill）与 PVT march 的摩擦是经验相关式，本构律进不去——这两行保持 OPEN 等混合规则决策，A5 由 2 增量收缩为「增量 1 ~0.5 + 待决策的剩余」 |
| 2026-11-12 | A5（增量 2） | §5.76 | 3135635 | 2884.71 s / 131 | 黑油/PVT 混合规则决策：**只重定价液相那一份黏度**。共享 `wallRateApparentViscosityPaS`（教材 K → Rabinowitsch `K'`，按段壁率 `γ̇ = 8\|v\|/D` 定价，Papanastasiou `m = 100` 正则化屈服）；BB 保留已发表 `lambda` 无滑移加权、持液率/流型/倾斜/`exp(S)`/`Ek` 一字不动，5 个 march/探针调用点统一转发声明；PVT 表仍供给 `rho(p)`/`mu(p)` 压力依赖；loader 撤掉黑油/PVT 对 `yieldStress`/`flowIndex` 的两条拒收、新增 `consistency` 无 `flowIndex` 一条。内核 ~90 行净增（4 文件 109 增 / 19 删）+ 72 行单元锚 + 2 case 目录 + 2 checker + 2 `add_test`，**121 判据**：`steadyPvtRheology`（平表 → `τ_w = 0.25 Pa` 是几何量，五律各自独立反解；未声明 = Hagen–Poiseuille `Q = 1.96349540849e-05`，`f·Re = 64` 全档，`K'/K = 1.118034` 排除教材 K 的 1.25× 速率）54 条、`steadyBlackOilRheology`（干油 + 水平 → `lambda = HL = 1`、`exp(S) = 1`、`Ek = 0`，checker 只从 dict 重算 `Bo(p)`/`rhoOil(p)`/`muL` 并逐段重推进 24 步，落点最差 0.081 Pa < 割线门 0.1005 Pa、出口速度 1e-11、`m/Q = 850`、`regime` 恰是已发表 `fr` 图）60 条、`networkTests.cpp` 7 条精确解析锚（`(4/D)K'γ̇^n` 1e-12、HB 两份额之和、声明黏度不能移动 holdup/regime/`exp(S)`）；两份 case 各四条拒收大声点名；零声明逐位不变由既有多相/PVT/屈服等 9 项 + 全量回归见证。**边界**：`tau_0` 在相关式链里只是有限表观黏度、**没有死支**（真闩锁是 5.49 / A5 增量 3）；湍流半区未独立验证（Re 1.649 ~ 518.6 刻意深层流）；环空与耦合边拒收不变。**case 设计教训**：BB 的经验相关式让严格闭式很难，靠「干油 + 水平 + 两端定压」把 `lambda/HL/S/Ek/ele` 全退化掉、混合物密度在层流下相消，只留 `Bo(p)` 一个 1.8e-5 的漂移，checker 因此能重推进而不是拟合；PVT 则用平表把这唯一的漂移也消掉，`τ_w` 变成纯几何量 |
| 2026-11-12 | A5（增量 3） | §5.77 | 6a388ef | 834.76 s / 132 | 两流体分布式行拿到**静止的墙**，A5 至此 CLOSED：阈值 `dp_y = Σ_i (4dx_i/d)·α_i·τ₀`（只有液相份额屈服，均匀 α 退化成 `4τ₀L/d·α`，**不认 n / K**），drive = `dpDrive − Σρ_m g sinθ dx`；步首 latch 的 held 行 `u = 0`、`staticHold()` 报 1、`nodePressure()` 报摩擦损失 0（静态壁扛全部压差）；流动支 `τ_i = α_i(μ_L γ̇ + τ₀) + (1−α_i)μ_G γ̇` 在 `v → 0+` 从上方趋于阈值。内核 ~310 行净增（`TwoFluid.h` 57 / `twoFluid.cpp` 241 / `twoFluidApp.cpp` 11，含注释）+ 7 份 dict + `check_two_fluid.py` 新模式 379 行 + 1 `add_test`，**100 判据**：主臂 `dp_y = 400 Pa > dp 300` 全程 held（`p = 400` 边界仍 held）；`case_open` 由独立全 blend Churchill + 隐式 Euler 递归逐帧复现 vm / dp，稳态 `(dp−dp_y)d²/(32μ_L L) = 0.0125`、`f·Re = 64`；**分布式见证** `case_graded` 半管液 ⇒ `dp_y = 200 Pa`，同一 300 Pa 开路，`p = 199/201` 把带边钉在 ±1 Pa（Godunov 扩散使 α 漂移 6.96e-3 ⇒ 阈值漂移 0.70 Pa，判据显式核算）；零声明 `case_newt` 逐帧复现旧层流启动、`yieldStress 0` 哨兵逐位相同；`case_cross` 从 `initial vm 0.01` 单调减速、t = 0.14 停住转 held、全程无反向流；`case_hb`（n 0.5 / K 0.1 / τ₀ 1）阈值仍 400、稳态 `vm* = 0.05`；`case_network` 两 span 全 held 且 `junction net = 0`，1200 Pa 探针验证**单向释放棘轮**让两管都流动并与 `case_open` 逐位吻合；三条拒收大声点名。**case 设计教训**：① 屈服壁律在 `u = 0` 不是函数，中心差分跨零把阈值跳跃塞进斜率，serial 速度在 `8.9e-05 ~ 1.8e-03` 永久振荡、network 直接报 singular Newton matrix → crossing 改判「运动支上是否还有根」（rest 残差变号）、`momentumResidual(dt,0,dp)` 取带符号阈值、两处 Newton 在 `\|u\| <= h` 改单边差分；② `churchillFriction` 的 Re 钳位 1e-12 会抹掉屈服阈值下限（`Re ~ v²`），需在 `Re < 1e-9` 邻域直接用层流恒等式 `(4dx/d)τ_w`；③ 两端定压是唯一能观察「不开 / 开路」的面（`phaseRates` 强加、`sealed` 恒 0），`initial vm` 让运动→静止的 crossing 不需要时间变化边界；④ junction 压力是 Newton 未知量，步首 latch 会过期，必须单向释放，否则超阈值网络永久冻结。零声明逐位不变由既有 twoFluid 全系列 + 全量回归见证 |
| 2026-11-12 | A8（增量 1） | §5.78 | ee791a6 | 974.68 s / 133 | 1D↔3D 伙伴块第一次认识两相：**接口层 `ICoupledElement::lossAt(q)` 零改动**、`src/3d` 求解器零改动，只在 `DetailBlock` 装一个冻结相场并复用两个现成钩子——`IEnergyCoupling::buoyancyForce()` 给逐 cell 加速度 `g·rho_m/rho_ref`（FVM 的 p 是 kinematic，`measure()` 乘回 `rho_ref` 后才是物理静压头 `−L_meas <rho_m> g_flow`，只有液体那份是重的），`IEddyViscosity::nuEffFace()` 让 `viscosityMix` 用同一份 alpha 逐面给壁律定价；逐面 hook 只有一个 owner，与非牛顿/RANS 同时占用即大声拒收（§5.63 规则延伸）。内核 `CoupledNetwork.h` 78 行 + `coupledNetwork.cpp` 169 行 + 4 app 26 行 + 1 dict 96 行 + checker 308 行 + 1 `add_test`，**70 判据**：单路径稳态网络固定 `q = 2.4e-4`，使 `dp − dp_plain` 只含相场贡献（第一版照 §5.18 抄并联旁路，声明臂与孪生臂 q 不同、压降差混进分支重分配，是首要 case 设计教训）；打印的 `alphaMean 0.5 / rhoMean 700 / muMean 0.1 / bodyDrop 686 Pa` 对 dict-only 闭式 1e-9，实测压降差 689.21 Pa（+0.47% 离散足迹，全液 1176 → 1167.36、全气 196 → 194.56 也在 2% 内显式核算）；重力线性 0.500005 / −0.249959、端点密度 200 < 700 < 1200 且 700 ≠ 零声明的 800（**判据锐化**：rhoG 若取 400，rhoMean 会恰好等于参考密度 800，「忽略相场直接用 rho_ref」的实现也能骗过密度判据）；冻结界面使 `probeSpread` 0.905 → 101.672 Pa（> 10×，证明体力逐 cell 加）；均匀混合臂发展段梯度 ×1.9955、分层臂 `muMean 0.125` 且 `dp` 落在 978.44 < 1207.44 < 1855.56；零声明孪生与声明-空块对 `dp/q/sweeps/probeSlope/dpPa` 逐位相同（0 容差）。**测量锚教训**：`bodyDrop` 的跨度是首末 **cell centre**（50 cell → 0.49 L）而非 bbox L，用 0.5 会算出 800 Pa 而实测 784 Pa，正好差一个 cell。7 条大声拒收（负密度/黏度、alpha 越界或不足两带、alpha 与 alpha0 互斥、gravity 三分量、viscosityMix 词、抢 hook）；9 臂与 7 探针全写系统临时目录，case 目录不留临时 dict；`verify.coupledMultiphase`（Test #53）。**边界**：alpha 冻结、动量仍是常 `rho_ref` 不可压、可压 / 变密度伙伴块留待增量 2，A8 当时保持 `IN_PROGRESS`，增量 2 见下一行 |
| 2026-11-12 | A8（增量 2） | §5.79 | 842dca3 | 1077.29 s / 134 | 伙伴块第一次可压，**A8 至此 CLOSED**：`detail{ compressible{ pRef; RsT; mu; } }` 让 `lossAt(q)` 整条短路到沿 mesh 截面 march 的等温理想气体 duct（不是把闭盒 Sod 管接进来——`chflow-3d-compressible` 无库级入口、无 inlet/outlet，`src/3d` 的 `phi` 是体积通量且 `fvm::div` 逐字假定常密度）。**接口层 `ICoupledElement::lossAt(scalar qMag)` 与 `src/3d` 零改动**；`p_in^2 − p_out^2 = 24 μ L p_in q_in/(W H^3)`（L 为 cell-centre 0.49 m），`qMag` 定为入口体积流量、`q_out` 是结果；march 逐站由 `ρ=p/RsT`、`ṁ` 守恒重算 `u` 与局部壁摩阻 `dp/dx=−12μu/H²`，再以**局部** `2p(dp/dx)` 推进 `p²`（不是读闭式直线）——因 `p·u=ṁRsT/A` 沿程恒定、局部梯度逐站相同，故 `p(x)` 弯而 `p^2(x)` 严格线性（落在解析线上到舍入）、`rho·q = mDot` 逐站守恒，`p2Slope` 只留作 checker 参考、出口态取 march 末站。内核 header 45 行 + cpp 170 行 + 4 app 50 行 + 1 dict 98 行 + checker 377 行 + 1 `add_test`，**126 判据**：主臂 16 个打印量对 dict-only 闭式 1e-9（`dp 924.758956369 Pa`、`rho 0.0333 → 0.03025`、`q_in 2.4e-4 → q_out 2.6446e-4`、`mDot 8e-6`、`d(p^2)/dx −3.6e7`），5 站 × 4 恒等式各 1e-9，七臂各 6 量对闭式；`RsT` 加倍时 `dp` **逐位相同**而 `rho_in`/`mDot` 减半（RsT 在 dp 里相消）；三档 `pRef` 977.56 > 924.76 > 902.36 Pa（**pRef 刻意取低压 1e4**，取 1e6 会让密度只变 0.09%、常密度假实现藏进容差，与增量 1「rhoG 别取 400」同一条教训）；三档 rate 的 `dp/q` 严格递增（非线性）、双 seed 同固定点；`probeSpread` 11.21 Pa 对零声明 0.91 Pa 且随速率单调增大；零声明 twin 把同名块搬到顶层（不属于 `detail{}`）对 11 量 **0 容差**；8 条大声拒收（负值 / 壅塞 / 与 phases、非牛顿、RANS、非稳态抢本构）。**边界**：`pRef` 是 dict 声明的绝对入口压强、不随网络节点压强变（要变必须扩 `lossAt(q,pUp)` 并重做准 Newton 论证）；等温 + 局部充分发展 Poiseuille + 忽略轴向惯性，不是可压 Navier–Stokes；与 `phases{}` 互斥故「既分层又可压」未建。`verify.coupledCompressible`（Test #54） |
| 2026-11-13 | C3（增量 1） | §5.80 | 28ff301 | 1013.58 s / 135 | 矩形 duct 进稳态管网摩阻：可选键 `ductHeight`（`diameter` 保持宽度）， bore 成 W×H 矩形。层流锚是**自推精确 Fourier 级数**（截面 Poisson `grad²u = −G/μ`、no-slip）：`bracket(aspect) = 1 − (192α/π⁵)Σ_{n odd} tanh(nπ/2α)/n⁵`（α = short/long），`addPipe` 解析一次（4001 奇项、尾 < 1e-15）；`Q = (dp/L)·long·short³·bracket/(12μ)`、`C = f·Re_Dh = 96/((1+α)²bracket)`（平行板 96、正方 56.91）。`Dh = 2WH/(W+H)` 进 Re/粗糙度/L-Dh；湍流 = Churchill 换层流成员（generic `churchillAnnulusFriction`，C = 64 时圆管 IEEE 逐位还原）。内核 ~148 行（`Network.h` 20 / `network.cpp` 118 / `transient.cpp` 4 / `twoFluid.cpp` 6）+ 1 case + 1 checker + 1 `add_test`，**26 判据**：五档 `f·Re = C_rect` 与 `Q = conductance·dp` 双 1e-9、单调阶梯 94.705 > 84.676 > 72.931 > 62.192 > 56.908、薄档→96（< 1.5）、正方 56.908（< 0.01）、**series-free 2D FD Poisson 交叉验证**（单位正方 SOR、面积分 = bracket(1)/12、Richardson(48,96) 5.2e-7）、旋转孪生 Q/v/Re/f 1e-9、圆管孪生 `f·Re = 64` + Poiseuille 1e-9、Dh 接线 1e-9 全档、湍流档双 Re > 4000 且 Q 按面积分流、五条拒收大声点名（环空组合/屈服/幂律/耦合边/负 ductHeight）。零声明逐位不变（既有 `steadyEquipmentPump` 等全套 + 全量回归 135/135 为证）；拒收范围另含两相/PVT march、transient/twoFluid loader。`verify.steadyRectDuct`（Test #94，1.68 s）。**case 设计教训**：checker 阶梯顺序必须 thin→square（`consts[0]` = 薄档≈96、`consts[-1]` = 正方 56.91），第一版 square→thin 让单调/极限三断言取错臂假红——内核正确（per-arm `f·Re = C` 与 FD 交叉验证全过），仅断言索引错。**边界**：牛顿单相稳态；偏心环空留待 C3 增量 2（数值 bipolar-Poisson 锚） |
| 2026-11-13 | D1（增量 1） | §5.81 | 9b54720 | 696.52 s / 136 | 多级 ESP + 亲和律进稳态泵曲线：泵管两个可选键 `pumpStages N`（≥1，N 级串联在同一流量上叠加扬程，能量守恒）与 `pumpSpeed s`（>0，Euler 相似律 Q∝s、H∝s²）。组合曲线 `H(q) = stages·(s²h0 − (h0−h1)(q/q1)²)` 截断在 0；`pumpRate 0` 时常扬程 `stages·s²·h0`。**两条物理自推**（非转录经验相关式）：级数叠加 = 能量守恒；亲和律 `H_s(q) = s²H₁(q/s) = s²h0 − (h0−h1)(q/q1)²`——二次项的 s² 恰好相消、只有截止扬程缩放。`stages1/speed1`（或缺省）与 §5.33 泵律 **IEEE 逐位相同**（`x*1.0==x`），零声明一位不动。内核 ~30 行（`Network.h` 10 含注释 / `network.cpp` `pumpHeadOf` 6 + loader 键 2 + 校验 9）+ 1 case（22 节点 11 管）+ 1 checker + 1 `add_test`，**63 判据**：速率注入臂在已知 q=2e-4 上逐点读曲线（级数阶梯 FedS1:FedS2:FedS3 = 1:2:3 1e-9、变速臂 FedV12−FedV08 = −ρg(1.2²−0.8²)h0 二次项相消见证、自由源节点 p_src = p_term+Rq−ρgH 闭式 1e-6、终端被 honoring）、闭合回路工作点 = ρgH(q)=Rq 正二次根（LoopBase 复现 §5.33 锚 3.085629e-4 1e-9、2Stage/Speed0.8/Comb 各对自己根）、常扬程 stages3/speed0.9 给 H=97.2 m 精确 + q=ρgH/R 线性、三条大声拒收（无泵声明 stages/speed、speed 0、stages 0）。零声明逐位不变由 `steadyEquipmentPump`（Test #41）+ 全量回归 136/136 见证。`verify.steadyEquipmentEsp`（Test #43，0.16 s）。**case 设计教训**：速率注入终端取 **10 bar 而非 1 bar**——3 级 ESP 在 q=2e-4 处扬程 100 m（≈8.3 bar），1 bar 终端下自由源节点 p_src=p_term+Rq−ρgH 被压成**负压**（FedS1≈−0.48 bar、FedS3≈−6 bar），而管网 Newton 可行性守卫（`network.cpp` globalization `pTry>0` 否则跳过）拒收非正压强 → 第一版 6 条注入臂全不收敛（Q≠2e-4、源节点卡正压侧），**闭合回路臂却全精确**，从而证明内核律正确、隔离为 case 设计而非物理；抬高终端到 10 bar 后所有源节点为正（最小≈2.97 bar），Q/dpPump/dpFriction/Re 全不变（只依赖强加流量与曲线）。另回路 Re 带下界取 40（speed0.8 臂 Re≈48）。**边界**：稳态单相牛顿泵曲线（黑油/PVT march 经验相关式不进）；气体干涉（head 随 GVF 退化）需已发表相关式、留待增量 2；亲和律是几何相似+不可压推论，未含效率/Re 随转速偏移；级数叠加假定各级几何相同同一流量（未含级间泄漏）。D1 保持 IN_PROGRESS |
| 2026-11-13 | D3（增量 1） | §5.82 | 24c0d6a | 1158.91 s / 137 | 径向地层导热网格耦合两流体温度方程：`thermal{}` 之外可选 `formation{}` 块（`conductivity k` / `outerRadius r_e` / `rhoC` / `shells M`），把固定壁（无限等温热库）换成**有限热容径向岩石网格**。声明后 `wallTemperature(i)`（constant/geothermal）改读为**远场 T_e**（r_e Dirichlet），壁温由网格解出。每段 M+1 几何节点 `ρ_j=r_w(r_e/r_w)^{j/M}`（r_w=D/2），node0=壁经膜 `kw=k'dx` 耦合流体、nodeM=远场。**自推柱坐标 Fourier**：面导纳 `g_j=2πk·dx/ln(ρ_{j+1}/ρ_j)`（柱壳精确热导，非平板近似），节点容量 `C_j=rhoC·dx·π(upper²−lower²)`（总=`rhoC·dx·π(r_e²−r_w²)`），隐式 Euler 径向 march 与流体温度步**算子分裂**（稳态不动点精确）。**稳态锚**：常热流逐壳守恒→`T_j=T_w+(T_e−T_w)ln(ρ_j/r_w)/ln(r_e/r_w)`，几何节点上 = j/M **线性**（与网格无关）；膜阻 `R_film=1/kw` 与地层阻 `R_form=ln(r_e/r_w)/(2πk·dx)` **串联**→流体见 `g_eff=1/(R_film+R_form)`，稳态流体=换热递推 `r=mcd/(mcd+g_eff)`；壁节点分压 `T_w=T_e+(T_i−T_e)R_form/(R_film+R_form)`。内核+app ~190 行（`TwoFluid.h` ~48 / `twoFluid.cpp` 构造几何 ~50 + `formationStep` ~55 + `thermalStep` 接线 ~10 + loader 8 / app ~14）+ 1 case + 1 checker + 1 `add_test`，**132 判据**：稳态流体递推逐 cell 1e-9（20）、壁分压逐段 1e-9（20）、对数剖面 j 线性逐段 max 1e-9（20）、远场=T_e（20）、径向单调（20）、轴向单调（1）、稳态能量 qStore→0 / qWall=−qFar / qIntake=qFar / 流体焓损=远场排热 / `qFar=G_FORM·Σ(T_w−T_e)`（5）、**瞬态 intake=storage+far 逐帧精确**（20 帧）+ 首帧 storage>1W 有限热容见证（1）、五拒收（无 thermal / k≤0 / r_e≤r_w / rhoC≤0 / shells<1）。实测 T_out **318.463084 K**、r 0.986969、g_eff 0.346749、R_form 2.383928 K/W、qFar 303.0023 W；t=100 帧 300.575=78.253+222.322。零声明逐位不变由 `twoFluidThermalRelax/Exchanger/Geothermal`（Test #129-131）+ 全量回归 137/137 见证。`verify.twoFluidFormationRadial`（Test #132，0.31 s）。**边界**：逐段独立无轴向岩石导热；算子分裂 O(dt) 瞬态误差（稳态精确，判据锁稳态+逐帧地层自身能量平衡）；远场 Dirichlet，线源 E1 瞬态留待增量 2；rhoC 测试值 1e4（真实岩石 ~2e6）只验离散与守恒；formation 只进两流体温度方程。D3 保持 IN_PROGRESS |
