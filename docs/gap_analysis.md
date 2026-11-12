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

> **落点编号说明（2026-11-11 审计）**：带「（预估）」的 §5.x 落点只是**占位编号**，实际按增量关闭顺序分配，不预留。架构文档当前已用到 **§5.72**（§5.63 = C1/A4，§5.64 = B6，§5.65 = C2，§5.66 = C4，§5.67 = B5，§5.68 = D5，§5.69 = A7，§5.70 = A3，§5.71 = A9，§5.72 = A10）；因此 A2 原写 §5.63 已被占用，C3 原写 §5.66 与 A6 的 §5.67 不一致。这几条在本轮审计里改为「待分配 / 顺延」，其余预估编号维持不动，等该项真正启动时再定。

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
- **当前状态**：`OPEN`
- **阻塞项**：无（需要一条本构律，不是文献问题）
- **内核改动**：`src/1d/src/transient.cpp::holdsStatic` 当前是 `qPrev == 0 && |drive| <= threshold` 的刚性判据；亚阈值蠕变需要把它改成 `qPrev <= q_creep(t_rest, …)` 的时变判据，并在节点行加一笔凝胶累积项。
- **成本估计**：1.5 增量（本构 + 一个启停算例）
- **落点**：architecture.md 待分配（原写 §5.63，已由 C1/A4 占用）
- **备注**：当前内核 §5.49 分支 `Q` 严格为 0、无蠕变、无触变性、无"静止 N 小时后重启需要的压力"（OLGA/Pipesim restart-pressure 类）。§5.64 的注入臂把这条缺席**量化**了：全死网络一旦有注入，只能靠把自由节点抬到 33978.9 Pa 让一根管越过阈值开路，而不是靠任何亚阈值流动。

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
- **当前状态**：`OPEN`
- **阻塞项**：无（表面比 A4 大）
- **内核改动**：`src/1d/src/twoFluid.cpp` + `src/pvt/*` 调度面
- **成本估计**：2 增量
- **落点**：architecture.md §5.66（预估）
- **§5.63 后的范围收缩**：耦合边那一路的拒收已由 C1/A4 解除；本条现在只剩 PVT 表 / 黑油 / 两流体分布式求解那三行，`steadyNetwork` 里的 `hasBlackOil || hasPvtFluid` 拒收仍在。

### A6. 环空 / 矩形 duct（Dodge–Metzner）

- **参考工具**：Pipesim annulus / OLGA eccentric annulus
- **当前状态**：`OPEN`
- **阻塞项**：无
- **内核改动**：`src/1d/src/network.cpp` 增加一个截面分支 `annulus{ dInner … }`；层流锚点取 Fredrickson–Bird 1958 闭式；`Re_GP` 用 Dodge–Metzner。
- **成本估计**：1 增量
- **落点**：architecture.md §5.67（预估）
- **备注**：README 已写"圆管 only（5.45）"。

### A7. 求解器退出门槛已测未立

- **参考工具**：所有商业软件都有用户可调收敛门槛
- **当前状态**：`CLOSED`（2026-11-11，§5.69）
- **阻塞项**：无
- **内核改动**：`src/1d/include/cf/1d/SolverGates.h`（新建，`cf::gates` 五个 `inline constexpr`：steadyPipeResidualGate 1e-13 / steadyPipeBracketGate 1e-9 / steadyPipeMaxIter 100 / transientNewtonGate 1e-11 / transientNewtonMaxIter 60）；三处站点改用常量名，零行为变化（`network.cpp` 残差门 `:343`、括号门 `:360`、`binghamWallStress` 帽 `:145` 也并入同一管解帽政策；`transient.cpp` 门 `:465`、帽 `:365`）。两相族（`solvePipeFlowTwoPhase`，残差门 1e-7/1e-12·pUp）按声明范围留在 header 之外。
- **成本估计**：0.5 增量（文档 + 一个 policy assertion）；实测 ~0.5
- **落点**：architecture.md §5.69
- **已测量**：§5.47 稳态括号门先于残差门触发，§5.49 瞬态牛顿出口不收敛到不动点而累积缓慢松弛；checker 实测六管闭式足迹 2.4e-13~4.9e-12（括号门名义 1e-9 的实际余量）、瞬态 march 漂移非零 < 1e-9（worst 1.46e-12）、tol 1e-16 逐位不变、maxIter 0 大声失败。

### A8. 1D↔3D 界面多相 / 可压伙伴块

- **当前状态**：`OPEN`
- **阻塞项**：无
- **成本估计**：2 增量
- **落点**：architecture.md §5.69（预估）

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
- **新边界**：`turbulence{}` 与非牛顿流体互斥（两者争同一个逐面黏度钩子）；多边同时带流变已转发但未端到端验证；PVT/黑油/两流体仍拒收（A5）

### C2. HB × 阀 `K`（无日程）

- **当前状态**：`CLOSED`（由 §5.65 关闭）
- **内核改动**：**无**。parser 自 §5.47 / §5.49 / §5.50 起就同时收一根管上的 `flowIndex` / `consistency` / `yieldStress` / `valve`，没有 guard 拒收这个组合——只是从未被问过。成本估计 0.5 → 实际 ~0.3 增量（2 case 文件 + 1 checker + 1 `add_test`）
- **落点**：architecture.md §5.65
- **实测**：56 判据 / 0 失败，16 管 × 2001 状态全部被「HB 剖面积分 + `K ρ|Q|Q/2A²`」递推复现（最差 4.347e-10）；持停管五个数量级 K 逐 token 相同；同一 `K = 500` 在 `n = 0.5` 管上拿走 2.33 Pa、在 `n = 1.5` 管上拿走 308.76 Pa（133×）——**阀的分成由壁律指数决定**；三条停止迹（第 9 / 9 / 5 步）逐状态复现；稳态核对与行进行在同一个两损失根上会合。`verify.transientYieldIndexValveK`
- **边界**：阀 K 仍按几何项计，不随流变变稠变稀（阀有自己的 K，流体有自己的 K_cons）；带日程的组合由 §5.58–5.62 覆盖

### C3. 稳态环空

- **当前状态**：`OPEN`（与 A6 同体）
- **内核改动**：`network.cpp` 新截面分支
- **成本估计**：1 增量
- **落点**：architecture.md §5.67（与 A6 合并；原写 §5.66，与 A6 条目不一致）

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
- **当前状态**：`OPEN`
- **成本估计**：3 增量（ESP 多级曲线 + 气体干涉；气举阀动力学；有杆泵井下示功图）
- **落点**：architecture.md §6.x
- **README 现状**：只有离心泵。

### D2. 蜡 / 沥青质 / 结垢沉积

- **参考工具**：OLGA / Pipesim wax & scale 模块
- **当前状态**：`OPEN`
- **成本估计**：3 增量
- **落点**：architecture.md §6.x
- **依赖**：需要在 `src/pvt` 里增加固相质量平衡，并与热力场耦合。

### D3. 瞬态地层热耦合（径向地层网格）

- **参考工具**：OLGA `HEATFORM`
- **当前状态**：`OPEN`
- **成本估计**：2 增量
- **落点**：architecture.md §6.x
- **README 现状**：只有"恒壁温或地温线"，无径向地层离散。

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
7. ~~**B6**~~ — **CLOSED**（§5.64）；~~**C2**~~ — **CLOSED**（§5.65，HB × 无日程阀 K，阀的分成由壁律指数决定）；~~**C4**~~ — **CLOSED**（§5.66，ceil + 末步裁剪 + 末态恒打印，endTime 为 dt 整倍数时逐位兼容）；~~**B5**~~ — **CLOSED**（§5.67，三张逐边表 × 逐边 profile 同台，quasi+slab+charge 三表叠一边）；~~**A7**~~ — **CLOSED**（§5.69，`SolverGates.h` 五常量单一来源 + policy assertion，零行为变化）；**C3** — 按需要穿插。
8. ~~**A9**~~ — **CLOSED**（§5.71，两条耦合边串联：一个速率逐位绑定、行可加、账按 M1/(M1+M2) 分配、斜率阶梯反比于 M 和；实测 ~0.8 增量而非 1）。~~**A10**~~ — **CLOSED**（§5.72，双重计费准则成文：首层按档阶梯界 + 塌缩日程 + M_dyn 排他界，`slab` 成为推荐选项；实测 ~0.8 增量而非 1）。**A2 / A5 / A6 / A8 / A11** — 按引用热度排。
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
