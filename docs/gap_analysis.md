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

> **落点编号说明（2026-11-11 审计）**：带「（预估）」的 §5.x 落点只是**占位编号**，实际按增量关闭顺序分配，不预留。架构文档当前已用到 **§5.64**（§5.63 = C1/A4，§5.64 = B6）；因此 A2 原写 §5.63、A3 原写 §5.64 均已被占用，C3 原写 §5.66 与 A6 的 §5.67 也不一致。这三条在本轮审计里改为「待分配 / 顺延」，其余预估编号维持不动，等该项真正启动时再定。

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
- **当前状态**：`OPEN`
- **阻塞项**：无；但该量在 1D 定截面积稳态管道里恒为 0，只能在**瞬态行**或**截面积突变节点**处被观测。
- **内核改动**：`src/1d/src/transient.cpp` 节点行（degree-generic §:395-:416）需要在 `Σ s_i A_i v_i = 0` 之外再带一笔 `Σ s_i β_i A_i v_i²`；`transient.cpp::pipeFrictionPa` 不动。
- **成本估计**：1 增量
- **落点**：architecture.md 待分配（原写 §5.64，已由 B6 占用）
- **被命名 5 次**：§5.34 边界① / §5.35 边界② / §5.36 边界② / §5.37 边界② / §5.38 边界③。
- **现状接受数**：第一层 χ 间隙 5–7%（ramp 类瞬态用例）。

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
- **当前状态**：`OPEN`（策略决定，不是 bug）
- **阻塞项**：需要一条文档化的策略
- **内核改动**：文档为主；`src/1d/src/network.cpp:342`（1e-13·dp 残差门）与 `:359`（1e-9 相对括号门）以及 `src/1d/src/transient.cpp:464`（1e-11 maxRel）当前都是硬编码（**2026-11-11 审计修正行号**：原文写 `:347` / `:364` / `:317`，分别漂移 5 / 5 / 147 行；网络 Newton 的外部门 `ctl.tol` 可调，但这三处内部门不受它控制）。
- **成本估计**：0.5 增量（文档 + 一个 policy assertion）
- **落点**：architecture.md §5.68（预估）
- **已测量**：§5.47 稳态括号门先于残差门触发，§5.49 瞬态牛顿出口不收敛到不动点而累积缓慢松弛。

### A8. 1D↔3D 界面多相 / 可压伙伴块

- **当前状态**：`OPEN`
- **阻塞项**：无
- **成本估计**：2 增量
- **落点**：architecture.md §5.69（预估）

### A9. 两条耦合边同处串联位置

- **当前状态**：`OPEN`
- **阻塞项**：无
- **成本估计**：1 增量
- **落点**：architecture.md §5.70（预估）

### A10. 行 `M·dq/dt` 双重计费已测未消

- **当前状态**：`OPEN`（已测：层 1 `rowError` 44–63%）
- **阻塞项**：需要准则"多少层才允许 `rowError` 低于 X"
- **成本估计**：1 增量
- **落点**：architecture.md §5.71（预估）
- **相关**：§5.38/§5.39 测得，§5.40 分解，§5.44 仍未消。

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
- **当前状态**：`OPEN`（**2026-11-11 审计收窄口径**；原标题把 `profile` 也算作"表"，与 §5.41 已落地的算例冲突）
- **口径**：`profile` 是 `starts.<arm>` 里的**字段**，不是独立表块。真正的第三张表是 §5.43 的 `convectiveToRow{ <臂>{ <边> charge|report } }`。本格要的是 `ddtFrom` + `inertiaFrom` + `convectiveToRow` **三张独立逐边表**连同逐边 `profile` 进**同一份 case**。
- **已落地但仍不够**：`tests/cases/transientCoupledMultiTurbulentProfileInertia/case.dict:114/125/134`（§5.41）已有 `ddtFrom` + `inertiaFrom` + 逐边 `profile`；`tests/cases/transientCoupledMultiTurbulentShapeAndInertiaSplit/case.dict:137`（§5.44）已有 `inertiaFrom` + `convectiveToRow` + `profile`。全仓**没有**一份 case 同时含三张表（审计核对了所有 `case*.dict`）。
- **内核改动**：无（三张表两两已交叉：§5.40 / §5.41 / §5.44）
- **成本估计**：1 增量
- **落点**：architecture.md §5.76（预估）

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

- **当前状态**：`OPEN`
- **内核改动**：与 §5.50 相同路径，只是 key 换成 HB
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.78（预估）

### C3. 稳态环空

- **当前状态**：`OPEN`（与 A6 同体）
- **内核改动**：`network.cpp` 新截面分支
- **成本估计**：1 增量
- **落点**：architecture.md §5.67（与 A6 合并；原写 §5.66，与 A6 条目不一致）

### C4. `reportEvery` 与 `endTime` 解耦

- **当前状态**：`OPEN`
- **内核改动**：parser；`transient.cpp` 报告路径
- **成本估计**：0.25 增量
- **落点**：architecture.md §5.79（预估）
- **备注**：§5.49 以来 1D 瞬态核心是固定 `dt`（`src/1d/src/transient.cpp:739` 由 `endTime/dt` 定步数、`:740-744` 每 `reportEvery` 步采样），隐藏了"求解器自己的 `dt` 能不能跑在打印步之下"的问题。**审计补充（2026-11-11）**：部分解耦其实已有——1D app `src/1d/apps/transientApp.cpp:100` 用 `|| i == nSteps` 强制末步打印，3D PISO `src/3d/apps/pisoFoamApp.cpp:148` 用 `dt = min(dt, endTime−t)` 裁剪末步（§5.34 已记）。这些不是 CFL 自适应，C4 主体仍 `OPEN`。

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
- **当前状态**：`OPEN`
- **成本估计**：1.5 增量
- **落点**：architecture.md §6.x
- **备注**：每个瞬态用例都用固定 `dt`；在最慢管上放一个 CFL 探针能直接消掉 §5.38 记录的 `dt 0.05` 发散墙。

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
6. **D5** — 自适应时间步。A1 阻塞期间**这是 D 类里钩子最小的第一项**；按 §2 的标准，零内核改动的 B 类项（B5 / C2 / C4）风险更低，应优先穿插，不要因为这条编号就跳过它们。（2026-11-11 审计澄清：原文"可执行的第一项"与 §2「首选从 B 类开始」相互冲突。）
7. ~~**B6**~~ — **CLOSED**（§5.64，本轮关闭）；**B5 / C2 / C3 / C4** — 按需要穿插。
8. **A2 / A3 / A5 / A6 / A7 / A8 / A9 / A10 / A11** — 按引用热度排。
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
| 2026-11-11 | 文档审计 | — | 待补 hash | — | **按文档排查代码完成度**：架构 §5.1–5.64（`###` 小节 63 个，另 5.17 是 `##` 标题）、审计时 CTest 117 项（本轮新增 B6 后 118 项）；30 条能力差异中审计当时 **6 条 CLOSED**（A4 / B1 / B2 / B3 / B4 / C1）、**2 条 BLOCKED**（A1 / A11）、其余 OPEN —— 本轮关闭 B6 后为 7 条 CLOSED、21 条 OPEN。7 处文档漂移：B6 孤立草稿虚假声明 + Jacobian 注释反了、B5 标题口径与 §5.41 冲突、A7 三处行号漂移 5/5/147 行、A2 落点 §5.63 与 C1/A4 撞号、A3 落点 §5.64 预留、C3 与 A6 落点不一致、C4 漏记 1D/3D 已有的末步解耦；优先级第 6 项「D5 第一项」与 §2「B 类优先」互斥。B5 收窄为「三张独立逐边表 + 逐边 profile」（全仓无此 case），仍 OPEN |
| 2026-11-11 | B6 | §5.64 | 待补 hash | 714.41 s / 118 | 全网络全死支的 0 迭代不动点：0 内核改动，2 份入库 dict + 4 个临时派生臂、45 判据；E 随定压均值 25000→25500、阈值刀刃 Q 仍为 0、单管放开时 iterations 仍为 0、注入臂 Newton 19 步且 E = 33978.9475076 Pa 与 Buckingham 闭式反解一致；`verify.steadyAllBelowThreshold`（Test #45，全量里 0.88 s）、全量 **118/118 绿** |
