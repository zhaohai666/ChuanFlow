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

## 1. A 类 —— 已命名且已记录的 open cell

> 这些差异在之前的 §5.x 增量里**已经被发现**、**已经被量化**、**已经被写进 architecture.md 的"边界①/②/…"小节**；它们不是新发现，只是尚未关闭。除非阻塞项解除，不要重新从零推导。

### A1. 屈服 + 湍流摩擦闭合（Darby–Metzner / Wilson–Thomas 类）

- **参考工具**：Pipesim / OLGA 非牛顿湍流管段默认行为
- **当前状态**：`BLOCKED`
- **阻塞项**：Wilson & Thomas 2006 (*Can. J. Chem. Eng.* 84:520–525) 与 Thomas & Wilson 2007 (Hydrotransport 17) 两份文献对 `β = τ_w/τ_B` 的区间收缩写法与 `11.6/χ` 钝化项的陈述不一致； adjudicator 已推得 `D U* ρ/μ_B = He^0.5·β^0.5`（对数项必须在 `τ_B → 0` 时塌回 `Re*`），但 β 自身的解方程形式仍待第三方核对。Task #84 保留已失败路径与两份已提取文本（`/tmp/wilson_thomas.pdf`、`/tmp/wilson_rough.pdf.txt`，本箱无 poppler，可用 stdlib zlib 重抽）。
- **内核改动**：`src/1d/src/network.cpp` 增加 `yieldedTurbulentFriction(Re, He, n)`；`src/1d/src/transient.cpp::pipeFrictionPa` 同路径；§5.51 的 `yieldedTransitionOf` 已报告 locus，此处只需写摩擦侧。
- **成本估计**：1 增量（设计已 settled，只差转录核对）
- **落点**：architecture.md §5.61
- **判据**：§5.51 的 `f·Re_B·g(xi) = 64` 在 4× 与 151× 超过 locus 的两格**必须由新公式复现**（替换而非删除）；"rate is the laminar root" 家族（最差 3.2e-12）同样替换。

### A2. 亚阈值蠕变 / 凝胶强度 / 重启压力

- **参考工具**：OLGA restart-pressure / Pipesim gel-strength
- **当前状态**：`OPEN`
- **阻塞项**：无（需要一条本构律，不是文献问题）
- **内核改动**：`src/1d/src/transient.cpp::holdsStatic` 当前是 `qPrev == 0 && |drive| <= threshold` 的刚性判据；亚阈值蠕变需要把它改成 `qPrev <= q_creep(t_rest, …)` 的时变判据，并在节点行加一笔凝胶累积项。
- **成本估计**：1.5 增量（本构 + 一个启停算例）
- **落点**：architecture.md §5.62（预估）
- **备注**：当前内核 §5.49 分支 `Q` 严格为 0、无蠕变、无触变性、无"静止 N 小时后重启需要的压力"（OLGA/Pipesim restart-pressure 类）。

### A3. β（对流动量项）与静止支同行

- **参考工具**：Fluent 截面动量通量因子
- **当前状态**：`OPEN`
- **阻塞项**：无；但该量在 1D 定截面积稳态管道里恒为 0，只能在**瞬态行**或**截面积突变节点**处被观测。
- **内核改动**：`src/1d/src/transient.cpp` 节点行（degree-generic §:395-:416）需要在 `Σ s_i A_i v_i = 0` 之外再带一笔 `Σ s_i β_i A_i v_i²`；`transient.cpp::pipeFrictionPa` 不动。
- **成本估计**：1 增量
- **落点**：architecture.md §5.63（预估）
- **被命名 5 次**：§5.34 边界① / §5.35 边界② / §5.36 边界② / §5.37 边界② / §5.38 边界③。
- **现状接受数**：第一层 χ 间隙 5–7%（ramp 类瞬态用例）。

### A4. 耦合边拒绝流变与屈服

- **参考工具**：Fluent system coupling + 非牛顿伙伴块
- **当前状态**：`OPEN`
- **阻塞项**：无
- **内核改动**：`src/1d/src/transient.cpp:226–229` 当前显式拒绝非牛顿 key；耦合调度 `src/coupling/coupledNetwork.cpp` 需要把 `pipeFrictionPa` 用完整 rheology key 调用。
- **成本估计**：1 增量
- **落点**：architecture.md §5.64（预估）
- **备注**：与 C1 本质相同；A 类与 C 类的区别是 A 在架构文档里有命名边界、C 没有。

### A5. PVT / 黑油 / 两流体分布式拒绝流变与屈服

- **参考工具**：OLGA + Pipesim 全管线非牛顿
- **当前状态**：`OPEN`
- **阻塞项**：无（表面比 A4 大）
- **内核改动**：`src/1d/src/twoFluid.cpp` + `src/pvt/*` 调度面
- **成本估计**：2 增量
- **落点**：architecture.md §5.65（预估）

### A6. 环空 / 矩形 duct（Dodge–Metzner）

- **参考工具**：Pipesim annulus / OLGA eccentric annulus
- **当前状态**：`OPEN`
- **阻塞项**：无
- **内核改动**：`src/1d/src/network.cpp` 增加一个截面分支 `annulus{ dInner … }`；层流锚点取 Fredrickson–Bird 1958 闭式；`Re_GP` 用 Dodge–Metzner。
- **成本估计**：1 增量
- **落点**：architecture.md §5.66（预估）
- **备注**：README 已写"圆管 only（5.45）"。

### A7. 求解器退出门槛已测未立

- **参考工具**：所有商业软件都有用户可调收敛门槛
- **当前状态**：`OPEN`（策略决定，不是 bug）
- **阻塞项**：需要一条文档化的策略
- **内核改动**：文档为主；`network.cpp:347`（1e-13·dp 残差门）与 `:364`（1e-9 相对括号门）以及 `transient.cpp:317`（1e-11 maxRel）当前都是硬编码。
- **成本估计**：0.5 增量（文档 + 一个 policy assertion）
- **落点**：architecture.md §5.67（预估）
- **已测量**：§5.47 稳态括号门先于残差门触发，§5.49 瞬态牛顿出口不收敛到不动点而累积缓慢松弛。

### A8. 1D↔3D 界面多相 / 可压伙伴块

- **当前状态**：`OPEN`
- **阻塞项**：无
- **成本估计**：2 增量
- **落点**：architecture.md §5.68（预估）

### A9. 两条耦合边同处串联位置

- **当前状态**：`OPEN`
- **阻塞项**：无
- **成本估计**：1 增量
- **落点**：architecture.md §5.69（预估）

### A10. 行 `M·dq/dt` 双重计费已测未消

- **当前状态**：`OPEN`（已测：层 1 `rowError` 44–63%）
- **阻塞项**：需要准则"多少层才允许 `rowError` 低于 X"
- **成本估计**：1 增量
- **落点**：architecture.md §5.70（预估）
- **相关**：§5.38/§5.39 测得，§5.40 分解，§5.44 仍未消。

### A11. 真 MPI 传输 / 域耦合矩阵求解

- **当前状态**：`BLOCKED`
- **阻塞项**：用户 OK 之前不得安装 MPI 工具包（ standing 约束）
- **成本估计**：≥ 3 增量
- **落点**：architecture.md §5.71（预估）

## 2. B 类 —— 零内核改动增量（case + checker only）

> 与 §5.50 / §5.52 / §5.60 同一形态：内核已经能回答这个问题，只是没有一个新的拓扑或新的参数交叉去把它问出来。**首选从这一类开始**，因为风险最低、回归零成本。

### B1. 源 / 汇日程唤醒停泊臂

- **参考工具**：Pipesim scheduled injection / OLGA `MASSFLOW` 边界
- **当前状态**：`CLOSED`（由 §5.61 / task #95 关闭）
- **内核改动**：`src/1d/include/cf/1d/Transient.h::TransientSource` 加 ramp 字段；`src/1d/src/transient.cpp` 加 `sourceRateAt()` helper，`step()` 与 `commitStep()` 把 `nd.source` 改成按 `tNew` 求值；parser 加 `schedule{}` 子块。
- **成本估计**：1 增量（实质是 C 类小钩子，与阀日程同形态）
- **落点**：architecture.md §5.72（预估）
- **来源**：§5.56 边界④。
- **拓扑**：一根屈服壁 + 一根带 `valveClose/Open` 的泄放管 + 一个 `source{ schedule … }` 在同一自由节点；日程在 sealed 窗口内把节点抬过该壁的放开压力，让"本来一直 held"的壁被叫醒一次。

### B2. 节点初值 +600 Pa，使 n=0.5 臂被第一时间问

- **参考工具**：OLGA `INITIAL` 扫描
- **当前状态**：`OPEN`
- **内核改动**：无
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.73（预估）
- **来源**：§5.59 边界③ / §5.60 边界④。
- **做法**：把 `Fa..Fd` 的 `pressure 2.7e4` 改成 `2.76e4`（释放压力 23000/22000/21000/20000，+600 后 n=0.5 臂 c、d 的起始 drive 都跨过自家阈值）。

### B3. HB × 阀日程 × 双向闩锁 三路同台

- **参考工具**：OLGA 带 Herschel-Bulkley 泥浆的阀行程
- **当前状态**：`CLOSED`（由 §5.59 / task #93 关闭）
- **内核改动**：无（§5.52 是 Bingham，§5.54 也是 Bingham；把 `flowIndex` 声明进同一份 dict 即成新算例）
- **成本估计**：1 增量
- **落点**：architecture.md §5.59
- **判据**：§5.52 的"阶梯只动翻步不动是否"必须仍在 HB 下成立；§5.54 的"谁能重睡由邻居阈值次序决定"必须仍在 HB 下成立。

### B4. 分支节点 + 双向闩锁（一节点 ≥ 三根管）

- **参考工具**：OLGA Tee + check-valve 重关
- **当前状态**：`OPEN`
- **内核改动**：无（§5.54 给每个自由节点**只**一进一出；§5.55 的分流节点已经在 5.55 里被实现但只出现在 HB 停泊带，没有与闩锁同台）
- **成本估计**：1 增量
- **落点**：architecture.md §5.75（预估）
- **关键**：重闩必须被"被 hold 的那根管自己的流回"触发，而不只是邻居的日程。

### B5. 逐边 `ddtFrom` × 逐边 `inertiaFrom` × 逐边 `profile` 三张表同台

- **参考工具**：Fluent system coupling 全特性
- **当前状态**：`OPEN`
- **内核改动**：无（三张表两两已交叉：§5.40 / §5.41 / §5.44）
- **成本估计**：1 增量
- **落点**：architecture.md §5.76（预估）

### B6. 全网络每条边都低于自家阈值（稳态）

- **参考工具**：Pipesim gelled pipeline 启动研究
- **当前状态**：`OPEN`（§5.49 边界⑧ 记录"HP 帽让 Jacobian 非奇异"，未立判据）
- **内核改动**：无
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.77（预估）

## 3. C 类 —— 小内核钩子 + 算例（与 §5.48 / §5.51 同形态）

> 与 B 类的差别：需要在某处加一条调度线或新增一个 key；与 D 类的差别：不动方程、不新增模块。

### C1. 流变参数（`flowIndex` / `consistency` / `yieldStress`）进入耦合边

- **当前状态**：`OPEN`（与 A4 同体）
- **内核改动**：`coupledNetwork.cpp` 调度 `pipeFrictionPa` 时带完整 rheology key；`transient.cpp:226–229` 的拒绝逻辑去掉。
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.64（与 A4 合并）

### C2. HB × 阀 `K`（无日程）

- **当前状态**：`OPEN`
- **内核改动**：与 §5.50 相同路径，只是 key 换成 HB
- **成本估计**：0.5 增量
- **落点**：architecture.md §5.78（预估）

### C3. 稳态环空

- **当前状态**：`OPEN`（与 A6 同体）
- **内核改动**：`network.cpp` 新截面分支
- **成本估计**：1 增量
- **落点**：architecture.md §5.66（与 A6 合并）

### C4. `reportEvery` 与 `endTime` 解耦

- **当前状态**：`OPEN`
- **内核改动**：parser；`transient.cpp` 报告路径
- **成本估计**：0.25 增量
- **落点**：architecture.md §5.79（预估）
- **备注**：§5.49 以来每个瞬态用例都是 `endTime/dt = reportEvery × states`，隐藏了"求解器自己的 `dt` 能不能跑在打印步之下"的问题。

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
2. **B4** — 分支节点 + 双向闩锁。B3 之后闩锁线的下一个自然跳。
3. **B2** — 节点初值扫描（+600 Pa 使 n=0.5 臂被问）。0.5 增量。
4. **C1** — 流变进入耦合边（与 A4 合并关闭）。
5. **A1** —— 任务 #84，β 项被第三方文献确认后启动。
6. **D5** — 自适应时间步。D 类里唯一钩子较小的一项。
7. **B5 / B6 / C2 / C3 / C4** — 按需要穿插。
8. **A2 / A3 / A5 / A6 / A7 / A8 / A9 / A10 / A11** — 按引用热度排。
9. **D1 / D2 / D3 / D4 / D6 / D7 / D8 / D9** — 大型能力块，按产品规划排期。

## 6. 实施日志

> 每完成一个增量，按 `| 日期 | ID | 增量 | commit | 实测 ctest | 备注 |` 追加一行。

| 日期 | ID | 增量 | commit | 实测 ctest | 备注 |
|------|----|------|--------|------------|------|
| 2026-11-11 | — | 文档初版 | 93b2c69 | — | A 类 11 条、B 类 6 条、C 类 4 条、D 类 9 条 |
| 2026-11-11 | B3 | — | — | — | **审计为 CLOSED**：§5.59（task #93, commit 已在闩锁线）= Herschel-Bulkley × 阀 `valveClose`/`valveOpen` 日程 × 双向闩锁三路同台，三档 (n,K)=(0.5,1.0)/(1.0,0.1)/(2.0,0.02)；本节 5.59 起即关闭，文档初版误记为 OPEN |
| 2026-11-11 | B1 | §5.61 | 6581c96 | 1311.16 s / 116 | 源日程唤醒停泊臂：`sources{ schedule{ openRate; close; open } }` 与阀日程同一 ramp helper；1 HB 壁 + 1 阀泄放 + 1 注入源同台；6423 判据、2 case 文件、Test #90 (2.23 s) |
