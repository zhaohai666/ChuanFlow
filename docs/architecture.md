# ChuanFlow 架构说明（当前实现）

> 本文档描述仓库中已经落地的部分：`src/common` 基础层、`src/pvt` 统一物性层
> （PR 状态方程 / 多组分闪蒸 / PVT 表）、`src/1d` 一维管网内核（M1 稳态 +
> M2 首增量：瞬态水力 / 持液率输运 / 水合物动力学）与 `src/3d` 三维 FVM
> 内核，以及 `src/bindings` Python 绑定与 `cli/` Rust 入口。
> README 中的 1D↔3D 耦合 / GUI 等模块属于路线图（M4–M6），尚未实现。
> 物性层与一维/三维内核的耦合（密度由 pvt 实时求取而非用户输入常数）也已
> 列入 M2/M3 后续增量。
>
> ChuanFlow 不是 OpenFOAM 的 Fork：以下设计参考了 OpenFOAM-dev 的工程思想
> （网格拓扑、场分层、fvc/fvm 分离、SIMPLE 算法回路），但所有数据结构、离散
> 格式与求解器均为独立自研实现。

## 1. 总体分层

```
apps (chflow-1d-steady / chflow-1d-transient / chflow-pvt / chflow-scalar-transport / chflow-simple-foam)
  ├── src/1d          一维管网内核（M1 稳态节点法 + M2 首增量：瞬态/持液率/水合物）
  │      └── src/common  基础层（字典解析、日志/异常）
  ├── src/pvt         统一物性层（PR 状态方程 / 多组分 TP 闪蒸 / PVT 表插值）
  │      └── src/common  基础层（字典解析、日志/异常）
  ├── src/bindings    pybind11 模块 chflow_core（Python 包 python/chflow）
  │      └── src/1d + src/pvt
  └── src/3d          三维 FVM 内核（mesh / fields / discrete / solvers / io）
        └── src/common  基础层（向量、Field 容器、稀疏矩阵+Krylov求解、字典解析、日志/异常）
cli/  Rust 自动化入口 chuanflow（run steady|transient|pvt <case.dict>，定位并转发引擎二进制）
```

内核零第三方依赖：纯 C++20 + 标准库，CMake ≥ 3.24 构建，测试用 CTest +
Python 标准库脚本校验（无外部测试框架）。

## 2. src/common 基础层

| 头文件 | 作用 | 对应 OpenFOAM 概念 |
|---|---|---|
| `Labels.h` | `label`(int32) / `scalar`(double) 别名 | `label` / `scalar` |
| `Vector.h` | 3 维向量，`dot/mag/cross` | `vector` |
| `Field.h` | `Field<T>` = `std::vector<T>` 薄封装 | `Field<Type>` |
| `Matrix.h/.cpp` | CSR 稀疏矩阵 + 求解器 | `lduMatrix` + solvers |
| `Dictionary.h/.cpp` | OpenFOAM 风格字典解析（嵌套 block、`key value;`） | `dictionary` / `entry` |
| `Log.h` / `Error.h` | 分级日志与 `cfFatal` 异常 | `Pout` / `FatalError` |

### 线性求解器
`SparseMatrix::solveBiCGStab`：van der Vorst (1992) 右预处理 BiCGStab，Jacobi
预处理。关键约定（改动时务必保持）：

- 搜索方向 `p` 跨迭代必须保持**未预处理**，`z = M⁻¹p` 只是瞬态量；
- 第二个半步同样要过预处理：`sh = M⁻¹s`、`t = A·sh`、
  `ω = (t,s)/(t,t)`、`x += ω·sh`、`r = s − ω·t`；
- `infNorm` 对 NaN 返回 NaN（否则全 NaN 残差会被误报为 0 并“收敛”）。

## 3. src/pvt 统一物性层

README 把 `pvt` 列为"整个项目最重要的共享模块"：一维、三维、耦合仿真必须
共用同一套相平衡与密度。首版落地 PR 状态方程、多组分等温闪蒸与 PVT 表插值。

### 3.1 组分与混合物 `Components.h` / `components.cpp`
- 内置 10 个天然气/凝析油组分（甲烷、乙烷、丙烷、异丁烷、正丁烷、正戊烷、
  正己烷、CO₂、H₂S、N₂）的 `M / Tc / Pc / ω` 表值，工厂函数构造时统一换算
  g/mol → kg/mol、bar → Pa；
- `Mixture`：组分 + 摩尔分数，`setFractions` 校验非负并归一，另有
  `molarMass()`、`massFraction(i)`、`describe()`；
- `componentByName`（返回 found 标志）与 `requireComponent`（字典里写错组分
  名直接 `cfFatal` 并列出可用名）。

### 3.2 Peng-Robinson 状态方程 `PengRobinson.h` / `pengRobinson.cpp`
自研实现，全部来自公开定义式：

- `a_i = 0.45724 (R·Tc)²/Pc · α(T)`、`b_i = 0.07780 R·Tc/Pc`——注意是
  **R·Tc 而不是 R·T**：写成 R·T 会让 `a/b` 量级崩塌，得到 Z≫1 的非物理根。
  这条由两个不变量测试钉住：临界点压缩因子 `Z_c = 0.30740`（±2%）和
  低压理想气体极限 `ρ → pM/RT`；
- van der Waals 单流体混合规则（`k_ij = 0`，接口在 `aIJ` 上预留）；
- 压缩立方程 `Z³ − (1−B)Z² + (A−3B²−2B)Z − (AB−B²−B³) = 0` 用 Cardano
  （单实根）/ 三角法（三实根）求全部实根，临界点三重根处有 `p ≥ 0` 保护分支；
- `RootPreference`：`LowestGibbs`（稳定根，最小化 `Σx_i ln φ_i`）/
  `Liquid`（大于 B 的最小根）/ `Vapor`（最大根）；`ln φ_i` 解析导出；
- `pureVaporPressure`：纯组分饱和压由逸度相等 `ln φ_L = ln φ_V` 求得，
  Wilson 式起步 + 带 ±0.5p 信任域的 Newton，`T ≥ Tc` 返回 0。

### 3.3 多组分 TP 闪蒸 `Flash.h` / `flash.cpp`
Wilson 初值 `K = (Pc/p)·exp[5.373(1+ω)(1−Tc/T)]` ⇒ Rachford-Rice ⇒ PR 逸度
系数 ⇒ 逐次代换 `K_i ← exp(ln φ_L,i − ln φ_V,i)`，收敛判据
`max|Δln K| < tol`（默认 1e-12）：

- Rachford-Rice 用 Newton + 二分括号保护（`f(β)` 在 (0,1) 上严格递减）；
- β 钳位：`f(0) ≤ 0` ⇒ 全液（β=0，气相为初生相），`f(1) ≥ 0` ⇒ 全气；
- **相态判定只看收敛后的 K**，不能沿用上一次迭代的钳位标志；同时要识别
  平凡解（`max|K_i − 1|` 过小）并要求 `x ≠ y`，否则单相区里 SS 收敛到
  `K = 1` 时残差只是数值噪声，会把两相态误标成单相；
- 单相时相标签由 PR 自身密度给出：`v/b < Z_c/B_c ≈ 3.951` 即液相；
- `bubblePressure`：纯组分直接退化为 `pureVaporPressure`（单元/绑定测试都
  断言两者在 1e-5 内一致）；混合物用"闪蒸是否找到两相"这一单调判据做二分，
  从 Wilson-Raoult 估计压力先向下走 ×0.6 找两相、再向上 ×1.25 找单相，
  然后 80 次二分收敛到 1e-10。这比 `ΣKz − 1` 残差稳健——后者在不存在液相根
  的压力区间会因平凡解产生伪根（曾把泡点误报为 20.8 MPa，实际包络在
  8.3 MPa 塌陷）。

### 3.4 PVT 表插值 `PvtTable.h` / `pvtTable.cpp`
- `BilinearTable`：要求 (p, T) 网格严格递增，采样在网格外夹边；
- `liquidDensityTable / vapourDensityTable / vapourFractionTable`：每个网格
  点跑的是同一个闪蒸，因此表路径与 EOS 直接求值路径热力学一致（回归测试在
  网格节点上比到 1e-14，半格插值误差 < 2%）。求解器每个单元/每步都要物性时
  可以预先制表，只付双线性查表代价。

### 3.5 算例字典与应用 `PvtCase.h` / `apps/pvtApp.cpp`
语法：`temperature`；`pressures (…)`（或单个 `pressure`）；
`composition { methane 0.34; … }`；`solve { tol; maxIter; bubblePoint; }`；
可选 `table { property liquidDensity|vapourDensity|vapourFraction;
pressures (…); temperatures (…); }`。
`chflow-pvt <case.dict>` 输出 `pvt flash T … p … beta … single … iterations …
converged … residual …`、`pvt liquid|vapour rho … M …`、逐组分
`pvt x/y <name> <value>`、`pvt bubble T … p …` 与 `pvt table …`
（12 位有效数字，供 `tests/check_pvt.py` 解析）。

尚未实现（pvt 后续增量）：水合物相平衡/生成条件（HH-2、kJ-i 类）、黑油/
挥发油分离器、黏度与界面张力模型、非等温（PH/PT 焓耦合）闪蒸，以及把
pvt 密度直供 1D/3D 求解器（当前一维的 ρ/μ 仍由 case.dict 给定常数）。

## 4. src/1d 一维管网内核

### 4.1 稳态节点法（M1）

单相节点法（node model）管网水力计算：未知量是各节点压力，每条管道给出
流量—压降本构，整体用阻尼 Newton 求解。

- **管道本构**：Darcy-Weisbach 压降 + 重位水头 `rho*g*dz`；摩阻系数用
  Churchill (1977) 全流态显式关联式（层流极限恰为 64/Re，湍流自动过渡），
  对给定 |dp| 以层流闭式解为初值做压缩定点迭代求 |v|（层流导数 ≤ 1/2，
  湍流超线性收敛）；
- **方程组**：自由节点质量守恒 `Σ Q_in − Σ Q_out − q_source = 0`，雅可比用
  中心差分的电导 `dQ/d(dp)` 解析装配（并以上/下限约束到层流电导
  `πD⁴/(128μL)` 防止 dp→0 时退化）；小规模稠密 LU + 部分主元求解；
  残差增长超 1.5 倍时对 Newton 步做对分线搜索。要求至少一个定压节点作参考；
- **算例字典**（`tests/cases/steadyNetwork/case.dict`）：
  `fluid{rho;mu}` / `nodes{ N { pressure …; } }` /
  `pipes{ P { from; to; length; diameter; roughness; elevation } }` /
  `sources{ s { node; rate } }` / `solve{ maxIter; tol; relax }`；
- **应用**：`chflow-1d-steady <case.dict>`，stdout 输出
  `node <name> p <value>` / `pipe <name> Q … v … Re … f … dp …`（12 位有效
  数字，供 `tests/check_steady_network.py` 解析），未收敛返回非零码。

M1 剩余的井筒/热力耦合与环路迭代（对标完整 Pipesim）尚未实现；一维的
`rho/mu` 目前仍由 case.dict 直接给定，改由 §3 的 pvt 层按 (p, T) 求取属于
后续增量。

### 4.2 瞬态水力 + 持液率输运 + 水合物动力学（M2 首增量）

`src/1d/include/cf/1d/Transient.h` + `src/1d/src/transient.cpp`。未知量向量
`x = [Q_1..Q_P, p_free]`（管流 + 自由节点压力）**同时** Newton 求解，
隐式 Euler：

- **管道动量**：`I·(Q − Q_n)/dt = p_from − p_to − ρg·dz − F(Q)`，
  惯性 `I = ρL/A`；摩阻 `F(Q)` 直接调用稳态内核同一套 Churchill 关联式，
  并做成 Q 的奇函数（反向流同号阻力）。雅可比对角 `I/dt + F'(Q)`，
  `F'` 用中心差分并下限约束到层流斜率的 0.2 倍，避免 Q→0 时雅可比退化；
- **节点压缩性**：`(p − p_n)·V/(K·dt) = Σ q_in + q_src`，`K` 为混合物体积
  模量；`K = 0` 退化为不可压端口，该行改为代数约束 `Σ q = 0`。定压节点是
  无穷 Dirichlet 水库，不参与体积累积；
- **收敛**：≤60 次 Newton，逐行按自身量级缩放残差（分母 +1 防零），
  `maxRel < 1e-11`；
- **持液率（drift-flux 段塞）输运**：对旧 α 做快照的显式上游差分
  `α_i^{n+1} = α_i + (Q·dt/V)(α_{i−1} − α_i)`，即一串全混釜（CSTR）级联——
  这个递推本身就是可复现的参考解，测试按 1e-9 对比；
- **水合物动力学**：`h = V_h/(V_L + V_h)`，`dh/dt = k(1 − h)`
  ⇒ `dV_L/dt = −k·V_L`，用积分因子精确积分
  `V_L^{n+1} = V_L·e^{−k·dt}`、`V_h = V_total − V_L`（体积守恒的
  液→水合物转移），因此**任意 dt 都无离散误差**；
- **应用** `chflow-1d-transient <case.dict>`：输出
  `t …` / `node <n> p … alpha … hydrate …` / `pipe <n> Q …`（12 位有效数字）。

字典语法见 `Transient.h` 末尾注释（`nodes{}` 增
`volume/bulkModulus/fixedPressure/liquidFraction/hydrateFraction/hydrateRate`，
`pipes{}` 增 `initialFlow`，新增 `time{ dt; endTime; reportEvery; }`）。

三个已验证的物理回归（`tests/check_transient.py`）：

| 算例 | 校验内容 |
|---|---|
| `transientBlowdown` | 罐排气：同时刻隐式 Euler 2×2 闭式递推（1e-9）+ 连续包络 `exp(−t/RC)`（2%）|
| `transientRlc` | 两柔性节点 + 惯性管 = RLC 欠阻尼振荡：首个流量峰的时间与幅值（3%）、≥4 次符号翻转、封闭系统液量守恒 |
| `transientSlug` | 段塞推进：串联管层流流量（1e-3）+ CSTR 级联递推（1e-9）|

尚未实现（M2 后续增量）：真两相动量交换（气/液各自动量方程）、非等温与
成核/生长动力学、液膜夹带与段塞合并破碎、管件/阀门模型。

### 4.3 Python 绑定 `src/bindings` + `python/chflow`
pybind11（FetchContent v2.13.6，`-DCHFLOW_WITH_PYTHON=ON`，缺少 Python
开发头文件时自动跳过）。模块 `chflow_core` 导出
`OneDNetwork / NetworkPipe / FluidProperties / SteadyControls /
SteadyNetworkResult / PipeFlowResult` 与 `solve_steady_network /
load_network / load_steady_controls / churchill_friction`；瞬态侧导出
`TransientCase / TransientNodeConfig / TransientPipeConfig / TransientSource /
TransientSample / TransientSeries / TransientTimeControls`、可单步推进的
`TransientNetwork(step/time/…)` 以及 `run_transient / hydrate_kinetics_step /
load_transient_case / load_transient_controls`；物性侧导出
`Component / Mixture / RootPreference / EosState / FlashControls / FlashResult /
BilinearTable / FlashCase` 与 `component / wilson_k / pure_vapor_pressure /
peng_robinson_state / tp_flash / bubble_pressure / liquid_density_table /
vapour_density_table / vapour_fraction_table / load_flash_case`。C++ 侧的
`cf::Error` 自动转成 Python `RuntimeError`。Python 包 `chflow` 只是
再导出层。求解与闪蒸调用释放 GIL。

### 4.4 CLI `cli/`（Rust）
零依赖 std-only 二进制 `chuanflow`：
`chuanflow run steady|transient|pvt <case.dict>` 通过 `$CHUANFLOW_BIN_DIR`（或 PATH）
定位 `chflow-1d-steady` / `chflow-1d-transient` / `chflow-pvt` 并转发
stdout/退出码；`chuanflow version`。
定位为自动化入口（HPC/CI），不做计算。

## 5. src/3d 三维内核

### 5.1 网格 `PolyMesh.h` / `blockMesh.cpp` / `gmshMesh.cpp` / `meshLoader.cpp`
`polyMesh` 是唯一的网格容器：面心 `faceCentres`、面法积 `faceSf`、
`owner/neighbour`（约定：`Sf` 由 owner 指向 neighbour，边界面上为外法向）、
单元体积/重心、面面积、`deltaCoeffVec = (cN−cO)/|cN−cO|²`。
容器与文件格式无关：只要给出点、面顶点环、owner/neighbour 与 patch，
几何量（含面的朝向）都由拓扑重算，因此非结构网格与结构化网格走同一条路径。
单元体积和重心用面扇形三角化的散度定理精确计算（一阶矩/体积），
对六面体、四面体、三棱柱都给出真实重心——早期按面面积加权的重心会让
斜单元的 `cO→面心` 方向偏离面法向，非正交补偿随之失真。

- `BlockMesh::fromDict`：六面体结构化网格（blockMesh 风格，
  `min/max/cells` + 6 个边界面 patch）；
- `GmshMesh::fromFile`：Gmsh **ASCII** `.msh`（2.2 与 4.1 两种写法），
  一阶四面体(4)/六面体(5)/三棱柱(6)；面由体单元连接关系重建并按顶点集
  去重（同一面被第三个单元引用 → 非流形报错），patch 名取表面 physical
  group（四边形边界面也能匹配被三角化 Split 后的两个三角形），文件不带
  表面分组时回退到包围框平面名 `xmin…zmax` 并告警，两者都不匹配则报错；
  金字塔单元与任何二阶单元直接拒绝而不是悄悄退化成角点；
- `MeshLoader::fromCase`：case.dict 选型——`block{}` 走结构化，
  `mesh { format gmsh; file "..."; patches { inlet { type inlet; } } }`
  走文件；相对路径先按 cwd、再按 case 文件所在目录解析，因此算例可以把
  网格和 `.dict` 放在一起、从任意目录运行。

### 5.2 场系统 `Fields.h` / `volFields.cpp`
- `volScalarField` / `volVectorField`：单元场 + 每 patch 一个边界条件
  （`fixedValue` / `zeroGradient` / `symmetry`…），`faceValue/linearise`
  提供边界面的 `Tf = c·Tp + k` 线性化；
- `surfaceScalarField`：面法积通量场（`phi`），`correctBoundaryFlux(U)`
  按 patch 上速度 BC 重算边界通量（壁面/对称 → 0）。

### 5.3 离散算子 `discrete/`
仿 OpenFOAM 的 fvc/fvm 双命名空间：

- `fvc::grad`：高斯（散度定理）梯度，`linear | corrected`。`linear` 对
  任何平面面多面体网格精确复现线性场（测试断言到 1e-14）；`corrected`
  再加两轮非正交补偿，正交网格上等价于 `linear`，斜网格上是一个收敛的
  近似——粗四面体网格（实测最坏非正交 54°、单元体积比 12:1）上补偿项会
  放大，需要靠网格质量控制而不是格式；
- `fvc::interpolate`：有界线性面插值（权重夹到 [0,1]）；
- `fvc::divPhi`：单元面通量求和（保持广度量纲，不除体积）；
- `fvm::div`：`upwind / linear / blend`，对 `phi` 组装隐式对流项；
- `fvm::laplacian`：调和平均面扩散系数；`fvm::laplacianFace` 直接接收
  面心 gamma（SIMPLE 的 `rAUf`），不做单元→面插值；
- `fvm::ddt`、`fvm::sp`：瞬态与隐式源。

非正交修正通过 `fvMatrix::registerCorrector` 挂到矩阵上：每次线性求解外
迭代用当前解重算修正面通量（OpenFOAM `correctedSnGrad` 的对应物）。

### 5.4 fvMatrix
单元方程容器：`aDiag/aOffDiag/rhs/内部与边界系数`，`solve()` 组装 CSR 后调
BiCGStab，支持外迭代（非正交修正次数）与欠松弛。矩阵为广度量纲
（对角 = 面系数之和，源项乘体积），这一约定见 5.5 的 rAU 讨论。

### 5.5 SIMPLE 求解器 `simpleSolver.cpp`
单相不可压定常回路（simpleFoam 结构，全部自研）：

1. 动量：`div(phi,U) − laplacian(nu,U) = −V grad(p)`，欠松弛通过
   **对旧场混合**实现：`U = Uold + α_U(U* − Uold)`（不是对矩阵零初值混合）；
2. 预测通量 `phi = interp(U) & Sf`；
3. 压力修正：`laplacian(rAUf, p') = div(phi)`，`p'` 继承 `p` 的 BC
   （fixedValue→0 Dirichlet，否则 zeroGradient）；
4. 校正：`U −= rAU grad(p')`，`p += α_P p'`，然后显式面通量校正
   `phi_f −= rAUf_f·SnGrad(p')·|Sf| + 非正交补偿`（公共件 `correctFluxes`，
   即 OpenFOAM 的 `phi -= pEqn.flux()`，见 5.6），校正后的 `phi` 在散度
   意义下按构造为零；
5. 收敛判据：预测通量散度 `contErr`；当动量欠松弛使 `contErr` 出现
   滞后地板时，用 `|p'|` 衰减作为真实的外迭代收敛信号（`pcorrTol`）。

两个刻在实现里的经验结论（回归历史，勿回退）：

- **rAU 必须与 p' 方程、通量校正、U 校正三处使用同一算子**，校正后的
  `phi` 才能严格散度自由；校正和 p' 方程用同一系数即可自洽；
- `rAU = α_U/aDiag`（不含体积）。由于 p' 方程与通量校正使用同一算子，
  `rAU` 整体缩放不破坏散度自由性，只缩放 `p'` 的幅值：使用物理响应
  `V/aDiag` 会把通量闭合所需的整段压力斜坡直接注入 `p`，SIMPLE 发散；
  Jacobi 缩放隐式地对压力做了欠松弛，是经典 SIMPLE 的稳定选择。
  代价是预测散度存在与滞后成正比的地板值，因此测试用 `pcorrTol` 判收敛。

两条 M3 期间实测出来的**已知限制**（不是 bug，后续要处理）：

- 压力场在 Jacobi 尺度上累加，因此物理量由**流向压力梯度**决定的算例，
  速度场先收敛、`p` 远慢于它；需要一致动量插值（Rhie-Chow）才能让压力与
  通量同速收敛；
- 收敛判据只看连续性，而动量平衡不在其中：流向均匀的算例（单列渠道，
  所有 `phi` 恒为 0）在第 1 次外迭代就读到 `continuity 0` 并退出，壁面剪
  切还没进入平衡。湍流渠道回归因此用 PISO 跑同一套物理——时间尺度会保
  住耦合；给 SIMPLE 加动量残差判据是后续的独立改动。

### 5.6 压力修正公共件 `pressureCorrection.cpp`
SIMPLE 与 PISO 共用四个算子，保证两种耦合闭合的是**同一个离散连续方程**：

- `updatePhiInternal`：内部面 `phi = interp(U) & Sf`（不动边界面）；
- `interpolateToFaces`：单元系数（`rAU`）按面求值，供 `fvm::laplacianFace`；
- `correctFluxes`：`phi_f −= rAUf_f·SnGrad(p)·|Sf| + 非正交补偿`，边界面用
  压力场的 `linearise`（fixedValue 给 Dirichlet，zeroGradient 不贡献）——
  即 OpenFOAM 的 `phi -= pEqn.flux()`；
- `continuityError`：`max|div(phi)|` 除以最大面面积，作为无量纲残差。

### 5.7 瞬态 PISO 求解器 `pisoSolver.cpp`
单相不可压瞬态回路（pisoFoam 结构，全部自研）：

1. 每步一次动量：`ddt(U)/dt + div(phi,U) − laplacian(nu,U) = V(f − grad pⁿ)`；
2. `HbyA = rAU·(H(U) + V grad pⁿ)/V`，即不含压力梯度的动量响应
   （OpenFOAM 的 `UEqn.H()*rAU()`），边界按 `constrainHbyA` 处理：速度
   给定的 patch 保留其面值，其余 zeroGradient；
3. 每个 corrector：`phi = interp(HbyA) & Sf` → 隐式求解**全压力**
   `laplacian(rAUf, p) = div(phi)`（以 pⁿ 为初值）→ `U = HbyA − rAU grad(p)`
   → `phi −= rAUf SnGrad(p)|Sf|`；
4. `rAU = V/aDiag`（物理动量响应）。因 `fvm::ddt` 把 `V/dt` 放在对角上，
   `rAU ≤ dt` 天然约束了每次校正，所以 PISO 不需要 SIMPLE 那套欠松弛
   （对比 5.5 的 Jacobi 缩放结论）；
5. 收敛判据：校正后的 `max|div(phi)|`，实测降到 1e-18 量级（机器精度），
   corrector 一次即退出。

两个刻在实现里的经验结论（回归历史，勿回退）：

- **必须解全压力 p，不能累加修正量 p'**：同位网格上 `interp(U)` 对
  压力模式是"盲"的（交错场的半格平均相消），累加 `p = Σp'` 会把一个
  钟形压力模态放大（实测每步 1.6~1.7 倍），最终污染动量对角导致发散。
  这里压力只通过面的紧凑 `SnGrad` 进入通量——正是压力方程反演的算子，
  该模态无法被泵起；
- **稳定性上限**：压力 Poisson 在拉长网格上各向异性，其最弱模态（流向
  压力斜坡 + 横向涡卷）只靠黏性衰减。启动渠道算例实测：扩散数
  `nu·dt/dy² ≤ 0.5` 时严格稳定且精确；接近 1 时该模态从舍入误差起振并
  每步放大，直至线性求解器报 NaN（app 会显式告警，不会静默继续）。
  瞬态算例请把 `nu·dt/min(cell edge)²` 控制在 0.5 以内。

### 5.8 应用与 IO
- `chflow-scalar-transport <case.dict>`：速度驱动的对流扩散标量求解
  （`fvm::ddt + fvm::div + fvm::laplacian` 一次求解）；
- `chflow-simple-foam <case.dict>`：不可压 SIMPLE；
- `chflow-piso-foam <case.dict>`：不可压瞬态 PISO，`piso{dt endTime
  correctors continuityTol linearTol linearMaxIter nonOrthCorrectors
  reportEvery force}`，步进循环 `dt = min(dt, endTime − t)`，逐步打印
  `t / uMax / cont / correctors` 与结束时第一条 x 截面的 `profile y u`
  供解析核对，`piso done: steps ... worstContinuity ...` 收尾；
- `VtkWriter`：ASCII `.vtu`（单元中心场），ParaView 直接可读。

### 5.9 湍流模型层 `turbulence/`
RANS 两方程模型 + 壁面函数，全部运动黏度形式（密度 1），与不可压求解器一致。

**接口方向**：`IEddyViscosity`（`include/cf/fvm/EddyViscosity.h`）声明在 FVM
层，只有一个方法 `nuEffFace()`；SIMPLE/PISO 因此以
`IEddyViscosity* turb = nullptr` 结尾，不引入 `turbulence → fvm` 的反向依赖。
`RansModel : IEddyViscosity` 持有 `U`/`phi` 的引用，每次 `correctTurbulence(dt)`
读当前场，app 只需把模型指针交给求解器。

**基类承担**（`turbulenceModel.cpp`）：k 的对流扩散衰减输运、壁面层次、
`nut`、以及面级 `nu + nut`（动量方程与两个湍流输运共用同一份）；子类只给
耗散变量方程、`nut` 公式、扩散系数与平衡对。`LaminarModel` 是"无模型"路径：
`nut` 恒为 0、`hasDissipation()` 为 false，但仍然报告网格的 `y+`。

**壁面链**（`WallFunction.h/.cpp`）：
1. 壁面积由 `polyMesh` 的 patch type == `"wall"` 给出，`wallDist` 取
   `deltaCoeffVec · n̂`（owner 中心到面的法向距离）；
2. Spalding 全范围对数率写成 `y+` 关于 `u+` 的显式式
   `y+ = u+ + e^{−κC}[e^{κu+} − 1 − κu+ − (κu+)²/2 − (κu+)³/6]`，
   κ=0.41、C=ln(E)/κ=5.5676。反演不做两次反函数，而是合并成**一个**关于
   `u+` 的单调方程 `u+·y+(u+) = u_par·y/ν` 二分求根；
3. `nut_wall = u_τ²y/u_par − ν`，于是壁面动量扩散系数
   `(ν+nut_wall)|Sf|/y = u_τ²|Sf|/u_par`——动量方程在壁上拿到的正是剪应力；
4. `k` 在壁上取 zeroGradient（它的值由剪切决定），耗散变量取面对值：
   ε = `u_τ³/(κy)`，ω = `√(ω_vis² + ω_log²)`，
   `ω_vis = 6ν/(β1y²)`、`ω_log = u_τ/(√β*κy)`——二次合成在次层与对数层
   各自给出正确极限，因此分辨次层与不分辨次层的网格都能用。
   面对值通过 `volScalarField::setFaceValue` 写入、由 `faceValue()/linearise()`
   读回。

**平衡壁面单元**（`wallEquilibrium()` + `imposeWallEquilibrium()`，M3 实测需求）：
壁面函数网格上，紧贴壁面的单元**不能**靠输运到达对数层平衡态——面的 Dirichlet
由剪切给定，而单元值只通过分子黏度向它扩散；更要害的是该单元的应变率是相对
无滑移面值算的 `U_P/(dy/2)`（`y+≈52` 时比对数层的 `u_τ/(κy)` 大数倍），于是
`P = nut|S|²` 经 `nut = C_μk²/ε` 反馈成 `k²` 增长。实测：三步之内
`k = 1.87e8 m²/s²`、`nut = 9.6e11`，随后每个线性求解器报 NaN。
标准做法是让该单元变成代数单元：输运照旧组装（它负责把剪切带向外层），但
求解后**替换**壁面单元的 `k`/耗散变量为对数层平衡对
`k = u_τ²/√C_μ`（SST 用 `√β*`）、`ε = u_τ³/(κy)`、`ω = u_τ/(√β*κy)`，
再重建 `nut`。这对平衡值自带 `nut = κu_τy`（对数层涡黏度），外层看到的就是
正确的混合。一个单元同时触到两面时由较强剪切的那面给定**成对**的值。

**壁面距离**：`cellWallDistance()` 以壁面相邻单元为源、在单元图上跑 Dijkstra，
代价为跨面的中心距；面存储为 owner<neighbour，所以两侧都要能到达（此前只走
owner→neighbour，导致上层行永远收不到距离）。这是正交网格上的精确值、非正交
网格的上界，对只需要"量级为一个单元尺寸"的 F1/F2 足够。

**常数**：k-ε `C_μ=0.09, C1=1.44, C2=1.92, σk=1, σε=1.3`；
SST `a1=0.31, β*=0.09, β1=0.075, β2=0.0828, σk1=0.85, σk2=1, σω1=0.5,
σω2=0.856, ζ*=0.3`，`γ = β/β* − ζ*/σ_ω` 保持自由流 `C_μ=0.09`，
F1/F2 为 `tanh(arg^4)` / `tanh(arg²)` 的标准混合，生产项用 Menter 限幅器
`Pk ≤ 10β*kω`；`kFloor = dissFloor = 1e-12`，湍流方程自带的线性求解设置
`1e-8 / 2000` 次。

**case 侧**：`turbulence{ model intensity lengthScale relaxation }`，缺块即
laminar；`lengthScale ≤ 0` 回落为包围盒对角线的 7%（水力直径经验值）；
`intensity < 0` 与 `relaxation ∉ (0,1]` 直接报错。`initial{ U (Ux Uy Uz) }`
给一个起步速度场（`caseInit.h`），壁面函数需要一个已有速度的边界层可谈。
vtu 额外写 `nut`/`k`/耗散变量，收尾打印 `u_tau / y+ mean / y+ max`。

## 6. 概念映射速查（OpenFOAM → ChuanFlow）

| OpenFOAM-dev | ChuanFlow |
|---|---|
| `polyMesh` / `primitiveMesh` | `polyMesh`（`src/3d/include/cf/fvm/PolyMesh.h`） |
| `blockMeshDict` | `BlockMesh::fromDict`（case.dict 的 `block{}`） |
| Gmsh `.msh` 导入 | `GmshMesh::fromFile` + `MeshLoader`（case.dict 的 `mesh{}`，ASCII 2.2/4.1，一阶 tet/hex/prism） |
| `volScalarField` / `surfaceScalarField` | 同名类（`Fields.h`） |
| `fvc::` / `fvm::` | 同名命名空间（`FvMatrix.h`） |
| `fvVectorMatrix` / `solve(UEqn == −grad(p))` | `fvMatrix`（多分量）+ apps 手写回路 |
| `lduMatrix` + PBiCGStab | `SparseMatrix` + 右预处理 BiCGStab |
| `simpleFoam` | `solveIncompressibleSimple` + `chflow-simple-foam` |
| `pisoFoam` | `solveTransientIncompressiblePiso` + `chflow-piso-foam` |
| `system/fvSchemes` 字典 | case.dict 的 `schemes{}`（`loadSchemes`） |
| `correctedSnGrad` 非正交修正 | `fvMatrix::registerCorrector` |
| `pEqn.flux()` | `correctFluxes`（`pressureCorrection.cpp`，SIMPLE/PISO 共用） |
| `constrainHbyA` / `UEqn.H()*rAU()` | pisoSolver 的 `HbyA` 组装与 patch 约束 |
| `basicThermo` / `pengRobinson`（reacting 系列的热力学底座） | `cf::pvt::Mixture` + `pengRobinsonState`（`src/pvt`，等温闪蒸由自研 `tpFlash` 提供） |
| `thermophysicalDict` / `species` | case.dict 的 `composition { name fraction; }`（`loadComposition`） |
| `incompressible::turbulenceModel` / `eddyViscosity` | `IEddyViscosity`（FVM 层）+ `RansModel`（`src/3d/src/turbulence`） |
| `kEpsilon` / `kOmegaSST` | `KEpsilon` / `KOmegaSST`（运动黏度形式，自研离散与线性化） |
| `nutWallFunction` / `epsilonWallFunction` / `omegaWallFunction` | `wallState()` + `wallDissipation()` + `diffusionFace()` 的壁面分支 |
| `kqRWallFunction`（壁面单元代数化） | `wallEquilibrium()` + `imposeWallEquilibrium()`（平衡壁面单元） |
| `wallDist` / `nearWallDist` | `RansModel::cellWallDistance()`（单元图 Dijkstra） |
| `fvOptions` / `turbulenceProperties` | case.dict 的 `turbulence{}`（`TurbulenceConfig::fromCase`） |

## 7. 测试体系
- `tests/unit/kernelTests.cpp`（CTest `kernel.unit`）：字典解析、网格几何
  恒等式（散度定理 Σ Sf = 0）、场 BC 线性化、矩阵装配与 BiCGStab 收敛；
- `tests/unit/networkTests.cpp`（CTest `network.unit`）：Churchill 层流极限
  64/Re、单层流管 Hagen-Poiseuille 正/反流、静水平衡（dp=ρg·dz→Q=0）、
  单井 Newton 闭式解、对称 Y 分流、case.dict 解析与求解；
- `tests/unit/transientTests.cpp`（CTest `transient.unit`）：水合物动力学
  精确指数（1e-12）、单步 2×2 隐式 Euler 闭式对照（1e-10）、静止平衡、
  无速率时相不变、带速率时液量指数衰减、slug case.dict 往返与级联递推；
- `tests/unit/pvtTests.cpp`（CTest `pvt.unit`）：全部是 EOS 内部不变量而非
  文献数值——临界点三重根 `Z = 0.30740`、`v_c/b = Z_c/B_c`、低压理想气体
  极限、纯组分饱和压处逸度相等（1e-9）与随 T 单调、纯组分泡点 ≡ 饱和压
  （1e-5）、超临界饱和压为 0、多组分闪蒸的组分物料衡算（1e-10）/相分和/
  组分逸度相等（1e-8）/`K_C1 > K_C5`/轻端在气相富集/升压降 β/高压塌成单相
  液/真空塌成单相气、混合物泡点被压力扫描包络并随 T 上升、双线性表节点
  与直接求值一致（1e-14）、气相密度随压升、两相区液相密度沿包络下降、
  单相区密度对 p 单调增对 T 单调减、case.dict 往返与未知组分名报错；
- `tests/unit/pisoTests.cpp`（CTest `piso.unit`）：非正时间步直接报错、
  无黏第一步每个单元 `u = f·dt`（1e-9）且压力停在选择基准、1500 步瞬态
  推进后剖面逐点等于**离散**稳态抛物线（含半格壁面偏移，1e-5 相对峰速）、
  三阶差分逐点满足 `nu∇²u = −a`（1e-9）、横向速度保持 0、1500 步内校正后
  `max|div(phi)| < 1e-12`；
- `tests/unit/gmshMeshTests.cpp`（CTest `gmshMesh.unit`）：单六面体（无
  physical group → 包围框平面命名 + 类型覆盖）、2×2×2 六面体与
  `BlockMesh` 几何指纹逐元素相等、同一四面体网格的 2.2/4.1 两种写法读出
  同一网格、三棱柱渠道 20 单元/28 内部面/6 patch 与真实重心
  `(2/3, 1/30, 0.05)`、逐单元 `Σ Sf = 0`（闭合）、`fvc::grad(z)` 精度、
  `MeshLoader` 的 `block{}`/`mesh{}` 选型与两个错误分支；
- `tests/unit/turbulenceTests.cpp`（CTest `turbulence.unit`）：Spalding 对数
  率的两端极限（次层 `u+=y+`、对数层趋近 `ln(Ey+)/κ`）、反演单调、F1/F2 混合
  梯级（含"外层 F1 随 y 递减"与"在 0.15 处正在切换"）、`nut` 定义式、
  壁面链恒等式 `u_par(ν+nut_wall)/y = u_τ²` 与壁面耗散 Dirichlet、
  **平衡壁面单元**（`k = u_τ²/√C_μ` 与配对的 ε/ω 逐面成立、`nut = κu_τy`
  对 k-ε 精确而对 SST 只允许被剪切分支压低、闭合只作用于触壁单元）、
  均匀衰减对闭式解（k-ε `k = k0[1+(C2−1)(ε0/k0)t]^{1/(1−C2)}` 与
  `ε = ε0(k/k0)^{C2}`、SST 在近壁分支 `ω = ω0/(1+β1ω0t)`、
  `k = k0(1+β1ω0t)^{−β*/β1}`）、壁面距离、laminar 路径、配置校验；
- `tests/cases/*`：十一个物理回归算例 + `tests/check_*.py` 断言：
  - `steadyNetwork`（两井集输层流管网：节点压力/流量/Re/压降对比
    Hagen-Poiseuille 串并联闭式解，含重位水头，1 次 Newton 收敛）
  - `transientBlowdown` / `transientRlc` / `transientSlug`（见 4.2）
  - `pvtFlash`（四组分天然气 320 K 压力扫描：闪蒸恒等式——相分和、
    `β·y + (1−β)·x = z`、ρ_L > ρ_V、泡点两侧相态塌陷、表节点 = 直接求值）
  - `advectionDiffusion`（Pe=10 一维对流扩散，解析解对比）
  - `pureDiffusion`（Laplace/Poisson，线性场精确再现）
  - `gmshDiffusion`（非结构网格上的同一道题：三棱柱渠道 `tests/meshes/
    prismChannel_22.msh`，20 单元、实测最坏非正交 78°，两端 Dirichlet
    0/1；`T = z` 仍是离散解的精确解（斜面的补偿项恒为 0），断言误差
    1e-10——网格来自 Gmsh 文件而不是 `block{}`）
  - `poiseuilleChannel`（平面渠道 Poiseuille：抛物线剖面、通量守恒、
    `v≈0`、中心线与解析压降对比；SIMPLE 收敛约 1500 外迭代 / 37 s）
  - `startedChannel`（渠道瞬时启动：10×20×1 网格、体积力驱动、
    `nu dt/dy² = 0.2`、6000 步 / 12 s 物理时间；速度历史对比解析级数
    （早期相对误差 8.7e-3，为 Euler 一阶时间离散；后半程 2.8e-7）、
    终态剖面逐点等于离散抛物线（1.4e-8）、壁面偏移实测等于
    `(a/2nu)(dy/2)²`、校正后散度 1.4e-18）
  - `turbulentChannel`（充分发展湍流平面渠道，1×16×1 单列、体积力驱动、
    `y+ ≈ 52`、750 步 / 0.5 s，两个模型各跑一遍）：壁面剪应力对积分平衡
    `u_τ² = f·δ`（k-ε 2.2e-5、SST 3.7e-4 相对误差）、首格 `u+` 对 Spalding
    对数率（2.7e-7 与 6.2e-7，说明壁面链在整条求解链上是自洽的）、剖面对称
    单调、`u_tau/y+` 落在壁面函数区间、终态漂移 <1%、壁面单元的平衡对与
    `nut = κu_τy` 从 vtu 复核；摩阻对 Pope 光滑平板律
    `1/√f = 2log10(Re_b√f) − 0.6` 的偏差为 k-ε −14.5%、SST +5.4%
    （`Ub/u_τ` 17.8 与 21.9，期望 20.8）——8 个单元跨过 `y+ = 52` 到中心线的
    网格只能把这个数核对到网格自己的精度，因此带宽取 20%，而近壁闭合严格核
    查。checker 把 case 复制到临时目录再改 `model`，产物不落进仓库。
- `tests/python/test_bindings.py`（CTest `python.bindings`）：Python 绑定
  构建/求解/读结果与闭式解对比、异常翻译、case.dict 往返；瞬态侧校验
  罐排气 RC 包络、逐步 `TransientNetwork.step` 与批量 `run_transient`
  结果一致、水合物动力学指数；pvt 侧校验临界点 Z、纯组分泡点 ≡ 饱和压、
  `tp_flash` 的物料衡算与逐组分逸度相等（由 `ln_phi` 复算）、case.dict
  加载 + 泡点两侧相态、`liquid_density_table` 与直接求值一致；
- `cli`（CTest `cli.unit` + `cli.runSteady` + `cli.runTransient` +
  `cli.runPvt`）：cargo 参数解析单测 + 通过 CLI 端到端跑稳态/瞬态管网与
  pvt 闪蒸冒烟测试（存在 cargo 时才注册）；
- 运行：`cmake -S . -B build && cmake --build build && ctest --test-dir build`
