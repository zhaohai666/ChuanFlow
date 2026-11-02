# ChuanFlow 架构说明（当前实现）

> 本文档描述仓库中已经落地的部分：`src/common` 基础层、`src/1d` 一维管网
> 内核（M1 稳态 + M2 首增量：瞬态水力 / 持液率输运 / 水合物动力学）与
> `src/3d` 三维 FVM 内核，以及 `src/bindings` Python 绑定与 `cli/` Rust 入口。
> README 中的 PVT / 1D↔3D 耦合 / GUI 等模块属于路线图（M3–M6），尚未实现。
>
> ChuanFlow 不是 OpenFOAM 的 Fork：以下设计参考了 OpenFOAM-dev 的工程思想
> （网格拓扑、场分层、fvc/fvm 分离、SIMPLE 算法回路），但所有数据结构、离散
> 格式与求解器均为独立自研实现。

## 1. 总体分层

```
apps (chflow-1d-steady / chflow-1d-transient / chflow-scalar-transport / chflow-simple-foam)
  ├── src/1d          一维管网内核（M1 稳态节点法 + M2 瞬态/持液率/水合物）
  │      ├── src/bindings  pybind11 模块 chflow_core（Python 包 python/chflow）
  │      └── src/common  基础层（字典解析、日志/异常）
  └── src/3d          三维 FVM 内核（mesh / fields / discrete / solvers / io）
        └── src/common  基础层（向量、Field 容器、稀疏矩阵+Krylov求解、字典解析、日志/异常）
cli/  Rust 自动化入口 chuanflow（run steady|transient <case.dict>，定位并转发引擎二进制）
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

## 3. src/1d 一维管网内核

### 3.1 稳态节点法（M1）

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

M1 剩余的井筒/热力耦合与环路迭代（对标完整 Pipesim）尚未实现。

### 3.2 瞬态水力 + 持液率输运 + 水合物动力学（M2 首增量）

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

### 3.3 Python 绑定 `src/bindings` + `python/chflow`
pybind11（FetchContent v2.13.6，`-DCHFLOW_WITH_PYTHON=ON`，缺少 Python
开发头文件时自动跳过）。模块 `chflow_core` 导出
`OneDNetwork / NetworkPipe / FluidProperties / SteadyControls /
SteadyNetworkResult / PipeFlowResult` 与 `solve_steady_network /
load_network / load_steady_controls / churchill_friction`；瞬态侧导出
`TransientCase / TransientNodeConfig / TransientPipeConfig / TransientSource /
TransientSample / TransientSeries / TransientTimeControls`、可单步推进的
`TransientNetwork(step/time/…)` 以及 `run_transient / hydrate_kinetics_step /
load_transient_case / load_transient_controls`。C++ 侧的
`cf::Error` 自动转成 Python `RuntimeError`。Python 包 `chflow` 只是
再导出层。求解调用释放 GIL。

### 3.4 CLI `cli/`（Rust）
零依赖 std-only 二进制 `chuanflow`：
`chuanflow run steady|transient <case.dict>` 通过 `$CHUANFLOW_BIN_DIR`（或 PATH）
定位 `chflow-1d-steady` / `chflow-1d-transient` 并转发 stdout/退出码；
`chuanflow version`。
定位为自动化入口（HPC/CI），不做计算。

## 4. src/3d 三维内核

### 4.1 网格 `PolyMesh.h` / `blockMesh.cpp`
`polyMesh` 是唯一的网格容器：面心 `faceCentres`、面法积 `faceSf`、
`owner/neighbour`（约定：`Sf` 由 owner 指向 neighbour，边界面上为外法向）、
单元体积/重心、面面积、`deltaCoeffVec = (cN−cO)/|cN−cO|²`。
`BlockMesh::fromDict` 生成六面体结构化网格（blockMesh 风格，
`min/max/cells` + 6 个边界面 patch）。

### 4.2 场系统 `Fields.h` / `volFields.cpp`
- `volScalarField` / `volVectorField`：单元场 + 每 patch 一个边界条件
  （`fixedValue` / `zeroGradient` / `symmetry`…），`faceValue/linearise`
  提供边界面的 `Tf = c·Tp + k` 线性化；
- `surfaceScalarField`：面法积通量场（`phi`），`correctBoundaryFlux(U)`
  按 patch 上速度 BC 重算边界通量（壁面/对称 → 0）。

### 4.3 离散算子 `discrete/`
仿 OpenFOAM 的 fvc/fvm 双命名空间：

- `fvc::grad`：最小二乘梯度，`linear | corrected`（corrected = 精确对偶
  + 非正交补偿）；
- `fvc::interpolate`：有界线性面插值（权重夹到 [0,1]）；
- `fvc::divPhi`：单元面通量求和（保持广度量纲，不除体积）；
- `fvm::div`：`upwind / linear / blend`，对 `phi` 组装隐式对流项；
- `fvm::laplacian`：调和平均面扩散系数；`fvm::laplacianFace` 直接接收
  面心 gamma（SIMPLE 的 `rAUf`），不做单元→面插值；
- `fvm::ddt`、`fvm::sp`：瞬态与隐式源。

非正交修正通过 `fvMatrix::registerCorrector` 挂到矩阵上：每次线性求解外
迭代用当前解重算修正面通量（OpenFOAM `correctedSnGrad` 的对应物）。

### 4.4 fvMatrix
单元方程容器：`aDiag/aOffDiag/rhs/内部与边界系数`，`solve()` 组装 CSR 后调
BiCGStab，支持外迭代（非正交修正次数）与欠松弛。矩阵为广度量纲
（对角 = 面系数之和，源项乘体积），这一约定见 4.5 的 rAU 讨论。

### 4.5 SIMPLE 求解器 `simpleSolver.cpp`
单相不可压定常回路（simpleFoam 结构，全部自研）：

1. 动量：`div(phi,U) − laplacian(nu,U) = −V grad(p)`，欠松弛通过
   **对旧场混合**实现：`U = Uold + α_U(U* − Uold)`（不是对矩阵零初值混合）；
2. 预测通量 `phi = interp(U) & Sf`；
3. 压力修正：`laplacian(rAUf, p') = div(phi)`，`p'` 继承 `p` 的 BC
   （fixedValue→0 Dirichlet，否则 zeroGradient）；
4. 校正：`U −= rAU grad(p')`，`p += α_P p'`，然后显式面通量校正
   `phi_f −= rAUf_f·SnGrad(p')·|Sf| + 非正交补偿`（OpenFOAM 的
   `phi -= pEqn.flux()`），校正后的 `phi` 在散度意义下按构造为零；
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

### 4.6 应用与 IO
- `chflow-scalar-transport <case.dict>`：速度驱动的对流扩散标量求解
  （`fvm::ddt + fvm::div + fvm::laplacian` 一次求解）；
- `chflow-simple-foam <case.dict>`：不可压 SIMPLE；
- `VtkWriter`：ASCII `.vtu`（单元中心场），ParaView 直接可读。

## 5. 概念映射速查（OpenFOAM → ChuanFlow）

| OpenFOAM-dev | ChuanFlow |
|---|---|
| `polyMesh` / `primitiveMesh` | `polyMesh`（`src/3d/include/cf/fvm/PolyMesh.h`） |
| `blockMeshDict` | `BlockMesh::fromDict`（case.dict 的 `block{}`） |
| `volScalarField` / `surfaceScalarField` | 同名类（`Fields.h`） |
| `fvc::` / `fvm::` | 同名命名空间（`FvMatrix.h`） |
| `fvVectorMatrix` / `solve(UEqn == −grad(p))` | `fvMatrix`（多分量）+ apps 手写回路 |
| `lduMatrix` + PBiCGStab | `SparseMatrix` + 右预处理 BiCGStab |
| `simpleFoam` | `solveIncompressibleSimple` + `chflow-simple-foam` |
| `system/fvSchemes` 字典 | case.dict 的 `schemes{}`（`loadSchemes`） |
| `correctedSnGrad` 非正交修正 | `fvMatrix::registerCorrector` |
| `pEqn.flux()` | simpleSolver 显式面通量校正块 |

## 6. 测试体系
- `tests/unit/kernelTests.cpp`（CTest `kernel.unit`）：字典解析、网格几何
  恒等式（散度定理 Σ Sf = 0）、场 BC 线性化、矩阵装配与 BiCGStab 收敛；
- `tests/unit/networkTests.cpp`（CTest `network.unit`）：Churchill 层流极限
  64/Re、单层流管 Hagen-Poiseuille 正/反流、静水平衡（dp=ρg·dz→Q=0）、
  单井 Newton 闭式解、对称 Y 分流、case.dict 解析与求解；
- `tests/unit/transientTests.cpp`（CTest `transient.unit`）：水合物动力学
  精确指数（1e-12）、单步 2×2 隐式 Euler 闭式对照（1e-10）、静止平衡、
  无速率时相不变、带速率时液量指数衰减、slug case.dict 往返与级联递推；
- `tests/cases/*`：七个物理回归算例 + `tests/check_*.py` 断言：
  - `steadyNetwork`（两井集输层流管网：节点压力/流量/Re/压降对比
    Hagen-Poiseuille 串并联闭式解，含重位水头，1 次 Newton 收敛）
  - `transientBlowdown` / `transientRlc` / `transientSlug`（见 3.2）
  - `advectionDiffusion`（Pe=10 一维对流扩散，解析解对比）
  - `pureDiffusion`（Laplace/Poisson，线性场精确再现）
  - `poiseuilleChannel`（平面渠道 Poiseuille：抛物线剖面、通量守恒、
    `v≈0`、中心线与解析压降对比；SIMPLE 收敛约 1500 外迭代 / 37 s）
- `tests/python/test_bindings.py`（CTest `python.bindings`）：Python 绑定
  构建/求解/读结果与闭式解对比、异常翻译、case.dict 往返；瞬态侧校验
  罐排气 RC 包络、逐步 `TransientNetwork.step` 与批量 `run_transient`
  结果一致、水合物动力学指数；
- `cli`（CTest `cli.unit` + `cli.runSteady` + `cli.runTransient`）：cargo
  参数解析单测 + 通过 CLI 端到端跑稳态/瞬态管网冒烟测试（存在 cargo 时才注册）；
- 运行：`cmake -S . -B build && cmake --build build && ctest --test-dir build`
