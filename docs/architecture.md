# ChuanFlow 架构说明（当前实现）

> 本文档描述仓库中已经落地的部分：`src/common` 基础层、`src/1d` 一维稳态管网
> 内核（M1）与 `src/3d` 三维 FVM 内核。
> README 中的 1D 瞬态 / PVT / 耦合 / GUI 等模块属于路线图（M1 后续与 M2–M6），尚未实现。
>
> ChuanFlow 不是 OpenFOAM 的 Fork：以下设计参考了 OpenFOAM-dev 的工程思想
> （网格拓扑、场分层、fvc/fvm 分离、SIMPLE 算法回路），但所有数据结构、离散
> 格式与求解器均为独立自研实现。

## 1. 总体分层

```
apps (chflow-1d-steady / chflow-scalar-transport / chflow-simple-foam)
  ├── src/1d          一维稳态管网内核（节点法水力计算，对标 Pipesim 的基础能力）
  │      └── src/common  基础层（字典解析、日志/异常）
  └── src/3d          三维 FVM 内核（mesh / fields / discrete / solvers / io）
        └── src/common  基础层（向量、Field 容器、稀疏矩阵+Krylov求解、字典解析、日志/异常）
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

## 3. src/1d 一维稳态管网内核（M1）

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

M1 剩余的 Python 绑定与 CLI 封装、以及井筒/热力耦合与环路迭代（对标完整
Pipesim）尚未实现。

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
- `tests/cases/*`：四个物理回归算例 + `tests/check_*.py` 断言：
  - `steadyNetwork`（两井集输层流管网：节点压力/流量/Re/压降对比
    Hagen-Poiseuille 串并联闭式解，含重位水头，1 次 Newton 收敛）
  - `advectionDiffusion`（Pe=10 一维对流扩散，解析解对比）
  - `pureDiffusion`（Laplace/Poisson，线性场精确再现）
  - `poiseuilleChannel`（平面渠道 Poiseuille：抛物线剖面、通量守恒、
    `v≈0`、中心线与解析压降对比；SIMPLE 收敛约 1500 外迭代 / 37 s）
- 运行：`cmake -S . -B build && cmake --build build && ctest --test-dir build`
