# 不可压 Navier–Stokes Python 绑定设计说明（格 1i）

本说明先于实现落盘，固定 `chflow_core` 不可压 NS 绑定的暴露面与四件核心
契约（面通量 `phi`、边界条件、场读写、时间推进）。实现只驱动既有内核求解器
（`solveTransientIncompressiblePiso` / `solveIncompressibleSimple`），不改内核
数值。

## 1. 范围与排除

- **本格覆盖**：层流（无湍流模型、无能量方程）的瞬态 PISO 与稳态 SIMPLE
  两条不可压压力–速度耦合链路；`turb` 与 `energy` 钩子均传 `nullptr`。
- **本格排除**（后续增量再接，不在本格占位）：
  - RANS/LES 涡粘模型作为 `IEddyViscosity` 钩子（已由 LES/RSM 各自绑定覆盖，
    不走本格的通用 NS 入口）；
  - 温度场 / 浮力 Boussinesq 耦合（`IEnergyCoupling`）；
  - `p_rgh` 静水压力吸收模式（app 侧编排，本格不还原 `pTrue` 重构）；
  - 多相、可压、变密度与非牛顿流体。

## 2. 四件契约

### 2.1 面通量 `phi`（surfaceScalarField）

- `phi` 为面心标量通量，内部面由速度场重建，边界面由 `correctBoundaryFlux(U)`
  按边界条件闭合；绑定内部持有，不向 Python 暴露 NumPy 视图（面通量的
  owner/neighbour 寻址与场字段不同，留待后续面场增量）。
- 瞬态 PISO：每步 `solveTransientIncompressiblePiso` 闭合的通量即下一步动量
  方程的对流通量；稳态 SIMPLE：`solveIncompressibleSimple` 返回的通量是
  Rhie–Chow 一致插值后的闭合通量。
- 连续性判据沿用内核结果：PISO 步 `continuityError`、SIMPLE
  `finalContinuityError`，绑定只读数，不自造容差。

### 2.2 边界条件（`patches{}`）

复刻 `pisoFoamApp` / `simpleFoamApp` 的补丁读法，逐补丁：

- 速度 `U`：补丁声明 `U (ux uy uz)` 时按 `UType`（缺省 `fixedValue`）落固定值；
  未声明 `U` 时，`wall` 类型补丁按**无滑移** `fixedValue 0`，其余补丁按
  `zeroGradient`（避免旧默认把壁面静默变成自由滑移）。
- 压力 `p`：补丁声明 `p` 时按 `pType`（缺省 `fixedValue`）；未声明按
  `zeroGradient 0`。
- 补丁名经 `patchIdByName` 解析，解析失败抛 `RuntimeError`。

### 2.3 场读写

- 内核场：`volVectorField U`、`volScalarField p`、`surfaceScalarField phi`。
- 读出：经 `field_numpy.h` 转 NumPy 拷贝，float64、C-order；`U` 形状 `(n,3)`，
  `p` 形状 `(n,)`；Python 持有拷贝，不持内核场引用。
- 写入：本格不开放 `from_numpy` 回写初场；初场只接受 case dict 的
  `initial { U (…); }` 播种（缺省静止启动）。场回写入口沿用 `ScalarField3d` /
  `VectorField3d`（fields3d 绑定），不在本格重复。

### 2.4 时间推进

- **瞬态 PISO**：`piso{ dt; endTime; correctors; continuityTol; linearTol;
  linearMaxIter; nonOrthCorrectors; reportEvery; force; }`；`schemes.ddt` 固定
  `Euler`；步长 `dt = min(ctl.dt, endTime - t)`；每步记录 `(t, uMax,
  continuityError, correctorsUsed)`；末态汇总 `worstContinuity`、步数。
- **稳态 SIMPLE**：`simple{ maxIterations; continuityTol; pcorrTol; momentumTol;
  relaxU; relaxP; pRelaxStart; nonOrthCorrectors; force; }`；一次外迭代闭合，
  返回 `iterations / finalContinuityError / finalMomentumResidual / converged`。
- 体作用力 `force (fx fy fz)` 进 `bodyForce`；本格不做 `p_rgh` 吸收。

## 3. 暴露面（拟定）

- `solve_piso_ns_case(case_dict_path) -> NsPisoCase`
  - `U (n,3)`、`p (n,)`、`history`（`t / u_max / continuity / correctors`）、
    `steps`、`worst_continuity`、`nu`、`end_time`、`dt`、`n_cells`。
- `solve_simple_ns_case(case_dict_path) -> NsSimpleCase`
  - `U (n,3)`、`p (n,)`、`iterations`、`continuity`、`momentum_residual`、
    `converged`、`nu`、`n_cells`。

## 4. 对账锚点（复用既有校验器，不自创阈值）

- 瞬态：`tests/check_piso.py` 的启动通道级数解（半离散步长壁面偏移显式建模），
  晚段历史相对误差、终态剖面、连续性 `1e-12` 量级、单调启动无超调。
- 稳态：`tests/check_poiseuille.py` 的充分发展抛物线剖面、横向速度近零、
  截面通量守恒；容差沿用校验器（剖面 `0.012`、中心线 `8%`、通量 `2%`）。
- 拒绝路径：缺 `piso{}`/`simple{}` 控制块、非正 `nu`、未知补丁等由内核校验
  抛出，绑定透传原信息。

## 5. 边界

只改 `src/bindings/**`、`tests/python/**`、`tests/CMakeLists.txt`；不触碰
`gui/**`、`cli/**`、`docs/**`、`src/coupling/apps/**`、`src/3d/apps/**`、
`src/1d/apps/**` 与内核数值。
