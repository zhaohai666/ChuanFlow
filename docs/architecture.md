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
apps (chflow-1d-steady / chflow-1d-transient / chflow-pvt / chflow-scalar-transport / chflow-simple-foam / chflow-piso-foam / chflow-vof / chflow-3d-compressible / chflow-3d-parallel / chflow-coupled-1d3d / chflow-coupled-1d3d-transient / chflow-coupled-1d3d-multi)
  ├── src/1d          一维管网内核（M1 稳态节点法 + M2 首增量：瞬态/持液率/水合物）
  │      └── src/common  基础层（字典解析、日志/异常）
  ├── src/pvt         统一物性层（PR 状态方程 / 多组分 TP 闪蒸 / PVT 表插值）
  │      └── src/common  基础层（字典解析、日志/异常）
  ├── src/bindings    pybind11 模块 chflow_core（Python 包 python/chflow）
  │      └── src/1d + src/pvt
  ├── src/3d          三维 FVM 内核（mesh / fields / discrete / solvers / turbulence / multiphase / parallel / io）
  │     └── src/common  基础层（向量、Field 容器、稀疏矩阵+Krylov求解、字典解析、日志/异常）
  └── src/coupling    分区耦合层（稳态管网的一条边换成三维详细块；M4 首增量）
        ├── src/1d + src/3d（各自独立，互不知道对方存在）
        └── src/common  基础层
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

尚未实现（pvt 后续增量）：水合物相平衡/生成条件（HH-2、kJ-i 类）、挥发油/
黑油分离器闪蒸（相关式形式的黑油物性已落地，见 §3.6）、EOS 黏度与界面
张力模型、非等温（PH/PT 焓耦合）闪蒸，以及把
pvt 密度直供 1D/3D 求解器（当前一维的 ρ/μ 仍由 case.dict 给定常数）。

### 3.6 黑油物性相关式 `BlackOil.h` / `blackOil.cpp`
油藏/集输工程标准的黑油（black-oil）物性栈，全部 SI 接口、内部用现场单位
计算并锚定文献值到 12 位：

- **Standing 泡点压力**（whitson wiki bopvt / Ahmed eq. 2-76；T 用 °F 直进，
  不做 460/459.67 的舍入差）；
- **Vazquez-Beggs** `Rs`、`Bo` 与欠饱和压缩系数 `co`（`co` 分母对退化泡点
  ——死油 Rs=0 给出 pb=0——下限约束到关联式拟合范围，防 bo 指数坍缩）；
- **饱和油密度**由单 STB 质量守恒闭式（Ahmed eq. 2-63，350γo + 0.0764γgRs
  除以 5.6146·Bo，ft³/bbl 换算）；
- **Beggs-Robinson** 死油/活油黏度、**Lee-Gonzalez-Eakin** 气体黏度；
- **Dranchuk-Abou-Kassem Z**（11 系数）：残差在 Z∈[0.05,3] 连续变号，用
  区间二分求解（阻尼 Newton 在 ρr≈1 附近发散）；
- `evaluateBlackOil(fluid, p)` 一次给出 rs/bo/freeGas/bg/z/ρg/μo/μg/μ 混合，
  `blackOilMassPerOilSm3` 给出地面油单位的管线质量流量换算；
  单元测试（`tests/unit/pvtTests.cpp`）对上述每条相关式锚定独立 Python
  转写值，并用 Ahmed 恒等式做交叉验证。

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

**两相扩展：黑油 PVT + Beggs-Brill（对标 Pipesim 两相稳态）**。case.dict 的
`fluid{}` 里写 `blackoil{gasGravity; oilDensity; gor; waterCut; temperature;
surfaceTension;}` 即切换两相模式，`sources` 变为地面油产量（sm³/s）：

- **管道本构**：Beggs-Brill (1973) 流型判别（segregated/transition/
  intermittent/distributed，过渡区按 A 权重混合）、倾斜修正 ψ 与两相摩阻
  倍数 S；沿程 24 段 march 从上游压力积分 BB 梯度（每段用局部 p 的黑油
  物性就地算 λ、HL、ρs），`marchBeggsBrillPipe`（`BeggsBrill.h`）对
  fluids 库 doctest 黄金锚点（686.9724506803469 Pa/m）吻合到 1e-6；
- **流量—压降反解**：BB 梯度对流量是 U 形的（持液率释压先降、摩阻后升），
  `g(q)=p_out(q)−p_low` 有多个根；根查找从高流量端向下括住**最大根**（稳定
  支），避免从 q=0 括到静液柱附近的近零根导致管 law 在分支间跳变——
  这是两相管网 Newton 收敛的关键；
- **管网方程组**：同一套节点法 Newton；雅可比电导取单向差分（探针伸向
  活流分支深处，不跨越静液柱拐点），死支（q=0）用拐点弦电导并回退
  1e-6 保非奇异；全局化线搜索接受死区平台上的等残差步（方向一致地推向
  拐点），活流分支要求严格下降，出物理域（压力非正/非有限）的试步跳过；
- **回归**：`tests/cases/steadyNetworkMultiphase/case.dict` 两井集输 + 
  `tests/check_steady_network_multiphase.py`（流量/质量平衡到 1e-9、HL∈(0,1]、
  流型合法、与 gor=0 干油对照：气体抬升水平管压降）；单元测试另含纯液/
  纯气退化（精确 Darcy-Weisbach + 重位）、`solvePipeFlowTwoPhase` 闭合差
  与两端管网一致（<5 Pa）。

- **节点分析：IPR–VLP 交汇（对标 Pipesim nodal analysis）**。`src/1d/
  Nodal.h` 的 `solveNodalAnalysis` 把流入动态曲线与井筒垂直管流曲线交叉，
  工作点是 `f(q) = p_vlp(q) − p_ipr(q)` 的最大根；`chflow-1d-nodal` 从
  case.dict 读 `ipr{form; reservoirPressure; …}` + `well{length; diameter;
  roughness; elevation;}` + `wellheadPressure`：

  - **IPR 三种闭式**（SI，Pa、sm³/s 地面油）：linear 达西生产指数
    `q = J(p_r − p_wf)`、Vogel 溶解气驱 `q/qmax = 1 − 0.2x − 0.8x²`（x =
    p_wf/p_r）、Jones 产能方程 `q = C(p_r² − p_wf²)^n`；各带自己的闭式
    反函数 `pWf(q)`/`rate(p_wf)` 并按定义域 [0, rateMax] 钳位成全函数；
  - **VLP 内层反解**：producer 流向上行，march 必须自底向上——
    `wellboreInletForRate` 括号扩张 + Illinois secant 反解"哪个井底流压
    march 后正好落在井口回压"（pOutlet 随 pIn 次线性增长——压力越高
    自由气体积越小、摩阻越省，残差单调可括根）；内层闸门 1e-12·p_wh
    远紧于外层 1e-7·p_r，嵌套不漏水；下坡增压剖面直接报错；
  - **最大根原则**：VLP 对 q 是 U 形（持液率释压先降、摩阻后升）而 IPR
    单调下降，f 可多次过零——仍从高流量端（IPR 容量上界）向下括住最大
    根（稳定支），与管网管 law 同一原则；油藏顶不起静液柱时如实报
    flows=false、q=0，井底取静液柱本身；
  - **回归**：三算例（`nodalVogel` 活井 gor=90：q=0.00290785579167、
    pwf=15120468.4282、静液柱 18642604.4027、HL=0.40 distributed；
    `nodalDeadOil` 孪生 q 降约 5 倍、pwf=20796742.5522——气体举升；
    `nodalDeadWell` p_r=5 MPa 顶不起静液柱）+ `tests/check_nodal.py`
    23 条断言（Vogel 恒等式在 Python 端由 pwf 重算、流动井底低于静液柱、
    golden 锚点 1e-9）+ 单元测试（闭式曲线往返 1e-12、无摩阻极限
    q* = J(p_r − p_static) 到 1e-6、过号稳定性 f(0.9q)<0<f(1.1q)、静液柱
    march 回落井口 <1e-4 Pa、死井井底=静液柱）。

- **设备单元：阀 + 泵（对标 Pipesim 稳态设备）**。pipe 字典可带
  `valve <K>`（kq² 损失系数，`dp = K·ρ|v|v/2`，恒与流向相反）与
  `pumpHead <h1> / pumpRate <q1> / pumpShutoffHead <h0>`（离心泵扬程曲线
  `H(q) = h0 − (h0−h1)(q/q1)²` 钳位 ≥ 0、对 q 偶），管 law 变为
  `dpTotal = ρg·dz + (fL/D + K)ρ|v|v/2 − ρg·H(q)`。流量—压降反解用**对
  流量括根的 Illinois secant**（残差 `G(q) = loss(q) + ρg·dz − ρg·H(q) −
  dpTotal`）而不是泵扬程定点迭代——后者在扬程曲线下降快于系统曲线上升
  （`ρg|H'| > d(loss)/dq`，本算例比值 ≈ 1.32）时不收缩，振荡后锁死在
  H=0 钳位支，流量与管网 Newton 全错；按 `G(0) = ρg·dz − ρg·h0 −
  dpTotal` 的符号选支后残差严格单调（loss 无界增长、−ρg·H 有界），每支
  恰一根。输出增列 `dpValve` / `dpPump`；
- **设备回归**：三算例 `tests/check_steady_equipment.py`（闭式二次根
  1e-9）：

  | 算例 | 校验内容 |
  |---|---|
  | `steadyEquipmentValve` | 定压差 + K=500 阀：q 对层流二次根 1e-9，dpValve = Kρq²/2A²（占压降 2.5%）|
  | `steadyEquipmentPump` | 闭环路（两端等压）：工作点 = 系统曲线∩扬程曲线二次根（q*=3.085629e-4、H=24.1315 m），dpPump = −ρg·H(q*) 与摩阻逐 Pa 相等（回路恒等式），净 dp = 0 |
  | `steadyEquipmentNetwork` | 源率驱动串联网（阀→自由节点→泵 + 10 m 高程→定压终端）：p_J1=72750 / p_W1=173511 对闭式链 1e-6；泵几乎顶平高程（p_J1 < 终端压力）——泵头符号错误的符号测试（差 ~4 bar）|

- **耦合边 `coupledResistance`（把一条管交给三维求解器）**：pipe 字典可带
  `coupledResistance <R>`（Pa·s/m³）。它为正的瞬间，本管的 Churchill 摩阻、重位、阀、
  泵在 `pipeLossAt` 里全部被 `R·q` 取代——这条边不再有自己的 1D 物理，它的压降由伙伴
  报来；Newton 的电导仍是对本管 law 的中心差分， law 此时恰为线性所以斜率就是 `1/R`，
  而那个"层流上限"也从 HP 几何值换成 `1/R`（HP 是几何的性质，不是这条边的物理）。
  加载时校验三条：R ≥ 0、该管不得是黑油/PVT 管（`march` 逐段算压降，不是一个系数）、
  不得同时带阀或泵（否则"这条边的压降"语义不再唯一）。回路见 5.18；
  瞬态内核的同一钩子是 `TransientPipeConfig::coupledResistance`，但它只替换本管摩阻——
  重位与区域惯性 `ρL/A` 仍留在 1D 的隐式动量行里，伙伴因此只需报"此刻抵抗的压降"（见 5.19）。

- **表驱动单相物性（PVT 供给 ρ(p)/μ(p)）**：`PvtFluid`（严格递增压力
  网格 + ρ/μ 表，节点间分段线性、表端钳位；NaN/p≤0 取首条，Newton 早期
  迭代永不产生坏物性）+ `solvePipeFlowPvt`——与黑油 march 同构：24 段
  局部机械能平衡 `dp/ds = −ρ(p)g·sinFlow − f(Re(p))·ρ(p)v|v|/(2D)
  + ρ(p)g·H/L`，**质量流量守恒**（`w = ρ(pUp)·q`，局部 `v = w/(ρA)`，
  气体沿程膨胀、出口速度 > 进口 q/A）、阀 K 均摊到段、Re/f 取出口值；
  流量由**单调律的扩括 Illinois secant** 反解（loss 无界增长、静水头有界
  → 恰一根，无需黑油的取最大根技巧；残差 1e-12·pUp）；电导沿用两相的
  单侧深探（静水头膝下 Q=0 平台，中心差分会把 Newton 指回死区），死侧
  探过膝。**平表严格退化为常物性管律**（每段复现同一隐式 law，单位测试
  对 `solvePipeFlow` 1e-9）；等温气柱（ρ=p/c）静水闭式 `p(h)=p₀exp(−gH/c)`
  在 Euler 截断 `x²/(2N)` 内成立。网络分派支路 `hasBlackOil > hasPvtFluid
  > 常物性`；字典 `fluid { pvt { pressure(...); rho(...); mu(...); } }`
  （等长 ≥2、正值、严格递增，blackoil/pvt 共存即 fatal）。回归
  `tests/cases/steadyPvtColumn` + `tests/check_steady_pvt.py`：checker 从
  case 自身的表出发复算 24 段推进，闭合报告的节点压力至 2.1e-12 rel、
  出口速度/Churchill f 至 ~1e-12 rel，气柱膨胀（出口速度 > q/A）与
  可压缩性可见性（dp 偏离进口密度常物性律 2.8%/0.6%）双断言。

环网/环路不需要专门的环路迭代（对标 Pipesim 环路）：Newton 节点压力
表述对任意环状拓扑直接适用——闭环回归在 `networkTests`（层流环，注入在
环上分流 72.2%/27.8%，对 2×2 Cramer 闭式 6.5e-11；残差地板 ~2e-10 相对
由稠密求解噪声决定，tol 取 1e-9）。黑油模式下物性由 §3.6
按 (p) 实时求取，单相模式物性由 case.dict 常量或上表 `fluid { pvt }`
供给。井筒热力耦合见 §4.3。

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
- **管件/阀门**：每管可带局部损耗 `dp = K·ρ|v|v/2`（与稳态内核同一 K 约定），
  并支持线性关阀/开阀时间表——`valveClose{ start; end; closedK; }` 把 K 从
  开启值线性升到 closedK，`valveOpen{ start; end; }` 再线性降回。时间表在
  **隐式新时刻 t+dt** 取值、牛顿步内保持常数；残差与雅可比分别加
  `K·ρ|Q|Q/(2A²)` 与 `K·ρ|Q|/A²`，残差缩放分母计入 `|dp_valve|`。时间表对
  t 连续，界点两侧 K 相等不跳变；
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
`pipes{}` 增 `initialFlow/valve/valveClose/valveOpen`，
新增 `time{ dt; endTime; reportEvery; }`）。

四个已验证的物理回归（`tests/check_transient.py`）：

| 算例 | 校验内容 |
|---|---|
| `transientBlowdown` | 罐排气：同时刻隐式 Euler 2×2 闭式递推（1e-9）+ 连续包络 `exp(−t/RC)`（2%）|
| `transientRlc` | 两柔性节点 + 惯性管 = RLC 欠阻尼振荡：首个流量峰的时间与幅值（3%）、≥4 次符号翻转、封闭系统液量守恒 |
| `transientSlug` | 段塞推进：串联管层流流量（1e-3）+ CSTR 级联递推（1e-9）|
| `transientValve` | 关阀/开阀时间表：代入节点柔度后每步只剩单个单调标量方程，二分到机器精度精确复现 (Q, p, α)（1e-9，时间表按 `t += dt` 与内核锁步累计——打印时刻仅 12 位有效数字，斜坡上会引入 dK~1e2 的舍入误差）+ 关阀泄漏 6.2e-8 m³/s、罐压冻结（漂移 ~1.0 Pa）、重开恢复排液 q(20)=3.385e-3、全程层流 Re=1229 |

尚未实现（M2 后续增量）：机理型段塞追踪（液膜/液塞双区制模型）——段塞
首增量（入口速率时间表驱动的段塞体检测/统计与合并算例）已由 §4.3 落地；
液膜夹带已由 §4.3 的
entrainment 增量落地；非等温与成核/生长动力学已由 §4.2 热力耦合与
§4.3 水合物动力学落地）。

### 4.3 分布式两流体瞬态内核（M2 第二增量，OLGA 主干）

`src/1d/include/cf/1d/TwoFluid.h` + `src/1d/src/twoFluid.cpp`。与 §4.2 的
集总节点容腔模型互补，这是**沿管长分布**的漂移流两流体模型（不可压相的
三方程核心）：状态是逐网格元的持液率 α（气相 1−α）与随液相输运的
水相分数 ϕ：

- **质量守恒**：`∂α/∂t + ∂(α·v_L)/∂x = 0`，气相由体积守恒
  `α_g = 1 − α` 给出；不可压 ⇒ 混合体积通量 `v_m` 沿程均匀（等截面约束
  现为**管链级**——`pipes{}` 串接允许逐段倾角但需等径，见串接管链条；
  气相 EOS 压缩性已落地——逐段气量输运开启后通量不再均匀，见下）；
- **漂移流闭合**：Wallis 常滑移 `v_G = v_m + α·u_s` ⇒ 液相通量
  `F(α) = u_s·α² + (v_m − u_s)·α`（凸二次）；或 Zuber-Findlay
  `v_G = C0·v_m + v_d` ⇒ `F(α)` 退化为线性；或**动态滑移**（第三种，
  携带真两相动量交换的状态变量）——见下一条。相内速度由
  `v_m = α·v_L + (1−α)·v_G` 反解；
- **动态滑移闭合（model dynamic）**：逐段输运滑移态 `s = v_G − v_L`，
  由双相动量方程消去共享 `dp/dx` 后的降阶单压力形式驱动：
  `K_s ds/dt = (1/ρ_G − 1/ρ_L)(ρ_m·dv_m/dt + loss/Δx) − D·σ`，其中
  `K_s = 1 + α·α_g·Δρ²/(ρ_L ρ_G)` 是滑移惯性因子、
  `D = 1/(ρ_L α) + 1/(ρ_G α_g)`、`σ = C·s|s|` 是界面阻力（C =
  interfacialDrag，kg/m⁴），`loss` 为该元摩阻 + 重位压降。阻力项半隐式
  （过原点割线），无条件稳定且不动点恰为终端滑移
  `s* = sign·√(|drive|/(C·D))`；混合动量感受滞后一步的滑移惯性
  `α·α_g·Δρ·ds/dt`（稳态时消失，故代数闭合的稳态持液率被精确复现）；
  输运的通量函数变为逐面 `F(α) = s̄·α² + (v_m − s̄)·α`（s̄ 取相邻元平均，
  开边界取 ghost 滑移），CFL 逐面取 max；
- **全气相动量闭合（model full，分布式压力动量）**：第四种闭合，
  逐段求解气相动量而非降阶滑移近似——气相动量衡算
  `α_g ρ_G Δx·dv_G/dt = α_g Δp − fr_G − α_g ρ_G g sinθ Δx − σ Δx`，其中
  Δp 取混合动量行的逐元压力梯度（`fr + ρ̄ g sinθ Δx + ρ̄ Δx·dv_m/dt`），
  壁面摩阻按相分数比例分摊 `fr_G = α_g·fr` ⇒ 摩阻与 `α_g Δp` 里的
  `α_g fr` 部分恰好相消。代入 `v_G = v_m + α·s` 并展开气相物质加速度
  `dv_G/dt = dv_m/dt + α·ds/dt + v_G·dv_G/dx` 得逐段滑移 ODE
  `ds/dt = [(ρ̄ − ρ_G)(g sinθ + dv_m/dt) − ρ_G v_G dv_G/dx − σ/α_g]/(α ρ_G)`
  （对流项前差取流向、出口面零梯度——与出口通量同约定），σ 过原点
  割线半隐式，终端滑移 `s* = sign·√(α_g|R|/C)`、
  `R = (ρ̄ − ρ_G)(g sinθ + dv_m/dt)`。均匀稳态下对流项恒为零（稳态通量
  均匀 ⇒ ρ_G 与 v_G 均匀）——稳态锚同时是对流项消失的恒等式校验，
  瞬态前沿推进期对流项活跃（实测峰值 0.199 m/s²）；含对流惯性 BVP 的
  非均匀稳态（压差驱动三元组、串接管链水平段剖面）是后续增量。与
  model dynamic 共享同一套逐面通量机器（faceSlip、Godunov 选型、CFL、
  输运 clamp——`slipClosure()` 同时覆盖两者，model 0/1/2 路径逐位不变），
  仅滑移 ODE 不同；rate 驱动时 v_m 自构造即钉在注入率上 ⇒ dv_m/dt = 0，
  滑移 ODE 纯浮力 + 对流驱动；
- **Godunov 通量**：对二次通量函数取精确黎曼解——凸 F 落差为激波
  （Rankine-Hugoniot 速度 `s = ΔF/Δα`）、上升为稀疏波（声速点
  `−q_b/2q_a` 处取值），线性时按符号迎风。CFL 条件
  `dt·max|F'|/Δx ≤ 1` 违反即 fatal，格式单调 ⇒ α 严格保持在 [0,1]；
- **混合动量**：隐式 Euler `ρ̄L·(v_m − v_m^n)/dt = Δp − Σ fr_i −
  Σ ρ_m(α_i)g·sinθ·Δx`（柱惯性是单位面积质量 `ρ̄L`——未知量是速度而非
  体积流量，不除以截面积），摩阻用同一条 Churchill 曲线的均匀流模型
  （混合物 `ρ_m/μ_m` 体积加权，Re 逐元），固定压差端用标量 Newton
  （残差对 v_m 单调），速率边界直接施加 `v_m = vsl + vsg`，双封闭端
  强制 `v_m ≡ 0`（不可压约束）+ 静水压报告；动态闭合时上式左侧另加
  `−Σ α·α_g·Δρ·(ds/dt)·Δx`（滞后一步）；
- **水相分数**：`∂(αϕ)/∂t + ∂(ϕ·F)/∂x = 0`，按通量方向迎风，ϕ ∈ [0, α]；
- **温度方程（thermal{} 块出现时）**：各段混合温度（两相热平衡），
  `C_i dT_i/dt = H_in,i − H_out,i + k'·Δx·(T_wall,i − T_i)`，`C_i` 是相加权
  `ρc` 容量（`α ρ_L c_L + (1−α) ρ_G c_G`，气密度取本步输运后 EOS 值）；
  焓通量**骑各相自己的体积通量**（与本步输运共用的面液/气通量，滑移下
  液相不再以混合速度运动——混合 `ρc·v_m` 加权会对液相容量过度拖曳，
  滑移为零时两种形式逐位相同），各相按**自己的通量方向**迎风；气密度
  迎风与质量通量同侧（rate 边界取注入参考密度、定压边界取 ghost
  `p/(gasR·T)`、管内取迎风元密度）；`phaseRates` 入口的组成按
  `vsl/vsg` 构造，开放边界的 ghost 温度两侧都必须给——出口那侧只在
  回流时被用到，而回流是定压端的正常结果。壁温 constant 或 geothermal
  （`T = 地表温 + 梯度·垂深`，垂深 = inletDepth − s·sinθ，出地即 fatal）。
  温度隐式，不加新的 CFL；每步与输运共用 DenseSolve 稠密求解（段数
  ≤ 数百）。每帧报告 `thermal: qWall … qInlet …`（壁面换热与入口焓流
  之和，稳态能量衡算用）；
- **水合物动力学（hydrate{} 块，需 EOS 与 thermal 同时开启）**：逐段输运
  水合物体积分数 `h = V_h/(V_L + V_h)`（油为 1 − 水 − h 的惰性余量），
  与水相分数同一迎风规则、同一 clamp 语义骑液相输运。每步输运后按元
  求过冷度 `ΔT = T_eq(p_i) − T_i`，平衡曲线 `ln p = curveA + curveB/T_eq`
  （p 取本步输运后 EOS 压力、T 取当前段温）；生长（ΔT ≥ 0）对
  `dh/dt = k_g·ΔT·(ϕ_w + h)` 做**隐式 Euler（对自由水一阶、对水汇隐式）**
  `h⁺ = (h/dt + k_g·ΔT·(ϕ_w + h))/(1/dt + k_g·ΔT)`，并以气源封顶
  `h⁺ ≤ h + max((m_G − trace)/(gasContent·α·Δx), 0)`；分解（ΔT < 0）为
  `h⁺ = h/(1 + dt·k_d·|ΔT|)`；随后 1:1 体积移动 `ϕ_w⁺ = ϕ_w + h − h⁺`
  （保体积约定）、按化学计量从 EOS 气量扣气 `m_G −= gasContent·α·Δh·Δx`
  （分解放回）。潜热反馈进温度方程：`rhs += latentHeat·α·Δh·A·Δx/dt`
  （生长放热、分解吸热）。应用输出逐元追加 `h …`；
- **气相 EOS 压缩性（fluid gasR > 0 开启）**：气相密度成为逐段输运状态
  ——气相质量连续 `d(ρ_G α_g)/dt + d(ρ_G α_g v_G)/dx = 0` 以单位面积气量
  `m_G = ρ_G(1−α)Δx`（kg/m²）为状态变量，气相体积通量
  `G(α) = v_m + (s̄ − v_m)α − s̄α²` 对三种闭合都是 α 的二次/线性函数且
  `G(1) = 0`，复用同一 Godunov 精确黎曼解与同一 CFL 界（凸凹与液相互为
  镜像，激波/稀疏波选择相反）；质量通量按体积通量方向迎风——
  `phaseRates` 注入字典参考密度 ρ_G，`fixedPressure` 边界的 ghost 密度由
  **绝对**压力 `p/(gasR·T)` 构造（EOS 关闭时既有 gauge 语义不变），管内取
  迎风元密度。状态方程 `p_i = ρ_G,i·gasR·T_i` 逐段报告（rate 驱动管线的
  热力学压力由注入质量决定，与水力压降剖面无关——稳态体积守恒把逐段
  ρ_G 钉在注入参考密度上，水力压降剖面由动量行经 pnode 另行报告；
  集总闭合不携带局部压力梯度，**分布式压力动量首增量已由 §4.3 的
  model full 落地**——逐段滑移 ODE 的驱动项携带混合行压力梯度）。逐段
  密度反馈进混合 ρ̄（摩阻/重位/柱惯性/滑移驱动的 inv(Δρ)、K_s、D 全部
  逐元化，滞后一步）；
  EOS 关闭（gasR = 0）时所有路径逐位回到常数 ρ_G，既有 7 个回归不动。
  封闭元保护：α→1 时气量呈指数泄空而气相体积线性关闭，数值保护为
  残气冻结密度、非残气 fatal；
- **串接管链（pipes{} 块，井轨迹首增量）**：`pipe{}` 单段之外，
  `pipes{}` 具名子块按出现次序串成管链，逐段 length/diameter/roughness/
  angleDeg/sections；逐元 Δx、sinθ、相对粗糙度与累计垂降（geothermal
  壁温按元心垂深取值）随链生成，段交界就是普通 Godunov 内部面（面滑移
  取相邻元平均，无需特殊语法）。等径校验：集总混合动量只携带一个速度，
  链必须等径（ctor 违反即 fatal）；总长 = Σ 段长，`initial fills` 的
  from/to 按链入口弧长定位，CFL 取最小元长。斜井/水平混合下各段独立
  稳态根由该段逐元 drive（摩阻 + 重位）决定，交界由**驻激波**衔接——
  交界通量恰为 vsl，晕圈只有几个元；
- **分支管网（network{} 块，Y/T 节点 + 双节点跨管）**：`pipes{}`/`junctions{}`/
  `boundaries{}` 三个子块把多根 `TwoFluidPipe` 拼成管网——每根管接一或两个
  节点（单节点管：该端无边界，另一端携带边界：定压/速率；速率驱动的
  远端必须节点在出口侧且总注入率 ≥ 0；**双节点管 = 跨管 span**，两端
  均无边界，由两端节点压力夹持驱动），节点端以 `pipe.inlet`/
  `pipe.outlet` 引用。每个时间步一步解全网再推进各管：**管网级全局
  Newton**（solveNetwork）的未知量是 [自由管（定压远端或跨管）混合速度
  …, 各节点压力 p_J …]，行 = 逐管动量行（残差 = 集总动量
  `ρ̄L(u − uⁿ)/dt + fr(u) + 重位 − (p_endIn − p_endOut)`，p_end 取定压
  远端或所接节点压力；跨管行同时耦合两个 p_J，解析导数 ∓1）+ 逐节点
  体积衡算行 `Σ s_i A_i v_i = 0`（DenseSolve 求解，收敛 |r| ≤ 1e-12·scale，
  60 次迭代）；单节点网络的矩阵与逐节点解逐位相同（twoFluidJunction
  回归作位兼容证明）。动态闭合与负注入暂不支持（fatal）。收敛后**来流
  混合**：持液率按体积流量加权、水相/水合物按液相体积流量加权、温度按
  热容流量 `A(|F|ρ_L c_L + |G|ρ_G c_G)` 加权、气密度按气相质量流量加权
  （EOS 下节点气质量衡算精确闭合）；相向流动的支路端 fatal。随后改写
  节点端边界：**供给支**（s·v > 0）得到混合态定压 ghost（面通量自然从
  边界元迎风），**被供支**得到 `phaseRates`——`vsl/vsg` 取该管自身通量
  函数在混合持液率处的分配（F + G = v_m 对一切闭合精确保住体积衡算，
  均相闭合下相分配也精确）；节点混合密度经 `rhoGInjected` 注入 EOS 质量
  通量与导热系数。跨管由上游节点供给（入口 phaseRates）、向下游节点
  供给（出口 ghost）；混合只读元端状态、改写只写边界，无次序风险。
  串接语义不受影响：逐位复用同一 momentumResidual/momentumScale
  （serial Newton 与管网 Newton 共享算术），serial 7 个边界组合校验仅对
  单管生效；
- **液膜夹带（entrainment{} 块）**：液相按元分为壁膜与液滴两部分——
  液滴分数 `d`（输运层语义，闭合/动量/物性不动）：液滴含量
  `qd = α·d` **骑气相**，面通量 `qd·v_G` 按气相体积通量方向迎风
  （gvol = 0 的面零通量），壁膜通量 `F − qd·v_G` 补全液相总通量 F——
  α 的输运完全不变（膜+滴+气仍整体携带 v_m），分裂是线性一致的。每步
  输运后按元弛豫 `dd/dt = (E_eq − d)/relaxTime`（隐式 Euler，
  `E_eq = clamp((v_G − vCrit)/vRef, 0, 1)`，v_G 取本步输运后状态的
  gasVelocity），α ≤ 1e-9 的元跳过；CFL 另加 max|v_G| 覆盖。边界语义：
  `phaseRates` 端注入 `entrainedFraction ×` 施加液相率（与率一致的组成
  下恰等于含量×气速），`fixedPressure` ghost 端按 gvol 符号注入
  `α_ghost·E·v_G_ghost`（入流）或带出元内容（出流）。管网节点混合的
  液滴分数按液相体积流量加权并改写被供支的 entrainedFraction。E_eq 弛豫
  与入口注入共同作用下稳态**不是均匀 `d*`** 而是沿程发展剖面——面衡算
  与隐式交换合成的线性递推
  `d_i = (v_G* d_{i−1} + (Δx/τ)E_eq*)/(v_G* + Δx/τ)`、
  `d_0 = (E_in·vsl + α*(Δx/τ)E_eq*)/(α*(v_G* + Δx/τ))`；
- **入口速率时间表 + 段塞体统计（slug{} 块 + boundaries.inlet.schedule）**：
  `schedule` 按时间表逐点改写入口的标量 vsl/vsg（每步在区间中点
  `t + dt/2` 处判定——t 按步累加会有 ~1e-12 漂移，中点语义对漂移鲁棒，
  表点落在整步边界时触发步与机器精度一致）；每步输运后
  `slugBodies()` 扫描 `α ≥ alphaThreshold` 的连续元段、汇出总长
  `≥ minLength` 的段塞体 `{rear, front, length, inv}`（inv = ΣαΔx 为
  段塞体持液量）并计数——这给合并/破碎提供了逐帧可观测量，合并算例
  直接断言两个段塞并成一个的时间；
- **应用** `chflow-1d-twofluid <case.dict>`：输出 `t … vm … dp …` /
  `node <j> p …` / `cell <i> alpha … water … vsl … vsg …`（thermal 时
  追加 `T …`，动态/full 闭合时追加 `s …`，EOS 开启时追加 `p …`，水合物开启时
  追加 `h …`，夹带开启时追加 `entr …`）+ `thermal: qWall … qInlet …`；
  段塞统计开启时追加 `slug count N` 与
  `slugbody k rear … front … length … inv …`；
  `network{}` 算例改发
  `t … net` / `junction <名> p … net …` / `pipe <名> vm … dp …` /
  `pnode <名> <j> p …` / `pcell <名> <i> alpha … water … vsl … vsg … [T/p/h/entr]`。

十七个解析回归（`tests/check_two_fluid.py`，均把 Godunov 递推按 1e-9 复现）：

| 算例 | 校验内容 |
|---|---|
| `twoFluidFront` | 液柱推气：激波位置 = 液量积分 = `40 + v_m·t`（1e-9，与 u_s 无关）、水相前沿同速、剖面 TVD 单调；注：Wallis 闭合下 `F'(0.5) = v_m` 与激波速结构简并，一阶格式前沿按 √t 扩散 |
| `twoFluidSteadyHoldupWallis` | 稳态分层持液率 = 二次方程 `F(α*) = vsl` 的物理根（1e-9）+ 压降逐项重算（Churchill 摩阻 + 30° 重位）|
| `twoFluidSteadyHoldupZuber` | 线性闭合的稳态根 `α* = (vsl − v_m(1−C0) + v_d)/(C0·v_m + v_d)` + 水平摩阻压降 |
| `twoFluidSettling` | 封闭水平管重力分异：双激波各以 `u_s·α0` 精确推进（±2.5Δx）、液量守恒 1e-9、壁面平台 α→1/气袋 α→0 |
| `twoFluidPressureLaminar` | 定压启动：独立二分解隐式 Euler 动量 ODE（1e-8）+ 指数包络 `u* = ΔpD²/32μL`、时间常数 `ρD²/32μ` = 0.125 s（3%）|
| `twoFluidDynamicSlipSettling` | 动态闭合重力分异（封闭竖管）：滑移映射 + 输运全耦合逐元复现（1e-9）、液量与重位压降精确守恒（`dp = ρ_m(ᾱ)gL`，实测每帧 100027.83 Pa）、终态下半 α > 0.99 / 上半 α < 0.01、界面钉在中点、终端滑移 `s* = √(drive/(C·D))`（实测滞后 2.0e-4 相对，容差 2e-3）|
| `twoFluidDynamicSlipHoldup` | 动态闭合稳态持液：耦合根 `(α*, s*)` 独立求出（F(α, s*(α)) = vsl + Churchill 摩阻驱动，实测 0.5529126041 / 0.6185778539），终场逐元对齐 1e-9、相分配 vsl = 0.4 / vsg = 0.6、压降 4437.1846 Pa 重算、末两帧全等（真稳态）|
| `twoFluidCompressibleHoldup` | 气相 EOS 压缩性：`twoFluidDynamicSlipHoldup` 同管线 + gasR 287。滑移 + 持液 + 逐段气量 `m_G = ρ_G(1−α)Δx` 三态全耦合逐元复现（1e-9）、逐段 `p = ρ_G·gasR·T_0` 恒等式、全局气量衡算 `Σm_G(t) = Σm_G(0) + ∫(入−出)`（实测终值 89.417 kg/m²）、稳态气量平衡钉住 ρ_G* = 20 精确（EOS 压力 1.722e6 Pa）且 `(α*, s*)` 与不可压回归同根 0.5529126/0.6185779、压降 4437.1846 Pa 不变 |
| `twoFluidWellChain` | 串接管链：竖直 riser（20 m/10 元，sinθ=1）+ 水平 flowline（20 m/10 元）动态闭合，交界共享 Godunov 面。全耦合逐元复现（1e-9）；两段独立稳态根——riser (α*, s*) = (0.7580622562, 1.9523156870)、flowline (0.5529126041, 0.6185778539)——bulk 元分别 1e-3 / 1e-9 对齐（交界晕圈 riser 第 7–9 元为驻激波边界层，豁免）；逐元相分配恒等式 `vsl + vsg = v_m`（1e-12）；压降按复现稳态像逐项重算（Churchill 摩阻 + 重位，实测 169374.682443 Pa）；末两帧全等（真稳态）|
| `twoFluidHydrateSealed` | 封闭静置水合物：wallis 滑移 0 使一切输运通量恒为零，各元跑纯逐元递推（生长对自由水隐式 + 气源封顶 + 1:1 体积移动 + 化学计量气汇 + 带潜热的对角热步）逐帧复现（1e-9）；闭式终态 `h* = water0 = 0.5`、`water* → 0`（实测 2.2e-16）、`T* = 壁温 280`、`p* = m_G*·gasR·T_wall`，`m_G* = 200 − 122·0.5·0.5·2 = 139 kg/m²` ⇒ `p* = 1.117004e7 Pa` 精确；末两帧全等（转换壁面受限，尾部指数衰减后弛豫到壁温）|
| `twoFluidHydrateFlow` | 流动水合物（rate 驱动 + EOS + thermal + 动力学全开）：整步全耦合逐元复现（1e-9）——Godunov 输运 α（wallis F 凸迎风取左）、水/水合物骑同一液相通量同 clamp 语义、气量按 wallis G 迎风 + 化学计量气汇、相通量导热的隐式温度解；恒等式：逐元组成 `w + h = 0.25`（1e-12）、相分配 `vsl + vsg = v_m`（1e-12）、发射 EOS 气量对复现 m_G（1e-9）、稳态气衡算 `ρ_G 注入·vsg − ρ_G,out·vsg,out = gasContent·vsl_out·h_out`（实测 8.4795689623 对 8.4795689622 kg/m²/s）、Churchill 压降逐项重算（4774.6772 Pa）、壁面换热求和、稳态能量衡算 `入口焓 + 壁面 + 潜热 = 出口焓`（1.7736709e5 W，潜热 = gasContent 消耗率×面积×latent）；水沿程单调耗尽/水合物单调累积；末两帧全等 |
| `twoFluidJunction` | 对称 Y 节点合流（branchA/branchB 定压 200/300 Pa → 节点 J → outlet 0 Pa，均匀持液 0.5 + wallis 滑移 0 ⇒ 混合物性处处常数）：校验器独立积分消元后的 2×2 层流隐式 Euler 系统（每步，rtol 1e-4 覆盖 Churchill 与层流 ~1e-6 偏差），并逐元复现节点混合 + Godunov 水相输运（1e-9）+ 隐式迎风温度解（1e-7——轨迹速度携带层流偏差）。逐帧：节点体积衡算 net ≈ 0（Newton 衡算行，1e-12）、p_J/三管 vm 对齐轨迹、发射 dp = Churchill 摩阻和逐项重算（1e-9）、α 恒 0.5（1e-9）、`vsl + vsg = v_m`（1e-9）。稳态闭式锚：`p_J* = Σ(A_i/K_i)p_i / Σ(A_i/K_i) = 60.0375234522 Pa`（K_i = 32μ̄L_i/D_i²）、`u_A* = 0.0349556635`、`u_B* = 0.0599306884`、`u_C* = 0.0421717120` m/s（Re 71–129 层流）、三管损失行进节点压力收于 p_J*；来流混合闭式：`φ_J* = 0.2473580247`（液相体积流加权，branchB 注入 0.1 前沿冲刷）、`T_J* = 312.6320988 K`（热容流加权，branchB 注入 320 K），前沿 ~713 s 清空几何，末端帧对齐 1e-10 |
| `twoFluidNetworkSpan` | 双节点跨管 + 管网级全局 Newton：srcA/srcB 定压 200/300 Pa → 节点 J1 → 跨管 mid → 节点 J2 → sink1/sink2 定压 0 Pa（5 管 7×7 系统：3 动量行 + 2 衡算行，mid 两端均接节点、动量行被两端 p_J 夹持）。校验器独立积分 7×7 层流隐式 Euler 系统（每步，rtol 1e-4 覆盖 Churchill 与层流偏差），并逐元复现两级节点混合 + Godunov 水相输运（1e-9）+ 隐式迎风温度解（1e-7）。稳态闭式锚（节点导纳 `c_i = A_i/K_i`，K_i = 32μ̄L_i/D_i²）：`p1* = 141.9280795715 Pa`、`p2* = β·p1* = 12.2417750574 Pa`（β = c_m/(c_m + c_1 + c_2)，p1* = 500·c_s/(2c_s + c_m(1−β))）、五管 `u* = 0.0145035/0.0394785/0.0539820/0.0171978/0.0152869` m/s（Re 30–110 层流）、pnode 恒等式（srcA/srcB 出口 = mid 入口 = p1*，mid 出口 = p2*，sink 入口 = p2*）；来流混合闭式：`φ1* = 0.2074690`（液相体积流加权，branchB 注入 0.1 前沿冲刷）、`T1* = 314.62655 K`（热容流加权），J2 由跨管单独供给（混合态 = mid 末端元态）；前沿 ~700 s 清空几何、级联瞬态尾（srcB → J1 混合 → mid → J2 → sink 逐级上风翻转 + 隐式热尾，e 折叠 ≈ Δx/u 秒/级）到 endTime 2500 s 才衰减到位；发射 dp 从内核 vm 做纯损失行进重算（内核 Churchill vs 校验器层流 ~1e-9 相对偏差经 dp ∝ u²f 放大，sink2 dp = 全额 p2* 最敏感，vm 物理由 1e-4 轨迹行独立覆盖）|
| `twoFluidEntrained` | 液膜夹带首增量（rate 驱动水平管 + wallis 滑移 0.3 + entrainment vCrit 0.2/vRef 0.6/τ 8 s/初值 0.9，入口注入 10% 夹带）：整步逐元复现（1e-9）——α 的 Godunov 输运不变、液滴含量 `qd = αd` 骑气相按 gvol 符号迎风（gvol = 0 面零通量）、rate 入口注入 `E_in·vsl`、隐式 E_eq 弛豫（初值 0.9 → 下游沉到 E_eq，双向交换都被触发）；恒等式 `vsl + vsg = v_m`（1e-9）；稳态闭式锚：`α* = F⁻¹(vsl) = 0.5598164906`（1e-9）、发展剖面递推 `d_0 = 0.1492899147 → d_19 = 0.6026486661`（`E_eq* = 0.6132415786` at `vG* = 0.5679449472`，1e-9——入口注入 0.015 与内部输运 `α*d*v_G*` 不平衡，均匀 `d*` 不是不动点）、水相冲刷到注入 0.25、Churchill 压降逐项重算（实测 3585.64 Pa）；末两帧全等 |
| `twoFluidFullMomentum` | 分布式压力动量（含气相对流惯性；rate 驱动竖直 riser 20 m/20 元 + EOS，model full C = 1000，重气 ρ_G 500/gasR 1000 使滑移 ODE 真惯性 rate·dt ≈ 0.5）：滑移 + Godunov 持液 + 水相 + 逐段气量四态全耦合逐元复现（1e-9，瞬态前沿期对流项 vG dvG/dx 活跃、峰值 0.199 m/s²）；耦合不动点 `s*(α) = √(α_g α Δρ g sinθ/C)`、`F(α*, s*) = vsl` 独立迭代（实测 α* = 0.7760652210、s* = 0.9231137603——均匀稳态对流项恒为零，锚兼作该恒等式校验）、体积守恒钉住 ρ_G* = 500 精确（EOS 压力 1.5e8 Pa 均匀）、相分配 vsl = 0.15/vsg = 0.25 逐元、气量衡算（实测终值 2239.35 kg/m²）、压降 = Churchill 摩阻 + 重位逐项重算（实测 179140.75 Pa）、pnode 自出口锚 1.5e8 逐元行进、水相冲刷到注入 0.3（浮力回落使瞬态大半程液相通量为负——前沿需 ~350 s 全冲刷，endTime 400）；末两帧对齐 1e-9 |
| `twoFluidSlugPulse` | 入口速率时间表 + 段塞合并（Wallis 闭合，三段速率表 0.05→0.15@30→0.05@57→1.2@65，dt 0.05/endTime 120）：每步按区间中点判定调度点（对 ~1e-12 的 t 累加漂移鲁棒，表点 30/57/65 均落在整步边界——先做 1140 步 `t += dt` 证明 `t_ = 56.999…999` 漂在 57 之下、逐点 `t <= t_` 判定会晚一步触发）；vm 逐帧钉在调度率（1e-12）；α + 水相逐帧复现（1e-9）；基线窗均匀 `α = F⁻¹(0.05; v_m 0.6) = 0.23166247903554`、段塞数 0；脉冲 1 前沿 = 激波速 `S1 = 0.50739027146868·(t−30)`、段塞数 1、TVD 单调；凹陷期尾部按特征速 `F'(0.3) = 0.4` 后退、前沿 `S1P = 0.40739027146868`；脉冲 2 入口激波 `S2 = 1.73608724052384` 追上前段塞——65.5 s 数到 2、66.0 s 合并为 1，合并后尾部以 `S1PP = 1.5574545` 侵蚀、前沿 `16.9587 + 1.5575(t−65)`，到达出口前整场均匀至 `A_HI2 = 0.7406029237394383`；液量积分 E 的闭式检查点（30/56.5/64.5/65.5 s，斜率取离散伸缩窗精确值 SLOPE1/SLOPE2A）；段塞体统计 rear/front/length/inv 逐帧复现（1e-9）、终态整管为一段塞体；pnode 逐元 Churchill 行进 + 压降 |
| `twoFluidPressureFull` | 定压驱动 model-full 稳态边值问题：竖直不可压 riser 两端定压 3e5/1e5 Pa，**不给任何速率**——混合行 Newton（隐式 Euler、持液冻结、滞后滑移惯性）从驱动压差解出 v_m，相速率在边界**涌现**。校验器逐步复现 Newton（与内核同初始点/同中心差分/同 1e-12·scale 终止）+ 滑移 ODE（对流项 + 原点割线半隐式阻力）+ Godunov 水相（1e-9，含瞬态对流项峰值诊断）；稳态锚：滑移闭式 `s* = √(α_g α Δρ g sinθ/C)`（对 ρ_G 不敏感，实测 1.4646517）、`v_m*` = 稳态行二分根（实测 4.92）、涌现入口速率 `vsl = F(0.5, s*) = 2.0938` / `vsg = G(0.5, s*) = 2.8261`、ghost 组成即稳态持液 `α* = 0.5`（第二通量根需 v_m < s，此处排除）；BVP 闭合恒等式：入口锚起的损失行进恰好落在出口锚 `p[N] = 3e5 − 20·1e4 = 1e5`、`dp = Σ损失 = Δp驱动`（1e-9）；末两帧 = 精确离散不动点（1e-12）。参数注记：半隐式割线阻力需 |s_old| > 0 才接上，首步会冲到 `(ρL−ρG)g dt/ρ_G`——ρ_G 300/dt 0.05 把该过冲压在 `s*(0.8) = 1.31` 之下，启动平滑 |

三个温度回归（`tests/check_two_fluid_thermal.py`，断言离散闭式而非
连续指数）：

| 算例 | 校验内容 |
|---|---|
| `twoFluidThermalRelax` | 封闭管、恒壁温：每段独立集总容量，逐元对齐离散闭式 `T^n = T_w + (T_0−T_w)/(1 + dt·k'Δx/C_i)^n`（1e-9）——液段（C = 656.6 J/K → 316.511874 K）与气段（C = 6.75 J/K → 280.016706 K）弛豫速率相差四个量级，直接考相加权容量 |
| `twoFluidThermalExchanger` | 稳态单相流 + 恒壁温：离散稳态是精确递推 `T_i − T_w = (T_{i−1} − T_w)·r`，`r = ṁc/(ṁc + k'Δx)`（1e-9，逐元），且稳态恒等式 `Σ k'Δx(T_i − T_w) = ṁc(T_in − T_out)`（实测双侧 1023.5341 W）；连续 Ramey 指数 290.9027 K 只作参照打印（O(Δx) 壁面离散差）|
| `twoFluidThermalGeothermal` | 竖直关井弛豫到地温线：稳态逐元 `T_i = 288 + 0.03·d_i`（1e-9）、沿深度单调、终态 qWall → 0（2.3e-10 W）——关井井温测井的极限形态 |

注意瞬态寿命由**驻留时间**决定而非弛豫时间：扰动每穿过管程衰减
`r^N`，换热器算例要跑 20 个驻留时间（2000 s）才到稳态。

### 4.4 Python 绑定 `src/bindings` + `python/chflow`
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
vapour_density_table / vapour_fraction_table / load_flash_case`；节点分析侧
导出 `BlackOilFluid / IpRForm / IpR / WellLift / NodalControls / NodalResult /
NodalCase` 与 `solve_nodal_analysis / vlp_bottom_hole_pressure /
load_nodal_case`。C++ 侧的
`cf::Error` 自动转成 Python `RuntimeError`。Python 包 `chflow` 只是
再导出层。求解与闪蒸调用释放 GIL。

### 4.5 CLI `cli/`（Rust）
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
单相不可压定常回路（simpleFoam 的**全压力形式**，全部自研）：

1. 动量：矩阵**不含压力**，`div(phi,U) − laplacian(nuEff,U) = V(bodyForce −
   grad(pⁿ))`——`−V·grad(pⁿ)` 走源项（OpenFOAM 的 `solve(UEqn == −fvc::grad(p))`），
   解出的 `x` 是对当前压力迭代的物理响应。欠松弛通过对旧场混合实现：
   `Hfree = (1−α_U)·Uold + α_U·x`（不是对矩阵零初值混合）；
2. `HbyA = constrainHbyA(Hfree)`：速度给定的 patch 保留其法向通量，其余
   zeroGradient；
3. 预测通量 `phi = flux(HbyA)`（OpenFOAM 的 `phiHbyA`）——单元压力梯度
   不进预测通量，面速度才与 `SnGrad` 闭合成 Rhie-Chow 一致（实测：预测
   通量里带单元梯度时，渠道出口柱 Ux 塌陷 27%）。`adjustPhi` 仅在压力
   无任何 fixedValue（`needReference`）时做全局平衡，否则边界物理吸收
   失衡，强行缩放反而把自由边界通量反转；
4. **全压力** Poisson：`laplacian(rAUf, pFull) = div(phi)`，`pFull` 用算例
   自己的 p 边界条件——修正量钉在 0 的 `p'` 形式会把 `pⁿ` 的边界误差
   冻结进定点；
5. 校正：面通量按**未松弛**的解出压力闭合 `phi −= rAUf·SnGrad(pFull)|Sf|`
   （OpenFOAM 的 `pEqn.flux()`，公共件 `correctFluxes`），然后
   `p += aP·(pFullNew − p)`、`U = Hfree − aP·rAUr·(grad(pFullNew) −
   grad(pⁿ))`——deferred 拆分的两半用同一个 `aP`，定点就是离散动量平衡
   `aDiag·U − offdiag·U = V·bodyForce − V·grad(p)` 本身。

**可选的时间项（同一回路，不做 PISO 的分步投影）**：`SimpleControls` 带 `dt` 与上一层的
`uTimeOld` 时，动量方程额外装配 `fvm::ddt(UEqn, *uTimeOld, dt)`（隐式 Euler：`V/dt` 进对角、
`V·u_old/dt` 进右端），于是 SIMPLE 可以在**一个时间层内**迭代到该层的不动点，而不必把时间精度交给
PISO 的算子分裂。三者必须同时成立，缺一即 fatal：给了 `uTimeOld` 却没给 `dt`、给了 `dt` 却没给旧场、
以及 `schemes{ ddt … }` 不是 `Euler`（`src/3d/src/solvers/simpleSolver.cpp`）。默认值
（`dt = 0`、指针为空）对所有既有调用者为惰性，稳态回路与全部既有算例逐位不变。实际用法见 5.33：
耦合接口按时间层调用 `beginTimeLevel(dt, qPrev)` 快照一次 `u_old`（**层**拥有历史，不是每次 sweep），
`V/dt` 落在对角上同时把 `rAU` 兜住不超过 `dt`——实测每次块解的扫数因此比准稳态伙伴**更省**
（渠道 11.8 → 7.8~8.4 扫），代价转移到界面那一侧的趟数。

**压力欠松弛 `aP`（自适应步长限制器）**：`rAUr = V/aDiag` 是局部响应，
而通量闭合所需的整段压力斜坡由全域壁面介导——细长算例第一个 Poisson
把启动失衡解到 K 倍物理压力（100×20 渠道实测 K≈70，即拉长比），该步长
经 `grad(p)` 反馈回动量会以 ~K 的增益振荡。限制器把每次压力步长压到场
尺度的 15% 以内（线性区），随步长收缩向 `relaxP` 打开：
`aP = min(relaxP, max(pRelaxStart, 0.15·pScale/pStep))`。线性区内
sweep 间是速率 `|1−2aP|` 的定点二循环（`aP=0.02` 实测衰减 0.955 精确吻合）；
`pRelaxStart=0` 退回常数 `relaxP`。固定大 `aP`（≥0.05）会经入口角点局部
模态进入有界极限环，`relaxP` 取 0.02~0.03（`relaxU=0.9`）。

收敛判据：预测通量由 `HbyA` 单独构成、压力通量那一半由校正补上，因此
预测散度收敛到 HbyA/压力通量的固有失衡量（验证渠道 ~5e-3；这个量还随
边界上的**速率分布**变，同一渠道铺平入口时实测 4.8e-2、铺抛物入口时
2.0e-3，见 §5.34），
`continuityTol` 只在失衡较小的算例可达；**`|dp|` 分支（`pcorrTol`）才是
该结构下的真实外收敛信号**。两分支都过动量残差门 `momentumTol`
（`max|U−Uold|/max|U|`，含动量欠松弛与压力响应的总固定点变化）——
只看连续性判不了动量平衡：流向均匀的单列渠道所有 `phi` 恒 0，
第 1 次外迭代就读到 `continuity 0`，壁面剪切还没进入平衡（实测关掉
`momentumTol` 时 1 次迭代即"收敛"，`u_τ` 只有 0.554，对积分平衡差 11%）。
加门后体积力渠道 ~30 次外迭代收敛（k-ε 27 次、SST 39 次），`u_τ` 与
拖力平衡符合到 3e-7，与 PISO 750 步的定常解一致（Ub 10.96 vs 11.10）。

刻在实现里的经验结论（回归历史，勿回退）：

- **边界上的速度 prescribe 是逐面的，`HbyA` 必须逐面继承它**（`adoptFaceValuesFrom`）：
  `constrainHbyA` 原来只搬 patch 的那一个数，于是耦合界面把入口速率铺成一个形状时，
  方程被形状强迫、预测通量却仍是平的——边界那一列留下一个形状自身大小的失衡，任何
  迭代都消不掉（抛物入口实测地板 6.46e-2 → 3.6e-3，SIMPLE 与 PISO 同修）；
- **外停机的 `settled` 比较的是两次扫之前的连续性，不是历史最小值**：冷启动的前两个
  二循环会先跌破它最终停住的地板，从别的速率热启动也会从更低的位置起步，拿历史最小
  值当参照就在此后几百扫都拒绝停机（同一算例 400 扫不收敛 → 21~44 扫）；二循环正是
  这个循环自己的周期，所以是"同相对同相"的比较，而增长仍然表现为增长；
- **rAU 必须与压力方程、通量校正、U 校正三处使用同一算子**，校正后的
  `phi` 才能严格散度自由（公共件 `correctFluxes` 保证）；
- 旧 `p'` 修正形式靠 Jacobi 缩放（`rAU = α_U/aDiag`）隐式欠松弛压力来
  稳定，代价是 `p` 在 Jacobi 尺度上累加、物理压力远慢于速度收敛
  （`poiseuilleChannel` 停在 `|p'|=1.4e-4` 时约 1500 外迭代）。全压力
  形式直接用物理响应 `V/aDiag`，压力收敛到真实物理值（入口角点
  实测 1.226 ≈ 解析 1.2 Pa），同一算例 47 次外迭代收敛。

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
- `chflow-simple-foam <case.dict>`：不可压定常 SIMPLE，`simple{maxIterations
  continuityTol pcorrTol momentumTol nonOrthCorrectors relaxU relaxP
  pRelaxStart force}`，
  逐次外迭代打印 `continuity / corrected / |p'| / dU`，收尾 `SIMPLE done:
  iterations ... continuity ... momentum ... (converged|NOT converged)`；
  退出码即收敛标志；case 带 `energy{}` 时每扫描多解一个温度场、vtu 多写单元量
  `T`、收尾多报一次热平衡（见 5.11）；
- `chflow-piso-foam <case.dict>`：不可压瞬态 PISO，`piso{dt endTime
  correctors continuityTol linearTol linearMaxIter nonOrthCorrectors
  reportEvery force}`，步进循环 `dt = min(dt, endTime − t)`，逐步打印
  `t / uMax / cont / correctors` 与结束时第一条 x 截面的 `profile y u`
  供解析核对，`piso done: steps ... worstContinuity ...` 收尾；带 `energy{}`
  时报告行追加 `energy t .. Tmin .. Tmax ..`，收尾报热平衡（见 5.11）；
- `chflow-vof <case.dict>`：两模式——相分数输运与界面几何（速度由 case 给定
  而不是求解），或 CSF 动量模式（`surfaceTension{}` 冻结界面、驱动 PISO）
  （见 5.10），vtu 写 `alpha`/`nHat`/`kappa`，可选 `vtkInitial` 写同一组量的
  初始场；case 带 `reconstruction{}` 时每次报告多一行 `plic ...` 并额外写
  `nHatPlic`/`offsetPlic`/`areaPlic`/`centrePlic`/`kappaPlic`（见 5.12），
  `alpha{ circle{ ... fill exact } }` 用圆盘与单元盒的解析交面积初始化 α；
- `VtkWriter`：ASCII `.vtu`（单元中心场），ParaView 直接可读；它**保存**加进来的
  场而不是拷贝，所以每个场要活到 `write()` 之后。

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

### 5.10 多相流层 `multiphase/`（M3 首增量：VOF 相分数输运 + 界面几何；M4：CSF 与动量耦合、解出通量输运 α）
`include/cf/multiphase/VoF.h` + `src/multiphase/vof.cpp`，命名空间 `cf::vof`。
只依赖 fvc/fvm/fvMatrix：`chflow-vof` 是两模式 app——case 带 `velocity{}` 用**给定**
无散速度输运 α（M3），带 `surfaceTension{}` 冻结界面、驱动 PISO 读出它隐含的
压力-速度响应（M4）——两条路径各自单独验收。

**输运是 M-矩阵**：`ddt(α) + div(phi, α) = 0` 用隐式 Euler + 迎风，作用在**离散
无散**的 phi 上时，每行对角占优且行和 = V/dt、非对角全 ≤ 0，于是守恒与 α ∈ [0,1]
是构造性质——不裁剪、不做 redistribution、不加限制器。实测旋转 case 32 步后
`min α = 1.58e-39`、`max α = 0.999901`、相体积漂移 −5.1e-14 相对。
单元测试拿**一步的离散闭式解**对表：阶跃场在 θ = U·dt/dx 下，锋后不变、锋前第 k
格是 `(θ/(1+θ))^k`，这条尾部正是迎风写在下游干区的东西。

**线性速度 ⇒ 逐单元精确无散**：面心取值的 Σ_f U(x_f)·S_f 是 ∫div U dV 的精确值
（面是平面、U 线性），所以守恒不是运气。均匀平移实测 max|div(phi)| = 0（真的是
0），刚体旋转 U = ω ẑ×(x−c) 实测 1.08e-19。旋转有穿壁速度，因此闭合盒只能在
"相永远到不了壁"的前提下成立：入口给 0、出口带走 ~0。

**界面法向** `interfaceNormal()`：n̂ = ∇α/|∇α|，指向 α 增大方向（对 phase-1 液滴
朝内）；只有 |∇α| 达到全场最大值 `interfaceThreshold`（case 里 0.1）比例的单元
带界面，其余给**精确 0**——这个 0 就是告诉曲率算子"别读我"。阈值取相对值是因为
|∇α| 的量级由初始化方式决定（同一个界面的单元填充阶跃与 tanh 能差一个量级），
而界面**位置**不由它决定。

**曲率** `interfaceCurvature()`：κ_o = (1/V_o)·Σ_f (n̄_f − n̂_o)·S_f，且**只累加两侧
单元都带界面的面**。减去自身法向在连续介质里恒等（闭合单元 Σ_f S_f = 0），它带来
两件事：一致法向读 0，以及带外邻面直接退出模板。不做这个减法时，把带外法向置 0
会让边界面的插值法向减半，两侧各注入 ±|S_f|/(2V) = ±1/(2dx) 的假象——dx = 0.0125
上实测 ±40，足以淹没真值 −1/R = −6.25。修好后平面界面在**任意**网格上读
κ = 0（单元测试：带内带外 max|κ| 都恰好是 0.0），手工构造的径向单位法向满足
`|κ·r − 1| / (dx/r)² < 1.85` 的二阶律（模板是单元中心法向插值到面，误差是 r̂ 的
二阶导 1/r² 对比散度本身 1/r）。
两条推论：κ 只由**带内面**承载，所以界面沿自身方向要跨 ≥2 单元——单元填充的阶跃横跨
只有 2 格且沿面方向无变化，读不出曲率（80×80 上 140 个带单元、法向偏径向 15.4° 均值
/35.5° 最大、带均值 κ = −4.57 对 −1/R = −6.25、逐单元 |κ·r+1| = 2.9 均值）；tanh 平滑
过的阶跃能读，而对流不在乎阶跃被怎么抹开。
**误差由带的厚度（以单元计）决定，不由分辨率决定**：80×80 上 smooth = 0/1/2/3/4/6 dx
给出法向偏离 15.4/9.5/4.0/2.4/1.4/0.7°（均值）与 35.5/22.9/12.2/6.9/4.5/2.4°（最大），
带均值 κ −4.57/−5.04/−5.21/−5.82/−6.00/−6.22 对 −1/R = −6.25；而固定 smooth = 4 dx 在
40、80、160 格上给出 1.39/1.40/1.40° 与 κ 均值 −6.14/−6.00/−5.97。抹开几个单元就是**梯度
路径**的几何精度设定；sharp 的 α 由 §5.12 的几何界面重构接住——它只需要 α，不需要抹开，
其曲率由 §5.15 的截面几何直接给出并可用 `surfaceTension{ kappaPlic true; }` 喂给动量。

**对流后的各向异性扩散**：`|U| = omega r` 沿流向把界面带拉宽（412 → 699 单元），逐单元
法向相对**相质心**从初始的 1.40° 均值/4.49° 最大涨到 3.98°/19.6°，逐单元 |κ·r+1| 从
0.137 涨到 0.178 均值（最大 1.93），而带的**均值** κ 仍是 −6.40 对 −1/R = −6.25。checker
因此分两处断言：解析初始场用 tight 带宽（方向 2.5°/6°、逐单元曲率 0.2/1.0），对流后只断言
"仍然是圆的"（8° 均值/30° 最大）。参照系必须跟着液滴走——用**初始**圆心量对流后的法向，
会把整个 0.4 rad 的转动算到界面头上（读出 19.7° 均值/76° 最大，那是轨道角速度不是格式
耗散）。逐单元几何之所以放在初始场上严格断言，是因为那是唯一有解析真值的场。

**守恒的符号边界**：平移 case 的质量平衡是有界非零——入口喂 U·A·t = 0.002，场
只涨 0.00199999973，缺的 2.72e-10 是出口带着 α = 1.19e-6 的干尾流出去的，
上界 U·A·t·α_outlet = 2.38e-9。它不随 solverTol 变化（1e-13 也一样），是物理泄漏
而不是求解误差，checker 因此断言**带符号的界**而不是 0。
另一个"不是最后一位"的实测：α 的 y 不变性最大跨行差 1.00004e-16 = α=1 的 1 ulp，
因为迭代求解的扫描顺序让各行独立收敛。

**连续表面力（M4 首增量）** `surfaceTensionForce(sigmaOverRho, alpha, kappa, grad)`：
`f = −(σ/ρ)·κ·∇α`，符号从法向约定整个推出来：n̂ 指向 α 增大方向（液滴**朝内**），
κ = div n̂ = −1/R 为负；物理上表面张力把流体拉向曲率中心（也朝内），于是
`0 = −∇p + f` ⇒ `p_in − p_out = σ/(ρR)`，**内侧更高**。直观的 σκ∇α 符号正好相反
（会把液滴吹爆）。力天然局限在带上：κ 在法向为 0 处为 0，厚度由 ∇α 承载——同一
界面用 n 单元厚的 tanh 把同一个积分摊到 n 个单元上、峰值低 n 倍。
**进入 PISO 的通路**：`PisoControls::cellForce` 把每单元加速度加进动量右端。它由
`explicitOperator`（即 H(U)）随 **HbyA** 一起携带，因此通量看到的是与压力同一套
单元梯度/面法向梯度配对——表面力需要的正是这种配对：压力在动量写出的地方（单元）
与连续性闭合的地方（面）同时平衡它，只到达其一时另一个就会加速流体。`correctFluxes`
因此一行不改。
**静态液滴实测**（`vofDroplet`，σ=0.08、ρ=1、R=0.16 ⇒ 真值 σ/(ρR) = 0.5，闭盒
全无滑移 + 一块 p=0 锚定）：解出的 Laplace 压差 0.49511（−0.98%，与带自身曲率误差
同阶）、远场 −3.1e-6、两侧体相压强平到跃变的 0.63%/0.38%（斜坡只属于带）、解场的
镜像对称 |Δp| 0.006% 与 |Ux+Ux'| 0.17% of max|U|、离散表面力净和 = ∫|f| 的
2.9e-13（内力）、闭盒净动量 = ∫|U| 的 3.4e-7、解出通量 max|div(phi)| = 1.1e-14
（1 个 corrector 即收敛）、α 逐位不变、伪速度 max|U| = 0.0123 = √(σ/(ρR)) 的
1.7% 且增量递减（0.0020→0.0014，是旋涡在建立而不是失稳）。
单元测试给恒等式：**均匀 κ 上离散力成对抵消**——Σ V·f = −(σ/ρ)κ·Σ V·∇α，而离散
梯度是 Gauss 和，闭合域上 α=0 边界使 Σ V∇α 精确为 0，实测 3.5e-17 相对；同时断言
力恰好在带上的单元非零（带外恰好 0）、液滴受压向内（outward 精确为 0）、κ=0 ⇒ f=0。
**范围与两相耦合**：vofDroplet 是**单向**的（α 冻结，只读出压力-速度响应）。第二增量
把回路闭合：`surfaceTension{ advect true; }` 时每步先按 αⁿ 的力跑 PISO，再用**解出的**
φ 输运 α（`advectAlpha`），再从 αⁿ⁺¹ 重算 n̂/κ/力。PISO 通量只在线性求解器容差内无散，
所以两相耦合下的质量守恒不再是闭式——它变成**边界交换**：闭盒壁面 α=0 时交换为零，
`Σ(αⁿ⁺¹−αⁿ)V` 精确到舍入（实测 1.0e-13 相对）；有界性仍是 M-矩阵性质（CFL 内）。
**椭圆液滴振荡**（`vofOscillating`）把周期当作内核没被给的数来测，并且**裁决是哪条
定律**：n̂ 朝内约定下释放 a≠b 的椭圆，n=2 模振荡。2D 自由液滴公式 ω² = n(n²−1)σ/(ρR³)
（势流 + 运动学 + Bernoulli 对 Laplace 跃变，κ′ = (n²−1)ζ/R²）在这里**差 35%**——它假设
环境不携带惯性；本内核动量是单流体式（ν、σ/ρ 为常数），环境就是同一种流体。给外域
补上 φ_out = B r⁻ⁿ、两侧同配运动学条件后得 **ω² = n(n²−1)σ/((ρ_in+ρ_out)R³)**：环境把
自己的惯性加进分母（同一陈述的能量形式：内外两域动能 ρπ ω²ζ²R²/n 对表面能
σπ(n²−1)ζ²/2R）。匹配密度下 n=2 读 **ω² = 3σ/(ρR³)**；3D 自由球是 l(l−1)(l+2)=8、
3D 匹配球是 4（单位 σ/(ρR³)）——算例同时裁决几何与惯性。R = √(ab) = 0.15819
（解出的 Laplace 跃变落在它上面：dp 均值 0.502 对 σ/(ρR) = 0.5057），线性周期
T = 2π√(ρR³/(3σ)) = 0.8071 s，实测（aspect 二阶矩穿越钟）**0.8602 s，+6.6%**：振幅
(a−b)/(a+b)=0.15 的有限振幅软化（首四分之一周期 0.243 对 T/4=0.202）、迎风让界面
厚度在程内翻倍、阈值带低估端部曲率（t=0 带均值 κ −4.7 对弧长均值 −6.2），三者同向。
变形每周期衰减 ~×0.5——黏性自己只给 ×0.8（2D n=2 模 γ = 6ν/R² = 0.24/s），主导阻尼
是迎风加厚的带。首步 |U| 冲到 0.62（≈0.87×毛细速度 √(σ/(ρR)) = 0.711）：物理释放 +
`interfaceNormal` 阈值截断的带缘力尖峰（截断处离散 κ 跨 ~2/dx，少数带缘单元 |f| ~141），
t~0.2 内 ring down、不触及形状（aspect 从首个报告起就光滑）——与 vofDroplet 量到的
伪电流同源，此处叠在物理释放速度上。

**case 侧**：`velocity{ type uniform|rotation; U | omega center }`、
`alpha{ plane{ point normal } | circle{ center radius [smooth | fill exact] } | ellipse{ center radiusA radiusB smooth } }`、
`reconstruction{ [normal youngs|centroid|mixed|fit] [transport alpha|geometric] [bandTolerance] }`（见 §5.12 与 §5.14，缺省不重构）、
`patches{ name{ alphaType alpha } }`、`time{ dt endTime reportEvery interfaceThreshold }`、
`schemes{ ddt Euler; div; grad }`（ddt 只接受 Euler）、`solverTol`/`solverMaxIter`、
`output{ vtk vtkInitial }`；逐步报告行
`vof t .. mass .. massErr .. min .. max .. cx .. cy .. x50 .. kappa .. iface N .. solver I .. residual R`，
其中 `x50` 是相质心所在行上 α 的 0.5 穿越位置（checker 与 vtu 独立复算互校）。开了
`reconstruction{}` 时每次报告前多一行
`plic <stamp>: band N cells area A readBack R kappa K`（法向取 `fit` 时行尾再加
`scatter S rounds N`：邻域平面拟合后仍不一致的距离，与最多用掉的轮数），vtu 多写
`nHatPlic`/`offsetPlic`/`areaPlic`/`centrePlic`/`kappaPlic`（动量模式的 vtu 同样带上这组）。

**动量模式（M4）**：`surfaceTension{ sigma rho nu [advect] [kappaPlic] correctors continuityTol }` 给出时
驱动 PISO 而不是输运——case 必须恰好含两者之一；`advect true` 接通两相耦合（解出的
φ 输运 α、每步重建力），缺省 false 冻结界面；`kappaPlic true` 把力的 κ 换成 §5.15 的
截面几何算子（需 `reconstruction{}`，两相耦合下每步随场重建）；`patches{}` 在同一块里加读
`U (0 0 0)`（无滑移）与 `p`（一块 patch 的 Dirichlet 锚定闭合盒的常数），pType
可写 `zeroGradient`；报告行换成
`vof t .. mass .. massErr .. min .. max .. cx .. cy .. dp .. pIn .. pOut .. umax .. [aspect] .. corr .. cont .. dpMax`，
其中 dp/pIn/pOut 分别是跃变与按 α>0.99/α<0.01 分离的两侧体相均值，advect 时多报
`aspect`（α 二阶矩轴比 √(Ixx/Iyy)，振荡的形状信号；checker 从 vtu 独立复算互校）；
vtu 额外写 `U`、`p`、`csf`。

### 5.11 能量方程 `fvm/EnergyCoupling.h` + `fvm/EnergySolver.h` / `solvers/energySolver.cpp`（M3：对流扩散 + 壁面热流 + Boussinesq 浮力）
case 带 `energy{}` 块时，SIMPLE 与 PISO 各自多解一个温度场；不带则这两个求解器
与 §5.5–5.7 完全一致（`IEnergyCoupling*` 是新的可选形参，缺省 nullptr）。

**方程与单位**：`dT/dt + div(u,T) − div(α∇T) = S`，
`α = κ/(ρcp) + nut/Pr_t` [m²/s]，
`S = volumetricSource/(ρcp) + ν|S_ij|²/cp + ε/cp` [K/s]
（分子耗散由 `viscousDissipation` 开、湍流释热由 `turbulentDissipation` 开，ε 由湍流层给）。
求解器本身保持**运动学**（压力方程里 ρ=1、`force` 是加速度），所以 `rho`、`cp`
只用来把 κ 与 W/m³ 的源折算成扩散系数和升温速率，浮力也以加速度进入动量。

**装配只复用已有 fvm 项**：`fvm::ddt`（隐式 Euler）、`fvm::div(eqn, phi.values(), div)`
——这里的 `phi` 就是求解器**闭合连续性那一个**通量，以及
`fvm::laplacian(eqn, alphaCell(), scheme, −1)`（取负号让 LHS 正定）+ `addRhs(V·S)`。
正因为 T 看到的是同一个 φ，离散焓平衡才不是近似而是恒等式：
`Σ Q_patches + ∫q dV − ρcp·Σ_{边界面} φ_f T_f = 0`，定常 case 实测残差
4.0e-12 / 2.3e-10 / 1.1e-10 / 9.9e-15 W（对 280–301 W 的量级）。面扩散系数用
`harmonicFaceGamma`，边界面取 owner 单元的 γ——这样 `patchHeatRate()` 报的就是
算子自己装配的那套系数，不是另一套。

**热边界来自同一个 `patches{}` 块**：`T 320` → fixedValue；`heatFlux q` →
**符号约定：正 = W/m² 进入流体**，app 按 Fourier 定律换成 `fixedGradient = q/κ`
（沿外法向）；`TType ...` 直接给类型；三样都不写的 patch 保持标量缺省
`zeroGradient`，即绝热——所以 case 只列出真正带热的面。
`output` 的 vtu 多写一个单元量 `T`，日志多报
`energy: Tmin .. Tmax ..`、逐面 `heatPatch <name> Q <v>`、以及
`energy: in .. source .. enthalpyOut .. balance ..`（这些数用 12 位精度打印，
`log::info` 缺省的 6 位会把平衡残差截成 0）。

**耦合次序**（这一条决定浮力能不能被压力回答）：SIMPLE 在 p/U 松弛**之后**跑
`correctEnergy(0)`，让它搭本次扫描 Rhie-Chow 一致的 φ，而浮力在扫描**顶部**读取，
于是力与温度属于同一次扫描；PISO 把 `buoyancyForce()` 并进 `force`，因此它随
HbyA 一起到达通量（与 §5.10 的表面力同一条通路），`correctEnergy(ctl.dt)` 放在
corrector 循环**之后**。定常欠松弛是手动的——先解（`eqn.schemes().relaxation = 1.0`
避免矩阵把自己混到零初值上）、再对着求解前的快照混合；瞬态直接取解（隐式 Euler
的步长就是解），所以 tight 的定常回归一律写 `relaxation 1`。

**两条实测出来的坑**（都写在算例头注释里）：
① **全壁绝热 + 无体源的定常热矩阵是奇异的**：Krylov 落在哪一级由迭代子空间决定，
不是由 case 决定（实测 T 整体塌到 0）。定常热 case 至少要有一个温度条件或一个源。
② **闭合盒 + 均匀力的压力-速度地板**：均匀力是投影求解能吃下的力，理论上盒中
什么都不该动、`dp/dz = f_z = β(T−tRef)|g|`。SIMPLE 会停在 continuity 1.05e-4、
uMax 1.5e-4（压力锁住的模态）；换成 PISO 跑 1000 步后 u_x 6.9e-12、
worstContinuity 2.5e-10，但 max|U| 仍是 **2.706e-4**、`dp/dz` 0.099972 对 0.098100。
关键是**对照实验**：把 `energy{}` 整块删掉、改由等价的均匀
`force (0 0 0.0981)` 驱动同一个盒，读到同一个 2.706e-4 与同样的两端压力——
这个地板属于闭合盒耦合本身，不属于能量方程。**给定均匀体力的那份已由 `p_rgh`
表述解决**（`closedBoxGravity` 算例，PISO/SIMPLE 两 app 共用 `pRgh` 开关）：
`p_rgh = p − force·x` 把均匀力精确吸收进变量本身——动量变成 `a = −∇p_rgh +
（只含 Boussinesq 偏差）`，闭盒的逐 patch 常数 BC 对 p_rgh 恰好是正确 BC
（真压在求解后按 `p = p_rgh + force·x` 重建）。实测同一只盒子：plain-p 蠕升
到 uMax 1.82663170862 m/s 且持续增长；p_rgh 模式 `uMaxAll = 0`、
`prghSpread = 0`、重建压程 9.56475 = 9.81·0.975（格心 y 跨度上的静水头），
maxCorrectors 1、worstContinuity 0——第一遍 U = 0 就满足动量。注意这条不覆盖
上面 Boussinesq 偏差驱动的盒子：偏差力就是物理本身，不能被吸收，压住它要
buoyantFoam 式的浮力修正压力/通量处理，仍待做；`energyBuoyantRest` 的断言
继续按"回答不了这个力的耦合会给出什么"来定标：
Stokes 速度 `f H²/ν = 2.45 m/s`（比残差高四个量级）与压降的**符号**。

**涡扩散用梯度输运类比 + Pr_t**（`alphaCell = kappa/(rho cp) + nut/Prt`）。壁面
处理照搬动量壁函数的结构：`faceGamma()` 在调和平均上**只覆盖** type 为
fixedValue 的 wall 面（heatFlux/绝热墙保留自己的通量定义），给出温度壁面律

```
T+ = Pr·y+                       (y+ ≤ y+_T，导热副层)
T+ = Prt/κ·ln(E·y+) + Prt·P      (Jayatilleke P = 9.24[(Pr/Prt)^0.75−1][1+0.28e^(−0.007Pr/Prt)])
α_face = u_tau·y / T+(y+)   ⇒   q_w = ρcp·u_tau·(T_p − T_w)/T+
```

y+_T 是两支的交点（从最小值点 `Prt/(κPr)` 向上括根二分）。副层极限 `T+ = Pr·y+`
恰好回到 `α_face = ν/Pr` = 分子扩散率——解析副层的网格自动退化为纯导热，所以
同一个代码路径对两类网格都成立。分子 Pr < 0.5（液态金属，两支律失效）直接
`cfFatal`。`patchHeatRate` 与装配矩阵共用同一个 `faceGamma()`，报告的墙速率和
方程实际用的通量恒等；`thermalWall` 逐 patch 报 12 位的 u_tau/y+/T+/q。
黏性耗散两项都可开：分子应力 `ν|S_ij|²/cp`（`strainRate` 由湍流层给）与
**湍流释热** `ε/cp`（`turbulentDissipation`）。ε 存在湍流层上
（`RansModel::eps_`，`turbDissipation()` 访问器；每次 `correctTurbulence` 与
`initializeUniform` 末尾经虚函数 `updateTurbDissipation()` 刷新）：k-ε 就是
输运变量本身，SST 按 `ε = β* k ω` 换算，层流模型报零；能量层在
`turbulentDissipation true` 而没有湍流模型、或场尺寸不符时直接 `cfFatal`。
`volumetricHeatRate()` 报**全部已开源项的体积分**（ρcp·ΣV·S），不再只是
volumetricSource 一项。守恒语义要读对：总平均内能加热 = 分子 ν|S̄|² + ε——
`nut|S|²` 把平均动能抽进 k，ε 再把它放出来，两项无重复计算；但**壁面函数的
壁单元不是能量守恒的**：那里的 (k, ε) 对被 `imposeWallEquilibrium` 代数钉在
对数层平衡上，不是输运出来的。湍流渠道实测 `source/pump`（源 vs 泵功
`force·Ub·V`）：SIMPLE 0.9548 / 1.1082、PISO（endTime 8）0.7623 / 1.1774
（k-ε / SST）；把 ε 源整个关掉读数只有 0.004 / 0.03。所以 checker
（`tests/check_turbulent_channel.py` 第 8 组断言）对泵功只断言可见性带
[0.5, 2]，tight 的断言放在离散能量平衡（SIMPLE 8.7e-16 / 2.6e-16 W、
PISO 2.9e-9 / 1.1e-7 W——瞬态残余就是储能率）、`∫ε dV` 对报告 source
（vtu 复核 1e-8 相对）与壁单元 `ε = u_τ³/(κ y₁) = 475.23`（两模型都对，
SST 经 `β* k ω` 读出）上。

**七个回归算例**（`tests/check_energy.py` + `tests/check_thermal_wall.py`，每个模式都拿内核**没被给的**数对表；
实测值即 ctest 打印的 notes）：

| 算例 | 核对什么 | 实测 |
|---|---|---|
| `energyConductionFixed` | 无流动 + 两端 fixedValue：逐点直线、壁面速率 `kA dT/dx` | 剖面误差 **0.000e+00** K；0.2 W 对 0.2 W；balance −3.975e-12 W |
| `energyConductionFlux` | `heatFlux` 的**符号**与 fixedGradient：Fourier 定斜率不定水平 | 剖面离 Fourier 解 1.600e-08 K；受热面 319.75 K，面上真值 320 K；balance −2.287e-10 |
| `energyConductionSource` | 体源：镜像对称 + 两端各带走 `S·Vol/2` | 对称 0.000e+00 K（25 K 升温上）；峰值 325 K 对 `SL²/8k` = 325 K；剖面离连续抛物线 1.563e-02 K（升温的 0.06%，即半格壁面偏移）；balance −1.112e-10 |
| `energyDecay` | 隐式 Euler + 绝热均值守恒 + **扩散算子自己的离散误差** | 200 步；离散特征值 9.864532 对连续 k² 9.869604；离散衰减因子 0.373802218 对连续 0.372707839（比真 PDE 慢 0.294%）；幅值历史离离散本征模 4.640e-13 相对；离散均值 0.000e+00；balance 0 |
| `energyBuoyantRest` | 浮力进入**压力能回答**的地方（含上面那条地板） | max\|U\| 2.706e-04 对无平衡时的 Stokes 尺度 2.453e+00；`dp/dz` 0.099972 对 `f_z` 0.098100，直线残差 3.475e-04 比 `f·H` = 0.04905；T 不漂移；balance 2.254e-11 |
| `energyChannelFlux` | 对流↔连续性通量的耦合 + `Nu_Dh = 140/17` + 剖面形状不变量 | 23 次扫描收敛；`Umax` 0.100000000833 对 `aH²/8ν`；流量 0.00134 m³/s（沿程散布 9.915e-09）；混杯升温 **0.0189393944 W** 对面供给 **0.0189393945 W**；`Tw−Tb` 2.42409318 K 对 `17qh/35k` = 2.42857143（`Nu_Dh` 8.2505 对 8.2353）；两壁不对称 0.000e+00 K；形状比 `(Tb−T_centre)/(Tw−Tb)` = 0.285116 对 `39/136` = 0.286765；balance −9.898e-15 |
| `turbulentChannel`（能量块） | 湍流释热 ε/cp：离散焓平衡闭合、source = ∫ε、壁单元 ε 对数层平衡、source/pump 可见性带（与摩阻断言同一个 checker，`tests/check_turbulent_channel.py` 第 8 组） | balance 8.7e-16 / 2.6e-16 W（SIMPLE）、2.9e-9 / 1.1e-7 W（PISO t=8）；∫ε(vtu) 对 source 1e-8 相对；壁单元 ε 475.23 = u_τ³/(κ y₁)（两模型，SST 经 β* k ω）；source/pump 0.76~1.18，关掉 ε 源读 0.004 / 0.03 |
| `thermalWallChannel` | 温度壁面律接入能量方程（湍流传热的墙面闭合） | 力驱动湍流渠道列（u_tau 0.6245 是输入，首格 y+ 52.04 = 4.2 个热副层），底/顶壁 310/300 K：两模型报告的 T+ 对 Python 重写律 **11.8956036 全同**（1e-8）；墙速率 = ρcp·u_tau·A·ΔT/T+（1e-8）；vtu 场对 checker 重述的离散问题（内部调和 γ + 墙面 u_tau·y/T+）最大偏差 **3.5e-5 K**（kEpsilon）/ **1.7e-5 K**（SST），重述系统通量恒定 1.5e-12 相对；壁函数通量是分子首格通量的 **3.08 倍**（闭合确实在起作用）；镜像对称 0.0；balance ~1e-11 相对 |

`energyChannelFlux` 的两条设计值得记：几何是 y ∈ [0, 0.2] 的平行板渠道，所以
`H = 0.2`、半通道高 `h = 0.1`、`Dh = 4h = 0.4`——`Umax = a H²/(8ν)` 用**全通道高**，
`Tw − Tb = 17 q h/(35 k)` 用**半通道高**，两者混用壁面偏移就差一倍。κ 从 0.01 降到
0.002 是为了把轴向佩莱数 `Ub L/α` 推到 33（Pr = ν/α = 5）——0.01 时等温入口沿轴向
倒灌 −0.0035 W，把 `Tb` 斜率压到 13.4（解析 15）并留下 0.65 K 的线性拟合残差；
形状比 `39/136` 与 q、κ、流量都无关，纯导热的场在那里给 0，
所以它断言的是对流-扩散平衡本身而不是某个系数。

**case 侧 keys**：`energy{ rho cp kappa turbulentPrandtl volumetricSource
viscousDissipation expansion referenceT gravity relaxation initialT initialMode
initialAmplitude initialAxis }`（除 `volumetricSource` 之类缺省外，`fromCase` 对
非法值直接 `cfFatal`：ρ/cp/κ 必须为正、`relaxation ∈ (0,1]`、`initialMode ≥ 0`）；
`patches{ name{ T | heatFlux | TType/TValue } }`；`initialMode > 0` 时初始场是
`initialT + A·cos(mode·π·ξ)`——它在零通量墙上正是中心差分 Neumann Laplacian 的
**精确本征向量**，这就是 `energyDecay` 能有闭式对照的原因。

### 5.12 几何界面重构 `multiphase/Plic.h` + `src/multiphase/plic.cpp`（M4：PLIC——法向 + 精确截体积 + 界面几何）
`include/cf/multiphase/Plic.h` + `src/multiphase/plic.cpp`，与 VoF 同处 `cf::vof`。它接住 §5.10 结尾那条
限制：**梯度路径**要 α 被抹开几个单元才量得出几何，而 PLIC 只需要 α 本身——每格切一个平面，
平面的位置由"它切下的体积 = 这格的 α"决定，所以一格厚的界面能被完整带走。

**几何契约**：`cellBoundary()` 把单元还原成带绕序的面环（owner 用 `faceSf`、neighbour 取负，环自身
面积向量与面法向点积为负时整环反转）；`clipConvexCell()` 对每个环做 Sutherland–Hodgman 裁剪，再按
绕 +n 的 atan2 递减把截面色排序（绕序相对保留侧的外法向 −n 为 CCW）；体积与一阶矩用带符号四面体
对一个固定原点求和（原点取顶点均值以 conditioning，矩存在平移坐标里，`centroid = origin + moment/V`）。
**只处理凸单元**（hex/tet/prism），非凸直接 `cfFatal`。

**反解偏移** `offsetForFraction()`：带括号的 Newton，`dV/d(offset) = −A_cut`，每轮同时收紧括号，
容差 `1e-15·V`，≤80 轮。

**带由 α 容差定，不用梯度阈值**：`bandTolerance < α < 1 − bandTolerance`（缺省 1e-6）就是界面带。
实测 40×40 上一格厚圆盘：α 带 60 格，PLIC 重构 60/60，而 `interfaceNormal` 的梯度阈值标 124 格、
其中 64 格根本没有界面——梯度阈值是**过包含**的，它按 |∇α| 的峰值比例选格，而一格带里 |∇α| 沿界面
方向本来就抖。法向四种取法：`youngs`（∇α 归一）、`centroid`（phase1 截体心 − phase0 截体心，在
亚单元几何上量方向而非在场上量斜率）、`mixed`（两者归一均值，缺省）——后两者都拿自己切出的那一刀
再精化一次；`fit` 不在自己这一格里读，读邻域各格的 α 把它切成一个平面（5.14）。

**精确初始化** `exactDiskAlpha()`：α = 圆盘与单元盒的解析交面积 / 盒面积，按 x 分断点（盒缘与
±√(R²−|y|²)）逐段积分 `∫√(R²−x²)dx`。它是"重构必须读回来的那个 α"，因此把几何误差和初始化误差
分开了。Σα 守恒到 2.1e-10 相对——极限由 `asin(±1)` 端点与 √ 处的放大决定，不是重构的锅。

**实测精度**（CTest `plic.unit`，全部对闭式而非存档值）：单位立方体轴向截 0.35 给出体积 0.35、
形心 (0.825, 0.5, 0.5)、截面 1、截面形心 (0.65, 0.5, 0.5)，whole/void 标志与 `offsetForFraction`
的逆都精确；对角截割对 `{x+y+z ≤ s}` 的闭式 `(1/6)[s³ − 3(s−1)³₊ + 3(s−2)₊³ − (s−3)₊³]` 与其
√3 倍导数，9 个截平面上 |dV| = 4.16e-16、|dA| = 9.67e-16；hex/tet（msh 4.1 与 2.2）/prism 网格上
`cellBoundary` 闭合到 `mesh.cellVolumes()` 与 `cellCentres()`，|dV|/V = 1.1e-16…1.7e-15、|dc| ≤ 2.0e-15，
斜剖两侧互补 ≤ 5.2e-16；每个平面把自己那格的 α 读回来 ≤ 3.4e-15（40×40）/ 9.3e-15（80×80）。

**两条老实地写下来的限制**：① 三种**单元内**取法的法向**方向**在一格带上不随分辨率收敛——单元中心上量
40→80 是 16.84°→16.49°（均值），截面色心上量是 8.34° 均值/14.11° 最大（case `vofPlic` 断言的口径），
40→80→160 也只从 8.34° 走到 7.59°：一格 α 里的方向信息本来就只有一格。这条由 5.14 的 `fit` 解掉。
② 带均值 κ = −3.696 对 −1/R = −5（−0.74），因为 κ 仍走 §5.10 的单元中心散度模板，一格带撑不起它。
所以**本增量是测量，不是动力学**：这一节的 α 仍由 §5.10 的迎风方程输运、CSF 仍读梯度路径的 κ；
动力学在 5.13（α 交给平面带走），法向换成邻域面积加权平面在 5.14，κ 换成截面几何在
5.15（`interfaceCurvaturePlic`，`surfaceTension{ kappaPlic true; }` 接入动量路径）。
`vofPlic` case（40×40、R/dx = 8、`fill exact` + `reconstruction{ normal mixed; }`、刚体旋转 0.1 rad）
实测：α 带 = 平面数 = 报告值 = 60，|n|−1 ≤ 6.4e-13，截面色心躺在自己平面上 ≤ 9.2e-13，一格带跨的
界面面积 0.111878 对圆周拉伸面积 2πR·dz = 0.125664（89%），相体积 π R² dz 守到 5.9e-15 且 40 步内漂移
≤ 6.7e-15；对流把锋摊厚后带 60→421 格、Σ截面 0.368——那是锋的厚度而不是圆的周长，checker 明说这点。

### 5.13 几何输运 `src/multiphase/plicAdvect.cpp`（M4：α 由平面带走，而不是由 α 的通量带走）
`vof::advectGeometric(alpha, phi, dt, reconstruction)`，声明在 `multiphase/Plic.h`。它换掉的正是 5.12
结尾那条限制的另一半：α 方程的迎风输运是一阶格式，一个格厚的台阶过它必被摊平；把 α 交给平面，
锋就只是平面在动，不会摊。

**可再裁剪的凸区域** 是这一步的类型支点：`ConvexRegion`（若干闭合外法向面环）+ `cellRegion(mesh, c)`
+ `clipRegion(region, n̂, offset)` 与 `clipRegion(region, halfSpaces)`。`clipConvexCell` 与它共用同一套
`clipLoops`（whole/void 早退、截帽从共面点回收、容差随**当前**区域尺度走，所以区域能被反复切）。
体积仍由带符号四面体对自身边界积分——`clipRegion` 保持闭合，因此"区域有多少落在盒里/扫过的棱柱里"
是一次积分而不是对形状的假设。单元测试里两刀切出的斜 slab 对上 `cubeCap` 之差（≤1e-12），
box∩box 的棱柱体积恰为 面积×距离。

**每面一个通量**：`sweptPrism(faceLoop, intoDonor, |φ|dt/|S|)` 是面的平移扫过的棱柱，写成半空间：
面平面与到达平面之间的 slab，加上每条边一个侧面（法向 `cross(n̂, edge)`，朝向棱柱内部一点取向，
而不是跟环的绕序——这让它对任意边数的多边形与任意存储方向都成立）。 donor 是 φ 的上游格，
其相区域再被这个棱柱切一刀就是这一面交出的相体积；full 格直接给 `|φ|dt`，空格给 0。
**同一个数在 donor 处减去、在 receiver 处加上**，所以在离散闭合的 φ 上守恒是构造性的，不靠求解器容差。

**两道限、迭代到不动点**：一个格的多个出流棱柱在角上重叠，因此可能一起要超过它持有的相；接收侧
同理受"空体积 + 自己正在交出的量"限制（否则相中间的满格会被判为无处接收而饿死整步）。比值沿面共享，
所以限流不破坏守恒。两侧假设互相作废（缩小一个入流同时缩小它的 donor 出流、也就缩小了那个入流的
额度），单趟过一遍会留下 α > 1 的格，因此两道限交替重算直到最大改动 ≤ 1e-13。

**边界是水库**：patch 的 α 取 `alpha.faceValue(f)`（其界面不重构），所以走几何输运的算例必须让界面
远离承载通量的 patch；`reconstruction{ transport alpha|geometric }`（缺省 `alpha`，其他值 `cfFatal`）
选路径，几何路径下每步先重构再输运，报告行 `geo t <t> flux N limited N rounds N minFactor C
maxCourant C clamp V volErr V min A max A`。

**实测精度**：`plic.unit` 上，C=0.4 的台阶锋走 3 步后与"平移后的锋"逐格差 3.33e-16、带仍只有 1 格、
进口交出的体积精确到 1e-12、clamp 恰为 0；40×40 圆盘绕中心旋转 40 步，逐步相体积漂移 ≤ 6.9e-16 相对、
clamp ≤ 4.42e-17、α ∈ [0,1]、带 60→70 格、限流最多 336 面且最小因子 0.5、最多 35 轮。
case `vofZalesak`（Zalesak 转盘：40×40、R/dx = 8、`fill exact`、刚体旋转一整圈 6.28 rad、dt 0.01、
Courant 0.195）实测：每面通量之争 ≤ 3.3e-14 相对、整圈漂移 9.4e-15、clamp ≤ 3.1e-17、α 精确落在 [0,1]、
限流最多 314/448 面（最小因子 0.5，正是同格两个出流棱柱在角上重合处）、最多 35 轮；平面读回自己那格的
α ≤ 2.55e-15（14 次报告），带 60→76 格、Σ截面 0.13372 对 2πR·dz = 0.125664（+6.4%）——迎风的 α 方程
在 40 步内就把同一道锋摊成 421 格、Σ截面 0.368。形状误差 checker 自己写闭式（圆盘与轴对齐单元的解析交面积，
与求解器的初始化独立，二者在初始场互验到 5.0e-13/格）：一整圈后与解析圆盘差其面积的 **3.85%**，等价于边界
法向位移 0.154 dx = 1.92% R，76 个界面格里 2 格的前后位置翻了一侧，质心偏 0.198 dx（2.5% R）。

**限制照旧要说清**：① 输运不改 κ 的精度——CSF 缺省仍走 §5.10 的单元中心散度（法向方向一项已由
5.14 的 `fit` 解决：同一个转盘 `vofFit` 的带更窄、周长更接近圆、形状误差更小）；κ 由截面直接给出在
§5.15 补上，`surfaceTension{ kappaPlic true; }` 可把动量路径的 κ 也换掉——但胶合压力能否吃下
截面力的分布，见 §5.15 的实测；② 3.85% 的形状误差是 PLIC 在粗网格刚体旋转上的典型残余变形，不是守恒破口；
③ 守恒前提是 φ 逐格闭合：给定速度场上到 round-off（max |div φ| = 2.2e-19），若 φ 由压力修正解出则只到
连续性残差；④ 只处理凸单元；⑤ 界面不能压在有通量的 patch 上（见上）。

### 5.14 邻域最小二乘平面法向 `NormalChoice::Fit`（M4：法向的方向随网格收敛）

§5.12 的限制 ① 说得很硬：`youngs`/`centroid`/`mixed` 都在**自己那一格**里读 α 的楼梯，40→80→160 加密
时均值只从 8.34° 走到 7.59°——方向不收敛，因为一格厚的带在单元中心上本来就只装得下一格的信息。
`fit` 换的是**读什么**：不读这格的楼梯，读这格周围每一格界面上各自的那一刀。

**契约**：带里每格 j 都有一个把自己 α 切准的偏移 o_j（`offsetForFraction` 的唯一解）。界面若真是一个
平面，这些平面本就重合；令它们重合的**倾斜**就是界面的方向。对候选法向 n̂ 把邻域每格重切，得
(o_j, A_j, p_j)——截面积与截面色心。绕自己截面上一点转平面，体积只按二阶变（面环对自己形心的面积矩
为零），所以一阶关系是 o′_j = o_j + δ·p_j。在 n̂ 的切平面基里最小化 Σ A_j (e_j + δ·q_j)²
（q_j = p_j − 面积加权形心）给出 2×2 正规方程；解 δ、转、重切，最多 4 轮，tilt ≤ 1e-16 即停，
|δ| ≤ 0.35 rad 是信任半径。α 是体积分数数据，不是位置数据：与"把重构出的界面**点**拿来拟平面"不同
（那些亚单元位置沿一格带跳约一整格，比它们本该描出的弧还大），这样拟合是收敛的。

**邻域按顶点取**（不是共面）：一格带沿界面方向只角角相连，共面邻域会让一格除了自己谁也不用同意。

**阻尼不是防奇异**：邻域的截面色心几乎共线时——锋斜穿网格，对角带正是这样——切平面里有一个方向几乎
没有数据，正规方程在那条方向上读到的只是另一个方向斜率的噪声，无阻尼的拟合顺着走过去，停在自己的
邻域都不认的地方：同一片圆盘上最坏格 23°、邻域散差 0.1 格（带均值的 10 倍），三次加密都不退。对角加
1e-2·tr(G) 的 Levenberg–Marquardt 阻尼后最坏格 3.1°；1e-3 与 1e-2 给到位数相同的答案——阻尼只改走法，
不改它停下的地方（停下时 b = 0，阻尼动不了那个点）。

**轮数是取舍**：平面上 4 轮把方向收敛干净（0°，平面常数差 5e-8 格，邻域散差 6e-9 格），2 轮还留 0.004°
与 6e-5 格；但弯面上多给的轮数开始用来解释楼梯——16 轮时 R/dx = 8 的圆盘均值从 1.16° 退回 2.90°，报出的
界面面积从圆周的 97% 退到 91%。

**实测**（`plic.unit` 组 7 与 7b，全对闭式）：圆盘 40×40→80×80→160×160，`fit` 对解析径向的均值
1.16°→0.488°→0.304°（每加密一次至少掉四分之一，两次掉到 26%），最坏格 3.10°→1.40°→2.44°；同一组网格上
`mixed` 是 8.34°→7.94°→7.59°、最坏 14.11°→14.78°→15.23°，即不收敛。邻域散差按格计 0.0369→0.0196→
0.0093，每加密一次约减半。界面面积 0.122191→0.123813→0.118976 对周长 0.125664（97%/99%/95%），每格 α
读回 ≤ 3.0e-14。一个真的平面（30° 倾斜、16×16）：`fit` 的法向对解析法向最坏 <1e-6°、平面常数差 3.2e-9
（= 5e-8 格）、邻域散差 6.3e-9 格，而同一场上 `mixed` 均值 8.17°。带均值 κ 仍是 −3.75/−3.62/−3.56 对
−1/R = −5：限制还是 §5.10 的单元中心散度模板，不是法向——κ 这一块由 §5.15 的
`interfaceCurvaturePlic` 接走。

**接入**：case 里 `reconstruction{ normal fit; }`（`youngs|centroid|mixed|fit`，其他值 `cfFatal`），
求解器的 `plic` 报告行尾因此多出 `scatter <均值> rounds <最多>`（只有 `fit` 会给）。case `vofFit` 与
`vofZalesak` 是同一个转盘、同一套几何输运，只把法向换成 `fit`：初始平面 1.16° 均值/3.10° 最坏（`mixed`
8.34°/14.11°），一整圈后带 60→65 格而不是 76、Σ截面 0.128212 而不是 0.133723（周长 0.125664）、形状误差
3.3292% 的圆盘面积而不是 3.8453%、边界位移 0.133 dx 而不是 0.154、质心偏 0.141 dx 而不是 0.198、
0 格翻到锋的另一侧而不是 2 格；守恒、读回、α ∈ [0,1] 那些数不变，因为截割与通量是同一套。代价是每格
四轮邻域重切：ctest 上 51 s 对 11 s。

### 5.15 曲率由截面几何直接给出 `interfaceCurvaturePlic`（M4：κ 不再走单元中心散度）

§5.10 的散度模板 κ_o = (1/V_o)Σ(n̄_f − n̂_o)·S_f 从单元中心法向里读几何，而单元中心法向在
一格带上只有楼梯可读（§5.12/5.14）：`fit` 把**方向**换成邻域拟合后，κ 还欠着。`vof::interfaceCurvaturePlic(mesh, r)`
把这一块也换成截面几何——直接对重构出的界面点拟合二次曲面：

**做法**：对带里每格取**两环角邻域**（共享顶点的格再取一圈；单层网格上一环只有 3–5 格，且全部
截面点共高程——切向基把与法向最不对齐的轴取成 e1，挤出方向的列全为 0）持有界面的格，收它们的
截面形心 `r.areaCentroid`、按截面面积 `r.area` 加权。在**该格自己的法向标架**里（z 轴 = n̂）拟
z = L(c₀x² + c₁xy + c₂y² + c₃x + c₄y + c₅)，坐标按邻域切向跨度 span 归一；加权正规方程先做
**列剪枝**（加权列范 ≤ 1e-6·colMax 的基方向丢掉——单层网格上挤出方向只有 5e-7 的舍入数据，
不剪会被均衡放大成 ±5e5 的毒 κ），列均衡 + 1e-8 岭 + Cholesky；κ = div n̂ 在原点取值
（z_xx(1+z_y²)+z_yy(1+z_x²)−2z_xy z_x z_y 除以 (1+z_x²+z_y²)^{3/2}，系数按 span 反缩放）。
**欠定如实报 0**：点数少于存活系数、无跨度、Cholesky 非正定，都报 0——"这里没量出曲率"，
而不是编一个。

**实测**（`plic.unit`，全对闭式）：精确充填圆盘 R/dx = 8/16/32，带均值 κ **−5.02…−5.12 对
−1/R = −5**（同一组网格上 §5.10 的模板是 −3.25…−3.90），最坏格 5.87（1.17 倍 1/R），无一爆表；
40×40 精确盘上 `mixed` 60/60 全报，`fit` 在 160×160 上 220/252——`fit` 的法向在单层挤出场上有一个
面外退化的倾斜自由度，倾斜出去的格的截面点是 5·dx 的高程离群点，把邻域推成欠定，如实报 0
（`mixed` 的梯度/质心差在 z 向无倾斜分量，全部网格 100% 覆盖——这也是动量路径在挤出场上取
`mixed` 的理由，已写进 Plic.h）。**平面锋**（30°、16×16）：κ 最坏 **2.77e-14**、24/24 全报——
二次曲面过一组共面点就是平的，舍入之外一无所有。

**喂给 CSF**：`surfaceTension{ kappaPlic true; }`（需 `reconstruction{}` 块，否则 `cfFatal`）时
动量的 κ 取 `kappaPlicNow`（每步 `advect` 后先 `rebuildPlic` 再取），`vof csf:` 日志行的 kappa
因此改报**力实际读到的那张 κ 的带均值**（两条路径同算术，非 PLIC 时与原 kappa0 一致）。

**实测的边界（`vofDropletPlic` case，`fill exact` 单向静态液滴）**：κ 带均值 **−6.2918 对
−1/R = −6.25（0.67%）**，力内力净和 1.2e-12、镜像对称 0.14%、伪电流 8.4% 且增量递减——曲率与力
都对；但一格带上的力是 **delta**，胶式压力（Rhie–Chow 面梯度配单元力）无法在 delta 上立起斜坡：
跃变被摊到带的邻居格上，体相均值读 dp **−30.6%**（α>0.999 的"体相"里还有 23% 的斜坡没走完）。
抹开的 tanh 上试过另一头：二次曲面读到的是重构点描出的**水平集环带**（带宽从 0.05 加到 0.35 的
bandTolerance 扫描，κ −3.7→−9.1，毒恰好在 ∇α 最重的中带格——那里两侧的邻域把环带读进来），
dp −28% 且伪电流 32%。所以**这个开关今天是锐利场工具**：`vofDropletPlic` 按实测把 dp（40%）与
体相平直（30%）作为"分布位置"契约断言，κ（8%）照旧按 −1/R 断；要闭式 3% 的 Laplace
跃变仍走 `vofDroplet` 的抹开场 + §5.10 模板。闭合这最后一环需要力带上重新抹一层权重、或压力
格式能看见界面——下一个增量的题目。

### 5.16 密度基可压缩 Euler `apps/compressibleApp.cpp`（M3：激波捕捉首增量）

前面的 3D 层全是不可压的：压力把速度投回无散，声速只是"根本不存在的数"。`chflow-3d-compressible`
是第一条走**密度基**的回路——没有压力投影、没有线性求解器，守恒量自己带压缩性，激波由面通量
自己涌现。

**状态与格式**：每单元存**广延**总量 `(ρV, ρu V, ρE V)`，理想气体 `e = p/((γ−1)ρ)`、
`p = (γ−1)ρe`、`c = √(γp/ρ)` 每次 `primitives()` 从总量重算；时间推进是一阶 Godunov
`U_i += dt·Σ_f Φ_f`（广延量，所以**不除 V**——写成 `dt/V·ΣΦ` 时第 2 步 ρ 就跳到 ±4.3e4，
这条单位恒等式现在由 checker 的逐格复现钉住）。面通量 `Φ = F(U_L,U_R)·|S_f|` 取 **HLLC
近似 Riemann 通量**（Toro）：`S_M = (p_R − p_L + ρ_L u_L(S_L−u_L) − ρ_R u_R(S_R−u_R)) /
(ρ_L(S_L−u_L) − ρ_R(S_R−u_R))`，星态 `ρ*_K = ρ_K(S_K−u_K)/(S_K−S_M)`、
`U*_K = ρ*_K(1, S_M n̂ + u_tK, (ρE)_K/ρ_K + (S_M−u_nK)(S_M + p_K/(ρ_K(S_K−u_nK))))`，
通量 `F*_K = F_K + S_K(U*_K − U_K)`；切向速度原样穿过接触面。几何全部来自 `polyMesh`
（`faceSf()` 的带符号面积、`cellVolumes()`），所以同一份代码对块网格和非结构网格都成立；
`y/z` 两向的墙通量在单截面列上逐分量相消，离散系统**就是** 1D，横向速度精确保持 0。
边界一律取不可透滑移壁：`Φ = (0, p·S_f, 0)`——质量和总能量一格都不出门，动量交给墙。
于是守恒语义要按闭盒来读：`mom(t) + impulse(t) = mom(0)`，`impulse(t) = ∫dt Σ_b (p S_f)_x`
由 app 自己累加上报；锋未触壁的两段里墙压就是初值，闭式 `dt(p_R−p_L)A` 每步，checker 两者
都断言（实测 1e-9 相对）。

**波速是这个格式的全部难点**（也是本增量唯一一次"格式看着对、物理全错"）：声学包络
`S_L = min(u_L−c_L, u_R−c_R)`、`S_R = max(u_L+c_L, u_R+c_R)` 在隔膜处给 `S_R = 1.183`，
而 Sod 的真激波跑 **1.752**（RH：`p* − p_R = ρ_R S_sh u*`，`p* = 0.30313`、`u* = 0.92745`）。
估计的扇形**没包住它要携带的那道波**，于是隔膜面的动量通量读 0.059 对真值 0.532、激波被摊成
一条谁也认不出的斜坡（单格最陡 Δρ 只有跃变的 18%）。Einfeldt/Roe 形式
`S_R = max(0, u_R+c_R, ũ+ã)` 也挡不住：Roe 平均 `H̃ = 3.3172 ⇒ c̃ = 1.1519`，`ũ+ã = 1.152`。
改用 Toro 的**激波自适应**估计：先取两稀疏近似解的压力
`p*_0 = [((γ−1)/2)(u_L−u_R) + c_L + c_R] / [c_L p_L^{−(γ−1)/2γ} + c_R p_R^{−(γ−1)/2γ}]`
的 `2γ/(γ−1)` 次幂，再 `q_K = 1`（`p*_0 ≤ p_K`，那侧是稀疏波，头速就是 `u_K − c_K`）否则
`q_K = √((γ+1)/(2γ)·p*_0/p_K + (γ−1)/(2γ))`，最后 `S_K = u_nK ± c_K q_K`。`q` 就是
Rankine–Hugoniot 激波速度关系本身，所以 `p*_0` 估准时 `S_R` 就是真激波速——Sod 上
`p*_0 = 0.30674`、`q_R = 1.66536`、`S_R = 1.7625` 对 1.7517。换掉估计后什么都不用再调：
星态读回 `p*` 到 0.02%、激波位置到 0.76 dx。

**checker（`tests/check_euler.py`，case `sodShockTube`）五组断言**，前四组各用一种不同的办法：
① **它就是那个离散解**——Python 逐格重放同一套 HLLC 递推（含同一个 `q` 估计、同一个散射与
更新顺序），400 格上写出的 ρ/u/p 与复现值差 **5.0e-13**（契约 1e-9）；下面任何物理断言只有在
这条成立时才有意义。② **闭盒守恒**：8 个报告点上质量与总能量 1e-12 相对、
`mom + impulse − mom(0)` 恰为 0、且锋触壁前 `impulse = dt(p_R−p_L)A·steps` 闭式（1e-9）。
③ **精确 Riemann 解**（`p*, u*` 由压力函数的 Newton + 二分、星态密度由 RH/等熵关系）——求解器
**没有**被告知其中任何一个。一阶格式的数值黏性 `|λ| dx/2` 把每道锋摊成几个
`σ_λ = √(|λ| dx t/2)` 宽（本网格：头 0.017、尾 0.004、接触 0.015、激波 0.021），窗口的豁免
因此按**物理宽度**而不是格数定，同一组容差在任何网格上说的是同一件事：远离锋的平台逐点对
精确采样比对（稀疏扇区不豁免，直接对自己那条自相似剖面比）——far left ρ 0.16%/p 0.22%、
fan 4.58%/6.35%、star left 0.76%/0.96%、star right 0.02%/0.02%、far right ~0；两道锋的位置读
密度剖面的**半跃变交叉点**（格心线性内插，不是"最陡那一格"——抹开的锋照样带着自己的位置），
接触 0.205σ（1.25 dx）、激波 0.091σ（0.76 dx）；`max u = 0.92961` 对 `u* = 0.92745`（0.23%）。
④ **收敛率**：全场 L1 误差在本格与 checker 自己写的加密格（x 加倍、dt 减半，CFL 不变）上
各算一次，`log2(coarse/fine)` 落在 [0.4, 1.5]——实测 ρ 0.66、u 0.83、p 0.77；带三道间断的解
里锋的扩散贡献 √dx，所以这不是"一阶"的漂亮数字，而是它真实的样子。⑤  positivity 与横向速度
恒 0。精确解自己先被钉在 Toro 公布的 Sod 星态上（4 个数 5e-4 内）才被用来当参照。

**本增量的边界**（往下补的顺序）：无黏——没有应力张量与热通量，所以还不是可压缩
Navier–Stokes；一阶——没有 MUSCL/TVD 重构，激波宽度就是 `σ`；单截面列——横向只用来验证
散度定理逐格成立，真正的 multidim 激波（斜激波、马赫盘）还没算例；也没和不可压层或 1D
两流体共享状态（`pvt` 的真实气体 `c(p,T)` 还没接进 `primitives()`）。

## 5.17 区域分解与 halo 交换 `parallel/DomainDecomposition.h` + `src/parallel/`（M3：并行首增量）

对标 OpenFOAM 的 `decomposePar` + `processorMesh` + `Pstream`，但**不绑定任何传输后端**：
本增量交付的是分解拓扑、耦合面消息计划与全局归约——MPI（或任何别的后端）只负责把
已经按 Isend/Recv 形状备好的缓冲区搬到线上。本机没有 mpirun/mpi.h，所以所有断言都在
**同进程内**按域循环完成：这正是"后端无关"能被验证的方式。

**分解**：递归坐标中位数**排序**（`orderBlock`）再连续切片 `d = (i·nDomains)/nCells`，
负载不均衡度天生 ≤ 1。切片是纯数组区间，因此 checker 能一模一样的复现它。排序有两个
必须钉死的细节，都是踩出来的：

- **质心有 ~1e-16 的舍入噪声**。一层格（沿轴向坐标相同）在精确算术下并列，实际却不：
  若直接按 double 排，切片的分界点由噪声决定，分解就不再"由几何可复现"。所以先按
  (坐标, 格号) 排，再把**噪声带**（1e-9·该轴向范围）内的连续一段整体按格号重排。
- **对称块的三个轴向范围在精确算术下相等**（单位立方体四面体格就是这种），于是"选最长
  轴"同样由噪声决定，第一次切的方向都能变（实测 hex 上 48 格、tet 上 76 格跑错域）。
  规则改成：在 1e-9·max 带内取最长轴，**并列时优先 x**。

**寻址**：`DomainFace{face, cell, peer, sign}`。同域内部面在本域列两条（±1），耦合面在
两侧各列一条，边界面只由 owner 列一条——三类相加恰好重组全局网格，app 每次运行都断言
这点（`boundary == mesh 边界数`、`internal/2 + coupled == mesh 内部面数`、`coupled 条目
== 2·coupled 面`）。halo 槽位**惰性分配**（`slotOf`：不在本域 map 里就追加到 owned 之后），
所以一个 peer 跨几个耦合面也只占一格，随之而来的消息不变式是
`exchange 搬动值数 == 2 · 耦合面数`（不是 halo 槽位数——这条一开始就写错了）。
每个 peer 一条连续、索引互相对齐的收发计划：发送只写自己的 owned 槽、接收只写 halo 槽，
因此没有次序冲突，天然就是 MPI 点对点该搬的东西。

**算子**：`DomainMesh` 带本域的逐面几何（`faceSf`/`faceCentres`/`faceDelta`/`faceOwnerWeight`）
与 owned 格的体积/质心。这里 `faceDelta` 必须是 `面心 − 全局 owner 质心`，不能是
`deltaCoeffVec()`（那是 C_n − C_o）——串行 `gaussGrad` 的非正交修正用的就是前者。
`assembledGrad`/`assembledDivPhi` 逐域扫**本域面表**（全局面号升序），对每个格做与串行
**完全相同**的浮点加序，所以 `gradDev`/`divDev` 是 0 而不是"某个小容差"。corrected
格式额外需要两次修正pass之间插一次 `exchange(grad)`——串行同意的做法。

**归约**：`reduceSum/Min/Max` 就是 Allreduce 的替身。app 打印各域 partial 的和/最小/最大
并与串行同序累加的结果对比（`sumDev/minDev/maxDev`），checker 再从 .vtu 独立重算
`Σ t·V`（体积由 checker 自己按散度定理在网格面上积出来）。

**测试**（`tests/check_parallel.py`，两个 case：`case.dict` 12×8×5 hex 切 1/2/3/5/7/8/16 域、
`caseTet.dict` 95 个 gmsh 四面体切 1/3/7/13 域）：

| 断言 | hex | tet |
|---|---|---|
| `gradDev`（分解算子 vs 串行算子，两格式两场） | 恰 0 | 恰 0 |
| `divDev`（分解净出流 vs `fvc::divPhi`） | 恰 0 | 恰 0 |
| `haloDrift`、`sumDev/minDev/maxDev` | 恰 0 | 恰 0 |
| 16 域与 1 域写出的场逐格相同 | 恰 0 | 恰 0 |
| 线性锚场梯度对解析值（checker 自己从 .vtu 算） | 4.8e-14 linear / 7.7e-13 corrected | 不断言，见下 |
| 无散通量 `U=(y,z,0)` 在内区格（**含切缝穿的格**）净出流 | 2.2e-19 | 4.2e-17 |
| checker 复现的 `domain` 场 + 逐域 owned/halo/internal/coupled/boundary/messages | 全对 | 全对 |

两个场各有分工：**线性**场是"与网格无关的精确性锚"（有界线性插值对线性面值是精确的，
所以散度和梯度该给出解析答案，切几刀都一样）；**二次**场不是——面插值对二次带 O(dx²)
误差，边界格上留下 `c₄·dx/4 = 0.0625` 的残差，它的作用是检验"分解逐位复现串行"。
（一开始把这两个角色混了：串行参照在**所有**面用解析值、域路径只在边界用解析值内部插值，
于是 1 域也读出 0.0625 的假偏差。）

tet 案例**不**断言精度，只断言"切割不改变精度"：在这 95 个粗四面体上串行 corrected 格式
连线性场的梯度都复现不了（实测 gradAnalytic ≈ 217，linear 也 5.7），因为斜面让有界权重
饱和截断——这是格式在该网格上的性质，1 域时就存在。checker 把它作为 note 如实打印。

**本增量的边界**：① 没有传输后端——收发计划、槽位与缓冲区都按 MPI 点对点形状备好了，
但 `MPI_Isend/Recv`/`Waitall` 与 `MPI_Allreduce` 还没接（需要装 MPI 工具链，属于系统级变更，
要先确认）。② 只有**显式**算子被并行化；`fvMatrix` 的装配与求解仍是整网格单进程，
域耦合的线性求解（半对角 + 界面 Schur/多重网格）是下一条。③ 没有重叠/多级 halo，
没有网格重排（scotch/metis），分解只有 RCD 这一种。

### 5.18 分区 1D ↔ 3D 耦合 `src/coupling/`（M4 首增量：稳态管网里嵌入一段三维详细区）

对标 Fluent 的 **system coupling** 与 OLGA 的 embedded detailed section：管网的一条边的
**内部**换成一块三维 SIMPLE 网格，两侧在界面上交换"该速率对应的压降"直到一致。
没有联合隐式矩阵、没有 Jacobian——是界面上的 Gauss–Seidel。新目录 `src/coupling`
（`cfcoupling`，`PUBLIC` 依赖 `cf1d` 与 `cf3d`）与 app `chflow-coupled-1d3d`。

**1D 侧的钩子**：`NetworkPipe::coupledResistance`（Pa·s/m³）。置为正数后，这条管的
全部摩阻/重位/局部损失被 `R·Q` **替换**（`pipeLossAt` 直接短路），1D 只带着一个接口
系数，细节物理全在三维侧。字典里写 `coupledResistance <value>;`，`loadNetworkFromDict`
校验三条硬约束：非负、该管不得是黑油/PVT 管、不得同时带阀或泵——三种混合语义下
"这条边的压降"不再是单一数。

**伙伴接口**：`ICoupledElement::lossAt(qMag)`（恒正，流向两侧同电阻）加
`lastSweeps()/lastConverged()`（默认 0/true，让解析伙伴也能被回归）。纯虚是为了
**接口本身可以先于三维实现被回归**：`tests/unit/couplingTests.cpp` 用一个
`dp = R0·q + k·q²` 的解析块测回路，三维块只是另一个实现。

**三维详细块** `DetailBlock`：网格由 `MeshLoader::fromCase` 读（`detail{ block{…} }`
或 `mesh{…}`），入口/出口 patch 由 `polyMesh` 自己的 patch 类型（`inlet`/`outlet`）识别，
流轴向取入口面积向量 `|ΣSf|` 的主导分量、`flowSign_` 修正朝向——**不**再维护第二份
接口描述，两套配置迟早对不上。构造时还断言两个 patch 的面积向量互相抵消到 1e-6
（否则"同一个截面"不成立）。

`lossAt(qMag)` 每次把入口速度重设为该速率隐含的塞流 `qMag/inletArea_`、出口 `p` 固定 0、
壁面 `U` 无滑移、其余面 `U` 与 `p` 一律 `zeroGradient`，再跑一次 SIMPLE，返回测得的压降。
**为什么不是 `symmetry`/`wall`**：SIMPLE 把三个动量分量装配在同一矩阵、共享对角
（`simpleSolver.cpp:110` 断言），而这两个 BC 只固定法向分量——第一次跑就被那条断言
撞崩。`zeroGradient` 对三分量一致，是这里唯一安全的侧条件。

`BlockMeasure` 带出最后一次解的全部可核验量：格数、SIMPLE 外迭代数、是否收敛、块测压降
`dpPa`（= ρ × 体积加权截面均值的首−末；三维核的 p 是运动压力 m²/s²，管网用 Pa，
换算在这里做）、被要求的速率、入口/出口 patch 通量、**每一横截面平面净通量的 min/max**
（含 patch 平面——这是离散的连续性陈述）、探针窗内的最小二乘梯度 `probeSlope` 与到直线的
最大偏离 `probeSpread`、入口面积、塞流速度、`uxMax` 与只限探针截面的 `uxCentre`。
数全从三维场自己量出来，回归因此分得出"块报了个不同的压降"与"块不再守恒"。

**外层回路** `solveCoupled1d3d`：

```
R ← couple.seeds 给定的初值
每趟：
  ① 以 R 作 coupledResistance 解网络  → q、dpNetwork = R·q、节点压力
  ② 块在 q 上解一次                  → dpBlock（块自己的测值）
  ③ residual = |dpNetwork − dpBlock| / dpBlock
     ——两个求解器之间的差，不是拿一个数和自己比
  ④ 记下这一趟**开始时**的 R 作为报告状态，再刷新 R ← (1−ω)R + ω·dpBlock/q
```

不动点就是"网络所收的压降 = 块在所交付速率上报的压降"，即界面上的连续+动量陈述本身。
层流伙伴 dp∝q ⇒ 一次刷新即命中；`relax` 与残差判据是为非线性（湍流）伙伴准备的。

**"报告的状态只有一致"是踩出来的**：曾有循环结束后再解一次取最终态的写法，于是
`out.resistance` 描述的系数与 `out.q`、块的 `fluxIn` 来自不同迭代——一个跨求解器的
不一致。现在的结构上不可能：接受的那一趟原样报告，循环后不再解。

**残差地板归谁管**：接口残差的下限是块自己解的误差（本例 `pcorrTol 1e-5` ⇒ ~1e-5），
`couple{ tolerance }` 设得比伙伴分辨率更紧只会白烧外迭代；网络侧 `solve{ tol }` 同理
受限于节点残差地板（对注入速率 ~1e-12 相对），本例因此取 `1e-10`。这两条都在 case 注释里。

**物理锚与实测**（`tests/cases/coupled1d3d`：`Src → A → { X(三维块) ∥ By(一维旁通) } → B → Snk`，
0.5×0.02×0.02 m 渠道、μ=0.1、ρ=800、注入 2.4e-4 m³/s ⇒ 塞流 ~0.39 m/s、Re_Dh≈64 层流，
入口长度 `≈0.06·Re_Dh·Dh≈0.077 m`、探针窗 `0.12–0.49` 落在其外）：

| 断言 | 数值 |
|---|---|
| 已发展段 `dp/dx` 对平面 Poiseuille `12μU/H²` | 1163 对 1183（−1.7%，带 3%） |
| 块实测 `dp/Q` 对层流 duct 电阻 `12μL/(H³W)` | 3.933e6 对 3.75e6（+4.9%，多的是入口段） |
| 端到端压差对闭式管网（解析电阻） | 718.2 对 699 Pa（+2.7%，带 4%） |
| 分支流量对闭式分流 | 1.5777e-4 对 1.6032e-4 m³/s（−1.6%） |
| 网络分支流量 == 块 `fluxIn` | 1e-12 相对（机器精度） |
| `fluxOut − fluxIn`、逐截面平面通量极差 | ≤1e-7 相对 |
| 剖面逐格对抛物线（从 .vtu 独立重算） | 最大偏差 0.0088 m/s（峰速 1.5%） |
| 横向速度 / 中心线速 | ≤1.9e-5 m/s / = 1.5 倍均值 |
| 压力链闭合、`dp_X == dp_By`、`pJ1−pJ2 == dp_X`、节点衡算 | 1e-9 相对 |
| 三个种子（1.0e6/3.75e6/1.2e7）的 q / dp / R 互差 | 1.6e-5 / 1.0e-5 / 4.7e-5 相对 |

**单元侧**（`coupling.unit`，解析伙伴）：线性伙伴恰好 2 趟命中且网络最终带的就是报告
电阻；`relax 0.5` 与 `1.0` 落同一不动点、只是趟数更多；二次伙伴收敛到 `dp(q)=Δp` 的
解析根（1e-9）；分支网络复现闭式分流与节点压力（1.8e5/1.6e5/1.4e5/1e5）；被 `maxOuter`
掐停的回路如实报 `!converged` 且残差恒等式仍成立。

**本增量的边界**：① 只打通**稳态、单相、不可压、单向**（网络→块给速率，块→网络给压降）；
瞬态耦合见 5.19，可压块、多相块未做。② 一次只详细化一条边，没有多块接口、没有块间的次结构。
③ 外层只有 Gauss–Seidel + 常数松弛——Aitken/IQN 拟 Newton 当时未做，强非线性伙伴靠 `relax`（§5.21 补上）。
④ 界面只在"压降↔速率"这一个标量上耦合；速度剖面、湍流量、相分数的界面映射未做——
伙伴今天可以是湍流块（§5.22），但过界面的仍然只有那一个标量。（速度剖面那一格由 §5.34 补上、
§5.35 带进瞬态、§5.36 换上湍流伙伴重新量过；相分数与湍流量 `k/omega` 过界面仍未做。）

### 5.19 分区 1D ↔ 3D 瞬态耦合 `solveTransientCoupled1d3d`（M4：把同一段三维详细区放进时间步进里）

同一张块网格、同一个网络，改成**从零速起步**：`time{ dt 0.02; endTime 1.2; }` 里每一步
界面上再套一层子迭代，交换的还是"压降↔速率"，app 是 `chflow-coupled-1d3d-transient`。

**伙伴为什么必须保持准稳态**（这一条决定了本增量的形状）：让块报**含惯性的**压降
（时程精确的 PISO 块，报 `ρL/A·(q−q_prev)/dt + 摩阻`）在分区 Gauss–Seidel 下不收敛——
网络解出来的系数会依赖上一步的速率，接口映射变成 `R → a − c/R`，其在 stiff 区的
线性化导数 `|c/R²|` 在本算例是 **2.7e8 ≫ 1**。所以三维侧继续只报"此刻该速率抵抗的压降"
（就是 5.18 的 `lossAt`），而区域的惯性 `ρL/A` **留在 1D 侧**，与网络其余部分在同一
隐式动量行里被解。这个切分是可反证的：时间常数 `τ = (M_X + M_By)/(R_X + R_By)`，把区域
惯量算两遍会把 τ 抬高 `M_X/(M_X+M_By)`——本算例是 **43%**。

**1D 侧的钩子**：`TransientPipeConfig::coupledResistance` + `TransientNetwork` 上四个新
入口。`step(dt)` 被拆成 `solveStep(dt)`（解出这一步、把状态挂在 `xPending_` 上，**什么都不
推进**）与 `commitStep()`（采纳挂起的步：流量、压力、持液输运和时钟一起走过去）。接口因此
能在**同一个时间层**上反复重解，而不是沿着轨迹漂：刷新 `R` 走 `setCoupledResistance`，它
显式作废挂起的步——那一步是用旧系数解的，已经不再是系统的状态（有单测钉住"重解逐位相同"
与"刷新后提交被拒绝"）。约束与稳态侧一致：`coupledResistance` 为正替换本管摩阻（重位和
惯量仍归 1D），且与阀互斥。

**每步的回路**：

```
每步 n：重复
  ① net.solveStep(dt)（不提交）→ q^n = |pendingPipeFlow(edge)|
  ② 块在 q^n 上解一次          → dpBlock
  ③ residual = |R·q^n − dpBlock| / dpBlock，≤ tolerance 就 commitStep()
     否则刷新 R 再来一遍（最多 couple{ maxSubIterations } 趟）
```

不接受"没谈拢的一步"：谈不拢直接 fatal 而不是提交一个两侧不一致的状态再往下滚。

**离散闭式（不是连续解）**：网络从静止起步，`t=0` 的连续方程并不成立，于是第一步特殊。
隐式欧拉加不可压端口给

```
q_X^1 = Q·(R_By + M_By/dt) / ((M_X + M_By)/dt + R_X + R_By)      // 突加速率按惯量分
q_X^n = q_∞ + (q_X^1 − q_∞) / (1 + dt/τ)^(n−1),  q_∞ = R_By·Q/(R_X+R_By)
```

若丢掉区域的惯量，序列会**直接**从 `q_∞` 起步——这就是上面那条反证在第一步的形式。

**实测**（`tests/cases/transientCoupled1d3d`，网格、流体、拓扑与 `coupled1d3d` 逐字相同，
所以末端就是稳态算例那个工况点）：

| 断言 | 数值 |
|---|---|
| 每步节点压降 == `R·q + M·Δq/dt`（X，用该步实测 R） | 最差 1.75e-10 相对（Newton 地板） |
| 同上 By / A / B（旁通用 Hagen-Poiseuille，lead 含首步惯性冲击） | 1.15e-13 / 2.0e-12 / 3.9e-13 |
| 每步接口残差恒等式、`dpNetwork == R·q`、`q == q_X`、`t == n·dt` | 1e-9 / 逐位 / 逐位 / 6.7e-16 绝对 |
| 每步三段衡算：`q_A == Q`、`q_A == q_X + q_By`、`q_B == q_X + q_By` | ≤1.4e-20 绝对 |
| 第一步速率对惯量分压闭式 | 1.38651e-4 对 1.38623e-4（2.07e-4 相对，带 1e-3） |
| 起步深度 `q^1/q_∞` | 0.879（带 <0.95：序列真的在爬） |
| 由序列回归 `q^n = A + B q^(n−1)` 的 τ 对闭式 | 0.26314 对 0.25885 s（1.66%，带 5%；算两遍惯量读 0.3698） |
| 同一回归的 `q_∞` 对闭式 | 1.57774e-4 对 1.57775e-4（1e-6 相对，带 1%） |
| 末端距工况点 | 0.14%（带 0.5%，且必须 >0：末端不是被假定的） |
| 收敛系数对层流 duct 电阻 | 3.93214e6 对 3.75e6（+4.86%，与稳态算例同一个数） |
| 每步块解趟数：首步 2 趟修正种子、其后 1~2 趟 | 60 步共 95 次块解（带 ≤120） |
| 块内守恒、已发展梯度、抛物剖面、横向速度 | 同 5.18 各带（`dp/dx` 1162 对 1182） |

**层流伙伴为何每步只需 ~1.6 次块解**：区域的 `dp/q` 与速率无关，接口把种子修正一次之后
系数就基本正确，其后每步的失配只剩块自己的解误差（`pcorrTol 1e-5` ⇒ 残差 ~1e-5…1e-4，
在 `tolerance 1e-4` 下多数步一趟即定、少数两趟）。

**顺带修的一处核问题**：SIMPLE 外收敛的 `settled` 判据原写作 `contErr <= bestContErr`
严格单调。瞬态接口每步都**热启动**，而热启动的首个 predictor 连续性可以低于该 case
最终落到的 HbyA/压力通量地板（实测 0.00768788 对 0.00769587），于是"best"再也打不破、
`|dp|` 分支永远走不到——首步的块解烧满 300 迭代并报未收敛。判据改成留 1% 余量
（`simpleSolver.cpp:343`），防发散的初衷不变；改后每步块解都收敛，本算例总时长 28 s。

**单元侧**（`coupling.unit` 新增 20 条，解析伙伴）：12 步线性伙伴逐步对上离散闭式（含
首步惯量分压与 `Δp_nodes = dpBlock + M_X·Δq/dt`）到 1e-9；每步 `block.lastQ() == res.q`
到 1e-15（一步之内两侧说的是同一个状态）；块解总次数恰好 `nSteps + 1`；`relax 0.5` 与
`1.0` 落同一条爬升曲线、只是趟数更多；`maxSubIterations 1` 时 fatal 且网络**仍停在起点**
（`time()==0`、流量为 0）。`transient.unit` 另钉住 `solveStep`/`commitStep` 的语义、重复
提交被拒、刷新系数作废挂起步、耦合管带阀被拒。

**本增量的边界**：① 伙伴仍是准稳态不可压层流块——区域的非物理"瞬时惯性"由 1D 侧承担，
块内的真实非定常（涡脱落、可压波）未进界面。稳态侧的湍流伙伴见 §5.22，瞬态侧的湍流伙伴见
§5.23——那条把"伙伴准稳态"这件事变成一条可反证的物理：分工的时间常数由伙伴的**斜率**定。
区域惯量的**来源**在 §5.24 变成一个可切换的归属：本算例的 `M = ρL/A` 用的是圆孔径，对方渠道
多报 27.3%，所以下面那张表的 τ 0.2631 s 是"孔径臂"的读数（区域自测给 0.2396 s，两臂末端速率差
0.107%）。
② 界面上仍是单个标量（压降↔速率）。
③ 步内子迭代当时是 Gauss–Seidel + 常数松弛，惯性效应带来的 added-mass 靠"把惯量留在 1D"
回避，而不是靠 Aitken/IQN 硬扛（界面的准 Newton 在 §5.21，瞬态臂已可用）。④ 块场在步内被解多次，只有**最后**一次留在 `measure()`
里，中间趟的场不可回放（记录的只有各趟的 sweeps/收敛位）。

### 5.20 多边同时详细化 `solveCoupled1d3d(net, edges, blocks, ...)`（M4：一张管网同时挂两块三维详细区）

5.18/5.19 只把一条边交给一个伙伴。真实站场常常同时有几处需要详细化的地方（一个节流器
加一段复杂几何），而核里"边"原本是标量参数。这一增量把它变成列表：

```cpp
CoupledResult        solveCoupled1d3d(OneDNetwork&, const std::vector<label>& edges,
                                      std::vector<ICoupledElement*> blocks, ...);
TransientCoupledResult solveTransientCoupled1d3d(TransientNetwork&, const std::vector<label>&,
                                                 std::vector<ICoupledElement*>, ...);
```

`blocks` 借用不拥有，长度必须等于 `edges`；`CoupledResult::edges/state` 与
`TransientCoupledStep::edges` 都按调用方给的**边序**排列，每趟的全部逐边记录留在
`OuterIteration::edges` / `TransientInterfacePass::edges` 里。原先的 `dpNetwork` 字段
去掉了：它就是 `resistance * q`，检查器按逐位恒等式自己重算，反而更硬。

**一趟里怎么做**（block Gauss–Seidel，稳态与瞬态同构）：

```
解一次网络（隐式，带所有边的当前系数）
对每条边 i：q_i = |网络给出的该边流量|，dpBlock_i = blocks[i]->lossAt(q_i)
                        残差_i = |R_i·q_i − dpBlock_i| / dpBlock_i
若 max_i 残差_i <= tolerance 就采纳这一趟；否则所有 R_i 一起按 dpBlock_i/q_i 松弛刷新
```

三条设计决定都值得说明：

1. **所有伙伴对着同一次网络解求值**。逐边各自迭代到收敛会让边 i 看到边 j 的半成品状态，
   两侧就不再是同一时间层上的同一个状态——5.18 的"一次交换说一件事"会被破坏。
2. **停止条件取最差边**（`res.residual == max_i 残差_i`），不是最后一条边。单元侧
   `coupling.unit` 与算例检查器都按这条断言（实测 `worstResidual` 与接受趟里较大者逐位相同）。
3. **成本是 1 + N**：一趟 = 一次网络解 + 每伙伴一次解，绝不更多。算例侧直接数日志里
   `SIMPLE iter 1:` 的条验证：3 组初值共 8 趟、块解恰 16 次 = 2 × 8。

**实测**（`tests/cases/coupled1d3dMulti`：`Src → A → { X1(3D) ∥ By(1D) } → X2(3D) → B → Snk`，
注入 3.0e-5 m³/s，Snk 定压 1.0e6 Pa，块 1000 + 4800 单元，全跑约 64 s）：

| 断言 | 数值 |
|---|---|
| 每边 `R·q ==` 网络自己的压降；每趟残差 == 两解之差 | 1e-9 相对带 / 1e-9 |
| 接受趟 `q ==` 网络流量 == 块入口积分流量 | 逐位 / 最差 5e-16 相对 |
| 串联边 X2 每趟 `q == Q`（拓扑钉死） | 1e-16 相对 |
| 两条 lead `q == Q`、`q_X1 + q_By == Q`、`dp_X1 == R_By·q_By` | ≤1e-8 / ≤1e-8 / 1e-6 带内，实测 ~1e-13 |
| 压力链 Src→J1→J2→J3→Snk 跨两个界面闭合 | 1e-9 相对 |
| 首趟真的远：X2 三组同为 0.66007，X1 最高 2.0649 | 界面不是被初值喂出来的 |
| 趟数 lo 3 / mix 2 / high 3，全部收敛 | 带 ≤8 |
| 三组初值逐边落同一点（q/dp/R 跨度） | 最差 2.4e-5 相对（带 1e-3） |
| 收敛系数对各自 duct 的平面 Poiseuille | X1 3.6856e6 对 3.75e6（0.983）、X2 1.55916e7 对 1.6e7（0.9745），带 [0.94, 1.06] |
| 端到端压降对串并联闭式 | 554.3 对 567.4 Pa（−2.3%，带 5%） |
| 分流对闭式 `q_X1 = Q·R_By/(R_X1+R_By)` | 2.01547e-5 对 2.00399e-5（+0.57%，带 3%） |
| 块内守恒（进 == 出）、逐平面同速率 | X1 ~1e-9、X2 5.4e-5 相对（热启动只解 9 趟的解精度，带 1e-3） |
| `probeSlope·L` 对该块报告的 `dpPa` | 0.16% / 0.68%（带 2%） |
| 中心线速度对平均（平面 Poiseuille 为 1.5） | 1.472 / 1.472（带 ±0.04） |
| Re（按各自 `q`） | 10.3 / 20.4，层流——上面那些线性闭式的前提 |

算例故意让两条边是**不同**的问题：并联边的速率由分压决定、串联边的速率被拓扑钉死，
于是有些恒等式可以完全不依赖块就验证。可伪证性靠两处：三组初值给 X2 的起点与流量都
相同 ⇒ 它的首趟 `dpBlock` 必须逐位相同（伙伴之间没有串信息），而 X1 的首趟速率在三组
之间差 2 倍以上（两边确实在各自移动）。

**单元侧**（`coupling.unit` 新增 35 条，两个解析伙伴）：稳态两边——并联边对上闭式分流
`Q·R_By/(R_By+R_1)`、串联边对上全流量、两条边各落在自己伙伴的系数上、一次共享刷新即定
（`iterations == 2`）、每个伙伴每趟恰好被调一次、`res.residual` 等于接受趟的较大者、
远初值（一边 40× 高、一边 8× 低）落同一点；一边线性一边二次（`dp = R0 q + k q²`）时，
分流对上二次方程的正根到 1e-6，且需要多于一次刷新。畸形输入被拒：同一边出现两次、
伙伴列表长度不符、空伙伴。瞬态两边——每步两个残差都收敛、串联边每步都等于注入速率、
并联边仍对上 5.19 那条离散闭式爬升、每边的动量行用**自己**伙伴的压降闭合
（`Δp_nodes == dpPa_i + M_i·Δq_i/dt`）、块解次数 == 趟数（每边一次）、fatal 的消息
**点名最差那条边**且时钟停在起点。

**顺手测到的一件事（原始结论里"发散"与"−730"两处已被 §5.21 的实测替换）**：非线性强的伙伴
会让未加速的 block Gauss–Seidel **走得极慢**。二次伙伴挂在并联边时，接口映射是
`R_2 → R_0 + k·Q·R_1/(R_1+R_2)`，不动点处导数 `−k·Q·R_1/(R_1+R_2)²`；用分流关系
`R_1(Q − q_2) = (R_0 + k q_2) q_2` 把 `k` 与 `R_2` 消掉，这个导数就是
`−(Q − q_2)/Q + R_0 q_2/(R_1 Q)`，**随 `k` 单调地从 0 走到 −1**：`k = 0`（线性伙伴）是 0，
一趟就中；`k` 越大越接近 −1，即锯齿形慢收敛——它不会发散，但也没有任何固定趟数预算兜得住。
实测：单边 `k = 1.75e10`（导数 −0.875）要 **183 趟**到 `1e-10`，40 趟收不下；多边（一线性一二次，
同 `k`）要 **211 趟**；瞬态同一问题每步只给 8 趟时未加速臂**谈不拢**（fatal），加速臂 20 步全过。
所以缺的不是"能不能收敛"，而是"多少趟"以及"固定预算下过不过"。

**本增量的边界**：① 多边当时只在稳态算例里跑过，瞬态多边只有单元覆盖（伙伴解太贵）——成本
压下来之后由 §5.25 合上：两条边同台的瞬态算例已发布，一行内核没改。
② 界面上仍是每边一个标量；把上面那笔"多少趟"的账还上的 Aitken/IQN 在 §5.21。
③ 伙伴仍是准稳态不可压层流块——多边算例到今天也只挂层流块；`details.X{}` 已经接受
`turbulence{}`，但两块同时湍流的回归成本没进预算（稳态单边的湍流伙伴见 §5.22）。
④ 一趟里各边的系数独立刷新，边与边之间不共享 Jacobian：
这是"一次网络解 + N 次独立块解"，不是块对角的隐式耦合。

### 5.21 界面准 Newton（Aitken / IQN-ILS）：把慢收敛的接口拉回几趟

§5.20 测到的那笔账——非线性伙伴让未加速的接口用 183 趟、211 趟才到机器精度——就是这一块。
界面本身就是"接口变量 ↔ 压力"的映射，正是 Aitken/IQN 的标准对象，而 ChuanFlow 此前在
`couple{}` 里只有常数松弛 `relax`。**它只是把同一批求解
排得更聪明**：一趟仍然是一次网络解 + 每伙伴一次解，`lossAt` 调用次数 == 趟数，一条断言
就能证伪"加速是靠多解几次换来的"。

**接口变量与残差**。每边归一化到 `x_i = R_i / R_i(种子)`，伙伴报的电阻与网络带的电阻之差
给残差 `r_i = (dp_i/q_i − R_i)/R_i(种子)`：归一化让 3.75e6 与 1.6e7 两条边在同一次最小二乘里
不被量纲大的那条主导（多边算例里两者差 4.3 倍，不归一化时 IQN 的边耦合方向会被 X2 吃掉）。

**一步怎么算**（`InterfaceQuasiNewton`，`coupledNetwork.cpp` 文件内）：

```
每趟：dR = x_now − x_prev, dQ = r_now − r_prev  →  放进滑动窗口（≤ qnHistory）
基：从**最新**的 pair 开始做 Gram–Schmidt，正交化后模长 < 1e-6·原始 的丢弃，
    列数上限 == 边数（接口算子的秩不可能更高）     →  {w_j, q_j}, q_j 两两正交归一
步：u = relax·(r − Σ_j q_j (q_j·r)) − Σ_j w_j (q_j·r)
    |u| > qnLimit·|r| 时整步缩到该上限（Venable 保护）
    x_next = x + u, 逐分量地板 1e-6·种子（电阻非正会给网络记一个反着自己速率的压降）
```

三个设计决定：
① **基从最新 pair 重建，而不是只追加**。单边接口所有 `dQ` 都平行（秩 1），"只追加 +
正交化才决定留不留"会把窗口冻在第一个 pair 上——第一版就这么错过：单边要 14 趟；改成
最新优先重建基后 11 趟，残差从 `8.5e-12` 到 `1.4e-13`。多边（秩 2）同理：40 趟 → 9 趟。
② **单边单 pair 恰好退化成 Aitken Δ²**：`r` 完全落在 `q` 上时垂直分量为零，剩下
`u = −(q·r) w = −Δx·r/Δr`，正是割线/Newton 步。这不是"近似 Aitken"，是它本身。
③ **瞬态每步重建历史**（`reset()`）：接口问题含本步的惯性项，上一步的 pair 量的是另一个映射。

**实测**（`coupling.unit`，解析伙伴 `dp = R_0 q + k q²`）：

| 场景 | 未加速 | 加速 | 界面上是否同一个解 |
|---|---|---|---|
| 单边并联，`k=1.75e10`，`tol 1e-10` | 183 趟（40 趟不够） | **11 趟** | 速率对闭式 `6.787e-3` 到 1e-6 ✓ |
| 多边：X1 线性 ∥ X2 二次，同 `k` | 211 趟 | **9 趟** | 两条边分别落在闭式分裂与 `q_1+q_2 = Q` ✓ |
| 线性伙伴（层流块就是这类） | 2 趟 | **2 趟** | 首趟系数逐位相同（还没有 pair 可用）✓ |
| 瞬态单边 20 步，`tol 1e-6`，每步预算 8 趟 | **fatal（谈不拢）** | 20 步全过，67 趟 | 末端速率与宽松预算的未加速臂同值 `1.25153e-3` ✓ |
| 同上，每步预算 120 趟 | 122 趟 | 67 趟 | 两条臂末步速率一致 ✓ |

最后一行是这件事最实用的形态：瞬态的每步预算本来就是有限的（步内"解而不采纳"不能无限滚），
未加速时首步即 fatal，加速后 6 趟以内每步都谈拢。

**接口映射导数的闭式**（顺手把 §5.20 的错数改掉）：`R_2 → R_0 + k·Q·R_1/(R_1+R_2)` 在不动点处
导数 `−k Q R_1/(R_1+R_2)²`，用分流关系消元后是 `−(Q−q_2)/Q + R_0 q_2/(R_1 Q)`：`k=0` 时 0，
`k→∞` 时趋于 −1，**中间不越过 −1**。所以未加速的 block Gauss–Seidel 对这一族伙伴不会发散，
但会任意慢——先前文档写的"−730、60 趟也收不下"是把导数公式里的 `R_1+R_2` 平方算漏了一项，
现已按实测改正（单元里那条断言现在测的是 `slope < −0.8` 与两条臂同解，而不是那个错数）。

**本增量的边界**：① 写这一节时，非线性伙伴只能由解析伙伴提供——3D 层流块报的 `dp/q` 本质是线性的，
所以 `verify.coupled1d3d*` 三个算例里加速臂与未加速臂同解（这是回归，不是收益证明）；收益证明
在 `coupling.unit` 的闭式上。下一个增量（§5.22）把 RANS 块接进接口，这条限制随之解除：真实三维
伙伴的 8 趟 → 5 趟就是同一笔收益在伙伴侧的实测。多边算例把 `accelerate yes` 打开后实测仍是 8 趟（`lo` 3 / `mix` 2 /
`high` 3）、块解 16 次，落到与未加速臂同一个工况点（X1 `q 2.01548e-5 / R 3.6856e6`，
X2 `q 3.0e-5 / R 1.5592e7`，压降 467.75 Pa）。② 界面仍是每边一个标量（电阻），不是 profile（场到场的 IQN）。
③ `qnHistory`/`qnLimit` 有默认值（8 / 100）且未在 case 里暴露；只 `accelerate` 一个开关可配。
④ 加速只作用于系数向量，块自己的 SIMPLE 收敛不受影响。

### 5.22 湍流伙伴块：接进接口的第一个真三维非线性伙伴（§5.21 的那笔账在 3D 侧兑现）

§5.21 的加速在真实算例上只能算回归：层流块报的 `dp/q` 几乎不随速率变——把同一套最小二乘指数
拟合跑在 `verify.coupled1d3d` 的 11 趟 `(q, dpBlock)` 上只得 `m = 1.06`（那 6% 是对流惯性叠入口段），
电阻只变 1.06 倍，所以接口靠种子修正一两次就到位；真正的非线性当时只存在于 `coupling.unit` 的
解析伙伴 `dp = r0·q + k·q²` 里。这一增量把
RANS 块（§5.9 的 `kOmegaSST`）接到 `couple{}` 的 `detail{}` 上，让**伙伴本身**成为非线性来源：
湍流渠道的 `dp ∝ q^1.74`，它的 `dp/q` 随速率变 2.1 倍——这正是 §5.20 那条接口映射的三维版本，
而且这张图里两头的电阻都在动：1D 侧的 A/B 落在 Re≈12700、旁通 By 落在 Re≈20400，
用的也是各自的湍流摩阻律。

**接进去的地方只有一处排序**。`DetailBlock::lossAt` 现在是
`applyBoundary(qMag)` → `RansModel::create(...)`（仅一次）→ `solveIncompressibleSimple(..., turbModel_.get())`：
模型的入口湍流由 `uReference()` 读 **inlet patch 上的 `fixedValue` 速度**算出
`k0 = 1.5(U·I)²`、`ω0 = √k0/(C_μ^0.25·L)`，所以先把界面这一趟送来的速率刷到边界上、再建模型，
就不需要任何场播种。模型对象在块的生命周期里只建一次、逐趟暖启动，
`nuEfficient = nu + nut` 进动量矩阵、`nut` 进 `k`/`ω` 的扩散，壁面按壁函数给 `nut` 与 `u_τ`。
层流侧没有多解一次、也没有多一条方程：`turbModel_` 为空时传给 SIMPLE 的指针是 `nullptr`，
矩阵装配与 §5.18 完全同一条代码路径（三个层流耦合算例原样通过；日志只是多了几个恒为 0 的键）。

**实测**（`tests/cases/coupledTurbulent`，同一个种子 `R = 3.0e6` 跑两臂）：

| 臂 | 趟数 | 接口残差 | 落到 |
|---|---|---|---|
| 未加速 block Gauss–Seidel | 8 | 4.30e-4 | `q 5.204663e-4`，`R 5.686581e5` |
| 准 Newton | **5** | 1.41e-5 | `q 5.203993e-4`，`R 5.688702e5` |

两臂同工况点：`R` 差 3.7e-4 相对、`q` 差 1.3e-4 相对——都在接口容差之内，所以省下的 3 趟是净收益，
不是换了个答案。块解次数 == 两臂趟数之和（13 == 13，数日志里 `SIMPLE iter 1:` 的条数），
加速没有靠多解几次换来的成分。

**块侧自己也要能对上账**，checker 里最硬的三条都是纯三维的：动量闭合
`|dp/dx|_发展 · H/2 == ρ·u_τ²` 实测 4.4651 vs 4.5134 Pa（1.1%）——注意左边必须用探针窗的已发展梯度，
用 `dpPa/L` 会把入口段算进去而差到 6%；涡粘对分子粘 152.8 倍；`y+` 96.7 在壁函数带（30–300）内。
剖面平坦度 `u_中线/U_块 = 1.094`（抛物线是 1.5），Pope 光壁律由块自己的 `u_τ` 反推来流速度得
1.3609 m/s、实测 1.3012（−4.4%，湍流模型与经验律之间该期望的符合度）。

**两条从实测里长出来的配置约束**（写进了 case.dict 头注释）：
① RANS 块的相对连续性残差停在 ~2e-3：修正后的压力残差能到 1e-10，但通量带着不平衡走——
所以 `couple.tolerance` 必须松于此（这里取 1e-3），`simple.continuityTol` 取 1e-2 而不是层流算例的 1e-6，
否则块每趟都烧满 100 次扫还不认自己收敛（第一版就是这样：一个臂跑 6 分钟没结束）。
② 种子电阻（3.0e6）是收敛值（5.687e5）的 5.3 倍，第一趟的速率只有正确值的 0.44 倍——
入口 `k0` 是按那个速率定的，之后暖启动；接口要几趟把速率找回来，这趟数就是 §5.21 加速的对象。

**本增量的边界**：① 本节只覆盖稳态单边。瞬态单边**已跑通**（§5.23：块 2560→320 单元、每步
热启动，一趟从 ~4 s 降到 ~0.1 s，27 步 × 两臂 88 趟共 9 s），所以"湍流块太贵、不进时间层"
只是当时没缩网格；多边 + 湍流仍无算例。② `turbulence{}` 的入口条件由第一趟的速率定，
不做"逐趟刷新入口 k/ω"（那是另一种迭代格式，会破坏 `lossAt` 只做一次求解的语义）。
③ 壁面仍靠壁函数，网格必须把第一格放在 log 层；没有低雷诺数近壁解析的选项。
④ 界面变量仍是每边一个电阻；`nut` 或剖面参与接口（场到场 IQN）在 §5.21 边界②里已经记下，仍未做。

### 5.23 瞬态湍流伙伴：非线性伙伴进时间层内，分工的时间常数由斜率定（§5.19 的接口 × §5.22 的伙伴）

§5.19 的瞬态臂带的伙伴是层流块，§5.22 的湍流伙伴只落在稳态单边。这一增量把同一条 RANS 渠道
放进时间层内的接口，代价与收益都不是"多跑一个臂"那么简单：**伙伴的 `dp/q` 一旦随速率变，
网络自己的分工动力学就跟着变**，而 §5.19 那条被验证到 1.6% 的公式正好是错的形状。

**为什么形状会错**。注入速率固定、总流量守恒，所以状态是一个标量：`q_By = Q_IN − q_X`，
两条支的动量行等于同一个节点压降，相减消掉 `dp12`：

```
(M_X + M_By)·dq/dt = f_By(Q_IN − q) − f_X(q),   f_i(q) = R_i(q)·q
```

在不动点附近线性化，`f_i′ = m_i·R_i`（`m_i = dln f_i/dln q` 是伙伴在该点的**割线斜率**，
`m ≡ 1` 就是纯电阻），于是连续衰减率 `S = m_X·R_X + m_By·R_By`；隐式欧拉逐步递推给出

```
B = 1/(1 + dt·S/M),   τ_eff = dt/ln(1 + dt·S/M),   M = M_X + M_By
```

层流时 `m ≡ 1`，`S` 退化成弦 `R_X + R_By`——所以 §5.19 拿弦验到 1.6% 并没有区分这两者。
这里两个伙伴都 `m ≈ 1.7`，弦形式把 τ 读高 62%。这条区分是可反证的，不是文字：

| 对照 | 值 | 备注 |
|---|---|---|
| 27 步爬升回归 | τ 0.2468 s，`B 0.7842` | 只取内部点最小二乘，一个时间常数跨 4.6 步 |
| 斜率闭式 `S = Σ m_i R_i` | τ 0.2460 s（**−0.34%**），`B 0.7835` | `m_X 1.684`（块自己跨 2.3% 速率变化回应的割线）、`m_By 1.729`（Re 11175） |
| 弦闭式 `S = Σ R_i` | τ 0.3990 s（**+62%**） | 层流算例验到 1.6% 的那条公式 |

`τ_chord/τ_slope = 1.622` 对阻力加权的 `m̄ = S/ΣR_i = 1.703` 差 −4.7%：差在 `τ_eff` 对 `S`
是对数关系而不是正比（`ln(1+x)/x` 那一步），checker 把这 6% 余量单独断言出来，没把它算成斜率
公式自己的精度。

**没有改一行源码**：瞬态 app 早就透传 `TurbulenceConfig` 与 `seeds`，`lossAt` 在时间层内被调
和被调在稳态循环里是同一处。控成本靠三件事——块从 2560 单元缩到 320（20×8×2）、步数缩到 27、
每步只给 6 趟；小网格叠每步热启动把一趟块解从 ~4 s 压到 ~0.1 s，两臂 88 次块解共 9.1 s。
`dt 0.06` 对 τ 0.247 s 是每常数 4.6 步、跑到 6.5 个常数之外，所以末端距渐近点还剩可测的尾巴。

**动量行的账，以及为什么旁路不能拿网络自己反证**。两侧同一状态仍由逐步行闭合：
详细支 `dp == R·q + M_X·Δq/dt` 到 **7.34e-11**，前后引管到 **2.47e-12**。旁路支就麻烦了——
它的"摩阻"就是 1D 内核自己算的那份，拿 `frictionDp` 去核对等于让网络跟自己比。checker 因此
在 Python 里重写一份 **Churchill 1977 公开式子**作外部锚，旁路行 `dp == R_churchill·q + M_By·Δq/dt`
闭合到 **7.33e-14**（机器精度），顺带把 1D 摩阻律在湍流分支上独立验了一次。

**接口与首步**。两臂首步都恰好 3 趟谈拢，全程最多 3 趟（预算 6），每步平均 1.63 趟；每趟块解
做 10–13 次 SIMPLE 扫（首步冷场 13 与 12，checker 只按 `>= 8` 断言它是冷场，因为热启动的那几步
也在同一量级）。逐步界面残差落在 4.4e-5…2.0e-4，最大在第 23 步，正压在 `tolerance 2e-4` 下面
——这个地板是 §5.22 那两条配置约束之一在时间层里的版本，稳态侧 ~2e-3 的通量不平衡一条没减。首步速率
5.09491e-5 对"按惯量分"闭式 5.09976e-5 差 9.5e-4（闭式里用第 1 趟的 `R`，不是收敛后的 `R`），
比它弛豫到的电阻工况点高 2.3%；末端 4.97828e-5 距拟合渐近线 0.0038%，正是 τ 预测的 0.0033%
指数尾。两臂（种子 1.0e7 在不动点外一倍、2.2e7 紧贴）落同一状态：速率 2.8e-5、幂指数 6.2e-5、
τ 3.0e-5 相对；块解次数 == 界面趟数之和（88 == 88）。

**块侧第一性原理照旧独立闭合**（这些不借 1D 的任何数）：`Re_bulk 19913`、`y+ 41.6` 在壁函数带、
动量 `|dp/dx|·H/2 == ρ·u_τ²` 到 0.35%、`nut/ν 86.6`、剖面平坦度 1.094、Pope 光壁律由块自己的
`u_τ` 反推来流速度 −5.2%、总压降是已发展梯度的 1.0279 倍（入口段占 2.8%）。

**本增量的边界**：① 伙伴仍是**准稳态**——`lossAt` 只报"此刻该速率抵抗的压降"，块内没有
`du/dt`，所以区域的真实非定常（涡脱落、启动压力波）不进界面；由此 τ 里只有电阻的斜率修正，
没有伙伴的相位滞后。区域的**整体动量**在 §5.24 已能由块自测并过界（进动量行，不进这条被
否决的系数形式），但那是惯量的归属、不是块内的时间导数。② `m_X` 是跨 2.3% 速率变化的**两点割线**，不是局部导数：弦/斜率差 62%
远大于这个不确定度，但 `m` 自身在两臂间有 1.684/1.687 的离散。③ 旁路也是湍流的（Re 1.1e4），
所以这 62% 不该记在"三维伙伴"头上——本增量钉的是"界面的 τ 由**伙伴们的斜率**定，无论伙伴是
三维块还是 1D 管"。④ 界面仍是每边一个标量；`accelerate yes` 在这里是默认配置而非收益证明
（两臂都 3 趟首步，未加速臂没跑）。⑤ 湍流 + 多边 + 瞬态三者同台的算例已由
**§5.27** 补上——两条边同时交给 RANS，一臂平谈、一臂准 Newton。

### 5.24 区域惯量的归属：块在自己的网格上量 `ρV/A²`，接口把它装进 1D 的动量行（§5.19 的接口 + 一条被实测否决的系数形式）

§5.19 把区域惯量留在 1D 侧，靠的是一个没说出口的假设：**代理管的名义孔径就是区域的截面**。
`pipes{}` 里给耦合管写的 `length/diameter` 因此同时决定了区域的电阻归属和它的 `M = ρL/A`。
不规则区域没有"那个孔径"可写（节流器、异形管、阀体），所以这一增量把假设换成**可对比的一个数**：
块在自己的网格上量 `ρ·Volume/Area²`，接口把那个数装进 1D 的动量行。

**形状误差是量出来的**（`tests/cases/transientCoupledInertia`：区域是 0.5×0.02×0.02 m 方渠道，
代理管取能塞进它两壁之间的圆管 `d = 0.02`，两臂只差这个数的来源）：

| 归属 | 值 kg/m4 | 出处 |
|---|---|---|
| 1D 名义孔径 `ρL/A` | 1273239.5447351625 | 圆截面 3.14159e-4 m2 |
| 块自测 `ρV/A²` | 1000000.0000000163 | 方截面 4.0e-4 m2，逐格体积和 |
| 比值 | 1.273240 | 对 `A_square/A_circle = 1/0.785398` 到 1e-12 |

同一种流体、同一段长度，圆孔径多报 **27.3%** 惯量——内切圆把四角的 21.5% 体积漏在门外。

**为什么惯量只能进动量行，不能进界面系数**。§5.19 从推导里否决过"块报含惯量的压降"
（映射 `R → a − c/R`，stiff 区导数 `|c/R²| = 2.7e8`）；这一条把那个否决**跑成了实测**：给伙伴
装上它自己的 `M`、仍按稳态形状把阻抗交给接口（`Z(q) = R(q)·q + M·(q − qPrev)/(dt·q)`），
分流的爬升能被接口分辨的条件正是 `dt < τ = M/R`——即"每个在迭代自己正在分辨的那条斜坡"的瞬态，
本算例 `dt 0.02` 对 `τ 0.24`。于是在任何 `q < qPrev` 的那一步 `c` 压过 `a`，系数过零。两次实测：

| 现象 | 数值 |
|---|---|
| 伙伴报非正阻抗，接口 fatal | `Z = −2180.209163 Pa`（该步提出 `1.00e-4`、上一步 `1.52e-4`，`M = 1.0e6`、`dt 0.02`）|
| 放宽 fatal 后块自己谈不拢 | 烧满 300 趟 SIMPLE，连续性钉在 **0.00522332**（corrected 4.74214e-10）|

（这两行来自被否决的那个原型，代码没留在仓库里——留下的是它决定的形状：`ICoupledElement` 只报
`lossAt` 与一个惯量数，没有"报含惯量压降"的那条口子。）

系数过零等于要求块去解一条"负阻 + 反向日惯量"的通道：块自己也不再认收敛，实测就停在通量不平衡上
（第二行的那组数）。**动量只进行，不进映射**——这不是调参能救的，是形式错了。

**现在的约定**：`couple{ inertias ( network partner ) }` 按**种子臂**给归属（长度 1 摊平到全部，
`network` 是缺省）。`partner` 时每条耦合边取它的块量的 `ρV/A²`，走
`TransientNetwork::setCoupledInertia`（src/1d/src/transient.cpp:399）替换该管的 `ρL/A`，动量行
仍是 `M/dt·(q − qPrev)`、隐式程度一字未改，只是 `M` 换了来源。测量是网格的属性、与时间层无关，
所以在任何一步被解之前取一次，且**全边校验完才装**——报不出区域的伙伴不能让网络半声明
（src/coupling/src/coupledNetwork.cpp:793）。面积按流向平面分组、每组 `Σ|Sf_axis|`；平面之间面积
差过 1e-9 相对就 fatal，因为截面在变的区域没有单一惯量可报（src/coupling/src/coupledNetwork.cpp:358）。
接口每条边回读装上的那个数（`EdgeExchange::inertance`，0 = 不在解任何时间层），块摘要行报它自测的
`inertance`，所以归属可核对而不只是可断言。

**可反证的断言**（同一张管网、同一个种子 3.75e6，只换一个数的来源；checker 七族，`tests/check_transient_coupled_inertia.py`）：

| 断言 | 实测 |
|---|---|
| 分区物理等价：两臂末端速率 | 0.000157557 / 0.000157726 m3/s，差 **0.107%**；外推工况点差 **0.0014%** |
| 惯量只改路径、不改去向：两臂速率差随步衰减 | 1.292e-5 → 1.693e-7 m3/s（峰值的 1.31%）|
| 首步分工纯由惯量定（各臂用自己的 M 对闭式）| 1.386514e-4（孔径）对 1.515730e-4（区域），1e-3 内；差 **8.2%** 的工况点；过报惯量的臂只得另一臂的 **0.9147** |
| `τ = (M_X + M_By)/(R_X + R_By)` | 拟合 0.263142 / 0.239584 s 对闭式 0.25885 / 0.235038 s（+1.66% / +1.93%）|
| τ 之比 = 惯量之比（被旁路 `M_By = 1.69765e6` 稀释后）| 实测 **1.09833** 对预测 **1.10129**（−0.27%）|
| 归属核对：日志 `M` == 该臂来源，逐位常数 | 孔径臂 == `ρL/A`、区域臂 == `ρV/A²`，各 1e-12；比值 == 面积反比 1e-12 |
| 每臂动量行用自己的 M 闭合 | 最差 1.75e-10 / 1.42e-10；旁路与前后引管 1.2e-13 / 2.0e-12 / 3.9e-13 |
| 界面在每个时间层谈拢 | **120** 层最差残差 **9.851e-5** < 该 case 的 1e-4 |
| 块侧不借 1D 的数 | 已发展梯度 −1161.58 / −1162.83 Pa/m 对面 Poiseuille 1182 / 1183；入口面 4.0e-4 m2 正是惯量所参考的截面 |

**单元侧**（`coupling.unit`，解析伙伴 + 一个覆写 `regionInertance()` 的局部块）：伙伴报的数
**等于**孔径给的数时两臂速率与压力历史逐位相同（1e-15）——这条把"分区物理等价"从算例观察升成
算术恒等；报的数不同则轨迹必须不同；每步 `edges[0].inertance` 就是装进去的那个数；护栏顺序是
**先验数、后碰行**——无伙伴的管、`0`、`NaN` 各自 fatal 且网络一个系数都没被改（`pipeInertia`
读回原值），只报不出区域的伙伴（`regionInertance() < 0`）也在装载前 fatal。

**代价**：两臂 120 个谈拢的时间层共 **2137** 次块扫、**175** 趟界面，checker 全程 53 s。两臂都不
需要准 Newton（§5.21）：系数仍是伙伴在它被给的速率上抵抗的电阻，惯量没进映射——收益正是上面
那条否决给出的分工。

**本增量的边界**：① 过界的仍只有**一个标量惯量**，它假设区域有一条等截面、且区域动量按
`V/A²` 折算成 1D 速率；速度剖面带来的动量修正系数 β 没有进界面（β 该由块自己的剖面算，这条渠道的
日志里 `uxCentre/plugSpeed = 1.474` 就是它的原料，但没被用）。变截面区域直接 fatal 而不是报一个错数。② 伙伴的
`lossAt` 仍准稳态，块内没有 `du/dt`：本增量搬过去的是**区域整体的动量**，不是块的时间导数项，
所以 §5.23 边界①那条"没有相位滞后"一字未改。③ 归属只有 network / partner 两种，`partner` +
湍流伙伴目前只在单元里覆盖；`partner` 同时挂两条边的**瞬态算例**已在 §5.26——两条边被名义孔径
错估的幅度不同（+27.3% 对 +69.8%），所以那不是把一个数复制两遍。**按边混合**归属（一条
`network`、一条 `partner`）当时还表达不出来——`inertiaFromPartner` 每次运行一个 bool
（src/coupling/include/cf/coupling/CoupledNetwork.h:314），**已由 §5.29 补上**：同一张动量系统里
X1 吃块测、X2 吃图纸，且逐边严格可分。④ 已发布的
`transientCoupled1d3d` 仍用缺省的 `network`，所以它记录的 τ 0.2631 s 带着这 27.3%——现在它是
**可量的**（区域自测给 0.2396 s），但没有顺手改那个算例：那会重写 §5.19 整张实测表。

### 5.25 瞬态多边同台 `chflow-coupled-1d3d-transient-multi`（§5.20 的多边 × §5.19 的时间层，一行内核没改）

§5.20 让两条边在**稳态**里同台，§5.19/§5.23/§5.24 的瞬态接口却始终只挂一条边。核里那个
`edges` 循环一直是真的（`solveTransientCoupled1d3d` 从 §5.20 起就收 `vector<label>`），
缺的是有人把两条边同时放进时间层里跑一遍、把账记下来、并让成本进得了预算。
**本增量的 `git diff` 里没有任何 `src/` 内核文件**：一个新 app、一个算例、一个 checker、两处
CMakeLists——它合上的是 §5.20 边界①那条"瞬态多边只有单元覆盖"。

**为什么挑这一对边**。拓扑照搬稳态多边算例 `Src → A → { X1(3D) ∥ By(1D) } → X2(3D) → B → Snk`，
`time{ dt 0.05 }` 注入 3.0e-5 m³/s、Snk 定压 1.0e6 Pa、四个节点都是不可压端口。两条边在同一趟里
被**同一次网络解**评估，却各自被要求服务完全不同的界面：

* **X2 在串联位置上**：不可压网络把 `q == Q_IN` 每步钉死（1e-16 相对），所以它的界面从头到尾
  只谈压降、没有任何速率可谈；但它自己的惯量仍要在动量行里付启动的账——首步节点压降
  1554.2153 Pa 对上 `R·q + M·q/dt`，其中 **1086.5 Pa 纯粹是那个突然的起步**。
* **X1 在并联支里**：速率是动量平衡的结果，首步按惯量分得 1.76228e-05，再按电阻分弛豫向
  2.01566e-05。它的界面每步都在谈速率。

一条"没什么可谈"、一条"每步都在谈"，而停止判据只有一个数：`worst = max_i 残差_i`
（src/coupling/src/coupledNetwork.cpp:849，测试 :855，全部谈拢才 `commitStep()` :903，
谈不拢 fatal 且**点名最差那条边** :889）。

**一趟的账**（隐式欧拉的时间层内 block Gauss–Seidel，:822）：

```
第 pass 趟：net.solveStep(dt)   ← 整张网络一次，带所有边的当前系数
对每条边 i：q_i = |网络此刻给出的该边流量|，dpBlock_i = blocks[i]->lossAt(q_i)
            残差_i = |R_i·q_i − dpBlock_i| / dpBlock_i，worst = max_i 残差_i
worst <= tolerance ⇒ 采纳这一步、推进时钟；否则所有 R_i 一起按 dpBlock_i/q_i 刷新，重来一趟
```

app 侧把**每趟**的逐边残差/系数/速率连同 `worst` 的**下标**一起打进日志
（src/coupling/apps/transientCoupledMultiApp.cpp:290），checker 从这些行重算出调度：
"停在最差边"因此是可从日志复算的事实，不是文档里的假设。

**实测**（`tests/cases/transientCoupledMulti1d3d`，块 1000 + 800 单元，两组初值 × 12 步，
checker 十族 76 条全过）：

| 断言 | 实测 |
|---|---|
| 每个时间层两条边同时 `residual <= 1e-4` | **24** 层谈拢；全程最差 8.5309e-5 / 8.5303e-5（两组都落在第 7 步）|
| 采纳的恰是每步最后那一趟：`lo` 21 趟（9 趟超差被重跑）、`mix` 20 趟（8 趟），每步 1 趟被采纳 | 预算 `maxSubIterations 8`，每步实际 1–3 趟 |
| 其中 **15** 趟只有一条边还在超差——先谈拢的那条跟着没谈拢的一起刷新 | 调度是全网一个数，不是每边自己的循环 |
| 块解次数 == 趟数 × 边数 | 42 == 21×2、40 == 20×2（合计 **82** 次块解；48 个被采纳的边态各记 9–22 次 SIMPLE 扫，`blockConverged` 全 1）|
| 每边动量行用**自己**装上的 `M` 闭合 `Δp_nodes == R·q + M·Δq/dt` | X1 3.04e-11 / X2 1.82e-13 相对；旁路与前后引管 1.3e-12…1.4e-11 |
| 串联边每步速率 == 注入速率，且首步压降含惯量项 | 0.000%；1554.2153 Pa 对闭式逐位，其中 1086.5 Pa 是启动 |
| 首步分工纯由惯量定（用该步的 `R`，不是收敛后的）| 1.76228e-05 对闭式 1.76228e-05 m³/s，**4.6e-12 相对** |
| 逐步衰减对上**隐式欧拉自己的** `1/(1 + dt/τ)` | 0.840671 对 0.841045（**−0.044%**），一个常数跨 6.3 步 |
| τ 的两种写法必须分清 | 连续闭式 `(M_X1+M_By)/(R_X1+R_By)` = 0.264555 s；但隐式欧拉在同一个 τ 上的**格式等效**常数是 `dt/ln(1 + dt/τ)` = 0.288834 s（**+9.18%**），整条轨迹拟合出来 0.288094 s 与它差 −0.26%——`dt/τ = 0.19` 时把这半步当成连续式会读错 9%，而它跟耦合没关系 |
| 弛豫到的去向 == 电阻分流 `Q·R_By/(R_X1+R_By)` | 预测 2.01566e-05 / 拟合 2.01554e-05（0.006%）；12 步 = 2.27 个 τ，末端距渐近点仍差 **1.87%**——故意的：checker 按离散闭式预测**整条轨迹**（每步 0.008% 内），不需要等它衰减完 |
| 两组初值只差 X1 的系数（1.25e6 对 3.75e6，3×），X2 同为 5.3e6 | 第 1 趟：`lo` 两边各差 **0.6605 / 0.6601**，`mix` 的 X1 已在 **1.877e-2** 而 X2 仍差 **6.601e-1**——一条边不能替界面报状态；`lo` 的 `worst` 下标是 0、`mix` 是 1，停在最差边因此是可复算的 |
| 逐边落同一状态（速率/系数/压降跨度） | 最差 1.02e-7 相对（带 2e-4）；X2 的速率两组逐位相同 |
| 邻居的初值不侵入自己的轨迹 | X2 第 1 步 `R` 两组 15590584.01 / 15590530.21，差 3.5e-6 相对——剩下的是它自己被热启动过的趟数，不是被喂的速率 |
| 收敛系数 vs 各自 duct 的平面 Poiseuille | X1 3.68466e6 / 3.75e6 = **0.9826**，X2 1.55906e7 / 1.6e7 = **0.9744**——与 §5.20 稳态多边记录的 0.983 / 0.9745 互锚 |
| 块侧不借 1D 的数 | 逐平面同速率、进 == 出；已发展梯度 145.6 对 148.3（X1）、1177 对 1200（X2）Pa/m；`uxCentre/plug = 1.472`；入口面 4.0e-4 / 3.0e-4 m² 正是自测 `ρV/A²`（1.0e6 / 1.06667e6 kg/m⁴）所参考的截面；Re 7.9 / 12.0 层流 |

**可反证性靠三处**（前两处借自 §5.20，第三处是这个算例独有的）：① X2 两组起在同一个系数上、
被同一个速率喂着 ⇒ 它的第 1 趟残差必须是同一个数（实测两组同为 **0.6601**），第 1 步谈拢后的 `R`
也只剩 3.5e-6 相对——伙伴之间没有串信息的通道；② X1 两组的系数与首趟残差都不同 ⇒ 两边确实在
各自移动，不是同一个问题跑两遍；③ **先谈拢的边不能自己收工**——那 15 趟单边超差把这条变成了观测：
若调度允许每边各自的循环，这些趟里邻居的块解根本不会发生，而实测块解次数严格等于
`趟数 × 2`。这条断言只有多边能立。

**成本怎么压下来的**（三处，都能核对）：
* X2 跨 0.03 m 的对称方向从 6 格降到 1 格：两侧都是 `symmetry` 的直渠道本质是 2D 问题，实测
  该块弦阻只因此移动 6.7e-5 相对，一趟便宜 6 倍。
* 步数 30 → 12：不再等爬升衰减完再验渐近点，而是用离散闭式预测整条轨迹（见上表第 10 行）。
* 初值 3 组 → 2 组：原先记录的"高种子那组每步烧满 300 趟 SIMPLE"用专门探针（同拓扑、同块、
  种子 1.125e7、跑 2 步）**复测不成立**（首步扫 11 / 10 趟），于是那一组删掉，文档只留能复现的数。

全算例实测 **235.3 s**：82 次块解、41 趟界面、24 个谈拢的时间层、两组初值。这些时间花在每个
时间层重谈一次界面（12 层 × 每层 1–3 趟 × 两条边）上，而不是花在边数上——稳态多边只谈 8 趟就
收工，这里 41 趟是因为 24 个时间层各要谈一次。

**那笔被删掉的账其实只是搬了家**：300 趟的 plateau 仍在，出现在 **`lo`** 组（日志第 496 行：
连续性钉在 **0.000683438**，corrected 6.54177e-11、dp 9.92704e-10、dU 1.49808e-07）——低系数那一侧
的块在冷场首趟会撞上上限，而接口照样谈拢，因为接受趟的残差地板由伙伴自己的分辨率给（本算例
最差 8.531e-5，正压在 `tolerance 1e-4` 下面）。这坐实了 §5.20/§5.24 反复出现的那条配置约束：
`tolerance` 问得比块答得更细，只会多烧扫，不会多要精度。

**本增量的边界**：① 多边 + 瞬态 + **湍流**三者同台仍无（§5.23 边界⑤一字未改）——这两块是层流，
所以 §5.21 的加速、§5.23 的斜率 τ 都还没在多边时间层里被量过；算例里 `accelerate no` 是**刻意**的，
测的就是两边同台的裸 GS。② 两条边的惯量都走缺省 `network`（代理管 `ρL/A`）——多边 ×
`inertias ( partner )` 的瞬态算例已在 §5.26，同样一行内核没改。③ 每边仍是**一个标量**系数、各边独立刷新，边与边之间不共享 Jacobian——这是
"一次网络解 + N 次独立块解"，不是块对角的隐式耦合（同 §5.20 边界④）。④ `tolerance` 全网一个数，
于是最钝的那个伙伴决定整张网络的趟数：那 15 趟单边超差就是它的代价形状，按边给容差没做。
⑤ 伙伴仍准稳态，块内没有 `du/dt`——多边不改变这一点，只是让它同时发生在两块网格上。

### 5.26 瞬态多边 × 区域惯量归属同台：`inertias ( partner )` 一次挂满两条边（§5.25 的多条边 + §5.24 的归属，一行内核没改）

§5.24 让**一条**边在两种惯量来源之间做选择，§5.25 让**两条**边进同一个时间层，但那个算例
两组初值都走缺省的 `network`。这两条边界在核里其实不存在：`solveTransientCoupled1d3d` 的归属
循环从一开始就是逐边的（先量完所有边、全过校验才装，src/coupling/src/coupledNetwork.cpp:793），
app 侧的 `couple.inertias` 也一开始就是**按种子臂**解析并透传的
（src/coupling/apps/transientCoupledMultiApp.cpp:176）。缺的是有人把两条边同时交给块自测的那个数、
把两边的账各自记下来，并证明它们**不共用一个误差**。本增量的 `git diff` 里同样没有任何 `src/`
文件：一个新算例、一个 checker、一处 `add_test` 注册与文档——用的还是 §5.25 那个 app 二进制。

**为什么这一对边值得单独记**。§5.24 的形状误差只有一个数（内切圆对 4.0e-4 m² 方截面多报 27.3%），
所以"归属影响多大"在它那里是一条斜坡。这里两条边的截面/名义孔径之比**不同**：X1 是 0.02×0.02 方渠道
配 `d 0.02`，X2 是 0.01×0.03 矩形配 `d 0.015`——

| 边 | 归属 `network` = ρL/A kg/m⁴ | 归属 `partner` = ρV/A² kg/m⁴ | 比值 | 比值 = 实截面 / 内切圆 |
|---|---|---|---|---|
| X1 | 1273239.5447351625 | 1000000.0000000163 | **1.273240** | 4.0e-4 / 3.14159e-4 |
| X2 | 1810829.5747344538 | 1066666.6666666756 | **1.697653** | 3.0e-4 / 1.76715e-4 |

两条边被同一个"名义孔径"错估的**幅度差 0.4244**，所以这不是 §5.24 跑两遍：一个全局的归属选择
现在同时改变两行动量方程，而它们各自的界面仍要分别谈拢。两臂只差这一个来源——种子逐位相同
（X1 3.75e6、X2 1.2e7），拓扑、网格、`dt`、容差、趟数上限全部照搬 §5.25。

**归属只动动量行，这一条在多边下重新被实测钉住**（三处，都是可从日志复算的）：

* 每条边日志里的 `M` == 该臂该边的来源，**1e-12 相对**；而块自测的 `inertance`（`X1`
  1000000.0000000163、`X2` 1066666.6666666756）**两臂逐位相同**且 == `ρV/A²`（1e-9）——
  被换掉的只有装进网络的那一行，块侧一个数没动。
* 界面回读的 `dpNetwork` 恒等于 `R·q`：全部边 × 全部步 × 两臂，最差残差 **0.0 Pa**。
  §5.24 那条"动量只进行、不进映射"的否决在多边下成立，不是单边的巧合。
* 串联边 X2 的速率两臂逐位相同（1e-12）——拓扑钉死 `q == Q_IN`，归属碰不到它；但它的**节点压降**
  在两臂差 **446.497744 Pa**，而闭式 `ΔM·Q/dt` 给 **446.497745 Pa**。首步那笔启动账因此是
  1086.4977 Pa（孔径）对 640.0000 Pa（区域），节点压降 1554.2137 对 1107.7159 Pa。

**并联边 X1 的分工与斜坡各自跟着自己的 `M`**：首步纯按惯量分流
1.7622833840174925e-05（孔径）对 1.9099580818343594e-05（区域），各对该臂自己的 `M` 的闭式
差 **4.6e-12 / 5.9e-13 相对**；`τ = (M_X1 + M_By)/(R_X1 + R_By)` 拟合 0.264555 对 0.240207 s，
比值 0.907968 对闭式 0.908028（−0.006%）；逐步衰减对上隐式欧拉自己的 `1/(1 + dt/τ)`，实测
0.840671 对 0.841045（−0.044%）与 0.827248 对 0.827709（−0.056%），整条轨迹按离散闭式预测到
0.008% / 0.005%。**去向不变**：两臂拟合的渐近速率 2.015535554e-05 对 2.015510611e-05，差
**0.0012%**；而 12 步（2.3 个 τ）截断处的末端速率仍差 **+1.23%**（1.9780002e-05 对 2.0024037e-05），
方向正是"多报惯量的那臂爬得慢"——checker 比的是拟合出来的目的地，不是被截断的尾巴。

**接口日程**（每步 1 趟被采纳，`maxSubIterations 8` 全程未触及）：

| 量 | `bore` 臂 | `region` 臂 |
|---|---|---|
| 界面趟数（12 步）| **20**（每步 [2,2,2,2,2,2,1,2,1,2,1,1]）| **17**（[2,2,1,2,1,1,2,1,1,1,1,2]）|
| 块解次数 == 趟数 × 2 | 40 | 34 |
| 谈拢的时间层 | 12，最差被采纳残差 8.5303e-5 | 12，9.5498e-5 |
| 首趟残差（X1 / X2）| 1.8770e-2 / 2.3030e-1 | 1.8016e-2 / 2.3030e-1 |

两臂的趟数差 3 不是块变慢——被采纳状态的 SIMPLE 扫两臂同落在 9–22 趟——而是**同一个接口在两条
不同的斜坡上重谈**：任何一步两边都必须同时到 `1e-4` 以下才采纳，`worst` 的下标记账停在哪个伙伴上
（首步两边都超差、下标 1，第 2 趟起一路是 0）。全算例合计 **74 次块解、37 趟界面、24 个谈拢的
时间层**。

**分区物理仍然等价**：两臂收敛后的系数只差 **7.40e-4**（X1）/ **1.27e-6**（X2）相对；末端系数
X1 3684657.459 对 3685402.007、X2 15590640.83 对 15590640.75，对平面 Poiseuille 是
**0.9826 / 0.9744**（孔径臂）与 **0.9828 / 0.9744**（区域臂），与 §5.25 那组的
0.9826 / 0.9744 对上——归属不移动伙伴对"电阻"的答案，它移动的是到达那里的路径。

**成本**：算例本身 364.7 s，checker 466.7 s（机器负载 10–13，与 §5.25 同一台机器同一时段）。
合计 **74 次块解 / 37 趟**，低于 §5.25 的 82 / 41——因为这里两臂都把 X1 种在它自己的系数上，
§5.25 那组把 X1 种低 3 倍的 `lo`（单独 21 趟）在本算例不存在，而 X2 的种子从 5.3e6 抬到 1.2e7
也省下几趟超差。归属自身的代价是两臂之间那 **3 趟**（20 对 17），没有别的。

**本增量的边界**：① 一个**按臂**的归属意味着两臂只能整体二选一：`inertiaFromPartner` 是每次运行
一个 bool（src/coupling/include/cf/coupling/CoupledNetwork.h:314），"X1 用孔径、X2 用区域"在**同一趟**
里表达不出来——两臂同台是本增量能给出的最接近的观测。该边界**已由 §5.29 补上**：`couple.inertiaFrom`
把同一个决定交到每条边手上，一条 `network`、一条 `partner`，一张动量系统里同时装上两个来源，
且逐边严格可分（被交出去的那条边与整体 `partner` 臂逐位相同，没交出去的那条与整体 `network` 臂相同）。② 归属仍是**每边一个标量**，等截面假设与
动量修正系数 β 一字未改（§5.24 边界①）。③ 湍流 + 多边 + 瞬态三者同台已由 **§5.27** 补上（同一拓扑、两块 k-ω SST 伙伴）；
那里为了把变量剩一个，用等效孔径把归属**中和**了，湍流 × 归属的叠加仍未测。
④ 伙伴仍准稳态，块内没有 `du/dt`。⑤ `tolerance` 仍全网一个数。

### 5.27 湍流 + 多边 + 瞬态三者同台：两条边都交给 RANS 伙伴，接口一趟平谈一趟准 Newton（§5.25 的两条边 × §5.23 的非线性伙伴，一行内核没改）

§5.23 把**一个** RANS 伙伴放进时间层内的接口（单边），§5.25 把**两条**边放进同一个时间层（两块
层流渠道），§5.26 让惯量归属同时挂满两条边（仍是层流）。§5.25 与 §5.26 各自记下过同一条边界——
"湍流 + 多边 + 瞬态三者同台仍无"。这条边界在核里也不存在：`turbulence{}` 从一开始就是**逐块**
声明的（`details.<edge>.turbulence`，两块各建一套 k-ω SST），`accelerates` 从一开始就是**逐臂**
解析并透传的（src/coupling/apps/transientCoupledMultiApp.cpp:161）。缺的是有人把两条边同时交给
非线性伙伴、把两个接口算子放在同一张管网上了同台，并且把"非线性改了哪些结论、没改哪些"逐条
记下来。本增量的 `git diff` 里没有任何 `src/` 文件：一个新算例、一个 checker、一处 `add_test`
注册与文档——用的还是 §5.25 那个 app 二进制。

**拓扑与两臂**

```
Src -> A -> { X1 (3D k-ω SST 0.2x0.005x0.005) || By (1D 0.3x0.008) } -> X2 (3D k-ω SST 0.12x0.004x0.008) -> B -> Snk
fluid { rho 1000; mu 1e-3; }   sources { feed rate 1.2e-4; }   time { dt 0.05; endTime 0.6; }
couple { edges ( X1 X2 ); tolerance 2e-4; maxSubIterations 6; accelerates ( 0 1 );
         starts { gs { X1 1.0e7; X2 1.5e7; }  iqn { X1 1.0e7; X2 1.5e7; } } }
```

两臂的种子逐位相同（X1 1.0e7、X2 1.5e7），只差 `accelerates` 这一个开关——所以臂间任何差别
都只能归给接口算子。水而不是 §5.25 那 800 Pa·s 的重油：1.2e-4 m³/s 进 5 mm 缝隙是
Re_bulk 19963、进 4×8 mm 是 30000，两条边都是壁面函数分辨得动的湍流；旁通 `By` 自己也在
Re 11156，一维侧同样不在层流区。

**惯量归属在这里被"恰好等效的孔径"中和掉了**。§5.24/§5.26 那条自由度在本算例不构成混淆：
两根代理管都取 `d = 2√(H·W/π)` 的 16 位数值，于是装进动量行的 `ρL/A` 与区域的 `ρL/(H·W)`
差 **2.2e-16**，块在自己网格上量的 `ρV/A²` 与同一个数差 **4.6e-15**（X1 两条都是
8000000 kg/m⁴，X2 都是 3750000 kg/m⁴）。`inertias` 保持缺省 `network`，两臂在这个算例里
**选不出差别**——剩下能测的只有伙伴的本构。

**非线性改写了什么**

* **斜坡的恢复项必须换成斜率和**（本增量最重要的一条）。X1 在并联位置上，速率从惯量分向电阻分
  弛豫。§5.25 那对层流伙伴用弦 `R_X1 + R_By` 就把时间常数算到 1.6% 以内；这里同一个式子给
  0.394059 s，**读高 64%**（实测拟合 0.240828 s）。换成伙伴的斜率和
  `m_X·R_X + m_By·R_By`（m_X = 1.682 由块自己两头的 `dp(q)` 答复拟合，m_By = 1.7285 由
  Churchill 弦在同一两点上算出）给 **0.241385 s，对拟合 +0.23%**。§5.23 在单边湍流上量到的
  那个"电阻加权指数"解释在这里一字不动地成立：`τ_chord/τ_slope = 1.632` 对
  `m̄ = 1.7012`（差 4.1%，其余是隐式步自己的 `ln(1+x)`）。串联伙伴的非线性**不进入**这个式子——
  它的速率不是未知量。
* **离散闭式仍然逐步管住整条轨迹**：每步衰减 B = 0.812518 对 `1/(1+dt·m̄R̄/M)` = 0.812907
  （−0.048%）；首步 X1 得 5.09933e-5 对闭式 5.10318e-5（−0.075%），离它弛豫去向的电阻分
  **2.6%**；拟合目的地 4.97833e-5 对两弦电阻分 4.96942e-5（+0.179%）；12 步（2.5 个 τ）截断处
  末端还差 0.2475%，正是它自己衰减尾巴该烧掉的 0.2477%；整条轨迹按斜率式离散闭式预测到
  **0.159%**（层流那组 0.16%——非线性没有把轨迹预测拖垮，只是把恢复项的算法换了）。
* **伙伴的律**：X1 的块在 2.2% 的速率行程上拟合出 `dp = c·q^1.682`（§5.25 的层流伙伴同法
  1.06）。X2 **读不出指数**：它的速率锁在注入值上（最差 1.1e-16 相对），它自己的 `dpBlock`
  整条斜坡只动 3.4e-8（iqn 臂 2.1e-7）——网络没有给它任何被测到的 `q` 变化。于是它的非线性只能
  对解析参照读：末端系数对同几何的平面 Poiseuille 是 **5.7664×**（X1）与 **7.8174×**（X2），
  §5.25 那两块同几何层流是 0.9826 / 0.9744。两个倍数不同（差 0.42×）正是"两条边不是同一道题
  做两遍"。

**非线性没改写的什么**

* 界面**仍是弦**：24 个边×步上 `dpNetwork == R·q` 的最差残差是 **0.0 Pa**，残差仍恰等于两个
  求解器答数的相对差（0.0）。湍流伙伴把 `R` 变成一个随工作点走的数，但网络充电用的仍是
  `R·Q`；动量行的闭合与接口映射仍然是两件事。
* 一维侧照旧用**公开式**独立判卷：`By` 的动量行对 Churchill 弦 + 自身惯量最差 **7.6e-14**，
  两条 lead 2.5e-12 / 7.2e-13；X1 的行 8.8e-11、X2 的行 2.2e-14（对块刚报的那条弦）。
* 串联边首步照样吃惯量踢：节点压降 11638.381828 Pa 对 `R·Q + M·Q/dt` 到 **1.9e-15**，其中
  **9000.0 Pa 是 `M·Q/dt`**，2638.4 Pa 才是摩阻。

**串联边在湍流伙伴下更便宜，不是更贵**（本增量第二条新事实）

| 量 | X1（并联，被谈的那条）| X2（串联，被钉的那条）|
|---|---|---|
| 速率 | 5.099e-5 → 4.991e-5，逐步下降 | 每步每趟恒等 `Q_IN`，1.1e-16 |
| 系数跨度（第 2 步起）| 随工作点走（2.2143e7 一路）| **8.8e-9** 相对（iqn 3.4e-8）|
| 被采纳残差（第 2 步起）| 最高 1.369e-4（贴着容差）| **7.74e-9**（容差的 1/25840）|
| 块 SIMPLE 扫 | 10–13 趟（iqn 10–14）| **每步恰好 6 趟** |

于是 `worst` 的下标在**两臂全部 53 趟**里都是 0：日程是 X1 一个人的。29 趟超差重跑里有 **25 趟**
X2 早已在容差内、仍被陪着刷一遍（多边日程是全网一个数的代价，§5.25 记过，非线性不改变它）。
层流那组里串联块要 9–22 趟 SIMPLE，这里它 6 趟就到位——**伙伴越非线性，被钉住的那条边越不需要
重谈，因为它连速率都没有被问；被谈的那条边则把整个预算吃掉**。

**两臂同台：一个目的地，一条稍短的路**

| 量 | `gs`（平谈）| `iqn`（Aitken/准 Newton）|
|---|---|---|
| 界面趟数（12 步）| **27**（首步 4 趟）| **26**（首步 3 趟）|
| 块解次数 == 趟数 × 2 | 54 | 52 |
| 超差重跑 / 被采纳 | 15 / 12 | 14 / 12 |
| 最差被采纳残差 | 1.3693e-4 | 1.3725e-4 |
| 每步最大趟数对预算 6 | 4 | 3 |

末端逐边落同一点：X1 速率差 **0.0002%**、系数 0.0001%、压降 0.0003%；X2 速率**逐位相同**、
系数差 5.4e-10、压降差 8.5e-11。两臂拟合出的斜坡时间常数 0.240828 对 0.241034 s（0.085%）、
伙伴指数 1.6820 对 1.6816。**加速在这个算例只买到 1 趟**——§5.24 那个锯齿问题在这里不存在
（接口本来就不慢），本增量的价值是三者同台而不是提速，这一点要按原样记下来。
一个诚实的副作用：首步 iqn 让串联边以 **2.001e-5** 被采纳，平谈是 2.313e-7（86.5 倍），
仍在容差的 1/10 内——准 Newton 作用在接口**整体**上，它把两条边的修正混在一起，于是在它把并联边
推得更狠的那一步，邻居带着更大的间隙被收下。

**块侧仍是自己的第一性**（每块、每臂都判）

| 块 | y+ | ν_t/ν | Re_bulk | `\|dp/dx\|·H/2` 对 `ρu_τ²` | Pope 给的 `ū` 对实测 | `u_max/ū` |
|---|---|---|---|---|---|---|
| X1（0.005 缝隙）| 41.7 | 86.8 | 19963 | **+0.37%** | 2.1045 对 1.9963（−5.1%）| 1.094 |
| X2（0.004 缝隙）| 58.4 | 96.8 | 30000 | **+6.5%** | 3.9546 对 3.7500（−5.2%）| 1.090 |

X2 那 6.5% 是网格：0.12 m 长只有 12 个 x 单元，插件入口剖面之上的壁面剪应力要走到
约 5 个水力直径（≈0.027 m）才饱和，探头窗里剩下的已发展段不足以把整块平均压降压到发展线上
——整块压降对发展线斜率是 0.9604（X1 同式为 1.0279）。两块都仍守恒到 1e-11、截面通量一致到
1e-11、界面末端压降与块自己场的压降逐位相同。

**成本**：106 次块解 / 53 趟界面 / 24 个谈拢的时间层（320 + 192 单元，两臂各自暖启动）。
checker 十族 141 条全过，含求解器在内一次运行约 90 s（机器负载 16）。

**本增量的边界**：① 伙伴仍**准稳态**——`lossAt` 只报"此刻该速率抵抗的压降"，块内没有 `du/dt`；
两条边都湍流不改变这一点，只是让它同时发生在两块网格上（§5.23 边界①、§5.25 边界④）。
② 界面仍是**每边一个标量电阻**，`nut`、剖面、界面法向应力到场到场的接口仍未做（§5.21 边界②）。
③ 归属仍是每边一个标量，且本算例被等效孔径**中和**成同一个数——湍流 × 归属同台没有测到，
不能因为两者各自做过就认为它们同台过。④ 多边湍流只测了"一并联 + 一串联"这一种组合；两条边
都在并联位置（各自有速率自由度、`worst` 会在两边之间换手）或都在串联位置（谁都谈不出速率）
未测。⑤ `tolerance` 与日程仍全网一个数。⑥ 日志头一行 `interface plain Gauss-Seidel` 取的是
`couple.accelerate` 的缺省值（src/coupling/apps/transientCoupledMultiApp.cpp:216），
`accelerates ( 0 1 )` 逐臂生效时它只描述第一臂；checker 按每臂自己的 `transient start` 行分段，
不读这一行。

### 5.28 湍流 × 区域惯量归属同台：两条边的孔径都写错、还错得不一样——归属只挪惯量账，不挪电阻（§5.27 的两块 RANS 伙伴 × §5.26 的 `inertias ( partner )`，一行内核没改）

§5.27 把三样东西摆上台，却**故意抽掉了一格自由度**：两根代理管都取区域自己的等效孔径
`d = 2√(H·W/π)`，于是 `network` 与 `partner` 装进动量行的是同一个数，归属在那个算例里**选不出
差别**（它的边界③原话：湍流 × 归属同台没有测到，不能因为两者各自做过就认为同台过）。§5.26 测了
归属，但两条边都是层流渠道。核里同样没有缺代码：`inertias` 逐臂解析并透传、`turbulence{}` 逐块
声明（src/coupling/apps/transientCoupledMultiApp.cpp:161）。本增量的 `git diff` 里没有任何 `src/`
文件：一个新算例、一个 checker、一处 `add_test` 注册与文档，用的还是同一个 app 二进制。

**把抽掉的那格放回去，并且放歪**——两条边的孔径都写成图纸上的名义圆，且歪的倍数不同：

```
Src -> A -> { X1 (3D k-ω SST 0.2x0.005x0.005) || By (1D 0.3x0.008) } -> X2 (3D k-ω SST 0.12x0.004x0.008) -> B -> Snk
fluid { rho 1000; mu 1e-3; }   sources { feed rate 1.2e-4; }   time { dt 0.05; endTime 0.6; }
couple { edges ( X1 X2 ); tolerance 2e-4; maxSubIterations 6; inertias ( network partner );
         starts { bore { X1 1.0e7; X2 1.5e7; }   region { X1 1.0e7; X2 1.5e7; } } }
X1 { diameter 0.005 }   X2 { diameter 0.0049 }      // 都不是 2√(HW/π)
```

| 边 | 截面 `H·W` | 名义孔径 | 区域自测 `ρL/(H·W)` | 孔径臂装的 `ρL/A` | 歪多少 |
|---|---|---|---|---|---|
| X1（1:1）| 2.5e-5 m² | 0.005 m | 8000000 kg/m⁴ | 10185916.3578813 | **1.273240×** |
| X2（2:1）| 3.2e-5 m² | 0.0049 m | 3750000 kg/m⁴ | 6363546.2460733 | **1.696946×** |

两臂种子逐位相同、都不加速，所以臂间任何差别只能归给惯量从哪来。两个歪倍数不相等这一条本身
是被断言的（`1.273240 != 1.696946`）：**不存在**一个全局修正因子能同时复现两臂的斜坡——§5.24
"2:1 矩形比 1:1 更经不起名义圆"这一点第一次落在非线性伙伴上。块自己那侧的自测值在两臂**逐位
相同**（8000000 / 3750000，就是 `ρL/(H·W)`），归属是臂的选择，不是块的选择。

**四条账，逐条分清**

* **界面仍是弦**：24 个边×步上 `dpNetwork == R·q` 的最差残差**两臂都是 0.0 Pa**，残差仍恰等于
  两个求解器答数的相对差（0.0）。孔径写错多少都不进这一行。
* **每边的动量行用本臂那份 `M`** 闭合：并联 4.0e-11 / 8.8e-11，串联 2.8e-14 / 3.1e-14；一维侧
  照旧由公开式判卷——旁通对 Churchill 弦 7.8e-14 / 7.2e-14，两条 lead 2.0e-16。
* **斜坡的时间常数是 `(M_X1 + M_By) / 斜率和`**——`M` 只在分子里，所以换归属**只动 τ、不动 `m`**：

| 臂 | 装的 `M_X1` | 拟合 τ | 斜率式预测 | 弦式预测 | 每步 `B` 预测/实测 | `c/s` 对 `m̄` |
|---|---|---|---|---|---|---|
| `bore`（孔径）| 1.01859e7 | **0.277231 s** | 0.276586 s（**−0.232%**）| 0.453205 s（+63.5%）| 0.834623 / 0.834974（+0.042%）| 1.6386 对 1.6988 |
| `region`（块测）| 8.00000e6 | **0.240829 s** | 0.241385 s（**+0.231%**）| 0.394059 s（+63.6%）| 0.812907 / 0.812519（−0.048%）| 1.6325 对 1.7012 |

  伙伴指数在两臂读出来是同一条律：X1 的块 **1.6762**（孔径臂，速率行程 8.0%）对 **1.6820**
  （区域臂，2.2%），旁通 1.7301（Re 11274）对 1.7285（Re 11156）——**惯量碰不到斜率**。
  整条轨迹仍各自被本臂的离散斜率式闭式管住：孔径臂到 **0.0086%**、区域臂到 **0.0023%**；
  首步也仍按**本臂的** `M` 分那股突然注入的流（4.52323e-5 对闭式 4.5021e-5，离它弛豫去向的
  电阻分 **9.9%**；区域臂 5.09933e-5 对 5.10603e-5，离 **2.6%**——同一道题，两臂的踢法不同）。
* **τ 之比就是惯量之比，被旁通稀释**：0.277231/0.240829 = **1.15115** 对
  `(M_X1+M_By)/(M'_X1+M_By)` = 1.15649（**−0.46%**）。而**目的地不含 `M`**：两臂拟合出的渐近速率
  4.9789e-5 / 4.97833e-5 m³/s 只差 **0.0113%**，末端速率却还差 **1.492%**（4.9162e-5 / 4.99065e-5）
  ——12 步只走完 2.16 / 2.49 个时间常数，`dt` 就是为了让这笔账在最后一步仍未结清。
  并且**把惯量写大那一臂落得更远**（−1.2592% 对 +0.2475%）：重的区域更慢到达电阻谈好的速率。

**串联边是归属的照妖镜**：它的速率不是未知量，于是两臂末端**逐位相同**（速率差 1.4e-20 m³/s、
系数差 8.8e-11），界面残差连首步那一下都两臂相同（**2.3003e-7**，容差的 1/869；第 2 步起
1.77e-9 / 7.21e-9）——归属**根本进不了界面这一行**。它只进启动的账：首步节点压降孔径臂
**17910.892818 Pa** = 摩阻 2638.3818 + `M·Q/dt` 15272.5110，区域臂 **11638.381827 Pa**
= 2638.3818 + 9000.0000；两臂之差 6272.510991 Pa 与 `ΔM·Q/dt` 对上到 1e-6 相对，闭式本身
到 3.8e-15。反过来说，被谈的那条边的系数在两臂间可以差到 **8.51%**——那不是归属动了电阻，
而是它先动了速率、电阻再随工作点走。

**日程与成本**（两臂都是平谈块 Gauss–Seidel）

| 量 | `bore` | `region` |
|---|---|---|
| 界面趟数（12 步）/ 超差重跑 | **30** / 18 | **27** / 15 |
| 块解次数 == 趟数 × 2 | 60 | 54 |
| 每步最大趟数对预算 6 | 4 | 4 |
| 最差被采纳残差（对容差 2e-4）| 1.617e-4 | 1.369e-4 |
| `worst` 下标 | **全部 30 趟 = 0** | **全部 27 趟 = 0** |
| 超差趟里串联边已在容差内 | 16 / 30 | 13 / 27 |
| 钉住的块 SIMPLE 扫 / 被谈的块 | 每步恰 6 / 10–13 | 每步恰 6 / 10–13 |

装更重的惯量要多花 3 趟——这是本算例里归属唯一买到的东西（接口更"僵"一点，每步多谈一轮）。
末端系数对同几何平面 Poiseuille：X1 **5.7064×**（孔径臂）/ **5.7664×**（区域臂），X2 **7.8174×**
（两臂逐位相同）——两条边仍是两道不同的题。

**块侧仍是自己的第一性**（每块、每臂都判）：X1 在 `Re_bulk` 19665（孔径臂）/ 19963（区域臂）、
`y+` 41.2 / 41.7、`ν_t/ν` 85.7 / 86.8、`|dp/dx|·H/2` 对 `ρu_τ²` +0.32% / +0.36%、Pope 光壁律给的
`ū` 比块实测高 5.5% / 5.4%（与 §5.27 记的 −5.1% 同一个间隙，只是反过来读）、`u_max/ū` 1.094；
X2 两臂同一状态：`Re` 30000、`y+` 58.4、`ν_t/ν` 96.8、
闭合 +6.48%（还是那 12 个 x 单元的入口段问题，整块压降对发展线斜率 0.9604）、`u_max/ū` 1.090。
每块都守恒到 1e-9 相对、各截面通量一致到 1e-9、界面末端压降与块自己场的压降逐位相同。

**成本**：114 次块解 / 57 趟界面 / 24 个谈拢的时间层（320 + 192 单元，两臂各自暖启动）。
checker 十族 135 条全过，串行跑完两臂约 45 s（ctest `-j 8` 下 142 s）。

**本增量的边界**：① 伙伴仍**准稳态**，块内没有 `du/dt`（§5.23 边界①）；归属说的是这块区域的
质量记在谁账上，不是块内多了一时间层。② 界面仍是**每边一个标量电阻**，剖面/`nut`/法向应力
到场到场未做（§5.21 边界②）。③ 归属仍是**每次运行一个全局选择**，只是这次两臂选了不同的值来
对照；按边混合（同一臂里 X1 用块测、X2 用孔径）**已由 §5.29 补上**——`inertiaFromPartnerPerEdge`
让一张表在同一趟里两种来源并存，并且逐边可分。
④ 湍流 × 归属只测了"一并联 + 一串联"，且两臂都不加速；湍流 × 归属 × 准 Newton 三者同台**已由
§5.30 补上**（完整 2×2 四臂 ⇒ 算子进不了归属的表，归属也不动不动点），三轴全交叉（四张按边表 ×
两个算子 = 八臂）在 §5.31；两边都在并联位置（`worst`
会换手）仍未测。⑤ `tolerance` 与日程仍全网一个数：30 趟里 16 趟
（区域臂 27 里 13）是早已谈拢的串联边陪着刷的。⑥ 本算例故意把孔径写错来制造归属的差别，
这不等于推荐这么建模——它测的是"错了以后，哪些结论仍然只对惯量敏感"。

### 5.29 按边混合惯量归属：一次运行里两个来源同时装上——每条边只回答自己那份（§5.28 的两块 RANS 伙伴 × 新的一格 `couple.inertiaFrom`，本增量动了 3 个 `src/` 文件）

§5.28 的边界③原话："归属仍是**每次运行一个全局选择**……按边混合（同一臂里 X1 用块测、X2
用孔径）仍未做——`inertiaFromPartner` 只有一个 bool"。一个 bool 不够，是因为区域描述的**质量
天生不齐**：一块是机加工渠道，块能在它自己的网格上量出截面；下一块是图纸上有人敲进去的名义
圆，谁也没法从它"量"出一个数。整网二选一要么把量得到的扔掉、要么把猜的当真，而这两种都不是
工程上想回答的问题。内核其实早就按边写：安装块逐边向 `blocks[i]` 要 `regionInertance()`、逐边
`setCoupledInertia`（src/coupling/src/coupledNetwork.cpp:796），缺的只是"每边可以不同意"这一格
控制。所以本增量**动了 3 个 `src/` 文件**，且只动这一格：`TransientCoupledControls` 加
`inertiaFromPartnerPerEdge`（src/coupling/include/cf/coupling/CoupledNetwork.h:330）、安装块把原来
的 bool 降级成**默认值**、表逐边覆盖它（仍是先全部检查再全部安装）、app 解析 `couple.inertiaFrom`
并把每臂每边的选择逐条打进日志（src/coupling/apps/transientCoupledMultiApp.cpp:201）。界面循环、
准 Newton、伙伴接口、块内求解一字未改。

逐边选择是**按臂**的一张表，所以"整体二选一"的两次运行仍是参照物：

```
Src -> A -> { X1 (3D k-ω SST) || By (1D) } -> X2 (3D k-ω SST) -> B -> Snk
couple { edges ( X1 X2 ); inertias ( network partner network );
         inertiaFrom { mixed { X1 partner; X2 network; } }
         starts { bore { X1 1e7; X2 1.5e7 }  region { 同 }  mixed { 同 } } }
```

| 装的 `M`（kg/m⁴）| X1 | X2 | X1 的来源 | X2 的来源 |
|---|---|---|---|---|
| `bore`（都吃图纸）| 10185916.36 | 6363546.246 | network | network |
| `region`（都吃网格）| 8000000 | 3750000 | partner | partner |
| `mixed`（一边一个）| **8000000** | **6363546.246** | **partner** | **network** |

三臂种子逐位相同、都不加速，块自测值三臂**逐位相同**（8000000 / 3750000 == `ρL/(H·W)`），
两条边被孔径错估的幅度仍不同（1.273240× / 1.696946×）。于是 mixed 那一行**本身就是一张表里
装着两个来源**，checker 十族 132 条断言全部围绕"这一格是不是严格按边兑现"来判：

* **叠加，不是平均**：mixed 的 X1 在**全部 12 步**上与 region 臂的 X1 **逐位相同**（速率、系数、
  块报压降、残差），而它的 X2 **不是** region 臂的 X2——被换掉的那条边跟着换，没被换的那条形同
  复制。反过来 mixed 的 X1 与 bore 臂不同：首步 5.09933e-05 对 4.52323e-05 m³/s（差注入速率的
  **4.80%**）。
* **界面仍是弦**：三臂各 24 个边×步上 `dpNetwork == R·q` 最差 **0.0 Pa**（不管这一边的 M 从哪来、
  孔径写错多少）。每边的动量行用**本边那份** `M` 闭合：并联 4.0e-11 / 8.8e-11 / 8.8e-11、串联
  2.8e-14 / 3.1e-14 / 3.1e-14（bore / region / mixed）；一维侧照旧对发表的 Churchill 弦
  7.8e-14 / 7.2e-14 / 7.2e-14，两条 lead 2.0e-16。
* **启动的账按边走**：串联边首步节点压降 bore **17910.892818 Pa**（= 摩阻 2638.3818 + `M·Q/dt`
  15272.5110）、region **11638.381827 Pa**（= 2638.3818 + 9000.0000）、mixed **17910.892818 Pa**
  ——与 bore 臂**同一个数**，因为它那一条边被交给了图纸；与 region 之差 6272.510991 Pa 恰等于
  `ΔM·Q/dt`。
* **斜坡也按边走**：拟合 τ bore 0.277231 s、region 0.240829 s、mixed **0.240829 s（与 region
  逐位相同）**——做斜坡的是 X1 那条边，它拿到哪个来源就按哪个弛豫，X2 的选择对斜坡**完全不可见**。
  预测仍是 `(M_X1 + M_By) / 斜率和`：0.276586（−0.232%）/ 0.241385（+0.231%）/ 0.241385（+0.231%）；
  τ 之比 1.15115 对惯量之比 1.15649（−0.46%）。伙伴的律不被归属碰：X1 的指数 1.6762 / 1.6820 /
  1.6820。**目的地不含 `M`**：三臂渐近速率 4.9789e-05 / 4.97833e-05 / 4.97833e-05 m³/s（差 0.0113%），
  末端速率却仍差 0.620%。
* **日程也按边走**：bore 30 趟（18 趟超差）、region 27（15）、mixed **27（15）与 region 完全同**——
  它谈的是同一个 X1，就花同样的趟数；`worst` 下标三臂全部 84 趟都是 0；块解 168 == 趟数 × 2；
  每步最多 4 趟（预算 6 从未尽用）；最差被采纳残差 1.617e-4 / 1.369e-4 / 1.369e-4（容差 2e-4）。
* **单元侧另有一条与 RANS 无关的判据**（`coupling_tests`）：同一个两条边的线性瞬态接口，
  `inertiaFromPartnerPerEdge = {1, 0}` 与 `{1,1}`、`{0,0}` 两次整体运行逐边对照——被改的那条边
  跟着 `{1,1}`（1e-12），没改的那条跟着 `{0,0}`，于是混合运行**既不是**任何一个整体选择；首步串联
  边节点压降 == `{0,0}` 那份，与 `{1,1}` 差 `ΔM·Q/dt` 到 1e-9。表少给一条边是 fatal；被问的那块
  不报区域也是 fatal，且**一条都不装**——每条边仍停在图纸的 `ρL/A`、时钟仍为 0。

**成本**：168 次块解 / 84 趟界面 / 36 个谈拢的时间层（三臂，两臂 24 个的时间照旧 ×1.5）。
checker 十族 132 条全过，串行约 147 s（ctest `-j 8` 下 194 s）；逐边选择对每条臂都打进日志
（`transient arm <臂> edge <边> inertia <来源>`），
用 `transient arm` 而不是 `transient start` 开头，因为后者是各 checker 用来切臂的表头。

**本增量的边界**：① 表是**每臂一张**，跨臂比较才有参照；单臂运行时它就是唯一的选择，没有第二
条腿可对照。② "两个来源"仍是**两个整体选择按边分配**；同一条边内部同时吃网格与图纸（比如把
区域切成两块伙伴）没有做，那需要块报一个分布而不是一个数。③ 伙伴仍**准稳态**（§5.23 边界①）。
④ 界面仍是**每边一个标量电阻**（§5.21 边界②）。⑤ 湍流 × 归属 × **准 Newton** 三者同台**已由 §5.30 补上**（完整 2×2 四臂，两条轴各自独立）。
⑥ app 头部的 `region inertia <词>` 行仍只描述第一臂，per-edge 时要靠逐边日志行读——头部那一行
如何同时如实描述三臂仍未解决（§5.26 边界里就记过）。

### 5.30 湍流 × 归属 × 准 Newton 三者同台：2×2 四臂，两条轴各自独立（§5.29 的拓扑 × §5.21 的接口算子，一行内核没改）

§5.28 边界④和 §5.29 边界⑤停在同一句话上——"湍流 × 归属 × 准 Newton 三者同台未测"。两条轴此前各
只被单独扳过：§5.27 拆算子但用等效孔径把归属**中和**，§5.28/§5.29 扳归属而全链都不加速。它们没有
理由互相泄漏：准 Newton 改的是界面的**更新**，而更新与"这块区域的质量记在谁账上"之间只隔着同一个
`EdgeExchange`。可"没有理由"不是证据，所以本增量把两条轴**交叉**成完整 2×2——
`couple{ inertias ( network network partner partner ); accelerates ( 0 1 0 1 );
starts{ bore_plain bore_newton region_plain region_newton } }`——拓扑、两块 k-ω SST 伙伴、两个写歪
的名义孔径、四组逐位相同的种子全部照旧，**一行内核没改**。十族 172 条断言从两个方向证明独立：

* **算子进不了归属的表**：同来源的两臂装上**逐位相同**的 `M`（bore 对 10185916.36 / 6363546.246，
  mesh 对 8000000 / 3750000 kg/m⁴），块自测的惯量四臂只有两个数（8000000 / 3750000 == `ρL/(H·W)`）。
* **界面在两种算子下仍是弦**：四臂各 24 个边×层上 `dpNetwork == R·q` 最差 **0.0 Pa**，残差仍是两个
  解之间的差（0.0）；每条边的动量行仍按**本边**那份 `M` 闭合——并联 4.0e-11 / 4.0e-11 / 8.8e-11 /
  8.8e-11、串联 2.8e-14 / 3.2e-14 / 3.1e-14 / 3.0e-14（bore_plain / bore_newton / region_plain /
  region_newton）——**准 Newton 采纳的状态仍是同一张动量系统的状态**。
* **算子轴买到的是日程**：同来源两臂末端速率 4.9162e-05 对 4.9162e-05（3.2e-07 相对）与 4.99065e-05
  对 4.99066e-05（2.0e-06），拟合渐近 4.9789e-05 对 4.97891e-05（2.1e-06）与 4.97833e-05 对
  4.97831e-05（4.0e-06），τ 0.277231 对 0.277255 s（8.9e-05）与 0.240829 对 0.241034 s（8.5e-04），
  串联边首步踢出 17910.892818 对 17910.848757 Pa（差 0.044 Pa）与 11638.381827 对 11638.329025 Pa
  （0.053 Pa）——**界面只有一个不动点**，两对里加速臂各少 1 趟。
* **归属轴照旧只管惯量账**：τ 是 `(M_X1+M_By)/斜率和`，四臂分别 −0.232% / −0.247% / +0.231% /
  +0.158%；τ 之比 1.15115 对惯量之比 1.15649 在**两对里都成立**；目的地四臂同一速率（渐近 spread
  0.0119%）；伙伴指数 1.6762 / 1.6763 / 1.6820 / 1.6816（极差 0.0058）——两条轴谁也不动谁的律。
* **准 Newton 那处已知不整齐被量化**：它作用在界面**整体**，于是那条速率被拓扑钉死、本来不需要谈的
  串联边陪着刷新残差：平谈 2.300e-07 对加速 1.670e-05 / 2.001e-05（**73× / 87×**），但都仍在
  `tolerance 2e-4` 之内，而真正被谈那条边不受影响（1.617e-04 对 1.617e-04、1.369e-04 对 1.372e-04）。
* **日程与代价**：30 / 29 / 27 / 26 趟（合计 112）、`worst` 下标四臂全部 112 趟都是 0、每层最多 3 或
  4 趟（预算 6 从未尽用）、块解 224 == 趟数 × 2；每块场独立闭合（X1 的 `y+` 随**本边来源**落在
  41.2（孔径臂）/ 41.7（区域臂）两点、`Re` 19665 / 19963、`|dp/dx|·H/2 == ρu_τ²` +0.32% / +0.36%，
  X2 恒为 58.4 / 30000 / +6.48%）。

**成本**：224 次块解 / 112 趟界面 / 48 个谈拢的时间层。app 四臂串行 93 s，checker 串行约 99 s。

**本增量的边界**：① 2×2 **不含按边混合**——那是 §5.29 的第三格自由度，三轴全交叉要 8 臂；本增量
证明的是"算子不改变归属的账"，于是按边表在准 Newton 下也应逐边可分——当时那是**推论**，
**已由 §5.31 做成实测**（四张表 × 两个算子 = 八臂，同来源臂在该边的历史逐位相同，分离度
1.6e8×..4.8e14×）。
② 伙伴仍**准稳态**（块内无 `du/dt`，§5.23 边界①）。③ 加速在这里只买到 1 趟：本算例的接口映射不
病态（平谈 2–4 趟就收敛），准 Newton 的价值在近 1 映射上（§5.21、§5.27）。④ 界面仍是**每边一个标量
电阻**（§5.21 边界②）。⑤ app 头部的 `region inertia` 行仍只描述第一臂（§5.29 边界⑥）。

### 5.31 三条轴同台：四张按边表 × 两个界面算子 = 八臂，按边归属在准 Newton 下是实测而不是推论（§5.30 的 2×2 × §5.29 的按边表，一行内核没改）

§5.30 把归属与算子交叉成 2×2，但它的边界①明说了那是**半张桌子**：四臂只用了两个"整界面"来源，
"这条边吃块测、那条边吃图纸"在准 Newton 下**没有被跑过**——它当时是两条已证结论的推论。本增量把
推论变成实测：`couple{ inertias ( network partner network network network partner network
network ); accelerates ( 0 0 0 0 1 1 1 1 ); inertiaFrom{ mixX1_gs{ X1 partner; X2 network; }
mixX2_gs{ X1 network; X2 partner; } mixX1_qn{ … } mixX2_qn{ … } } }`，四张表（两边都图纸、两边都
网格、X1 网格 × X2 图纸、镜像那张）各自在块 Gauss–Seidel 与准 Newton 下跑一遍 = **八臂**，拓扑、
两块 k-ω SST 伙伴、两个写歪的名义孔径、八组逐位相同的种子全部照旧，**`git diff` 里没有任何 `src/`**。
十族 **373** 条断言（`tests/check_transient_coupled_multi_three_axis.py`，RC=0）：

* **算子进不了归属的表**：八臂、两种算子，`network` 永远装 X1 **10185916.36** / X2 **6363546.246**
  kg/m⁴，`partner` 永远装 **8000000** / **3750000**（== `ρV/A²`），每步常数；孔径把 X1 多报 **27.3%**、
  X2 多报 **69.7%**，所以按边点名不是把一个数复制两遍。
* **界面在两种算子、四张表下仍是弦**：八臂 192 个边×层 `dpNetwork == R·q` 最差 **0.0 Pa**，残差仍
  是两个解之间的差；每条动量行仍按**本边**那份 `M` 闭合（并联 ≤8.8e-11、串联 ≤3.2e-14）。
* **可分性实测到几个量级**：判据按"这条边真正会动的量"来量——X1 量速率、被拓扑钉死速率的 X2 量节点
  压降。同一算子里，给某条边**同一个来源**的两臂在该边的历史上：准 Newton 四对全部**逐位相同**
  （0.0），平谈里 `partner` 源那对逐位相同、`network` 源那对差 **2.9e-10**；而给**另一个来源**的最近
  那臂差 **4.801e-02**（X1 速率）/ **3.502e-01**（X2 压降）。**分离度 1.6e8× 到 4.8e14×**。
* **混合表在加速下照样拼装**：`mixX1_qn` 的 X1 逐位等于 `region_qn` 的 X1、X2 的节点压降等于
  `bore_qn` 的 X2 到 **4.9e-07** 相对（平谈臂是 **1.1e-09**）；`mixX2_qn` 镜像成立（X1 逐位 == `bore_qn`，
  X2 差 **7.5e-07**）。也就是说准 Newton 让"另一条边的表"在**钉死那条边**上进到了 1e-6 量级（平谈是
  1e-9），但比来源本身的效应小五个数量级。
* **启动的账仍按边走**：平谈四臂的首步串联压降精确等于本边来源的踢——孔径 **17910.892818** Pa、区域
  **11638.381827** Pa，两混合臂之差 **6272.510991** == `ΔM·Q/dt`；加速臂付同一份踢到 **0.009 Pa**
  （17910.840016 / 11638.337767，两混合臂之差 6272.502249 对 `ΔM·Q/dt` 差 8.7e-3 Pa，1.4e-6 相对）。
* **斜坡仍只做斜坡那条边的决定**：`τ` 平谈 **0.277231**（X1 = 孔径）/ **0.240829**（X1 = 区域），加速
  **0.277255** / **0.241034** s；混合臂的 `τ` 与"做斜坡那条边同源"的整臂一致到 **0.0 / 4.0e-09** 相对，
  而来源那条轴本身值 **0.151**（15%）——另一边值 4e-9。**目的地不含惯量**：八臂拟合渐近速率
  4.97831e-05..4.97891e-05 m³/s，spread **0.0119%**。
* **算子轴仍只买到日程**：同一张表的两臂末端速率差 **9.7e-07**（`bore`/`mixX2` 群）与 **4.4e-06**
  （`region`/`mixX1` 群）相对注入速率；趟数 30 / 27 / 27 / 30（平谈）对 29 / 26 / 26 / 29（加速），合计
  **224** 趟界面、**448** 次块解 == 趟数 × 2、**96** 个谈拢的时间层，`worst` 下标 224 趟全是 0。
  准 Newton 那处不整齐按表复现：钉死那条边的残差地板 2.300e-07 → **1.670e-05 / 2.001e-05**
  （**73× / 87× / 87× / 72×**），真正被谈那条边 1.617e-04 / 1.369e-04 分毫不动，全部仍在 2e-4 内。
* **块侧仍自己闭合**，且 X1 的终点**按本边来源**分成两点：孔径臂 4.9162e-05 m³/s（`Re` 19665、
  `y+` 41.2、壁面剪闭合 +0.32%）、区域臂 4.99065e-05 m³/s（`Re` 19963、`y+` 41.7、+0.36%）；X2 恒
  58.4 / 30000 / +6.48%。§5.30 里那两个 `y+` 数值到此有了出处——**它们不是同一臂的数**。

**成本**：app 八臂串行 **275 s**（ctest 单条 262.7 s），checker 串行 **267 s**（含跑 app）。

**本增量的边界**：① 三条轴的交叉只做在**一并联一串联**这个拓扑上；两边都在并联位置（`worst` 会换手）
或都在串联位置仍未测。② 平谈下 `network` 源的 X1 被另一条边的表碰到 **2.9e-10**、`partner` 源的 X1
一点没碰到，而准 Newton 下两者都是 0.0——这个不对称是**第一次量到**，机理（网络 Newton 收敛地板与
逐边更新顺序的交互）没有在本增量里定死，所以只记成实测而不写成结论。③ 伙伴仍准稳态（块内无
`du/dt`）。④ 界面仍是每边一个标量电阻。⑤ app 头部的 `region inertia` / `interface` 两行仍只描述
第一臂，混合臂要看它自己的 `transient arm … edge … inertia …` 行。

### 5.32 两条耦合边都在并联位置：`worst` 会换手、一个节点压降要同时对两份动量行、四张表把斜坡排成严格有序的阶梯（§5.31 的八臂 × 新拓扑，一行内核没改）

§5.30 与 §5.31 的边界①写的是同一句话：**两边都在并联位置没测**。那里 X2 串在旁管之后，它的速率被
拓扑钉死，于是三件事在结构上就测不到——(a) `worst` 下标永远只会指向能谈的那条边；(b) 被钉住
那条边的残差地板其实是另一条边串话的读数，所以它会在两种算子之间单向漂 **73× ~ 87×**；(c) 斜坡 `τ`
只由"做斜坡那条边"决定，混合表在 `τ` 上和同源整臂重合到 4e-9，按边归属**动不了时间常数**。本增量把
拓扑换过来：`Src → A → { X1 ‖ X2 } → B → Snk`，两条耦合边接在同一对节点 `J1/J2` 上、中间没有任何旁
管，于是接口里**有两个都要谈的速率**。改的只有拓扑与两个名义孔径（X1 0.2×0.005×0.005 m、名义圆
5 mm；X2 0.12×0.004×0.008 m、名义圆 4.9 mm，真实截面是名义的 1.70 倍），八臂的按边表声明、两个 k-ω SST
伙伴、`starts`、容差 2e-4、6 趟 sub-iteration 预算全部照抄 §5.31；**`git diff` 里没有任何 `src/`**。
十族 **332** 条断言（`tests/check_transient_coupled_multi_both_parallel.py`，RC=0）：

* **两条边都在动**：八臂里每条边的速率都待在注入量之内且都走完一段路——离 `Q_IN` 最远
  **0.397 ~ 0.706**，跨层行程 **0.0178 ~ 0.1026**（`bore` 臂动得最少也有 1.78% 注入量）。X1 末层
  4.8743e-05 / 4.8473e-05 / 4.9500e-05 / 4.7591e-05 m³/s（四张表），X2 = `Q_IN` − X1，两条支路之和
  逐层精确等于注入速率，**没有第三条支路可以藏不平衡**。
* **`worst` 真的换手**：279 趟里打印的 `worst` 下标与两条残差的 argmax **逐趟一致**（0 处不符）；边 1
  当最差的次数 **11 / 7 / 7 / 4 / 13 / 8 / 8 / 9**（八臂），§5.31 的 224 趟这个下标全是 0。
* **一个节点压降、两份动量行**：同一臂里 `J1 − J2` 既要对上 X1 那行（用 **X1 那份** `M`）也要对上 X2
  那行（用 **X2 那份** `M`）——八臂最差 **1.1e-10 / 7.6e-11** 相对。`dpNetwork == R·q` 在 192 个
  边×层仍最差 **0.0 Pa**；两条主渠道仍对公开 Churchill 弦 + 自己的惯性闭合到 **2.0e-16**。
* **分配可以从外面重算**：`q_e = (dp + M_e q_prev,e/dt)/(R_e + M_e/dt)`、`Σq_e = Q_IN` 用**每边自己
  那本账**重建每个时间层的分流，八臂 96 层最差 **4.7e-12** 相对注入量；把任一条边的 `M` 换成另一本
  账，重建立刻偏 **0.047 ~ 0.108**（比真值差 10 个数量级）——"这条边吃哪本账"是这次运行的可识别属
  性，不是贴标签。
* **斜坡阶梯（这个拓扑才有的判据）**：四张表的分子 `M_X1 + M_X2` = **1.65495e7 / 1.43635e7 /
  1.39359e7 / 1.175e7** kg/m⁴，拟合 `τ` = **0.29115 / 0.25287 / 0.25075 / 0.21374** s（平谈），两种
  算子各自严格单调，且**每一对 `τ` 之比等于对应分子之比到最差 3.40%**；同一张表换算子 `τ` 不动
  （≤0.25%）。两条混合表**严格落在两张整表之间**——这正是 §5.31 边界①预言"两边都并联时混合表不再是
  任何一个整表"，现在是实测。单一速率判据：一条指数穿过全部 12 层对历史的偏离 ≤**0.258%** 行程
  （尾部四窗的带宽 0.2 ~ 10.1% 只作诊断报告，因为每个窗是两层小差之比，会放大末位数字）。
* **目的地仍然不含惯量**：八臂都弛豫到注入量的 **0.409319 ~ 0.409364**（spread **0.0206%**），且各臂
  都等于**自己那两条边报出的幂律** `dp ∝ q^m`（X1 **1.6601 ~ 1.6861**、X2 **1.6703 ~ 1.6842**，两条
  边都可读）的电阻分流，最差 **0.070%**。
* **没有乘客**：八对（4 表 × 2 边）的残差地板，两种算子之比落在 **0.044 ~ 25.6**，而且**方向不一**
  致（8 对里 5 对平谈更差）；每对最差地板不超过容差的 **21.1%**。对照 §5.31 那条被钉住的边：**单
  向** 73× ~ 87×，因为它无所谈、只是在读另一条边的串话。这里两条边都无所遁形。
* **算子轴仍只买到日程，块侧仍自己闭合**：趟数 **32 / 38 / 34 / 37 / 32 / 37 / 32 / 37** = **279**
  趟界面、**558** 次块解 == 趟数 × 2、**96** 个谈拢的时间层、每层 1D 自己解到 `<1e-11`；X1 `y+`
  **40.1 ~ 41.4** / `Re` 9518 ~ 9900 / 壁面剪闭合 **+0.20% ~ +0.33%**，X2 `y+` **37.3 ~ 38.1** /
  `Re` 11750 ~ 12068 / **+4.35% ~ +4.38%**，全部在 log-law 带内。

**成本**：checker 含跑 app，ctest 单条 **377.5 s**（八臂串行）。

**本增量的边界**：① 伙伴仍**准稳态**——块内 `du/dt` 没有进界面，斜坡是网络侧的 1/ dt 项给的。
② 界面仍是**每边一个标量电阻**，profile 级映射没做（§5.34/§5.35/§5.36 已补上）。③ `τ` 阶梯按**实测比**定判据（3.40%）而不是解析式
 `τ = ΣM/Σ(m·R)`（把两支路线性化成 DAE、层内不迭代时的应有值）：实测弛豫率系统性**偏快**，解析值
 是实测的 **1.479 ~ 1.533** 倍（实测比它小 32.4% ~ 34.8%），而这个倍数在四张表、两个算子下几乎不
 变——所以它像是一层内界面松弛（每层 1 ~ 6 趟）本身的离散映射效应，而不是归属表的效应；**机理没有
 在本增量里定死**，只记成实测。④ §5.31 边界②（平谈下 `network` 源那条边被另一条边的表碰到 2.9e-10）在这个
拓扑里**没有对应物**——两条边都能谈，那条被钉死的边不存在，其机理仍未定。⑤ app 头部的 `interface
( … )` 与 `inertia …` 两行仍只描述第一臂，混合臂要看它自己的 `transient arm … edge … inertia …` 行。
⑥ §5.31 边界①的另一半——**两条边都在串联位置**——仍未测；那里两个速率都被拓扑钉死，接口只剩谈压力，
本增量的判据（换手、分流重算、`τ` 阶梯）都没有对应物。

### 5.33 块内真实非定常 `du/dt` 进界面：界面仍在摩擦弦上谈、加速度记账给动量行，归属第一次带上可测误差（§5.32 的拓扑 + 六臂，含内核改动）

§5.25 到 §5.32 十一段瞬态耦合共用同一句前提：**伙伴按瞬时求**。块在被交接的速率下解到定常，于是它报的
压降只是速率的函数，刷新后的系数永远是两条正分支的交点——这正是接口能收敛的原因，也是**行里那份惯量
从来没有对手**的原因：`network`、`partner` 和逐边混合只是同一句界面条件的三种写法，三者之间唯一的证据
是各自走的路（趟数、`τ`、首增量踢差）。本增量把这块前提去掉：块在自己的网格上解 `du/dt`（隐式 Euler
一步对着上一层留下的场，同一个 SIMPLE 回路迭代到不动点），它报的压降从此来自一条**带加速度项的动量
平衡**。界面**继续在摩擦部分上谈**——`R ← dpFriction/q` 同时喂 GS 刷新与 QN 残差向量，因为用整段压降
刷新就是 §5.19 那条 `a − c/R` 映射，`dt < M/R` 时交点直接消失——而把块实测到的加速度**记进 1D 动量
行、只评分不迭代**：

    rowError = (R*q + M_row*(q - qPrev)/dt - dpBlock) / dpBlock

这条量**不可迭代**：界面迭代只会把摩擦残差压下去，`rowError` 恰好等于行里那份惯量错的数量。拓扑、两个
网格、两个 k-ω SST 伙伴、容差 2e-4 全部照抄 §5.32，六臂 = 三个来源（新增第三本账 `slab` =
`rho*V_probe/A^2`，就是被报压降真正跨过的那块流体）× 准稳态/时间准确 × 两个算子；控制臂的块字典照样写
`ddt Euler`、由 case 级 `couple{ ddts ( … ) }` 覆盖回 `steady`，所以覆盖两个方向都在跑，控制臂就是 §5.32
的 `bore_gs` 那一套计算。十族断言（`tests/check_transient_coupled_unsteady_partner.py`，RC=0）：

* **界面确实只在弦上谈**：24 个边×层里 `R*q` 落在摩擦弦上最差 **1.987e-04**（容差 2e-4 的 1.0×），而离
  **整段**压降差 **22× ~ 3039×** 容差、最坏一次差到压降的 **61%**——这就是"拿加速度刷新系数"那条路，量化
  了它为什么不能用。准稳态控制臂里弦**就是**整段压降（两者都在容差内，1.672e-04 / 1.655e-04）。
* **拆分是块自己场里的恒等式**：`dpFriction + dpInertia == dpBlock` 在 5 个非定常臂 120 个边×层
  ≤**1e-12**；`rowDrop == R*q + M*dq/dt`、`rowError == (rowDrop - dpBlock)/dpBlock` 最差 **0.0**（就是
  界面自己写下的那几个数）；`rowError*dpBlock` 对 `(M_row - inertanceDynamic)*dq/dt` 最差
  **1.94e-04** 相对压降，即容差之内——评分里剩的松紧是求解器的，不是归属的。控制臂 `dpInertia`、
  `rowError` 都恰为 **0.0**。
* **时间层拥有自己的历史**：`qPrev` 逐层等于**上一个被提交的时间层**的速率，最差 **1e-15**，首层两条边
  都是 0（从静止起步）；`dq` 就是两者之差。层内多趟因此都在同一份历史上加速，而不是沿轨迹走步。
* **加速度不是装饰**：它占报出压降的 **0.4% ~ 60.8%**，最大就是起步那一层，且随斜坡沉降而衰减（末层
  1.1%/0.4%）。符号在每臂**只有一处**例外——X2 第 2 层：那条支路的通量已经翻负（−3.2e-07 ~
  −1.9e-06 m³/s）而它的场仍在加速（+61 ~ +123 Pa）。这就是"块不是一个速率函数"最直接的读数：边界层还
  在长，几乎不变的通量可以带着仍在上升的动量。
* **块给自己的加速度评分（`chi`）**：`inertanceDynamic / inertanceProbe`——起步层 **92.8% ~ 94.6%**（少
  的那几个百分点是剖面刚出生）；第 4 层起 **|chi − 1| ≤ 2.8%**，五臂、两边、九层全都一样，而且同一
  (边,层) 上五臂之比 <1.05（第 4 层 X1 **1.016 ~ 1.028**、X2 **0.978 ~ 0.997**）；中间第 2~3 层
  **|chi − 1| > 5%**（实测 37% ~ 965%），那两层集总块什么都说不了。`betaProfile` **1.0163 ~ 1.0185**，
  `dtLevel` 逐位等于声明的 0.05 s。
* **归属第一次带上误差**：在沉降层（4~12）`rowError` 恰等于 `(M_row - M_probe)*(dq/dt)/dpBlock` ——绘图
  那份比 **0.968 ~ 1.015**（它比 slab 错 **+96%**（X1）/ **+191%**（X2）），区域那份 **0.971 ~ 1.012**
  （错 **+54%** / **+71%**）；`slab` 源同一式预测 0，实测残留 ≤**0.23%** 压降（X1 0.22/0.23%、X2
  0.01%），也就是 `chi` 带剩下的形状而不是错惯量。中位数阶梯（层 4~12）：X1 绘图 **1.84%** / 区域
  **2.53%** / slab **0.0388%**，X2 **1.50%** / **1.37%** / **0.0115%**——前两者**不按错的倍数排序**，
  因为读孔径同时也把轨迹挪了、被乘的那个 `dq` 也变了；两者都比 slab 高一个数量级以上。
* **评分与界面映射无关**：同一格的平谈/准 Newton 中位数之比 **1.01× / 1.00×**（绘图）、**1.29× /
  1.56×**（slab）——`rowError` 是物理性质，不是算子产物。
* **首层的弦是起始壁面剪**：非定常伙伴第 1 层的摩擦弦是准稳态同位置弦的 **2.46 ~ 3.05** 倍（10 个
  边×臂，跨度 <1.6 倍，所以这是一件事不是某臂的怪脾气），第 3 层起回到准稳态律 ±**3.8%** 内——被谈的
  系数重新成为速率的函数，这才是这套格式能收敛的原因；代价由下一层用趟数付。
* **目的地仍由弦决定，但位置第一次随归属变**：六臂末端 X1 分流 **0.400006 ~ 0.406190**，对照控制臂自己
  幂律的电阻分流 **0.409367** 偏 **0.78% ~ 2.29%**（<5%）；三个来源末端三条不同速率（X1
  **4.8423e-05 / 4.80008e-05 / 4.85344e-05** m³/s，两两差 >1e-9）。准稳态链里归属只能挪"路"，在这里它
  连"位置"一起挪——而上面的阶梯说哪个位置更对。
* **成本与块侧闭合**：趟数 **32 / 41 / 49 / 56 / 42 / 39** = **259**、块解 **518** == 趟数 × 2、**72** 个
  谈拢的时间层，每层最多 **4 / 6 / 6 / 9 / 5 / 5** 趟（预算 10；§5.32 的 6 趟预算在这里会让 `slab_gs`
  第 2 层以残差 2.2e-3 谈不拢——首层那条被起始壁面剪撑大的弦需要 5~9 趟才解得开）；每次块解的扫数反而
  **更省**：X1 控制 **11.8** vs 非定常 **7.8 ~ 8.4**、X2 **7.4** vs **5.5 ~ 6.3**（隐式 Euler 把 `V/dt`
  放上对角，`rAU` 被 `dt` 兜住）；1D 每层自己解到 `<1.7e-11`；两块场闭合——`y+` **37.6 ~ 40.9**、
  `Re_bulk > 4000`、`|dp/dx|*H/2 == rho*u_tau^2` 8% 内、`fluxIn == fluxOut == 末端 q`、`dpPa ==` 界面
  末端 `dpBlock`（1e-9）。

**成本**：checker 含跑 app，单条 **约 36 s**（六臂串行）。

**本增量的边界**：① 界面仍是**每边一个标量电阻**，profile 级映射没做（§5.34 已把这一格补上：界面现在
交得出形状，但 1D 那一行仍不认识 `beta`）；加速度整笔进 1D 那一行，没有回
传成块内的第二边界条件（一次 pass 内块的 `u_old` 固定）。② `slab` 是第三个**几何**账，不是形状账：
首层那 5% ~ 7% 的 `chi` 亏与第 2~3 层"通量与动量不成比例"没有任何归属表能表达，那需要 `beta` 本身进
行——记在这里但没有本增量的判据。③ 阶梯判据要**人工挑沉降层**（4~12）：第 2~3 层的 `rowError` 与
`(M_row - M_probe)*(dq/dt)/dp` 不同阶，混进中位数会把绘图与区域压成同一档（§5.32 那种"整段 12 层取
中位"的写法在这里给出 X2 上 1.99% 对 1.96%，不可分辨）。④ 首层弦膨胀 **2.46 ~ 3.05** 倍只记为实测，
没有做起始层显式模型；其机理（壁面剪还是 Stokes 的那份）与 §5.32 边界③的 `tau` 系统性偏快是两个独立
的、尚未定死的观测。⑤ 伙伴仍**单相不可压**、无热力学、无自由面。⑥ `ddts` 表是**逐臂**的（长度 == 臂
数），不是逐边——同一臂里两块要么都步进要么都定常；逐边 `du/dt` 与逐边 `ddt` 方案还没做。⑦ `slab`
没有静默回退：probe 窗内没有流体时直接 fatal（src/coupling/src/coupledNetwork.cpp:408），所以"第三本
账"要么真跨过被测量的那块流体、要么这次运行不成立。

### 5.34 profile 级界面映射：界面对象交的是一个速率，交的是这个速率**长成什么形状**（§5.20 的多边稳态接口，含内核改动）

从 §5.18 的第一段耦合到 §5.33 的十一段瞬态耦合，界面上走的都是**一个数**：网络要 `q`，块把 `q/A` 均匀
铺满入口面。对一个长到忘掉了入口的管段这是真话；对一块**正因为不够长才要被详细化**的区域是假话——同
一个速率铺成塞状和铺成 Poiseuille 要在不同长度上发展，因此代价不同压降，于是两侧谈拢的那个不动点本身也
不同。此前接口没有任何方式说得出某块区域正被喂成哪一种。

这一格现在补上：`detail{}` 的每条边可选 `profile`，三个词 `flat` / `parabolic` / `power1over7`，而且
`starts{}` 的每一**臂**各自命名（缺省即 `flat`，一个词等于两条边同形状，一个列表就是逐边）。六臂 =
两种塞状起点（系数相差三倍）+ 抛物 + 1/7 + `mixed`（一次运行里两条边喂不同形状）+ 抛物走准 Newton。
拓扑、两个网格、`relax`、容差 1e-4 与 `maxOuterIterations 12` 全部照抄 §5.20。

内核侧三处（`src/3d/src/fields/volFields.cpp:209`、`src/3d/src/solvers/simpleSolver.cpp:185`、
`src/3d/src/solvers/pisoSolver.cpp:128`）与耦合侧一处（`src/coupling/src/coupledNetwork.cpp:418`
`buildInletProfile`）：

* **速度边界可以逐面 prescribe**。`volVectorField` 的 `fixedValue` 原来只有一个 patch 级的数，现在多
  一张按面的覆盖表；`linearise` 读它，`surfaceScalarField::correctBoundaryFlux` 用
  `dot(U.faceValue(f), Sf)` 积分通量。`buildInletProfile` 按声明的形状给每个入口面定权、再按面积分
  归一化，使 `∫u dA = q` **精确**成立——权重只沿被无滑移壁面夹住的那根轴取，对称面那根轴没有可对抗的
  梯度所以保持均匀。
* **`HbyA` 必须逐面继承这份 prescribe**（`adoptFaceValuesFrom`）。第一版只让方程被形状强迫、预测通量
  却仍是塞状，边界那一列留下一个形状自身大小的失衡：抛物入口的 X2 连续性地板钉在 6.46e-2，400 扫不收
  敛。修好后同一算例每臂每次块解都自己停机（21~44 扫，首趟冷解 188~222），全案从 206 s 降到 ~170 s，
  顺带把 `betaIn 1.2209484` 那处 X2 抛物异常也消成了 1.2009419——它从来不是物理，是这个 bug。
* **外停机的 `settled` 参照改成"两次扫之前"**（`simpleSolver.cpp:378`，原为"历史最小值的 1.01 倍"）。
  两个理由都被实测撞上：冷启动的前两个二循环会先跌破它最终停住的地板（iter 3 的 0.0625734 对地板
  0.0682594），而从另一个速率热启动会**从地板下方**起步（`flat_lo` 第 2 趟 X1：0.0543008 一路降到
  0.0470024 而 `dp` 已到 1e-7）。历史最小值被这两件事同时污染，"iter≥5 才开始记"只救得了前者。二循环
  是这个循环自己的周期，所以新写法是相位对相位，而增长仍然表现为增长。§5.5 关于"预测地板是 HbyA/压力
  通量的固有失衡量、`|dp|` 分支才是真实外收敛信号"的结论照旧成立。

十族断言（`tests/check_coupled_1d3d_profile_map.py`，296 条，RC=0）：

* **交出去的是什么**：六臂各自按臂、按边报出 `inletShape`，`mixed` 一臂里 X1 是 `parabolic`、X2 是
  `flat`——形状是逐边的选择，不是全局的。
* **速率活过映射**：`fluxImposed == 网络要的那个数` 最差 **1.2e-15** 相对、`fluxImposed == fluxIn`
  **逐位相同**，而且每个横截面通量的 min/max 之差 ≤**2e-5**——面积分归一化不是近似。注入量与参考压力各回到
  1e-9。
* **形状就是声明的那个形状**：块从**自己的边界面上**读回的 `betaImposed` 对 checker 独立重算的中点求和
  （`parabolic`/20 面板 → 1.1970121617640874、/16 → 1.1953421565609941；`power1over7` → 1.0134305810983724
  / 1.0128704004797027）一致到 **1e-9** 相对；离连续极限 6/5 与 64/63 各差 0.25% ~ 0.39%（离散面板的
  固有偏差，不是误差）；`flat` 离 1.0 的绝对偏差 ≤**1.5e-15**。
* **动量通量台账**：`momentumFluxPa == rho*beta*q^2/A` 在 12 个边×臂上**逐位相同**——这就是 §5.33 边界②
  那句"形状本身还没写进行里"目前唯一可读的那个形状量。
* **还剩多少路要走**：20 个站从 `dx/2` 排到 `L - dx/2`；塞状入口 X1 从 **1.0245788** 起步、单调爬到平台
  **1.1955544**（差 0.171），1/7 从 **1.0558** 起步走同一条路但更短；抛物入口在承载边界的那一平面之外
  **漂移 ≤0.16%**（X1 0.154%、X2 0.162%）——Poiseuille 入口是平移不变的，这一条本身就是形状被真喂进去
  的证据。首平面读 1.2009419 对 `betaImposed` 1.1953422（0.47%）：那是入口单元自己的偏置。
* **平台是区域的，不是入口的**：三种入口下的发展 `betaOut` X1 **1.1955544 / 1.1955896 / 1.1955600**
  （相对极差 **2.9e-5**）、X2 **1.1930892 / 1.1930909 / 1.1930896**（**1.4e-6**），都贴着连续 Poiseuille 1.2。
* **阻力与它自己的物理**：对 `12 mu L_meas/(H^3 W)`（X1 3.70e6、X2 1.568e7，`L_meas` 就是块量压降的那
  一段）——抛物入口 **−1.81% / −2.48%**（壁面那一面上的梯度只跨半个单元，偏低是离散性质）、塞状入口
  **+22.97% / +14.08%**、1/7 **+19.03% / +11.27%**，且 `R_flat > R_17 > R_parab` 两条边都成立。每米压降对
  发展斜率：抛物 ≤**0.61%**，塞状 +24.3% / +16.6%。probe 窗内散射 X1 **0.4387 → 0.00422 Pa**（**104×**）、
  X2 **0.4385 → 0.0192 Pa**（**23×**），1/7 落在中间。
* **工作点因此会动**：只换形状，并联分流 X1 从 **4.9905859e-4** 到 **5.3998308e-4**（**+8.20%**），而串联
  边恒等于注入量（最差 **7.9e-15** 相对）；J1 的体积衡算 `q_X1 + q_By == q_A` 最差 **4.3e-14** 相对；旁路那条 1D 管用自己
  的 Churchill 1977 弦落在块的压降上最差 **5.0e-14**，前管 `Src − J1` 对同一根弦 **3.8e-13**；两域
  `Re` 落在 399 ~ 432（X1）与 427（X2），层流，所以臂间差异只可能是入口。
* **什么动不了它**：起始系数相差三倍的两臂落在同一个 R（X1 4549846.57 / 4550046.96，差 **4.4e-5** 相对）；
  `mixed` 的 X1 与 `parab` 的 X1 **逐位相同**（`dpPa` 1961.8431819560024 两边一样）——一条边的形状不影响
  它的伙伴；平谈与准 Newton 对抛物压降差 **2.0e-7** 相对，只买到趟数（3 == 3）。
* **谁为映射付账**：六臂 `iterations` **5 / 5 / 3 / 5 / 3 / 3**、`blockSolves == 2 * iterations`、SIMPLE
  解的次数与之一致、`worstResidual` ≤**3.4e-5** < 容差 1e-4，`settled` 行与块报告逐位相符且每块都在自己
  的门上停机（21~44 扫，无一次撞 400）。最后把**预测器地板**按边归位：抛物入口 **1.997e-3 / 3.597e-3**，
  比塞状入口同一边的 4.813e-2 / 5.430e-2 低 **24× / 15×**——形状进了 `HbyA`，而不是只进了校正。

**成本**：单臂串行整个 case ~166 s（六臂、两网格），含跑 app 的 ctest 单条 **166.17 s**；两网格各
1500 / 3200 单元。

**本增量的边界**：① 形状目前只到**块的入口边界条件**这一层，还没有成为界面自身的未知量：1D 那一行仍不
认识 `beta`（§5.33 边界②那格没动），也就是说接口交的是形状、但还不为形状的动量后果负责。② 三个词是
**解析**形状；块自己解出来的形状不会被回写成网络的入口条件，也没有从伙伴场读回的重映射通道。③ 只做了
**稳态**多边（瞬态侧的形状见 §5.35）。④ 湍流伙伴下的形状没验证——`power1over7`
是按层流连续极限选的对照物，RANS 渠道的 `beta` 与它的关系另是一格（§5.36 已量：那一格里 RANS 的平台是
**1.0206** 而不是 6/5，抛物入口因此从"欠满"翻成"过满"，判据变成 `(beta_in − beta_out)` 的符号）。⑤ 壁面那一面的梯度偏置（−1.8% /
−2.5%）记为离散性质、没做 `correctedSnGrad` 之外的补偿。⑥ 逐面的 prescribe 现在只有 `fixedValue` 读；
其他 patch 型仍走 patch 级的那一个数。

### 5.35 profile 级界面映射进瞬态接口：形状在每个时间层里重新交一次，于是它同时决定分流的**去处**和**快慢**（§5.34 的形状 × §5.25 的瞬态多边，含耦合侧改动，三维内核一行没动）

§5.34 把形状接在稳态多边上，那里形状只有一件事可说：区域要跑多远才算发展完全，因此它的阻力值多少。
瞬态里同一个数有**两**件事要说——除了去处，阻力还是斜坡**松弛到**的那个分母，也就是
`tau = (M_X1 + M_By)/(R_X1 + R_By)` 里的两项之一，所以形状还决定这处分流**多快**落位。而在此之前，
没有任何一层能说出某块区域在这一步被喂成了哪一种：瞬态接口逐层交的仍然只是一个 `q`。

拓扑、两个网格、驱动（`t=0` 开一个恒定注入率 `3.0e-4 m3/s`）、`dt 0.05` × 12 层、容差 1e-4、
`maxSubIterations 8` 全部照抄 §5.25 的 `transientCoupledMulti1d3d`，两条 duct 的几何也照抄 §5.34 那
两块，只是沿管向更粗（X1 30×20×1 / X2 25×20×1，probe 窗各 0.12–0.49 与 0.15–0.38）。五臂 =
`lo`（两边塞状，对照）+ `parab` + `p17` + `mixed`（X1 抛物、X2 塞状）+ `parab_qn`（抛物走准 Newton），
全部用未加形状那一版引用的同一对起点系数（`X1 1.25e6`、`X2 5.3e6`，各低于自己 duct 值三分之一）。

内核侧三处（`src/coupling/apps/transientCoupledMultiApp.cpp:170`、
`src/coupling/include/cf/coupling/CoupledNetwork.h` 的 `ProfileLedger`/`EdgeExchange`、
`src/coupling/src/coupledNetwork.cpp:1203`）：

* **瞬态 app 认得 `profile`**：与稳态侧逐字同款的解析（一个词 = 两条边同形状，一个列表 = 逐边），逐臂
  逐边在 `transient start` 行报出，缺省仍是塞状。
* **每层的边状态带一份形状台账**：`DetailBlock::profileLedger()` 把 `measure_` 里的
  `inletShape`/`fluxImposed`/`betaImposed`/`betaIn`/`betaOut`/`inletArea` 交给 `EdgeExchange`，瞬态的
  逐趟循环在块解之后填好，于是 `transient step n: edge i` 一行里读得到"这一层要的形状、该形状在入口面
  上被读回的值、以及区域自己的首末截面读到的值"。稳态那条路径直接读 `measure()`，一行没动。
* **`fluxImposed` 是新加的那个数**：它是界面写下的那个条件自己的面积分，因此它等于这一层的 `q`——这一
  条就是把"交了形状而不是交了另一份流体"和"形状顺带改了流量"分开的判据。

十族断言（`tests/check_transient_coupled_profile_map.py`，131 条，RC=0）：

* **形状是逐层的声明**：`120 = 5 臂 × 12 层 × 2 边` 每一条都报出该臂该边的形状——界面在**每一层的每一
  趟**重写它，而不是一次运行开头交一次。
* **速率活过映射，逐层都活**：`fluxImposed == 该层解出的 q` 最差 **6.66e-16** 相对；`betaImposed` 每层
  都回到自己声明形状的解析锚点（`flat` 1、`parabolic` 1.197012161764087、`power1over7`
  1.013430581098372），最差 **1.44e-15**，而且第 1 层与第 12 层是同一个数（层间极差 ≤**1.55e-15**）——
  形状是声明，跟着速率走的只有它的归一化。
* **交进去的形状不等于区域内读到的形状**：塞状 X1 以 1 到达、首个求解截面读 **1.09592**、末截面
  **1.19557**（区域还得自己发展）；1/7 以 1.01343 到达、首截面 **1.12133**（比塞状更满——它的加权取在
  发展截面上而不是边界面上）；抛物以 1.19701 到达、首截面 **1.20274**、末截面 1.19556（区域被**往下修**
  而不是往上长）。X2 同构（1.1253 / 1.14453 / 1.20223 对 1.19555）。
* **两端不对称**：末截面在五臂之间几乎不可分辨（X1 极差 **6.86e-6**、X2 **1.01e-8**），首截面却分开
  0.107（**9.7%**， outlet 极差的 **15566×**）与 0.077（**6.8%**，**7.6e6×**）；那一个平台值 **1.19557**
  对连续平面 Poiseuille 的 6/5 差 **0.37%**，是有限长管在这套网格上留下的量。
* **台账是每层的数，不是每案的常数**：塞状臂 X1 的首截面因子在第 1 层到第 12 层之间移动 **0.004177**
  （同一份几何上斜坡把更多速率推过去），抛物臂只移动 **4.984e-05**（**84×** 更少）——到货即成形的东西
  没有可发展的余地。
* **去处随形状动**：斜坡末端 X1 分流 `lo 1.937475e-4` → `p17 1.94691e-4`（**+0.4868%**）→
  `parab 1.991039e-4`（**+2.7647%**），三种形状是三个不动点而不是同一不动点的三种写法；串联边 X2 每层
  恒等于注入量（拓扑定的），它换的是阻力：`1.63217e7` → `1.61353e7`（**−1.14%**）→ `1.5317e7`
  （**−6.16%**），`mixed` 的 X2 回到塞状值。
* **形状是逐边的，而且逐边到逐位**：`mixed` 的 X1 与 `parab` 的 X1 在 12 层 × 6 个关键量上**逐位相同**
  （X1/By 组成电流分配器，下游任何区域动不了它的分流）；这一非平凡性由同一臂的 X1 与 `lo` 在同层相差
  **2.76%** 撑住。`mixed` 的 X2 跟着**塞状**臂到 **2.58e-6** 相对，而离"把自己的 X1 塑形成抛物"那一臂
  差 6.6%——不到逐位是因为界面在容差处停下，而每层的热启动取决于此前花了几趟。
* **网络那侧一行方程都不为形状改账**：五条管在每个节点压降上都以 `R*q + M*dq/dt` 闭合，最差
  **4.46e-11**（每臂都测）；`M` 在每层每条边都还是 `inertias ( network )` 放的 `rho*L/A`
  （X1 1.27324e6、X2 1.81083e6 kg/m4）。`unsteady` 在对照臂的 12 层与每臂的末层都是 0（两块都是
  `ddt steady`）：这一格测的是形状，不是区域自己的历史，`rowError` 在此不评分。
* **快慢也随形状动**：第 1 层严格等于纯惯量分流（≤1e-9）；把每臂自己的轨迹按离散惯量解
  `1/(1 + dt/tau)` 回归——`lo` 拟合衰减 0.833970 对预测 0.837734（**−0.449%**）、`parab` 0.842224 对
  0.842082（**+0.017%**），全 12 层对斜坡的偏离最差 **0.0747% / 0.0033%**；`tau` 由该臂落定的阻力构成：
  **0.258135 s → 0.266619 s（+3.29%）**（`p17` 0.259611 s 落在中间）——更满的入口落得更慢，因为它要走的
  路也更长。到第 12 层为止抛物臂离终点 **1.99%**、塞状臂 **1.49%**。`parab_qn` 落在 `parab` 的轨迹上
  （末层 1.991039066e-4 对 1.991039023e-4），算子没动形状定下的任何事。
* **谁付账**：五臂的界面趟数 **25 / 17 / 25 / 17 / 17**，块解次数恰为 `2 ×` 趟数（每趟刷新两条边，包括
  已经 agree 的那条），全案 **101 趟 / 202 次块解 / 60 个采纳的时间层**。形状改变的不只是答案，还有
  时间表：抛物臂 12 层只花 17 趟；而准 Newton 在这个容差下**什么都没买到**（同样 17 趟），因为塑形后的
  界面每层本来就一两趟谈拢。

**落定的阻力对 duct 律**（`12 mu L/(H^3 W)`）：X1 塞状 **1.0570×**、1/7 **1.0396×**、抛物 **0.9594×**；
X2 **1.0201×** / **1.0085×** / **0.9573×**——塞状高（塞子还得被发展，弦为之付账），抛物低（壁面那一面
的梯度只跨半个单元，§5.34 已记为离散性质）。X1 probe 窗内的每米压降散布从 **0.4045 Pa** 塌到
**0.006182 Pa**（**65×**）；X2 只从 0.02243 Pa 塌到 0.01842 Pa（**1.22×**，因为它本来就短，窗内那一段
与入口发展长度无关的斜率占主导），所以这里承担结论的不是散布之比而是**每个块自己的窗宽**：它对**自己的**
压降始终 ≤**0.053%**（`lo` X1 768.01 Pa 的 0.0527%，X2 4896.5 Pa 的 0.0005%），界面读到的弦不是弯线上的
割线。两域 `Re`（按水力直径）**155.0 / 160.0**，全部层流
（`nutMax`、`yPlus` 为 0），这正是"形状之间可比"的前提。

**动量通量台账（可读，但在这个 `Re` 下不是压降的主体）**：`momentumFluxPa` 对
`rho*betaImposed*q^2/A` 每块每层一致到 **1e-9**，而它对区域报告压降的比值是 `lo` X1 **9.78e-5**、X2
**4.90e-5**、`parab` X1 **1.33e-4**、X2 **6.25e-5**。X2 在两臂承载**同一个**速率，所以它的入口动量通量
之差（0.287283 对 0.24 Pa）完全是声明的形状：比值 **1.197012 = beta(parabolic)** 到 **1e-9**。

**成本**：五臂、两网格、60 个时间层整案 **35.5 s**（§5.34 的稳态六臂版 ~166 s）。

**本增量的边界**：① 形状仍然只到块的入口条件；1D 那一行动量方程还不认识 `beta`（§5.33 边界②、§5.34
边界①同一格），所以动量通量只被**记账**，没进任何一行被求解的方程——上面那个 ≤2e-4 的比值就是它还不需要
进账的理由，也是它进账时会最先显形的地方。② `tau` 里的 `R_X1` 是**块落定后的割线**，不是界面正在谈判的
那个数，因此"形状决定快慢"目前是逐臂落定后的比较，不是对未收敛过程的预测。③ 湍流伙伴 + 形状 + 瞬态三者
同台没测（§5.34 边界④那一格在瞬态侧同样没动）——§5.37 已同台。④ 逐层的 `betaIn` 只报首末截面；窗内的完整序列在稳态版
已按站输出，瞬态版此处只在 checker 里比对。⑤ `mixed` 的 X2 到 **2.58e-6** 而不是 0，其因由（每层热启动
取决于趟数）已记在上面，没有把它做成一条可关闭的判据。

### 5.36 RANS 伙伴下的 profile 级界面映射：形状第一次交给的不是一个 duct 律，而是一块自己会决定"发展完成"长什么样的区域（§5.34 的形状 × §5.22 的湍流伙伴，含内核改动）

§5.34/§5.35 把形状交给的伙伴都是层流 duct：那里区域的-developed 剖面就是连续 Poiseuille，动量通量因子
6/5，所以"入口越满、区域要走的路上欠的账越少"是一句关于**剖面**的判断。RANS 渠道把这个前提抽走了——
它发展到的不是 1.2，而是本例实测的 **1.0206**（平坦核心区 + 薄壁面层）。于是两根轴第一次同台，能问出
两根轴各自问不出的那个问题：**"满入口更便宜"是剖面的性质，还是那块层流 duct 恰好有的那个发展终态的
性质。**答案是两者都不是：判据是 `(beta_in − beta_out)` 的**符号**，而这个符号在湍流伙伴下会翻。

算例 `tests/cases/coupled1d3dTurbulentProfile/case.dict`：拓扑照抄 §5.20（`Src → A → { X1 ‖ By } → X2 →
B → Snk`），两块网格与湍流设置照抄 §5.27 那两条 `kOmegaSST`（intensity 0.05、`relaxation 0.7`、
`div upwind`、`laplacian/grad corrected`），只是 X2 从 0.12 m 加长到 **0.24 m** 好让抛物形的入口段落在块
内。驱动 `1.2e-4 → 6.0e-5 m3/s`（原因见下面的边界①），`rho 1000 / mu 1.0e-3`，两条耦合边的 proxy 孔径
仍是各自截面的等效直径 `d = 2 sqrt(H W / pi)`（1.4e7 / 1.36e7 起点系数），`couple tolerance 1e-3`、
`maxOuterIterations 12`、`accelerates ( 0 0 0 0 1 )`，块侧 `simple maxIterations 300 / continuityTol 3e-3`。
五臂 = `lo`（两边塞状）+ `parab` + `p17` + `mixed`（X1 抛物、X2 塞状）+ `parab_qn`。

内核侧一处新工具（`src/coupling/include/cf/coupling/CoupledNetwork.h:209` 的 `pressureStations`、
`src/coupling/src/coupledNetwork.cpp:798`/`:814` 逐站填写、
`src/coupling/apps/multiCoupledApp.cpp:358` 逐臂逐边输出）：

* **压强与形状报在同一组截面上**。等截面 duct 里壁面剪切只能把压强往下拉，所以"某一站压强上升而上面
  那站下降"就是入口段在**把压降还回来**——这是 §5.34 那版唯一不能直接读到的那件事，因为层流下满入口
  的账恰好被"少发展一段路"抵掉，曲线从未抬头。β 序列本来就在这些站里，加一列 `p` 不改任何求解。

十族断言（`tests/check_coupled_1d3d_turbulent_profile.py`，272 条，RC=0，整案 **29.7 s**）：

* **形状是逐边的声明，伙伴是非线性的**：五臂 × 两边各自报出 `inletShape`，且 `turb kOmegaSST` 在场——
  这一格里被喂形状的是真三维 RANS 块，不是任何弦。
* **速率活过映射**：`fluxImposed == 网络要的那个数` 最差 **1.0e-12** 相对、与 `fluxIn` 逐位同、每个横截
  面通量极差 ≤**2e-5** 相对。注入量与参考压力回到 1e-9。
* **形状就是声明的那个形状，与谁 owns 壁面无关**：`betaImposed` 对 8 面板中点求和重算一致到 **1e-9**
  （`parabolic` 1.181719848566793、`power1over7` 1.0103203595676944、`flat` 与 1 差 ≤1.5e-15）；离连续
  极限 6/5 与 64/63 各 ≤**2%**。动量通量台账 `rho*beta*q^2/A` 每块逐位复现（X1 抛物 **0.0670 Pa**、X2
  **0.1329 Pa**——注意这个数在湍流压降面前是 4e-4 量级，它不进账没关系，见边界②）。
* **符号决定区域怎么处理它**：十个"臂×边"里 **5 个过满 / 5 个欠满**（`over = parab×2、mixed X1、
  parab_qn×2`）——层流那版抛物是**欠**满的。过满三者的 β 沿程**单调下降**（1.181720 声明 → 1.119795 首
  站 → 1.020266 平台，后四站离平台 ≤**2e-4** 相对），欠满者单调**上升**（塞状 1.00110 → 1.02104；1/7
  1.01384 → 1.02095，连第一个内部截面都已经比交给它的 1.01032 更满）。
* **平台属于区域**：十对 `betaOut` 落在 **1.020266 .. 1.021045**，极差 **7.79e-4**；它离层流 duct 的
  **1.2 差 14.9%**——所以"形状值多少"这件事在换伙伴时必须重新量，不能换算。
* **形状的账用压强支付，而且看得见**：`parab` X1 压强从入口面 **171.39 Pa** 升到 **391.10 Pa**
  （x = 0.045 m，**+219.7 Pa**）再一路降到 15.53 Pa；`parab` X2 从 750.29 升到 **1083.71**
  （x = 0.035 m，+333.4）。`lo` 与 `p17` 的四条曲线**逐站不升**（首站即最大）。这是等截面块里唯一能
  把"入口比区域自己的终态更满"读成数字的方式。
* **摩擦没有被形状动过，被动的是割线**：developed 窗内的每米压降对模型自己的壁面平衡
  `2 rho u_tau^2 / H` 落在 **0.943~1.107**（十个"臂×边"全测；`lo` 0.9449/0.9434、`p17` 0.9842/0.9870、
  `parab` 1.1067/1.0778）。而整块报出的 `dpPa` 对"按同一平衡铺满整条 L"的外推值是：塞状 **1.0118 /
  1.0104**、1/7 **0.9800 / 0.9835**、抛物 **0.2835 / 0.5370**（准 Newton 0.2833 / 0.5368）。抛物 X1 还给
  回的 **393.96 Pa** 是 `rho (beta_in − beta_out) U^2 = 225.7 Pa` 的 **1.75×**，X2 的 621.0 Pa 是 283.7 Pa
  的 **2.19×**——剩下的那一份来自入口段自己的壁面剪切低于发展值，本增量只把它记为差额，没有拆开。
  `R_flat > R_17 > R_parab` 的顺序在两条边都**保持**（X1 1.4078e7 > 1.3094e7 > 4.1397e6），但抛物那一档
  从层流的"便宜 8%"变成"便宜 **68%**"。
* **工作点跟着割线走**：只换形状，并联分流 X1 `2.439311e-5 → 3.765199e-5`（**+54.4%**）、1/7 **+3.32%**
  （比抛物小 **16.4×**，与它离平台只有 0.0071 相称）；串联边 X2 每臂恒等于注入量（**1e-12** 相对），
  J1 体积衡算最差 **1.2e-12**，旁路 1D 管以自己的 Churchill 1977 弦落在块压降上（parab：**155.92 Pa**）、
  前管 `Src − J1` 对同一根弦闭合；`J3` 五臂同值（极差 <**1e-6 Pa**），因为串联边下游什么都没有动过。
* **逐边归属与算子各自仍然成立**：`mixed` 的 X1 与 `parab` 的 X1 分流差 **7.5e-13** 相对、压降差
  **1.3e-4**（那一份是块自己的扫数 29 对 32，不是形状串味）；`mixed` 的 X2 回到塞状臂（压降差 **1.5e-7**）
  且 `betaImposed == 1` 到 1e-12；平谈与准 Newton 对抛物压降差 **4.2e-4 / 3.0e-4**，只买到趟数（5 → 9）。
* **时间表，以及容差为什么是 1e-3**：五臂 `iterations` **3 / 5 / 5 / 5 / 9**、`blockSolves` 与 SIMPLE 解
  次数都恰为 `2 ×` 趟数、`worstResidual` **3.94e-4 ~ 9.24e-4** 全部在容差内、`settled` 行与块报告逐位相符。
  唯一不在自己门上停机的是 **抛物 X2：每趟都花满 300 扫**（`blockConverged 0`），它在被拓扑钉死的速率下
  逐趟报出的 `dpBlock` 在 **719.35 ~ 720.30 Pa** 之间摆（**1.3e-3** 相对）——界面容差低于伙伴自己的分辨率
  只会烧趟数，这就是 `1e-3` 的实测地板（§5.22 因同一理由选了同一个数）。塞状/1/7/`mixed` 的六条边全部
  7~29 扫停机；每次求解首扫的连续性地板最差 **1.08e-2**。
* **壁面函数仍然合法**：`y+` 首格 **21.3 ~ 32.6**（>11，落在对数律带内，未进尾迹区），`nutMax/nu`
  **46 ~ 146×**，水力直径 `Re_bulk` **4879 ~ 10000**——两条边的湍流是真的，这正是"形状之间可比"的前提。

**本增量的边界**：① **操作点有一条实测的墙**：`1.2e-4 m3/s` 下塞状臂与 1/7 臂**照样收敛**（**7** 趟 /
**6** 趟），只有抛物臂死在**第 1 趟**里——串联块 X2 要把整个 1.2e-4 推进 4 mm 缝隙（3.75 m/s），它的连续性
残差到第 40 扫已经是 `1.9e48`（`dp 7.4e92`），第 **47** 扫场变成非有限值，块随即触发 "the detailed region
reports a non-finite pressure drop"；同一次运行里的并联块 X1 只是不停在自己门上（300 扫用满，残差平台
`0.0163498`）。把 `relaxP` 从 0.3 改到 0.1：X1 的前 **32** 扫逐位不变（`aP = min(relaxP, max(pRelaxStart,
0.15*pScale/pStep))`，该式不超过 0.1 时两个 relaxP 给同一个数，见 `src/3d/src/solvers/simpleSolver.cpp:308`），
第 33 扫之后 X1 反而在 101 扫停下，而死因、死在哪一块、第几扫全都一样。把 `model` 换成 `laminar` 死得
更早：并联块第 **8** 扫（`77 -> 1.03e6 -> 5.5e16`）。加密也救不了：法向 8 -> 16 格后 13 次块求解**全部**
花满 300 扫、残差在 `0.059 ~ 1.5` 之间摆，界面把 X1 的速率推到 `1.13e-4` 时块报出**负的**压降
`-0.4367 Pa` 而停机。所以这堵墙不是压力松弛、不是湍流模型、也不是网格能旋钮出来的，而是**往壁函数渠道
上写层流形入口**这件事在 1.2e-4 处的实测边界。本例因此落在 **6.0e-5**。② 1D 那一行动量方程仍然不认识 `beta`（§5.34 边界①、§5.35 边界②同一格）：这里
`rho beta q^2/A` 只有 4e-4 量级的压降占比，形状的全部后果目前都由块的**割线**代劳，所以它进账的时机是
"界面要预测而未收敛的过程"，不是本例。③ 入口段还回的 393.96 Pa 与动量项 225.7 Pa 之间那一份差额没有
拆开归因。④ 形状 + 湍流 + **瞬态**三者仍未同台（§5.35 边界③那格一字未改）——§5.37 已同台。⑤ 两条耦合边仍是一条并联
一条串联；抛物形状在串联位置会直接把注入量改掉，那是另一格。⑥ 发散边界（①）目前是**记录**，没有做成
可断言的判据——checker 只跑 6.0e-5 那一档。

### 5.37 形状 × 湍流伙伴 × 瞬态三者同台：稳态那条符号律在每个被采纳的时间层里重新成立一次，而形状第一次改掉的是斜坡的**快慢**本身（§5.36 的形状 × §5.35 的瞬态，含内核改动）

§5.36 的判据是 `(beta_in − beta_out)` 的**符号**，但那是在一次稳态求解里读的；§5.35 把形状放进了时间层，
它的伙伴却是层流 duct——developed 剖面就是连续 Poiseuille 的 6/5。两边各留下同一格疑问：**符号律是不是只
属于"落定后的那个操作点"**。瞬态里平行边的速率逐层在涨、块给回的是逐层的割线，如果符号会在斜坡中途翻，
§5.36 那条结论就只是端点巧合。三轴第一次同台，能问出第三件两根轴各自问不出的事：**入口段还压那笔账有没有
时间签名**。稳态里它是一条沿 x 的压强曲线；这里它必须在每一个被采纳的层上重读一次，才能说"还压"是层 `n`
自己要的东西，而不是这一臂结束时剩下的东西。

算例 `tests/cases/transientCoupledMultiTurbulentProfile/case.dict`：拓扑照抄 §5.34/§5.35
（`Src → A → { X1 ‖ By } → X2 → B → Snk`），两块网格与湍流设置照抄 §5.36 那两条 `kOmegaSST` 渠道（X1
20×8×2 over 0.2×0.005×0.005，X2 24×8×2 over 0.24×0.004×0.008，intensity 0.05、`relaxation 0.7`、probe 窗
0.12–0.19 / 0.16–0.23、`lengthScale` 3.0e-4 / 2.4e-4），操作点仍是 §5.36 的 **6.0e-5 m3/s**（1.2e-4 是那格
边界①里层流形入口死掉的地方）。时间层 `dt 0.05 / endTime 0.4` → **8 层**，注入速率恒定，所以斜坡是管网自己
对并联分流的**惯性弛豫**，不是时间表驱动的源。两条边的 proxy 孔径取各自截面的等效直径（A 2.5e-5 / 3.2e-5，
`coupledResistance` 起点系数 1.4e7 / 1.36e7），`couple tolerance 1e-3 / maxSubIterations 12 / accelerates
( 0 0 0 0 1 )`，块侧 `simple maxIterations 300 / continuityTol 3e-3`、`schemes { ddt steady; }`（块内非定常
是 §5.33 那一格，本例不动）。五臂 = `lo` + `parab` + `p17` + `mixed`（X1 抛物、X2 塞状）+ `parab_qn`，
`starts` 五臂共用同一对种子，臂间唯一允许的差别是声明的形状（和最后一臂的界面算子）。

内核侧一处改动，是 §5.36 那件工具的时间化（`src/coupling/include/cf/coupling/CoupledNetwork.h:40` 的
`ProfileLedger::pRisePa/pRiseX`、`:296` 在 `profileLedger()` 里把 `measure_.pressureStations` 归约成"最高站
压强减首站压强、以及它落在哪一站"、`:461` 的 `EdgeExchange` 字段、
`src/coupling/src/coupledNetwork.cpp:1224` 在瞬态接口每层每边的台账里填写、
`src/coupling/apps/transientCoupledMultiApp.cpp:486` 逐层逐边输出）：

* **上升量跟着求解走，不跟着运行走**。此前 app 是在整臂跑完后再读一次块的
  `measure().pressureStations`，那只能给出这一臂最后一层的那条曲线；把归约放进 `ProfileLedger`，层 `n`
  报的就是层 `n` 那次解自己的还压，于是 checker 能对 **80 个（层 × 边）状态**各自的高峰下断言，而不是对
  一条端点曲线。三维内核、界面谈判式、1D 那一行都没改：`pRise` 是**读**出来的，不进任何被解的方程。

十一族断言（`tests/check_transient_coupled_multi_turbulent_profile.py`，652 条，RC=0，整案 **64.5 s**）：

* **符号律逐层成立，不是端点巧合**：80 个状态的 `betaOut` 全部落在 **(1.02, 1.0213)**，实测
  **1.02027438 .. 1.02105813**——平台是**区域的**，不随层也不随形状动（比层流 duct 交回的 6/5 低 **14.9%**，
  1.020716 对 1.2）；入口比平台满的恰好 **40** 个状态、比它差的恰好 **40** 个，且与 `pRisePa > 0` 的 40/40
  **是同一批**，没有一个状态骑在边界上。`lo` / `p17` 在每一层都不比区域自己发展到的更满，`parab` / `mixed`
  的抛物边 / `parab_qn` 在每一层都更满。
* **仪器不会自己造出上升**：`pRisePa` 为 0 的那些状态报的站就是区域自己的第一站（0.005 m）——单调下降的
  压强，其最大值就是起点；非 0 的那些永远落在同一站（X1 **0.045 m**、X2 **0.035 m**），在每一层速率被抬上去
  的过程中都不换站。末层的 `pRisePa` / `pRiseX` 与该次求解自己的压强表**逐位相符**（10 条边 × 臂）。
* **还回来的就是动量那一笔，而且逐层**：抛物 X1 的上升是 `rho*(beta_in-beta_out)*U^2` 的 **0.9888 ..
  0.9980**（预测 136.19 Pa 对实测 134.66 Pa），抛物 X2 落在 **1.1748 .. 1.1754**（283.74 对 333.43 Pa）；
  比值不随斜坡漂，说明它是每一层的同一件物理，不是起动瞬态。上升量对速率的增长指数 **2.1659**（101.7514 →
  134.6606 Pa，而速率只涨 1.1381×），`parab` / `mixed` / `parab_qn` 三臂给到同一个数。被钉死的串联边则坐在
  **一个**比值上到 **6.3e-4**（层 1 与层 8 都是 333.4297 Pa）——它的速率从未动过，所以它的一切也不许动。
* **形状是逐边的声明，且在时间里仍然是**：每臂两边的 `inletShape` 在 8 层全部与声明一致，80 个状态报的正是
  该臂所名的形状（界面每层重交一次，不是在开头谈拢一次）；`betaImposed` 等于声明剖面在 8 面板上的重算值到
  **1.1e-15**、跨 8 层散布 ≤**9.4e-16**（`parabolic` 1.181719848566793、`power1over7` 1.010320359567694、
  `flat` 与 1 差 ≤1.5e-15，分别距连续极限 6/5 与 64/63 ≤2%）；`fluxImposed` 对截面积分重算最差
  **6.66e-16**——形状是用形状买的，不是用另一份流体买的；`momentumFluxPa` 对 `rho*beta*q^2/A` 到
  **2.2e-16**。
* **形状不进 1D 那一行动量方程**（§5.34 边界① / §5.35 边界① / §5.36 边界②同一格）：`rowDrop == dpNetwork
  + M*dq/dt` 在 5 臂 × 8 层逐位成立，`M` 就是 `rho L/A`（**8e+06 / 7.5e+06**，与块的 `regionInertance`
  相同；惯量归属这根轴本例固定在管网侧，是 §5.26 那一格），`dpInertia` 恒 0 / `dpFriction` 恒等于伙伴报的
  全量压降（准稳态伙伴不认识 `dq`）；残差是两份压降的相对不同意、且在容差内；先导管 Churchill 压降与
  `Src − J1` 到 **1.8e-12**，`Snk` 恒等于参考压力。
* **第一层是纯惯性分流**：`M_By/(M_X1+M_By)` = **0.427275**，实测 `lo` **0.424732**（−0.595%）、`parab`
  **0.438354**（+2.593%）、`p17` **0.426449**（−0.193%），次序仍是抛物 > 1/7 > 塞状——同样一记惯性踢下，
  更满的入口欠更少压降，于是分流给它更多。
* **τ 第一次被形状按**倍数**改掉**：`lo` **0.3869 s**、`parab` **0.7942 s**（比值 **2.0531**）、`p17`
  **0.3897 s**（+0.72%，即抛物那份移动的 **0.7%**）。三臂的逐层 `|dq|` 衰减都是几何的（0.858487 ..
  0.940144，层内散布 <5%），所以一个时间常数能描述整条序列。§5.35 对层流 duct 只移动 **3.29%**，因为那里的
  割线本来就贴着落定斜率；这里 X1 的落定割线被形状改到 **3.1×**（塞状 1.422271e7 > 1/7 1.314546e7 > 抛物
  4.614610e6），斜坡的快慢就跟着改。教材式 `(M_X1+M_By)/(R_X1+R_By)` 给 **0.5823 / 1.2163 / 0.6160 s**，
  对绝对值一律高估 **1.5053 ~ 1.5807×**，但抛物/塞状的**比值** 2.0887 对实测 2.0531 只差 **1.2%**——同一个
  1.5× 在两臂之间被除掉了。
* **去处也被形状改掉**：斜坡弛豫到的分支份额 `lo` 0.413869、`parab` 0.498902、`p17` 0.422758（`mixed` /
  `parab_qn` 与 `parab` 同）。抛物形在第一层只值 **+3.2071%**，在末层值 **+20.5460%**——形状的价格是在斜坡
  上付的，不在斜坡起点；1/7 是 +0.4041% → +2.1479%，与它离平台那 **0.0072** 的 β 差距相称。
* **可分性在时间里保持**：`parab` 与 `mixed` 的并联边 **7 个字段 × 8 层逐位相同（0.0e+00）**——串联边被喂
  什么，到不了并联边，即使两边伙伴都是 RANS；`mixed` 的串联边就是 `lo` 的那条（速率与压降相对差 <**5e-5**、
  `pRisePa` 在 8 层恒为 **0.0**）；串联边的速率在全部 40 个（臂 × 层）状态上都等于注入量（≤**1e-11**），
  所以形状在那里只能移动压降，永远移动不了速率；`parab` 与 `parab_qn` 差 <**1e-4**，算子改的是路径不是物理。
  `mixed` 一层之内 X1 抛物、X2 塞状，`betaImposed` 分别是 1.181719848566793 与 1——按边表在时间循环里是真表。
* **时间表，以及那堵诚实的墙**：五臂 **16 / 17 / 13 / 17 / 18** 趟，SIMPLE 解次数恰为 `2 ×` 趟数
  （32/34/26/34/36），最坏残差 **2.673e-06 ~ 6.040e-04** 全在 1e-3 内；每层的状态就是它最后一趟的状态，
  `worst` 下标与同行里残差最大的那条边一致（`mixed` 每一趟都停在被塑形的那条边上，`[0]`——平的串联边即使
  在自己谈压降也只是乘客）；管网线性代数本身每层到机器精度。唯一不在自己门上停机的是**抛物形的串联块**：
  `parab` 8 层里 **4** 层、`parab_qn` **6** 层花满 **300 扫**，且花满与 `blockConverged 0` 是同一批状态；
  同一个被钉死速率上塞状/1/7 **7 扫**停机、塑形并联边 **24 扫**。`mixed` 的串联边是平的，所以邻居抛物也
  不影响它——墙属于"喂给钉死速率的那个形状"，不属于串联拓扑。网络握着那个速率，所以未收敛的伙伴表现为
  **摆动的压降**而不是失败，这正是容差取 1e-3 而不再小的理由。
* **代价仍然落在割线与壁面那侧**：末层 X2 塞状 **1538.4003 Pa**、1/7 **1428.4563 Pa**、抛物 **719.9844 Pa**
  ——形状占掉这条边压降的 **53.2%**，而它的速率从未变过；抛物 X1 报 **138.11 Pa**，是发展外推 389.47 Pa 的
  **0.3546**，缺的 251.35 Pa 是入口弛豫项 136.19 Pa 的 **1.85×**（X2：0.5368 与 **2.19×**，缺 621.18 Pa 对
  283.74 Pa）；塞状/1/7 则报自己 fully-developed 身的 **0.9801 ~ 1.0124**（3% 内）。`pRise` 相对该臂报告的
  压降是 X1 **97.5%**、X2 **46.3%**。壁面函数仍然合法：`y+` 首格 **21.4 ~ 32.6**、`nutMax/nu` **47 ~ 134×**、
  按水力直径 `Re_Dh` **4966 ~ 10000**，发展窗每米压降对模型自己的壁面剪力 `2 rho u_tau^2/H` 闭合到
  **−5.7% ~ +10.1%**；每块入口进、出口出到 **≤5.2e-12**——形状是在截面上重新分配，不是被造出来的。

**成本**：五臂、两块 RANS 网格、8 个时间层整案 **64.5 s**（§5.36 的稳态五臂 29.7 s、§5.35 的层流瞬态 35.5 s）。

**本增量的边界**：① 块内仍解 `ddt steady`：本例的"瞬态"是管网瞬态 + 逐层准稳态伙伴；§5.33 那块真 `du/dt`
与形状还没在同一格碰过面，所以这里的还压完全由准稳态割线代劳，"层 `n` 的上升量"里不含块自己的时间导数。
② 1D 那一行动量方程仍不认识 `beta`（§5.34①/§5.35①/§5.36②同一格），但这一格的理由变了：`pRise` 已经是
134.66 / 333.43 Pa 量级，分别是抛物臂 X1、X2 所报压降的 **97.5%** 与 **46.3%**——形状后果里最大的一笔仍然只
在块的割线里体现，没进被求解的那一行。③ `τ` 的判据是**实测比**（2.0531），教材式对绝对值仍系统性高估
1.5×；其来源（"发展中的区域的落定斜率 ≠ 它的割线"）没有拆开归因。④ 抛物串联块花满 300 扫目前是**记录 +
一致性断言**（花满的那批与 `blockConverged 0` 那批同一），还没有做成"界面应当在几层内把未收敛伙伴的压降摆动
压到多少"的可关闭判据。⑤ 两条耦合边仍是一条并联一条串联；形状 × 双边并联（§5.32 的拓扑）与形状 × 非定常伙伴
（§5.33）这两格都没与本例交叉。⑥ 1.2e-4 那堵发散墙仍是记录（§5.36 边界①同一格，本例未重跑那一档）。

### 5.38 形状 × 块内真 `du/dt` 同台：区域的加速度对形状是盲的、摩擦半区不是，于是形状的价格整个落在界面谈判的那半区，而行侧的惯量声明被记了两次（§5.37 边界① × §5.33 的块内非定常，零内核改动）

§5.37 自己写下边界①：块内仍解 `ddt steady`，"层 `n` 的上升量"里不含块自己的时间导数，还压完全由准稳态割线
代劳；§5.33 把一块真 `du/dt` 的伙伴放进过时间层，但从未递给它一个形状。这两格碰面的地方，能问出一根轴各自
问不出的问题：**区域的压降现在是两笔之和 `dpBlock = dpFriction + dpInertia`，而界面只迭代第一笔——那么形状
到底住在哪一项里**。如果区域自己流体的加速度不知道入口怎么铺开，形状的全部价格就留在电阻被指派的那条割线
里，行侧的惯量声明塑形与不塑形是同一个数；如果它知道，形状就有一部分成本落在动量行里，§5.33 做的归属就
成了界面剖面的函数。

算例 `tests/cases/transientCoupledMultiTurbulentProfileUnsteady/case.dict`：拓扑、两块网格、壁面函数、操作点
与 8 个时间层全部照抄 §5.37（`Src → A → { X1 ‖ By } → X2 → B → Snk`，6.0e-5 m3/s，X1 20×8×2 over
0.2×0.005×0.005 / probe 0.12–0.19，X2 24×8×2 over 0.24×0.004×0.008 / probe 0.16–0.23，proxy 孔径取截面等效
直径所以行侧 `M` 就是块自己在网格上量到的 `rho*V/A^2`）。唯一改的是时间层与表述：`time { dt 0.2; endTime 1.6; }`
（§5.37 是 0.05/0.4），块 dict 两边都写 `schemes { ddt Euler; }`，`couple{ ddts ( steady steady Euler Euler
Euler ) }` 把前两臂**反向**覆盖回 `steady`，于是这个 override 两个方向都被走一遍。五臂就是一张 2×2 加一格
第三形状：`lo_steady` / `parab_steady` 是 §5.37 的那两臂，`lo_euler` / `parab_euler` 是交叉本身，`p17_euler`
让"加速度不认识形状"不只是一次巧合。`accelerates ( 0 0 0 0 0 )`——界面算子不是本例的变量（§5.33 交叉过
du/dt × 准 Newton，§5.37 交叉过形状 × 准 Newton），惯量归属也固定在管网侧。

**内核一行没改**，本例读到的只有既成语义：`src/coupling/src/coupledNetwork.cpp:1182` 的
`chord = ex.unsteady ? ex.dpFrictionPa : ex.dpPa` 与 `:1194` 的 `residual = |R*q - chord|/chord`——界面停机
判的**只是摩擦那半区**，`rowDrop = R*q + M*dq/dt` 与 `dpBlock` 的差（`rowError`）从来不被迭代。这个增量做的
是把这句话变成可断言的数字。

十一族断言（`tests/check_transient_coupled_multi_turbulent_profile_unsteady.py`，490 条，RC=0）：

* **加速度是形状盲的，摩擦半区不是**：斜坡真正在动的那层（层 1，X2 的 `dq` = 6e-5），实测惯性对
  `rho*V_probe/A^2*dq/dt` 的比在 X1 是三形状 **0.9623 / 0.9658 / 0.9630**（塞状/抛物/1/7，散布 **0.34 个
  百分点**）、X2 是 0.9792 .. 0.9794（**0.07 pp**）；而同样这三个形状把**摩擦半区**移动 X1 **−14.8%**（层 1）
  / **−51.2%**（层 8）、X2 **−53.2%**。斜坡进入自尾段后（层 3-8）X1 三臂全在探针 slab 预测的 **2.31%** 内，
  且抛物在**每一个**这种层上都最接近 1（L3 1.0084 对 1.0231 … L8 1.0004 对 1.0049）；末解的
  `inertanceDynamic/inertanceProbe` = 1.00491 / 1.00036 / 1.00443（0.6% 内，§5.33 的准稳读法允许 2.8%）。
  形状值这条边压降的一半以上，却值它惯量的零。
* **界面照常停机，因为它停的那半边不含惯量**：80 个状态的 `residual` 全部 <1e-3，而 `residual == |R*q −
  chord|/chord` 这个恒等式在同一批状态上到 **0.0e+00**——同处它们的惯量那半区最多 disagreement **62.7%**。
  循环收敛在动量的一半上，从不看另一半。
* **行侧的惯量在块也开始付费后被记了两次**：`rowError` 层 1 = X1 **53.9% / 62.7% / 54.7%**、X2 **44.1% /
  56.0% / 45.3%**（lo/parab/p17）；末层 X1 仍差 **12.3%（抛物）/ 1.4%（1/7）/ 1.0%（塞状）**，而被钉死的串联边
  降到 **0.0051%**——把它翻倍的那一项已经停止运动了。形状只在它本身很大时才改变双重计费的多少：抛物 62.7%
  对塞状 53.9%，而 1/7 与塞状差不到一个百分点，因为被重复计的是行的 `M*dq/dt`，形状只移动这一层穿过它的速率。
  准稳两臂的 `rowError` 恒 **0**，且 `rowDrop − dpBlock == M*dq/dt − (dpBlock − dpNetwork)` 到 **2.3e-14 Pa**。
* **行误差恰好就是归属那笔算术，不含别的**：恒等式 `rowError*dpBlock == (R*q − dpFriction) + (M*dq/dt −
  dpInertia)` 在 48 个步进状态上到 **2.3e-13 Pa**。于是形状的还压**不在**行误差里：`parab_euler` 的并联边末层
  有 **152.22 Pa** 还压、占它自己 177.43 Pa 压降的 **85.8%**，而行与块只差 −0.11 + 21.98 = **21.87 Pa**。割线把
  还压吞进它谈判的摩擦半区，行从未看见这个形状。
* **符号律与平台第一次是时间的函数**：`beta_out` 在准稳两臂的 32 个状态**全部**落在 (1.02, 1.0211)
  （1.02026815 .. 1.02105374），而且每臂每边跨 8 层的散布 ≤ **8.0e-5**——那是一个数，不是一条轨迹（伙伴自己没有
  时间层就永远不会落后）。步进三臂的**第一层**在带之下（1.01427360 .. 1.01750843），即平台间隙的 **68.4 / 71.0 /
  68.6%**（X1）与 **83.2%**（X2）；但这个未完成的状态与被交的是哪个入口无关：X1 层 1 三臂是 1.014351 /
  1.014274 / 1.014325，彼此差 <1e-4，而强加的 β 相差 **0.18**——块回答 `du/dt` 时，是**层**拥有发展，不是形状。
  到末层步进区域落回与准稳孪生**同一个**平台（1.02104568 对 1.02105374，不高于）：伙伴的表述决定**何时**到达
  终态，从不决定**哪个**终态。入口比平台满的恰好 **32** 个（全部抛物供料）状态、比它差 **48** 个，且与
  `pRisePa > 0` 的那批**同一批**——符号律被问了 80 次，每层两次、每种表述一次。
* **还压的时间签名**：`parab_euler` X1 的 `pRise / (rho*(beta_in−beta_out)*U^2)` 从层 1 的 **0.0704** 爬到末层
  的 **0.9490**，是一条曲线；同一把尺子在准稳孪生上读 **0.9782 .. 0.9928**，正是 §5.36 那条带——差的那一截是
  斜坡，不是律。峰站随区域发展下移：**0.015 m**（层 1）→ **0.035 m**（层 2-3）→ **0.045 m**（层 4 起），准稳臂
  八层恒在 0.045 m。串联边从层 2 起坐在一个比值上（1.1735 .. 1.1800，散布 <0.01）；它的层 1（**0.4811**）是这一族
  唯一的例外，因为一块刚开始回答 `du/dt` 的区域还没发展出形状所对照的那个平台。
* **τ 被 `du/dt` 加长，而形状用同一个乘子**：`lo` 0.4799 → **0.5849 s**（**+21.9%**）、`parab` 0.9786 →
  **1.1683 s**（**+19.4%**）；两种表述都是几乎同样的加价，而形状本身在步进侧 ×**1.9973**、准稳侧 ×**2.0391**。
  1/7 只动 +6.40%（0.6224 s），与它离平台只差 0.0071 的 β 相称。五臂尾段 `|dq|` 衰减都是几何的（0.65920 ..
  0.84267，层内散布 <5.6%），所以一个时间常数能描述整条序列。教材式 `(M_X1+M_By)/(R_X1+R_By)` 是本例最干净的
  一格发现：对准稳两臂它高估 **+22.6% / +27.5%**，对步进三臂它到 **−0.59% / −2.93% / +0.75%**——行侧那笔被
  重复计费的惯量，正是它在准稳侧缺的那一笔。
* **第一层的分流与末层的去处**：`M_By/(M_X1+M_By)` = **0.427275**。准稳两臂仍在它附近（0.419889 / 0.464563，
  9% 内），步进三臂系统性**低于**它 **−18.92% / −14.01% / −18.40%**：在 `dt` = 0.2 s（短于这些臂自己的 τ）时，
  区域惯量不是第一脚唯一要推动的质量，旁路拿走的比行侧惯量声明肯给它的更多。斜坡弛豫到的份额在每种表述内部
  保持形状次序（准稳 0.585000 > 0.407092；步进 0.537958 > 0.412185 > 0.400581），而每个步进臂到末层都还落在
  自己的准稳孪生之后（0.537958 < 0.585000、0.400581 < 0.407092）——1.6 s 大约是三个 τ，不是斜坡的尽头。层 1 由
  伙伴表述移动的量大于由形状移动的量：`du/dt` 使塞状 **−17.5%**、抛物 **−20.9%**，而抛物在一种表述之内只值
  **+10.6% / +6.1%**——双重计费已经在掌舵，形状还没有。
* **被钉死的边是干净实验室**：串联边的速率在全部 40 个（臂 × 层）状态等于 6.0e-5（≤1e-11），层 2 起
  `|dq|` ≤ **1.36e-20**；然而它照报 **8.27 / 8.02 / 8.21 Pa** 的 `dpInertia`——质量加权速度动过，是因为**剖面**
  动了。它的动态惯量（`q/0`）因此不作任何判据，只记录为 5.5e13 / 6.4e16 / 1.4e14。同一个开关在串联边的压降上
  只留 ≤ **0.168%**（层 3+），在并联边的速率上留 **−1.60%（塞状）/ −8.04%（抛物）**，后者是前者的 **5.0 倍**：
  速率自由的地方，双重计费不是可忽略项，它就是轨迹本身。X2 末层压降在两种表述下都保持抛物 < 1/7 < 塞状。
* **时间表、成本与那堵诚实的墙**：五臂 **16 / 18 / 27 / 37 / 31** 趟，SIMPLE 解次数恰为 `2 ×` 趟数
  （32/36/54/74/62），最坏残差 **7.35e-4 ~ 9.58e-4** 全在 1e-3 内，管网线性代数每层到机器精度（newton ≤9.1e-12）
  ——迭代的是界面，不是线性代数。`du/dt` 不免费：×**1.69**（lo）、×**2.06**（parab），因为界面谈判的那个系数
  现在随层一起漂。但**步进不花块内扫数**：平/1-7 的两条 Euler 并联边最差 **17** 扫。300 扫花满与
  `blockConverged 0` 仍是同一批状态（抛物串联边 3/8 准稳、1/8 步进，其余 0/8），同一被钉死速率上塞状/1-7
  **7 扫**停机。
* **伙伴仍然合法，网络仍然闭合**：`y+` 首格 **21.0 ~ 32.6**、`nutMax/nu` **46 ~ 134×**、`Re_Dh` **4807 ~
  10000**、每块进/出到 **≤5.2e-12**（形状在截面上重新分配，步进它也不造质量）。发展窗梯度对模型自己的壁面剪力
  `2 rho u_tau^2/H` 闭合：不带该层加速度时 **−5.67% ~ +18.13%**；把 `rho*du/dt` 记账后抛物两臂朝零走
  （+10.53% → +6.07%、+18.13% → +9.63%，对应 **+104.1 / +169.1 Pa/m**，因为它们并联边的速率到末层还在涨），
  塞状/1-7 的步进臂本来就 <4%（−3.79%、+0.62%）——这项是真的、是剪力的几个百分点，既不能拿来救闭合也不需要
  藏起来。先导段 Churchill 压降与 `Src − J1`、`J3 − Snk` 到 **1e-9**（11.264128 Pa），每层节点对逐边等于那一行
  动量式（含惯量项），`Snk` 恒等于参考压力。

**成本**：五臂、两块 RANS 网格、8 个时间层整案 **70.6 s**（§5.37 同样五臂八层的 `dt 0.05` 是 64.5 s；多出来的
是步进臂的交换次数：16 → 27、18 → 37 趟）。

**本增量的边界**：① `ddts` 是**逐臂**的，不是逐边的：还没有"X1 步进、X2 准稳"这一格，所以"同一层里一半用真
`du/dt`、一半用割线"仍未测——这是本例唯一没能交叉的那根轴，也是下一步的第一格。② 行的惯量在块也开始付费后
**没有被摘掉**：这正是本例测到的 44-63% 双重计费。"`inertias ( partner )` 归属在步进块上是否应当把行侧那一项
挪走"、以及"界面应当在几层内把 `rowError` 压到多少"目前都只有实测，没有做成可关闭判据；归属这根轴本例固定在
管网侧。③ 1D 那一行动量方程仍不认识 `beta`（§5.34①/§5.35①/§5.36②/§5.37② 同一格），而理由更强了：还压现在
是一条时间曲线（4.68 → 152.22 Pa），它那 85.8% 的最大一笔仍然只体现在块被谈判的那半区里。④ `dt` 墙是**记录**：
0.05 那一档下 X1 区域被要求比一层前更慢、于是真实地报出压强上升，撞上"区域以正压降回答"的内核护栏——四种设法
留下 −22.75 / −541.94 / −50.12 / −128.13 Pa 四个数，割线欠松到 0.4 的那档活下来但 12 趟后停在 1.7035e-3 对
1e-3。本例把 `dt` 提到 0.2 s 使整案闭合，但那几档仍未升格为判据。⑤ 300 扫花满还是"记录 + 与 `blockConverged
0` 同批"的一致性断言。⑥ 两条耦合边一条并联一条串联；形状 × 双边并联（§5.32 的拓扑）未与本例交叉。⑦ 1.2e-4
那堵发散墙仍是记录（§5.36① / §5.37⑥ 同一格）。⑧ 界面条件里惯量那半区被**度量**而从未被**迭代**——本例把这件
事变成了可断言的数，但没有改它。

### 5.39 逐边 `ddt` 表：同一层里一半真 `du/dt`、一半割线——你没动的那条边逐字就是它自己控制臂的算数，混合改变的只有界面停哪半区、双重计费住在哪、以及谈判落到的价（§5.38 边界① × §5.38 的算例，含 app 声明面改动，零内核改动）

§5.38 自己写下边界①：`ddts` 是**逐臂**一张表，所以一次运行里所有耦合边要么一起步进、要么一起准稳，"同一层
一半用真 `du/dt`、一半用割线"这一格从未被测。这一格能问出一根轴各自问不出的问题，而且问的对象是内核早就
存在的一段代码：`src/coupling/src/coupledNetwork.cpp:1182` 的 `chord` 是**逐块**取的
（`ex.unsteady ? ex.dpFrictionPa : ex.dpPa`），所以一个混合臂在一趟交换里谈判的是**两种不同的表述**——一条边
的电阻被摩擦半区计费，另一条的电阻被它的整压降计费。于是要问的是：**分流是否等于它各条边的和，还是这一层会长
出全定常与全 Euler 都没有的性格**。

算例 `tests/cases/transientCoupledMultiTurbulentEdgeDdt/case.dict`：拓扑、两块网格、壁面函数、操作点、`dt` 与
8 个时间层全部照抄 §5.38（`Src → A → { X1 ‖ By } → X2 → B → Snk`，6.0e-5 m3/s 在 t=0 接入，X1 20×8×2 over
0.2×0.005×0.005 / probe 0.12–0.19，X2 24×8×2 over 0.24×0.004×0.008 / probe 0.16–0.23，proxy 孔径取各自截面的
等效直径，所以每条边的行侧 `M` 就是那条边的块在自己的网格上量到的 `ρV/A²`），唯一新增的是**表述可以逐边
声明**：`couple{ ddtFrom { lo_x1free { X1 Euler; X2 steady; } lo_x2pin { X1 steady; X2 Euler; } parab_x1free
{ X1 Euler; X2 steady; } } }`。优先级是 `ddtFrom` > `ddts` > 块 dict 自己的 `schemes{ ddt }`，而且两个覆盖方向
都走一遍：三臂的 per-arm 词是 `Euler`，表把其中**一条**边收回 `steady`（表压过逐臂词，也压过块自己的声明）；
`lo_steady` 的 per-arm 词是 `steady`，它压过两块 dict 都写着的 `Euler`。`accelerates ( 0 0 0 0 0 0 0 )`、惯量归属
固定在管网侧、`tolerance 1e-3` / `maxSubIterations 12`，所以臂与臂之间只剩"哪条边被表点了名"。

七臂是两组"镜像 + 各自的两个对照"：`lo_steady`（S/S，§5.38 的控制）、`lo_both`（E/E，§5.38 的交叉，这里由
逐臂词而非表达到）、`lo_x1free`（E/S：速率**弛豫**的那条边步进）、`lo_x2pin`（S/E：镜像，被拓扑钉死的那条边
步进），以及把这三行按抛物入口重做一遍的 `parab_steady` / `parab_both` / `parab_x1free`。塑形族存在的理由很
具体：混合臂的每条边都要有一个**同时固定该边表述与该边形状**的对照，否则"这条边没变"只能和一条平入口的边
比，差异无法归给邻居。X1 = 并联/"自由"边（它的速率是分流在弛豫的量），X2 = 串联/"钉死"边（它的速率被拓扑钉在
注入量上，层 2 起 `|dq|` ≤ 1.36e-20）。

**内核一行没改**。改的是 app 的声明面（`src/coupling/apps/transientCoupledMultiApp.cpp:272`（`ddtFrom` 逐臂逐边
解析，缺边或非法词直接 `cfFatal`）、`:361`（`schemes{ ddt }` → `ddts[臂]` → `ddtFrom[臂][边]` 三级覆盖）、`:440`
（打印**生效的**逐边词，而不是声明的那个）），本例读到的接口语义仍全是既成的。

十族断言（`tests/check_transient_coupled_multi_turbulent_edge_ddt.py`，202 条，RC=0）：

* **声明到达"边"，且生效词与块自己的回报一致**：112 个（层 × 边）状态各带**自己那条边**的标志——56 个说区域
  在该层答了 `du/dt`、56 个说没有，且在两个混合臂里这两种同时坐在同一层里；三个混合臂的两个词确实不同。块的
  回报逐边一致并交回它解过的那一层：步进边 `dtLevel` = 0.2 s、准稳边 = 0。每条边的行惯量是该边网格等效孔径的
  `ρL/A`（X1 **8e6**、X2 **7.5e6**，到 1e-14）且块在自己的网格上量到同一个数；每条边计费加速度的 slab 也就是
  它开场宣布的那一块（X1 2.8e6、X2 2.1875e6）。
* **半区是逐边切的，不是整臂混的**：`dpFriction + dpInertia == dpBlock` 在 112 个状态到 **2.2e-16**；56 个准稳
  状态的惯性半区**恰为 0**、摩擦半区就是整压降、`rowError` 恰为 0；56 个步进状态全带非零惯性半区，**包括速率被
  钉死的那些**；混合臂的八个层里每层恰有一条边有惯性半区——这个臂是逐层混合的，不是被时间平均成一种表述。而行
  那笔 `ρL/A·dq/dt` 在**两条**边上都留在账上（节点对 == 摩擦 + 行惯量项，112 个状态到 **0.0 Pa**），不管那条边
  的块是否自己在付费。
* **界面停在哪半边也是逐边读的**：按边重算 `|R·q − chord|/chord`（步进边取摩擦半区、准稳边取整压降）复现打印
  残差到 **0.0e+00**，且 112 个状态全部落在 1e-3 内——混合的一层仍然是收敛的一层。这个选择不是装饰：层 1 的 7
  条步进边按摩擦半区读 **0.0000 ~ 0.0137%**，同一批状态按整压降读 **17.62 ~ 32.00%**（高两到三个数量级）；而每
  条准稳边这两种读法是**同一个数**（最差 0.0e+00），因为在那里摩擦半区就是整压降。
* **本例最强的一格——你留下的那条边就是它自己控制臂的算数**：混合臂在 X1 上与"同表述同形状"的对照**逐位相同**
  （15 个字段 × 8 层全部差 **0.0**：速率、压降两半区、系数、残差、扫数、行误差、β、还压、标志）。它在 X2 上
  复现自己那条边的控制：平入口 **6.210e-4 Pa** of 1538.4（4e-7 相对）、塑形 **4.476e-1 Pa** of 719.8（6.2e-4
  相对），速率到 6.8e-21，行误差 0.0 / 2.3e-5。镜像臂同样成立：`lo_x2pin` 的步进钉死边距全 Euler 臂的同一条边
  **8.6e-4 Pa**，且每一层扫数相同——所以相等是对"该边自己的表述"成立，不是对某个特定的臂成立。作为尺度：同
  一条边把**自己的**表述翻过来会动 **662.0 / 2109.5 / 799.7 / 2148.9 Pa**（平 X1 / 平 X2 / 塑 X1 / 塑 X2），比
  一个混合邻居在它所没点名的边上留下的痕迹高三到六个数量级（最坏 6.2e-4 相对）。
* **两样东西不成加和，且都落在"价"与"预算"上**：一是界面最后留在被留下的那条边上的系数——它移动 **27.2 /
  588.8 / 10138.7 Pa·s/m³**（相对 1.1e-6 / 2.3e-5 / 8.5e-4）而该边的压降只动 4e-7 ~ 6.2e-4：**能相加的是那条边
  的物理，不是谈判落到的那个价**。二是塑形钉死块的预算：`parab_x1free` 与 `parab_steady` 的 X2 差 **290** 扫，
  混合臂有 6/8 层花满 300 扫、它的控制只有 3/8——这是 §5.37 那堵容差地板透出来：速率被管网握着时，塑形入口是
  **解到接近预算**的，所以某一层够不够预算是邻居的轨迹可以翻的一枚硬币，而它交回的压降仍按 6.2e-4 相对一致。
* **双重计费在步进的那条边上出生，也止于那条边**：层 1 的 `rowError` 是逐边的——X1 **53.89%**（平）/ **62.66%**
  （塑）、X2 **44.06%** / **56.01%**；同一批混合臂那 8 层里被留准稳的那条边报**恰 0.0**（双重计费不越过它共享的
  节点对）。末层自由边仍差 **1.021%**（平）/ **12.328%**（塑），被钉死的步进边降到 **0.000005 ~ 0.005102%**——
  把它翻倍的那一项已经停止运动。而步进钉死边仍是干净实验室：层 2 起 49 个状态的速率都是 6.0e-5（`|dq|` ≤
  1.36e-20），却照报 **8.271 / 8.271 / 8.020 Pa** 的 `dpInertia`——层 1 从静止跳出去留下一笔块至今还带着的账，
  而塑形的块付的是另一个数，因为它的场是另一个场。
* **斜坡属于速率自由的那条边**：混合臂层 1 的分流与"自由边同表述同形状"的对照**逐位相同**（`lo_steady` ==
  `lo_x2pin` 都是 0.419889、`lo_both` == `lo_x1free` 都是 0.346431、`parab_both` == `parab_x1free` 都是
  0.367420）——表没点名的那条边把分流改变了 **0.0**。末层的去处同理（0.407092 / 0.400581 / 0.537958），尾段衰减
  比也是同一个数而不是相近的数（0.64159 / 0.71416 / 0.84483），于是 τ 也是同一个数（自由边层 3→8 的几何衰减）：
  **0.4506 / 0.5941 / 0.9897 / 1.1861 s**。逐边看：`du/dt` 加在**弛豫**那条边上使斜坡长 **+31.8%**（平）/ **+19.8%**
  （塑），加在**钉死**那条边上使斜坡加长 **0**（那两臂的 τ 逐位相等）；混合臂内部形状仍按 ×**2.00** 乘它（1.1861 s
  对 0.5941 s）——§5.37/§5.38 量的那根杠杆，现在第二条伙伴从未离开割线的臂上也成立。交换时间表也逐层相同
  （[3,2,2,2,2,2,1,2] / [4,7,5,3,3,2,1,2] / [3,7,5,5,4,4,5,4]），七臂总趟数 **16 / 27 / 27 / 16 / 18 / 37 / 37**：
  `du/dt` 放在弛豫边上花 **×1.69**（平）/ **×2.06**（塑），放在钉死边上花 **×1.00**。SIMPLE 解次数恰为 `2 ×` 趟数
  （32/54/54/32/36/74/74，共 356 次）——时间表是界面的工作，不是块的工作。对教材式 `M_By/(M_X1+M_By)` =
  **0.427275**：准稳自由边 −1.73%（平）/ +8.73%（塑），步进自由边 −18.92 / −18.92 / −14.01 / −14.01%。
* **行被记的那笔加速度仍然不认识形状**：混合臂层 1 自由边的实测惯性对 `ρV_probe/A²·dq/dt` 的比是塑形
  **0.965755**、塞状 **0.962319**（差 **0.34 个百分点**），且每种形状内部与 §5.38 的同形状臂同一个数——把第二个
  伙伴从时间层里拿出来，并没有教会第一个认识形状；而**同一个开关在同一条边上**把摩擦半区移动 **−14.8%**（层 1）/
  **−51.2%**（末层）。交进去的 β 仍是 8 面板离散声明值（到 1.1e-15）、入口条件的截面积分仍是该层的速率（到
  6.7e-16），步进边与准稳边一样。
* **平台与还压也是逐边的事实**：准稳的边从层 1 起就**坐在**平台上（1.02028287 ~ 1.02105374，跨 8 层散布
  **3.0e-9 ~ 8.0e-5**，是一个数不是一条轨迹），步进的边从带下出发再落到带上（1.01427360 / 1.01435119 / 1.01728347
  / 1.01750844 → 1.02009350 / 1.02098341 / 1.02104568）——区域落后自己的伙伴，差的量是**层**，不是形状。符号律
  同样逐边逐层：报出压强上升的恰好是 **48** 个（112）状态，即三个塑形臂的两条边 × 八层——两种表述、两个位置都
  算——且与"交进去的 β 高于区域读回的 β"那批是同一批。同一个塑形臂的两条边只差"自己的 `du/dt` 在不在产生压降
  的那条方程里"：准稳的钉死边八层恒在**一站**（333.4297 Pa @ 0.035 m，层间动 **1.2e-6 Pa**，速率被钉住且形状固定
  就没有斜坡可签名），步进的自由边是一条曲线（**4.6838 → 152.2237 Pa**，峰站随区域发展走
  **0.015 → 0.035 → 0.045 m**）。这条边的镜像是 `parab_both` 的钉死边：当它**也**步进时，峰站走 [0.025, 0.035]、
  摆 **200.9 Pa**——还压峰在哪一站是这个边自己的表述在说话，不是形状、也不是位置。可分离性也成立：混合臂与它
  的塑形控制在准稳钉死边的还压上差 **1.0e-6 Pa**，对照步进边上那 **147.5 Pa** 的曲线。
* **块内预算、那堵墙，以及墙与失败第一次不完全重合**：逐扫也是逐边的——两条把它留着定常的臂里准稳钉死边每层
  都 **7** 扫，两条步进它的臂也各 **7** 扫，两条步进塑形自由边的臂序列都是 [17,19,20,21,20,21,23,22]。墙属于形状
  且在每种表述下都在：`parab_steady` **3/8**、`parab_both` **1/8**、`parab_x1free` **6/8** 层花满 300 扫，而每一条
  平臂与每一个 X1 块从不碰它（其余十一列最差 **28** 扫，平的混合臂不超过 **17**）。新东西在这里：`parab_x1free`
  那 6 层里有 **5** 层报 `blockConverged 0`——花满预算与没收敛**不是同一事件**，有一层到了预算仍然满足自己的容差。
  全部 112 个状态最坏残差 **9.578e-4**（在 1e-3 内），所以花满是块预算的事实，不是耦合失败；管网线性代数每层到
  **9.09e-12**——迭代的是界面，从来不是矩阵。
* **伙伴仍然合法，网络仍然闭合**：首格 `y+` 23.0 ~ 32.6（两种表述、两条边、七个臂都在对数律带里）、`nut/ν`
  **46 ~ 134×**、末层 `Re_Dh` **4807 ~ 10000**、每块进/出 ≤ **3.8e-12**（形状与时间层都不造质量）。发展窗梯度对
  模型自己的壁面剪力 `2ρu_τ²/H` 在不带该层加速度时 **−5.67% ~ +18.13%**；把该层加速度记账后有三种读法，而且它们
  分别说出这笔账归谁：钉死边**是同一个数**（−5.66 / −5.67 / +7.73 / +7.78 / +7.79%，速率停了就没得加，块带不带
  `du/dt` 都一样）；混合臂里**被留着定常的自由边仍然移动**（−5.51 → −5.29%、+10.53 → +6.07%），因为那项是**行**
  的、来自管网自己的速率变化，块的表述不拥有它；步进的自由边按形状分方向（平 −3.79 → −5.41% 离零，+28.5 Pa/m；
  塑 +18.13 → +9.63% 朝零，+169.1 Pa/m）——修正的符号是区域自己的历史，不是表述的性质。先导段 Churchill 压降与
  `Src − J1`、`J3 − Snk` 到 1e-9（**11.264128 Pa**），每层节点对逐边等于那一行动量式（含惯量项），`Snk` 恒等于参考
  压力，`A` 恒等于注入量（<1e-18）且 X1 + By 分流到 1e-11。

**成本**：七臂、两块 RANS 网格、8 个时间层、178 趟界面 / 356 次 SIMPLE 解；单独跑整案 **129 s**，在 `-j8` 的
ctest 里两次分别记 **137.6 s** 与 **214.4 s**（墙钟由同批负载决定）。§5.38 的五臂是 70.6 s；多出来的是两臂 ×
两形状的四条步进并联边。

**本增量的边界**：① 逐边 `ddt` 与逐边惯量归属（`couple.inertiaFrom`，§5.29 的机制）**没有同时**出现——本例把归属
固定在管网侧，所以"这条边自己测的惯量 + 那条边步进"仍是一格空白。② 形状侧只写了三臂：塑形臂只有"X1 步进 /
X2 准稳"这一种混合（`parab_x2pin` 未做），因此塑形混合行的镜像仍未测；塞状侧的镜像是本例可分离性的主要支柱。
③ 逐边 `ddt` × 两条边都在并联位置（§5.32 的拓扑）未交叉，"worst 换手在两种表述之间怎么发生"还没被问过。④ 界面
最终落在哪个价（27.2 ~ 10138.7 Pa·s/m³）本例**测到而不解释**：它不进任何守恒律，只是谈判的终点。⑤ 44 ~ 63% 的
双重计费仍然只被**度量**：行侧那项 `M·dq/dt` 在块也开始付费后没有被摘掉，"几层内把 `rowError` 压到多少"与
`inertias ( partner )` 是否应当在步进块上摘掉行侧那一项都仍是实测而非判据。⑥ 1D 那一行动量方程仍不认识 `beta`
（§5.34①/§5.35①/§5.36②/§5.37②/§5.38③ 同一格）。⑦ 300 扫花满仍是记录 + 一致性断言，而本例第一次出现"花满与
`blockConverged 0` 不重合"（6 对 5），这个不重合本身没有升格为判据。⑧ `dt 0.05` 那堵墙与 1.2e-4 那堵发散墙仍是
记录（§5.36① / §5.37⑥ / §5.38④⑦ 同一格）。⑨ 界面上的伙伴仍是单相不可压块；多相/可压伙伴未进接口。

### 5.40 逐边 `ddtFrom` × 逐边 `inertiaFrom` 同台：行侧那笔双重计费第一次被**分解**成两项——一项是界面永不迭代的弦差，一项是"这条边的惯量记在谁账上"的账单，而 `slab` 只消掉后者（§5.39 边界① × §5.29 的归属表，零内核改动）

§5.39 自己写下边界①：逐边 `ddt` 与逐边惯量归属（`couple.inertiaFrom`，§5.29 的机制）**没有同时**出现。这一格
不只是"少跑了一个组合"：§5.39 还把归属固定在管网侧，并且把两根 proxy 孔径写成各自截面的**等效直径**，于是
`ρL/A_proxy`（图纸）与 `ρV/A²`（区域）本来就是**同一个数**——归属那张表在任何结果里都看不见。要把"记了谁的惯量"
变成可测量的东西，得先故意把孔径写歪；要问出它**为什么**可测量，得先知道内核那一行到底在数什么：
`src/coupling/src/coupledNetwork.cpp:1160` 起的每一行是 `rowDrop = R·q + M·dq/dt`，而评分是
`rowError = (rowDrop − dpPa)/dpPa`（`:1201`，准稳边为 0，`:1202`）。于是在一条步进边上恒有

```
rowError · dpBlock = (R·q − dpFriction) + (M·dq/dt − dpInertia)
                       \________ 弦差 _______/   \____ 归属账单 ____/
```

第一项是界面**从不迭代**的那半（`chord = unsteady ? dpFrictionPa : dpPa`，`:1182`），第二项才是"这一行的 `M` 从
哪本账来"。本例要问的是：**这个分解是否逐边成立（两张表在同一层里交叉时也一样），而 `slab`——那个数值正好就是
块量自己压降所跨的那块流体的来源——是否把第二项消掉，只剩第一项。**

算例 `tests/cases/transientCoupledMultiTurbulentEdgeDdtInertia/case.dict`：拓扑、两块 `kOmegaSST` 网格、probe、
壁面函数、操作点、`dt 0.2` 与 8 个时间层全部照抄 §5.39（`Src → A → { X1 ‖ By } → X2 → B → Snk`，6.0e-5 m³/s 在
t=0 接入，X1 20×8×2 over 0.2×0.005×0.005 / probe 0.12–0.19，X2 24×8×2 over 0.24×0.004×0.008 / probe 0.16–0.23）。
**唯一改动是每根耦合管各自的图纸孔径**：X1 写 0.005（等效 0.0056418958354775631）、X2 写 0.006（等效
0.00638307648642292245）。于是三本账真的成为三个数——X1 图纸 **1.0186e7** / 区域 **8.0e6** / slab **2.8e6**，X2
图纸 **8.4883e6** / 区域 **7.5e6** / slab **2.1875e6** kg/m⁴；图纸把区域高估 **1.273×**（X1）与 **1.132×**（X2），
slab 是区域的**长度分数** 0.350000（0.07 m of 0.20）与 0.291667（0.07 of 0.24）。八臂：

* `nn_steady` / `pp_steady` —— 两条边都准稳，图纸 / 区域。这是 §5.39 无法评分的那一格（准稳边的 `rowError` 按定义
  为 0），用来问"评分看不见时，归属改不改轨迹"。
* `nn_both` / `pp_both` / `ss_both` —— 两条边都步进，图纸 / 区域 / slab：纯归属阶梯。
* `es_nnn` —— X1 步进、X2 准稳，归属处处图纸：即 §5.39 的混合臂。
* `es_mix` / `se_mix` —— 同一格混合，但**步进那条边改吃 slab**，互为镜像。这两臂是本例的新格：在**同一条边**上
  `ddtFrom` 与 `inertiaFrom` 各自点名，而逐臂开关说不出"X1 是 Euler+slab 且 X2 是 steady+network"这种话。

`ddtFrom` 逐臂逐边解析在 `src/coupling/apps/transientCoupledMultiApp.cpp:273`，`inertiaFrom` 在 `:303`（非法词、
缺边直接 `cfFatal`），两份表在 `:361` / `:395`–`:409` 各写成逐边向量，归属在安装时一次决定、**在时间循环之前**
（`coupledNetwork.cpp:1092`），所以它对准稳边同样生效；`:420`/`:432` 打印每臂每条边**生效的**来源，`:441` 打印块
自己宣布的 `regionInertance` 与 `probeInertance`。**内核一行没改**：两张表都是既成机制，本例只是把它们放到同一
条边上，并把让它们可分辨的那个数（孔径）写歪。

十族断言（`tests/check_transient_coupled_multi_turbulent_edge_ddt_inertia.py`，203 条，RC=0）：

* **两张表到达同一条边，且块能核对接口装上去的数**：128 个（层 × 边）状态各带**自己**那条边的两个词——72 个步进
  / 56 个准稳，混合臂里两种坐在同一层；`es_mix` 与 `se_mix` 的两份表在**同一条边**上互相不同意。装上行的 `M` 与该
  来源宣布的数到 **1e-14** 一致，而三个来源给出三个不同的数；块的回报里同时带着 `inertance`（区域）与
  `inertanceProbe`（slab），所以接口安装的是伙伴可以自己核对的数，不是它猜的。
* **本例的主判据：双重计费确实是两项之和，逐边成立**——`rowError·dpBlock == (R·q − dpFriction) + (M·dq/dt −
  dpInertia)` 在 72 个步进状态全部成立，最差 **6.4e-17**（按它由那四项的绝对值之和归一；`es_mix` 层 1 X1）。半区
  本身仍是逐块切的：`dpFriction + dpInertia == dpBlock` 到 **2.2e-16**，56 个准稳状态的惯性半区与 `rowError` **恰为
  0**。层 1 的两项根本不是一个量级：X1 弦差 **+0.0472 / +0.1015 / +0.6006 Pa** 对归属账单 **+718.87 / +551.40 /
  +11.79 Pa**（15236× / 5430× / **20×**），X2 弦差 ≤**2.5e-4 Pa** 对 **+1903.86 / +1607.38 / +13.63 Pa**。到平台
  顺序翻转：`ss_both` 层 8 X1 是弦差 −0.3224 Pa 对归属 **+0.0066 Pa**——能被来源移掉的只有第二项。
* **阶梯单调，且第二项就是图纸多收的那笔**：层 1 `rowError` 在 X1 上 **77.0836% → 53.8941% → 0.9364%**（图纸 /
  区域 / slab），在 X2 上 **52.1912% → 44.0637% → 0.3737%**；`slab` 消掉图纸账单的 **98.36%**（X1）/ **99.28%**
  （X2）。而 `(M_图纸 − M_slab)·dq/dt` 预测出 **708.19 Pa**（X1）与 **1890.23 Pa**（X2），对实测 **+718.87 /
  +1903.86**——差额就是块自己那项 `dpInertia` 的量级，不是别的东西。钉死边的阶梯在**绝对价**上更高（1903.9 对
  718.9、1607.4 对 551.4、13.6 对 11.8），因为它的 `dq/dt` 就是整步。
* **归属只挪行侧那半，区域那半几乎不动**：层 1 X2 的 `dpInertia` 在三种归属下都是 **642.6188 Pa**，最差 **1.1e-8**
  相对；更强的是**两次归属之差恰是 `(M_a − M_b)·dq/dt`**——**296.4791** 对 296.4791（4.3e-8 相对）、**1593.7496**
  对 1593.7500（2.7e-7）、**1890.2286** 对 1890.2291（2.4e-7），而每一档那点残差（1.3e-5 / 4.3e-4 / 4.5e-4 Pa）
  恰好被块自己那两项之差（弦差 + `dpInertia`，同量级 1.3e-5 / 4.5e-4 / 4.6e-4 Pa）吸收。所以 `rowError` 对装上的
  `M` 是**线性的，斜率就是该层的 `dq/dt`**。反面：块的惯性半区对 `M_slab·dq/dt` 的比在 9 个层 1 步进状态里是
  **96.02 ~ 97.92%**——slab 是**对的板、错的 `du/dt`**，这正是它留下的那一点。
* **消掉之后剩下的就是界面本来在追的那个数**：装了 slab 的 24 个（层 3 起 × 步进）状态里 `|rowError|` 最大
  **9.646e-4**，**低于界面自己停机的 1e-3 容差**；在并联（自由）边上它还落在同层残差的 **1.73×** 以内。对照图纸的
  自由边：同样这些层上 `|rowError|` 到 **1.176e-1**，是容差的 **118×**、自身残差的 **153.2×**——这一项界面既迭代
  不掉、容差也没要求它迭代。
* **准稳边照付，却一分不评**：`nn_steady` 层 1 X1 的行装上 **1165.3 Pa** 的 `M·dq/dt`，而那块区域报 **310.8 Pa**
  且没有加速度——若评分就会是 **375%**，但 `rowError` 打印 **0.0**；整个准稳族（56 个状态）都是几百到上千帕的
  `M·dq/dt` 对上恒零的评分。这不是记账漏洞，而是 `:1202` 那句定义的直接后果。
* **两条局部性定律，且不是同一条定律**：`es_nnn` 的 X1 与 `nn_both` 的 X1 **15 字段 × 8 层逐位 0.0**（邻居从 Euler
  改成 steady，图纸孔径照旧），`se_mix` 的 X1 与 `nn_steady` 的 X1 也**逐位 0.0**（邻居同时改了**两张表**）——表述
  跨边可交换。但 `es_mix` 的 X1 对 `ss_both` 的 X1（该边**同表述同来源**，只有邻居不同）只 **2/15** 字段逐位相同，
  速率差 **7.41e-15**（2.9e-10 相对）、压降差 **4.54e-7 Pa**（3.4e-10 相对）：**表述可交换，归属不可交换**，因为
  装上去的 `M` 是全局线性系统的一项。同一规律在准稳一侧也成立：三条准稳钉死边相对 `nn_steady`/X2 只差 **2.232e-4
  ~ 5.900e-4 Pa**（对 1538.4 Pa，约 4e-7 相对）。而"评分为 0"绝不等于"来源无所谓"：`nn_steady` 与 `pp_steady` 的
  自由边速率差到 **9.18%**、压降差 **50.5 Pa**，两者的 `rowError` 却是**同一个 0.0**；`es_mix`/X1 与 `es_nnn`/X1
  在同一条边上只差来源，速率最差相差 **6.735e-6 m³/s（28.3%）**。钉死边的**压降**沿阶梯不动（4.4e-8 相对），
  它的 `rowError` 却跨 **51.8 个百分点**——账的两侧各走一头。
* **斜坡跟着自由边那两条表走，且阶梯是单调的**：层 1 分流 slab **0.431868** > 区域 **0.346431** > 图纸
  **0.319615**（步进）、区域 **0.419889** > 图纸 **0.381355**（准稳）——行越抵抗加速，第一步越少走它；τ 同序
  **0.2997 / 0.6014 / 0.6134 s**（slab → 图纸 **2.0×**），准稳侧 **0.4329 / 0.4773 s**（那是管网自己的惯性松弛，
  不是伙伴的）。三个混合臂与"自由边同表述同来源"的对照**逐位相同**（`es_nnn` == `nn_both` 0.319615 / 0.6134、
  `es_mix` == `ss_both` 0.431868 / 0.2997、`se_mix` == `nn_steady` 0.381355 / 0.4773）——**表没点名的那条边把分流
  改变 0.0**。去处仍聚在一起（0.396201 / 0.400581 / 0.407166 / 0.404869 / 0.407092）：归属在层 1 分开注入量的
  **11.2%**，到层 8 只剩 **1.1%**，所以平台散布是容差，不是某个来源选定的物理。
* **成本由自由边的两张表定，且单调**：趟数 **16 / 16 / 26 / 27 / 30 / 26 / 30 / 16**（共 **187**），SIMPLE 解次数
  恰为 `2 ×` 趟数（共 **374**），每臂最差扫数 **15 ~ 17**、**没有一个块碰 300 扫预算**、`blockConverged` 全 1、最坏
  残差 **1.861e-5 ~ 9.668e-4** 全在 1e-3 内。准稳 16 趟 → 图纸 26 → 区域 27 → slab 30；而 `se_mix`（步进的那条边是
  **钉死**的）只花 **16** 趟，与一条完全准稳的臂同价——界面 effort 住在自由边。
* **伙伴仍合法、网络仍闭合**：首格 `y+` **22.8 ~ 32.6**（对数律带里）、`nut/ν` **45.2 ~ 74.9×**、`Re_Dh` **5215 ~
  10938**、`β_out` **1.020947 ~ 1.021052**、每块进/出 ≤ **2.6e-16**，且所有边的入口形状都是 flat——下面任何差异都不
  是形状效应。每层节点对逐边等于**那一行的动量式**（`R·q` 加上装上去的 `M·dq/dt`，1e-9），不是块的压降：所以案件
  命名的那本账**进了网络采纳的压力场**；先导段 Churchill **11.264128 Pa** 到 1e-9，`Snk` 恒等于参考压力，`A` 恒等于
  注入量（<1e-18）、X1 + By 分流到 1e-11；层 1 行与区域最大 disagreement **2546 Pa**，两项都有界，没有发散。

**成本**：八臂、两块 RANS 网格、8 个时间层、187 趟界面 / 374 次 SIMPLE 解；单独跑整案约 **28 s**（§5.39 的七臂是
129 s——本例的八臂全部是平入口，且没有任何块碰预算）。

**本增量的边界**：① 准稳边 × `slab` 归属这一格**没跑**——表允许它（归属在时间循环之前安装，与表述无关），但本例
的准稳边一律吃图纸或区域，所以"评分恒零而账单被消掉"那一格仍只有推理、没有实测。② 形状侧完全没交叉：八臂全是
flat 入口，§5.39 已经证明塑形块会把扫数预算翻出来，而"塑形 × 归属"这一格还空着（`parab_x2pin` 也仍缺）。③ 界面
算子固定在 GS（`accelerates ( 0 0 0 0 0 0 0 0 )`），所以"归属改变收敛"在准 Newton 下会怎么表现（§5.33 交叉过
`du/dt` × QN）未测。④ 逐边两表 × 两条边都在并联位置（§5.32 的拓扑）未交叉。⑤ `rowError` 仍然**不是判据**：本例
把它分解、量化、并证明了哪一项能被来源消掉，但没有升格成"几层内压到多少"的停机条件，也没有把 `slab` 变成推荐选
项——它是一次声明，不是推导。⑥ `es_mix`/X1 与 `ss_both`/X1 那 **2.9e-10** 的速率残差只被量化、机制未定（与 §5.31
边界② 同一类观测：归属进了全局解，究竟哪一步在最后一位分岔没查出来）。⑦ 1D 那一行动量方程仍不认识 `beta`
（§5.34①/…/§5.39⑥ 同一格）。⑧ `dt 0.05` 那堵墙、1.2e-4 那堵发散墙、"花满 300 扫 ≠ `blockConverged 0`"仍是记录
（§5.36①/§5.38④⑦/§5.39⑦⑧）。⑨ 界面上的伙伴仍是单相不可压块；多相/可压伙伴未进接口。

### 5.41 塑形入口 × 逐边两表同台：§5.40 分解出的两笔账单里，形状只动**块那一半**——弦差被形状放大近百倍却仍是账单的 1e-5，而账单只跟着 `dq/dt` 走（§5.40 边界② × §5.39 缺的 `parab_x2pin` × §5.40 边界①，零内核改动）

§5.40 把一条步进边的行误差写成两笔——

```
rowError · dpBlock = (R·q − dpFriction) + (M·dq/dt − dpInertia)
                       \____ 弦差 ____/   \____ 归属账单 ____/
```

并留下边界②：八臂全是平入口，"形状 × 归属"这一格还空着，`parab_x2pin`（§5.39 记下的另一格）也仍缺；同时留下边界①：准稳边 × `slab` 归属表允许、但没跑。这三格合起来其实是一句话——**形状到底作用在哪一笔上**。这不是修辞：账单那笔是 `(M_装上 − M_块自己付)·dq/dt`
的算式，所以在**钉死边**（`dq` 由拓扑给定，与轨迹无关）上，如果连形状都能移动它，那么"记了谁的惯量"和"块做了什么"就不可分；如果它不动而块的 `dpInertia`/`dpBlock` 大动，那么归属阶梯就是所装 `M` 的性质，入口形状只是给块那一半重新定价。

算例 `tests/cases/transientCoupledMultiTurbulentProfileInertia/case.dict`：拓扑、两块 `kOmegaSST` 网格、probe、壁面函数、操作点、**故意写歪的 proxy 孔径**、`dt 0.2`、8 个时间层与 1e-3 容差全部照抄 §5.40（`Src → A → { X1 ‖ By } → X2 → B → Snk`，6.0e-5 m³/s 在 t=0 接入；X1 20×8×2 over 0.2×0.005×0.005 / probe 0.12–0.19 / 图纸孔径 0.005 对等效 0.00564190；X2 24×8×2 over 0.24×0.004×0.008 / probe 0.16–0.23 / 0.006 对 0.00638308），**唯一加的东西是 `starts.<arm>.profile`**——逐边的入口形状列表。三本账仍是三个数：X1 图纸 **1.018592e7** / 区域 **8.0e6** / slab **2.8e6**，X2 图纸 **8.488264e6** / 区域 **7.5e6** / slab **2.187500e6** kg/m⁴；图纸高估区域 **1.2732×** / **1.1318×**，slab 是区域的长度分数 **0.350000** / **0.291667**。

七臂：`f_net` / `f_region` / `f_slab` 是平口的纯归属阶梯（就是 §5.40 的 `nn_both` / `pp_both` / `ss_both` 换了名字），`p_net` / `p_region` 是把同一阶梯的前两档塑形（形状单轴的对照），`p_mix`（X1 步进 + 图纸 + 塑形，X2 准稳 + slab）与 `p_pin`（X1 准稳 + 图纸，X2 步进 + slab）把两张逐边表与形状放到同一条边上——`p_pin` 正是 §5.39 缺的 `parab_x2pin`，而 `p_mix`/X2 正是 §5.40 边界① 的准稳 × slab。

**一堵墙，先记下来**：`slab` 归属 × 步进 × **自由**边一旦同时遇上塑形入口，层 1 内就发散——`coupledNetwork.cpp:628 detailBlock: the detailed region reports a non-positive drop (-89.475517 Pa) at 0.000022 m3/s`。所以本例的 `slab` 只出现在准稳边与钉死边上，"塑形 × slab × 自由步进边"仍是推理而不是结果。

声明面：`profile` 的逐边列表在 `src/coupling/apps/transientCoupledMultiApp.cpp:170`–`:193` 解析（长度不等于耦合边数直接 `cfFatal`），每臂每条边**生效的**形状在 `:433` 打印；内核词表 `flat` / `parabolic` / `power1over7` 在 `src/coupling/src/coupledNetwork.cpp:423`–`:426`，写入在 `:491`。`ddtFrom` 在 `:273`、`inertiaFrom` 在 `:303`——两张表和形状列表都是既成机制，**内核一行没改**，本例只是把它们放到一起。

十一族断言（`tests/check_transient_coupled_multi_turbulent_profile_inertia.py`，158 条，RC=0）：

* **两张表 + 形状逐边到达，且都能读回**：112 个（层 × 边）状态各带自己那条边的三个词——96 个步进 / 16 个准稳，`p_mix` 与 `p_pin` 的两份表在**同一条边**上互相不同意，而逐臂开关说不出"X1 是 Euler+图纸+parabolic"这种话；装上的 `M` 与该来源宣布的数逐字一致，块的回报同时带 `regionInertance` 与 `probeInertance`。
* **恒等式对形状免疫**：`rowError·dpBlock == (R·q − dpFriction) + (M·dq/dt − dpInertia)` 在两种形状共 96 个步进状态全部成立，最差 **9.3e-17**（按那四项绝对值之和归一，`p_region` 层 1 X1）；块的自身切分 `dpFriction + dpInertia == dpBlock` 到 **1.7e-16**；16 个准稳状态 `rowError` **恰为 0.0**。
* **本例的主判据（钉死边）：形状不动账单，只重定价块那一半**——X2 层 1 账单 **1903.8603 Pa 平 → 1903.4177 Pa 塑（0.023%）**，而块的 `dpInertia` 从 **642.6188 → 643.0614 Pa**、`dpBlock` 掉 **21.4%**（沿 8 层的最大逐层相对变化 **113.7%**）。形状**能**动的是弦差：同一边 **2.47e-4 → 2.38e-2 Pa（97×）**，可它即使塑形也只占账单的 **1.3e-5**。全案 112 状态的最大弦差 **0.6423 Pa**，而最小的层 1 账单 **11.79 Pa**（`f_slab`/X1）仍是它的 **18×**——两笔账始终在不同人手里。
* **自由边的账单跟着它被收费的那个速率走**：层 1 速率 **1.917688e-5 → 2.014387e-5 m³/s（+5.04%）**、账单 **718.8665 → 754.2383 Pa（+4.92%）**——形状进不了 `M`，只进得了轨迹。
* **阶梯对形状盲（在它削掉的比例上）**：`图纸 → 区域` 削掉 X1 层 1 行误差的 **30.08%**（平）对 **29.30%**（塑），差 **0.78 pp**；X2 削掉 **15.57%** 对 **15.58%**，差 **0.00 pp**。平口三档逐字复现 §5.40（`rowError` **77.0836 / 53.8941 / 0.9364%**，账单 **718.8665 / 551.3981 / 11.7934 Pa**；X2 **52.1912 / 44.0637 / 0.3737%** 对 **1903.8603 / 1607.3812 / 13.6312 Pa**），塑形的 region 档逐字复现 §5.39 记下的 `parab_both`（**62.6559 / 56.0164%** 对 §5.39 写的 62.66 / 56.01%）——因为 §5.39 的 `network` 数就是本例的 `partner` 数，两个算例是同一臂的两种叫法（跨例回归）。
* **边界① 测到了**：`p_mix`/X2 的准稳边吃 `slab`，付 **656.25 Pa** 的 `M·dq/dt`，而那块报 `dpInertia` **0.0** 且评分恒 **0.0**；同样的 `dq/dt` 下图纸会收 **2546.48 Pa**，所以**来源藏掉 74.2% 一笔评分看不见的账**。`p_pin`/X1 的准稳边吃图纸，层 1 藏 **1264.30 Pa**、层 8 仍藏 **35.06 Pa**。而准稳**钉死**边的隐藏账单只在层 1 出现一次（**656.25 然后 0.00**），因为注入以后 `dq` 就是 0。块的回报逐边给出 `dtLevel`：说 `steady` 的那条边恰为 **0 s**，说 `Euler` 的那条 **0.2 s**——形状仍然每层重写在这条看不见时间层的边上。
* **斜率律只在拓扑钉住速率的地方精确**：X2 上两档归属的账单差恰是 `(M_a − M_b)·dq/dt`——**296.4791 对 296.4791 Pa**（`p_net`/`f_net` vs `p_region`，0.0000%）、**1890.2276 对 1890.2291 Pa**（`p_net` vs `p_pin`，0.0001%，残差 0.0015 Pa 落在块自己那 0.0560 Pa 里）。反面是自由边：同一层两臂的 `dq` 本来就差 **7.74%**（平）/ **8.62%**（塑），所以拿 `(M_a − M_b)·dq/dt` 去预测会偏 **22.6%**。同时归属**进不到块里面**：X2 的 `dpInertia` 沿 8 层在两种归属下不动到 **8.1e-6 Pa**（平）/ **1.8e-2 Pa**（塑）。
* **局部性多了一个触发**：`p_mix`/X1 与 `p_net`/X1（自己三样全同，邻居两张表都改）在**层 1..5 的 15 个字段逐字相同**，只在层 6..8 分开（`rowError` 最差 **3.1e-4**、`dpBlock` **6.8e-3 Pa**），而 `p_net` 的钉死块第一次填满预算就在**层 5**——共享前缀恰好在邻居的块放弃的那一层结束（实测巧合，不是推出来的因果）；`τ` 1.2824 对 1.2826 s，终点 0.518647 对 0.518634。反面：只改**自己的**表述（`p_pin`/X1 对 `p_net`/X1，`M` 逐字相同）就把自由边速率动 **4.681e-6 m³/s**、层 1 分流从 **0.335731** 到 **0.413742**，而那本该报出 88.6217% 差额的评分在准稳边恰是 **0.0**。
* **墙是形状的，预算由归属定价，且本例"墙"与"未收敛报告"是同一批**：平口处处不碰 300（最差 **17 扫**），塑形只在**钉死边**碰——`p_net`/X2 层 5、`p_mix`/X2 层 1,2,4,5,6,8、`p_pin`/X2 层 4,7；同一条塑形钉死边在图纸下要 **300**、在区域下峰值 **162**。这 **9** 个填满预算的状态与报 `blockConverged 0` 的 9 个**完全重合**（与 §5.39 的"6 个碰墙 / 5 个报未收敛"相反），而 9 个的界面残差最差 **8.912e-4** 仍在 1e-3 容差内、那一层照常走完——"花满扫数"仍是记录，不是界面执行的停机规则。
* **成本由两个旋钮一起定，形状是贵的那个**：趟数 **26 / 27 / 30 / 35 / 37 / 36 / 17**（共 **208**），SIMPLE 解恰为 `2 ×`（**416**），每臂最差残差 **7.279e-4 ~ 9.555e-4** 全在 1e-3 内；同归属下塑形比平口贵（**26 → 35**、**27 → 37**）。
* **斜坡由入口形状定尺度、由惯量定次序**：`τ` 平口 **0.6494 / 0.5932 / 0.4235 s**（图纸 / 区域 / slab），塑形 **1.2824 / 1.2023 / 1.0907 s**——同归属下形状把松弛时间**翻倍**（**×1.97 / ×2.03**）；层 1 分流平口 **0.319615 / 0.346431 / 0.431868**、塑形 **0.335731 / 0.367420 / 0.413742**，终点平口 **0.396201 / 0.400581 / 0.407166**、塑形 **0.518647 / 0.537958 / 0.518634 / 0.563791**。
* **伙伴仍合法、网络仍闭合**：首格 `y+` **22.8 ~ 32.6**、`nut/ν` **45.2 ~ 133.7×**、`Re_Dh` **5365 ~ 11968**、`β_out` **1.020056 ~ 1.021052**、每块进/出 ≤ **2.3e-12**；`β_in` 平口 **1.0011** 档、塑形 **1.1162 / 1.1012**（Poiseuille 的 1.18172 被壁面函数管道自己的发展段松掉），所以上面每一条差异都能归到形状或归属其一，不是第三种东西。每层节点对逐边等于**那一行的动量式**（`R·q` 加上装上的 `M·dq/dt`，1e-9），两段先导管仍是管网自己的 Churchill **11.2641281367 Pa**（1e-9），`Snk` 恒等于参考压力。

**成本**：七臂、两块 RANS 网格、8 个时间层、208 趟界面 / 416 次 SIMPLE 解；单独跑整案 **143.4 s**（`ctest -j8` 记
**153.4 s**）。同样拓扑、八臂全平口的 §5.40 是 28 s（ctest 34.58 s），同样七臂但 4 平 / 3 塑的 §5.39 是 129 s（ctest
137.6 ~ 214.4 s）——本案是 3 平 / 4 塑、208 趟，比 §5.40 那 187 趟还少，墙钟却差五倍，所以在这两案之间趟数不是墙钟时间的界。

**本增量的边界**：① 塑形 × `slab` × 步进**自由**边是发散墙，本例只把 `-89.475517 Pa @ 2.2e-5 m³/s` 记成事实，没有升格成判据，也没有解释为什么 slab 的低惯量在自由边上会把区域压成负压降。② **平口** × `slab` × 准稳边仍未跑（本例的准稳 × slab 只出现在塑形的 `p_mix`/X2 上），§5.40 边界① 只填了一半。③ 界面算子固定在 GS（`accelerates ( 0 0 0 0 0 0 0 )`），"形状 × 归属"在准 Newton 下会怎么表现（§5.33 交叉过 `du/dt` × QN）未测。④ 逐边两表 × 两条边都在并联位置（§5.32 拓扑）仍未交叉。⑤ `power1over7` 入口从未与归属同台——本例只有 flat / parabolic。⑥ `rowError` 仍然**不是判据**：本例证明了哪一笔被哪个旋钮移动，但没有把"层 1 弦差 ≤ 账单的 1e-5"变成停机条件，也没把 `slab` 变成推荐选项。⑦ `p_mix`/X1 与 `p_net`/X1 的尾部（层 6..8）分岔只量化到"恰好在邻居块填满预算之后"，究竟哪一步、哪一位开始分岔未查出来（§5.40 边界⑥ 那 2.9e-10 的同一类观测）。⑧ 1D 那一行动量方程仍不认识 `beta`（§5.34①/…/§5.40⑦ 同一格）——§5.42 把这个数逐层**报出**了，行本身仍然没有它。⑨ `dt 0.05` 那堵墙、1.2e-4 那堵发散墙、"花满 300 扫 ≠ `blockConverged 0`"仍是记录（§5.36①/§5.38④⑦/§5.39⑦⑧/§5.40⑧）。⑩ 界面上的伙伴仍是单相不可压块；多相/可压伙伴未进接口。

### 5.42 对流动量账第一次逐层报出：这笔账既不在归属那张表里、也不在行里——`rho*(betaOut - betaIn)*q|q|/A²` 在 112 个状态上恒等、符号由形状定、大小从块自身压降的 5.7%（平口）到 121%（塑形准稳边 `p_pin`/X1），而 1D 那一行仍是 `R·q + M·dq/dt`（§5.41 边界⑧ / §5.34 边界① 那一格的第一半；内核只多一个被打印的数）

§5.34 / §5.35 已经在界面上记过一本动量通量的账，但那是**入口那一半**：`momentumFluxPa == rho*beta*q^2/A`，伙伴是层流 duct，它对区域自报压降的比值是 **4.90e-5 ~ 1.33e-4**，所以当时的结论是"形状只到块的入口条件，1D 那一行还不认识 `beta`，而这个比值就是它还**不必**进账的理由"。§5.36 起伙伴换成壁面函数 RANS 管道，§5.41 又把入口形状放上逐边归属表，于是"进出口形状不同"这件事第一次有了量级：区域的 `beta_out` 在 112 个状态上是 **1.0141 ~ 1.0211**（平口档 1.0141 ~ 1.0211、塑形档 1.0141 ~ 1.0206，它跟着发展段走，几乎不看入口是什么），而 `beta_in` 平口 **1.0010 ~ 1.0012**、塑形 **1.0962 ~ 1.1182**——两者之差平口 **+0.0130 ~ +0.0200**、塑形 **-0.0789 ~ -0.0980**，不再是小数点后的离散性质，它是一个和块自己那笔压降同量级、甚至可以超过它的数。这不改变行的写法，只是把"还不必进账"那句理由**测掉**。

算例、算例目录、七臂、两块 `kOmegaSST` 网格、`dt 0.2`、8 个时间层、1e-3 容差、GS、逐边 `ddtFrom` / `inertiaFrom` / `profile` 三张表全部照抄 §5.41（`Src → A → { X1 ‖ By } → X2 → B → Snk`，6.0e-5 m³/s 在 t=0 接入），**没有新增 case 关键字**：新增的只有内核里每个块每次求解多算的一个减法，和 `transient step` 行里 `pRiseX` 之后、`unsteady` 之前的那一个字段 `dpConv`。

声明面：字段 `ProfileLedger::dpConvectivePa` 在 `src/coupling/include/cf/coupling/CoupledNetwork.h:58`，取值在该文件 `:314`–`:319`（用同一次求解的 `betaIn` / `betaOut`、该层的速率与实测入口面积），随 `EdgeExchange`（`:490`）在 `src/coupling/src/coupledNetwork.cpp:1226` 交给边界，打印在 `src/coupling/apps/transientCoupledMultiApp.cpp:527`。**行算术一处没动**：`coupledNetwork.cpp:1200` 仍是 `rowDrop = R*q + M*dq/dt`，弦仍取 `unsteady ? dpFriction : dpBlock`，`rowError` 仍是 §5.40 那两笔。判据面是 §5.41 那两个文件（case + checker），老 158 条**一条未改就全绿**——这就是"行为惰性"的实测形式；只有局部性那条的字段集从 15 扩到 16（把新字段纳入逐位比较，见下面最后一条）；新增第 12 族 **11** 条，共 **169** 条（RC=0）。

* **恒等与符号，逐状态**：112 个（层 × 边）状态每一个都满足 `dpConv == rho*(betaOut - betaIn)*q|q|/A^2`，最差相对偏差 **1.44e-15**；符号恰好跟着形状走——48 个平口状态全 **正**（区域得长出比塞状入口更饱满的剖面，`f_net`/X1 层 8 是 `betaOut` 1.0209 对 `betaIn` 1.0011）、64 个塑形状态全 **负**（它卸掉 Poiseuille 的过量，同一格是 1.1162 降到 1.0201；`p_net`/X2 层 8 则是 1.1012 降到 1.0206），无一例外。
* **它有多大**：平口 `f_net` X1 **7.632 → 17.924 Pa**、X2 **57.615 → 69.844 Pa**（沿 8 层；X2 那一列三个平口臂一致到 **5.7e-8**，X1 在层 1 是 **7.632 / 9.157 / 15.044 Pa**，见下面归属那条）；塑形 `p_net` X1 **-55.239 → -148.974 Pa**、X2 **-277.295 → -283.272 Pa**。占那块自己报的压降：平口最坏 **5.73%**（`f_slab`/X1 层 3 的 19.715 Pa 对 344.36 Pa），塑形 `p_net`/X2 平台 **39.3%**，`p_pin`/X1 层 8 **121.1%**（**-179.391 Pa** 对 **148.15 Pa** 的整块压降——那还是条准稳边）。§5.35 那个"≤2e-4 所以不必进账"的比值，在这个伙伴下翻了三个数量级。
* **归属那张表进不了这笔账**：钉死边三个平口臂的 `dpConv/q²` 在全部八层上一致到 **5.7e-8** 相对（同一批状态里块自报的 `dpInertia` 也是同一个 **642.6188 Pa**），而那张表在该边装上的惯量相差 **3.64×**；自由边同一层里 `dpConv` 本身差到 **97.1%**（层 1：**7.632 / 9.157 / 15.044 Pa**），但那里速率本身就差 **35.1%**，除掉 `q²` 之后只剩 **7.96%**，且逐层单调塌到 **0.69%**——和 §5.41 的账单同构：形状的价格只通过轨迹进来，不通过"记了谁的惯量"。
* **形状才是那个旋钮**：同一个系数平口 **1.9401e10** 对塑形 **7.8687e10 Pa·s²/m⁶**，**4.06×**。
* **与压强回升同号、同量级、但不相等**：56 个层 ≥2 且回升 >50 Pa 的塑形状态上 `|dpConv| / pRise` 落在 **[0.8475, 1.2989]**，X2 平台（层 4..8 × 四个塑形臂）收窄到 **[0.84760, 0.85110]**（`p_net`/X2 层 8：**-283.27 Pa** 对 **334.20 Pa** 的回升）；层 1 在带外（**2.0786**，回升只发育到 133.40 Pa）。平口处处报 `pRise == 0.0` 却仍欠一笔正的账（`f_net`/X2 层 8 **+69.844 Pa**），所以回升是这笔失衡的**症状之一**，不是它的金额。
* **它也不是 §5.40 那笔双重计费**：账单超过 100 Pa 的 **9** 个状态（全部在层 1）上 `dpConv` 至多占 **17.3%**，而且在塑形臂上符号与账单**相反**——把对流项收进行里既关不掉那 9 个状态上 **44.06% ~ 88.62%** 的 `rowError`，在那几条边上还会放大它；反过来 **16** 个步进状态的账单已经消失（`|rowDrop - dpBlock| < 0.01 Pa`）却仍欠一笔实账：**15** 个是平口 X2 平台各 **+69.84 Pa**（`f_net`/`f_region`/`f_slab` 层 4..8），**1** 个是塑形 `p_region`/X2 层 7 的 **-283.27 Pa**。两笔账，两个人。
* **准稳边同样欠**：`p_mix`/X2 每层 **-283.670 ~ -283.783 Pa**，而它的行按定义不装惯性项（`unsteady 0`、`rowError` 恰 0.0）；`p_pin`/X1 **-90.637 → -179.391 Pa**。§5.41 那句"来源藏掉一笔评分看不见的账"现在有了第二笔，而且**行也看不见它**。
* **行仍是那一行，且新字段不引入新的非局部**：96 个步进状态上 `rowDrop == R*q + M*dq/dt` 到 **1e-6**；§5.41 的局部性前缀（`p_mix`/X1 对 `p_net`/X1，层 1..5）把判据字段集从 15 扩到 **16**（加入 `dpConv`）后仍然**逐位相同**，尾部照旧只在层 6..8 分开。

**成本**：七臂、两块 RANS 网格、8 个时间层、**208** 趟界面 / **416** 次 SIMPLE 解、112 个（层 × 边）状态里 96 个步进 +
16 个准稳、每臂最差残差 **7.279e-4 ~ 9.555e-4**——四组数与 §5.41 逐字相同（新增的是一个减法与一次打印）；单独跑整案
**199.4 s**（`ctest -j8` 记 **284.17 s**，全量 97 项 **297.95 s**）。这两条墙钟都要打折扣读：本机跑本案时 `load avg` 在
**7 ~ 11** 之间，§5.41 记那对 **143.4 / 153.4 s** 时是 **4.7**，而同批里 §5.34 的稳态 profile 映射也从 **255.9 s** 涨到
**297.9 s**——同一台机器上一起涨的那部分不是本案的成本。趟数、块解数与残差这三样才是本案自己说的数。

**本增量的边界**：① **报出 ≠ 收进**。把 `dpConv` 写进 1D 那一行动量方程需要（a）逐臂声明谁为这条边付这一项，以及（b）让块的弦半区停止重复计它——否则它就是 §5.40 那笔账单的 β 版本，一行里两处收同一笔钱。本增量两者都没做，行仍是 `R·q + M·dq/dt`。② 那一格 121% 只是**记录**，不是判据：界面照旧在 1e-3 容差内走完 8 层，没有"层内 `|dpConv|` 超过块压降多少就停"的规则。③ 为什么 `|dpConv| / pRise` 是 **0.85** 而不是 1、层 1 为什么落在 **2.08**，没有解释；两者在物理上同源于通量松弛，机制未查。④ 平口 × `slab` × 准稳边、塑形 × `slab` × 步进**自由**边（那堵 `-89.475517 Pa @ 2.2e-5 m³/s` 的墙）仍是 §5.41①② 那两格，β 账在其中只有推理、没有数。⑤ `power1over7` 入口从未与这笔账同台（本例只有 flat / parabolic）。⑥ 准 Newton（§5.41③）与双边并联拓扑（§5.41④）未测，`dpConv` 在它们下面会怎么动未知。⑦ 尾部（层 6..8）分岔、`es_mix` 那 2.9e-10 一类观测继承未查。⑧ 界面上的伙伴仍是单相不可压块，`dt 0.05` / 1.2e-4 / "花满 300 扫 ≠ `blockConverged 0`" 三堵墙仍是记录。

### 5.43 对流动量账收进 1D 那一行：`rowDrop = R·q + M·dq/dt + dpConv`，而块的弦同时变成 `dpFriction − dpConv`——两侧相加仍是块那半区摩擦，所以这是一次改名而不是第二笔账，差别只从行的导数里进来（§5.42 边界① 的第二半）

§5.42 把 `dpConv = rho*(betaOut - betaIn)*q|q|/A^2` 逐层逐边**报**出来了，并记下它占块自身压降从 **5.73%** 到 **121.1%**，同时把"为什么还不能直接写进行里"说清楚：那一行的阻力项 `R·q` 与块的弦是**同一个不动点**的两半——界面对 `R` 迭代，直到 `R·q` 追上块的弦。如果只把 `dpConv` 加到行上而弦照旧，同一笔钱就在一次平衡里收了两次（正是 §5.40 那笔账单的 β 版本）。所以这一增量的两半必须**同时**发生：

```
行：   rowDrop = R·q + M·dq/dt + dpConv          （dpConv 由系数 C = dpConv/q² 装上，行里写 C·|x|·x）
弦：   chord   = dpFriction − dpConv             （界面把 R 追上的目标少掉这一项）
```

于是 `R·q + chord` 恒等于块的摩擦半区：**改名，不是新账**。但改名不是无事发生——行的导数变了：`∂(C|x|x)/∂x = 2C|x|` 进了 Jacobian、`|C·|x|·x|` 进了 Newton 的 `scale`，而弦对 `q` 线性、形状项对 `q` 二次。不动点因此只从**轨迹**里进来，这正是 `c_free`（装上自由边）与 `c_pin`（装上钉死边）分岔的地方：拓扑钉住速率的那条边，`q` 无处可去，改名就只是改名。

算例 `transientCoupledMultiTurbulentConvectiveRow` 把 §5.42 的拓扑、两块湍流伙伴（X1 `kOmegaSST`、X2 RANS）、`dt 0.2`、8 层、1e-3 容差、GS、逐边 `ddtFrom` / `inertiaFrom` / `profile` 三张表全部照抄，**唯一新增的声明面是第四张逐边表** `convectiveToRow { <臂> { X1 charge|report; X2 …; } }`，四臂 = `{c_off: report/report, c_both: charge/charge, c_pin: report/charge, c_free: charge/report}`。`c_off` 是行为惰性对照（什么都没人付，就是 §5.42 的世界），`c_pin` 是改名对照，`c_free` 才是那句真实的陈述。

声明面：`PipeDyn::convective`（Pa·s²/m⁶，`src/1d/include/cf/1d/Transient.h:189`）由 `TransientNetwork::setCoupledConvective`（声明 `:157`、实现 `src/1d/src/transient.cpp:429`，与 `setCoupledResistance` 同样要求该管已是耦合管，且因为这项**带符号**只查 `isfinite`）写入，动量行在 `src/1d/src/transient.cpp:248`–`:252` 收下它、`:257`–`:259` 把它计入 `scale`；界面上 `ProfileLedger::dpConvectivePa`（`src/coupling/include/cf/coupling/CoupledNetwork.h:60`）经由 `EdgeExchange::convectiveRowPa`（`:508`，报出那一侧恒 0）与新的 `chordPa`（`:466`）在 `src/coupling/src/coupledNetwork.cpp:1232`–`:1235` 分成两半，残差改用 `chordPa`（`:1255`–`:1256`）、`rowDropPa` 加回它（`:1266`）、QN 与 GS 两处刷新都读 `chordPa`（`:1289`、`:1303`），最后在 `:1317`–`:1318` 把 `C = convectiveRowPa/q²` 装回网络；开关 `TransientCoupledControls::convectiveIntoRowPerEdge` 在 `CoupledNetwork.h:660`，校验与解析在 `coupledNetwork.cpp:1133`–`:1145`（长度必须等于耦合边数，词表只有 `report`/`charge`），case 侧解析在 `src/coupling/apps/transientCoupledMultiApp.cpp:332`–`:356`，逐边生效的词打印在 `:483`–`:484`，每层的 `dpConv` 与 `rowBeta` 打印在 `:580`、`:588`。装上的次序不是随便排的：`net.pendingPipeFlow()` 一旦 `pending_` 为假就 fatal，所以这一项和 `setCoupledResistance` 一样**只谈不拢的那几趟**才装。判据是新算例的一对（case + checker），**35** 条全绿（RC=0）。

* **惰性到逐位**：`c_off` 与 §5.41 的 `p_net` 在 16 个字段 × 16 个（层 × 边）状态上**差 0.0**（首层与末层各按出厂数字复核：X1 层 1 `dpBlock` 851.08258956984685、层 8 177.93872006380388），每层界面趟数 **[3, 6, 5, 5, 4, 4, 4, 4]** 一字不差，`rowBeta` 处处 0.0——一个装上去却没人读的开关改不了"一层需要几趟"。同一句话的另一种证明是 §5.42 那 **169** 条判据：本案一行判据未改，重跑仍然全绿。
* **不收第三笔**：32 个装上状态 `rowBeta == dpConv`、32 个报出状态 `rowBeta == 0`，而全部 64 个状态上 `rowDrop − (dpNetwork + M·dq/dt + rowBeta) == 0.0`——行的总额仍然是节点对它的那纸账，只是多了一项署名。
* **弦确实停止重复计它**：内核自己的 `residual` 字段现在**逐字等于** `|dpNetwork − (dpFriction − rowBeta)| / chord`，最坏相对弦差 **9.495e-4**（`c_both`/X1，容差 1e-3 之内）、`c_off` **7.958e-4**，绝对剩余 ≤ **0.595 Pa** 对 **148.97 Pa**（X1）与 **283.27 Pa**（X2）两本账。这不是推论，是被界面自己那条判据读出来的。
* **改名的价格是算得出的**：`R_charged / R_off == 1 − dpConv / dpFriction_off` 在每层每条边上到 **1.722e-3**（X1）/ **2.683e-4**（X2），于是图纸那个系数不再是壁的数：X1 涨 **9.5% ~ 90.5%**、X2 涨 **12.5% ~ 39.3%**——形状账占摩擦半区的比例就是涨幅本身。
* **钉死边上只是改名**：`q` 在 4 臂 × 8 层逐字仍是 **6.0000000000000002e-05**；`c_pin` 对 `c_off` 在这条边上只有 `R` / `dpNetwork` / `rowBeta` 三项动，其余字段最大移动 **3.521e-2 Pa**（`dpBlock`）、行总额 **4.074e-2 Pa**、`rowError` **1.054e-4**，对 283.3 Pa 的账——那是界面停机的噪声，不是账。
* **声明只到自己那条边**：`c_pin`/X1 对 `c_off`/X1 **17/17 字段逐位 0.0**；`c_both`/X1 对 `c_free`/X1 同样 **17/17 逐位 0.0**。和 §5.35 为归属 bought 的那条可反证性同一条：邻居的声明改不了我这行的状态。
* **自由边上这是一条真轨迹**：X1 的 `q` 相对移动最大 **1.423e-4**（层 8 **3.111440e-5** 对 **3.111883e-5**，因为行里多了一项对速率二次的项），`dpBlock` **6.052e-4**，而行的**总额**只移动 **0.125 Pa**，对搬走的 **148.9 Pa** 是 **8.4e-4**——改名买的是路径，不是账单。多出来的那份由 1D 支管自己吃掉：层 8 `By` **2.88856e-5** 对 **2.88812e-5**（差 **4.427e-9**），而 `X1 + By` 在 4 臂 × 8 层上逐字等于注入速率（最差 **6.8e-21**）。
* **§5.40 那两只括号结构不变**：`rowError·dpBlock == (dpNetwork − chord) + (M·dq/dt − dpInertia)` 在 64 个状态上最差 **6.821e-13**，弦那半仍被容差压在 **0.5949 Pa** 以内、归属那半仍是 **754.2383 Pa**（X1）/ **1903.4177 Pa**（X2），两者相差 **2604 ~ 6664×**——把形状账搬走没有把它变成第三笔。
* **谁付这一项，改不了这一项值多少**：`|dpConv|/q²` 跨四臂在 X1 上一致到 **2.393e-5**、X2 到 **1.008e-4** 相对，每层都成立（X1 在层 1 与层 8 之间是 **1.3613e11 ~ 1.5384e11 Pa·s²/m⁶**，X2 是 **7.7026e10 ~ 7.8687e10**，后者与 §5.42 记下的塑形系数逐字相同）；`rowError` 对 split 也不敏感（钉死边 **1.054e-4**、自由边 **1.4e-3**）。这笔账是**区域**的性质，不是那行的。

**成本**：四臂、8 个时间层，`c_off` 与 `c_pin` 各 **35** 趟界面 / **70** 次块解，`c_both` 与 `c_free` 各 **34** 趟 / **68** 次——装上形状项没有让界面更贵，最坏单级 **7** 趟（`c_free` 层 2，上限 12），而 `c_free` 的层 7..8 反而从 4 趟塌到 3 趟（逐层序列 `c_off` **[3, 6, 5, 5, 4, 4, 4, 4]** 对 `c_free` **[3, 7, 5, 5, 4, 4, 3, 3]**）。

**本增量的边界**：① 新的第一堵墙在代码里而不是在结果里：`chord = dpFriction − dpConv` 一旦非正就 `cfFatal`（`coupledNetwork.cpp:1235`）。§5.42 记下过 **121.1%** 那一格（准稳边、塑形、`p_pin`/X1 层 8），那种边上弦付不起形状账；本案四臂都在层流剖面的发展段上取 `dpFriction`，`dpConv` 最大只占 **90.5%**，所以这条 fatal 从未触发——它是**记录**，还没有算例把它变成判据。② 只有 `report`/`charge` 两档，没有"部分付"；`C` 每层重算，因此带轨迹，但没有把它当阻尼项处理。③ 伙伴仍是那两块单相不可压 RANS 管道，`profile` 仍是 flat / parabolic 两档，`power1over7` 未与这项同台。④ 准 Newton 与双边并联拓扑（§5.41③④ 那两格）未测；在 QN 下弦的刷新走 `chordPa`（`:1289`）这条路只被本案的 GS 臂之外验证过一次，形状账与加速算子如何相互作用未知。⑤ 尾部（层 6..8）分岔、`es_mix` 那 2.9e-10、`-89.475517 Pa @ 2.2e-5` 那堵 slab 墙、`dt 0.05` / 1.2e-4 / "花满 300 扫 ≠ `blockConverged 0`" 全部继承未查。⑥ 装了 `dpConv` 并不等于消掉 §5.40 那笔 `M·dq/dt` 双重计费——本案归属表全档 `network`，账单那半仍在括号里（**754 / 1903 Pa**）。

### 5.44 两张"谁付"表落在**同一条边**上：归属那三条阶梯不被改名移动，形状那笔账不被阶梯移动，而 §5.40 分解出的两只括号第一次**同时**收缩（§5.43 边界⑥）

§5.43 把形状账收进行里，然后在自己的边界⑥里记下：这不等于消掉 §5.40 那笔 `M·dq/dt` 双重计费——因为本案的归属表全档停在图纸上，账单那半 **754.2383 Pa**（X1）/ **1903.4177 Pa**（X2）原地未动。这一增量不新增任何声明面，只是把已有的两张表放到同一条边上问一句：**那是两根杆子还是一根？**

```
行：   rowDrop = R·q + M·dq/dt + rowBeta        rowBeta = dpConv（charge）或 0（report）
弦：   chord   = dpFriction − rowBeta
恒等： rowError·dpBlock = (dpNetwork − chord)  +  (M·dq/dt − dpInertia)
       \________ 归 convectiveToRow 管 ________/   \____ 归 inertiaFrom 管 ____/
惯量三档（X2，闭合式）：network ρL/A = 8.488264e6 · partner ρV/A² = 7.500000e6 · slab ρV_probe/A² = 2.187500e6 kg/m4
```

算例 `transientCoupledMultiTurbulentShapeAndInertiaSplit` 是 **3 档 × 2 位置 = 六臂**的因式表，两张表都只写在 **X2**（拓扑钉死速率的那条边）上，X1 在六臂里始终是"图纸 + 报出"。为什么必须落在 X2：① 归属只在 `dq ≠ 0` 处可见，而 X2 除层 1 外 `dq` 恒 0——偏偏层 1 那笔 **1903 Pa** 是 §5.40 那个分解里最大的数；② `slab` 落在**自由**边上会发散。这一点是本案写作过程中重新量到的：把 `slab` 同时落到两条边的第一版跑到 `s_rep` 臂就停在

```
the detailed region reports a non-positive drop (-89.475517 Pa) at 0.000022 m3/s
```

——与 §5.41 记下那两个数**逐位相同**（148.47 s 处 fatal，前四臂完整）。因此本案的落位不是省事，是唯一既能把阶梯看清、又不会撞上那堵墙的落位。

声明面：一行内核代码也没有改。`inertiaFrom`（§5.29 的机制、§5.40 逐边化）与 `convectiveToRow`（§5.43）各自早就存在，本案只是同时打开它们并把两者都指向同一条边；`inertias ( network … )` 保持全档图纸，让逐边表成为归属轴的**全部**陈述。判据是新算例的一对（case + checker），**34** 条全绿（RC=0）。

* **阶梯的第一级是出厂数字**：`n_rep` 层 1 的 X2 账单 **1903.4177 Pa** 由 §5.41 发表的四个字段（`rowError` × `dpBlock` + `dpFriction` − `dpNetwork`）重新算出并对上（**<1e-9**），第二、三级是 **1606.9386**（区域）与 **13.1886**（薄块）Pa，一次 **144.3×** 的下台阶；三档装上的 `M` 逐层逐边等于各自闭合式（相对误差 <1e-9），X1 六臂都是图纸的 **1.018592e7**。
* **改名搬不动阶梯**：同一档内 `M` 在 8 层上**逐位相同**（三对差 **0.0 / 0.0 / 0.0**），块自己报的 `dpInertia` 在区域档与薄块档逐位相同、图纸档差 **2.273e-3 Pa**（对 643.1 Pa 是 **3.5e-6** 相对），三档的层 1 账单在 `report`/`charge` 两臂之间差 **0.0 Pa**。形状的署名与惯量的署名是两个地址。
* **阶梯也搬不动形状账**：X2 的 `dpConv` 在**六臂全部**层 1..7 上**逐位相同**，层 8 一致到 **2.8e-6** 相对，而它旁边的 `M` 相差 **3.8803×**；除掉 `q²` 的系数在 **7.7026e10 ~ 7.8687e10 Pa·s²/m⁶** 之间，与 §5.42/§5.43 记下的塑形系数逐字相同。§5.43 证明"谁付这一项改不了它值多少"，本案证明"装哪个惯量也改不了它值多少"。
* **两只括号第一次同时收缩**：层 1 的 X2 括号之和从图纸档报出的 **1903.44 Pa** 到薄块档装上的 **13.21 Pa**，**144.1×**；行自己的误差随之 **0.6635 → 0.5602 → 0.0046**。而**同一档内**装上与否只差 **1.4e-7**（`0.663516/0.560164/0.004589` 对 `0.663516/0.560165/0.004589`）——付掉形状那一半并不付掉惯量那一半，§5.43 边界⑥那句"仍未消掉"第一次是被读的而不是被推的。
* **归属不会把弦谈松**：六臂的最坏相对弦差**同一个数 7.9585e-4**（1e-3 容差之内；报出臂绝对 **0.5658 Pa**、装上臂 **0.5949 Pa**），换档买到的只是行里那个隐式项更小，不是界面更宽。
* **两个声明都只到自己那条边**：其余五臂的 X1 对 `n_rep` 的 X1 在 **17 字段 × 8 层逐位 0.0**，且 X1 的账单六臂恒 **754.2383 Pa**——§5.35 为归属、§5.43 为形状各自买到的那条可反证性，在这里是同一条定律的两次同时兑现。
* **仍不收第三笔**：24 个装上状态 `rowBeta == dpConv`、72 个报出状态恒 0（后者要么在 X1、要么是 `_rep` 臂的 X2），96 个状态上 `rowDrop − (dpNetwork + M·dq/dt + rowBeta) == 0.0`，§5.40 的恒等式在 96 个状态上最差 **6.821e-13**。
* **对照与闭合**：`n_rep` 逐字复现 §5.41 的 `p_net`（层 1 / 层 8 各按出厂字段复核，趟数 **[3, 6, 5, 5, 4, 4, 4, 4]** 一字不差），`X1 + By` 在 6 臂 × 8 层上恒等于注入速率（最差 **6.78e-21**）。

**成本**：六臂全部**做一样多的工**——每臂 **35** 趟界面（逐层 **[3, 6, 5, 5, 4, 4, 4, 4]**，最坏单级 **6** 趟，上限 12）、每臂 **70** 次块解（合计 420）。两张表都不是迭代，而是写进行里的系数，所以交叉不花钱。墙钟 **249.85 s**（user 247.78 s）。

**本增量的边界**：① `slab` 仍**没有**被判据落在自由边上——把它落在那里会发散，这一点本案又逐位重演了一遍（`-89.475517 Pa @ 0.000022 m3/s`），但"重演"仍是记录而不是判据，§5.43① 那堵 `chord ≤ 0` 的墙同样仍未成为判据。② 阶梯只在**层 1** 有意义（钉死边其余层 `dq ≡ 0`），所以"两只括号同时收缩"目前是一条**首层**陈述，平台尾没有这笔账可收。③ 薄块档仍留 **0.004589** 的 `rowError`（= **13.19 Pa**）——那是界面容差的地板而不是缺陷，但"把两只括号都压到 1e-3 以下"这条判据还不存在。④ 交叉只做了**一条边**、GS 一种算子、parabolic 一种形状、两块单相不可压 RANS 伙伴；X1 侧逐边混档（§5.40 的 `se_mix` 那一格）× 形状装上未测，`power1over7`、准 Newton、双边并联拓扑（§5.41③④）同样未同台。⑤ 尾部（层 6..8）分岔、`es_mix` 那 2.9e-10、`dt 0.05` / 1.2e-4 / "花满 300 扫 ≠ `blockConverged 0`" 全部继承未查。⑥ X1 的 **754.2383 Pa** 仍欠着——本案只在算例敢声明档位的边上收账。

### 5.45 非牛顿管流第一次进稳态管网：Metzner-Reed 广义雷诺数**替代** `rho·u·d/mu` 那一个位置，于是 Churchill 的层流半区从"关联式"变成"精确解"——三个结构指数、一条与 K 无关的壁应力，和两条被读出来的拒收

对标登记里的缺口 12（非牛顿流变）与 Products 侧的事实是同一句话：Pipesim / OLGA 的流体模型库里幂律与 H-B 是**默认能力**，而 ChuanFlow 从 5.1 到 5.44 的每一条摩阻陈述都写在 Newtonian 上。这一增量只做**一次替换**——不是新增一条摩阻关联式，而是把 Darcy-Weisbach 里问 `Re = rho u d / mu` 的那一处换成 Metzner-Reed 广义数，然后**照旧**调用同一个 Churchill：

```
Re_MR = rho u^(2-n) d^n / (K' 8^(n-1))        K' = K ((3n+1)/(4n))^n
f     = churchillFriction(Re_MR, eps/d)       dp = f (L/d) rho u^2/2 + K_valve rho|u|u/2 + rho g dz

n = 1 且 K = mu 时  K' = mu·4/4、8^0 = 1、u^1 = u，逐项就是 rho u d / mu
tau = K gdot^n 的层流精确解：dp = (L/d) K (3n+1)/n (8u/d)^n
                          Q  = (pi d^3/32) [ n dp d / ((3n+1) K L) ]^(1/n)
```

关键在于最后那两行**不是拟合**：Churchill 的层流半区是 `f = 64/Re`，把 `Re` 换成 `Re_MR` 之后它就是幂律圆管流动的解析解（n=1、K=mu 时退回 `128 mu L Q = pi d^4 dp`）。所以判据可以逐条打在闭式上，而不是打在"和关联式差不多"上。**这句"就是解析解"依赖 `K'` 上的那个指数 `n`**：Rabinowitsch 的精确壁面剪切率是 `gdot_w = (3n+1)/(4n) · 8u/d`，把它代进 `tau = K gdot^n` 就得到 `K' = K ((3n+1)/(4n))^n`；5.45 当时写的是教科书里那个**指数为 1** 的塌缩约定 `K (3n+1)/(4n)`，两者差 `(4n/(3n+1))^(1-1/n)`——`n = 1` 处为 1、`n = 0.5` 处是 **1.25**（25% 的流量、20% 的压差），而本案自己的 `PowN1 == Newt` 逐位对照恰好看不见它。这一格由 5.47 在 HB 的 `tau_0 -> 0` 极限上量出并改正，下面所有数字都是**改后重测**的。算例 `steadyPowerLawRheology` 把九根管的两端都接在**定压节点**上（A 3e4 / B 2e4，`case_dpDouble.dict` 把 A 抬到 4e4），于是每根管恰好承担 `dp = 1.0e4 / 2.0e4`，牛顿解不夹在断言和物理之间，九根管是九条彼此独立的陈述：`Newt`（未声明）、`PowN1`（n=1, K=mu）、`Thin05`/`Thin07`（n=0.5/0.7）、`Thick13`（n=1.3，同号反面同台）、`D2xNewt`/`D2xThin`（孔径翻倍）、`LongNewt`/`LongThin`（长度翻倍）。

声明面：`NetworkPipe::flowIndex` / `consistency`（`src/1d/include/cf/1d/Network.h:73`–`:74`，缺省 0 = Newtonian），解析与四条拒收在 `src/1d/src/network.cpp:1118`–`:1132`，广义数在 `:101`–`:117`，两个使用点是求解用的 `pipeLossAt`（`:239`–`:261`）与**打印**用的 `generalizedReynolds`（调用点 `equipmentFlow`，`:369`）——同一处代数，报出的那个数就是实际用了的那个数。判据是新算例的一对（四个 case 文件 + checker），**16** 条全绿（RC=0）。

* **精确到解，不是到容差**：九根管 × 两个驱动差 = 18 个速率对闭式最差 **3.249e-12** 相对（两个驱动差上同值）；报出的 `Re` 对 `Re_MR` 闭合式差 **2.204e-12**（12 位打印本身就在 5e-13 量级，改正后的常数又多进几位 `pow()`）；`f · Re_MR == 64` 在 18 个状态上最差 **1.563e-10**——这条把"上面那个闭式确实是内核用的那条律"从推理变成读数，而 `v == Q/A` 只差 **4.376e-12**，那是**打印精度**不是物理。
* ** declared-but-Newtonian 是逐位的**：`PowN1` 对 `Newt` 在 Q / v / Re / f / dp 五个字段、两个文件里**全部相等**（不是接近，是 `==`）——广义数在 n=1、K=mu 处逐项塌回原式，所以这一步替换对 98 条既有测试**结构上不可能**有影响。
* **三个结构指数**：孔径翻倍，Newtonian 涨 **16.000000** 倍、n=0.5 涨 **32.000000** 倍（指数是 `(3n+1)/n`，4 变成 5——剪切稀化管线的**尺寸律**和 Newtonian 不是同一个幂）；长度翻倍两者各除以 **2.000000** 与 **4.000000**（对 `dp` 的指数变了，对 `L` 的指数**没变**）；驱动差翻倍九根管各涨 **2^(1/n)**（实测 2.0000 / 4.0000 / 2.6918 / 1.7044 / 4.0000 / 4.0000，最差 **3.979e-12**）——`n>1` 与 `n<1` 在同一次运行里读那两个符号的反面。
* **壁应力不是材料属性**：共享孔径与长度的五根管 `tau_w = dp·d/(4L)` **逐位相同**（2.5 Pa），而 `K'(8u/d)^n` 在五根流变管上各自复现它自己的 `tau_w` 到 **5.413e-12**——流体可以自选**剪切率**，不能自选**壁应力**；`Thin05` 2.5 Pa 对 `D2xThin` 5.0 Pa 的差只来自几何。
* **两条墙第一次是判据而不是记录**：声明 `flowIndex` 而不给 `consistency` 被**拒收**（n≠1 时 `mu` 与 `K` 不是同一个量纲，默认值就是一个静默的错数）；耦合元带上流变被**拒收**（它的损失**就是**伙伴模型的陈述，5.43 为形状账立的那条不重复计费规则在这里没有可以承载流变的摩阻项）。checker 断言的是**非零退出码 + 消息**，不是任何物理量。
* **层流半区确实全程在位**：18 个状态里最大的广义雷诺数是 **170**（翻倍驱动差下的宽孔径 Newtonian 管），比 Churchil 的换段低一个量级——这正是上面那些断言敢写成"精确"而不是"容差"的原因。

**成本**：四个 case 文件、16 条判据、墙钟 **0.14 s**（`ctest -j8` 记 102 项那一轮；单次稳态求解是毫秒级，本案把 4 次运行打包进一条 CTest）。

**本增量的边界**：① 幂律**只有幂律**——Bingham / Herschel-Bulkley 的屈服应力在此架构里不是"再换一次 Re"：`loss(q→0) → 4 tau_0 L/d ≠ 0`，`equipmentFlow` 那条在 `|q|` 上展开的 Illinois 括号会在起始处拿到非零残差并落入错误的分支选择，**屈服是一个死分支问题，不是一次代入**，所以它单独成格：5.46 落了 Bingham（`n = 1`），5.47 落了 Herschel-Bulkley（`n ≠ 1` 的那条单调形），而 5.47 在 `tau_0 -> 0` 极限上量出本节 `K'` 的指数本该是 `n` 不是 `1`——那一格的常数错在这里被自己儿子的极限抓住。② 湍流半区被广义数**够到了但没有被验证**：本案全部状态停在层流侧，`Re_MR > 2100` 时把 Newtonian 的 Churchill 湍流项套在幂律上是行业惯例（Metzner-Reed 原始陈述），此处**只是使用，不是判据**。③ 只接在**常物性单相稳态律**上：PVT 表驱动的黑油/气柱（`network.cpp:656` 那一侧的逐步 `march`）与两流体分布式内核都明确**拒收**流变声明，瞬态 1D 那一行（`transient.cpp:72`–`:80`）当时仍是 Newtonian 摩阻——**这一句已被 5.48 解除**：那一行现在调用的就是这里同一个 `pipeFrictionPa`，声明的幂律在瞬态侧被读；而屈服在瞬态侧当时改成**拒收**，**这一半已由 5.49 解除**（那一行现在有自己的静止支：`transient.cpp:40`–`:44` 的入口闩锁与 `:353`–`:362` 的穿过钳制，阈值就是这里这个 `4 tau_0 L/d`）。④ 局部损失（`valve K`）仍按 Newtonian 的 `rho u|u|/2` 计，未与流变同台；`roughness` 与广义数的组合未测。⑤ 圆管：本案全是 `d`，矩形/环空的水力直径与 `K'` 修正（Dodge-Metzner）不在内。⑥ `n→0`（近似完美塑性）与 `n` 极大未扫，`flowIndex 0` 是"未声明"的哨兵而不是一个可声明的极限。

### 5.46 屈服应力第一次进稳态管网：死分支在括号**之前**判掉，流量落在 Buckingham 1921 的精确解上——`dp < 4 tau_0 L/d` 时 `Q` 打印出来就是 0，而报出的摩阻系数换成解自己含的那个，于是 `f · Re · g(xi) == 64` 成了一条可失败的断言（§5.45 边界① 的那一格）

5.45 把"屈服不是又一次代入"写成了边界，这一增量就是那一格的兑现。差别在第一行代数上：Newtonian 与幂律的 `loss(q)` 都过原点，所以 `equipmentFlow` 可以在 `|q|` 上开一个 Illinois 括号慢慢展开；带上 `tau_0` 之后

```
loss(q -> 0) = dpStatic = 4 tau_0 L / d  !=  0
```

被问到的压差如果不到这个数，**精确答案就是零流量**，而 `g(0)` 此时已经落在根的远侧——括号会朝着同一个符号的采样展开 60 次，然后拿两点割线算出一个不存在的根。所以这条分支必须在括号之前显式判掉（`src/1d/src/network.cpp:315`–`:327`，只在 `yieldStress > 0` 时生效，牛顿与幂律那两条的 `g(0)` 本来就是送它们来这边的符号，一行未动）。流量侧用的是圆管 Bingham 层流的**精确解**，反解给的是网络要的 `loss(q)` 方向：

```
Q = (pi d^3 / (32 mu)) tw [1 - (4/3) xi + (1/3) xi^4]      tw = dp d / (4L),  xi = tau_0 / tw
  ->  F(tw) = tw - (4/3) tau_0 + tau_0^4 / (3 tw^3) = 32 mu Q / (pi d^3) = B
      F(tau_0) = 0、F' = 1 - (tau0/tw)^4 > 0、F(B + 4 tau_0/3) >= B
      => 括号 [tau_0, B + (4/3) tau_0] 是解析的，Newton + 二分落到舍入位（:131–:160）
f = 8 tw / (rho v^2)   =>   f Re g(xi) = 64,   Q = Q_Hagen-Poiseuille · g(xi),   g in (0,1]
```

报出的 `f` **必须**跟着换（`:272`–`:282`，使用点 `:370`）：留着 Churchill 的数，日志上那一行的 `f` 就再也不能复现同一行的 `dp`，而"重复计费"式的自我矛盾正是 5.43–5.44 学过的形状。换来之后 `f · Re · g(xi) == 64` 不再是一句推理，而是"打印的系数就是解里的系数"。算例 `steadyBinghamYield` 是 `mu = 0.1`、`rho = 850`、`d = 0.02`、`L = 20`、两端定压 `dp = 1.0e4`，于是 `tw = 2.5 Pa` **与 `tau_0` 无关**，`xi = tau_0/2.5` 扫出 0.2 / 0.4 / 0.6 / 0.8 / 0.96 五根流动的加一根卡住的（`tau_0 = 2.6 → dpStatic = 1.04e4 > dp`）；另外四根管把这条律放到**网络 Newton** 上：`Sr1a/Sr1b` 是两条全同屈服管的串接（中点 C 是自由节点，`tau_0 = 1.0`），`Sr2a/Sr2b` 是一条卡住的上游串一条牛顿的（节点 D 自由）。

声明面：`NetworkPipe::yieldStress`（`src/1d/include/cf/1d/Network.h:99`，缺省 0 = 无阈值），解析与拒收在 `src/1d/src/network.cpp:1107`–`:1117`，内层括号解在 `:131`–`:160`，两个使用点是求解用的 `pipeLossAt`（`:252`–`:256`）与**打印**用的 `pipeFrictionReport`（`:272`–`:282`，调用点 `:370`）。判据是五个 case 文件 + checker（`case.dict` 12 根管、`case_driven.dict` 反过来给定速率、两个拒收文件），**26** 条全绿（RC=0）。

* **落在解上，不是落在容差上**：五根流动屈服管对 Buckingham 闭式最差 **4.920e-12** 相对；串接那一对在**半压**（`dp = 5.0e3`，而这个 `dp` 本身是 Newton 解出来的未知数）上对同一条闭式差 **3.323e-13**；两根无阈值管仍是 Hagen-Poiseuille 到 **1.844e-12**，也就是上面每个 `g(xi)` 比值所除的那条基线。
* **死分支报的是 0，不是"很小"**：卡住那根管 `Q / v / Re / f` 四个字段**打印就是 0**，同时 `dp` 仍是**被问到的** 1.0e4 Pa 而不是阈值 1.04e4 Pa——静摩擦吃掉要它吃的那笔，而不是被封顶的动摩擦。六根屈服管上 `Q > 0` 与 `dp > dpStatic` **同真同假**（阈值 2000 / 4000 / 6000 / 8000 / 9600 / 10400 Pa）。
* **近阈值没有被扫到左端点**：`xi = 0.96` 那根管仍然带着 **6.11731e-08 m³/s**——`Q` 随 `tau_0` 单调下降且这一段是 `Q ∝ (dp - dpStatic)²` 的软着陆，内层解如果松在括号左端，这里会先看起来"对"再在别处漏。
* **系数是解自己含的**：`f · Re · g(xi)` 在五根流动屈服管上离 64 最差 **1.843e-10**（打印精度就是这个量级），而 `dp == f (L/d) rho v²/2` 差 **2.904e-12**、`8 tw/(rho v²)` 复现报出的 `f` 差 **2.904e-12**——同一行日志里的系数与压差是一句话，不是两句碰巧接近。
* **阈值只删流量，不改壁应力**：八根并联管的 `tw = dp·d/(4L)` **逐位**都是 2.500 Pa（`tau_0` 不在式子里），而 `Q / Q_Hagen-Poiseuille` 逐根**就是** `g(xi)`：0.7339 / 0.4752 / 0.2432 / 0.0699 / 0.0031，最差 **4.920e-12**——屈服不弯抛物线，它挖掉抛物线的芯，挖掉的份额只是 `xi` 的函数；全部比值落在 `(0,1]` 里，这条就是"阈值只能把流量拿走"。
* **哨兵是逐位的**：`yieldStress 0` 对不写这个关键字在五个字段上 `==`——`tau_0 → 0` 让 `F(tw) = tw` 逐项塌回牛顿斜率，所以这一步对 100 条既有测试**结构上不可能**有影响；`f` 报告走的是同一条 `pipeFrictionReport`，未声明时那一支就是原来那一行 `churchillFriction(Re, eps/d)`。
* **死分支穿过网络 Newton 之后还站得住**：`Sr1a/Sr1b` 把自由中点落在**恰好** 2.5e4 Pa（`iterations=1`、残差恰为 0）且两根的速率到最后一个打印位相同；`Sr2a` 卡住时整条线清空——D 被下游定压拽到**恰好** 2.0e4 Pa，两根都打印 0，而 `Sr2a` 自己扛着全部 1.0e4 Pa。这一步不需要新的 Jacobian：闭式给出 `dQ/d(dp) = k[1 - (dpStatic/dp)^4] <= k`，所以 `solvePipeFlow` 里那个 Hagen-Poiseuille 上限对屈服管仍是**上界**，不会截掉真斜率。
* **反问题与偶性**：`case_driven.dict` 把速率钉成 `5.0e-6 m³/s` 让压差成为未知数，解出的 `dp = 7692.34 Pa` 对 checker 自己那套**二分**反演差 **5.326e-13**；反向取流那根打印 `Q = -5.0e-6`、`dp = -7692.34`，逐位是正向的负数——死分支判的是被选中的那个方向的符号，不是正向运气好。
* **三条墙第一次是判据**：PVT / 黑油管网带上屈服被**拒收**（屈服解是常物性单相层流律算的，逐步 `march` 那一侧没有它）；耦合元带上屈服被**拒收**（它的损失**就是**伙伴模型的陈述，这里没有墙留给阈值去 hold，5.43 为形状账立的不重复计费规则在这里落地为一条解析期拒绝）；负 `tau_0` 被**拒收**（`g(xi) > 1` 是这条律唯一能给出比牛顿更小压降的方式，那是符号写错而不是减阻剂）。第四条墙——`yieldStress` 与 `flowIndex` 同时声明——当时也是拒收（`n != 1` 时括号里的幂次换了，需要自己的单调形），5.47 落的就是这一条：同一行声明现在**解得出**，本案 checker 从"断言非零退出码"改成"断言解得出且 `dp` 仍是被问到的那一个"。

**成本**：五个 case 文件、26 条判据、墙钟 **0.28 s**（`case.dict` 12 根管 1 步收敛，`case_driven.dict` 6 步）。

**本增量的边界**：① **Bingham 而非 Herschel-Bulkley**：`tau = tau_0 + K gdot^n` 的 `n != 1` 会把 `F` 从"tw 与 1/tw³"换成一条不同的单调形，本案在解析期就拒收这种组合，不做静默近似——**这一条已被 5.47 解除**：`case_rejectRheology.dict` 改名为 `case_rheologyNowSolves.dict`，同一行声明现在解出 **2.60124e-07 m³/s**，本增量的 checker 只断言"它解得出、`dp` 仍是被问到的那一个"，数值本身归 5.47。② **层流侧 only**：`tw` 解的是 Buckingham 的层流解；屈服 + 湍流的工程关联式（Darby-Metzner 一类）没有落地，也没有被断言——本案最大的雷诺数是 **10.6**（那根未声明的牛顿管），屈服管里最大 **7.80**，全部离换段很远。③ `xi = 1` 那一格（`dp` 恰等于阈值）是浮点刀刃：`4·2.5·20/0.02` 与 `1.0e4` 的相等由 `0.02` 的表示决定，所以本案把 `tau_0` 放在阈值两侧各留 400 Pa，**不断言恰好压线**的那根管。④ 只接在常物性单相稳态律上：PVT / 黑油与两流体分布式内核**拒收**屈服声明，瞬态 1D 那一行当时仍是 Newtonian 摩阻——**幂律那一半已由 5.48 解除**（那一行现在调用同一个 `pipeFrictionPa`），**屈服那一半已由 5.49 解除**（静止支：`transient.cpp:40`–`:44` 入口闩锁 + `:353`–`:362` 穿过钳制，阈值正是这里这个 `4 tau_0 L/d`；5.49 的跨内核对照里，同一拓扑的稳态括号与瞬态静止支对两根卡住的对都报恰好 0，这是同一句话的两种算法），1D↔3D 耦合边的损失属于伙伴模型（且带屈服**仍**被拒收，`transient.cpp:226`–`:229`）。⑤ 阀的 `K` 项在 `pipeLossAt` 里加在 `4L tw/d` 之外，因此 `q = 0` 时不改变阈值（本案没有阀，这个组合未单独验证）；`roughness` 与屈服完全不进同一条式子。⑥ 圆管 only：环空/矩形的水力直径与 Dodge-Metzner 修正不在内。⑦ **全卡住的管网不在内**，而且它不是"收敛慢"而是**没有稳态解**：若每条边都停在阈值之下，注入的质量无处可去，而 `solvePipeFlow` 在退化差商上自动退回 Hagen-Poiseuille 上限（`g > 0` 不成立时取 `lam`），于是 Newton 手里是一个非奇异的 Jacobian 配一条不可能满足的连续方程——它只会撞到 `maxIter` 并报 `NOT converged`。本案把卡住那两根的两端都钉成定压节点（`Stk`）或让下游定压节点吸收（`Sr2a/Sr2b`），正是为了让这个可能性落在判据之外；真正要接的是"闸门全关 + 源还在注入"那一类瞬态/操作问题。

### 5.47 5.46 拒收的那一支第一次有自己的律：一个**剖面积分**同时装着 5.45 与 5.46——`n = 1` 时逐位是 Buckingham，`tau_0 = 0` 时逐位是幂律，而第二个极限顺手量出 5.45 的层流常数写错了 **25%**

5.46 把"屈服不是又一次代入"落地成 Bingham，代价是把 `n != 1` 那一格关在解析期外面。这一增量只补那一格，但补法是**把两格合成一个式子**：Herschel-Bulkley 的圆管层流流量本来就是剖面 `u(r) = int_r^R ((tau - tau_0)/K)^(1/n) dr` 的积分，而它逐项可积——记 `tw = dp d/(4L)`、`D = tw - tau_0`，

```
Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n) [ tau0^2/(n+1) + 2 tau0 D/(2n+1) + D^2/(3n+1) ]

n = 1 且 K = mu   ->  括号乘 D^2 就是 5.46 的 [1 - (4/3)xi + (1/3)xi^4] × tw^4/4
D -> 0            ->  Q -> 0，而阈值仍是 dp_static = 4 tau_0 L / d：**与 n 无关**
tau0 = 0          ->  只剩最后一项，Q = pi R^3 n/(3n+1) (tw/K)^(1/n)：5.45 的精确幂律排放
```

三条都是**同一行代数**的极限，不是三个分支，所以判据可以全部打在闭式上：`n = 1` 那根管必须和 5.46 的 Buckingham 支**逐位**给出同一个速率（两条代码路径——三次式的 Newton 与单调律的二分——在打印位上会合），`tau_0 -> 0` 那根管必须和 5.45 的 M-R 支会合。第二个会合是这一格真正严肃的地方：**它不成立**。5.45 把层流常数的 `K'` 写成 `K(3n+1)/(4n)`——那是湍流塌缩到层流时的教科书约定，指数是 1——而剖面积分要的是 Rabinowitsch 的壁面剪切率 `gdot_w = (3n+1)/(4n) · 8u/d`，于是 `K' = K ((3n+1)/(4n))^n`，指数是 `n`。两者之差是 `(4n/(3n+1))^(1-1/n)`：`n = 1` 处恰好为 1（这就是 5.45 自己的 `PowN1 == Newt` 逐位对照**看不见**它的原因），`n = 0.5` 处是 **1.25**——25% 的流量、20% 的压差。内核改用带指数的常数后，两条支路在 `tau_0 = 0` 上会合，5.45 那句"不是关联式而是精确解"才第一次成立。

反解（网络要的 `loss(q)` 方向）不需要导数：`Q(tw)` 对 `tw` 单调（屈服只从剖面的慢端删流体，所以壁应力越高扛得越多），左端就是阈值本身，右端从 `tau_0 = 0` 的幂律应力起——那是根的一个下界——不够就加倍，因为速率按 `tw^(1/n)` 增长，总会越过需求；二分在括号宽度 2ε 处退出。

声明面：`NetworkPipe::flowIndex` / `consistency` 在 `src/1d/include/cf/1d/Network.h:73`–`:74`、`yieldStress` 在 `:99`，两者同时声明的含义写在 `:88`–`:98`；解析与六条拒收在 `src/1d/src/network.cpp:1104`–`:1132`（三条属屈服 `:1107`–`:1117`、三条属指数 `:1118`–`:1132`，**其中"指数 + 阈值 = 拒收"那一条已删除**）；正定律 `herschelBulkleyRate` 在 `:173`–`:185`、反解 `herschelBulkleyWallStress` 在 `:196`–`:221`、两支路的选择器 `yieldedWallStress` 在 `:227`–`:232`；广义数 `generalizedReynolds` 在 `:101`–`:117`（HB 那一支 `:104`–`:109`），两个使用点仍是求解用的 `pipeLossAt`（`:252`）与**打印**用的 `pipeFrictionReport`（`:272`–`:282`），死分支仍是 `equipmentFlow` 里括号之前的那一段（`:323`，条件只有 `yieldStress > 0`，与 `n` 无关）。算例 `steadyHerschelBulkley` 十三根管：`rho = 850`、`mu = 0.1`，十一根并联的两端都钉在定压节点（A 3e4 / B 2e4，`case_dp2.dict` 把 A 抬到 4e4），所以每根恰好承担 `dp = 1.0e4 / 2.0e4`；两根串接经过自由节点 C。`n` 扫 0.5 / 0.7 / 1.0 / 1.3，`tau_0` 扫 0（未声明）/ 1e-9 / 0.5 / 1.0 / 2.0 / 2.4 / 2.6，其中 `2.6 → dp_static = 1.04e4` 是卡住的那根。判据是三个 case 文件 + checker，**22** 条全绿（RC=0）。

* **两条极限都在打印位上会合**：`HbN1`（`n = 1, K = mu`）与 `Bing`（不写指数）在 `Q / v / f` 三个字段上**逐位相等**（`==`，不是接近）——三次式 Newton 与剖面积分二分是两条路径，会合才说明新支是 5.46 的**推广**而不是对手；`tau_0 = 0` 那根管对幂律剖面积分差 **7.023e-13**，这一条是改常数**之后**才成立的。
* **落不进极限的那个差被算出来了**：同一个 `n = 0.5` 声明成 HB 与声明成幂律，`tau_0 = 1e-9` 的扰动让两支差 **1.001e-09**，而它声明的阈值本身是 `tau_0/tw = 4.0e-10`——差的是被声明进去的那个数，不是两支在分歧。
* **六个陡状态 + 串接那一对落在闭式上**：并联六根（其中一根就是死分支，两边都恰好是 0）最差 **2.331e-12** 相对（`dp = 1.0e4`）/ **1.835e-12**（`2.0e4`），串接那一对在半压 **5.0e3**（这个 `dp` 是 Newton 解出来的未知数）差 **3.634e-13**、翻倍后 **6.313e-13**。
* **同一积分再算一遍，不逐项**：把 `int tau^2 gdot(tau) d tau` 直接数值积分（`s^(1/n)` 在屈服半径处导数无穷，所以 [0,1] 按 2 的幂向 0 加密成 21 个面板 × 12 节点高斯-勒让德），对内核仍是 **1.511e-09** 相对；而**不加密**的单个面板差 **3.296e-06**、加密一倍再加密差 **0.000e+00**——比的是收敛了的积分。这一条和闭式是同一个陈述的两次独立推导，`(3n+1)` 写错一处不会在另一处得到同样的错。
* **软着陆那两格的 1e-9 是求解器的，不是律的**：`HbThin2`（`D` 只有阈值上方的 0.5 Pa）与 `HbNear`（差 0.1 Pa）落在 **1.511e-09** / **1.078e-10**；把打印的速率再喂回同一条闭式反演，内核留下的压差分别是 **1.112e-06 Pa** 与 **1.822e-08 Pa**，乘上各自的斜率 `d ln Q / d ln tw = 13.593` 与 `59.151` 恰好复现 **1.511e-09** 与 **1.078e-10**（预测与观测差在第 4 位）。`equipmentFlow` 在速率括号宽到 **1e-9** 相对时退出（`src/1d/src/network.cpp:364`），而压差是算例钉住的那个量，所以留在压差里的余量要乘上着陆自己的斜率才变成速率的差。断言因此分组：陡的那六根打到 **1e-11**，软的那两根打到 **1e-8** 并**说明**差额从哪来——不是把律放宽。
* **阈值的几何性与 `n` 无关，`D` 的形状有关**：`4 tau_0 L/d` 是 2000 / 8000 / 10400 Pa，`tau_0` 从 0.5 抬到 2.0 让速率从 **2.22509e-06** 掉到 **4.73333e-08**；卡住那根在 1e4 Pa 下 `Q / v / Re / f` **打印就是 0** 而 `dp` 仍是被问到的 1e4，同一根管在 2e4 下流出 **6.04532e-07**（闭式差 **9.720e-14**）——5.46 那条"死分支由被问到的压差决定，不由管的属性决定"在 `n != 1` 上重演了一遍，而 `D` 的主导幂次是 `1 + 1/n`：Bingham 是平方，`n = 0.5` 是**立方**。
* **串接与对称**：自由节点 C 落在**恰好** 2.5e4 Pa（翻倍后 3.0e4），两根 halves 五个字段**逐位相同**；驱动差翻倍时这一对从 **1.17871e-08** 涨到 **3.58889e-07**，即每半边的 `D` 从阈值上方 0.25 Pa 走到 1.25 Pa——软着陆是在**两个点**上被读，不是在一个点上被断言。
* **报出的 Re 换成表观数，并且这条不再是证据**：带指数的屈服管没有材料黏度可报，所以 `Re = rho v d / mu_app = 8 rho v^2 / tw`，配上 `f = 8 tw/(rho v^2)` 之后 `f · Re == 64` 是**定义**（最差 **2.911e-10**，即打印精度），断言它的理由只是"两个数不同步说明日志在印另一个状态"。有物理内容的是两条：`HbN1` 与 `Bing` 的 Re 之比**就是** `g(xi) = 0.475200`（**2.39928** 对 **5.049**，黏度之比即屈服挖掉的速率份额），以及 `dp == f (L/d) rho v^2/2` 在九根流动屈服管上差 **1.029e-10**。
* **整案停在层流**：两个文件里最大的雷诺数是 **84.10**（翻倍驱动差下的宽孔径管，表观数），那根未声明的牛顿管是 **21.25**——离换段两个量级，所以"湍流侧的屈服"仍不是这里能断言的东西。
* **拒收的账**：删掉一条（指数 + 阈值，改由 `case_rheologyNowSolves.dict` 断言它解得出），留下五条，其中**新写的一条**是"声明指数与阈值而不给 `consistency`"（`case_rejectNoConsistency.dict`）：`K` 的单位是 `Pa.s^n`，不能被 `mu` 静默默认。耦合元带屈服 / 带指数、PVT 与黑油带流变、负 `tau_0`、负 `n` 全部照旧。

**成本**：三个 case 文件、22 条判据、墙钟 **0.21 s**。改常数同时动了 5.45 的判据：它的 checker 把 `K'` 与闭式一起换掉之后 **16** 条仍全绿，速率对剖面闭式 **3.249e-12**、`f · Re_MR == 64` **1.563e-10**、`v == Q/A` **4.376e-12**（打印位）。

**本增量的边界**：① **层流 only**，而且 `pipeFrictionReport` 在 `yieldStress > 0` 时**永远**报解自己含的那个数、不看 `Re` 多大，所以屈服 + 湍流（Darby-Metzner 一类工程关联）在这里既没落地也没被证伪——5.46 边界② 的那一格，现在 HB 也共享它。② **表观 Re 是约定，不是测量**：`f · Re == 64` 在这一支是恒等式，跨两支比较雷诺数必须换算黏度（`HbN1` 与 `Bing` 之比 **0.475200** 就是那次换算），不能把两个数当同一个数读。③ **外层求解的退出精度第一次被量出来而没有变成判据**：`network.cpp:347` 的残差门（`1e-13 × dp`）与 `:364` 的速率括号门（`1e-9` 相对）在软着陆上会先碰后者；本案记录了哪几格、余量多少、乘数是什么，但**没有**给出"多少位以内算合格"的标准，也没有改那两个门（它们是全部稳态管网测试共用的）。④ 只接在常物性单相稳态律上：PVT / 黑油与两流体分布式内核仍在解析期**拒收**流变与屈服（`:1110`、`:1125`），瞬态 1D 那一行当时仍是 Newtonian 摩阻——**幂律那一半已由 5.48 解除**（那一行现在调用同一个 `pipeFrictionPa`），**HB / 屈服那一半已由 5.49 解除**（静止支读的正是这里这条剖面积分：`ReleaseHB05` 落在 **2.2250853568e-06**、对本案同一个闭式差 **1.138e-12**，而阈值本身不认识 `n` 这件事由 `StuckHB` 断言——它与未声明指数的 `Stuck` 报同一个 10400 Pa），1D↔3D 耦合边的损失属于伙伴模型。⑤ 极限的两端没扫：`n -> 0`（近似完美塑性）、`n` 极大、`D` 极小时的括号宽度与加倍次数上限（`200` 次），`flowIndex 0` 仍是"未声明"的哨兵而不是一个可声明的极限；`tau_0` 与 `K` 使 HB 在 `n = 1` 处等价于 5.46，但**没有**任何一条断言说这条律在湍流侧退化成 5.45 的 M-R 惯用式。⑥ 圆管 only（矩形/环空的水力直径与 Dodge-Metzner 修正不在内）；阀的 `K` 项与 `roughness` 都不和屈服同台。

### 5.48 稳态那条验过三次的壁面律第一次被**瞬态那一行**读到：一个共享函数换掉手工写的 `rho u d / mu`，于是声明 `flowIndex` 不再是被忽略而是被读——定态落点是同一个打印字符串，整条 2001 步轨迹由律自己复现，阻尼走律的斜率 `n·dp/Q`（弦与 Newtonian 下限各差 23%~50%），而屈服那一支当时由四条拒收断言（第四条已在 5.49 变成被解出的物理）

5.45、5.46、5.47 三格证的都在稳态管网里，三格的边界句子里都有同一句"瞬态 1D 那一行仍是 Newtonian 摩阻"。这一增量只做一件事：**把那一行的摩阻换成稳态用的同一个函数**。`pipeFrictionPa`（`src/1d/src/network.cpp:390`，声明 `src/1d/include/cf/1d/Network.h:202`）在 `cf` 作用域而不是上面那个匿名作用域——瞬态内核是另一个翻译单元，放进匿名作用域它就只能被 `pipeLossAt` 看见。稳态侧 `pipeLossAt`（`:254`）仍是在它之上只加阀项 `K·rho·v²/2`；瞬态侧 `frictionDp`（`src/1d/src/transient.cpp:72`–`:80`）现在除了调用与一次奇延拓（`std::abs(q)` 进、符号回）**什么都不算**，而这正是重点：它之前自己写 `re = rho u d / mu`，所以一份声明 `flowIndex 0.5` 的瞬态算例既不报错也不警告，只是把声明**忽略**掉然后报出 Newtonian 的速率——最坏的一种静默数字。阻尼那一侧同时换了：`frictionSlope`（`:85`–`:105`）在声明指数时取律自己的中心差分并在 `:106`–`:114` 早退，Newtonian 下限 `0.2·128 mu L/(pi d^4)` 只留在未声明的路径上（`:115`–`:116`）。行本身一字未动：两端定压节点时每根管恰好是一行方程 `rho L/A·(Q − Q_prev)/dt = dp − dp_wall(Q)`，`dp` 由节点压力给定。

声明面：`TransientPipeConfig::flowIndex` / `consistency` / `yieldStress` 在 `src/1d/include/cf/1d/Transient.h:75`–`:77`（`:68`–`:74` 的注释当时写明"`yieldStress` 只被读出来是为了拒收"，5.49 把它改成静止支的阈值），解析在 `transient.cpp:743`–`:745`，赋给几何在 `:208`–`:210`，回读接口 `pipeFlowIndex` / `pipeConsistencyOf` 在 `:556` / `:563`（声明 `Transient.h:181`–`:182`）。app 只在**声明了指数**的管上打印 `pipe X wall-law flowIndex n consistency K`（`src/1d/apps/transientApp.cpp:42`–`:54`），所以"未声明的案子打印和以前逐字相同"是一条能在输出字符串上失败的断言，而不只是数字上的。算例 `transientPowerLawRheology`：A 3e4 / B 2e4，六根管**同一几何**（`L = 20 m`、`d = 0.02 m`）、**同一压差** 1.0e4 Pa，流体 850 kg/m³ / 0.1 Pa·s，`dt 0.002`、`endTime 4.0`（2001 个状态；`dt` 必须小 against 最快那条弛豫，`n = 0.7` 的 `tau` 只有 7.5 ms，本案因此走过它 37 个时间常数）。五份声明：未声明、`n = 1, K = mu`、`0.5 / 1.0`、`0.7 / 2.0`、`1.3 / 0.02`，外加 `Thin05Rev` 以 `initialFlow -2.0e-6` 起步。判据是 checker 从 case 文件把律、剖面积分闭式与隐式欧拉递推**各自重新写一遍**再逐步对照，不比对记忆里的数字：**45** 条全绿（RC=0）。

* **一条律、两个内核、同一个打印字符串**：`t = 4` 定态时瞬态 app 报的 `Q` 与稳态 app 同一声明下报的 `Q`，`Newt`/`Decl1` **1.96349540849e-05**、`Thin05` 与 `Thin05Rev` **3.92699081699e-06**、`Thin07` **9.75728733974e-07** 三格是逐位相同的字符串；只有 `Thick13` 差在第 10 位（稳态 **3.41904523976e-05** 对瞬态 **3.41904523728e-05**，相对 **7.26e-10**）。
* **对剖面积分闭式**（5.47 修正 5.45 之后那一条，`Q = (pi d^3/8)·n/(3n+1)·(dp·d/(4 K L))^(1/n)`）：五根并联管最差 **7.256e-10**（`Thick13`），其余 **1.844e-12**（`Newt`/`Decl1`）、**7.023e-13**（`Thin05`、`Thin05Rev`）、**3.885e-13**（`Thin07`）——稳态侧 5.47 量到的是同一个数。
* **那 7.256e-10 是求解器的停，不是律的**：把闭式的速率再喂回内核自己的 `pipeFrictionPa`，被问的 1.0e4 Pa 回到 **3.638e-16**——Churchill 的层流半区与剖面积分在这里是同一句话。余下的差额是每步 Newton 的退出（`src/1d/src/transient.cpp:418`，`maxRel < 1e-11`，残差按 `rho L/(A dt)` 定标，`:389`–`:391`）在整条弛豫里往下沉时的累加，不是定律的误差。这一条不是写在句子里的解释，而是**另一条断言**（根差 < 1e-12），所以"是求解器不是律"这个说法本身可以失败。
* **整条轨迹都是那一行的代数**：每根管 2001 个状态由 checker 自己解隐式欧拉递推逐步对照，最差 **7.246e-10**（`Thick13` 第 2000 步），`Thin05Rev` **6.128e-10**、`Thin05` **2.137e-11**、`Newt`/`Decl1` **4.595e-12**、`Thin07` **8.109e-13**——不是"形状碰巧落在对的终点"。
* **`n = 1, K = mu` 是惰性的**：与未声明那根管 **2001/2001** 个状态的打印**逐字**相同（字符串相等，不是数值接近）——逐项同一条代数，所以任何一处不同都是漏。这也是 5.45 当年 `PowN1 == Newt` 那条对照在瞬态侧的对应物。
* **指数不是装饰**：同一压差同一几何，定态速率是 Newtonian 的 **0.2000**（`n = 0.5`）、**0.0497**（`n = 0.7`）、**1.7413**（`n = 1.3`）倍，最邻近的一根与 Newtonian 也差 **74.1%**——旧行为会把这六条轨迹压成同一条，所以这一格是**一次运行就能证伪**的。
* **阻尼走律的斜率，不走弦也不走下限**：从尾部拟合 `s = (I/dt)(1/r − 1)`（`r` 是拟合窗口两端偏差的几何比）得 `Newt` **5.092958e+08**（`tau = I/s` = 0.1063 s）、`Thin05` **1.275283e+09**（0.0424 s）、`Thin07` **7.179707e+09**（**0.0075 s = 3.8 步**）、`Thick13` **3.798422e+08**（0.1425 s），对律自己的 `n·dp/Q` 差 **9.569e-09 ~ 1.603e-03**；而弦 `dp/Q` 分别差 **49.9% / 29.9% / 23.0%**，Newtonian 下限 **1.018592e+08** 对每一根都差 10% 以上——两个"看起来能用"的候选都被量出来"不是行实际用的那个"。
* **奇延拓只有一行 marching 的话说得出来**：`Thin05Rev` 从 −2.0e-6 起步，全程**恰好一次**换号（`dp_wall(−q) = −dp_wall(q)`）、无回弹、无过冲（没有任何报告状态超过律自己的定态速率），终点与同一条律从静止起步那根管逐位相同（**3.92699081699e-06**）——定态属于律，不属于初始条件。稳态那三格说不出这句话，因为它没有初始条件。
* **仍整案停在层流**：2001 × 6 个状态里最大的广义雷诺数 **32.2165**（`Thick13`），Newtonian 那根 **10.625**，离换段两个量级，所以上面那条闭式是**精确解**而不是关联式；湍流侧的广义数被读到了但没有在这里被验证。
* **拒收的账第一次是本增量的边界本身**：四条（`case_rejectNegativeIndex` / `rejectNoConsistency` / `rejectCoupledRheology` / `rejectYieldStress.dict`，判据在 `transient.cpp:212`–`:229`）——负指数不是律（0 才是"未声明"的哨兵）、`n != 1` 时流体的 `mu` 不是 `K`（`Pa.s^n`）这个量、耦合元的损失已经是伙伴模型的陈述所以不能带流变、以及**任何** `yieldStress`：瞬态行只有**动壁律**，阈值需要一个自己在时间里积分的静止支，而这里 `dp` 正是行在解的未知数（稳态 5.46 那次便宜，是因为问到的压差已知）。**第四条已被 5.49 解除**：`case_rejectYieldStress.dict` 换成了 `case_yieldNowSolves.dict`（`Bing` 只有阈值、`Hb` 阈值 + `n = 0.5, K = 1`、`Stuck` 阈值 2.6 即壁要 1.04e4 Pa），本案第 10 族现在断言它们**解得出**并且定态落在 5.46/5.47 的屈服闭式上（`Bing` **1.4409438297e-05** 对 **1.4409438304e-05**，最差 **5.528e-10**；`Hb` **2.2250853568e-06**，**1.138e-12**；`Stuck` 的闭式自己的死分支说 0，轨迹也说 0），而"只有 `Hb` 报出一条 wall law"仍然是声明面那条断言。剩下三条照旧。

**成本**：5 个 case 文件、45 条判据、墙钟 **6.64 s**（ctest -j8 记，Test #77；5.48 落地那一轮同一条测试记的是 2.58 s，差额是 5.49 把第四条拒收换成一个被解的算例时加进去的那一族）。其余既有稳态与瞬态测试未改一字通过：`steadyPowerLawRheology` / `steadyBinghamYield` / `steadyHerschelBulkley` 三格的判据打的仍是 `pipeLossAt`，而它现在只是"共享函数 + 阀项"。

**本增量的边界**：① **幂律 only**，屈服在此格曾被拒（见上第四条）——**这一条已被 5.49 解除**：同一行的 `yieldStress` 现在有自己的静止支，`case_rejectYieldStress.dict` 改名成 `case_yieldNowSolves.dict` 并被断言**解得出**，本案 checker 的第 9 族因此从四条拒收变成三条，另加第 10 族测三根屈服管的定态落点（最差 **5.528e-10**，其中 `Stuck` 打印的是 `Q 0`）。② **层流 only**：湍流侧走广义数（`generalizedReynolds`）但没有屈服 + 湍流的判据（Darby-Metzner 一类仍是 5.46 边界② 那一格）。③ 每步 Newton 的退出精度第一次量到它**在定态里留下多少**（7.256e-10），但**没有**改那个门（`transient.cpp:418` 是全部瞬态测试共用的），也没有给出"多少位以内算合格"的标准。④ 只接常物性单相瞬态行：PVT / 黑油与两流体分布式内核仍在解析期**拒收**流变与屈服。⑤ 耦合边不带流变（`transient.cpp:215`），所以 1D↔3D 详细区里的非牛顿伙伴模型仍是空白。⑥ 圆管 only（环空的 Dodge-Metzner 不在内）；阀 `K` 项与 `roughness` 都不和流变同台。⑦ 瞬态响应只与"同一律的隐式欧拉递推"对照，**没有**和任何解析的非定常解（Womersley 一类）对照——本案证的是行读的是哪条律，不是时间离散本身的阶。

### 5.49 屈服壁第一次在**瞬态那一行**有自己的静止支：奇延拓的律在 `|drive| <= 4 tau_0 L/d` 里根本没有根，所以"把行改成 `Q = 0`"不是可接受的近似而是唯一能步进的写法——四根管各自因不同理由打印同一个 `Q 0`、停止**穿过**零并且与递推同步落在第 17 步、闩锁由网络自己在第 253 步放开，而放开那一对的终点逐字就是 5.46 稳态内核算出的速率（差 **1.864e-08**）

5.48 把幂律交给瞬态那一行之后留下的正是这一格：屈服管的壁压降不随速率消失，`dp_wall(Q) -> dp_static = 4 tau_0 L/d`（`Q -> 0+`），而在 `|drive| <= dp_static` 之下**那条律没有正根**。稳态不必回答这件事——它问到的压差是已知的，括号可以直接停在 `Q = 0`（5.46）；瞬态的 `dp` 恰是行在解的未知数，所以它需要的是一条**在时间里积分**的支。这一格把它落成三件东西，缺一件都不是同一句话：(i) **步入口闩锁** `holdsStatic`（`src/1d/src/transient.cpp:40`–`:44`：`qPrev == 0 && |drive| <= dp_static`）在 `:310`–`:316` 的循环里对每根管**一次**决定并写进 `PipeDyn::held`（`src/1d/include/cf/1d/Transient.h:221`）——之所以要闩而不是每趟都问：被 hold 的管看到的 drive 依赖同一次求解正在移动的节点压力，不闩就会在一个时间步内让该行在 `Q = 0` 与动壁律之间来回翻（注释 `:342`–`:346`）；(ii) **求解内的穿过钳制**（`:353`–`:362`）：`yieldStress > 0 && qPrev != 0 && |drive| <= dp_static && x[k]*qPrev <= 0` 时该行写成 `x[k] = 0; r[k] = 0; J[k][k] = 1`，这是有限时间停住的那一半，分辨率是 dt；(iii) **静止处的摩阻值与斜率**：`frictionDp`（`:72`–`:81`）在 `q == 0` 时返回 `copysign(dp_static, drive)`（`drive == 0` 才是 0），而 `frictionSlope`（`:85`–`:117`）在屈服支不走中心差分，走"超额割线" `(dp_wall(|q|) - dp_static)/|q|`（`:89`–`:103`）——屈服律在静止点的斜率是**无穷**，跨 `q = 0` 差分把阈值跳变本身除进商里（`dp_wall(h) - dp_wall(-h) = 2 dp_wall(h)`，斜率量级 1e15），Newton 于是永远坐在静止点上，那正是与 5.48 拒收方向相反的那个错。drive 由 `drivePa`（`:52`–`:56`）给出 `pFrom - pTo - rho g dz`，与行自己的 `dp - grav`（`:352`）同一条约定，所以闩锁与行不可能对"哪端是上游"各有说法；阀项与对流项不在 drive 里，因为两者都随速率一起消失。阈值 `yieldThresholdPa`（`:28`–`:31`）只吃几何，不认识 `n`——5.47 已经把它写成极限，所以这里直接沿用。

声明面：`TransientPipeConfig::yieldStress` 在 `Transient.h:77`（`:68`–`:74` 的注释现在说的是"静止支顶着的那个阈值"而不是"只被读出来是为了拒收"），解析 `transient.cpp:745`，赋给几何 `:210`，回读接口 `pipeYieldThreshold` / `pipeHeld` 在 `:542` / `:549`（声明 `Transient.h:192`–`:193`）。app 只对声明了阈值的管打印 `pipe X static-hold threshold 4 tau_0 L/d`（`src/1d/apps/transientApp.cpp:55`–`:64`），并在**每个样本**打印 `pipe X static hold 0|1`（`:84`–`:92`；t = 0 还没解过任何步，所以第一格报 0），于是"这根管此刻是被 hold 的"是一条能在打印字符串上失败的断言。算例 `transientYieldStaticHold/case.dict`：两端定压 A 3e4 / B 2e4（外加 A2 2.4e4），十根管**同一几何**（`L = 20 m`、`d = 0.02 m`）、流体 850 kg/m³ / 0.1 Pa·s、`dt 0.002`、`endTime 4.0`（2001 个状态），`rho L/A = 5.411e7`；阈值分别是 `tau_0` 0.5 / 1.0 / 1.5 / 2.6 → 2000 / 4000 / 6000 / 10400 Pa。另有 `case_series.dict`（同一 stuck 管放在**自由**节点的上游与下游各一次，`C`/`D` 体积 1.0、`F` 体积 0.1，故意让 `F` 的弛豫时间是 `C`/`D` 的十分之一）与 `case_steadySeries.dict`（同拓扑、`solver { type steadyNetwork; }`），判据仍由 checker 把律、阈值、递推与节点体积账**各自重写一遍**再逐步对照：**47** 条全绿（RC=0）。

* **四根管打印同一个 token `Q 0`，而理由各不相同**：`Stuck`（壁要 1.04e4 Pa，桌上只有 1.0e4）、`StuckHB`（同一个阈值配上 `n = 0.7`，所以阈值**不认识**流变指数）、`Hang`（升 1 m 之后 drive 只剩 1664.35 Pa，对 2000 Pa 的墙）、`HeldBelow`（静止在 4000 Pa drive、6000 Pa 墙下——动壁律**没有**任何速率可要）。四根都是 2001/2001 状态 `Q 0`、闩锁 2000/2000 步；checker 读的是打印的**字符串**，不是"小于某个容差"，所以 1e-19 也过不了这一关。
* **重力在 drive 里，而且它决定生死**：`Hang` 与 `Roll` 是同一根管、同一个阈值（`tau_0 = 0.5` → 2000 Pa）、同一对节点，只差 `dz` 的符号；`rho g dz = 8335.6525 Pa` 把它们分成 1664.3475 Pa（阈值之下，被 hold）与 18335.6525 Pa（阈值的 9.2 倍之上，流动）。`Roll` 落在**重力增强**那笔压差的 Buckingham 闭式上：**3.0767680512e-05** 对 **3.0767680529e-05**（相对 **5.419e-10**），对应的壁应力是 `tw = 4.5839 Pa` 而不是只看节点对会得到的 2.5 Pa。一个只比较节点压力的支会把 `Hang` 移动——它没有，因为墙抵抗的是"水头减去抬升"。
* **放开之后落在 5.46 / 5.47 的屈服闭式上**：`ReleaseB05` **1.4409438297e-05** 对 **1.4409438304e-05**（**5.250e-10**，`tw = 2.5`、`xi = 0.2`）、`ReleaseB10` 与 `ReleaseN1` **9.3305301783e-06** 对 **9.3305301812e-06**（**3.046e-10**，`xi = 0.4`）、`ReleaseHB05` **2.2250853568e-06**（**1.138e-12**）——稳态侧 5.46/5.47 量的是同一批闭式，一条律、两个内核。`n = 1, K = mu` 在屈服行上同样惰性：`ReleaseB10` 与 `ReleaseN1` **2001/2001** 个状态的打印逐字相同，所以 5.47 那条"一个闭式包含 5.45 与 5.46"在时间里也成立。
* **有限时间停住，而且它是**穿过**零停的**：`StopBelow` 以 2.0e-6 m³/s 起步，4000 Pa drive 对 6000 Pa 墙，在 `t = 0.034`（第 17 步）达到**恰好** 0，其后 **1984** 个状态打印单个 token `0`；闩锁从下一步才报（1983 步），停住那一步报的是穿过钳制——所以 flag 是**入口闩锁**，不是速率。对照量：没有静止支时这一行每步的位移是 `(dp_static - drive) dt / (rho L/A) = 7.392e-08 m3/s`，而停住前最后一个速率是 **5.409e-08 m³/s**，恰在一个这样的步长之内——两种读法之间隔的是一条支，不是一个容差。
* **整条 10 × 2001 的状态是那一行（含 `Q = 0`）的代数**：checker 自己解带支的隐式欧拉递推逐步对照，最差 **5.419e-10**（`Roll` 第 2000 步），`Stuck` / `StuckHB` / `Hang` / `HeldBelow` 是逐位 **0.000e+00**，`StopBelow` **1.428e-11**、`ReleaseHB05` **9.669e-11**、`ReleaseB05` **5.250e-10**；并且递推**自己**也在第 17 步把 `StopBelow` 停住（17 对 17），所以那个停止是支的算术，不是内核额外加的阻尼。
* **被 hold 的管把节点封住，封到的是求解器自己的退出**：`case_series.dict` 里 stuck 管的一端是**自由**节点，节点行保留它在自己那一列的系数、hold 行保留被强制的 0，于是离散体积账 `V/(K dt)·(p_i - p_{i-1}) - (Q_in - Q_out) = 0` 关到 `4.874e-13`（节点 C）与 **4.960e-13**（节点 D），对求解器停下的 `1e-11`（`:414` 的定标带一个 `+1.0 m3/s`，所以这里留下的是一条**绝对**余量而不是相对余量）。这一条既是好消息也是要说明白的：`4.960e-13 m3/s` 是本案最大节点吞吐 `1.657e-05` 的 **3.0e-08**，作为"墙漏了"它需要 **1.008e+08 s** 才能排空节点 C 那 5000 Pa，而本案 4.0 s 就走完了——所以剩下的是求解器的数字，不是物理的数字。
* **闩锁可以由网络自己放开**：`ReleaseLate`（`tau_0 = 2.0` → 8000 Pa）在 5000 Pa 的 drive 下起始被 hold，前 **253** 个状态打印恰好 0，在 `t = 0.506`（第 253 步）放开。翻的是**入口** drive 跨过阈值那一步：节点 F 那时已经掉到 **21993.1149 Pa**，留下 **8006.8851 Pa** 对 8000 Pa，而上一个样本只有 **7996.7725 Pa**。之后不回退——它排空的那个节点终点在 20450.6031 Pa，drive 9549.3969 Pa 是它自己的墙必须抵抗的。"被 hold" 是**行的状态**，不是管的属性。
* **放开那一对的终点与稳态内核同一拓扑**：瞬态 **8.8475721067e-07** 对稳态 **8.8475719418e-07**（相对 **1.864e-08**），它的节点落在 **20450.603132** 对 **20450.603139 Pa**（**3.374e-10**）；两根卡住的对在稳态内核里同样报恰好 0，这是上面"四根管"那格的跨内核一半。而两个被封的节点**还没**到它们的稳态压力（**22291.3** 对 20000.0、**25875.7** 对 30000.0 Pa），各自在自己正去的那一侧：这个 RC 网络的弛豫时间 `(V/K)·128 mu L/(pi d^4) = 5.09 s`，4.0 s 只是 **0.79** 个——所以本案与 5.46 先在**速率**上 agree，后才在浮动压力上 agree，这一点本身写成断言而不是解释。
* **拒收剩下两条**（`:222`–`:229`）：负 `tau_0`（阈值是静止的墙 hold 的应力，负的会让 held 集为空，于是每根屈服管都"付得出"任何速率）与耦合元带阈值（它的损失已经是伙伴模型的陈述，没有墙留给阈值去 hold）。5.48 那第四条——"任何 `yieldStress`"——不再是拒收，而是本节。

**成本**：5 个 case 文件（10 根管的定压案、6 根管 3 自由节点的串联案、同拓扑的稳态对照案、两份拒收）、47 条判据、墙钟 **140.14 s**（ctest -j8 记，Test #78；单独跑 122 s，而 5 次 app 调用只占 **0.31 s**——`case.dict` 0.198 s、`case_series.dict` 0.107 s、`case_steadySeries.dict` 0.007 s，所以贵的是 checker 把带支递推在 Python 里重解，不是求解器）。既有测试里只有 `transientPowerLawRheology` 那一格被改：第 9 族从四条拒收变成三条并新增第 10 族"屈服现在被解"（共 45 条，ctest 记 **6.64 s**，5.48 落地那轮同一条记 2.58 s）；`transient.unit` 与全部稳态屈服测试未改一字通过，打的仍是同一个 `pipeFrictionPa`。本轮全量回归 104/104 绿，实测 **553.98 s**。

**本增量的边界**：① 静止支是**刚**的：`|drive| <= dp_static` 时 `Q` 恰好 0，既没有阈值以下的蠕变/弹性形变，也没有任何随时间变化的结构（触变、凝胶强度、"静置 N 小时后要多少压差才能重新启动"那一类 OLGA / Pipesim 的启动压力）——本案证的是**有没有**这条支，不是它的流变学时间史。② 闩锁只在**步入口**决定，放开最早发生在下一个入口，所以开关的时间分辨率是 `dt`（本案第 253 步 = 0.506 s）；这是有意的（不闩会在一步内来回翻），代价是细于 dt 的开关看不见，而分区接口重解同一步时也不会重估（`Transient.h:217`–`:221` 那段注释说的就是这个）。③ 阈值只与 drive 比：阀 `K`、对流动量账 `beta`（§5.42–§5.43 那一格）**没有**与静止支同台，"静止的管上还挂着一个局部损失"这一格未验证。**（此处已在 5.50 更正：阀那一半由 5.50 落成算例——静止的管确实没有阀降可报，而 `drivePa` 里没有 K 是逐字符串被证的；**常数 K** 那一半由 5.50 关掉，**时间表**那一半（`valveClose` / `valveOpen` 的斜坡与静止支同台）由 5.52 关掉，结论是同一条：闩锁不看 K，所以日程改不了"停不停"，只改"第几步停"。`beta` 与静止支同台仍是空白。）** ④ **层流 only，而且本案不断言雷诺数**：与 5.48 第八族不同，本 checker 没有"全部状态落在层流半区"的判据，所以"落在层流侧"在这里是算例设计的后果而不是被证的东西；屈服 + 湍流（Darby-Metzner 一类）仍是 5.46 边界② 那一格。⑤ 跨内核终点对的是**同一拓扑的稳态内核**，不是任何实测或场数据；比较还发生在 0.79 个弛豫时间上，节点压力未收敛（上面那条 bullet 把这点写成断言，所以它同时也是这条边界的证据）。⑥ 只接常物性单相瞬态行：两流体分布式内核与 PVT / 黑油侧根本不读 `yieldStress`；耦合边带阈值被拒收（`:226`–`:229`），所以 1D↔3D 详细区里的屈服伙伴模型仍是空白。⑦ 每步 Newton 的退出门（`:418`，节点定标 `:414`）没有改：本增量第一次把节点账的**绝对**余量（4.960e-13 m3/s，来自定标里的 `+1.0`）量出来，但"多少位以内算合格"仍未立为标准——与 5.48 边界③ 同一格。⑧ "源还在注入而所有出口都卡在阈值下"这一类操作问题没有测：本案把注入端钉成定压节点，所以全 hold 的拓扑里节点压力会一路走上去这件事**没有**被断言（5.46 边界⑦ 在瞬态侧的对应物）。

### 5.50 屈服壁与局部损失第一次落在**同一条行**上：静止的管没有阀降可报（K 差五个数量级仍是同一个字符串，稳态侧连 `dpValve` 字段都不出现），动的管付两笔账（同一阈值四个 K，速率 **9.33053017832e-06 → 7.20655693864e-06**，阀吃掉 drive 的 **0 → 11.18%**，每一格都压在 checker 自己二分出的两律根上，最差 **5.419e-10**），而停住那一步由惯性对墙的余量决定、与 K 无关（第 **17 / 17 / 5** 步，`q* = 7.391983e-08` m³/s 逐条复现）——**内核一字未改**

5.46 把 `valve K` 与屈服写在同一段边界里，5.49 又把"阀与 `beta` 没有和静止支同台"列为自己的边界③，于是这一格被两边各让了一次而谁也没量。落地的方式不是加分支：两笔账早就在同一条行里（阀项 `transient.cpp:370`–`:371`、残差 `:381`、Jacobian `:383`；稳态侧 `network.cpp:454`–`:455` 由 `steadyNetworkApp.cpp:65` 打成 `dpValve`），静止支比的是 `drivePa`（`:52`–`:56` 的 `pFrom - pTo - rho g dz`，里面没有 K），而 `K rho|Q|Q/(2A²)` 随速率一起消失——这两件事合起来就是本节的全部物理，剩下的是把它们**各自**证到能被失败的粒度。于是本案把两者放到 14 根同一几何的管上（`L = 20 m`、`d = 0.02 m`、850 kg/m³ / 0.1 Pa·s、`dt 0.002`、`endTime 4.0`，2001 状态；节点 A 3e4 / A2 2.4e4 / B 2e4 全定压），一次运行同时给出：静止侧的**恒等**（不是容差）、运动侧的**根**、以及停住的**步号**。checker 把两条律都在 Python 里重写一遍，并且把每步的隐式欧拉解改写在**壁应力**上（`Q(tw)` 单调、`dp_wall = 4L tw/d` 线性、阀随 `Q` 单调，所以合成 residual 单调，一次二分即可），因此整案 14 × 2001 全轨迹复现只花几秒而不是 5.49 那种重解量级。

声明面没有任何新键：`valve` 与 `yieldStress` 各自早已解析（`:741`、`:745`），二者同台的拒收也一条没有增加——反过来说，本增量能列出的"仍被拒收"只有耦合边那两条：`requireNoValve`（`:143`，用于 `:234`–`:235` 与 `:624`）说耦合元不能带阀，`:226`–`:229` 说耦合元不能带阈值，所以详细区里"壁 + 阀"两笔账都由伙伴模型自己陈述这件事仍未被打开。app 侧只有一条可读性差异：瞬态 app 打印 `Q` 与 `static hold`（`transientApp.cpp:84`–`:92`），**不打印** `dpValve`，所以瞬态侧"阀被记账"只能由速率反证，压降的拆分只在稳态内核那一份输出里被断言。

* **静止的管不回答 K，而且是以字符串的方式不回答**：`StuckRef`（无阀）、`StuckK500`、`StuckK1e8` 两两在 **2001/2001** 个状态打印同一个 token（全是单个 `0`），闩锁 **2000/2000** 步；`HangRef` 与 `HangK500` 的整条轨迹是同一个字符串。稳态侧更直接：8 根最终静止的管**都没有** `dpValve` 字段（`steadyNetworkApp.cpp:65` 只在非零时打），其中 `StopK1e8` 的阀若按它自己的初始速率 2.0e-4 m³/s 计价是 **1.7225e+10** Pa——也就是说这笔"如果被流动就要付"的账，在静止点上连一个数字都不留下。
* **阈值也不认识 K**：14 行 `static-hold threshold` 全部逐位是 `4 tau_0 L/d`（最差相对 **0.000e+00**），三根 Stuck 打出的都是同一个 **10400** Pa。
* **一动起来，阀就被收进同一行**：同一阈值（`tau_0 = 1.0`，4000 Pa 墙）配四个 K 给出 **9.330530e-06 > 9.259901e-06 > 8.997374e-06 > 7.206557e-06**，阀在终点占 1.0e4 Pa drive 的 **0 / 0.37 / 1.74 / 11.18%**；最大那一格把同一根管从 5.46 的 Buckingham 速率砍掉 **22.8%**。
* **每一格都落在 checker 自己的两律根上**：6 根动的管在瞬态侧最差 **5.419e-10**（`RollRef`），稳态侧最差 **9.091e-13**；`ReleaseK5000` 与 `RollK500` 的 12 位打印就是根本身（**7.2065569386e-06**、**2.7557710254e-05**）。5.46 那条"稳态括号"与 5.49 那条"瞬态静止支"在这里是同一条方程的两个解法。
* **稳态打印的拆分加起来是对的，而且阀那一项是它自己的公式**：`dp - rho g dz - dpValve` 折成壁应力是 **2.4908 / 2.4564 / 2.2205 / 4.1751** Pa，律在该应力下的速率与打印的 `Q` 相对差 ≤ **8.039e-13**；`dpValve` 对 `K rho v^2 / 2`（用打印的 `v`）相对差 ≤ **4.982e-12**——所以它不是外层拟合出来的份额。
* **停住的步号是惯性对墙的余量，不是阀**：`StopRef` / `StopK5000` / `StopK1e8` 在第 **17 / 17 / 5** 步达到恰好 `Q = 0`，其起点速率分别是 **5.408929e-08 / 4.647801e-08 / 2.607645e-08** m³/s，而上一步分别是 **1.390704e-07 / 1.306228e-07 / 1.184044e-07**——`q* = (dp_static - drive) dt / (rho L/A) = 7.391983e-08` 恰好夹在中间，三根管三次都是。而 5.49 那个"没有静止支时每步走 `q*`"的振幅在这里第一次被反用成**判据**。停住那一刻阀最大的一根也只有 **292.8** Pa，是 **2000** Pa 余量的 **15%**。
* **整条 14 × 2001 是那一行（两笔账 + 支）的代数**：checker 自己解带阀的隐式欧拉递推逐步对照，最差 **6.890e-10**（`StopK5000` 第 16 步），5 根被 hold 的管逐位 **0.000e+00**；递推自己给出的停住步号与打印的**完全相同**（`[17, 17, 5]` 对 `[17, 17, 5]`），所以那个停住不是内核额外加的阻尼。
* **跨内核终点对同一拓扑**：t = 4 对稳态内核最差 **5.428e-10**（`RollRef`），8 根静止的管在两侧都是 token `0`。
* **两笔账的量纲差别就是它们能被区分的原因**：同一个 K = 500 挂在两种墙下，阀降之比 **1635.1043 / 174.2973 = 9.3811** 恰好是速率之比 **3.0629** 的平方（相对 < 1e-6），而壁降之比 **16700.5482 / 9825.7027 = 1.6997** 对同样的速率差只给出指数 **0.474**——屈服律在近零处是 `sqrt` 型的（`Q ∝ (tw - tau_0)^2`），所以"阀值多少"取决于它挂在哪面墙后面，而不是一个可以合并进 `K` 的常数。

**成本**：2 个 case 文件（14 根管的定压瞬态案 + 同拓扑稳态对照案）、**81** 条判据、墙钟 **8.51 s**（ctest -j8 记，Test #79；单独跑 checker **5.0 s**，其中两次 app 调用只占 **0.166 s + 0.043 s**）。这一格便宜在 checker 把每步的隐式欧拉解改写在壁应力空间（单调，一次二分），而不是像 5.49 那样在速率空间重解带支递推——同族的 5.49 那一条在本轮**同一负载**下记 **338.38 s**（上一轮轻负载时记 140.14 s），两者差的不是求解器而是 checker 自己那 10 × 2001 步。既有测试一字未改：两笔账的加法与 5.49 的静止支都在同一条行里，本节只是第一次让它们同台；本轮全量回归 **105/105** 绿，实测 **736.20 s**（各用例时间之和 5630.04 s；本机当时另有前台应用与 `-j8` 争用，故绝对值高于上一轮的 553.98 s，相对次序不变）。

**本增量的边界**：① 只有 **Bingham**：14 根管全都不声明 `flowIndex`，checker 的剖面积分是按 `n = 1, k = mu` 写的，所以"HB × 阀"这一格没被证——5.47 那条"一个闭式包含 5.45 与 5.46"在带阀的方程上仍然只是逐条律成立。② **层流 only**：稳态打印的最大 `Re` 是 **16.65**（`RollRef`），屈服 + 湍流（Darby-Metzner 一类）本来就还没落地，阀也不可能替它作证。**（此处已在 5.51 部分更正：这一族的"湍流屈服"被拆成两半落地——转捩 locus 那一半（`He` 三段拟合 + locus 上的 `dp_t`）已在稳态那一行落成算例，摩擦闭包那一半仍未落地，所以"超过 locus 付多少"依然没有答案。）**③ **阀的时间表没有与静止支交叉**：`valveClose` / `valveOpen`（`:125`，装配 `:241`–`:245`，按 `t + dt` 取值 `:293`–`:300`）在这里全是常数 K；"阀正在关的那一步，被 hold 的管按哪个 K 决定它是否放开"这一格未验证。**（此处已在 5.52 落成算例，而答案是"哪一个都不按"：闩锁比的是 `drive` 与 `4 τ_0 L/d`，`drivePa` 里没有 K，所以斜坡只能移动停住的**步号**（第 54 / 72 / 184 步，单调于斜坡）而不能造成停住；晚于停住才开始的斜坡、以及已经 hold 在静止的管上的完整循环，都在 2001 个打印 token 与 2001 面旗标里逐字不可见；`t + dt` 与 `t` 的采样之差恰好读成一个 `dt`（54/72 对 55/73）。）** ④ **两端定压**：本案没有自由节点，所以 5.49 的节点封存（第 7 族）与网络自己放开（第 8 族）都没有带阀重做——带阀的 hold + 浮动节点三样同台是空白。**（5.52 沿用同一族判据，但桌面仍是全定压，这一格未动。）**（**此处已在 5.53 落成算例**：九对管，每对的自由节点夹在被 hold 的屈服管与牛顿泄放管之间，日程放在**邻居**那一行时它经由节点压力碰到静止支——放开的那一步走 **253 / 299 / 1830 / 从不**（`K = 0 / 5e3 / 1e6 / 1e8`），而墙回答的放开压力只有一个 `pA - 4 τ_0 L/d = 22000 Pa`（七次放开落在 **21990.9 – 21999.2 Pa**）。5.49 的节点封存那一族在同一张网上重做：被 hold 的管以它打印的那个**强制零**进入连续方程，九个节点最差闭合到 **4.966e-14** 对求解器自己的 **1e-11** 退出判据。docs/architecture.md 5.53。）⑤ 瞬态 app 不打印逐管压降拆分，瞬态侧的阀证据全部经由速率反证；`dpWall + dpValve == dp` 只在稳态侧被断言。⑥ 静止支依旧**刚**：阈值以下没有蠕变，阀也没有给静止点添加任何东西（这正是本案要证的），所以"静置 N 小时后重新启动"那一类还是空白。⑦ 耦合边两笔账都拒收（`:143` / `:234`–`:235` 不收阀，`:226`–`:229` 不收阈值），所以 1D↔3D 详细区里的"壁 + 阀"仍是空白。⑧ 每步 Newton 的退出门与节点定标里的 `+1.0`（`:414`）没有改：跨内核与对根的 **5e-10** 量级余量仍归求解器所有，与 5.48 边界③ / 5.49 边界⑦ 同一格。⑨ 对流动量账 `beta` 与静止支同台仍是空白（5.49 边界③ 的另一半，本节只关掉了阀那一半）。

### 5.51 屈服管的**层流-湍流转捩 locus** 第一次落在稳态那一行上：locus 属于流体和几何而不属于状态（六对孪生管在压差相差 750 倍时 `He / ReBt / Qt / dpQt / Vrelt / band` 六个量逐位相同，`0.000e+00`），而 locus 上那一点**这个管自己付多少**由 5.46/5.47 的壁律给出（`dp_t` 对 checker 重解的 Buckingham 括号最差 **2.12e-12**），同时**摩阻一个钱都没被重定价**（每个速率仍是所问压差在层流律下的根，最差 **3.2e-12**；`f Re_B g(ξ) = 64` 在超出 locus 四倍与 151 倍的管上仍然成立，最差 **3.4e-12**）

5.45 到 5.50 六条增量各自把"湍流屈服"列成自己的边界，理由是同一句：`pipeFrictionPa`（`src/1d/src/network.cpp:393`–`:395`）一旦看见声明了阈值就**不看 `Re`**，直接回 `4 L τ_w / d`。这件事在这一节被拆成两半来落地，因为文献只让其中一半落得下来。Wilson & Thomas（*Can. J. Chem. Eng.* 84:520, 2006）证明了 Bingham 塑性的**转捩条件只依赖 Hedström 数** `He = τ_0 d² ρ / μ²`，并给出三段拟合 `Re_B,t = 2100`（`He ≤ 1700`）、`80 He^0.40`（`1700 < He ≤ 1e5`）、`25 He^0.5`（`He > 1e5`），其中 `Re_B = ρ v d / μ` 是塑性 Reynolds 数、`Trel = V/(τ_0/ρ)^0.5 = Re_B,t/He^0.5` 是恒等式而不是相关式。**这一整套只需要分组和常数，不需要湍流摩擦律**；而摩擦律那一条（他们的 Eq 15/16，含 `β = τ_w/τ_0` 的括号写法与 `11.6/χ` 钝化项）在本机没有 poppler 的条件下用 stdlib zlib 抽取时上下标位置丢失，两篇原文抽出来的形式互相差一个 `β^0.5`，且按抽取的写法 `β→∞` 时 `11.6/χ` 发散、退不回牛顿极限——**未确认的转写不进内核**，所以这一节报 locus 而不改摩阻。落地面：内核新增一个导出函数 `yieldedTransitionOf`（`network.cpp:413`）与 `PipeFlowResult::transition`，`solvePipeFlow` 在每个管上填一次（`:501`），稳态 app 打七个新 token（`He / ReB / ReBt / Qt / dpQt / Vrelt / band` 加 `pastTransition`/`belowTransition`）；**case.dict 没有任何新键**，瞬态那一行与 `pipeFrictionPa` 一字未改。

关键的一点是 `dp_t`：把 locus 的速率 `Q_t = Re_B,t μ A/(ρ d)` 送回**已经装好的**壁律，就得到"这根管离开层流陈述要付多少压差"。这是 Pipesim/OLGA 的转捩/层流包络线那一类输出，而它不需要任何新物理，只需要有人把两笔账放在一起。

* **locus 不是状态的函数**：`Small01/Small05`（200 / 500 Pa）、`TableHe04/TableHe04Run`（静止 / 流动）、`BigHe06/BigHe06Run`、`Bore02/Bore02Run`、`Edge1700Over/Edge1700Flat` 五对同几何管在压差相差最多 **750 倍**时六个 locus 量**逐位相同（0.000e+00）**；被 hold 住的 `TableHe04`（`Q 0`，静态阈值 800 Pa 顶住 500 Pa 的桌面）照样报出 `Qt 0.00597521486964`。这条判据之所以能被失败，是因为如果 locus 是从速率算的，孪生管就会分开。
* **三段拟合就是打印出来的那三段**：`He` 命中已发表参数集（Csizmadia & Tóth 2013 表 1：`D = 0.1 m`、`ρ = 1370`、`τ_y = 1 Pa`、`μ_B = 2e-2 Pa·s` → `He = 3.425e4`）到 12 位里的 4 位（`342.5`、`34250`、`3425000`、`3.425e8`、`3.425e10`），`Re_B,t` 对 Python 重算的带公式最差相对 **3.04e-12**（这个余量就是 app 的 12 位打印本身），带 0 逐字是 `2100`。
* **恒等式被读回成第二个已发表形式**：`Trel = Re_B,t/He^0.5` 对 15 根管最差 **2.14e-12**；于是在 `He^0.5` 带里 `Trel` 必须**恰好是 25**（打印就是 `25`，5 根管），在中带里退成 `80 He^-0.10`（`TableHe04` 打 **28.1594345251**、`Edge1e5Under` 打 **25.2982282384**）。也就是说论文自己的 Eq (19)/(20) 两种写法在这里是**同一条恒等式的两个投影**，而不是两条独立拟合。
* **`dp_t` 是这根管自己的壁律，不是相关式**：checker 用 5.46 的括号 `F(tw) = tw − (4/3)τ_0 + τ_0⁴/(3 tw³) = b` 重解 `tw(Q_t)` 再乘 `4L/d`，15 根管最差相对 **2.12e-12**（`BigHe06`）；并且 15 根管全部满足 `dp_t > 4 τ_0 L/d`——**离开庭阈值不够，还要再付到转捩才算离开层流**。数值：`TableHe04` **2024.01 Pa**（阈值 800 Pa）、`Small01` **403.075 Pa**（阈值 8 Pa）、`BigHe06` **1.02712e5 Pa**、`BigHe08` **8.62566e6 Pa**。
* **判决有三个独立的入口，而且它们同意**：`pastLaminarLimit` 逐管等于 `Q > Q_t`、等于 `ReB > ReBt`、等于 `dp > dp_t`，15 管 15 条；本案真的跨在两侧（**7 根 past、8 根 below**），`Small01/Small05` 是同一根管在 `dp_t = 403.075 Pa` 两边各站一步（200 / 500 Pa）。
* **读取 locus 没有碰任何一笔账**：每个动的管的速率仍是所问压差在层流屈服律下的根（最差 **3.2e-12**，`Small01`）；`f Re_B g(ξ) = 64` 在 `Edge1e5Under`（`Re_B = 7.86068e5`，**超出自己的 locus 151 倍**）、`Edge1e5Over` 与 `BigHe06Run`（`Re_B = 2.53550e5`）上仍成立，最差 **3.4e-12**。这正是下一格必须停止复现的那个数字。
* **`Re_B` 就是打印的 `Re`**：15 根 Bingham 管的 `ReB` 与 `Re` 相对差 **0.000e+00**（同一句 `ρ v d / μ`，两条打印路径）；`τ_0 = 2e-4` 的 `Bore02` 与它的牛顿孪生 `Bore02Newt` 在 1e-3 内符合 `Q = Q_Newt · g`，`g = 0.994667`（`ξ = 1/250`）——带 0 与那个 2100 是**牛顿极限的另一种读法**，不是另一条陈述。
* **两个带的接缝是量出来的台阶，不做平滑**：`He = 1e5` 处打印 **7999.99119999 → 7906.49554955**，高带比中带**低 1.17%**；`He = 1700` 处 **2100 → 1568.29172954**，**低 25.3%**，因为论文用来跨过去的 Eq (17) 的常数没能转写。两个方向都是"向下"，所以在两个边界上取哪一段是**未定的**，本节只把它记成观测而不是判据。带内以及带 1→带 2 之间 locus 对 `He` 单调不减。
* **locus 不认识不是 Bingham 的律**：`NewtRef`（无阈值）与 `HbRef`（阈值 + `flowIndex 0.5`）都**不打** `He`；后者打印的 `Re` 是表观数 `8 ρ v²/τ_w`（与 checker 重算差 **< 1e-9**），也就是说 HB 的律里根本没有塑性黏度，`He` 与 `Re_B` 都不是关于它的陈述——这是 5.47"一条闭式包含两条律"在**分组层面**的反面。

**成本**：1 个 case 文件（17 根管、6 个定压节点、五个桌面 200/500/5000/3e4/1.5e5 Pa）、**60** 条判据、ctest **Test #80 记 0.29 s**（单独跑 checker 同量级：一次 app 调用；这是全表最便宜的一条 verify 用例，因为 locus 是纯读出，零步进、零重解网格）。既有测试一字未改：locus 是新增的读出，不在任何一条摩阻路径上。本轮全量回归 **106/106** 绿，实测 **939.07 s**（各用例时间之和 7071.13 s，最重五条耦合用例合计 3563.00 s = 50.4%）。

**本增量的边界**：① 只有**稳态**那一侧打印；瞬态行调用同一个 `pipeFrictionPa`，locus 函数就在同一层，但目前不逐层打（`transientApp.cpp:84`–`:92` 仍是 `Q` 与 `static hold`）。② **超过 locus 之后什么都没说**：湍流屈服**摩擦**律仍未落地，`He`、`Re_B,t`、`Q_t`、`dp_t` 是"到哪里为止这话算数"，不是"过了这一点付多少"。③ 带 0 在 `He < 1700` 内的**形状**未落地（论文用 Eq (17) 桥接，其常数转写失败），所以带 0 是平的 2100。④ 两个接缝各自向下台阶（1.17%、25.3%），本案**记录而未给准则**：设计者在边界上该取哪一段没有断言。⑤ 只有 Bingham；HB 与幂律拿不到 locus，理由与 5.47 的表观 `Re` 同源。⑥ 等温、常物性：`μ(p)` 表驱动管（3.35 那一条）的 locus 会随填充变化，本案未涉及。⑦ **无粗糙度**：Thomas & Wilson（Hydrotransport 17, 2007）把同一分析推到粗糙壁并用 `Re* = k U*/(ν)` 与 `k/δ_k` 分类，其转捩式同样没转写；本案全部管 `roughness 0`。⑧ `BigHe08`/`BigHe10` 的 `τ_0 = 1e4 / 1e6 Pa` 是 `He^0.5` 带的**数值探针**，不是泥浆；它们的速率全被静态阈值顶住（`Q 0`），这恰好也证明 locus 在静止管上照样可读。⑨ 本案**不做**"past 就拒收"：在损失函数里抛错会在 Illinois 括号展开的 60 步里必然触发（哪怕答案是层流的），所以这一格是**报告**；真正的湍流闭包到来时应当**替换**答案而不是标记它，届时第 6 条那批 64 的复现断言要一起改掉。⑩ 没有对**实测**转捩数据做验证：论文的图 5–8 是图像，本案只对齐了它正文写出的分组与常数。

### 5.52 阀的**时间表**第一次与屈服壁的**静止支**同处一条行：日程移动的是**停住的步号**而不是**停不停**（第 **54 / 72 / 184** 步，斜坡收尾于 `0.15 / 0.6 / 从不`，单调），**晚于停住才开始的斜坡在整个 trace 里逐字不存在**（`StopLate == StopPlain`，2001 个 token 加 2001 面闩锁旗标全同），而阈值只看 `drive` 不看 `K`，所以阈值以上的管在一个价值 **9.8 倍整条 drive** 的关闭阀下仍然 **0/2001** held——**内核一字未改**

5.50 把自己的第三条边界写成"阀的时间表没有与静止支交叉"：那一格 14 根管的 `K` 全是常数，所以"阀正在关的那一步，被 hold 的管按哪个 `K` 决定它是否放开"没有答案。这一格把交叉做出来，而答案是**这两个开关读的不是同一个量**：闩锁 `holdsStatic`（`src/1d/src/transient.cpp:40`–`:43`）比的是 `drive` 与 `4 τ_0 L/d`，而 `drivePa`（`:52`–`:56`）里根本没有 `K`；斜坡 `valveKAt`（`:125`，采样 `:293`–`:300`：取 `t + dt` 的值并在 Newton 期间保持）只给**在动的那一行**定价。所以日程能碰到的只有后者，静止支既不吃 `K` 也不看日程。落地面只有测试：`tests/cases/transientYieldValveSchedule/case.dict`（12 根管，`ρ 850 / μ 0.1`，`L 20 / d 0.02`，除两根牛顿孪生外 `τ_0 1.5` → 阈值 **6000 Pa**；四个桌面：阈值以上 `1e4 Pa` 的关闭/循环/关到底/写死关闭，阈值以下 `4000 Pa` 的四种斜坡 × 无阀对照 + 两根静止在阈值下 + 两根同日程的牛顿孪生）与 `case_steady.dict`（同拓扑、每个阀冻结在 `t = 4` 的 `K`）、`tests/check_transient_yield_valve_schedule.py`。

* **阈值以上，阀是折扣而不是拒绝**：`Creep*` 四根在 `drive 1.0e4 Pa` 对 **6000 Pa** 墙下 **held 0/2001**，一根都不被 hold，尽管 `K = 1e6` 在开阀速率上计价 **98192 Pa = 9.8 倍整条 drive**；最小的动态速率是 **1.302827e-07** m³/s（关闭期间）而不是 0，全轨迹只有初值那一个 0。`CreepShut`（关到底就不动）与 `CreepK1e6`（一开始就写死关闭）**轨迹不同**——斜坡确实在行里；但两者的尾巴**逐位**汇到同一个根（开/关两桌各 **4.77522083e-06**、**7.93552215e-07** m³/s，rel `0.000e+00`）：日程在斜坡结束之后就不留痕迹。
* **阈值以下，是墙在停而日程只说什么时候停**：`StopPlain`（无阀）第 **184** 步、`StopMid`（斜坡 0.05→0.6 s）第 **72** 步、`StopEarly`（0.05→0.15 s）第 **54** 步，单调于斜坡；之后的闩锁旗标分别 **1816 / 1928 / 1946** 面（对 2001）。停住前最后一个动态速率分别 **3.276219e-10 / 4.251551e-08 / 7.076481e-08**——**被阀提前判停的管在更大的速率上穿过零**，因为那次穿过是钳制给的（`:353`–`:356`）而不是墙的近因渐近。
* **晚到的斜坡在整个 trace 里不存在**：`StopLate` 的阀 0.5 s 才开始关，而墙在 0.368 s 已经停了它，于是它的 **2001 个打印 token 逐字**等于 `StopPlain`，**2001 面旗标**也 state-for-state 相同。这条判据之所以能被失败，是因为同一张网里 `StopEarly/StopMid` 用同一形状的日程却在第 54/72 步分开——如果静止支读 `K`，`StopLate` 就会分开而不是相同。
* **没有日程能叫醒墙 hold 住的管**：`StopEarly/StopMid/StopLate` 在各自阀门开完（`t = 1.65 s`）之后的 **1176** 个状态仍然恰好 0，且这 1176 面旗标都说 hold；同一日程的牛顿孪生 `NewtEarly`（无 `yieldStress`）从未静止（最低 **9.064733e-07** 对起步的 2.0e-4），并在开阀后回到自己无阀的速率（尾巴对 `NewtPlain` rel **5.462e-10**）——不对称来自阈值而不是阀。而 **5.462e-10 而不是逐字相同**：牛顿那一行到 `t = 4 s` 仍在弛豫。
* **已被 hold 在静止的管对日程完全不作答**：`HeldCyc` 跑完整 close/hold/open 循环，**2001 个 token 逐字**就是 `HeldPlain` 的（每一个都是 0），两条轨迹对递推在**绝对与相对两个单位下都是 0.0**。
* **12 × 2001 = 24012 个状态由本文件自己的两律递推复现，两个门一起立**：动态状态（`|Q| ≥ 1e-8`）相对最差 **2.140e-09**（`StopPlain` 第 182 步），而全部 24012 个状态的绝对差不超过 **2.666e-15** m³/s（最差是 `NewtEarly` 末态）。两个门都要，因为内核 Newton 的退出（`:418`，定标 `:389`–`:391` 用 `drive` + 壁降 ≈1e4 Pa 而不是速率）在停住前那几个 `Q ~ 1e-10` 的状态上给的是**绝对**精度：那 2.14e-09 属于退出判据而不属于律——与 5.49 边界⑦ / 5.50 边界⑧ 同一格，本案第一次把它写成"相对门只在有速率的地方成立"。递推算出的停住步号与打印的 **4/4** 相同。
* **`t + dt` 这个采样点是被两个方向证的**：斜坡与动态重叠的五根管，按 `t` 采样的递推对 app 分别差 **5.526e-02**（`CreepCyc`/`CreepShut`）、**3.055e-01**（`NewtEarly`）、**1.000e+00**（`StopEarly`/`StopMid`，即"该停的那一步还在动"），而按 `t + dt` 的递推差 ≤ **5.187e-10**；最干脆的是步号读数：app 停在 **54 / 72**，按 `t` 采样的递推停在 **55 / 73**，差的正是**一个 `dt`**。反过来，斜坡从不与动态重叠的两根（`StopLate`、`HeldCyc`）两种采样给的是**同一条轨迹**（`0.000e+00` 之差）——采样约定在那里不可观测，理由与上一条相同：静止的行里没有 `K`。
* **稳态那一侧的对应物**：每个日程在 `t = 3.0 s` 已经结束并停在 family 7 冻结的那个 `K` 上（逐管断言）；march 的尾巴对上稳态内核 **< 2e-9**（`CreepOpen` 逐字 `0.000e+00`、`CreepCyc` 2.304e-11），而那个根就是本文件对 `drive = dp_wall + K ρ Q²/2A²` 的二分；阈值以下的六根管稳态内核打 `Q 0`、`dp 4000` 并且**连 `dpValve` 字段都不出现**——5.50 稳态侧那条在带日程的网络上重做一遍；唯一一根"阀关着且在动"的 `CreepShut` 打 `dp 10000 dpValve 2711.69100797`，对 `K ρ v²/2` 差 **1.270e-12**，而留给壁的 **7288.309 Pa** 高于它自己的 6000 Pa 阈值——这既是它没被 hold 的理由，也被写成判据。
* **阈值对每根屈服管是同一个数，`K` 不在里面**：10 行 `static-hold threshold 6000`（`4 τ_0 L/d`），带日程与不带日程的完全一致，两根牛顿管一行都不打。

**成本**：2 个 case 文件（12 根管、3 个定压节点、`dt 0.002 / endTime 4.0`）、**121** 条判据、ctest **Test #81 记 11.14 s**（单独跑 **8.45 s**：2 次 app 调用可忽略，成本在 12 × 2001 次壁应力空间的二分）。既有测试一字未改，内核也没动：阀项（`:370`–`:371`/`:381`/`:383`）、穿过钳制（`:353`–`:356`）与 Newton（`:318`）全是 5.38 / 5.49 / 5.50 已经装好的，这一格需要的是**被量而不是被建**。本轮全量回归 **107/107** 绿，实测 **775.07 s**（各用例时间之和 5893.78 s，最重五条耦合用例合计 3038.87 s = 51.6%）。

**本增量的边界**：① 只到**常物性单相**瞬态行：对流动量账 `beta` 与静止支同台仍是空白（5.49 边界③ 的另一半，本节只关掉阀那一半）。② 两端全定压，所以 5.49 的节点封存与"网络自己放开"没有带日程重做——`hold` × 浮动节点 × 时间表三样同台仍空（5.50 边界④）。**（此处已在 5.53 落成算例，而答案是"日程拥有放开，墙拥有压力"：邻居那一行的日程把放开那一步推到 253 / 299 / 1830 / 从不，同一处方放在被 hold 的那一行则连放开的那一步都不移动（还是 253，节点历史到第 252 个状态逐位相同）；七次放开全在 22000 Pa 的 ±12 Pa 里；而"先封后开"的邻居在 1290 步叫醒墙 hold 住的管，比日程第一次动东西晚 0.578 s——叫醒它的是网络，不是阀。）**③ 静止支依旧**刚**：阈值以下没有蠕变、没有触变、没有"静置 N 小时后重新启动要多少压差"，日程也没给静止点添加任何东西（这正是本案要证的）。④ 只有 Bingham（`n = 1`）：HB × 阀 × 日程未证，而 5.47 那句"一个闭式包含 5.45 与 5.46"在带日程的行上仍只是逐条律成立。⑤ 日程只用到线性斜坡与 `closeStart < 0` / `openStart < 0` 的组合；斜坡的**形状**由 checker 逐点核对（`2.5e5 / 5e5 / 7.5e5`，端点恰好 `0` 与 `1e6`），但非线性日程与"同一步内跨过 `closeEnd`/`openStart`"没有断言。⑥ 停住的**步号**分辨率是 `dt`（0.002 s），细于它的开关看不见（5.49 边界② 同一格），而上面那条 `55 / 73` 对 `54 / 72` 正是这条边界的直接读数。⑦ 每步 Newton 的退出门与定标里的 `+1.0` 未改：本增量第一次把**动/静两个单位下的门**同时立起来（相对 2.14e-09 只在有速率处成立、绝对 2.67e-15 全域成立），但"多少位以内算合格"仍未立为标准。⑧ 层流 only：本案阈值以上那组的速率仍在 5.51 的 locus 以下，past-locus 付多少归 5.51 边界②。⑨ 耦合边两笔账都拒收（`:143`/`:234`–`:235` 不收阀、`:226`–`:229` 不收阈值），所以 1D↔3D 详细区里的"壁 + 阀 + 日程"仍是空白。⑩ 跨内核终点对的是**同一拓扑的稳态内核**，不是实测或场数据；`NewtEarly` 的回归还带着 5.462e-10 的弛豫余量。

### 5.53 **浮动节点**、屈服壁的**静止支**与阀的**时间表**第一次同处一题：日程放在**邻居**那一行时它拥有放开（第 **253 / 299 / 1830 / 从不** 步，单调于 `K = 0 / 5e3 / 1e6 / 1e8`），放在**被 hold 的那一行**时它连放开的那一步都不移动；而墙在任何 `K` 下只回答**一个压力**（七次放开全落在 `pA - 4 τ_0 L/d = 22000 Pa` 的 **21990.9 – 21999.2 Pa** 里）——**内核一字未改**

5.50 把自己的第四条边界写成"两端定压，所以没有自由节点"，5.52 沿用同一个桌面。这一格把桌面改掉而**不动内核**：每对管子是"一根 `τ_0 = 2.0` 的屈服管 `A → F`（阈值 `4 τ_0 L/d =` **8000 Pa**，入口 drive `3.0e4 - 2.5e4 =` **5000 Pa**，余量 **3000 Pa**）串一根牛顿泄放管 `F → B`"，而 `F` 是自由节点（`volume 0.1 / bulkModulus 1e8`，即 `Cc =` **1e-9 m³/Pa**）。于是闩锁比的那个 drive 不再是算例文件里的数，而是**解出来的**——每对的行系统是

```
ρL/A (Qh - Qh')/dt = (pA - pF) - 4L tw_h/d - Kh ρ|Qh|Qh/(2A²)
ρL/A (Qb - Qb')/dt = (pF - pB) - 128 μL Qb/(π d⁴) - Kb(t) ρ|Qb|Qb/(2A²)
V/BC (pF - pF')      = Qh - Qb,        Kb(t) = valveKAt(t + dt)
```

`Qh` 在墙 hold 住时以**强制零**进最后一行（`transient.cpp:356`–`:362` 去掉那一行的压力列，而节点行仍留着它在 `Qh` 列上的 ±1），这就是 5.49 的节点封存在带日程的网络上的对应物。落地面：`tests/cases/transientYieldValveNodeRelease/case.dict`（11 个节点、18 根管，九个自由节点，`dt 0.002 / endTime 4.0`）、`case_long.dict`（只有 `K8`/`R1` 两对，`endTime 30.0 / reportEvery 5`）、`case_steady.dict`（同拓扑，`solver { type steadyNetwork; }`、自由节点写成 `F {};`、每个阀冻结在 `t = 4` 的 `K`，并带一个只设定残差速率标度的 `sources` 块——没有它，稳态残差冲到 **3.8116e+279** 并报 `NOT converged`）与 `tests/check_transient_yield_valve_node_release.py`。

* **日程放在邻居那一行，它拥有放开这件事**：泄放管上写死的 `K = 0 / 5e3 / 1e6` 把屈服管的放开推到第 **253 / 299 / 1830** 步（`t = 0.506 / 0.598 / 3.660` s，单调于 `K`），`K = 1e8` 把它推出整个窗口。0.2→0.3 s 关到 `1e8` 的斜坡（`R1`）**取消**放开，1.0→1.1 s 的同形状斜坡（`R2`）在 `K0` 那一步之后才可见：`R2` 的管与节点在**前 501 个状态**上与 `K0` 逐位相同（节点 **20656.051341** Pa 两边一样），第一个不同 token 在 **501** = `int(1.0/dt) + 1`，也就是采样点第一次落进斜坡的那一步——**同一个 `closedK`、同一段斜坡长度，两个结局，差别只在到达时间**。
* **日程放在被 hold 的那一行，它什么都拥有不了**：`O` 的 `1e8` 直接装在屈服管上、`OW` 的整段 close/hold/open 循环装在屈服管上，两根都在 `K0` 的那一步 **253** 放开，节点历史到 **252** 逐位相同；区别只在第一个动态状态上少掉的那点速率（`9.949789e-11` 对 `9.954104e-11` m³/s，即阀从第一个动的状态起就开始收钱）与尾巴上 **2.60e-05** 的弛豫余量。5.52 那句"静止的行里没有 `K`"在这里仍然是对的，而**多出来的一半**是：静止的行一旦放开，`K` 立刻收钱——它只是收不到"要不要放开"这一票。
* **叫醒墙 hold 住的管的是网络，不是阀**：`W` 的泄放管 `0 → 0.002 s` 关到 `1e8`、`2.0 → 2.1 s` 开回来。泄放管与节点从第 **1001** 个状态（`t + dt` 第一次落进开阀斜坡）开始分岔，而**被 hold 的屈服管仍打印 `0` 到第 1289 个状态**，放开在第 **1290** 步 = `t = 2.580` s——比日程第一次动东西晚 **0.578 s**。这条 289 步的延迟是本案测量出来的而不是假设的：5.52 的"晚到的斜坡逐字不存在"在带节点的网络上变成了"晚到的斜坡先动节点，节点再动墙"。
* **墙只回答一个压力**：`main["thr"]` 里 9 行 `static-hold threshold 8000`（`4 τ_0 L/d`），九根牛顿泄放管一行都不打——`K` 不在阈值里，跨三个数量级与两种日程都是。七次放开发生的节点压力 **21990.85 – 21999.24 Pa**（跨度 **8.4 Pa**）对 `pA - 8000 = 22000 Pa`，且每次都满足"上一个采样 `≤ 8000` < 本步入口 drive"。
* **放开是单向的，而且量得出为什么**：没有一个节点在放开之后重新被 hold（放开后峰值 **21980.9 – 21983.0 Pa** < 22000），而**唯一**能把节点压回放开压力之上的流就是那根屈服管自己的速率——它恰在那个压力处为零。网络自己的落点也在下面：那三对在放开之后又被封上的桌面，稳态内核给的节点是 **21634.2**（`R2`）与 **20450.6**（`W`、`OW`）Pa，全部低于 22000。
* **`held for 4 s` 是日历而不是墙**：`K8`/`R1` 在 4 s 窗口里 **2001** 个状态全打印 `0`（节点最深 drive 分别 **5419.36** 与 **6395.88** Pa 对 8000，而它们的泄放管全程在动：**8.040e-08 – 1.071e-07** 与 **9.090e-08 – 7.353e-06** m³/s），所以不是"闩锁睡过了"），而**同一拓扑**的稳态内核没有闩锁可等，报 **6.101454e-08** m³/s。30 s 那一跑把 `R1` 在第 **2425** 个采样 = **24.25 s** 交回来（最后一个被 hold 的采样 drive **7999.4773** Pa，下一步越过 **8000**），`K8` 到 30 s 仍差 **304.35 Pa**（drive **7695.6525**），而它的 drive 从 5000 Pa 起**单调上升**——"从不"是一个读数，不是判决。日历本身是算术：hold 期间节点每步掉 **0.2051 Pa**，即 **102.55 Pa/s**，正是泄放管自己那 **1.02549e-07** m³/s 除以节点的 **1e-9** m³/Pa；剩下的 2581 Pa 按同一速度要 **25.2 s**。
* **`9 × 2001 = 18009` 个状态加两对长跑（各 15000 步进）由本文件自己的两条行、节点行与采样 `K` 复现**：动态状态（`|Q| ≥ 1e-8`）相对最差 **1.514e-09**，节点相对最差 **4.610e-12**，全域绝对差不超过 **1.104e-16** m³/s；递推算出的放开步号与打印的 **9/9** 相同。相对门只在有速率处成立、绝对门全域成立，与 5.52 的动/静两个单位同一格。
* **被 hold 的管以它打印的那个强制零进连续方程**：九个节点对 `V/BC·dp/dt = Qh - Qb` 的最差定标残差 **4.966e-14** 对内核 Newton 自己的退出判据 **1e-11**（节点行定标里那个 `+1.0`，`:414`），即 `6.705e-09` 倍的最大移动速率（**7.406e-06** m³/s）——把它当成漏，需要 **1.0e+08** s 才能移动一个节点，而它自己的解在 4 s 里移动了 5000 Pa。
* **与稳态内核的对应**：九个冻结桌面的速率全是本文件对 `1e4 Pa = 4L tw/d + 128μLQ/(πd⁴) + (Kh+Kb)ρQ²/(2A²)` 的二分，最差 **1.96e-08**（在 `1e8` 那一组），`K0/W/OW` 那一组（阀门全开）**6.14e-13**；`dp_h + dp_b` 对 **1e4 Pa** 差 **4.0e-09** Pa、同一对里两根的速率差 **7.4e-10**；六把冻结阀的 `dpValve` 对 `K ρ Q²/(2A²)` 差 **2.44e-12**。**阀坐在节点的哪一侧，落定的速率不知道，节点知道**：三个"阀门最后全开"的桌面（`K0`、`W` 的循环在邻居上、`OW` 的循环在屈服管自己上）汇到同一个 **8.84757194e-07** m³/s，四个"一把 `1e8`"的桌面（三次在泄放管、一次在屈服管）汇到 **6.10145e-08**（跨度 **1.14e-08**）——但 `Fk8` 停在 **21634.1567** Pa 而 `Fo` 停在 **20031.0744** Pa，两者之差 **1603.08229** Pa 就是那把阀自己打印的降 **1603.08222** Pa，差 **6.5e-05** Pa。
* **尾巴按各自剩下的时钟定标**：`K0` **1.84e-08** / `K3` **8.81e-07** / `OW` **2.60e-05** / `O` **3.28e-05** / `W` **1.24e-03**（节点各为 **3.4e-10 / 2.8e-08 / 1.1e-06 / 5.5e-06 / 2.9e-05**）。这五格里**速率都比驱动它的压力离落点更远**，倍数 **5.9 – 53.8**：屈服律在近壁处 `q ∝ (tw - τ_0)²`，所以 1e-5 的应力余量对应 1e-3 的速率余量——本案的门因此**逐对**写而不是一个数。两对明说是半路：`K6` 放开后只剩 **0.34 s**，而按最后五个采样的 **571.6 Pa/s** 还需要 **1.43 s**（尾巴 rel **0.955**，节点 **21793.69** 对 **20978.18**）；`R2` 从**另一侧**走向那个 `1e8` 涓流（rel **0.580**，节点 **21426.42** 上升向 **21634.16**，速率 **1.45527e-07 → 1.45185e-07** 下降向 **6.10145e-08**）。
* **放开之后，闩锁与没有闩锁落在同一个数上**：`K8` 的泄放速率从**上方**走向稳态涓流（**1.02549e-07** @4 s → **7.25646e-08** @30 s → **6.10145e-08** 稳态；7.5 倍时钟收掉 **61%** 的差），而 30 s 里放开的 `R1` 其屈服管速率 **3.68432e-08** 仍在上升，走向的正是稳态内核那**同一个** **6.10145e-08**（节点 drive **8282.29 Pa** 已在阈值之上）。

**成本**：3 个 case 文件（主桌面 11 个节点 / 18 根管 / `dt 0.002` × 2000 步，长桌面 2 对 × 15000 步，稳态桌面同拓扑 + `sources`）、**104** 条判据、ctest **Test #82 记 8.51 s**（单独跑 **4.61 s**：三次 app 调用合计只有 **0.25 s**——主桌面 0.071、长桌面 0.099、稳态 0.076 s——成本全在 checker 自己那 `9 × 2001` 加两对 `2 × 15000` 次的壁应力空间二分与节点定点）。既有测试一字未改，内核也没动：静止支（`:310`–`:316`/`:353`–`:362`）、阀项（`:370`–`:371`/`:381`/`:383`）、节点行（`:395`–`:416`）与采样约定（`:293`–`:300`）全是 5.38 / 5.49 / 5.50 / 5.52 已经装好的，这一格需要的是**被量而不是被建**。本轮全量回归 **108/108** 绿，实测 **838.99 s**（各用例时间之和 6384.61 s，最重五条耦合用例合计 3152.26 s = 49.4%）。

**本增量的边界**：① 只到**常物性单相**瞬态行的 **Bingham**（九根屈服管都不声明 `flowIndex`）：HB × 日程 × 浮动节点未证，对流动量账 `beta` 与静止支同台仍是空白（5.49 边界③ / 5.52 边界① 同一格）。② 每个自由节点只挂**一对**（一进一出）：分支节点（三根以上）与 5.20 那种"两节点跨管"里没有重做这一格，所以"日程经更宽的网络碰到静止支"仍是空白。③ 节点的压缩性只用液相那条 `V/BC` 行；两流体 / 气相 EOS 的可压节点（3.9、5.13 那一族）不带静止支。④ 静止支依旧**刚**：阈值以下没有蠕变、没有触变，也没有"静置 N 小时后重新启动要多少压差"——本案那句"仍差 304.35 Pa"是**日历**读出来的，不是蠕变模型给的。⑤ 日程形状仍是线性斜坡加 `closeStart` / `openStart` 的组合；新的是**日程在哪一行**，不是新的斜坡。⑥ 放开只发生**一次**且结构上单向：本案的拓扑做不到"网络把节点压回去再锁上"，因为唯一的回灌流就是被 hold 的那根自己的速率，它恰在放开压力处为 0——**双向闩锁仍空**。⑦ 长桌面只跑了 `never` 的两对中的两对（成本）：本案断言"30 s 放开 `R1`、30 s 仍 hold `K8`"，但没有断言"所有 `never` 在某个有限时钟内都会放开"。⑧ 每步 Newton 的退出门与定标里的 `+1.0` 未改：**4.966e-14** 属于求解器而不属于律（5.48③ / 5.49⑦ / 5.52⑦ 同一格）。⑨ 层流 only：最大速率 **7.406e-06** m³/s（起步那一瞬）给 **0.0236** m/s、`Re = ρ v d / μ ≈` **4.0**，稳态九个打印 `Re_B` 全在 **0.48** 以下，都离 5.51 的 locus（`Re_B,t` 带 **1700 / 2100**）差三个数量级；past-locus 付多少归 5.51 边界②。⑩ 耦合边的两笔账都拒收（`:143`/`:234`–`:235` 不收阀、`:226`–`:229` 不收阈值），1D↔3D 详细区里的"壁 + 阀 + 日程 + 节点"仍是空白；跨内核终点对的是**同一拓扑的稳态内核**，不是实测或场数据。

### 5.54 **分支节点**上的两根屈服壁第一次让闩锁**双向**：第二根墙给节点一个它原来没有的**天花板**，于是同一份日程在**放开之后还能把它重新锁上**（**11** 根壁放开、共 **14** 次放开事件，对 **3** 次重锁，全文件 **17** 次翻），而谁能回去睡觉由**阈值的次序**决定、不由日历决定——**内核一字未改**

5.53 把两条边界写成结构性的：② 它的分支节点只有一进一出，⑥ 放开是单向的、因为一根墙自己的流恰在它需要的那个压力上归零，节点只能从下方到达并停住。这一格两条同时让位，而动的只是桌面：每臂从"一节点一壁"改成"**两根屈服壁进、一根带日程的泄放管出**"，七臂 **13** 根屈服壁、**8** 根牛顿管、9 个节点。内核不需要改，因为要用的三件东西本来就在了——节点行对任意度数写（`transient.cpp:395`–`:416`），闩锁是一条**状态**测试而不是历史测试（`:40`–`:43`：`Q_prev == 0` 且起步 drive 在 `4 τ_0 L/d` 之下，所以速率真的回到 0 时重新锁上是合法的），穿过零的那条内步钳也在（`:353`–`:355`）。每臂的行是

```
ρL/A (Qk - Qk')/dt = (pA - pF) - 4L tw_k/d - Kk ρ|Qk|Qk/(2A²)    k = 该节点上每根屈服壁
ρL/A (Qd - Qd')/dt = (pF - pB) - R Qd - Kd(t) ρ|Qd|Qd/(2A²)
V/BC (pF - pF')    = Σk Qk - Qd,     Kd(t) = valveKAt(t + dt)
```

* **第二根墙给的是天花板，不是第二个开关**：只要**弱**墙还在流，网络的极限就是**弱墙自己的**放开压力，每根比它强的墙都在上行路上从它下面穿过——这正是 5.53 那根孤壁永远做不到的动作。`Sa` 因此有三步翻（**333 / 507 / 1530** = 醒·睡·再醒），`Wt`/`Ws` 同样，而 `Wa`（6000 在 8000 之下）只翻一次：它的放开压力就是本臂天花板，从下方逼近。
* **闩锁两个方向读同一个数，所以没有迟滞可拟合**：**14** 次放开事件的起步 drive 超出自家阈值 **1.00 – 16.54 Pa**，**3** 次重锁差 **3.10 – 3.55 Pa** 不到**同一个**数（`Sa` 在 22003.10、`Wt`/`Ws` 在 24003.55 Pa 对 22000 / 24000），两侧都恰好差一步的节点行程。
* **次序而非日历**：臂 `i` 两根同为 8000，谁也不严格更弱，于是**都在第 318 步放开、都不重锁**；臂 `u` 一根独壁同样不重锁。不等式还会点第三根与第五根，文件各给一个理由：`St` **从未被问**（节点最低 **22593.53 Pa**，比它需要的 22000 高 **593.53 Pa**，两根更弱的墙全程替它把着门），`Se` 放开后节点再没回到阈值下（放开之后 `Fe` 最大 **21979.46 Pa**，差 **20.54 Pa**；它的稳态落点 **21379.73 Pa** 差 **620.27 Pa**）。
* **天花板只有当最弱入流是壁时才是壁的压力**：臂 `p` 的邻居是牛顿管，于是天花板是 `pA` 本身，那根 8000 Pa 的壁全程最接近drive 只到 **5101.52 Pa** 而从未放开。
* **被 hold 的行静默到最强的一种**：臂 `t` 就是臂 `s` 加上一根永不醒的壁，两者节点在 **2001** 个状态里差 **0.0 Pa**、两根流动壁与泄放管的打印 token 逐字相同；稳态内核不带闩锁也说同一句话，`Ft = Fs =` **22657.9280** Pa。
* **睡眠是数出来的区间**：`Sa` 的 **1355** 个 held 状态 = 第一次放开前 **332** + 重锁到第二次放开 **1023**；18 s 那一跑给它 **7627** 个（其中一次睡 **7295** 步），而 `Wa` 只有开头 **179** 个。
* **日历买 7001 步密封，读出的是渐近线不是日程**：长跑里弱墙的缺口从 **1048.65 Pa**（第 1000 步）只降到 **394.36 Pa**（第 7000 步），每千步收缩 **0.55252 / 0.79305 / 0.91052 / 0.96289 / 0.98491 / 0.99392**——"还没到"在这里变成"到不了"。长跑与 4 s 跑在共享翻步上逐条一致（`Wa` 180、`Sa` 333 与 507），共享的 **1251** 个状态节点数值最大差 **2.00e-07 Pa**（相对 **6.7e-12**）：21 条边与 3 条边求和的次序不同，所以这里对齐的是数而不是 token。
* **提前到达的密封把上行整段换成下行**：臂 `e` 在两根壁都还被 hold 时（0.2→0.4 s）关上，节点随后**向下**走完整 **1051** 步、每步 **-0.2201 ~ -0.2251 Pa**，正是 `dt · Q_泄放 / Cc`（那段里穿过 `1e8` 的均值 **1.11284e-07** m³/s 除以 **1e-9** m³/Pa）；两根壁各只拿到一次放开（`We` **1457**、`Se` **1624**）而丢掉重锁与第二次放开，日程冻结后又汇回与臂 `a` 同一个稳态落点 **21379.7317** Pa（marched 尾 **21375.82 / 21370.01**）。所以 5.53 说"重锁是放开的下游事件"在这里被正面量出：**日程碰不到停着的行，但它能决定节点还有没有第二次机会**。
* **复现与守恒**：**23000** 个状态（7 臂 × 2000 + 9000）由本文件自己的三条行、节点行与采样 `K` 复现——速率最差 **1.360e-09** 相对（只在地速率以上 **1e-8** m³/s 处判）、节点 **7.830e-12**、**零** 个闩锁旗标不同、翻步列表逐步一致；**19658** 个被 hold 的状态全打印单个 token `0`，节点平衡只用打印出的 token 闭合到 **4.951e-14** m³/s。
* **与稳态内核的对应**：稳态不带闩锁却恰好拒付 4 s 窗口从未叫醒的那两根壁——`Sp` 落定 drive **5000.00**、`St` **7342.07** Pa 对阈值 8000，速率打印为 `0`；其余 11 根非零，最小是 `Sa` 的 **1.6701e-07** m³/s 对 **8620.27** Pa。

**成本**：3 个 case 文件（主桌面 9 节点 / 21 根管 / `dt 0.002` × 2000 步，实测 **0.453 s**；长桌面同臂 × 9000 步 **0.254 s**；稳态桌面同拓扑 **0.054 s**），checker 单跑约 **5.6 s**（含 23000 状态的逐状态复现）。

**本增量的边界**：① 分支只做到**多进一出**；一进两出的分流节点与两侧都自由的汇点没有算例，天花板那条律在分叉上没有被问。② 静止支依旧**刚**：阈值以下没有蠕变、没有触变，日程与邻居都不给静止点添加时间依赖（本案要的正是"重新锁上"而不是"慢慢松开"）。③ 只有 **Bingham**（13 根壁都不声明 `flowIndex`），HB × 分支 × 双向闩锁未证。④ 放开与重锁的**步号**分辨率仍是 `dt`（0.002 s），细于此的开关看不见。⑤ 天花板律只覆盖了它的两个端点：入流全为壁（`a/i/t/s`）与"壁 + 牛顿"（`p`，天花板是 `pA`）都被量到，两根壁**加**一根牛顿的中间情形没有。⑥ 层流 only：本案速率远在 5.51 的转捩 locus 以下，past-locus 付多少归 5.51。⑦ `beta`（对流动量账）与静止支同台仍是空白。⑧ 耦合边两笔账照旧拒收（`:143`/`:234`–`:235` 不收阀、`:226`–`:229` 不收阈值），1D↔3D 详细区里的"壁 + 阀 + 分支"仍是空白。⑨ 7001 步的渐近线证的是**这个桌面到不了**，不是"任何壁都到不了"——它是数值证据，不是定理。⑩ 跨内核对的是**同拓扑的稳态内核**与打印 token，不是实测或场数据。

### 5.55 屈服壁第一次挂在节点的**流出侧**：drive 变成 `pF - pB`，于是网络的极限从天花板换成**地板**，重锁律也在同一枚硬币上成立（**14** 根壁、**8** 次翻 = **5** 醒 + **3** 重锁，被 trace 之前预测的那三根逐字命中，其中 `Wb` 的载体是一根**牛顿**管），而**死带是整个解的不动点不是某一行**（`Fd` 用同一个 token **27000.0** 走完 **2001** 个状态，稳态内核那里 `iterations=0 residual=0 converged` 且四个节点打印求解器自己的猜测）——**内核一字未改**

5.54 的边界②写着"分支只做到多进一出，一进两出与两侧都自由的汇点没有算例"，它自己的边界列表还写着另一半：5.46 到 5.54 每一根屈服壁都从定压高压源**流进**自由节点，所以闩锁比较的那个 drive 有一端是**被构造钉住的**。这一格两端同时让位，桌面把 `b`/`c`/`z` 改成"**一根进、两根出**"的分流、把 `d`/`f` 做成"**一根进、一根出但壁在出侧**"、把 `m` 的两端节点都放开（一根带日程的行同时出现在两个节点行里），而 `a` 保留 5.54 的两进一出作为同文件对照。内核不需要改，因为闩锁本来就不分方向：`holdsStatic` 取 drive 的**绝对值**（`transient.cpp:40`–`:43`），`drivePa = pFrom - pTo - ρg dz` 只是拓扑给的顺序（`:52`–`:56`），`Q = 0` 时 `frictionDp` 把 drive 的符号抄回去（`:72`–`:81`），节点行照度数求和（`:395`–`:416`）。七臂 **14** 根屈服壁 + **5** 根牛顿行、**10** 个节点、**19** 条边，每臂的行是

```
ρL/A (Qw - Qw')/dt = (pA - pF) - 4L tw_w/d - Kw ρ|Qw|Qw/(2A²)     入流壁
ρL/A (Qe - Qe')/dt = (pF - pB) - 4L tw_e/d - Ke ρ|Qe|Qe/(2A²)     泄放壁
ρL/A (Qd - Qd')/dt = (pF - pB) - R Qd - Kd(t) ρ|Qd|Qd/(2A²)       带日程的泄放
V/BC (pF - pF')    = Σin Q - Σout Q
```

* **方向不是装饰**：入流壁只能抬升节点，停在自家放开压力 `pA - 4 τ_0 L/d`（**22000 Pa**）从下方逼近；泄放壁是镜像——只能压低节点，停在 `pB + 4 τ_0 L/d` 从**上方**逼近。起点 27000 Pa 让两侧看到不同的数：入流侧 drive **3000**（对 8000，睡着），泄放侧 **7000**（对 8000 睡着、对 6000 和 4000 醒着），所以第一次决策（第 2 个状态起）**9** 根睡、**5** 根醒。
* **闩锁两个方向读同一个数**：**5** 次放开的起步 drive 超出自家 `4 τ_0 L/d` **0.41 – 16.54 Pa**，**3** 次重锁里 `Sa` 差 **3.10**、`Wz` 差 **1.73 Pa** 不到**同一个**数，**8** 次翻的前一步都在另一侧（`|drive - 阈值| < 120 Pa`）——没有迟滞带可拟合。
* **重锁律在流出侧同样成立，而且它要的是"更弱的同一行侧"而不是"另一根壁"**：trace 之前按 `ROWS_AT`（逐节点、逐侧）写出的预测**逐字命中**那三根——`Sa`←`Wa`（6000 在 8000 之下，入流侧）、`Wz`←`Tz`（4000 在 6000 之下，泄放侧）、`Wb`←**牛顿**的 `Db`（根本没有阈值）。其余 **11** 根在自己那一侧没有更弱的行，一根都没有重锁。臂 `c` 是这条律的反面教材：它的节点度数 **3**，可 `Sc` 是**唯一**的供给者，于是它在第 **173** 步醒一次后再没睡回去——买第二方向的是"自己那侧的更弱行"，不是度数。
* **`Wb` 不是单步穿越，这正好是镜像律里载体的那一半**：把它锁住的 drive 比阈值低 **99.04 Pa**，因为牛顿泄放管在 drive 已经退回墙内之后**还在减速液柱**，速率要再过几步才真的归零。
* **地板是渐近线，不是日程**：没人跨过 6000 泄放壁的地板，`Ff` 于是停在地板之上 **303.12 Pa**（速率从未到 0，闩锁从未被问）；`Fz` 停在**最弱**那根（4000）的地板之上 **292.22 Pa** 而不是 6000 那根的 26000 之上——**这正是 `Wz` 能在第 193 步被重新锁上的原因**。18 s 长跑里两块地板仍然只被逼近：`Ff` **475.39 → 82.54 Pa**、`Fz` **585.51 → 60.44 Pa**（第 1000 → 9000 步），臂 `z` 每千步收缩 **0.43708 → 0.73507**（处处 <1）。`Sz` 这类入流壁全程最接近 drive **5707.78** 对 8000，**从未被问**（`Sf`/`Sm`/`Sd` 同）——一串壁可以留其中一根永远睡着。
* **死带是整个解的不动点，不是某一行的**：臂 `d`（8000 进、8000 出）两行在 **2001** 个状态都打印单个 token `0`、旗标从第 2 个状态起 **2000** 个全是 `1`，节点 `Fd` 用**同一个 token 27000.0** 走完全程，因为没有任何行在动，节点行就是 `0 = 0`；那根泄放壁停在自家阈值内 **1000.00 Pa**（drive 恒 7000 对 8000）。
* **分流 × 日程：对称先立后破**：臂 `c` 两根几何相同的泄放管在日程动手前 **401** 个状态逐 token 相同（split ratio 恰为 1，不是拟合），从 `t = 0.802 s` 起 **1600** 个状态**每个**都不同；`closedK = 1e8` 下比值是 **0.0595** 而不是 0（关着的那根仍走 **2.4002e-08** 对开着的 **4.0325e-07** m³/s），因为 `K` 住在动量行里、这一行要付它已经有的惯性；从重开第一步（**2.502 s**）起差距单调收缩（**0.964 → 7.231e-06** 相对，750 个状态），末态 **5.342737e-07 / 5.342776e-07**——剩下的是液柱记忆，不是日程。
* **长跑证明翻步属于拓扑不属于日历**：`Wb` 83、`Wz` 193、`Sb` 317 在 18 s 那跑落在**同一步**（尽管重开晚了 12 s），并且继续睡（`Wb` **8918**、`Wz` **8808** 个 held 状态）；把 `Db` 在 **15.0 s** 重新打开**叫不醒**它锁住的那根壁（第 7500 – 9000 个状态全读 `1`），因为密封正是把节点拉到地板之下才锁的，重开只让它继续往下。
* **复现与守恒**：**38019** 个打印速率（19 行 × 2001）由本文件自己的三条行、节点行与采样 `K` 复现到 **1.555e-09** 相对（只在地速率 ≥ **1e-8** m³/s 处判），**8** 个自由节点压力（含臂 `m` 链那两个）到 **4.967e-12** 相对，**28014** 个闩锁旗标**零**不一致、**14** 根壁的翻步列表逐步一致；**8** 个节点行 × 2000 步只用打印出的 token 闭合到 **4.961e-14** m³/s，也就是末位打印的一个单位（**5.0e-14** = 1e-7 Pa × V/BC/dt），相对失配只长到 **1.889e-03** 且只在 `Fm1` 第 **2** 步（那一步行两侧都是 **8.1e-12 / 8.2e-12** 量级）。
* **闩锁不是对 drive 的猜测，是对上一行速率的后果**：**15748** 个 held 状态对应 **15765** 个打印的 `0`——差的 **14** 个是 `t = 0` 那一列（旗标还没决定就打印 `0`），另有 **3** 个"旗标仍开着却打印 0"的恰是三根壁各自的**钳步**（`Sa` 506、`Wb` 82、`Wz` 192，每根都在自己重锁的前一个状态）；任何壁打印过的最小非零速率是 `Sa` 第 **505** 步的 **2.1608e-13** m³/s。
* **与稳态内核：同一枚硬币的三种答案**：死带三臂单独成文件时 Newton **第一眼看就过**（`iterations=0 residual=0 (converged)`），四个自由节点打印求解器自己的起点 **25000 = (pA + pB)/2**、六根壁 `dp 5000`、`Mm 0`；这 7 行 + 4 节点在 7 臂文件里**逐 token 相同**，所以秩亏是**局部的**——7 臂文件停住（`iterations=200 NOT converged`）是**会流的那两臂**（`a`、`z`）造成的，其残差 **0.0321364485 × 6e-5 = 1.928187e-06** m³/s 正是臂 `a` 的节点失配。预算换 50000：`Fa` **21898.4842 → 21379.7317**（5.54 发布的落点）、`Fz` → **24000.040119**（地板之上 **0.0401 Pa**），残差 **2.6336e-11** 仍高于本文件 tol **1e-12**（banner 比较的就是它，`network.cpp:689`）——所以这是**速率**不是拒绝，而 `Fb 20450.6031`、`Fc 20272.1672` 两个预算**逐 token 不变**。稳态还把 march 末态醒着的 `Ef`/`Em` 也归零（它们在最后一个瞬态状态的 drive 是 **6303.12**/**6520.43** 对 6000）：臂 `f`/`m` 对稳态是死带（天花板 22000 在地板 26000 之下），对 march 是普通流动臂，因为** march 只逼近地板，稳态要落点**。臂 `b` 的稳态落点是串联对 `Sb → 节点 → Db` 两边 **8.84757194168e-07** m³/s（5.53 的串联速率），睡着的壁不参与；臂 `c` 的两根相同泄放管稳态给到最后一位（**5.34398984438e-07** / **5.34398984438e-07**），也就是瞬态在阀开着时的 split。
* **5.54 的两个翻步在这里被穿过一张 19 条边的网络复现**：`Sa` 醒·睡·再醒 **333 / 507 / 1530**，末态节点 **21375.816193 Pa** = 5.54 印为同臂尾巴的那个数。

**成本**：5 个 case 文件（主桌面 10 节点 / 19 根管 / `dt 0.002` × 2000 步，实测 **0.272 s**；长桌面 3 臂 × 9000 步 **0.438 s**；稳态桌面 **0.076 s**；死带桌面 **0.003 s**；`maxIter 50000` 那趟 **22.145 s**），checker 单跑约 **32 s**（其中 22 s 是那 50000 次迭代），**70** 条判据。

**本增量的边界**：① 只有 **Bingham**：**14** 根壁都不声明 `flowIndex`，HB × 泄放侧壁 × 双向闩锁未证。② 静止支依旧**刚**：阈值以下没有蠕变、没有触变；`Wb` 那一支"驱动已退回墙内、速率还要几步才归零"靠的是惯性，不是松弛律。③ 层流 only：瞬态那行的摩阻仍是层流斜率 + 二次 `K`，本案最尖的打印速率是 `Da` 的 **1.03678403182e-05** m³/s（等效 `Re` ≈ **5.61**，远在 5.51 的 locus 之下），past-locus 付多少仍归 5.51。④ 重锁律只量到三种载体各一例（壁←更弱壁、壁←牛顿行）；壁←壁而两者阈值相等（5.54 的臂 `i`）在这里没有对应臂，`Wb` 那 **99.04 Pa** 的缺口也没扫过 `K` 斜坡。⑤ 放开与重锁的**步号**分辨率仍是 `dt`（0.002 s）。⑥ 分流做到一进两出、汇点做到"一根行跨两个自由节点"（臂 `m`），真正的**两侧都自由且各有两根行**的汇没有算例。⑦ 死带"march 把哪个元素停下来"是**一个桌面、一个初值 27000** 的数值证据，初值没扫过，不是定理。⑧ `beta`（对流动量账）× 静止支仍空白。⑨ 耦合边两笔账照旧拒收，1D↔3D 详细区里的"壁 + 阀 + 分支"仍是空白。⑩ 跨内核对的是**同拓扑的稳态内核**与打印 token，不是实测或场数据。

### 5.56 死带被**扫**开：自由节点自己的初值从 `pB` 走到 `pA`，于是"哪个初值会停泊"从 5.55 的一句话变成一条**有宽度和边界的带**（**12** 臂 / **24** 根壁 / **14** 节点，**6** 个停泊点含**两个端点**、**6** 个滑动点，边界两侧各夹到 **1 Pa**，滑动臂的逼近律由内核自己的 Bingham 剖面**预测**而非拟合：`1/gap` 的斜率最差差 **1.4468%**、40 s 落点最差差 **3.619%**）——**内核一字未改**

5.55 的边界⑦承认得很直白：死带让 march 停在哪个元素，只有**一个桌面、一个初值 27000** 的证据，"初值没扫过，不是定理"。这一节就是那次扫描。**同一副两行臂**（`A → F → B`，两根 `τ_0 = 2.0` 的壁，阈值都是 `4 τ_0 L/d = 8000`）复制 **12** 份，每份只改**自由节点自己的初值**，让它从 `pB = 20000` 走到 `pA = 30000`：`20000 / 21000 / 21999 / 22000 / 23000 / 25000 / 27000 / 28000 / 28001 / 29000 / 30000`，外加 `25000` 那份给两根行都挂上关-开日程（臂 `m`）。被扫的量决定于两行代码：`holdsStatic` 把 `|drive|` 和 `4 τ_0 L/d` 比的是 **`<=`**（`transient.cpp:40`–`:43`），`frictionDp` 在 `Q = 0` 只抄 drive 的符号（`:72`–`:81`），于是节点行的两条右端同时归零、`:395`–`:416` 那一行成了 `0 = 0`。

* **带是闭的，而且两端各被夹到 1 Pa**：`[pA - 8000, pB + 8000]` 写成数值是 **`[22000, 28000]`**（宽 **6000 Pa**）。停泊的 **6** 臂在 **2001** 个状态里各打印**自己那一个 token**（`Fd 22000`、`Fe 23000`、`Ff 25000`、`Fg 27000`、`Fh 28000`、`Fm 25000`），两根行各打印单 token `0`；滑动的 **6** 臂 **2001** 个 token **两两不同**，且第 **2** 步就动（`a 20000.12092`、`b 21000.05584`、`k 29999.87908`）。**`d` 和 `h` 正是两个端点**——醒着那根的 drive 恰为 8000 = 阈值，第一次决策就被 `HELD`，所以本案读的是 `<=` 而不是 `<`：换成严格不等号，这两个起点在第 1 步就会流。再往外面一格：`Fc 21999 → 21999.00196` 动、`Fi 28001 → 28000.99804` 动，所以天花板落在 **`(21999, 22000]`**、地板落在 **`[28000, 28001)`**。
* **滑动只逼近，永不越界**：`2001 × 3` 个状态里下方三臂的最大值 **21999.00196 Pa** 距天花板还差 **0.99804 Pa**，上方三臂的最小值 **28000.99804 Pa** 距地板同样差 **0.99804 Pa**；**40 s / 20001 个状态**那一跑把"到不了"数了出来——离自家边最近的是起步只差 1 Pa 的 `c`/`i`，仍有 **0.9807 Pa**，`b`/`j` 差 **49.9313 Pa**，`a`/`k` 差 **51.5302 Pa**。
* **醒过的壁不会再睡，因为速率从来不是恰好的 0**：**2000** 个已决策状态里 **6** 条滑动臂的醒壁全部报 `hold 0`，其打印出的最小非零速率是 **4.78156e-13 m³/s**；同期睡着的那根报 **18 × 2001** 个单 token `0`。**停泊只能靠起点落在带里**，靠不了"流着流着停下来"。
* **上下两侧是同一枚硬币**：`a`/`k`、`b`/`j`、`c`/`i` 三对**逐 token 镜像**——`Fa + Fk = 50000 Pa` 在 **2001** 个状态里最大偏差 **0 Pa**，两侧醒壁的速率序列 **2001/2001** 完全相同。这不是对称性假设，是 `p → pA + pB - p` 把两行方程整体换名的结果，而本案把它按 token 验了一遍。
* **逼近律是预测的，不是拟合的**：近边速率对超出量是**二次**的（`d ln Q / d ln gap = 1.9647 – 1.9974`，**6** 对、跨 gap **0.981 – 51.5 Pa**），因为内核 `q_of_tw` 那支方括号 `1 - 4x/3 + x⁴/3` 在 `x = 1` 处**常数项与一次项同时相消**；把它接进节点行（`Cc = V/BC = 1e-9`）给出 `d(gap)/dt = -k_Q gap²/Cc`，即 `1/gap` 随时间**直线**增长、斜率 `K = k_Q/Cc = 4.9087e-4 Pa⁻¹s⁻¹`。**六条臂全在这条线上**：`1/gap` 对步号的实测斜率最差差 **1.4468%**（`a`/`k`），最贴近预测的是那两条 1 Pa 起步的（`c`/`i` 差 **0.0199%**）；把直线积分到 **40 s** 预测落点，`c`/`i` 差 **0.000%**、`b`/`j` 差 **2.943%**、`a`/`k` 差 **3.619%**。每 **2 s** 的收缩比单调升向 1 而从不越界（`a` **0.22988 → 0.88896**、`b` **0.36267 → 0.89273**、`c` **0.99804 → 0.99807**）——**一条 1/t 尾巴，不是指数**。
* **日历碰不到停泊的行**：臂 `m` 就是臂 `f` 加一份 `K` 一路涨到 **1e8** 的关-开日程，两臂 **2001/2001** 个节点 token、两根行的速率 token、**2000** 个已决策旗标全部相同——一行被 `HOLD` 时它的行里已经没有 `K` 项可读（`x[k]=0; r[k]=0; J[k][k]=1.0`，`transient.cpp:353`–`:355`）。
* **停泊场是节点自己两行的性质，不是矩阵的性质**：只留 3 臂（`c` 21999、`a` 20000、`d` 22000）的 `case_single.dict` 里，`Fd` 依旧单 token **22000**、`Fc` 依旧 **21999 → 21999.0019592**，节点与两行 **2001/2001 逐 token 相同**；只有远端那条 `a` 对不上 token（**1306/2001** 相同），漂移上界是**末位打印的一个量子** **1.000e-07 Pa**（速率 **3.800e-16 m³/s**）——两种矩阵尺寸把同一批行按不同顺序求和，所以这里对的是**数**不是 token。
* **复现与守恒**：**12** 个节点 × **2001** 个状态由本文件自己的两条行 + 节点行复现到最差 **5.978e-12** 相对，**24** 条速率序列在 **1e-8 m³/s** 地板之上复现到 **4.393e-10** 相对，**48000** 个闩锁旗标**零**不一致；**48024** 个打印速率里 **36024** 个恰为 `0`，每个节点行只用打印数字闭合到 **1.189e-06 m³/s**（每步只折合 **1.189e-15 Pa** 的节点运动）。`t = 0` 那一列 **24** 个旗标全打印 `0`（还没决策），所以事件扫描必须从 `i = 1` 起，否则 **6** 条停泊臂各自报一次假释放。
* **稳态那一侧：一个盲猜对十一个起点**：`case_steady.dict` 把 **11** 份串联壁给稳态内核——**22** 根行全部 `Q 0`、`dp 5000`，**11** 个自由节点全部打印 **25000**，banner 是 `iterations=0 residual=0 (converged)`。稳态**没有**闩锁、也就没有带：未定节点被留在两个定压的平均值，而 `(pA - 4τ_0L/d + pB + 4τ_0L/d)/2 = (pA + pB)/2` 这个恒等式说明**只要带非空，稳态的默认起点必然是停泊态**。5.55 的"`Fd` 是整解不动点"在这里第一次有了宽度：**内核的两个打印值（25000 与 converged）恰好落在被扫出来的带的正中间**。

**成本**：4 个 case 文件（主扫描 **14** 节点 / **24** 根管 × 2001 状态，实测 **0.27 s**；长扫描同 **12** 臂 × **20001** 状态 **3.07 s**；3 臂局部性文件 **0.08 s**；稳态 **13** 节点 / **22** 根管 **0.003 s**），checker 单跑 **14.24 s**、ctest `-R` 记 **14.56 s**、ctest `-j8` 记 **26.09 s**（Test #85），**49** 条判据、**0** 次内核改动，全量回归 **111/111** 绿（本轮实测 **354.01 s**）。

**本增量的边界**：① 扫的是**初值**，不是**参数**：**12** 臂共用同一副 `τ_0 / L / d / dt`，带对 `τ_0`、对 `Cc = V/BC`、对 `dt` 的**导数**没测（预测律里 `K ∝ k_Q/Cc` 依赖这三者，本案只在一个点上验了它）。② 只有 **Bingham**：`flowIndex ≠ 1` 时那支方括号在 `x = 1` 处一次项是否仍相消未查——若不相消，近边律就从二次退回一次，`1/t` 尾巴换成对数。③ 每臂只有**一根**入流壁和**一根**泄放壁，5.55 那种"同侧更弱邻居"的重锁载体在这里被构造排除了，所以本案**没有**任何一次重锁可数。④ 停泊点**不可进出**：`[22000, 28000]` 是个只能靠起点选择的口袋，没有任何日程能把一个滑动臂**放**进带里（`m` 就是这条证据），也没有任何东西能把停泊臂**叫**出来——后者只测了 `K` 日程，源/汇日程未测。⑤ 局部性只到**一个**远端臂的精度（**1.000e-07 Pa** = 一个打印量子，**1306/2001** token 相同），"矩阵尺寸无关"是**数**上的而非 token 上的。⑥ 层的**位置**分辨率仍是 `dt`（0.002 s）：`c` 的第 2 步打印 `21999`，靠 `%.10g` 才看得见它在动。⑦ `beta`（对流动量账）× 静止支仍空白；耦合边两笔账照旧拒收。⑧ 跨内核对的仍是**同拓扑稳态**与打印 token，不是实测或场数据。

### 5.57 死带 × Herschel-Bulkley：**带和闩锁不认识指数，认识指数的是带后面的那一行**（**3** 副七臂矩阵 × **2** 个时钟 = **6** 个瞬态文件 + **1** 个稳态文件，**42** 根 HB 壁，指数在**第一次已决策那一步**就被打印速率复现（最差相对 **4.698e-07 / 3.612e-07 / 2.874e-07**），近阈幂 `m = (n+1)/n` 与远阈幂 `1/n` 各归各的地理）——**内核一字未改**

5.56 的边界②留了一句话："`flowIndex ≠ 1` 时那支方括号在 `x = 1` 处一次项是否仍相消未查——若不相消，近边律就从二次退回一次"。这一节把整个停泊扫描搬到 Herschel-Bulkley 壁上：壁律 `τ = τ_0 + K γ̇^n`，三副文件只在 (`flowIndex`, `consistency`) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02) 上不同，`τ_0 = 2.0` 处处一样，每副七份起点 `25000 / 22000 / 28000 / 21999 / 20000 / 28001 / 30000`（`a/b/c` 在带里、`d/f` 在外面一格、`e/g` 在外面两千格），每副跑两个时钟（**2001** 与 **20001** 个状态）。

* **桌面是指数盲的**：`yieldThresholdPa = 4 τ_0 L / d`（`transient.cpp:28`–`:31`）里既没有 `n` 也没有 `K`，于是三档都在 `[22000, 28000]` 里停成**各一个 token** 的 `25000 / 22000 / 28000`，2001 与 20001 个状态都是；端点 `b`/`c` 在**第一次决策**就 `HELD`，醒着那根的 drive 恰为 **8000 = 阈值**，所以 `holdsStatic` 的 `<=`（`:40`–`:43`）同样只读几何。换句话说：**5.56 那条带的宽度、闭端与两个吸引子在 HB 上一字不动**。
* **壁律被打印 token 复现，而且从第一步就复现**：行自己说的是 `(inertia/dt)·ΔQ + 4L/d·τ_w = drive`（`:381`–`:382`），所以 τ_w 可从三个打印 token 反解。每档取 **1312** 个采样速率、**起点是第一次已决策那一步**，最差相对误差 **4.698e-07**（`n=0.5`）/ **3.612e-07**（`n=1`）/ **2.874e-07**（`n=2`）——那是 `%.10g` 的打印量子，不是定律；也没有一个反解出的壁应力掉到 `τ_0 = 2.0` 之下。**`n = 1` 的通用括号就是 5.56 的 Bingham 括号**（**3999** 个壁应力上最差差 **6.3e-16**），所以这一节加的是**同一支公式的第二半**，不是第二支公式。
* **第一步是斜坡，斜坡上括号照样闭合**：醒着的行把驱动力的 **539.5 / 1635.8 / 1821.0 Pa**（占 8000 的 **6.7% / 20.4% / 22.8%**）留给流体自己的惯性，壁应力因此起步于 **2.36512 / 2.09103 / 2.04472 Pa** 而不是同一驱动在静止时的 **2.49999 / 2.49997 / 2.49997**，且这个数**对 `n` 单调**——增稠的壁在"把液柱加速起来"这一笔上花得更多。镜像的两臂 `e`/`g` 落在同一个数上（差 **0.0e+00**）。
* **指数有地理，不是某一个数**：括号自己在屈服点上方 **1 nanoPa** 处读 **3.000000 / 2.000000 / 1.500000** = `(n+1)/n`，在上方 **五千倍**屈服应力处读 **2.000072 / 1.000048 / 0.500029** = `1/n`。所以 `m` 是**近阈**幂而不是唯一的幂：本案的 gap 跨 **0.4 – 2000 Pa**，`δ/τ_0` 始终不到 0.5，远阈渐近线**从没被走到**。真实轨迹只在这两点之间穿行，**18 / 26 / 34** 条弦全部严格落在 `(1/n, m)` 内，且每条与括号在同两个 gap 上的弦互差 **0.027% / 3.922% / 2.023%**。
* **5.56 的相消是 `n = 1` 的特例，不是一条定律**：`n = 1` 时近阈幂 `m = 2`——那正是 5.56 依赖的二次尾巴，`1 - 4x/3 + x⁴/3` 的双重相消是这一个指数的性质。`n ≠ 1` 时一次项**不再**相消，直线因此从 `1/gap` 挪到 **`gap^(-1/n)`**：实测斜率 **1.66666e-08 / 4.90772e-04 / 1.46335e-02** Pa⁻¹s⁻¹，距首项预测 `A_n/n`（`A_n = π R³ (n/(n+1)) (d/4L)^m / (τ_0 K^{1/n} C_c)`）差 **1.859% / 0.021% / 0.011%**，而 `n=1` 那一档就是 5.56 自己的 **4.908739e-04**。`n = 0.5 → 1 → 2` 之间这条直线的**陡度**差 **2.9e4** 与 **30** 倍，带却一步没动。
* **同样在外面一格，买到的东西完全不同**：40 s 后 1 Pa 起点距自家边界 **1.000000 / 0.980747 / 0.397901 Pa**，2000 Pa 起点 **1145.336119 / 51.530219 / 2.774363 Pa**，两种起点的单调序一致（陡尾慢、缓尾快）。`n=0.5` 那档 1 Pa 起点**一点也没走近**（**0.000000 Pa**，20001 个状态只打印 **4** 个 token），近边两臂因此**一条弦也不贡献**——所以 5.56 用来夹边界的"1 Pa 括号"是 `n = 1` 的便利，不是闩锁的性质；而 `n=1`/`n=2` 的近边臂反而把自己的 `m` 读了出来：**1.99976 – 1.99983** 与 **1.49985 – 1.49995**，各差 `m` 的 **0.012% / 0.010%**。
* **其余的停泊性质在 HB 上原样保留**：`d`/`f`、`e`/`g` **逐 token 镜像**（节点压力之和距 50000 差 **0.0e+00**，醒壁速率序列 **2001/2001** 与 **20001/20001** 相同）；**醒过的壁不再睡**（**20000** 个已决策状态里 4 条滑动臂的醒壁全报 `hold 0`，最小打印速率 **8.179690e-18 / 4.720579e-13 / 7.346240e-12** m³/s，另一侧仍在阈值下的行则始终 `HELD` 且单 token `0`）；`t = 0` 那一列 **14** 个旗标全 `0`，事件扫描照旧必须从 `i = 1` 起。**`n=1` 还是 5.56 的数**：4 s 落点 **21999.00196 / 21540.24169**、40 s 距界 **0.980747 / 51.530219 Pa**，与 5.56 为同一起点发布的数差 **0.00e+00 – 3.00e-07 Pa**。
* **稳态内核同样盲**：`case_steady.dict` 把三副 HB 串联壁（**6** 根行）交给稳态内核——全部 `Q 0`、`dp 5000`，**3** 个未定节点在 `n=0.5 / 1 / 2` 三档都打印 **25000**（同一个数），banner `iterations=0`。5.56 那个 `(pA + pB)/2` = 带正中间的恒等式与指数无关，因为阈值本身与指数无关。

**成本**：7 个 case 文件（三副瞬态 **0.96 / 0.44 / 0.37 s**，三副长时钟 **10.31 / 5.45 / 3.79 s**，稳态 **< 0.01 s**），checker 单跑 **26.45 s**、ctest `-R` 记 **26.61 s**、ctest `-j8` 记 **46 s**（Test #86），**106** 条判据、**0** 次内核改动，全量回归 **112/112** 绿（本轮实测 **353.63 s**）。

**本增量的边界**：① 指数是**和 `K` 成对**扫的：三档同时改 (`n`, `K`)，所以"直线陡度差 2.9e4 / 30 倍"要归给这一对，只有 `m = (n+1)/n` 这半边是纯 `n` 的、并由括号自己证明。② `τ_0` × 指数**交叉**未测：三档共用一个 `τ_0 = 2.0`，带的宽度对 `n` 的导数是零这条**只用**在阈值公式里，没有扫 `τ_0`。③ **0 条日程**：本案没有阀、没有局部 `K`、没有源/汇时间表，所以"带是闭的"这条只在无日程的行上验过 HB；5.50 / 5.52 那种日程 × 静止支的同台在 HB 下仍是空白。④ HB × 湍流（Darby-Metzner 类闭合）与重锁律的"等阈值"那一端照旧未做。⑤ 节点行"只用打印 token 闭合"那一步（5.56 的 **1.189e-06** m³/s）没在 HB 矩阵上重做，本案对的是**壁律复现 + 镜像 + 旗标**。⑥ 局部性（3 臂文件）只到 5.56 的 Bingham，HB 下"停泊场是节点自己两行的性质"未重测。⑦ 远阈 `1/n` 只是**括号自己**在 `δ = 5000 τ_0` 处读出来的，本案的任何轨迹都没走到那里。⑧ `beta` × 静止支仍空白；耦合边两笔账照旧拒收。⑨ 跨内核对的仍是**同拓扑稳态**与打印 token，不是实测或场数据。

### 5.58 阀日程 × HB 停泊带同台：**日程改的是速率与账单，不是带**——10 臂 × 3 个指数共用一副 4 s 时钟（open → 关闭中 → 闭合 → 开启中 → open），`h` 那份带满日程的停泊臂对 `a` **逐 token 相同**，`b` 在 `K` 扫到 1e8 时仍在第 1 步被 `HELD`；而同一支阀在同一根壁上花掉的钱，最差到 **46.3%** 的阈值——**内核一字未改**

5.57 把带切成两半：阈值 `4 τ_0 L/d` 与 `holdsStatic` 的 `<=` 不认流变指数，带下游的一切认。它留下的格子是 5.56 的边界④——**日程能不能碰到静止支**。阀日程没有第三条路可走，它就写在同一行里（`transient.cpp:363`–`:382`：惯性 + 摩阻 + `rho·K(t)·|Q|·Q / (2A²)` = drive），而 `:353`–`:355` 在**读 `kv` 之前**就把被 `HOLD` 的行换成 `x[k]=0; r[k]=0; J[k][k]=1.0`。于是两个断言在动手之前就已经写在这两段代码里，本案把它们逐 token 量了一遍，并量出**日程真正能动的那一格有多大**。

* **桌面与闩锁不认识日历**：**10** 副两行臂（`A → F → B`，两根 `τ_0 = 2.0` 的壁、阈值都是 `8000`）× **3** 个指数，`a` 25000 无日程、`h` 25000 把同一份日程挂在**两根行**上：`25000 / 25000 / 22000 / 28000` 四个停泊起点在 **2001** 个状态里**各打印自己那一个 token**，且 **`h` 就是 `a`**——**2001** 个节点 token、**2001+2001** 个速率 token **逐 token 相同**，两根行在 **2000** 个已决策状态里全报 `hold 1`。边臂 `b`（drive 恰为 `8000`）在 `i = 1` 依旧 `HELD`，身后的 `K(t)` 一路扫到 **1e8**——被 `HOLD` 的行里没有 `K` 项可读。**20** 个 `t = 0` 旗标全为 `0`，事件扫描照旧从 `i = 1` 起。
* **日程的第一字节就是日程自己的算术**：`valve 0` + 线性关闭 `0.5→1.0 s`（`closedK 1e8`）+ 线性开启 `2.0→2.5 s`，因为开阀值是 **0**，挂日程的行在 `t < 0.5 s` **逐字就是**没挂日程的那一行——**3** 对滑臂在 **0–250** 步（**2001/2001** 里前 251 个 token）全同；远的 `g`/`i` 对在**三个指数上都在第 251 步**分开，那正是第一个越过 `closeStart` 的采样点（`K = 4.000e+05`）。
* **一个日程，三种延迟**：近边的 `d`/`e` 对（1 Pa 起步）分开在 **None / 757 / 268** 步（`n = 0.5 / 1 / 2`）——在 `m = 3`（即 `n = 0.5`，`m = (n+1)/n` 见 5.57）那一侧，**2001** 个打印状态里日程**一次都没有露面**；`n = 1` 要走到闭合窗第 **506** 步、`n = 2` 只要 **17** 步。同一份日历，延迟由指数决定。
* **重建必须带上阀那笔钱，否则报的是流变而非阀**：从**第一个已决策步**起，挂日程那行的打印速率就是本文件的 HB 截面积分——把 `K(t) rho|Q|Q/2A^2` 从 drive 上搬下来以后最差 **3.757e-10 / 1.001e-09 / 8.235e-10** 相对；**不搬**则是 **6.4e-01 / 1.4e+01 / 1.9e+01**，即错 **0.6～19 倍**。同一份算法跑在没挂阀的孪生行上是同一个数（**3.658e-10 / 9.430e-10 / 6.909e-10**）。
* **闭合窗口值多少，逐指数**：`t = 1.5 s`（`K = 1e8` 保持、斜坡早已走完）速率被压到孪生的 **1/1.59 / 1/7.22 / 1/8.20**，整条轨迹到 `t = 4.0 s` 少走了 **24.1797 / 279.2209 / 253.1707 Pa**——各自那段接近量的 **14.4% / 18.1% / 13.9%**（接近量本身是 **167.8 / 1540.2 / 1822.4 Pa**，差一个量级）。**比例是日历的，距离是指数的。**
* **阀的钱是真的，可它买不到带**：`K rho|Q|Q/2A^2` 那一项的峰值是 **323.0 / 2483.1 / 3701.8 Pa**，即阈值 8000 的 **4.0% / 31.0% / 46.3%**，**随 `n` 单调**（`n = 2` 那臂的峰值落在关闭爬坡**刚起步**的第 **253** 步，`t = 0.506 s`——`K` 才爬到 **1.2e6**，速率却还是关闭前的那个速率，于是 **39.7%** 的 drive 一口气花在阀上；`n = 1` 在第 **254** 步、`K = 1.6e6` 处同理花掉 **26.2%**，而 `n = 0.5` 的峰值才落在爬完的第 **500** 步）；而把这笔钱搬下去之后，挂日程那行**自己 implied 的壁应力**与孪生行的差全程不超过 **3.1822% / 10.2906% / 11.1758%**（`K > 0` 的 **999** 个状态里）——日程改的是**多少 drive 走到壁上**，不是**壁听哪条律**。
* **开阀之后先欠后还，还得快慢仍由指数决定**：`t ≥ 2.5 s`（`K` 恰为 **0**）剩下的 **750** 个状态里，**每一个**的挂日程臂都比孪生流得**更快**（**1.03..1.04 / 1.59..3.57 / 1.93..5.31** 倍）——它被留在离天花板更远处，同一条律于是拿到更大的 drive；欠账 `2.5→4.0 s` 严格单调缩小 **26.3683→24.1797 / 691.7019→279.2209 / 848.1692→253.1707 Pa**，还掉 **8.3% / 59.6% / 70.2%**。`m = 3` 那一侧的债几乎是永久的：**1.5 s** 自由流动只追回 **2.2 Pa**——5.56 那条 `1/t` 尾巴本来就不追人。
* **带仍不进不出，任何指数都不**：`n = 0.5` 近旁差 **1.00000 Pa**、`n = 1` 差 **0.99804 Pa**、`n = 2` 差 **0.89252 Pa**（天花板侧），地板侧镜像；**没有一次重锁**——被阀压到最低的那两个速率是 **8.1797e-18 / 4.8957e-13 / 2.6280e-11 m³/s**，从来不是恰好的 `0`，闩锁抓不到把手。镜像 `e`（21999，阀在 `S` 行）/`f`（28001，阀在 `E` 行）**逐 token 相同**（节点相加对 50000 偏差 **0.0e+00**、速率 **2001/2001** 相同），因为损耗项对 `Q` 是**奇**的。**10** 个自由节点 × **2001** 个状态的节点行只用**打印 token** 闭合到 **4.9163e-14 / 4.9640e-14 / 4.9382e-14 m³/s**（每步只折合 **9.8e-08 Pa** 的节点运动）——阀在一行之内搬钱，没在节点上创造或销毁通量。

**成本**：**3** 个 case 文件（**12** 节点 / **20** 根管 × **2001** 状态，三副各实测约 **0.9 s**），checker 单跑 **3.5 s**、ctest `-R` 记 **3.65 s**、`-j8` 全量里记 **6.91 s**（Test #87），**115** 条判据、**0** 次内核改动；全量 **113/113** 绿 **368.59 s**。

**本增量的边界**：① 日程只测**一种形状**：线性关 + 线性开 + 一个常数 `closedK`；非线性阀特性、快关、水锤那一类未证。② 阀只挂在**壁上**（同一根管既是屈服壁又带阀），"阀挂在别处、经网络间接碰到静止支"归 5.53，且 5.53 只有 **Bingham**。③ 本案构造上**没有一次重锁**可数（唤醒行的速率从来不是恰好的 `0`），"日程能否把一个停泊臂叫出来"仍是 5.56 边界④那一格。④ 时钟只到 **4 s**：**8.3% / 59.6% / 70.2%** 是 `1.5 s` 内的还款比例，不是渐近终点。⑤ 层的**时刻**分辨率仍是 `dt`（0.002 s）：第 **251** 步那个"第一字节"是 `dt` 的整数倍。⑥ 只扫了**三对** `(n, K)`，且 `K` 与 `n` 绑在同一支 HB 律里，`n` 固定、`K` 单独扫的二维表没有做。⑦ 本案**没有**跑稳态内核——阀日程在稳态那一侧没有对应物，跨内核的对照留在 5.56/5.57。⑧ `beta`（对流动量账）× 静止支照旧空白；耦合边两笔账照旧拒收。⑨ 对的是**同拓扑的两条内核与打印 token**，不是实测或场数据。

### 5.59 HB 壁 × 双向闩锁同台：**阈值不认指数，"能不能重新锁上"认**——5.54 那 17 次翻搬到 **15** 根声明了 `flowIndex` 的壁上，醒来在 **318 / 333 / 343** 步、重锁在 `n = 1 / 2` 各发生一次、在 `n = 0.5` **根本没被问起**——**内核一字未改**

5.54 证明了"醒来过的壁能被自己的节点叫回去睡觉"，并把它唯一的流变留白记在边界③：**13 根壁没有一根声明 `flowIndex`**。5.57 又把死带切成两半——`4 τ_0 L/d`（`transient.cpp:28`–`:31`）与 `holdsStatic` 的 `<=`（`:40`–`:43`）不认指数，带**下游**的一切认；5.58 证明日历只在闩锁**之后**碰到那一行。这两条合起来问出一个比"HB 改不改带"更尖的问题：**重锁要的是节点回来，而节点是被壁律推着走的**。于是本案把 5.54 的桌面原样重写成 HB，只看一件事：哪一个指数能拿到那一次翻回去的旗标。

内核里有**两种**静止支形状，只有第二种能让一根**正在动**的壁睡下（`:340`–`:360`）：

* `held(i) = (qPrev == 0) and |drive| <= 4 τ_0 L/d` ——桌面，`:40`–`:43`；
* `crossing(i) = (qPrev != 0) and |drive| <= 4 τ_0 L/d` **且** Newton 迭代把 `x[k]` 带过了零 ——于是 `x[k]=0; r[k]=0; J[k][k]=1`（`:353`–`:355`）。

缺的那一半是壁自己那根管柱的动量，`x[k] * pd.qPrev <= 0.0` 这一行里没有 `n`——但**何时**满足它，全看每条律此前把多少动量放进了柱里。臂 `a` 就是 5.54 的双向闩锁，两根壁都改写为 Herschel-Bulkley：

```
(inertia/dt) dQ_k + 4L/d · tw_k            = pA - pF        （Sa: τ_0 2.0，Wa: 1.5）
(inertia/dt) dQ_d + Kd(t) rho|Q_d|Q_d/2A²  = pF - pB
(V/BC) dpF                                 = Q_Sa + Q_Wa - Q_d
```

`(n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02)` 三个文件，同一副关闭 `0.8→1.0 s`（`closedK 1e8`）+ 开启 `2.5→2.7 s`，节点从 `2.7e4` 起步；再加 `s`（`a` 去掉强行）、`u`（只有强行）、`p`（邻居换成牛顿管）。**阈值、入口闩锁、步内 clamp、度数通用的节点行——这四段代码一字未改。**

* **桌面这一半照旧盲，而且盲两次**：`Wa` 的旗标在**第 1 与 180** 步翻、`u` 那根孤壁在 **1 与 318**、`p` 那根只在 **1**——**三个指数同一个整数**；`Fa` 前 **180** 个状态、`Fu` 前 **318** 个状态在三个文件里打印**同一个 token**，因为那些状态里每一根壁都恰是 `0`，节点只由牛顿排放那一行决定。分开的第一步就是速率不再是零的那一步，而那里分得极小：其后 **40** 个状态三档只差 **8.505 / 3.680 Pa**，即阈值 8000 的 **0.11% / 0.05%**——律的第一次打印是一个小数，不是一条新轨迹。
* **到达不盲**：臂 `a` 那根 8000 Pa 的壁醒在第 **318 / 333 / 343** 步（**0.636 / 0.666 / 0.686 s**，跨 **25** 步 = **0.050 s**），三档打印 **1.4662e-14 / 6.72243e-11 / 1.45108e-10 m³/s**。每一步都是自己**第一个**入口 drive 越过 `4 τ_0 L/d` 的步：`8001.7883 / 8006.1182 / 8000.6927` Pa（前一步 `7991.3683 / 7998.6378 / 7994.4938`），超出量 **0.022% / 0.076% / 0.009%**，三档都在同一个数外 **80 Pa** 以内——阈值读不到里面，到达的时刻却是指数给的，因为 180 步之后把节点腾空的是**弱壁那条 HB 律**。
* **重锁要两半都齐，中间那段滑行由律决定**：`n = 1` 的入口 drive 在第 **507** 步落到 **7996.9049** Pa，旗标**同步**在 507 翻（滞后 **0** 步）；`n = 2` 在**第 465** 步就已低于阈值（**7998.4507** Pa），那一行却还醒着、还打印 **7.28766e-09 m³/s**，旗标要到第 **472** 步才翻——滞后 **7** 步 = **1.4e-02 s**。同一段 C++（`x[k] * pd.qPrev <= 0.0`，`:353`–`:355`）里没有一个 `n`，两档却差七步。
* **`n = 0.5` 不是被拒绝，而是从未被问**：醒来之后入口 drive 的最小值仍是 **8001.7883** Pa，桌面一次都没拿到"低于阈值"的机会；旗标全轨迹只翻两次（`[1, 318]`），其后 **1682** 个状态一律 `hold 0`。原因是那笔补给：封闭窗第一步（第 **500** 步）两根壁供给 **1.30083e-07** m³/s 对关闭阀的 **5.43469e-08**（比值 **2.394**），节点从 **21292.23214** 只爬 **108.2611** Pa，整个 `[500, 1350]` 里最好 **21400.49321** Pa——距这根壁的释放压力 **22000 Pa** 差 **599.507** Pa。`n = 1 / 2` 分别**超出 1192.446 / 1544.178 Pa**（`n = 2` 超得太早，密封还没走完就已经睡了）。弱壁独自峰值 **1.27525e-07** 对 **1.83724e-06 / 2.0945e-06 m³/s**：**一个量级的供给就是"有重锁"与"没有重锁"的差**。
* **账单本身逐指数**：第 500 步供给/泄漏比值 **2.394 / 27.268 / 29.821**，节点爬升 **108.2611 / 1210.4094 / 1391.8863** Pa（"回升在每个指数上都存在，只有它的大小认识 `n`"）；第二次放开在 **1530 / 1568** 步，分别是 `K` 归零（第 1350 步）之后 **0.36 / 0.44 s**，入口 drive **8001.0644 / 8002.4517** Pa——买回这第二次醒来的是重新打开的排放，不是流变。睡眠窗 `[507, 1530)` 与 `[472, 1568)` 各是 **1023 / 1096** 个恰好的 `0`、一个连续块，中间没有第三次决定；`n = 1 / 2` 的完整阶梯是 `[1, 333, 507, 1530]` 与 `[1, 343, 472, 1568]`——剪切胀稠那根壁**晚 10 步醒、早 35 步睡、再晚 38 步醒**，四次全读同一个 8000 Pa。
* **行诚实到用完根为止，clamp 跟在后面一步**：每个**未 held 且有根**的步，打印速率就是本行 implied 壁应力（drive 减去 `(inertia/dt)·Δq`，这一行没有阀项）的 HB 截面积分——`Sa` 最差 **1.957e-08 / 4.561e-07 / 4.429e-08** 相对（第 **319 / 505 / 470** 步，**1682 / 642 / 559** 个已决策状态），`Wa`/`Ws` **3.463e-09 / 4.956e-09 / 2.818e-09**，`Su` **1.539e-08 / 7.103e-09 / 7.153e-09**；`p` 那根壁 **2000** 个已决策状态全报 `hold 1`，括号一次都没被 offered。而**整条轨迹里恰好只有一个**未 held 步 implied 出低于屈服的壁应力：第 **506 / 471** 步，正是 clamp 前那一步（`n = 0.5` 一个都没有）——**行先失去根，`:353`–`:355` 一步之后跟上；"held"不是一个残差**。桌面与打印机也一步不差：唯一打印恰好 `0` 的清醒状态就是 506 / 471（`n = 0.5` 没有，**1682** 个清醒状态里一个都没有），因为速率已掉到打印量子之下（别处最小非零 **9.243e-14 / 2.161e-13 / 1.651e-10** m³/s）。
* **`n = 1, K = mu` 就是 5.54 那面壁，不是长得像**：旗标 **2001** 个状态逐位相等，`Sa` 阶梯 `[1, 333, 507, 1530]`、`Su` 阶梯 `[1, 318]` 两份文件同整数；节点 `Fa` / `Fu` 差 ≤ **1.0e-07** Pa（**1920 / 1643** 个 token 逐字相同，即一个打印小数 = 阈值 **1.2e-09** 的比例），速率最差相对 **7.2e-09**；在壁从不醒来的 `p` 臂上两份文件 **2001/2001** token **逐位相同**——不碰任何律的算术精确复现，只有流动的状态付小数。`p` 的节点最低 **24898.48194** Pa、从不触及 22000：天花板是**邻居那条律**的性质，在任何指数上都成立。
* **被 held 的那半行真的贡献零**：`s` 就是 `a` 去掉强行，直到强行第一次流动（**318 / 333 / 343** 步）弱壁的速率 token 与节点 token **逐字相同**（在 318/333/343 那一步分开）。**12** 条节点行（**4** 自由节点 × **3** 档）只用打印 token 在每一次放开、睡眠、第二次放开上闭合到 **4.9451e-14 m³/s**，即每步 **9.89e-08 Pa** 的节点运动；**15** 个 `t = 0` 旗标全为 `0`，事件扫描照旧从 `i = 1` 起；判据直接读到的 **25** 次旗标翻动（`Sa` **10** = 三档各自 `2/4/4`、`Wa` **6**、`Su` **6**、`Sp` **3**）里只有 **2** 次是重锁、**2** 次是第二次放开，其余是初始 `HELD` 与各壁的第一次醒来；臂 `a` 的阶梯上任意两次决定之间至少相隔 **50** 步（**0.1 s**）——没有 `dt` 尺度的抖动。

**成本**：**3** 个 case 文件（**6** 节点 / **10** 根管 / **5** 根 HB 壁 × **2001** 状态，另加 5.54 那份对照文件），checker 直跑 **4.96 s**、ctest `-R` 记 **4.92 s**、`-j8` 全量里记 **5.43 s**（Test #88），**87** 条判据、**0** 次内核改动；全量 **114/114** 绿 **372.08 s**。

**本增量的边界**：① 只测一副拓扑、一份日程（关闭 `0.8→1.0 s`、开启 `2.5→2.7 s`、4 s 时钟）与非线性阀特性/快关无关。② 只扫**三对** `(n, K)` 且成对扫，`τ_0` × 指数的交叉（5.57 那一格）在这里仍共用 `τ_0 = 2.0 / 1.5` 两个值。③ `n = 0.5` 那一格证的是"**没被问起**"，不是"任何 `n` 都到不了"——它是数值证据，不是定理；把节点初值抬高 600 Pa 就能问，本案没做。④ 步号分辨率仍是 `dt`（0.002 s）：那 7 步滑行是 `dt` 的整数倍，细于此看不见。⑤ "哪个清醒状态打印 `0`"受制于**打印量子**，不是内核的一个判决。⑥ 本案**没有**跑稳态内核，`n = 1` 的跨增量对照对的是 5.54 的瞬态文件。⑦ `beta`（对流动量账）× 静止支、耦合边两笔账照旧拒收。⑧ 对的是**同一内核的两份字典与打印 token**，不是实测或场数据。

### 5.60 τ_0 × flowIndex 交叉扫双向闩锁：**重锁是 `τ_0` 上的一个窗口，不是半根轴**——12 格里 **4** 格睡下、**2** 格醒着却从未被问、**1** 格晚到、**5** 格从不醒；`τ_0` 买的是天花板（**3000 Pa**），`n` 买的是补给（**≈2140 Pa**，比值 **41.6**），而**正余量还不够**——**内核一字未改**

5.57 把死带切成"桌面盲 / 下游认指数"两半，5.59 把双向闩锁整张桌面重写成 HB，但两案都只用了 `τ_0 = 2.0 / 1.5` 这**一对**屈服应力——5.59 边界②自己记着这一格空白。5.59 边界③还留了一句更要紧的话：`n = 0.5` 那格证的是"**没被问起**"，不是"任何 `n` 都到不了"。本案把两条轴交叉起来，在同一张桌面上问：**重锁到底是 `τ_0` 的单调函数，还是一个窗口？** 四臂只改强壁的 `τ_0`：

| 臂 | `τ_0` | 阈值 `4 τ_0 L/d` | 释放压力 `pA − 阈值` |
|---|---|---|---|
| `a` | 1.75 | 7000 Pa | 23000 Pa |
| `b` | 2.00 | 8000 Pa | 22000 Pa |
| `c` | 2.25 | 9000 Pa | 21000 Pa |
| `d` | 2.50 | 10000 Pa | 20000 Pa |

弱壁固定 `τ_0 = 1.5`（6000 Pa），`(n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02)` 三副文件 ⇒ **12** 格。每臂仍是 5.59 的三行，只有 `τ_s` 在动：

```
(inertia/dt) dQ_k + 4L/d · tw_k            = pA - pF        （S: τ_0 1.75/2.0/2.25/2.5，W: 1.5）
(inertia/dt) dQ_d + Kd(t) rho|Q_d|Q_d/2A²  = pF - pB
(V/BC) dpF                                 = Q_S + Q_W - Q_d
```

同一副关闭 `0.8→1.0 s`（`closedK 1e8`）+ 开启 `2.5→2.7 s`，节点从 `2.7e4` 起步，`dt 0.002`、`endTime 4.0`（**2001** 状态）。**阈值（`:28`–`:31`）、入口闩锁（`:40`–`:43`）、步内 clamp（`:353`–`:355`）、度数通用的节点行（`:395`–`:416`）——四段代码一字未改。**

* **桌面第三次盲，这次盲的是邻居的屈服应力**：**12** 根弱壁（三档 × 四臂）的旗标全部翻在第 **1** 与 **180** 步，与自家强邻居的 `τ_0` 无关；而且四臂在**前 239 / 240 / 241 个状态**里是**同一批字符串**——弱壁速率与节点压力 alike，因为那些状态里每一根强壁都还被 held，而 held 的行就是 `x=0, r=0, J=1`。前缀结束的步号**恰好等于最早那次强壁醒来**（`[239, 318, 1479, 2001]` / `[240, 333, 2001, 2001]` / `[241, 343, 2001, 2001]`），节点 token 从不更早分开。
* **阈值是几何，唤醒时刻是律**：四根壁的阈值 **7000 / 8000 / 9000 / 10000 Pa** 纯由 `4 τ_0 L/d` 给出，里面没有一个 `n`；每一次醒来都发生在**第一个**入口 drive 越过自家阈值的步上，超出量 `n = 0.5` 为 **5.6287 / 1.7883 / 3.2571 Pa**（阈值的 **0.080% / 0.022% / 0.036%**）、`n = 1` 为 **7.0834 / 6.1182 Pa**（**0.101% / 0.076%**）、`n = 2` 为 **10.2370 / 0.6927 Pa**（**0.146% / 0.009%**）。同一个几何数，到达它的步号由指数给。
* **12 格普查：窗口，不是半轴**——`n = 0.5`：`a [1, 239]` **AWAKE**、`b [1, 318]` **AWAKE**、`c [1, 1479]` **LATE**、`d [1]` **NEVER**；`n = 1`：`a [1, 240, 963, 1403]`、`b [1, 333, 507, 1530]` **SLEPT**，`c`/`d` **NEVER**；`n = 2`：`a [1, 241, 709, 1439]`、`b [1, 343, 472, 1568]` **SLEPT**，`c`/`d` **NEVER**。合计 **4 SLEPT / 2 AWAKE / 1 LATE / 5 NEVER**。重锁在**两侧**都丢：`τ_0 ≥ 2.25` 那侧是**从不醒**（`n = 1 / 2` 的最大入口 drive 只有 **8674.0542 / 8577.5570 Pa**，距 9000 差 **325.9458 / 422.4430 Pa**；`n = 0.5` 的 `d` 最近，**9802.5079 Pa**，仍差 **197.4921 Pa**），`n = 0.5` 那侧是**从未被问**。一个从不被问起的参数在打印里**不存在**：`n = 1 / 2` 的 `c` 与 `d`（`τ_0` 2.25 对 2.5）**2001/2001** 个 token 逐字相同，`n = 0.5` 的 `c`/`d` 相同到第 **1479** 步。
* **分解是定量的：`τ_0` 买天花板，`n` 买补给**：`τ_0` 从 1.75 扫到 2.5 把释放压力挪了 **3000 Pa**，却只把封闭窗 `[500, 1350]` 的最好节点挪了 **51.36260 / 45.20672 / 44.18058 Pa**（`n = 0.5 / 1 / 2`）；反过来，固定释放压力沿 `n` 扫，最好节点挪 **2138.60589 / 2143.68512 / 2145.78792 / 2145.78792 Pa**（四臂），四个跨度彼此只差 **7.18202 Pa**——**强壁的屈服应力改的是节点必须到达哪里，不是律能送到哪里**。比值 **41.6**（`n` 沿臂 `a` 的 2138.60589 对 `τ_0` 在 `n = 0.5` 的 51.36260）。最好节点还**单调下降**于强壁 `τ_0`（**21448.40242 > 21400.49321 > 21397.03982 ≥ 21397.03982**，另两档同形）：更强的邻居供给更少，所以补给更少，尽管它自己要的也更低。第 **500** 步的供给/泄漏比值 **2.931 / 2.394 / 2.348**（`n = 0.5`）、**29.515 / 27.268 / 27.453**（`n = 1`）、**32.899 / 29.821 / 30.076**（`n = 2`）——`τ_0` 只在这个比值的小数第二位上说话，`n` 说的是量级。
* **正余量是必要的，不是充分的**：`n = 0.5` 臂 `c` 的封闭窗最好节点 **21397.03982 Pa** 比它自己的释放压力 **21000 Pa** 高出 **397.040 Pa**，却**仍然一次都没睡**——那 **851** 个封闭状态**全部**打印 `hold 1`，因为它醒在第 **1479** 步，比密封破掉（第 **1350** 步）**晚 129 步**。桌面要有一根**正在动**的壁才有东西可 hold；余量再大，窗口里没有醒着的行也是零。对照 `n = 0.5` 的 `a`/`b`：余量 **−1551.598 / −599.507 Pa**，醒来后入口 drive 最小 **7005.6287 / 8001.7883 Pa**，仍高出阈值 **5.6287 / 1.7883 Pa**——那两格是"够不着"，`c` 这格是"够着了但没人在家"。
* **滑行同时认两条轴**：`n = 1` 两臂的入口 drive 落到阈值下的那一步（**963** / **507**，**6999.2383 / 7996.9049 Pa**）旗标**同步**翻，滞后 **0 / 0** 步，且那一步打印的速率恰为 `0`；`n = 2` 臂 `a` 第 **704** 步就已低于阈值（**6999.9119 Pa**）却还打印 **1.42139e-09 m³/s**，旗标 **5** 步后（**709**）才翻，臂 `b` 第 **465** 步（**7998.4507 Pa**）打印 **7.28766e-09 m³/s**、**7** 步后（**472**）才翻。**同一段 `x[k] * pd.qPrev <= 0.0` 里既没有 `n` 也没有 `τ_0`，四格给出四个滞后**（0 / 0 / 5 / 7）。睡眠窗 `[963, 1403)` / `[507, 1530)` / `[709, 1439)` / `[472, 1568)` 各是 **440 / 1023 / 730 / 1096** 个恰好的 `0`、一个连续块；第二次放开在 **1403 / 1530 / 1439 / 1568** 步，即 `K` 归零（1350）之后 **53 / 180 / 89 / 218** 步（**0.106 / 0.360 / 0.178 / 0.436 s**），入口 drive **7007.0388 / 8001.0644 / 7005.9030 / 8002.4517 Pa**——买回它的是重新打开的排放，不是流变，也不是屈服应力。
* **行诚实到用完根为止，clamp 一步之后跟上**：每个未 held 且有根的步，打印速率就是本行 implied 壁应力的 HB 截面积分——`Sa` 最差 **8.500e-09 / 1.876e-06 / 8.512e-08** 相对（第 **241 / 961 / 707** 步，**1761 / 1318 / 1027** 个已决策状态）、`Sb` **1.968e-08 / 4.402e-07 / 8.066e-08**（**319 / 505 / 470**，**1682 / 642 / 559**）、`Sc` **4.255e-08**（**1481**，**521**，只在 `n = 0.5` 有）、弱壁 **3.463e-09 / 4.956e-09 / 2.818e-09**（各 **1820** 个状态）。**5** 格（`n = 0.5` 的 `d`、`n = 1 / 2` 的 `c`/`d`）**2000** 个已决策状态全 `hold 1`，括号一次都没被 offered。而 implied 出低于屈服壁应力的未 held 步，在四个 SLEPT 格里**恰好各一个**：第 **962 / 506 / 708 / 471** 步，正是各自 clamp 的前一步；其余八格**一个都没有**。
* **臂 `b` 就是 5.59 的臂 `a`，不是长得像**：三档的旗标 **2001** 个状态逐位相等，阶梯同为 `[1, 318]` / `[1, 333, 507, 1530]` / `[1, 343, 472, 1568]`；速率 token 逐字相同的个数 `Sb` **865 / 1684 / 1655**、`Wb` **1643 / 734 / 615**、`Db` **1536 / 1391 / 624**（共 **2001**），最差相对差 `Sb` **6.062e-10 / 1.780e-08 / 1.336e-08**、`Wb` **2.957e-10 / 1.425e-10 / 2.530e-10**、`Db` **1.333e-11 / 9.721e-12 / 1.094e-11**；节点 `Fb` 对 5.59 的 `Fa` 差 ≤ **1.000e-07 Pa**（**1989 / 1958 / 1663** 个 token 逐字相同）。同一段物理换了一份文件，付的小数只来自 **12** 根管对 **10** 根管的求和次序。
* **记账**：**12** 条节点行（4 自由节点 × 3 档）只用**打印 token** 闭合到 **4.9661e-14 / 4.8836e-14 / 4.9100e-14 m³/s**，即每步 **9.93e-08 / 9.77e-08 / 9.82e-08 Pa** 的节点运动——`τ_0` 扫描在行之间搬钱，没在节点上创造或销毁通量。**24** 个 `t = 0` 旗标（每副 8 个）全为 `0`，事件扫描照旧从 `i = 1` 起（否则会凭空多出 24 次翻动）；强壁阶梯上任意两次决定至少相隔 **238 / 174 / 129** 步（**0.476 / 0.348 / 0.258 s**，各档 **3 / 6 / 6** 个间隔）——12 格里没有 `dt` 尺度的抖动。

**成本**：**3** 个 case 文件（**6** 节点 / **12** 根管 / **8** 根 HB 壁 × **2001** 状态 × **4** 臂，另加 5.59 那三副对照文件），六次 app 调用合计 **9.41 s**（首次 **4.06 s** 含冷缓存），checker 直跑 **10.93 s**、ctest `-R` 记 **10.61 s**、`-j8` 全量里记 **23.35 s**（Test #89），**146** 条判据、**0** 次内核改动；全量 **115/115** 绿 **933.34 s**。

**本增量的边界**：① 只测**一副拓扑、一份日程**（关闭 `0.8→1.0 s`、开启 `2.5→2.7 s`、4 s 时钟），与非线性阀特性/快关/水锤无关。② `τ_0` 只取 **0.25 网格上的四个值**：窗口的两条边都**没有定位**——`2.0` 睡、`2.25` 从不醒之间的那根临界 `τ_0` 未扫，`n = 0.5` 侧"从未被问"到"被问"的临界指数同样未扫。③ `(n, K)` 仍是**三对且成对**扫，`n` 固定、`K` 单独扫的二维表没有做。④ `n = 0.5` 那两格证的是"**没被问起**"，不是定理；把节点初值抬 **600 Pa** 就能问（5.59 边界③），本案没做。⑤ **LATE** 这个分类依赖本案选的封闭窗 `[500, 1350]`：换一份日程，`n = 0.5` 臂 `c` 可能落进另一类。⑥ 步号分辨率仍是 `dt`（0.002 s）：那 **5 / 7** 步滑行是 `dt` 的整数倍，细于此看不见。⑦ "哪个清醒状态打印 `0`"受制于**打印量子**，不是内核的判决。⑧ 本案**没有**跑稳态内核，跨增量对照对的是 5.59 的瞬态文件。⑨ `beta`（对流动量账）× 静止支、耦合边两笔账照旧拒收。⑩ 对的是**同一内核的两份字典与打印 token**，不是实测或场数据。

### 5.61 源日程唤醒停泊臂：**阀日程的形状第一次出现在源上**——瞬态 `sources{}` 学会 `schedule{ openRate; close { start; end; } open { start; end; } }`，与管上 `valveClose / valveOpen` **同一份 ramp helper**（`sourceRateAt` 复用 `valveKAt`），在**同一隐式新时刻 `t + dt`** 求值；注入脉冲把本来一直 held 的 HB 壁叫醒一次——**释放步相同，重锁步被日程推后**

> 把 §5.56 边界④（"源 / 汇日程能否叫醒一根停泊臂"）拿出来做掉。新内核钩子（transient.cpp:125-:145 加 `sourceRateAt`，:329-:337 每步 reset+sum 进 `nd.source`，parser 加 `schedule{}` 子块；header 加 `SourceEntry` 与 `TransientSource` 上的 ramp 字段）把瞬态 `sources{}` 从**常数**升级到**分段线性 ramp**：`0 → openRate` over `[close.start, close.end]`，hold，`openRate → 0` over `[open.start, open.end]`。无日程时退化为常数 `rate`（与 §5.52-§5.60 逐 token 相同）；parser 拒绝**既无 `rate` 也无 `schedule`** 的 entry（非零退出码 + 消息）。节点衡算 `sumIn` 与 `commitStep` 的液态库存都看见同一个日程率（`transient.cpp:400` / `:480`），`sourceFraction` 保持为常数。

**拓扑**：1 根 HB 屈服壁（`Sa`：`L 20, d 0.02, (n, K, τ_0) = (0.5, 1.0, 2.0)`，阈值 `4 τ_0 L/d = 8000 Pa`，释放压力 `pA − threshold = 22000 Pa`）+ 1 根泄放管（`Da`：`valve 0`，`valveClose 0.8→1.0 s closedK 1e8`、`valveOpen 2.5→2.7 s`）+ 1 个**自由节点**（`F`：`V 0.1, B 1e8, p_start 2.7e4`）。两副字典只差一个 `sources` 块：`case_off.dict` 无源（控制），`case_on.dict` 加 `inj{ node F; schedule { openRate 2e-6; close { start 1.2; end 1.4; } open { start 2.0; end 2.2; } } }`——注入 `0 → 2e-6 m³/s` 在步骤 **600..700**，hold 到步骤 **1000**，ramp 回 `0` 在步骤 **1000..1100**。所有常量沿用 §5.59 / §5.60（`ρ 850, μ 0.1, dt 0.002, endTime 4.0` 即 **2001** 状态）。

**行 per arm**（日程不进入行，只进入 `sumIn`）：

```
ρ L/A (Q_s − Q_s')/dt = (pA − pF) − 4 L tw_s / d           (Sa, HB)
ρ L/A (Q_d − Q_d')/dt = (pF − pB) − R Q_d − K_d(t) ρ|Q_d|Q_d/(2 A²)   (Da)
V/B (pF − pF')        = Q_s + Q_inj(t) − Q_d                (node F)
```

**内核**一行未改：**阈值** `4 τ_0 L/d`（`transient.cpp:28-:31`）、**entry desk** `holdsStatic = qPrev == 0 && |drive| ≤ threshold`（`:40-:43`）、**in-step clamp**（`:353-:355`，在 `kv` 被读之前跑）、**度泛型节点行**（`:395-:416`）——这**四条**正是 §5.59 / §5.60 用的同四条；日程盲，因为**它们不看 `nd.source`**。新加的只在**求值侧**：每条 `SourceEntry` 的 rate 在每步 `solveStep` 开始时**重置然后求和**到 owning 节点（`:329-:337`），用与阀同一个 `t + dt` 求值（保证 Newton 残差里的 `sumIn` 连续）。

* **释放步相同，重锁步被日程动**：`case_off` 的旗标阶梯 `[1, 318]`（events()，**2** 次翻：`t = 0` 边界 + 第一次释放），从步骤 **318** 到结束**一直流**；`case_on` 阶梯 `[1, 318, 531, 1466]`（**4** 次翻：`t = 0` 边界 + 释放 + 日程推动的重锁 + 日程 ramp-off 后的再释放）。前 **600** 个状态（日程为 `0`）的 `p_F`、`Q_Sa`、`Q_Da`、旗标 `_Sa` **逐 token 相同**——日程**不能**倒推回到它还没动的时刻。
* **阈值与释放压力照旧是 8000 Pa / 22000 Pa**：两副字典启动 banner 都打 `pipe Sa static-hold threshold 8000`；决定释放那一步的驱动 `pA − pPrev` 在 `case_off` 为 **8006.5905 Pa**（> 阈值），再前一步 **7999.3816 Pa**（≤ 阈值）；`case_on` 逐 token 相同。日程**不能**改阈值，因为阈值是 `τ_0` 的函数，不是 `nd.source` 的函数。
* **日程推动重锁的物理**：步骤 **318..530**（`0.636..1.060 s`）节点从 `22000` 跌到 `21239`（壁在流、泄放在流、两者把节点夹在释放与下游之间）；步骤 **600..700**（`1.2..1.4 s`）日程把 `Q_inj` 从 `0` 抬到 `2e-6`，`dp_F/dt = (Q_s + Q_inj − Q_d)/C_c` 变正，节点开始**回升**；步骤 **700..1000**（`1.4..2.0 s`）`Q_inj = 2e-6` 恒定，节点继续上行直到穿过 `22000`，步骤 **531**（`1.062 s`，实际在 `1.2 s` 之后的积分累积）壁被**重新锁上**；步骤 **1000..1100**（`2.0..2.2 s`）日程 ramp off，节点失去补给开始回落，步骤 **1466**（`2.932 s`）再次穿过 `22000`，壁**再次释放**。日程**只动了两个整数**（`531`、`1466`），**没动释放本身**。
* **HB 截面积分照旧**：每个 `Q_Sa > 0` 的步，`Q_Sa` 就是本行隐含壁应力的 HB 截面积分（`q_hb(tw, τ_0, n, K)`），最差相对 **1.4e-8**（步骤 **319**，**1682** 个已决策状态）。held 步 `Q_Sa == 0` 到打印位。
* **节点衡算照旧**：`case_off` 的节点闭合 `V/B · dp_F/dt = Q_s − Q_d` 到最差 **9.9e-8 Pa/step**（**2000** 个已决策状态）；`case_on` 的 `V/B · dp_F/dt = Q_s + Q_inj(t) − Q_d` 到最差 **1.2e-7 Pa/step**，注入贡献 `Q_inj(t)` 用**同一份 `sourceRateAt` 在 checker 里再求一次**——日程在 checker 与内核是**同一份 ramp 表达式**（`transient.cpp:125-:145` verbatim）。
* **泄放阀与日程无关**：`Da` 的 `K(t)` 由 `valveKAt` 求值，`case_off` 与 `case_on` 逐 token 相同；源日程**只**影响节点压力，泄放管只**看见**那份压力。
* **源日程** `schedule{}` 块与 `valveClose / valveOpen` 块**形态相同**：`close` 块把率从 `0` ramp 到 `openRate`（"阀 close 把 K 从 openK 推到 closedK"的类比），`open` 块把率从 `openRate` ramp 回 `0`（"阀 open 把 K 从 closedK 拉回 openK"的类比）；ramp 在 `t` 上**连续**，跨 ramp 边界的样本**不**看见跳变。parser 强制 `close.start ≤ close.end`、`open.start ≥ close.end`、`open.start ≤ open.end`，否则非零退出码。

**成本**：**2** 个 case 文件（**3** 节点 / **2** 根管 / **1** 根 HB 壁 × **2001** 状态 × **2** 副），两次 app 调用合计 **0.11 s**，checker 直跑 **0.38 s**、ctest `-R` 记 **2.67 s**、`-j8` 全量里记 **2.23 s**（Test #90），**6423** 条判据、**1** 次内核改动（parser + 求值侧；行算术**一行未改**）；全量 **116/116** 绿 **1311.16 s**。

**本增量的边界**：① 只测**一副拓扑、一份源日程**（`openRate 2e-6`、`1.2→1.4 s` 上 ramp、`2.0→2.2 s` 下 ramp），与源的正弦/脉冲/多段日程无关。② 源只挂在**一个**自由节点上；多源共享同一节点时的 fraction 累计**未**扫（parser 接受，但多源节点的 `sourceFraction` 是最后一个源的值——沿用旧行为，未修正）。③ `schedule{}` 块只支持 `close` + `open` 两段线性 ramp，**不**支持点序（`times / values` 表）；要的话走 `TwoFluidPipe.inlet.schedule`（已经支持）。④ 源日程与阀日程**没有**放在同一根管上做过；本案源在节点、阀在另一根管的端点。⑤ `beta`（对流动量账）× 静止支、耦合边两笔账照旧拒收。⑥ 耦合边拒收流变与屈服（`transient.cpp:226-:229`）照旧。⑦ 对的是**同一内核的两份字典与打印 token**，不是实测或场数据。

### 5.62 停泊带的下边缘：**`pB` 往下挪 500 Pa，`n = 0.5` 的臂 `d` 第一次醒来**——停泊带的下边缘是渐近地板 `pB` 与释放天花板 `pA − 4 τ_0 L/d` 的**较小者**；`pB = 20000` 时臂 `d`（`τ_0 = 2.5`，阈值 `10000 Pa`，释放压力 `20000 Pa`）的节点**只能趋近**释放压力而**永远到不了**，`pB = 19500` 时地板**被移走**，节点第一次穿过去——**内核一字未改**

5.60 的 `NEVER` 分类里臂 `d` 是"阈值 `10000 Pa`、释放压力 `pA − 4 τ_0 L/d = 20000 Pa`、max entry drive `9802.5 Pa`，差 `197.5 Pa`"；5.60 的边界④说"`n = 0.5` 臂 `d` **从没被问**，但**没有**说"问不到"是定理还是本案的参数。本案回答：**不是定理**。节点压力的渐近地板是 `pB`（泄放管另一端的定压），不是释放压力；`pB = 20000` 时地板恰与天花板相切所以节点**只能从上方逼近**而**永远停在上面**（5.55 的"死带"就是这种"两个吸引子之间的不动点"），把 `pB` 往下挪**任意正数**就给节点**穿过**天花板的空间——问题只是"穿过要多久"。

**拓扑**：与 5.60 完全相同（四臂 `a / b / c / d`、每臂一根强壁 `S{a,b,c,d}` + 一根弱壁 `W{a,b,c,d}` + 一根泄放管 `D{a,b,c,d}`），**唯一的差别**是 `B { fixedPressure 1.95e4 }`（5.60 是 `2.0e4`）。三档 `(n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02)`，每副四个臂同时跑，与 5.60 共享同一份内核、同一份 checker 框架。

**`pB` 扫描**（`n = 0.5` 臂 `d`，`τ_0 = 2.5`，阈值 `10000 Pa`，释放压力 `20000 Pa`）：

| `pB` (Pa) | 臂 d 流动状态数（共 2000 已决策步） | 分类 |
|---|---|---|
| 20000 | 0 | NEVER（5.60 的原始结果） |
| 19800 | 0 | NEVER |
| 19780 | 4 | 刚穿过 |
| 19760 | 61 | — |
| 19740 | 103 | — |
| 19720 | 136 | — |
| 19700 | 163 | — |
| 19500 | 313 | 本案采用值 |

转捩发生在 `pB ≈ 19780`（4 个流动状态），即天花板下方约 `220 Pa`；每多给 `20 Pa` 约多 `40..50` 个流动状态（线性区间），到 `pB = 19500` 共 `313` 个——仍不到总步数的 `16%`，说明**节点穿过天花板的速度很慢**（泄放管把节点往下拉的速率受 `Q_D` 限制，`Q_D` 又被 `pF − pB` 限幅）。

**`n = 0.5` 的三件事**（case_5.dict）：

* **臂 `d` 第一次醒来**：`Sd` 在第 **1689** 步释放（`3.378 s`），`313` 个流动状态直到结束。节点 `Fd` 在第 **2000** 步停在 **19752.84 Pa**——**低于**释放天花板 `20000 Pa` **247.16 Pa**，这在 `pB = 20000` 时**构造上不可能**（节点不会比 `pB` 更低）。释放那一步的 entry drive **10003.26 Pa** > 阈值 **10000**，再前一步 **9999.14 Pa** ≤ 阈值——**差 4.12 Pa**，与 5.60 其他臂的释放余量同一量级。
* **臂 `c` 的重锁回来了**：5.60 的臂 `c`（`τ_0 = 2.25`，阈值 `9000`，释放 `21000`）只在第 **1479** 步释放一次、**没有**重锁（因为节点从未回到 `21000` 以上）；本案 `pB = 19500` 把节点拉得更低，臂 `c` 在 **388** 步释放（比 5.60 的 1479 **早 1091 步**）、**929** 步重锁、**1390** 步再次释放——**完整的双向闩锁周期回来了**。重锁那一步的 entry drive **8997.22 Pa** ≤ 阈值 **9000**，再前一步 **9001.83 Pa** > 阈值——差 **2.78 Pa**，与释放侧的余量对称。
* **弱壁照旧同步**：四臂的弱壁 `Wa / Wb / Wc / Wd` 都在第 **1** 步释放（`[1, 1, 1, 1]`），与 5.60 / 5.59 / 5.58 / 5.57 逐字相同——**下边缘移动对弱壁没有影响**，因为弱壁的释放由"起始 drive 超过阈值 6000 Pa"决定，那是节点 `27000` 与源 `30000` 之间的距离，与 `pB` 无关。

**`n = 1` 与 `n = 2`**（case_1.dict / case_2.dict）：臂 `c / d` 仍然 `NEVER`——`n = 1` 的补给不足以把节点拉到 `20000` 以下（max entry drive 臂 `d` **9802.51 Pa**，与 5.60 逐 token 相同到小数第 4 位），`n = 2` 同理（max entry drive **9588.94 Pa**）。下边缘移动**只**在 `n = 0.5` 买到了新格子。

**节点闭合**照旧：`V/B · dp_F/dt = Q_S + Q_W − Q_D` 到最差 **4.97e-14 m³/s**（**12** 条节点行 × **2000** 已决策步），与 5.60 的闭合数同一量级；`t = 0` 旗标全 `0`（**24** 个）。

**成本**：**3** 个 case 文件（**4** 臂 × **3** 档 `(n,K)` × **2001** 状态），三次 app 调用合计 **0.72 s**，checker 直跑 **3.14 s**、ctest `-R` 记 **5.52 s**（Test #91），**84** 条判据、**0** 内核行。

**本增量的边界**：① `pB` 只取**两个值**（`20000` 来自 5.60，`19500` 是本案），转捩点 `pB ≈ 19780` 由 6 次手扫定位、精度 `±20 Pa`；临界值的精确位置（如二分法到 `1 Pa`）**未**做。② `n = 0.5` 的臂 `c` 重锁周期只出现在 `pB = 19500`；`pB = 19800` 时臂 `c` 仍只释放一次（节点最低 `20713 Pa`，差释放 `21000` 还有 `287 Pa`），`pB` 在 `19500..19800` 之间**有没有**一个值让臂 `c` 也走完重锁周期**未**扫。③ 节点初值仍是 `27000 Pa`；5.60 边界③ 说的"节点初值抬 600 Pa 使 `n = 0.5` 臂被第一时间问"那一格**未**做（那是另一个独立轴）。④ `beta` × 静止支、耦合边两笔账照旧拒收。⑤ 本案**没有**跑稳态内核。⑥ 对的是**同一内核的两份字典与打印 token**，不是实测或场数据。

### 5.63 耦合边的流体第一次进入三维详细区：**1D 行算术一字未动，3D 扩散项换成逐面表观黏度**——声明但牛顿档与未声明档的求解日志**逐位相同**，`dp` 对同一本构的平面渠道闭式落到 **1.0203 ~ 1.0752**

**这一格从哪来**：5.45 / 5.46 / 5.47 把幂律、Bingham、Herschel-Bulkley 三支壁面律写进 1D 稳态内核，5.48 把瞬态行接到同一个 `pipeFrictionPa`；四条解析期拒收里有一条是"`coupledResistance` 的耦合边不许带流变"，`docs/gap_analysis.md` 把它记为 **A4 / C1**，估 **2.5 增量**（"3D 侧需要完整非牛顿黏度回路，与 5.45-5.47 的 1D 侧同量级"）。本案关闭它，**实测 1 增量**。

**为什么不是"把 `pipeFrictionPa` 用完整 rheology key 调用"**（原 gap 清单的落法）：耦合边的损失**就是**伙伴的陈述，`pipeLossAt` 在 `network.cpp` 的第一行就 `return coupledResistance * q`，瞬态动量行同样只装 `R·q`；再算一次 1D 摩阻是 §5.43 那条**重复计费**。而 §5.49 的 `holdsStatic` / `frictionDp` / 穿零闩锁又都会读 `yieldStress`——只在 parser 里删掉四个 fatal guard 会让耦合边凭空长出一条"屈服壁"，把本来属于 3D 区域的阈值记到 1D 行上。所以落法是**把声明搬走而不是放开**：`NetworkPipe::partnerRheology` 收下 `(n, K, τ_0)`，`flowIndex / consistency / yieldStress` 三个壁面律字段清零；1D 侧一个物理读取点都不用改，既有测试**结构上不可能被移动**（`verify.coupled1d3d` 实测逐位照旧）。

**成本估计为什么错了一个量级**：3D 侧**没有**改 SIMPLE。§5.23 的 RANS 伙伴块早已把逐面黏度做成求解器钩子 `IEddyViscosity::nuEffFace()`，`simpleSolver.cpp:97` 每趟用 `fvm::laplacianFace` 读它一次。新增一个流体只需实现这个接口：

| 文件 | 作用 |
|---|---|
| `src/common/include/cf/common/Rheology.h` | POD 本构 `τ = τ_0 + K γ̇^n` 与 `apparentViscosityPaS()`；三个模块共用一份声明 |
| `src/3d/include/cf/fvm/GeneralizedNewtonian.h`<br>`src/3d/src/solvers/generalizedNewtonian.cpp` | `GeneralizedNewtonianViscosity : IEddyViscosity`：`γ̇ = √(2S:S)` 取单元值，`nuEff = [K γ̇^(n-1) + τ_0(1−e^(−mγ̇))/γ̇]/ρ`，`harmonicFaceGamma` 插到面，**wall 面覆盖**为 `|U_c,∥|/y_wall` 处的黏度 |
| `src/coupling/include/cf/coupling/CoupledNetwork.h`<br>`src/coupling/src/coupledNetwork.cpp` | `ICoupledElement::rheology()`；`DetailBlock` 收 `(n,K,τ_0)` 与数值控制；`lossAt` 把非牛顿模型作为那个钩子交给 `solveIncompressibleSimple` |
| `src/1d/src/network.cpp` / `transient.cpp` | 删 4 处 guard，声明搬进 `partnerRheology`，壁面律字段清零 |
| 4 个 `src/coupling/apps/*App.cpp` | 逐边 `readPartnerRheology(caseDict, edge)` 转发；打印 `nuWall / gdotWall / tauWall` 供 checker 独立重算 |

壁面覆盖是物理需要，不是数值修饰：剪切稀化流体在**壁面** `γ̇` 最大、`ν_eff` 最小，把单元中心黏度外推到壁面会系统性高估阻力。离散动量方程实际收的正是 `ρ·ν_f·U_c/y_wall`，所以用**同一个** `γ̇_w` 定值才让 `τ_w` 与报出的压降闭合。

**算例**（`tests/cases/coupledRheology/case.dict` + `tests/check_coupled_rheology.py`）：5.20 的稳态拓扑（`Src → A → {X(3D) ‖ By} → B → Snk`）、0.5 × 0.02 × 0.02 m 平面渠道、`(50 20 1)` 网格，只在 X 这根耦合边上声明流体；checker 在临时目录把声明改写成五臂（未声明 / 声明但牛顿 / 剪切稀化 / 剪切增稠 / HB 屈服），再用 `transientCoupled1d3d` 的字典跑第六臂。**306 条判据、0 失败**。

| 臂 | `(n, K, τ_0)` | `q` (m³/s) | `dp` (Pa) | 同一本构的平面渠道闭式 | `dp/dp_closed` |
|---|---|---|---|---|---|
| `plain` | — | 1.577703868e-4 | 620.445938 | 591.638951 | **1.04869** |
| `newton` | (1, 0.1, 0) | **1.577703868e-4** | **620.445938** | 591.638951 | **1.04869** |
| `thin` | (0.5, 1.0, 0) | 1.556433553e-4 | 636.466090 | 623.785531 | **1.02033** |
| `thick` | (1.5, 0.00913, 0) | 1.649709723e-4 | 566.091800 | 526.520390 | **1.07516** |
| `hb` | (0.8, 0.3, 0.5) | 1.436554338e-4 | 726.910256 | 711.827917 | **1.02119** |

闭式不是拟合：`q = 2W·K^(−1/n)·G^(−2)·[ n/(2n+1)·T^((2n+1)/n) + n/(n+1)·τ_0·T^((n+1)/n) ]`，`T = G·H/2 − τ_0`、`G = dp/L`，第二项就是塞流区的积分；`τ_0 = 0` 逐项塌回 `2Wn/(2n+1)·(G/K)^(1/n)·(H/2)^((2n+1)/n)`。checker 在 `G` 上二分反解，五个臂各自的比值都在 `[0.95, 1.13]` 内；牛顿臂那 **4.87%** 是入口段与离散误差，不是新误差——它就是 5.20 那条算例本来就有的地板。

**八族判据**：① `plain` 与 `newton` 的全部求解日志（去掉只有声明流体才打印的诊断行）**逐位相同**，尽管只有 `newton` 走 `laplacianFace` 路径；`nuWall` 逐位是 `μ/ρ = 1.25e-4`、`tauWall` 逐位是 `μ·γ̇_w`。② 从日志的 `nuWall / gdotWall / tauWall` 独立重建 `K γ̇^(n−1)/ρ` 与 `K γ̇^n + τ_0`（6% / 10% 内，平均值之积的相关性误差）。③ `|probeSlope|·H/2 == τ_w` 到 **9%**（壁面平均含入口段，充分发展窗口只取 `0.12..0.49`）。④ 上表闭式。⑤ 每趟界面 `residual == |R·q − dpBlock|/dpBlock` 到打印精度（3 seeds × 6~7 趟），三 seeds 落同一点（`spread < 1e-4`），末趟残差 ≤ `1e-4` 容差。⑥ 1D 侧 `Re` 逐位仍是 `ρvd/μ`（`79.268509937331771`），`pipe final X dp` 与并联支 `By dp` 逐位同为 `636.4812793667661 Pa`——**声明没有在 1D 行里被收第二次费**。⑦ 瞬态臂 60 个时间层每步 `dpNetwork == R·q`，末态落在稳态 `thin` 档 2% 内，块报出的 `nuWall` 与稳态档 6% 内。⑧ 四条拒收仍在：缺 `K`、负 `n`、负 `τ_0`（由本构模型自身命名原因）、以及 `turbulence{} + rheology` 争同一个逐面黏度钩子。

**既有拒收判据的迁移**：A4 解除的不只是四处 guard。5.45 / 5.46 / 5.48 / 5.49 四个既有 checker 各有一个 `case_rejectCoupled*.dict`，断言"耦合边带流变 / 屈服必须非零退出"——这正是 A4 的字面内容，所以第一次全量回归跑出 4 条失败（97% / 117）。按新契约它们改为**接受 + 逐位不变**：checker 在临时目录把同一份 case 删掉 `flowIndex` / `consistency` / `yieldStress` 三个 key 生成孪生档，断言两份 `rc == 0` 且 stdout 字符串完全相同；声明只被搬进 `partnerRheology`，1D 行里既没多出摩阻项，也没多出静止支。四个 case 的头注释同步改写（文件名保留 `reject` 只为兼容旧引用），四个 checker 的"判据族"说明也相应改写。`check_transient_power_law.py` 在模块作用域把 `re` 当广义雷诺数用，新引入的正则模块因此别名成 `re_lib`。这是本增量唯一一处既有测试判据改动，`ctest -R` 复跑 **4/4 通过**。

**本增量的边界**：① `turbulence{}` 与非牛顿流体**互斥**，构造里显式 fatal；湍流 + 广义牛顿同时用需要把 `nut` 与 `ν_eff` 合成一个面场，未做。② `detail{ generalizedNewtonian{ shearFloor 20; regularization 100; relaxation 0.15; } }` 是数值地板而非拟合常数：剪切稀化流体在 `γ̇ → 0` 时黏度发散，缺省 `1e-6` 会让第一趟 SIMPLE 拿到比牛顿大 500× 的黏度并实测发散（`iter 10` 起 `continuity 1.9e10`）；屈服档在缺省松弛 `0.7` 下 300 sweeps 内不收敛，降到 `0.15` 后 **10.7 s** 收敛。③ 多边同台带流变：四个 apps 都已转发，但**未**端到端验证。④ PVT 表 / 黑油 / 两流体分布式仍拒收（gap A5）；耦合边与 `valve` / `pump` 的互斥也照旧。⑤ 3D 详细区仍是平面渠道 + `flat` 入口；入口剖面 × 流变、环空（A6）未做。⑥ 闭式比的是**充分发展**解，实测比值里的 2~7% 是入口段与网格，不是系数拟合。

**成本**：1 个新 case 文件、1 个新 checker、3 个新源文件 / 头（common `Rheology.h` + 3D `GeneralizedNewtonian.h/.cpp`）、24 个既有文件的局部改动（含 4 个旧负例 case 头注释与 4 个旧 checker 的判据迁移）；5 稳态臂约 **63 s**、瞬态臂约 **43 s**、4 个拒收 case 可忽略；`verify.coupledRheology`（Test #49）直跑约 **110 s**，`-j8` 全量里两次分别 **93.63 s / 160.54 s**（并发争用）。**最终全量 `ctest -j8`（`CHFLOW_WITH_PYTHON=OFF`）：117/117 绿、375.04 s**；四个迁移过的旧测试 0.23 / 0.34 / 3.98 / 152.57 s 全通过。详见 `docs/gap_analysis.md` §6 实施日志。

### 5.64 全网络每条边都低于自家阈值：**0 次 Newton 迭代、`residual = 0` 的死不动点**——那个零是自由节点的 `0 = 0` 而不是"网络死了"；在自由节点注入 1e-6 m³/s 后 Newton 跑 **19** 步、E 升到 **33978.9475076 Pa**，恰是 Buckingham 闭式反解

**这一格从哪来**：§5.49 边界⑧ 只记了一句「HP 帽让 Jacobian 非奇异」，从没立过判据；`docs/gap_analysis.md` 把它编成 **B6**，状态 `OPEN`、零内核改动、成本 0.5 增量。工作区里曾留下一份自称"已由 checker 验证 / ctest 注册"的 `case.dict` + `probe.dict`，实际上两个文件都未被 CMake 或任何脚本引用，头注释还把 Jacobian 自由节点行写成"全零"——**与 `network.cpp:519` 的 Hagen-Poiseuille 电导帽相反**。本案把这份草稿补成正式算例，并同时修掉这两处错误陈述。

**为什么是 0 次迭代**：`solveSteadyNetwork` 把自由节点初始化为**全部定压节点的平均值**（`network.cpp:538-547`，本例 `(30000+20000)/2 = 25000 Pa`），随后在 `network.cpp:633` 先 `evaluate()` 一次；所有管都 `Q = 0`、无外部速率，故自由节点残差逐位是 0，`relativeResidual = 0 / 1e-300 = 0`，`network.cpp:638` 的 `while (relativeResidual > tol && iter < maxIter)` 根本不进入，`iterations = 0`、`converged`。死支的电导在 `:519` 回落到 `laminarConductance`，所以 **Jacobian 自由节点行不是零行、矩阵非奇异**；它只是没被要求求逆。

| 臂 | 改法（同一份 dict 在临时目录派生） | 实测 |
|---|---|---|
| `allDead` | 原案：5 根 Bingham 管，`4τ₀L/d` = 12000 / 16000 / 20000 Pa | 全部 `Q = v = Re = f = 0`，`iterations=0 residual=0 converged`，E = 25000 |
| `topology` | `probe.dict`：4 根管（去掉 Stk3） | 同一个死不动点，E = 25000，`iterations=0` |
| `shifted` | p(A) 3.0e4 → 3.1e4 | E = **25500**（新平均值）、并行 dp 11000、串联 5500，仍 `iterations=0` |
| `knifeEdge` | Stk1 τ₀ 3.0 → 2.5，使 `dpStatic == dp == 10000` | Q **仍逐位为 0**：放开判据是 `dp > dpStatic` 的严格不等式 |
| `oneFlowing` | Stk1 τ₀ 3.0 → 2.0，阈值 8000 < 10000 | Stk1 `Q = 1.37182879207e-6`（Buckingham 闭式同值），其余四管仍 0，**`iterations` 仍是 0** |
| `injected` | 在 E 注入 1e-6 m³/s、A 取出同量 | `iterations = 19`、`residual = 4.23516473627e-15`；E = **33978.9475076**、仅 Sr1b 流出 1e-6、其 `dp = 13978.9475076` |

**闭式不是拟合**：`dpStatic = 4τ₀L/d`，`Q = Q_HP·g(ξ)`，`Q_HP = πd⁴dp/(128μL)`，`ξ = dpStatic/dp`，`g(ξ) = 1 − 4ξ/3 + ξ⁴/3`（§5.46）。checker 从 dict 自己读出 `L / d / τ₀ / 定压` 重算阈值、各管的受力压降与 `pInit`；`oneFlowing` 正向用闭式核对 `Q` 与 `f·Re·g(ξ) = 64`；`injected` 在 `dp` 上**二分反解** `Q(dp) = 1e-6` 得 `dp = 13978.94750765034 Pa`、`E = p(B) + dp = 33978.94750765034 Pa`，与打印值在 1e-9 相对内相同。注入臂因此不是"看它动没动"，而是把升高的自由节点压力读回本构自己的反函数；Sr1a 只被反向推 3978.9 Pa，仍低于 12000 Pa 阈值而保持 `Q = 0`。

**五族判据、45 条、0 失败**：① 从 dict 重算「每条边的 `4τ₀L/d > 被问的 dp`」，而不是记住一张表；② 三个全死档的 `iterations=0 / residual=0 / 全字段恰为 0 / 各管 dp 恰为受力压降 / E 恰为定压均值`，并用 `shifted` 排除"E 硬编码 25000"；③ 阈值相等刀刃仍 `Q = 0`；④ `oneFlowing` 证明 **0 迭代 ≠ 网络死了**（一根管在流、Newton 仍未进场，因为自由节点只连着两根死管），且该管的闭式与 `f·Re·g(ξ)` 同时成立；⑤ 注入臂的 19 步 Newton 不发散、唯一流动管与源速率相等、`dp/E` 为闭式反解、反向管仍被阈值锁住——这是 HP 帽那条记录第一次变成可失败的判据。

**本增量的边界**：① 全死时自由节点压力**不可辨识**——方程是 `0 = 0`，任意 p 都满足，报出的 25000 / 25500 只是**初始化值**，不是解出来的物理量；不要把它当"凝胶管线的静止压力分布"。② 0 迭代只在**零外部速率**下成立，注入臂立刻让循环跑起来。③ 死支的 HP 帽保证矩阵可逆，但**不等于**该网络能凭空造出流量：注入臂是靠把 E 抬到 33978.9 Pa 让 Sr1b 越过阈值才找到解，这正是 A2 说的"重启压力"缺席的证据，不是亚阈值蠕变。④ 本次只验 Bingham（`n = 1`）；HB 阈值与 n 无关（§5.47）但未在本拓扑上单列一臂。⑤ 稳态专用，瞬态行的 `holdsStatic` 是另一套判据（§5.49），不在这里。

**成本**：**0 内核改动**；2 个已入库 case dict（此前是未跟踪草稿）+ 1 个 checker + 1 条 `add_test`，45 判据；直跑约 **0.1 s**，`ctest -R` 单跑 **0.44 s**，`-j8` 全量里 **0.88 s**（Test #45）。**最终全量 `ctest -j8`（`CHFLOW_WITH_PYTHON=OFF`）：118/118 绿、714.41 s**；对照上一轮 117/117 / 375.04 s，本次机器负载较高，`verify.vofOscillating` 186.64 s、`verify.vofFit` 124.48 s。详见 `docs/gap_analysis.md` §6 实施日志。

### 5.65 HB 壁 × 无日程阀 K：**阀的分成由壁律的指数决定**——同一个 `K = 500` 从 `n = 0.5` 那根管的 10000 Pa 里拿走 **2.33 Pa**，从 `n = 1.5` 那根管里拿走 **308.76 Pa**（133×）；持停管在五个数量级的 K 下逐 token 相同

**这一格从哪来**：`(阀, 日程, 指数)` 表还剩一格没问过——§5.50 量的是 **Bingham 壁 + 无日程阀 K**，§5.58–5.62 量的是 **HB/屈服壁 + 带日程阀**；`docs/gap_analysis.md` 把「HB × 阀 K（无日程）」编成 **C2**，估 0.5 增量、内核改动"与 §5.50 相同路径，只是 key 换成 HB"。实测**零内核改动**：parser 自 §5.47 / §5.49 / §5.50 起就同时收一根管上的 `flowIndex` / `consistency` / `yieldStress` / `valve`，没有任何 guard 拒收这个组合——它只是从未被问过。这也让 C2 的成本估计从 0.5 降到实际 ~0.3（2 个 case 文件 + 1 个 checker + 1 条 `add_test`）。

**同一根管上两项损失**：`ρL/A·(Q_new − Q_prev)/dt = drive − dp_wall(Q_new) − K ρ|Q|Q/(2A²)`，`dp_wall` 是 §5.47 剖面积分的反函数 `Q = q_hb(tw, n, K_cons)`。三项各自的边界在 HB 下全部照旧：静止时阀为 0（`dp_valve(0) = 0`）、阈值 `4τ₀L/d` 不含 n / K_cons / K、阀只在有速率可计时收费。

| 臂组 | `(n, K_cons, τ₀)` | 阀 K | 实测 |
|---|---|---|---|
| `StuckRef/K500/K1e8` | (0.5, 1.0, 2.6)，阈值 10400 Pa | 0 / 500 / 1e8 | 2001 个状态**逐 token 相同**、全为 `Q 0` |
| `RelHB05` 家族 | (0.5, 1.0, 1.0)，阈值 4000 Pa | 0 / 100 / 500 / 5000 | 末速 1.04050 / 1.04031 / 1.03958 / 1.03147 e-6，**严格单调**、各自落在自己的两损失根上 |
| `RelHB15(/K500)` | (1.5, 0.02, 1.0)，阈值 4000 Pa | 0 / 500 | 末速 1.2544e-5（12.1×），K=500 的阀分成 **308.76 Pa 对 2.33 Pa** |
| `HangHB05(/K500)` | (0.5, 1.0, 0.5)，elevation +1 | — / 500 | drive 1664.35 Pa < 2000 Pa：逐 token 相同、恒 `Q 0` |
| `RollHB05(/K500)` | 同上，elevation −1 | — / 500 | drive 18335.65 Pa：各自落在含重力的两损失根上 |
| `StopRef/K5000/K1e8` | (0.5, 1.0, 1.5)，initialFlow | 0 / 5000 / 1e8 | 分别在第 **9 / 9 / 5** 步停住，递推逐状态复现三条迹 |

**指数才是这一格的内容**：阀项是二次的 `∝ Q²`，而 Q 由壁律定——同一根几何管、同一个 `τ₀`、同一个 `K = 500`，剪切增稠壁（n = 1.5, K_cons = 0.02）流 1.2544e-5 m³/s，阀拿走 308.76 Pa；剪切稀化壁（n = 0.5, K_cons = 1.0）只流 1.0405e-6 m³/s，阀只拿走 2.33 Pa。§5.50 的 Bingham 壁只有一条律，这个分成无从谈起。

**七族判据、56 条、0 失败**：① `(flowIndex, consistency)` 逐位读回 + 阈值 `4τ₀L/d` 复算，且 `RelHB05` 与 `RelHB15` 打印**同一个**阈值；② 持停管在任何 K 下逐 token 相同（稳态核对无 `dpValve` 字段）；③ 11 根流动 / 停止管 × 2001 状态全部被「HB 剖面积分 + K ρ|Q|Q/2A²」的隐式 Euler 递推复现，最差 **4.347e-10**（`StopK1e8` 第 4 步）；④ 每根流动管的末速等于两损失根（checker 自己的 tw 空间二分），K 家族严格单调；⑤ 阀分成 133× 的指数效应（上一段）；⑥ 三条停止迹连同停止步被复现，闩锁在停住后的下一步接管；⑦ 稳态核对与行进行在同一个根上会合，打印的 `dp − ρg·dz − dpValve` 经 §5.47 反演读回打印的速率。

**本增量的边界**：① 无日程是**声明的缺席**而非新分支——`valveClose`/`valveOpen` 与本格的组合由 §5.58–5.62 覆盖，本格补的是"常数 K"；② 阀 K 仍按几何项计（`K ρ|Q|Q/2A²`），不随流变变稠变稀——阀有自己的 K，流体有自己的 K_cons，两个 K 不混；③ 只验正向流（Hang 档证明反向即死）；④ 稳态核对里 `dp` 打印的只是压降部分，重力项要自己减回去（`dp − ρg·dz` 才是摩擦半区）——这是打印约定，不是新物理。

**成本**：0 内核改动；2 个 case 文件 + 1 个 checker + 1 条 `add_test`，56 判据；transient app 16 管 × 2001 状态直跑约 **0.6 s**，checker 全程约 **9 s**。全量回归数字见 `docs/gap_analysis.md` §6 实施日志。

### 5.66 `reportEvery` 与 `endTime` 解耦：步数是 `ceil` 不是 `round`，末步裁到余量、`t` 用步索引落点、末态恒打印——endTime 不是 dt 整倍数的 case 第一次存在，而整倍数的 case 逐位不变

**这一格从哪来**：§5.49 定下的 1D 瞬态循环是 `round(endTime/dt)` 步、固定 `dt`、只在 `i % reportEvery == 0` 采样——它悄悄要求 endTime 是 dt 的整倍数，也悄悄回避了"求解器自己的 `dt` 能不能跑在打印步之下"这个问题。`docs/gap_analysis.md` 编为 **C4**，估 0.25 增量。审计（2026-11-11）已记两处部分解耦：1D app 的 `|| i == nSteps` 末态打印、3D PISO 的末步 `dt = min(dt, endTime − t)` 裁剪（§5.34）；本项把二者补进 1D 核心的步进循环（`transient.cpp runTransient` + `transientApp.cpp`），parser 零改动。

**新簿记**：

- 步数 `nSteps = ceil(endTime/dt − 1e-9)`——容差是给近整除陷阱的：`0.12/0.03` 在 IEEE 双精度下是 **4.000000000000001**，没有容差就会跑出 5 步、末步裁成 0；
- 末步 `h = min(dt, endTime − (i−1)·dt)`，隐式欧拉行把它当作自己的 `dt` 读（`I_h = ρL/(A·h)`），所以裁剪改的是**簿记**不是物理；
- `tPrev = min((i−1)·dt, endTime)` 用**步索引**定位，不用累加器——末态的 `t` 精确落在 endTime 上（0.1 就是 0.1，不是 0.09999999999999999）；
- 采样条件 `i % reportEvery == 0 || i == nSteps`——`reportEvery` 不再拥有末态；
- `dt`/`endTime` 必须同时为正，否则 `cfFatal`（核心）/ 报错退出（app）。

**case 与判据**：`tests/cases/transientTimeDecoupled`（dt 0.03、endTime 0.1，即 **0.03/0.03/0.03/0.01** 四步；三管 Wide d=0.04 Newtonian、Narrow d=0.008 Newtonian、Yield τ₀=1.0 Bingham，阈值 4000 Pa < drive 10000 Pa）。checker `check_transient_time_decoupled.py` 六臂 + 拒收、**40 判据**：

| 臂 | 设定 | 断言 |
|---|---|---|
| `clamped` | endTime 0.1, report 1 | steps=4、5 态、t = 0/0.03/0.06/0.09/**0.1**（≤1e-12）、全 15 个报告态被逐步 `h` 的隐式欧拉行复现（Newtonian 闭式 `Q_new = (I_h·Q_prev + drive)/(R + I_h)`，Bingham 在 tw 空间二分；最差 **3.27e-11**） |
| `report3` | report 3 | 4 步运行只出 3 态（步 0、3、末 4），末态 t = 0.1 |
| `report5` | report 5 | 2 态（步 0 + 末态）——**4 % 5 ≠ 0 而末态仍在**，这就是解耦本身；且末态与 `clamped` 的第 4 步**逐位相同**（打印是同一次 march 的窗口） |
| `exact` | endTime 0.12 | 4 步**恰好全为 dt**——近整除容差挡住 4.000000000000001 |
| `compat` | endTime 0.12, report 2 | 3 态 t = 0/0.06/0.12，与旧固定步循环的采样集相同 |
| `shortRun` | endTime 0.012 | **单步** h = 0.012 的裁剪运行，末态等于闭式 `Q_inf(1−a)` |
| 拒收 | endTime 0 | 非零退出 + 指名消息 |

**物理没有动**：第 5 族把 `clamped`/`exact`/`shortRun` 的 Newtonian 末态与整条 march 的闭式积 `Q_inf·(1 − Π a_i)` 对到 **1e-9**——6.40198e-05 / 7.5015e-05 / 8.6268e-06 各自落在自己的积上。裁剪只改簿记，一步不多、一步不少。

**向后兼容**：endTime 为 dt 整倍数的旧 case 逐位不变——`transientYieldValve` 仍 2001 状态、13 管 × 每 14 行、`steps=2000`、末态 t=4。

**checker 自己的坑（值得记）**：首跑 21 条失败的根因不在内核也不在模型，而在 checker 的键提取——case 头注释里写着示例字样（`endTime 0.1`、`reportEvery 3 and 5`），而 `re.search(r"endTime ([0-9.]+)")` 扫的是**整个文件**，第一个命中是注释。教训：**从 dict 里读键必须先剥 `//` 注释再正则**；同理步长序列只能用容差比较（`0.1 − 0.09` 是 `0.010000000000000002`，`0.12 − 0.09` 是 `0.030000000000000027`），不能对字面量。

**本增量的边界**：① 这是**簿记解耦**，不是自适应时间步（D5 仍 OPEN）——`dt` 还是那一个 `dt`，只是末步允许更短；② 末步更短意味着 Newtonian 管的松弛因子 `a_i = I_h/(R+I_h)` 更接近 1，数值稳定性不受影响（隐式欧拉无条件稳定），但精度随 h 变化——1e-9 的闭式积判据实测不受扰动；③ 3D 侧 PISO 的末步裁剪（§5.34）语义与此一致但实现独立，本项不动 3D。

**成本**：核心 `transient.cpp` 循环 ~10 行 + app 同构 + 1 case 文件 + 1 checker + 1 条 `add_test`，40 判据；checker 全程 **0.75 s**（Test #83）。全量回归数字见 `docs/gap_analysis.md` §6 实施日志。

### 5.67 三张逐边表 × 逐边 `profile` 同台：**quasi + slab + charge 三表叠一边**——quasi 边的行用 slab 的 `M` 付 `M·dq/dt`、同时让出 `dpConv`、grade 恒 0，三个声明互不越界

**这一格从哪来**：耦合语法的每一格都被单独问过——§5.40/5.41 把 `ddtFrom` + `inertiaFrom` + 逐边 `profile` 放进一份 case，§5.43/5.44 把 `inertiaFrom` + `convectiveToRow` + `profile` 放进另一份——但**没有任何一份 case 同时携带三张表**。`docs/gap_analysis.md` 编为 **B5**，估 1 增量、内核改动"无（三张表两两已交叉）"。实测**零内核改动**，实际成本 ~0.5：三张表都在时间循环之前安装、互不知晓，合在一起只是从未被写过。

**case**：`tests/cases/transientCoupledMultiTurbulentTripleTable`，拓扑/网格/探针/`dt`/seeds（1.4e7 / 1.36e7）全部复用 §5.41——每个臂都有已发布的孪生。5 臂：

| 臂 | `ddtFrom` | `inertiaFrom` | `convectiveToRow` | 出处 |
|---|---|---|---|---|
| `base` | Euler, Euler | network, network | report, report | §5.41 `p_net` |
| `mix` | Euler, **steady** | network, **slab** | report, report | §5.41 `p_mix` |
| `mix_charge` | Euler, **steady** | network, **slab** | report, **charge** | **新格** |
| `pin` | **steady**, Euler | network, **slab** | report, report | §5.41 `p_pin` |
| `both_chr` | Euler, Euler | network, **slab** | report, **charge** | §5.44 `s_chr` |

全部臂 `profile parabolic`（逐边字段，保持不动）。唯一的新格 `mix_charge` 把三张表叠在同一根 quasi 边上——运行回答了「quasi 边能不能同时被 charge」：能。它的行用 slab 的 `M = ρL_probe/A = 2.1875e6 kg/m4` 付 `M·dq/dt`，`rowBeta = dpConv ≠ 0` 在全部 8 层成立（最小 283.67 Pa），而 `rowError` 恒 0（quasi 边没有可评的账）。

**21 判据、0 失败**（checker `check_transient_coupled_multi_turbulent_triple_table.py`）：① 三张表 + 形状逐臂逐边从日志读回（scheme 从每个状态的 `unsteady` 标志、inertia/profile 从臂行、charge 从 `rowBeta`），臂序与声明一致；② 每个 rung 在 80 个状态上装它声明的闭式 `M`（drawing `ρL/A`、slab `ρL_probe/A`）；③ 24 个 quasi 态是真 quasi（`rowError 0`、`dpInertia 0`、`dpFriction == dpBlock`）；④ 新格三表同边成立 + `mix` 孪生（差一个词）`rowBeta 0` 且 `dpBlock/dpConv/dpFriction/dpInertia/M` 逐位相同——**charge 只动一行的一项**；⑤ charge 不及邻居（`mix_charge` 的 X1 与 `mix` 的 X1 17 字段 × 8 层逐位相同）；⑥ 80 状态 row 闭合 `dpNetwork + M·dq/dt + rowBeta` **slack 恰 0**、`rowBeta` 账户逐态跟随表；⑦ chord 在容差内（最差 8.9123e-4 对 1e-3）；⑧ `base` 1e-12 复现 §5.41 p_net 的 level-1/level-8 状态与 pass 序列 `[3,6,5,5,4,4,4,4]`；⑨ X2 速率钉在 6.0e-5（1e-12，`pin` 的 level-1 分流差 1 ulp 是 quasi 自由边对 split 的正常扰动）；⑩ 支路闭合 6.78e-21；⑪ charge 不多买 block 功（mix 与 mix_charge 同为 72 次、base 与 both_chr 同为 70；pin 的 34 是 quasi X1 让谈判更早结束）。

**本增量的边界**：① 三张表**安装期互不知晓**是设计而非巧合——`ddtFrom` 定 scheme、`inertiaFrom` 定 M、`convectiveToRow` 定行内 `rowBeta`，三者分别写行的不同字段，所以叠加不需要新的交互代码；② quasi 边 + slab + charge 的组合此前从未运行，本案证明它可行但**不证明它最优化**——quasi 边 grade 恒 0，charged 的 dpConv 只进网络行不进谈判弦；③ 全部臂同形状（parabolic），形状轴在本案是**保持项**不是变量——形状 × 三表的完整交叉仍属各原型小节；④ 运行成本 ~5.5 分钟（5 臂 × 8 层 × RANS SIMPLE），是全量回归里最贵的单项之一。

**成本**：0 内核改动；1 case 文件 + 1 checker + 1 条 `add_test`，21 判据；checker 复跑（含求解）约 **5.5 分钟**。全量回归数字见 `docs/gap_analysis.md` §6 实施日志。

### 5.68 CFL 自适应时间步：**调度是物理的性质**——`dt_next = clamp(cfl·min_k(L_k A_k/|q_k|), dtMin, dtMax)` 从刚提交的流量读出，checker 把整条调度从状态列重推出来，而不是记住一张表

**这一格从哪来**：§5.49 以来每个瞬态 case 都 marched 一个固定 `dt`；§5.66 解耦了 `reportEvery` 与 `endTime` 但步长序列仍是常量。`docs/gap_analysis.md` 编为 **D5**，估 1.5 增量、备注"在最慢管上放一个 CFL 探针"。实测 ~0.5：`TransientNetwork` 里几何都在手边，规则本身 `cflStep` 是 ~15 行。

**语法**：`time { dt <初始>; endTime; reportEvery; cfl <target>; dtMin <下限>; dtMax <上限>; }`——后三个键可选，`cfl` **在场**即切换自适应行进（`cfl` 声明但 ≤ 0、`dtMin ≤ 0`、`dtMax < dtMin` 三条拒收）。`cfl` 缺席 = 5.66 的固定步长 march **逐位不变**。

**自适应 march 的簿记**：首步用 case 自己的 `dt`；此后每步 `h = min(dt_next, endTime − t)`，其中 `dt_next` 由 `cflStep` 从**刚提交的流量**读出（`min_k(L_k A_k/|q_k|)` 是最紧的 Courant 预算，静止管不强加约束，**全静止的网络取 dtMax**）；`t` 累加——变步长没有索引可定位一步——但裁剪触发时 `t` **精确吸附**到 endTime（`h = endTime − t` 赋值后判 `clampFires` 置 `t = endTime`），末态不会落在 1.9999999999999998 上。app 打一行 `transient adaptive: cfl … dtMin … dtMax …`，让日志读者知道 t 列记录的是哪一种 march。Python 绑定暴露 `cfl/dt_min/dt_max` 字段与 controls 版 `run_transient` 重载。

**case 与判据**：`tests/cases/transientAdaptiveCfl`（三管 Fast d=0.015/L=1、Slow d=0.04/L=20、Yield τ₀=1 d=0.02/L=20；drive 1.5e5、cfl 0.8、dt 0.004、dtMin 0.08、dtMax 0.1、endTime 2）。case 设计靠两条 identity 联立定 drive：`rule_inf = cfl·ρL²/(drive·τ)`（floor 要绑定就要 rule_inf < dtMin）与 `Re_inf = ρd³·drive/(32 μ²L)`（checker 的层流闭式要复现就要 Re_inf < 2000）——**drive 1.5e5 让两条同时成立**（τ = ρd²/32μ = 0.0598 s，rule_inf ≈ 7.6e-2，Re_inf ≈ 1345 / 1272 全管层流）。第一版 μ=0.01 时 Re 14280 进湍流区，内核走 Churchill 摩阻，checker 层流闭式全面失配 0.83 才暴露——湍流档的 CFL 判据留作边界。

checker `check_transient_adaptive_cfl.py`，**16 判据 0 失败**：

- 声明行命名规则与钳；`steps=` 等于 t 列隐含的 h 数；末态 t 精确 = endTime；
- **调度逐步复现**：`h_i = min(clamp(0.8·min_k(L A/|q_k|(state_{i−1})), 0.08, 0.1), endTime − t_{i−1})` 对全部 26 步成立（首步 = dt 单独断言）——调度是从物理读出来的，不是记住的表；
- 调度访问全部四个 regime（init ×1 / ceiling ×2 / rule ×3 / floor ×20）+ 末步余量钳；
- 全部 27 态被**逐步 h** 的隐式欧拉行复现（Newtonian 闭式、Bingham tw 空间二分；最差 **8.53e-10**）；
- 固定 march 完好：去掉 `cfl` 后 ceil(2/0.004) = 500 均匀步、501 态闭式复现，且与自适应末态**不同**（调度真的动了离散化）；
- 三条拒收（`cfl 0` / `dtMin 0` / `dtMax < dtMin`）。

**本增量的边界**：① 这是 1D 管网的 Courant 规则，**不触碰** §5.38 记录的耦合界面发散墙——那是 3D 块侧的预算，不在 1D 行里；② 显式液量输运（`commitStep` 的 upwind 漂移通量）另有"runs dry − reduce dt"的 fatal，`dtMin` 必须把它挡在门外——本 case 的 dtMin 由 Courant 规则自然覆盖，未单列探针；③ 湍流档的 CFL 判据未单列（case 全层流）；④ 规则只看管流量，不看节点压力变化率——压力驱动的刚性（bulkModulus 大的节点）不进这个规则。

**成本**：内核 ~40 行（`cflStep` + core/app 各一条分支 + loader 校验）+ 1 case 文件 + 1 checker + 1 条 `add_test` + 绑定 3 字段，16 判据；checker 全程 **0.81 s**（Test #85）。全量回归数字见 `docs/gap_analysis.md` §6 实施日志。

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
| `rhoCentralFoam` / 密度基可压缩 Euler 回路（`divScheme` 的 Riemann 通量） | `chflow-3d-compressible`（`src/3d/apps/compressibleApp.cpp`：多面体网格上的广延守恒量递推 + 一阶 Godunov HLLC 面通量 + 理想气体，见 5.16） |
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
| `interFoam` 的 `alphaEqn.H` | `vof::advectAlpha` + `chflow-vof`（速度给定；`surfaceTension{ advect true; }` 时用解出的通量输运 α 并每步重建力；`reconstruction{ transport geometric; }` 时 α 方程整个不再参与，改由平面输运） |
| `MULES::limit` / `limiter` / `localErosion` | 不需要：迎风 + 离散无散 phi 即 M-矩阵，α ∈ [0,1] 是构造性质；几何平流路径上的对应物是 `advectGeometric` 的出/入流双限（比值沿面共享，因此不破坏守恒） |
| `interfaceProperties::correctM()`（`fvc::grad(alpha)`） | `vof::interfaceNormal()`（相对阈值只让界面带持有法向，带外精确 0） |
| `interfaceProperties::k()`（`fvc::div(M)`） | `vof::interfaceCurvature()`（只累加两侧都带界面的面，减去自身法向）；`vof::interfaceCurvaturePlic()`（§5.15：带内两环角邻域截面点的加权二次曲面，κ = div n̂，欠定如实报 0） |
| `surfaceTensionForce`（CSF，`f = −σκ∇α` 的密度归一形式） | `vof::surfaceTensionForce`（`f = −(σ/ρ)κ∇α`，符号从 n̂ 朝内的约定推出；p_in − p_out = σ/(ρR) 内侧更高；`surfaceTension{ kappaPlic true; }` 时 κ 取截面几何） |
| `cutCell-stencil.C` / `PLICReconstruct.C`（`fvModels` 的几何截割） | `vof::cellBoundary` + `vof::clipConvexCell` + `vof::offsetForFraction` + `vof::ConvexRegion`/`clipRegion`（凸单元解析截体积/截面面积/一阶矩、可反复裁剪的凸区域、带括号 Newton 反解平面常数，见 5.12/5.13） |
| `fvex8::interface_normal`（`limitedLevelFour` 重构后取法向） | `vof::reconstruct` 的 `youngs`/`centroid`/`mixed` 与一次精化，加 `fit`：邻域 α 的面积加权最小二乘平面（阻尼 Gauss–Newton，见 5.14）；κ 由截面几何直接给出（`interfaceCurvaturePlic`，见 5.15） |
| `interFoam` 每步用重构面输运 α（几何平流） | `vof::advectGeometric`（面扫棱柱写成半空间列表，逐面一个通量在两侧共享记账；出/入流双限迭代到不动点；`reconstruction{ transport geometric; }`，见 5.13） |
| `interFoam` 把 `f = F/rho` 挂到 `UEqn` 右端 | `PisoControls::cellForce`（每单元加速度进动量右端，由 HbyA/`explicitOperator` 携带，通量经与压力同一套配对看到它） |
| `buoyantBoussinesqSimpleFoam` / `PisoFoam` 的 `TEqn` | `EnergyModel::correctEnergy`（`fvm::ddt + fvm::div(phi,T) + fvm::laplacian(alpha,T,-1)`，读求解器自己的连续性通量） |
| `alphat = nut/Prt`（`turbulentTransport` 的热扩散） | `EnergyModel::alphaCell()` = `kappa/(rho cp) + nut/turbulentPrandtl`（梯度输运类比，前提是网格解析副层） |
| `alphatWallFunction`（温度壁面律，`u_tau y/T+` 面扩散率） | `EnergyModel::faceGamma()` 对 fixedValue wall 面的覆盖：两支 Jayatilleke 律 `T+`（`turbulence/WallFunction.h` 的 `thermalTPlus`），副层极限退回分子导热；Pr < 0.5 拒绝 |
| `fvOptions` 的 `wallHeatFlux` / `fixedGradient` T | `patches{ name{ heatFlux q; } }` → `fixedGradient = q/kappa`（正 = 进入流体）；无键的 patch 缺省 `zeroGradient` = 绝热 |
| `T` 的 `buoyancyForce = beta*(T-Tn)*g`（`BoussinesqBaseForce`） | `EnergyModel::buoyancyForce()` = `−beta (T − tRef) gravity`（只含偏差项，参考静水由 `force` 自己带或被 `pRgh` 开关吸收进 `p_rgh = p − force·x`；偏差力本身不能被吸收，压住它要 buoyantFoam 式浮力修正通量，未做） |
| `decomposePar` + `decomposeParDict`（`hierarchical`/`simple`） | `cf::parallel::Decomposition`（`src/3d/src/parallel/domainDecomposition.cpp`：递归坐标中位数排序 + 连续切片，case.dict 的 `decompose{ domains (…); }`，见 5.17） |
| `processorPolyMesh` / `growPatch` 的 coupled 面与 `/neighbour` patch | `Decomposition::faces(d)` 里的 `cf::parallel::DomainFace{face, cell, peer, sign}` + 惰性 halo 槽（`slot = nCells(d) + idx`）；`cf::parallel::DomainMesh` 的 rank-local 几何 |
| `Pstream::gather/scatter`（`processorLduInterface` 的 send/recv） | `cf::parallel::exchange(dec, local)`（每 peer 一段连续、收发同序的槽位计划 `sendPlan/recvPlan`）；MPI 传输后端尚未接 |
| `Pstream::sum/min/max`（`gsum`/`gMin`） | `cf::parallel::reduceSum/reduceMin/reduceMax`（同 `Allreduce` 语义的进程内归约） |
| `PstreamOps` 作用在 `fvc::grad` / `fvc::div` 的 halo 上 | `cf::parallel::assembledGrad` / `assembledDivPhi`（逐域按全局面序重放串行表达式，因此与串行 `fvc::` **逐位相同**，见 5.17） |
| Fluent **system coupling** / OLGA **embedded detailed section**（管网一条边的内部换成 3D 区，界面上交换压降↔速率） | `cf::coupling::solveCoupled1d3d`（`src/coupling`：界面上的 Gauss–Seidel + 常数松弛，`couple{ edge seeds relax tolerance }`，见 5.18） |
| 3D 区域约化成管网的一个阻力系数（system-couple 的"black box 压降"端） | `cf::coupling::ICoupledElement::lossAt(qMag)` + 生产实现 `DetailBlock`（SIMPLE 块 + `BlockMeasure` 自测：压降、进出通量、逐截面平面通量、探针梯度） |
| 1D 侧承载 3D 区域的边界条件（`fanPressure`/`porousJump` 一类压降模型） | `NetworkPipe::coupledResistance`（正的 `R` 时 `pipeLossAt` 短路成 `R·Q`，与本管摩阻/重位/阀泵互斥，见 4 与 5.18）；瞬态侧 `TransientPipeConfig::coupledResistance` 只替换本管摩阻，重位与区域惯量留在 1D |
| OLGA / Fluent 的**瞬态** embedded detailed section（每个时间层内再把界面上迭代一次） | `cf::coupling::solveTransientCoupled1d3d`（`TransientCoupledControls{ maxSubIterations relax tolerance }` + 逐步 `TransientCoupledStep`，见 5.19） |
| Fluent 的**多** detailed section（一张管网同时挂 N 个 3D 区），以及它们在**时间层内**的版本 | 两个入口都收 `const std::vector<label>& edges` + `std::vector<ICoupledElement*> blocks`：一趟一次网络解、所有伙伴对同一次解求值、系数一起刷新、停止条件取最差边（`couple{ edges ( ... ) starts{ name{ <edge> R; } } }` + `details{ <edge>{ ... } }`，稳态见 5.20、瞬态见 5.25——同一个内核循环，只是后者在每个 `n·dt` 上重谈一次，全部边谈拢才 `commitStep`；每条边的区域惯量来源按**臂**整体二选一（`couple{ inertias ( network partner ) }`，5.24 的单边与 5.26 的多边），或**按边**给（`couple{ inertiaFrom{ arm{ edge network|partner; } } }`——同一趟里一条吃块测网格、一条吃图纸孔径，5.29）；两条边都交给 RANS 伙伴、两个接口算子在同一张管网上同台（逐块 `turbulence{}` + `couple{ accelerates ( 0 1 ) }`，5.27），以及**归属 × 算子的完整 2×2**（四臂只差 `inertias` 与 `accelerates`，5.30），再到**按边表 × 算子的八臂全交叉**（`inertias` 给整臂缺省、`inertiaFrom{}` 逐边覆盖、`accelerates ( 0 … 1 … )` 同臂给出，四张表 × 两个算子，5.31），以及**两条耦合边都在并联位置**的那张拓扑（`Src → A → { X1 ‖ X2 } → B → Snk`，无旁管：两个速率同时谈、`worst` 下标换手、一个节点压降同时对两份各用本边 `M` 的动量行、`τ` 由 `M_X1 + M_X2` 排序，八臂（四张按边表 × 两个算子）整张换拓扑重跑，5.32），以及**伙伴自己带 `du/dt`**（块内 SIMPLE 里装配
隐式 Euler 一步、界面**仍只在摩擦弦上谈**、把块实测的加速度记进 1D 动量行并给出可评分的 `rowError`；
归属多出第三本账 `slab = rho*V_probe/A^2`（`inertias ( … slab … )` / `inertiaFrom{ edge slab; }`）、逐臂
`couple{ ddts ( steady Euler … ) }` 覆盖块字典的 `schemes{ ddt … }`，5.33） |
| Fluent 的 **velocity-inlet 剖面**：system-coupling 的 1D→3D 映射交的不只是流量，还有它横向怎么铺 | `DetailBlock::buildInletProfile`（`starts{ arm{ <edge> R; profile ( flat \| parabolic \| power1over7 ); } }` 逐边逐臂；按面积分归一化到网络要的 `q`，`betaImposed = <ux²>/<ux>²` 由块**自己的面**读回；底层是 `volVectorField` 的逐面 `fixedValue` 覆盖 + `HbyA.adoptFaceValuesFrom(U)`，因此形状进的是预测通量而不只是校正，5.34） |
| 时间层内"已解出但尚未采纳"的状态（分区接口要在同一层重解） | `TransientNetwork::solveStep` / `commitStep` / `pendingPipeFlow` / `pendingPressure`（`setCoupledResistance` 即作废挂起步，见 5.19） |

## 7. 测试体系
- `tests/unit/kernelTests.cpp`（CTest `kernel.unit`）：字典解析、网格几何
  恒等式（散度定理 Σ Sf = 0）、场 BC 线性化、矩阵装配与 BiCGStab 收敛；
- `tests/unit/networkTests.cpp`（CTest `network.unit`）：Churchill 层流极限
  64/Re、单层流管 Hagen-Poiseuille 正/反流、静水平衡（dp=ρg·dz→Q=0）、
  单井 Newton 闭式解、对称 Y 分流、case.dict 解析与求解；
- `tests/unit/transientTests.cpp`（CTest `transient.unit`）：水合物动力学
  精确指数（1e-12）、单步 2×2 隐式 Euler 闭式对照（1e-10）、静止平衡、
  无速率时相不变、带速率时液量指数衰减、slug case.dict 往返与级联递推；
  耦合边另钉五条：`coupledResistance` 管的隐式闭式 `q = Δp/(ρL/A/dt + R)`（1e-10——
  同一根管若还带自己的 2.5e8 摩阻会读到 1.30e-4 而非 1.96e-4，忘掉惯量则是 3.9e-4）、
  耦合管带阀被拒、`solveStep` 不推进任何东西且逐位可重解、`commitStep` 采纳挂起步且
  重复提交被拒、刷新系数即作废挂起步（§5.19）；
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
- `tests/unit/vofTests.cpp`（CTest `vof.unit`）：全部是构造性恒等式而非存储值——
  **一步输运的离散闭式解**（θ = U·dt/dx 的阶跃：锋后不变、锋前第 k 格
  `(θ/(1+θ))^k`，1e-12；α ∈ [0,1] 无裁剪；相体积变化恰等于边界交换
  `U·A·dt·(1−α_outlet)`）、零通量对场与体积都不动（1e-14）、
  **手工径向单位法向**的 `|κ·r − 1| / (dx/r)² < 3`（二阶律，剔除近轴距与触边界
  的模板截断单元）、**平面带回归**（带内带外 max|κ| 恰好为 0.0：曾有的
  ±1/(2dx) 减半假象不得复发）、tanh 平面的法向（32/64 单元成带、
  |n_y|+|n_z| 与 |n_x+1| < 1e-12、`phaseVolume` 与质心位置）、
  **均匀 κ 上表面力成对抵消**（ΣV·f = 3.5e-17 相对——Gauss 梯度在闭合域上的
  恒等式；力恰好非零于带、带外恰好 0；正曲率液滴受力向内，outward 精确为 0；
  κ=0 ⇒ f 逐位为 0）；
- `tests/unit/plicTests.cpp`（CTest `plic.unit`）：几何全部对闭式，不存参照值——
  单位立方体轴向截 0.35 的体积/形心/截面面积/截面形心与 whole/void 标志、
  `offsetForFraction` 的互逆（含分数 1 与 0 落在单元自身面上）、对角截割对
  `{x+y+z ≤ s}` 闭式及其 √3 倍导数（9 个平面，|dV| 4.16e-16、|dA| 9.67e-16）、
  偏移恰过顶点时截面是三角/六边形而体积不跳（曾漏收顶点截面角，差 1/6）、
  `cellBoundary` 在 hex/tet（msh 4.1 与 2.2）/prism 四种网格上闭合到
  `cellVolumes()`/`cellCentres()`（|dV|/V 1.1e-16…1.7e-15、|dc| ≤ 2.0e-15）、
  斜剖两侧互补、`diskRectangleArea` 的六个特例（整圆/四分之一/半圆/弧下条带
  `½(x√(1−x²)+asin x)`/割线下的弓形/完全在外）与 1500² 中点网格一致（最坏 3.5e-5）、
  `exactDiskAlpha` 守恒（Σα 对 πR²/dx² 到 2.1e-10 相对，极限是 `asin(±1)` 端点）、
  圆内单元恰为 1 圆外恰为 0、四种法向取法各自把带 60/60 重构完并读回 α（≤ 4.8e-15）、
  截面面积合计 0.109–0.122 对拉伸周长 2πR·dz = 0.1257、梯度阈值标 124 格而其中
  64 格无界面（PLIC 的带由 α 定，不多不少）、三种**单元内**取法的方向在 40→80→160 上
  不收敛（8.34°→7.94°→7.59°，最坏格 14.11°→14.78°→15.23°），κ 也不收敛
  （−3.696→−3.4997→−3.2509 对 −5，散度模板的限制），因此那三者断言的是**有界**；
  `fit` 断言的是**收敛**——1.16°→0.488°→0.304° 均值、3.10°→1.40°→2.44° 最坏、
  邻域散差 0.0369→0.0196→0.0093 格（§5.14）——另有一组对"真的平面"的 16×16 场断言
  方向 <1e-6°、平面常数差 5e-8 格、邻域散差 6.3e-9 格（同一场上 `mixed` 均值 8.17°）；
  **截面曲率**（§5.15）：同三组圆盘上 `interfaceCurvaturePlic` 带均值 −5.02…−5.12 对 −5
  （最坏格 5.87 = 1.17 倍 1/R），逐组都比模板更近真值；覆盖率 ≥ 0.8·带（`fit`-160 的
  面外退化格如实报 0，`mixed` 全覆盖）；平面锋上 κ 最坏 2.77e-14 且 24/24 全报；
  几何输运三组：`ConvexRegion` 被两刀切出的斜 slab 对两个 cube cap 之差（≤1e-12）、
  box∩box 的扫过棱柱体积恰为 面积×距离（0.4）、切空返回 `empty` 而不是失败；
  C=0.4 的**台阶锋**走 3 步后逐格等于"平移后的锋"（最坏 3.33e-16）且带始终只有 1 格、
  进口交出的体积精确、clamp 恰为 0——这是 α 方程的一阶迎风做不到的；40×40 圆盘自转
  40 步逐步相体积漂移 ≤ 6.9e-16 相对、clamp ≤ 4.42e-17、α ∈ [0,1]、带 60→70 格、
  限流最多 336 面且最小因子 0.5（同格两个出流棱柱在角上重合处）；
- `tests/unit/couplingTests.cpp`（CTest `coupling.unit`）：**分区回路本身**用解析伙伴
  `dp = R0·q + k·q²` 测，三维块只是另一个实现——线性伙伴恰好 2 趟命中、q = 闭式
  0.025、R = 4e6，且网络最终带的就是报告的那个电阻、伙伴最后一次被问的速率就是
  报回来的 q（报告状态一致性）；`relax 0.5` 与 `1.0` 落同一不动点、只是趟数更多；
  二次伙伴收敛到 `dp(q) = Δp` 的解析根 `(-R0+√(R0²+4kΔp))/(2k)`（1e-9，残差 <1e-10）；
  分支管网复现闭式分流 q_X = 0.01 与节点压力 1.8e5/1.6e5/1.4e5/1e5、两支压降相等、
  链闭合；被 `maxOuterIterations` 掐停的回路如实报 `!converged`、恰好跑了 3 趟、
  残差恒等式 `|R·q − dpPa|/dpPa` 仍成立（拿一个精确线性伙伴 + `relax 0.5` 才测得出
  "停了"——`relax 1` 下它一趟就残差 0）；两端定压的两节点网络 `nFree == 0`、
  节点残差恰 0，用作"网络侧不带任何迭代误差"的参照。瞬态回路同一套：12 步线性解析
  伙伴逐步对上离散闭式（首步按惯量分压 + 其后每步 `(1+dt/τ)^{-1}` 衰减，1e-9）、每步
  `block.lastQ() == res.q`（1e-15，一步之内两侧同一状态）、块解总次数恰 `nSteps+1`、
  每步残差 ≤ tolerance、节点压降 == `dpBlock + M_X·Δq/dt`（1e-8）；`relax 0.5` 落同一条
  爬升曲线只是趟数更多；`maxSubIterations 1` 时 fatal 且网络仍停在起点（§5.19）；
- `tests/cases/*`：六十个物理回归算例 + `tests/check_*.py` 断言（1D 两流体
  `twoFluid*` 十七个见 §4.3 的表）：
  - `steadyNetwork`（两井集输层流管网：节点压力/流量/Re/压降对比
    Hagen-Poiseuille 串并联闭式解，含重位水头，1 次 Newton 收敛）
  - `transientBlowdown` / `transientRlc` / `transientSlug` / `transientValve`（见 4.2）
  - `pvtFlash`（四组分天然气 320 K 压力扫描：闪蒸恒等式——相分和、
    `β·y + (1−β)·x = z`、ρ_L > ρ_V、泡点两侧相态塌陷、表节点 = 直接求值）
  - `advectionDiffusion`（Pe=10 一维对流扩散，解析解对比）
  - `pureDiffusion`（Laplace/Poisson，线性场精确再现）
  - `gmshDiffusion`（非结构网格上的同一道题：三棱柱渠道 `tests/meshes/
    prismChannel_22.msh`，20 单元、实测最坏非正交 78°，两端 Dirichlet
    0/1；`T = z` 仍是离散解的精确解（斜面的补偿项恒为 0），断言误差
    1e-10——网格来自 Gmsh 文件而不是 `block{}`）
  - `poiseuilleChannel`（平面渠道 Poiseuille：抛物线剖面、通量守恒、
    `v≈0`、中心线与解析压降对比；全压力 SIMPLE + 自适应 `aP` 收敛 47 次
    外迭代，入口压力对解析 1.2 Pa 符合到 2.2%）
  - `startedChannel`（渠道瞬时启动：10×20×1 网格、体积力驱动、
    `nu dt/dy² = 0.2`、6000 步 / 12 s 物理时间；速度历史对比解析级数
    （早期相对误差 8.7e-3，为 Euler 一阶时间离散；后半程 2.8e-7）、
    终态剖面逐点等于离散抛物线（1.4e-8）、壁面偏移实测等于
    `(a/2nu)(dy/2)²`、校正后散度 1.4e-18）
  - `turbulentChannel` / `turbulentChannelSimple`（充分发展湍流平面渠道，
    1×16×1 单列、体积力驱动、`y+ ≈ 52`；同一个 case 带 `piso{}` 与 `simple{}`
    两个块，两个 app 各读自己的，两个模型各跑一遍）：壁面剪应力对积分平衡
    `u_τ² = f·δ`（PISO k-ε 2.2e-5、SST 3.7e-4；SIMPLE 两者 3.2e-7）、首格 `u+`
    对 Spalding 对数率（说明壁面链在整条求解链上是自洽的）、剖面对称单调、
    `u_tau/y+` 落在壁面函数区间、壁面单元的平衡对与 `nut = κu_τy` 从 vtu 复核；
    定常性判据按耦合分别读取——瞬态看末两次 `uMax` 漂移 <1%，SIMPLE 看
    `momentum`（动量固定点变化）；摩阻对 Pope 光滑平板律
    `1/√f = 2log10(Re_b√f) − 0.6` 的偏差为 PISO k-ε −14.5%、SST +5.4%，
    SIMPLE k-ε −15.5%、SST +6.2%（`Ub/u_τ` 17.6~22.1，期望 20.8）——8 个单元
    跨过 `y+ = 52` 到中心线的网格只能把这个数核对到网格自己的精度，因此带宽取
    20%，而近壁闭合严格复查。checker 把 case 复制到临时目录再改 `model`，剖面
    一律从 vtu 质心取（两种耦合因此比的是同一组数），产物不落进仓库。
    case 的能量块（`turbulentDissipation`，300 K 双壁是唯一热出口）在第 8 组
    断言里核对湍流释热：离散平衡（SIMPLE 1e-15、PISO ≤1.1e-7 W，瞬态残余即
    储能率）、source 对 ∫ε(vtu) 1e-8 相对、壁单元 `ε = u_τ³/(κ y₁) = 475.23`
    （两模型），source/pump 落在 [0.5, 2] 的可见性带（壁函数壁单元非能量
    守恒，实测 0.76~1.18；关掉 ε 源读 0.004/0.03——见 5.11）。
    piso 的 `endTime` 3 → 8（2000 步）让热平衡进入准定常。
  - `vofTranslation`（平面界面被匀速平移携带过域：160 单元、入口 α=1；
    `max|div(phi)|` 恰为 0、α ∈ [0,1]（最小 1.2e-6 就是迎风干尾）、0.5 穿越位置
    与速度对 `x0 + U t`（0.9986·U，落后 0.33 dx）、vtu 与日志两条独立路径互校
    同一穿越、逐行 y 不变（1 ulp）、单调下降、平面曲率 max|κ| = 9.6e-12、
    法向偏 x 轴 6.5e-13；质量平衡断言**带符号的界**，见 5.10）
  - `vofRotation`（液滴绕刚体旋转 0.4 rad：80×80 闭合盒、平滑 4 dx 的圆界面；
    `max|div(phi)|` = 1.1e-19、相体积 −5.1e-14 相对、质心角对 ωt 差 −0.005%、
    半径 0.1995 对 0.2（−0.25%）、带均值 κ 由 −6.00（t=0）到 −6.40（t=0.4）
    对 −1/R = −6.25；初始场的逐单元几何——法向偏径向 1.40° 均值/4.49° 最大、
    |κ·r+1| = 0.137 均值/0.613 最大——从 `vtkInitial` 以 tight 带宽复核，对流后的
    逐单元法向相对**相质心**只断言"仍然圆"（3.98° 均值/19.6° 最大，带由 412 拉宽到
    699 单元），因为迎风会沿流向各向异性地抹开界面）；
  - `vofDroplet`（M4：静止液滴——α 冻结，PISO 读出 CSF 隐含的压力-速度响应；
    80×80 闭合盒全无滑移 + 一块 p=0 锚定、平滑 4 dx 的圆界面、σ=0.08 ρ=1 R=0.16 ⇒
    真值 σ/(ρR)=0.5）：解出的 Laplace 跃变 0.49511（−0.98%，vtu 与日志两条路径
    互校）、远场 −3.1e-6、两侧体相压强平到跃变的 0.63%/0.38%、解场镜像对称
    （|Δp| 0.006%、|Ux+Ux'| 0.17% of max|U|）、离散力净和 ∫|f| 的 2.9e-13、
    闭盒净动量 ∫|U| 的 3.4e-7、解出通量 max|div(phi)| 1.1e-14、α 逐位不变、
    伪速度 max|U| 0.0123 = √(σ/(ρR)) 的 1.7% 且增量递减——见 5.10）；
  - `vofOscillating`（M4：两相耦合——`advect true`，解出的 φ 输运 α、每步重建力；
    80×80 闭合盒、平滑 4 dx 的椭圆 a=0.184 b=0.136、σ=0.08 ρ=1 ν=0.001）：质量守恒到
    1.0e-13 相对（闭盒壁面 α=0）、α∈[0,1] 无截断、质心漂移 ≤1.1e-4、aspect 穿越
    周期 **0.8602 s 对密度匹配 Lamb 值 2π√(ρR³/(3σ)) = 0.8069 s（+6.6%）**（自由液滴
    公式会给 0.571 s——差 35%，穿越钟裁决）、振幅逐周期衰减（0.359→0.180→0.165→0.126，
    双轴都衰减）、dp 呼吸 0.432..0.525 均值 0.502 对 σ/(ρR)=0.5057、释放在
    0.87×毛细速度内、t>0.25 后 ≤0.202、vtu 独立复算体积/aspect/法向/κ——见 5.10）；
  - `vofPlic`（§5.12：同一个旋转盒，界面改为 `fill exact` 的逐单元精确充填，外加
    `reconstruction{ normal mixed }`；40×40、R=0.2、40 步转 0.1 rad）：`max|div(phi)|`
    2.2e-19、相体积对闭式 πR²dz = 0.0125663706144 差 5.9e-15 且全程漂移 ≤6.7e-15、
    带 60 格 = 重构平面数 = 报告数（α 带，不是梯度带）、法向 |n|−1 ≤ 6.4e-13、
    平面度（同格 ∫n·x 常数）≤ 9.2e-13、5 次报告的读回 α 误差 ≤ 2.1e-15、截面积合计
    0.111878 对拉伸周长 0.125664（89%）、切面上的法向偏内径向 8.34° 均值/14.11° 最大、
    `kappaPlic` 带均值 −5.1036 对 −1/R = −5（1.021 倍，§5.15 的截面几何；模板只有 −3.696）、
    梯度阈值标 124 格其中 64 格无界面、
    带由 60 拉宽到 421 而 Σarea 随之到 0.368——那是锋面厚度不是周长；断言的平面几何
    直接从 vtu 的 `nHatPlic/offsetPlic/areaPlic/centrePlic/kappaPlic` 复算；
  - `vofZalesak`（§5.13：Zalesak 转盘——同一个 40×40 盒、同一枚 `fill exact` 圆盘，
    但 `reconstruction{ transport geometric }` 让平面自己带 α，转完整一圈 6.28 rad、
    628 步、Courant 0.195）：每步面通量之争 ≤ 3.3e-14 相对、整圈相体积漂移 9.4e-15、
    clamp ≤ 3.1e-17 且 α **精确**落在 [0,1]、双限最多 35 轮、最多 314/448 面被限
    （最小因子 0.5）、读回 α ≤ 2.55e-15（14 次报告）、带 60→76 格（迎风 40 步就到 421）、
    Σ截面 0.13372 对拉伸周长 0.125664（+6.4%）；形状误差由 checker 自己写闭式
    （圆盘×轴对齐单元的解析交面积；它与求解器的初始化在初始场上互验到 5.0e-13/格，
    所以参照不是求解器算术的复本）：整圈后差圆盘面积的 **3.85%**＝边界法向位移 0.154 dx
    ＝1.92% R，76 格里 2 格的前后翻了一侧，质心偏 0.198 dx（2.5% R）；
  - `vofFit`（§5.14：同一个转盘、同一套几何输运，只把法向换成 `normal fit`）：初始平面
    对解析径向 1.16° 均值/3.10° 最坏（`mixed` 是 8.34°/14.11°）、κPlic −5.1216 对 −5.0000；
    一整圈后带 60→65 格而不是 76、Σ截面 0.128212 而不是 0.133723（周长 0.125664）、
    形状误差 **3.3292%** 圆盘面积而不是 3.8453%、边界位移 0.133 dx 而不是 0.154、
    质心偏 0.141 dx 而不是 0.198、65 格里 0 格翻到锋的另一侧而不是 2 格；守恒、读回、
    α ∈ [0,1] 那些数不变（同一套截割与通量）；日志每次报告给出 `scatter 9.2e-4…1.1e-3`
    与 `rounds 4`，checker 断言散差 < 0.06 格、轮数在预算之内；ctest 上 51 s 对 11 s；
  - `vofDropletPlic`（§5.15：单向静态液滴 + `fill exact` + `reconstruction{ normal mixed }` +
    `surfaceTension{ kappaPlic true; }`——力的 κ 来自截面几何）：κ 带均值 **−6.2918 对
    −1/R = −6.25（0.67%）**、力内力净和 1.2e-12、闭盒净动量 1.2e-6、镜像对称 |Δp| 0.14%、
    伪电流 8.4% 且增量递减、α 逐位不变、max|div φ| 3.9e-14；**断言为实测契约的**是胶式压力
    对 delta 截面力的分布位置——体相均值 dp −30.6%（容差 40%）、α>0.999 体相内仍有 23% 斜坡
    （容差 30%）；闭式 3% 的 Laplace 跃变由抹开场的 `vofDroplet`（−0.98%）持有；
  - `energyConductionFixed` / `energyConductionFlux` / `energyConductionSource`
    （无流动的三道热传导题，逐点从 vtu 质心比对）：两端 fixedValue 的直线剖面误差
    **恰为 0.000e+00 K**、壁面速率 0.2 W 对 `kA dT/dx` = 0.2；`heatFlux` 的符号约定
    （剖面离 Fourier 解 1.6e-8 K、受热面读 319.75 而面上真值 320）；均匀体源的镜像
    对称恰为 0、两端各带走 `S·Vol/2` = 1.0000000000 W、峰值 325 K 对 `SL²/8k`，
    离连续抛物线的 1.56e-2 K 就是中心格式施加壁值的半格偏移。三个 case 的
    `in + source − enthalpyOut` 分别为 −3.975e-12 / −2.287e-10 / −1.112e-10 W；
  - `energyDecay`（PISO，绝热盒里的余弦模态衰减）：`cos(mπ(i+½)/N)` 是中心差分
    Neumann Laplacian 的精确本征向量，checker 重放求解器自己的 `dt` 序列断言
    `1/Π(1+dt·α·λ_d)`（逐报告点 4.64e-13 相对）、200 步、离散均值恰等于 `initialT`
    （1e-12）、`|enthalpyOut| ≤ 1e-12`，再拿这个因子对连续 `exp(−αk²t)`
    （0.373802218 对 0.372707839）——报的是扩散算子的离散误差而不是把它藏起来；
  - `energyBuoyantRest`（PISO，均匀受热闭盒）：`dp/dz` 0.099972 对
    `f_z = β(T−tRef)|g|` = 0.098100 且符号正确、沿程压到直线的残差 3.475e-04
    比 `f·H` = 0.04905、T 不漂移。残差 max|U| = 2.706e-4 是**闭盒耦合自身的地板**
    （无 `energy{}` 块、改用等价 `force` 的对照 case 复现同一个数），断言因此按
    无平衡时的 Stokes 尺度 `f H²/ν` = 2.45 m/s 定标——见 5.11；
  - `energyChannelFlux`（SIMPLE，双壁等热流的渠道 Poiseuille）：23 次扫描收敛、
    `Umax` 0.100000000833 对 `aH²/8ν`、混杯升温 0.0189393944 W 对面供给
    0.0189393945 W（1e-6 相对，不需要任何解析输入）、`Tw−Tb` 2.42409 K 对
    `17qh/35k` = 2.42857（`Nu_Dh` 8.2505 对 140/17 = 8.2353）、两壁偏移差恰为 0、
    形状比 `(Tb−T_centre)/(Tw−Tb)` = 0.285116 对与 q/κ/流量都无关的 39/136——见 5.11；
    - `thermalWallChannel`（SIMPLE，热壁函数渠道，kEpsilon 与 kOmegaSST 各跑一遍）：
    checker 用 Python 重写两支温度律与离散传导问题——报告的 T+ 11.8956036 全同
    （1e-8）、墙速率对 `ρcp·u_tau·A·ΔT/T+`（1e-8）、vtu 场对重述三对角解最大偏差
    3.5e-5 K（kEpsilon）/ 1.7e-5 K（SST，即外迭代的残余而非公式差）、重述系统通量
    恒定 1.5e-12 相对、壁通量是分子首格通量的 3.08 倍、镜像对称 0.0、balance
    ~1e-11 相对——见 5.11；
  - `sodShockTube`（§5.16：密度基可压缩 Euler 的第一道算例——1 m 管段、单格方形截面、
    400 格沿 x，x=0.5 隔膜两侧 `ρ 1 / p 1` 与 `ρ 0.125 / p 0.1` 皆静止，t=0.2 结束；
    全部边界是不可透滑移壁）：checker 用 Python 逐格重放同一套 HLLC 递推，写出的场与复现值差
    **5.0e-13**（契约 1e-9）；闭盒质量与总能量 1e-12、`mom + impulse − mom(0)` 恰 0、锋触壁前
    `impulse = dt(p_R−p_L)A·steps` 闭式 1e-9；星态**没给**求解器，平台逐点对精确采样：
    far left ρ 0.16%/p 0.22%、稀疏扇区（对自己那条自相似剖面）4.58%/6.35%、star left
    0.76%/0.96%、star right 0.02%、far right ~0，接触与激波位置取密度剖面的半跃变交叉点，
    分别偏 0.205σ（1.25 dx）与 0.091σ（0.76 dx），`max u` 0.92961 对 `u*` 0.92745（0.23%）；
    全场 L1 在 400/800 格上给出收敛率 ρ 0.66、u 0.83、p 0.77（断言带 [0.4, 1.5]）。
  - `parallelDecomposition/case.dict` 与 `caseTet.dict`（§5.17：区域分解 + 耦合面 halo +
    全局归约的第一道算例，hex 切 1/2/3/5/7/8/16 域、95 个 gmsh 四面体切 1/3/7/13 域）：
    checker 从 .vtu 几何**自己重跑**同一套坐标中位数排序，逐域核对 owned/halo/internal/
    coupled/boundary/messages 与写出的 `domain` 场，并要求三类面条目恰好重组全局网格；
    分解算子对串行 `fvc::grad`/`fvc::divPhi` 的偏差、halo 校验漂移、`sumDev/minDev/maxDev`
    以及 16 域与 1 域写出的场全部**恰为 0**（不是拟合容差）；线性锚场的梯度由 checker
    独立按散度定理从 .vtu 面几何算出体积后与解析值比（hex 4.8e-14 / corrected 7.7e-13），
    无散通量 `U=(y,z,0)` 在包括被切缝穿过在内的内区格上净出流 2.2e-19（hex）/ 4.2e-17（tet）。
    tet 案例不断言梯度精度：该网格上串行 corrected 格式自己就复现不了线性场（gradAnalytic
    ≈ 217，1 域时即为如此），只断言切割不改变它。
  - `coupled1d3d`（§5.18：分区 1D↔3D 稳态耦合的第一道算例——`Src → A → { X(0.5×0.02×0.02 m
    三维 SIMPLE 块，1000 格) ∥ By(一维旁通) } → B → Snk`，注入 2.4e-4 m³/s、远端定压 1 MPa，
    三个初值电阻 `1.0e6/3.75e6/1.2e7` 各跑一整趟回路）：七组断言全绿——接口身份
    （网络分支流量 == 块 `fluxIn` 到 1e-12 相对、`dpNetwork == R·q`、残差 ==
    `|dpNetwork−dpBlock|/dpBlock`、每趟块自己报收敛）、节点衡算（`Q_A == Q_X+Q_By == Q_B`、
    主导线 carries 注入速率）、压力链（`dp_A+dp_X+dp_B == pSrc−pSnk`、`dp_X == dp_By`、
    `pJ1−pJ2 == dp_X`、汇点守住参考压）、块内守恒（`|fluxOut−fluxIn| ≤ 1e-7·fluxIn`、
    逐横截面平面净通量极差 ≤1e-7·q）、已发展段物理（`dp/dx` 1163 对 `12μU/H²` = 1183
    即 −1.7%、截面均值贴直线散差 0.091 Pa = 压降的 1.5e-4、中心线速 = 1.5 倍均值、
    `入口面积×塞速 == fluxIn`）、从 .vtu 独立重算的剖面（逐格对抛物线最坏 0.0088 m/s =
    峰速 1.5%、横向速度 ≤1.9e-5、逐列通量对 fluxIn 2%）、闭式工作点（端到端 718.2 对 699 Pa
    = +2.7%、分流 1.5777e-4 对 1.6032e-4 = −1.6%、块实测电阻 3.933e6 落在层流 duct 值
    3.75e6 的 1.00–1.06 倍带内、Re_X 20–200 保证解析对照成立）；初值无关性：三种子的
    q/dp/R 互差 1.6e-5 / 1.0e-5 / 4.7e-5 相对。
  - `transientCoupled1d3d`（§5.19：同一个网络、同一块网格改成从静止被突加的注入速率驱动，
    `dt 0.02 × 60 步`）：七组断言全绿——每步的离散 1D 陈述（节点压降 == `R·q + M·Δq/dt`，
    用该步实测的块系数，最差 1.75e-10 相对；旁通支用 Hagen-Poiseuille，1.15e-13；两条主导线
    含首步惯性冲击，2.0e-12）、每步接口身份（`dpNetwork == R·q`、残差 == 两求解器之差、
    `q == q_X`、块自己报收敛，均逐位到 1e-12 以下）、每步三段衡算（`q_A == Q == q_X+q_By == q_B`
    到 1.4e-20）、界面日程（首步 2 趟修正种子、60 步共 95 次块解 ≤ 120 的上限、每步残差
    ≤1e-4）、爬升对离散闭式（首步 1.38651e-4 对惯量分压 1.38623e-4 即 2.1e-4 相对、起步深度
    0.879 而非 1、序列回归 τ 0.26314 对闭式 0.25885 即 1.66%（把区域惯量算两遍读 0.3698）、
    q_∞ 回归对上闭式 1e-6 相对）、末端工况点（距 q_∞ 0.14%、收敛系数 3.93214e6 对 duct 值
    3.75e6 = +4.86%、Re 80 层流、接口末值 == 块场自测 `dpPa` 逐位）、块内守恒与已发展剖面
    （同 §5.18 各带：`dp/dx` 1162 对 1182、剖面最坏 0.0088 m/s、横向 1.8e-5）。
  - `coupled1d3dMulti`（§5.20：一张管网同时挂两块三维详细区，`Src → A → { X1(3D) ∥ By(1D) }
    → X2(3D) → B → Snk`，注入 3.0e-5 m³/s，X1 是 0.5×0.02×0.02 m、X2 是 0.4×0.01×0.03 m
    的平面渠道）：八组断言全绿——每边接口身份（每趟 `残差 == |R·q − dpBlock|/dpBlock` 到 1e-9，
    接受趟 `R·q ==` 网络自己的压降到 1e-9、`q ==` 网络流量 == 块入口积分流量到 5e-16）、
    纯 1D 侧的陈述（两条主导线与串联边都等于 Q 到 ≤1e-8、`q_X1+q_By == Q`、并联两边共享压降
    `dp_X1 == R_By·q_By` 到 1e-13、跨两个界面的压力链 1e-9）、界面日程（首趟确实远：X2 三组
    同为 0.66007、X1 最高 2.0649；趟数 3/2/3 ≤8；块解 16 次 == 2 边 × 8 趟，按日志里
    `SIMPLE iter 1:` 计数）、边互不借力（X2 首趟 `dpBlock` 在三组之间逐位相同，而 X1 首趟速率
    相差 >2×）、三组初值逐边落同一点（跨度 ≤2.4e-5）、闭式工况点（端到端 554.3 对 567.4 Pa
    即 −2.3%、分流 2.01547e-5 对 2.00399e-5 即 +0.57%、收敛系数对各自 duct 平面 Poiseuille
    0.983 / 0.9745）、两块各自的场（守恒 X1 ~1e-9 / X2 5.4e-5、`probeSlope·L` 对 `dpPa`
    0.16% / 0.68%、中心线 1.472 对 1.5、逐列积分流量对块速率 2% 带内、Re 10.3 / 20.4 层流）。
  - `coupledTurbulent`（§5.22：`detail{}` 挂 `turbulence{}` 的三维伙伴，同一个种子电阻跑
    未加速与准 Newton 两臂）：接口恒等式（`dpNetwork == R·q`、报的压降 == 伙伴自己的、残差 ==
    两者之差）、伙伴幂律 `dp ∝ q^1.74` 与 `dp/q` 随速率变 2.12 倍、块内动量闭合 1.1%、`y+` 96.7
    在壁函数带、涡粘 152.8 倍、剖面平坦度 1.094、Pope 光壁律 −4.4%、两臂同工况点 3.7e-4 相对且
    加速臂 5 趟 ≤ 未加速臂 8 趟、块解次数 == 两臂趟数之和 13 == 13。约 42 s。
  - `transientCoupledTurbulent`（§5.23：同一个 RANS 伙伴放进时间层内，27 步 × dt 0.06、
    两个种子、块 320 单元）：每步动量行——详细支 `dp == R·q + M·Δq/dt` 到 7.34e-11、旁路支拿
    **checker 里重写的 Churchill 1977 公开式子**作外部锚闭合到 7.33e-14（不拿网络自己的摩阻反证
    网络）、两条引线含首步惯性冲击 2.47e-12；伙伴幂律由它自己跨 2.3% 速率变化的回应读出
    （块 `q^1.684`/`q^1.687`，旁路 `q^1.729` @ Re 11175，且速率与压降同步单调下降）；爬升对闭式
    ——斜率 τ 0.2460 s 对整条序列回归的 0.2468 s（−0.34%）、弦闭式 0.3990 s（+62%）、逐步 `B`
    0.7835 对实测 0.7842、`m̄` 1.703 与 `τ_chord/τ_slope` 1.622 相差 −4.7%（余下是对数步）；
    首步按惯量分对上闭式 9.5e-4、且离它弛豫到的电阻工况点 2.3%；末端距拟合渐近 3.8e-3%，正是
    τ 自己的 3.3e-3% 指数尾；界面日程（首步 3 趟、全程最多 3 趟 / 预算 6、88 次块解 == 88 趟）；
    两臂同解（速率 2.8e-5、系数 6.2e-5、τ 3e-5、幂指数互差 1.8e-3）；块侧第一性原理独立闭合
    （Re_bulk 19913、y+ 41.6、动量 0.35%、涡粘 86.6 倍、平坦度 1.094、Pope −5.2%、入口段占 2.8%）。
    约 9.1 s。
  - `transientCoupledInertia`（§5.24：同一条 `X ∥ By`、同一个种子，`inertias ( network partner )`
    两臂只差区域惯量的**来源**）：归属可核对——两臂各 60 步的 `M` 是同一个常数，孔径臂 == `ρL/A`
    = 1273239.5447351625、区域臂 == 块自测 `ρV/A²` = 1000000.0000000163（各 1e-12），比值 ==
    `A_square/A_circle` = 1.273240（1e-12，即 +27.3%）；每臂动量行拿**自己那臂的 M** 闭合到
    1.75e-10 / 1.42e-10；首步纯惯量分压各按闭式到 1e-3（1.386514e-4 对 1.515730e-4，差 8.2% 的
    工况点，过报惯量的臂只得 0.9147）；τ 拟合 0.263142 / 0.239584 s 对闭式 0.25885 / 0.235038 s
    （+1.66% / +1.93%），τ 之比 1.09833 对惯量之比稀释后的预测 1.10129（−0.27%）；分工谈拢后惯性差
    消失（1.292e-5 → 1.693e-7 m³/s），两臂末端速率差 0.107%、外推工况点差 0.0014%；120 个时间层
    最差界面残差 9.851e-5 < 1e-4；块侧 `dp/dx` −1161.58 / −1162.83 对面 Poiseuille 1182 / 1183、
    入口面 4.0e-4 m² 正是惯量所参考的截面。代价 2137 次块扫 / 175 趟，约 53 s。
  - `transientCoupledMulti1d3d`（§5.25：稳态多边算例的拓扑与两块网格，换成时间步进驱动，
    `couple{ edges ( X1 X2 ) starts{ lo / mix } inertias 缺省 network }`，12 步 × dt 0.05）：
    十族断言全绿（76 条）——逐边逐步身份（`dpNetwork == R·q`、残差 == 两解之差、`q ==` 网络流量、
    `newton < 1e-11`、`t == n·dt`）、每边动量行用自己装上的 `M`（X1 3.04e-11 / X2 1.82e-13，
    旁路与引线 1.3e-12…1.4e-11）、串联边每步速率 == Q 且首步压降 1554.2153 Pa 里 1086.5 Pa 是
    启动惯量、**界面日程**（`worst == argmax`、每步恰 1 趟被采纳、`lo` 9/21 与 `mix` 8/20 趟超差
    被重跑、其中 15 趟只有一条边还超差）、块解次数 == 趟数 × 2（42 与 40）、边互不借力（X2 首趟
    残差两组同为 0.6601、第 1 步 `R` 差 3.5e-6）、两组初值逐边落同一点（跨度 1.02e-7）、
    首步纯惯量分工 1.76228e-05 对闭式 4.6e-12、爬升对上隐式欧拉自己的 `1/(1+dt/τ)` = 0.841045
    （实测 0.840671，−0.044%；连续闭式 τ 0.264555 s 与格式等效 0.288094 s 差 +8.90%）、去向 ==
    电阻分流 2.01566e-05（拟合 2.01554e-05，逐步 0.008% 内）、收敛系数对 duct 0.9826 / 0.9744
    （与 §5.20 的 0.983 / 0.9745 互锚）、两块各自的场（`dp/dx` 145.6 对 148.3、1177 对 1200，
    中心线 1.472，入口面 == 自测 `ρV/A²` 所参考的截面，Re 7.9 / 12.0 层流）。代价 82 次块解 /
    41 趟 / 24 个谈拢的时间层，235 s。
  - `transientCoupledMultiInertia`（§5.26：同一张管网、同样两条边，只把**两臂的区域惯量来源**换成
    `couple{ inertias ( network partner ) }`，种子逐位相同，12 步 × dt 0.05）：九族断言全绿（180 条）——
    每边每臂装上的 `M` == 该来源（1e-12）而块自测 `inertance` 两臂逐位相同且 == `ρV/A²`（1e-9）、
    两条边被名义孔径错估的幅度不同（X1 +27.3% 比值 1.273240、X2 +69.8% 比值 1.697653，都 ==
    实截面/内切圆到 1e-12，差 0.4244）、界面系数两臂每步同一个数（X1 最差 7.40e-4 / X2 1.27e-6）且
    `dpNetwork == R·q` 到 **0.0 Pa**、串联边速率两臂逐位相同而首步节点压降差 446.497744 Pa ==
    闭式 `ΔM·Q/dt` 446.497745 Pa（启动账 1086.4977 对 640.0000 Pa）、并联边首步按各自 `M` 分工
    （1.76228e-05 对 1.90996e-05，4.6e-12 / 5.9e-13）、τ 0.264555 对 0.240207 s（比值 0.907968
    对闭式 0.908028）、衰减各自对上 `1/(1+dt/τ)`（0.840671 / 0.827248，−0.044% / −0.056%）、
    目的地与归属无关（拟合 2.015535554e-05 对 2.015510611e-05，0.0012%）而 12 步末端仍差 +1.23%
    （方向由惯量定）、末端系数对平面 Poiseuille 0.9826/0.9828 与 0.9744/0.9744 == §5.25 的数、
    块解 == 趟数 × 2（40 与 34）、每步恰 1 趟被采纳且 `maxSubIterations 8` 未触及（20 对 17 趟、
    24 个谈拢的时间层，最差被采纳残差 8.5303e-5 / 9.5498e-5）。代价 74 次块解 / 37 趟，约 365 s。
  - `transientCoupledMultiTurbulent1d3d`（§5.27：§5.25 的两条边同时换成 §5.23 那类 RANS 伙伴，
    `couple{ edges ( X1 X2 ) accelerates ( 0 1 ) starts{ gs / iqn } }`、两块 `turbulence{ kOmegaSST }`，
    代理管孔径取区域的等效截面 ⇒ 归属被中和成同一个数，12 步 × dt 0.05）：十族断言全绿（141 条）——
    两臂各自分段且只差 `accelerate` 一个开关、`dpNetwork == R·q` 两臂 24 层最差 **0.0 Pa**、
    装的 `M` == `ρL/A` == 块自测 `ρV/A²`（2.2e-16 / 4.6e-15）、四根 1D 管各自闭合（对发表的
    Churchill 弦最差 7.6e-14）、串联边首步节点压降 11638.381828 Pa 对 `R·Q + M·Q/dt` 到 1.9e-15
    （其中 9000.0 Pa 是启动惯量、2638.4 Pa 才是摩阻）、X1 的律 `dp ∝ q^1.682`（层流伙伴 1.06）
    而 X2 读不出指数（速率锁在注入值 1.1e-16、压降跨度 3.4e-8）、末端系数对平面 Poiseuille
    5.7664× / 7.8174×（层流同几何 0.9826 / 0.9744）、**斜坡的恢复项是斜率和**（τ 预测 0.241385
    对拟合 0.240828 s，+0.23%；弦式 0.394059 s 读高 64%；弦/斜率 1.632 对 `m̄` 1.7012）、离散闭式
    管住整条轨迹（衰减 −0.048%、首步 −0.075% 且离电阻分 2.6%、目的地 +0.179%、末端 0.2475%
    对尾巴 0.2477%、轨迹 0.159%）、`worst` 下标在两臂全部 53 趟里都是 0 且 29 趟超差中 25 趟是
    串联边陪跑、块解 == 趟数 × 2（54 / 52）、预算 6 从未尽用（4 / 3）、两臂末端同一状态（≤5e-6）
    与同一 τ（≤5e-3）、每块场独立闭合（`y+` 41.7 / 58.4、`ν_t/ν` 86.8 / 96.8、`Re` 19963 / 30000、
    `|dp/dx|·H/2 == ρu_τ²` +0.37% / +6.5%、Pope −5.1% / −5.2%）。代价 106 次块解 / 53 趟 /
    24 个谈拢的时间层，约 90 s。
  - `transientCoupledMultiTurbulentInertia`（§5.28：§5.27 的两块 RANS 伙伴 × §5.26 的归属，且把
    上个算例中和掉的那格自由度**放歪**——两根代理管都写成名义圆 `X1 0.005 / X2 0.0049`，
    `couple{ inertias ( network partner ) starts{ bore / region } }`、两臂种子逐位相同且都不加速，
    12 步 × dt 0.05）：十族断言全绿（135 条）——每边每臂装上的 `M` == 该来源（1e-9）而块自测值两臂
    逐位相同且 == `ρL/(H·W)`（8000000 / 3750000）、两条边被错估的幅度**不同**（1.273240 对 1.696946，
    ⇒ 没有全局修正因子能同时复现两臂）、`dpNetwork == R·q` 两臂 24 层最差 **0.0 Pa**、每边行用**本臂**
    那份 `M` 闭合（并联 4.0e-11 / 8.8e-11、串联 2.8e-14 / 3.1e-14；对发表的 Churchill 弦 7.8e-14 /
    7.2e-14、lead 2.0e-16）、τ 是 `(M_X1+M_By)/斜率和`：孔径臂 0.277231 对预测 0.276586 s（−0.232%）、
    区域臂 0.240829 对 0.241385 s（+0.231%），弦式两臂都读高 63.5% / 63.6%（`c/s` 1.6386 对 `m̄` 1.6988、
    1.6325 对 1.7012），离散闭式各自管住整条轨迹（0.0086% / 0.0023%）且首步按本臂 `M` 分工
    （4.52323e-05 对闭式 4.5021e-05，离电阻分 9.9%；5.09933e-05 对 5.10603e-05，离 2.6%）、
    **τ 之比 == 惯量之比**（1.15115 对 1.15649，−0.46%）、目的地不含 `M`（渐近速率 4.9789e-05 对
    4.97833e-05，0.0113%）而 12 步末端仍差 1.492% 且写重那臂落得更远（−1.2592% 对 +0.2475%）、
    伙伴指数两臂读出同一条律（1.6762 / 1.6820 ⇒ 归属动 τ 不动 `m`）、串联边两臂逐位相同
    （速率 1.4e-20、系数 8.8e-11、连首步界面残差都是 2.3003e-07）而首步节点压降 17910.892818 对
    11638.381827 Pa，两臂之差 6272.510991 Pa == 闭式 `ΔM·Q/dt` 且各自闭合到 3.8e-15、`worst` 下标在
    两臂全部 57 趟里都是 0（30 趟里 16 趟、27 趟里 13 趟是已谈拢的串联边陪跑）、块解 == 趟数 × 2
    （60 / 54）、预算 6 从未尽用（4 / 4）、末端系数对平面 Poiseuille 5.7064× / 5.7664× 与 7.8174×、
    每块场独立闭合（`y+` 41.2 / 41.7 与 58.4、`ν_t/ν` 85.7 / 86.8 与 96.8、`Re` 19665 / 19963 与 30000、
    `|dp/dx|·H/2 == ρu_τ²` +0.32% / +0.36% 与 +6.48%）。代价 114 次块解 / 57 趟 / 24 个谈拢的时间层，
    约 45 s（串行）。
  - `transientCoupledMultiMixedInertia`（§5.29：§5.28 的两块 RANS 伙伴与两个写歪的名义孔径，但归属
    从**臂**交到**边**——`couple{ inertias ( network partner network ); inertiaFrom{ mixed{ X1 partner;
    X2 network; } } }`，三臂 `bore / region / mixed` 只差那张表，种子逐位相同、都不加速，12 步 × dt
    0.05）：十族断言全绿（132 条）——mixed 那一行**同时装两个来源**（X1 8000000 吃块测、X2
    6363546.246 吃图纸，1e-9）而块自测值三臂逐位相同（8000000 / 3750000 == `ρL/(H·W)`）、**叠加不是
    平均**（mixed 的 X1 全部 12 步与 region 臂逐位相同：速率/系数/块报压降/残差；它的 X2 不是；mixed
    的 X1 首步 5.09933e-05 对 bore 臂 4.52323e-05 m³/s，差注入速率 4.80%）、`dpNetwork == R·q` 三臂
    36 个边×步最差 **0.0 Pa**、每边行用**本边**那份 `M` 闭合（并联 4.0e-11 / 8.8e-11 / 8.8e-11，串联
    2.8e-14 / 3.1e-14 / 3.1e-14；对发表的 Churchill 弦 7.8e-14 / 7.2e-14 / 7.2e-14、lead 2.0e-16）、
    启动的账按边走（串联边首步节点压降 bore 17910.892818 对 region 11638.381827 Pa，mixed **== bore**
    那份，与 region 之差 6272.510991 Pa == `ΔM·Q/dt`，各自闭合到 2.2e-15）、τ 跟着会爬那条边（bore
    0.277231 / region 0.240829 / mixed **0.240829 s 与 region 逐位相同**），预测 `(M_X1+M_By)/斜率和`
    0.276586（−0.232%）对 0.241385（+0.231%）、τ 之比 1.15115 对惯量之比 1.15649、伙伴指数三臂同一条
    律（1.6762 / 1.6820 / 1.6820 ⇒ 归属不动 `m`）、目的地不含 `M`（渐近速率 4.9789e-05 / 4.97833e-05 /
    4.97833e-05，差 0.0113%）而末端速率仍差 0.620%、日程 30 / 27 / **27**（与 region 同趟数，因为它谈
    的是同一个 X1）且超差 18 / 15 / 15、`worst` 下标三臂 84 趟全 0、块解 == 趟数 × 2（168）、预算 6 从未
    尽用（4 / 4 / 4）、最差被采纳残差 1.617e-4 / 1.369e-4 / 1.369e-4、每块场独立闭合（`y+` 41.2 / 41.7
    与 58.4、`Re` 19665 / 19963 与 30000、`|dp/dx|·H/2 == ρu_τ²` +0.32% / +0.36% 与 +6.48%）。代价
    168 次块解 / 84 趟界面 / 36 个谈拢的时间层，串行约 147 s（ctest `-j 8` 下 194 s）。
  - `transientCoupledMultiTurbulentInertiaNewton`（§5.30：把归属与接口算子**交叉**成完整 2×2——同一
    张管网、同样两块 k-ω SST 伙伴、同样两个写歪的名义孔径，四臂只差
    `inertias ( network network partner partner )` 与 `accelerates ( 0 1 0 1 )`，种子逐位相同，12 步 ×
    dt 0.05）：十族断言全绿（172 条）——同来源两臂装上**逐位相同**的 `M` 表（bore 对 10185916.36 /
    6363546.246，mesh 对 8000000 / 3750000）而块自测值四臂只有两个数、`dpNetwork == R·q` 四臂 96 个
    边×层最差 **0.0 Pa**、每边行按**本边**那份 `M` 闭合（并联 4.0e-11 / 4.0e-11 / 8.8e-11 / 8.8e-11，
    串联 2.8e-14 / 3.2e-14 / 3.1e-14 / 3.0e-14）、算子轴只买到日程（同对两臂末端速率差 3.2e-07 /
    2.0e-06 相对、渐近差 2.1e-06 / 4.0e-06、τ 差 8.9e-05 / 8.5e-04 相对、串联首步踢出差 0.044 / 0.053 Pa，
    加速臂各少 1 趟）、归属轴照旧（τ 预测 −0.232% / −0.247% / +0.231% / +0.158%，τ 之比 1.15115 对惯量
    之比 1.15649 在两对里都成立，目的地四臂同一渐近速率 spread 0.0119%，伙伴指数 1.6762 / 1.6763 /
    1.6820 / 1.6816 极差 0.0058）、准 Newton 那处不整齐被量化（被拓扑钉死那条边的残差地板 2.300e-07 对
    1.670e-05 / 2.001e-05，即 73× / 87×，但都仍在 2e-4 之内，真正被谈那条边 1.617e-04 对 1.617e-04 与
    1.369e-04 对 1.372e-04）、`worst` 下标四臂 112 趟全 0、预算 6 从未尽用（4 / 3 / 4 / 3）、块解
    224 == 趟数 × 2、每块场独立闭合（X1 的 `y+` 随本边来源落在 41.2（孔径臂）/ 41.7（区域臂），X2 恒
    58.4；`Re` 19665 / 19963 / 30000、`|dp/dx|·H/2 == ρu_τ²` +0.32% / +0.36% / +6.48%）。代价 224 次块解 /
    112 趟界面 / 48 个谈拢的时间层，app 四臂串行 93 s、checker 串行约 99 s。
  - `transientCoupledMultiThreeAxis`（§5.31：三条轴同台——`inertias` 给整臂缺省、`inertiaFrom{}` 逐边
    覆盖、`accelerates ( 0 0 0 0 1 1 1 1 )` 给算子，四张表（两边图纸 / 两边网格 / X1 网格 × X2 图纸 /
    镜像）× 两个算子 = **八臂**，同一张管网、同样两块 k-ω SST 伙伴、八组逐位相同种子，12 步 × dt 0.05）：
    十族 **373** 条断言全绿——同一来源在两种算子下装**同一份** `M`（`network` 恒 10185916.36 /
    6363546.246，`partner` 恒 8000000 / 3750000 kg/m⁴）、八臂 192 个边×层 `dpNetwork == R·q` 最差
    **0.0 Pa**、每边行按**本边**那份闭合（并联 ≤8.8e-11、串联 ≤3.2e-14）；**可分性按"这条边真正会动的
    量"量出来**：X1 的速率历史在同来源臂之间准 Newton 四对全部**逐位相同**、平谈 `partner` 源逐位相同而
    `network` 源差 **2.9e-10**，异来源 **4.801e-02** ⇒ 分离度 **1.6e8×..4.8e14×**；被钉死的 X2 用节点压降
    量，同来源平谈 1.1e-09 / 1.7e-09、加速 4.9e-07 / 7.5e-07，异来源 3.502e-01；两镜像表在**两条边**上都
    分开（4.801e-02 与 3.502e-01）；首步踢出按边走（平谈 17910.892818 对 11638.381827、差 6272.510991 ==
    `ΔM·Q/dt`；加速臂付同一份到 0.009 Pa，混合臂之差 6272.502249 对 `ΔM·Q/dt` 1.4e-6 相对）；τ 恒等于
    "做斜坡那条边"的来源（混合臂与同源整臂差 0.0 / 4.0e-09 相对，而来源轴本身值 0.151）；目的地八臂同一
    速率（渐近 spread **0.0119%**）；算子只买到日程（同表两臂速率差 9.7e-07 / 4.4e-06，趟数 30/27/27/30 对
    29/26/26/29）；钉死那条边的残差地板 **73× / 87× / 87× / 72×** 而被谈那条边 1.617e-04 / 1.369e-04 不动；
    `worst` 下标 224 趟全 0、块解 448 == 趟数 × 2；X1 的终点按本边来源分两点（4.9162e-05 m³/s、`y+` 41.2、
    `Re` 19665、闭合 +0.32% 对 4.99065e-05、41.7、19963、+0.36%）。代价 448 次块解 / 224 趟界面 / 96 个
    谈拢的时间层，app 八臂串行 275 s（ctest 单条 262.7 s）、checker 串行 267 s。
  - `transientCoupledMultiBothParallel`（§5.32：把 §5.31 的八臂按边表整张搬到**两条耦合边都在并联位置**的
    拓扑上——`Src → A → { X1 ‖ X2 } → B → Snk`、无旁管，两个速率都要谈；只改拓扑与两个名义孔径，
    四张表 × 两个算子、种子 / 容差 2e-4 / 6 趟预算照旧，12 步 × dt 0.05）：十族 **332** 条断言全绿——两边都不
    被钉死（离注入最远 0.397 ~ 0.706，跨层行程 0.0178 ~ 0.1026，两支之和逐层等于注入）、`worst` 下标
    279 趟逐趟等于两残差 argmax 且边 1 当最差 **4 ~ 13** 次（§5.31 全 0）；**同一个 `J1 − J2` 同时对
    两份动量行、各用本边那份 `M`** 闭合到 **1.1e-10 / 7.6e-11**，弦 192 个边×层仍 **0.0 Pa**；每个时间层
    的分流可由外面重算（`q_e=(dp+M_e q_prev/dt)/(R_e+M_e/dt)`、`Σq_e=Q_IN`，最差 **4.7e-12**），换错账
    立刻偏 **0.047 ~ 0.108**；**`τ` 阶梯**：四张表分子 1.65495e7 / 1.43635e7 / 1.39359e7 / 1.175e7
    kg/m⁴ 对 `τ` 0.29115 / 0.25287 / 0.25075 / 0.21374 s，两算子内严格单调、**`τ` 比 == 分子比到 3.40%**，
    两条混合表**严格落在两张整表之间**；单一指数穿 12 层偏离 ≤**0.258%** 行程；目的地八臂收敛到
    0.409319 ~ 0.409364（spread **0.0206%**）且等于本臂两边幂律（X1 1.6601 ~ 1.6861、X2 1.6703 ~ 1.6842）
    的电阻分流 ≤0.070%；**没有乘客**：八对地板比 0.044 ~ 25.6 且方向不一致（5/8 平谈更差），最差地板
    ≤ 容差 21.1%，对照 §5.31 钉死那条边的单向 73× ~ 87×；趟数 32/38/34/37/32/37/32/37 = **279**、块解
    **558** == 趟数 × 2、96 层谈拢；X1 `y+` 40.1 ~ 41.4 / 闭合 +0.20% ~ +0.33%，X2 `y+` 37.3 ~ 38.1 /
    +4.35% ~ +4.38%。代价 558 次块解 / 279 趟界面 / 96 个谈拢的时间层，含跑 app 的 ctest 单条 377.5 s。
- `tests/python/test_bindings.py`（CTest `python.bindings`）：Python 绑定
  构建/求解/读结果与闭式解对比、异常翻译、case.dict 往返；瞬态侧校验
  罐排气 RC 包络、逐步 `TransientNetwork.step` 与批量 `run_transient`
  结果一致、水合物动力学指数；pvt 侧校验临界点 Z、纯组分泡点 ≡ 饱和压、
  `tp_flash` 的物料衡算与逐组分逸度相等（由 `ln_phi` 复算）、case.dict
  加载 + 泡点两侧相态、`liquid_density_table` 与直接求值一致；
- `cli`（CTest `cli.unit` + `cli.runSteady` + `cli.runTransient` +
  `cli.runPvt`）：cargo 参数解析单测 + 通过 CLI 端到端跑稳态/瞬态管网与
  pvt 闪蒸冒烟测试（存在 cargo 时才注册）；
- `tests/check_coupled_rheology.py`（CTest `verify.coupledRheology`，§5.63）：
  一份 case 派生 5 个稳态臂（未声明 / 声明但牛顿 / `n=0.5` / `n=1.5` / HB 屈服）
  与 1 个瞬态臂，另加 4 个拒收 case，**306** 条判据——声明但牛顿档与未声明档的
  求解日志逐位相同、从 `nuWall/gdotWall/tauWall` 独立重建 `K γ̇^(n−1)+τ_0/γ̇`、
  `|probeSlope|·H/2 == τ_w`、`dp` 对同一本构的平面渠道闭式（含塞流积分）、
  逐趟界面残差 `|R·q − dpBlock|/dpBlock`、1D `Re` 仍逐位是 `ρvd/μ`、
  瞬态每步动量行与末态落点、四条拒收边界；
- `tests/check_steady_all_below_threshold.py`（CTest `verify.steadyAllBelowThreshold`，§5.64 / B6）：
  2 份入库 dict + 4 个临时派生臂（移位定压 / 阈值相等 / 单管放开 / 自由节点注入），
  **45** 条判据、0 内核改动——从 dict 重算 `4τ₀L/d` 与被问的 dp、全死档
  `iterations=0 residual=0`、E 随定压均值移动、严格阈值不等式、0 迭代不等于网络死了
  （单管仍按 Buckingham 闭式流动）、注入臂 19 步 Newton 后 E 与 dp 落在闭式反解上；
- `tests/check_transient_yield_index_valve_k.py`（CTest `verify.transientYieldIndexValveK`，§5.65 / C2）：
  16 管（持停三管 / 四阀家族 / n=0.5 对 n=1.5 / 重力对 / 三停止管）+ 稳态孪生，
  **56** 条判据、0 内核改动——全部 2001 状态 × 11 根管被「HB 剖面积分 +
  K ρ|Q|Q/2A²」的隐式 Euler 递推复现（最差 4.347e-10）、阀分成 133× 的指数效应、
  三停止迹逐状态复现、稳态同根会合（打印分裂经 §5.47 反演读回速率）；
- `tests/check_transient_time_decoupled.py`（CTest `verify.transientTimeDecoupled`，§5.66 / C4）：
  6 臂 + 拒收，**40** 条判据——全部报告态被**逐步 h** 的隐式 Euler 行复现
  （Newtonian 闭式、Bingham 在 tw 空间二分，最差 3.27e-11）、末步裁剪落在
  endTime 精确、reportEvery 3/5 仍落末态且与 report 1 逐位相同、近整除
  0.12/0.03 不误计、endTime 整倍数向后逐位兼容、endTime 0 拒收；
- `tests/check_transient_coupled_multi_turbulent_triple_table.py`
  （CTest `verify.transientCoupledMultiTurbulentTripleTable`，§5.67 / B5）：
  5 臂，**21** 条判据、0 内核改动——三张逐边表（`ddtFrom`/`inertiaFrom`/
  `convectiveToRow`）+ 逐边 profile 首次同台；quasi+slab+charge 三表叠一边的
  新格可行（row 用 slab 的 M 付 M·dq/dt 且 `rowBeta == dpConv`）；base 臂
  1e-12 复现 §5.41 p_net 的 level-1/8 状态与 pass 序列；80 状态 row 闭合
  slack 恰 0、charge 只动一行的一项（邻居逐位相同）、chord 在容差内；
- `tests/check_transient_adaptive_cfl.py`（CTest `verify.transientAdaptiveCfl`，
  §5.68 / D5）：**16** 条判据——调度逐步被上一态流量的 Courant 规则复现
  （`dt_next = clamp(cfl·min_k(L A/|q_k|), dtMin, dtMax)`）、全部 27 态被
  逐步 h 的隐式 Euler 行复现（最差 8.53e-10，全管层流）、固定 march
  （`cfl` 缺席）500 步闭式复现且与自适应末态不同、四个 regime 全访问、
  `cfl 0` / `dtMin 0` / `dtMax < dtMin` 三条拒收；
- 运行：`cmake -S . -B build && cmake --build build && ctest --test-dir build`
