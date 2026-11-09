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
预测散度收敛到 HbyA/压力通量的固有失衡量（验证渠道 ~5e-3），
`continuityTol` 只在失衡较小的算例可达；**`|dp|` 分支（`pcorrTol`）才是
该结构下的真实外收敛信号**。两分支都过动量残差门 `momentumTol`
（`max|U−Uold|/max|U|`，含动量欠松弛与压力响应的总固定点变化）——
只看连续性判不了动量平衡：流向均匀的单列渠道所有 `phi` 恒 0，
第 1 次外迭代就读到 `continuity 0`，壁面剪切还没进入平衡（实测关掉
`momentumTol` 时 1 次迭代即"收敛"，`u_τ` 只有 0.554，对积分平衡差 11%）。
加门后体积力渠道 ~30 次外迭代收敛（k-ε 27 次、SST 39 次），`u_τ` 与
拖力平衡符合到 3e-7，与 PISO 750 步的定常解一致（Ub 10.96 vs 11.10）。

刻在实现里的经验结论（回归历史，勿回退）：

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
伙伴今天可以是湍流块（§5.22），但过界面的仍然只有那一个标量。

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
④ 湍流 × 归属只测了"一并联 + 一串联"，且两臂都不加速；湍流 × 归属 × 准 Newton 三者同台、以及
两边都在并联位置（`worst` 会换手）未测。⑤ `tolerance` 与日程仍全网一个数：30 趟里 16 趟
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
④ 界面仍是**每边一个标量电阻**（§5.21 边界②）。⑤ 湍流 × 归属 × **准 Newton** 三者同台仍未测。
⑥ app 头部的 `region inertia <词>` 行仍只描述第一臂，per-edge 时要靠逐边日志行读——头部那一行
如何同时如实描述三臂仍未解决（§5.26 边界里就记过）。

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
| Fluent 的**多** detailed section（一张管网同时挂 N 个 3D 区），以及它们在**时间层内**的版本 | 两个入口都收 `const std::vector<label>& edges` + `std::vector<ICoupledElement*> blocks`：一趟一次网络解、所有伙伴对同一次解求值、系数一起刷新、停止条件取最差边（`couple{ edges ( ... ) starts{ name{ <edge> R; } } }` + `details{ <edge>{ ... } }`，稳态见 5.20、瞬态见 5.25——同一个内核循环，只是后者在每个 `n·dt` 上重谈一次，全部边谈拢才 `commitStep`；每条边的区域惯量来源按**臂**整体二选一（`couple{ inertias ( network partner ) }`，5.24 的单边与 5.26 的多边），或**按边**给（`couple{ inertiaFrom{ arm{ edge network|partner; } } }`——同一趟里一条吃块测网格、一条吃图纸孔径，5.29）；两条边都交给 RANS 伙伴、两个接口算子在同一张管网上同台（逐块 `turbulence{}` + `couple{ accelerates ( 0 1 ) }`，5.27） |
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
