# 贡献指南（CONTRIBUTING）

## 1. 开发环境与构建

- C++20（Apple clang 14+ / GCC 12+ / Clang 14+ 实测通过），CMake ≥ 3.24；
- 内核层（`src/common`、`src/1d`、`src/3d`、`src/coupling`）**零第三方依赖**：只允许 C++ 标准库。
  Eigen/PETSc/spdlog 等只能出现在未来的可选后端模块，不进入内核；
- 测试需要 Python 3.10+（只用标准库）；
- 可选层：Python 绑定 `src/bindings` 经 FetchContent 拉取 pybind11
  （首次 configure 需要网络；离线或缺 Python 开发头文件时用
  `-DCHFLOW_WITH_PYTHON=OFF` 跳过，C++ 构建不受影响）；CLI `cli/` 为
  零依赖 Rust，未安装 cargo 时对应 CTest 用例自动不注册。

```bash
cmake -S . -B build
cmake --build build -j8
ctest --test-dir build --output-on-failure -j8  # 全量回归（117 项）约 260 ~ 1320 s ≈ 4.3 ~ 22.0 分钟（本轮实测 1311.16 s，116/116 绿，
                                             # 注：本轮与两个并发 ctest 共享 CPU，clean 基线约 404 s），
                                             # 墙钟的大头是最慢那几条（同批负载会把它们彼此推高；本轮各用例时间之和
                                             # 7371.58 s，五条最重的耦合用例合计 3568.0 s，占 48.4%）：
                                             # 形状账 x 归属账 875.73 s、稳态 profile 级映射 799.96 s、
                                             # 多边 x 惯量归属 699.54 s、逐边 ddt 635.59 s、对流动量账收进行 557.16 s、
                                             # VOF 振荡 233.76 s；非牛顿那几格里跨过这个量级的只有瞬态静止支那一条：
                                             # 稳态 HB 0.36 s、瞬态幂律 8.77 s、瞬态屈服 x 阀 12.95 s、
                                             # 稳态转捩 locus 0.24 s（本轮最便宜的 verify 用例之一：locus 是纯读出，一次 app 调用、
                                             # 零步进、零重解网格）、瞬态阀日程 x 静止支 17.38 s（12 x 2001 全轨迹重解，
                                             # 两次 app 调用可忽略）、瞬态浮动节点 x 日程 x 静止支 11.95 s（9 x 2001 加两对
                                             # 30 s 长跑全轨迹重解，三次 app 调用合计 0.25 s）、瞬态分支节点 x 双向闩锁
                                             # 20.94 s（7 臂 x 2001 加一段 18 s 长跑，三次 app 调用合计 0.76 s）、
                                             # 瞬态分流节点 x 双向闩锁 109.54 s、瞬态死带初值扫描 111.28 s（12 臂 x 2001
                                             # 加同一批 12 臂 x 20001 的长跑，四次 app 调用合计 3.43 s——贵的还是 checker：
                                             # 48024 个速率的壁应力空间二分复现 + 20001 状态上的 1/gap 直线）、
                                             # HB 停泊带 72.18 s、HB 停泊带 x 阀日程 13.08 s（三副各一次 app 调用约 0.9 s，
                                             # 贵的是逐状态的截面积分复现，但只有 2001 状态且没有长跑，所以 checker 3.5 s）、
                                             # HB 壁 x 双向闩锁 10.33 s（4 臂 x 2001 加 5.54 那份对照文件，四次 app 调用
                                             # 合计 4.45 s，其中首次 2.51 s 含冷缓存；checker 直跑 4.96 s）、
                                             # HB 壁 x tau_0 x 指数 x 双向闩锁 23.35 s（4 臂 x 4 档 tau_0 x 3 副加
                                             # 5.59 那三副对照，六次 app 调用合计 9.41 s；checker 直跑 10.93 s）、
                                             # 瞬态源日程 x HB 壁 2.67 s（2 x 2001 加 5.59 对照，两次 app 调用 0.11 s；
                                             # checker 直跑 0.38 s）、
                                             # 停泊带下边缘 5.52 s（4 臂 x 3 档 x 2001，三次 app 调用 0.72 s；
                                             # checker 直跑 3.14 s）、
                                             # 瞬态屈服静止支 399.55 s——最后一条贵的不是求解器（5 次 app 调用共 0.31 s），
                                             # 是 checker 把 10 x 2001 步的带支递推在**速率空间**自己重解了一遍；
                                             # 5.50 把同一种重解改写到**壁应力空间**（单调，一次二分），于是同规模的
                                             # 14 x 2001 全轨迹复现只剩 5.12 ~ 9.69 s
```

单个验证算例可以直接手动跑，便于调试：

```bash
./build/src/3d/apps/chflow-simple-foam tests/cases/poiseuilleChannel/case.dict
python3 tests/check_poiseuille.py \
    "$PWD/build/src/3d/apps/chflow-simple-foam" \
    "$PWD/tests/cases/poiseuilleChannel/case.dict"

# 三维瞬态 PISO：app 逐步打印 t/uMax/cont/correctors 与终态剖面，
# check_piso.py 对比启动渠道解析级数（含离散壁面偏移的显式核对）
./build/src/3d/apps/chflow-piso-foam tests/cases/startedChannel/case.dict
python3 tests/check_piso.py \
    "$PWD/build/src/3d/apps/chflow-piso-foam" \
    "$PWD/tests/cases/startedChannel/case.dict"

# 湍流渠道摩阻回归：同一个 case.dict 在临时目录里被复制改写成两个模型分别跑，
# 断言壁面剪应力对体积力的积分平衡、首格 u+ 对 Spalding 对数率、剖面对称单调、
# 平衡壁面单元与 Ub/u_τ 对 Pope 光滑平板律
python3 tests/check_turbulent_channel.py \
    "$PWD/build/src/3d/apps/chflow-piso-foam" \
    "$PWD/tests/cases/turbulentChannel/case.dict"

# 同一套断言换一种压力-速度耦合再跑一遍：case.dict 同时带 piso{} 和 simple{}，
# app 各读自己的块，剖面一律从 vtu 质心取，所以两个耦合比的是同一组数
python3 tests/check_turbulent_channel.py \
    "$PWD/build/src/3d/apps/chflow-simple-foam" \
    "$PWD/tests/cases/turbulentChannel/case.dict"

# VOF 多相首增量：速度场是**给定**的（uniform 平移 / rotation 刚体旋转），
# app 逐步打印相体积、α 上下界、质心、界面 0.5 穿越位置与界面带的均值 κ；
# check_vof.py 把 case 复制到临时目录运行、从 vtu 独立复算同一组量，
# 并按各自的闭式性质断言（守恒界、平移速度、曲率 −1/R、Laplace 跃变 σ/(ρR)）
./build/src/3d/apps/chflow-vof tests/cases/vofRotation/case.dict
python3 tests/check_vof.py \
    "$PWD/build/src/3d/apps/chflow-vof" \
    "$PWD/tests/cases/vofRotation/case.dict"

# 静止液滴（动量模式）：解出的拉普拉斯压差对 σ/(ρR)、伪速度幅值、镜像对称
./build/src/3d/apps/chflow-vof tests/cases/vofDroplet/case.dict
python3 tests/check_vof.py \
    "$PWD/build/src/3d/apps/chflow-vof" \
    "$PWD/tests/cases/vofDroplet/case.dict"

# 椭圆液滴振荡（两相耦合）：`surfaceTension{ advect true; }` 让解出的 φ 输运 α、
# 每步从新形状重建力；checker 测 aspect 穿越周期对密度匹配的 2D Lamb 值
# 2π√(ρR³/(3σ))（R = √(ab)），并断言质量守恒、振幅双轴衰减（约 80 s）
./build/src/3d/apps/chflow-vof tests/cases/vofOscillating/case.dict
python3 tests/check_vof.py \
    "$PWD/build/src/3d/apps/chflow-vof" \
    "$PWD/tests/cases/vofOscillating/case.dict"

# 密度基可压缩 Euler（Sod 激波管）：app 逐步报质量/动量/能量/壁冲量/压力上下界与 CFL，
# check_euler.py 用 Python 重放同一套 HLLC 递推（含同一个波速估计）与写出的场逐格比对，
# 再对精确 Riemann 解的平台、锋的半跃变交叉位置与全场 L1 的加密收敛率断言
./build/src/3d/apps/chflow-3d-compressible tests/cases/sodShockTube/case.dict
python3 tests/check_euler.py \
    "$PWD/build/src/3d/apps/chflow-3d-compressible" \
    "$PWD/tests/cases/sodShockTube/case.dict"

# 三维区域分解首增量：同一个块按 1/2/3/5/7/8/16 域各切一次、逐域跑显式算子，
# 每次域数写一份**全局单元序**的 vtu。check_parallel.py 从 vtu 的几何自己重跑
# 同一套坐标中位数分解，断言分区性与面数重装配、分解算子与串行逐位相等、
# halo 值不漂移、归约等于串行的 sum/min/max，以及不同域数写出的场一模一样
./build/src/3d/apps/chflow-3d-parallel tests/cases/parallelDecomposition/case.dict
python3 tests/check_parallel.py \
    "$PWD/build/src/3d/apps/chflow-3d-parallel" \
    "$PWD/tests/cases/parallelDecomposition/case.dict"
# 同一套契约跑在 Gmsh 四面体网格上（三角形面、无非正交对）：只断言"切割不改变
# 它"，因为这张网格自身的格式误差就吃不下线性场（见 caseTet.dict 头注释）
python3 tests/check_parallel.py \
    "$PWD/build/src/3d/apps/chflow-3d-parallel" \
    "$PWD/tests/cases/parallelDecomposition/caseTet.dict"

# 分区 1D↔3D 耦合首增量：管网的一条边（X）由三维 SIMPLE 块承担。app 按 couple.seeds
# 里的每个初值电阻各跑一整趟外层回路，逐趟打印 R/q/dp/残差与块的自测量
# （进出通量、逐截面平面通量极差、探针窗梯度与散布、剖面峰值），最后一趟的场写 vtu。
# check_coupled_1d3d.py 断言接口恒等式（网络分支流量 == 块 fluxIn 到 1e-12、
# 残差 == |dpNetwork−dpBlock|/dpBlock）、节点衡算与压力链闭合、块内守恒、
# 已发展段 dp/dx 对平面 Poiseuille、vtu 剖面逐格对抛物线，以及闭式工作点与初值无关性。
# 约 20 s（三个种子各带一次 1000 格 SIMPLE 扫）
./build/src/coupling/apps/chflow-coupled-1d3d tests/cases/coupled1d3d/case.dict
python3 tests/check_coupled_1d3d.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d" \
    "$PWD/tests/cases/coupled1d3d/case.dict"

# 耦合边流变（§5.63）：同一根 coupledResistance 管上声明 flowIndex/consistency/
# yieldStress。1D parser 把声明搬进 partnerRheology 并清零自己的壁面律字段（损失仍是
# 伙伴给的 R·q，不重复计费），四个 coupling apps 转发给 DetailBlock；3D 侧新增
# GeneralizedNewtonianViscosity 实现现成的 IEddyViscosity 逐面钩子，SIMPLE 方程装配未改。
# check_coupled_rheology.py 从一份 case 派生 5 个稳态臂（未声明 / 声明但牛顿 / 剪切稀化 /
# 剪切增稠 / HB 屈服）+ 1 个瞬态臂 + 4 个拒收 case，306 条判据：声明但牛顿档与未声明档
# 求解日志逐位相同、nuWall/gdotWall/tauWall 独立重建 K*gdot^(n-1)+tau0/gdot、
# |probeSlope|·H/2 对 tau_w、dp 对同一本构的平面渠道闭式（含塞流积分）、逐趟界面残差、
# 1D Re 仍逐位是 rho*v*d/mu、瞬态每步动量行、四条拒收边界。
# §5.63 还把 4 个既有 checker 的耦合边负例（case_rejectCoupled*.dict）从
# "必须拒收"改成"接受且与无声明孪生档逐位相同"：verify.steadyPowerLawRheology /
# verify.steadyBinghamYield / verify.transientPowerLawRheology /
# verify.transientYieldStaticHold。只改判据与 case 头注释，不改物理。
# 约 110 s（5 稳态臂 63 s + 瞬态臂 43 s）
python3 tests/check_coupled_rheology.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d" \
    "$PWD/tests/cases/coupledRheology/case.dict" \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient" \
    "$PWD/tests/cases/transientCoupled1d3d/case.dict"

# 瞬态分区耦合：同一个网络、同一块网格改成从静止被突加的注入速率驱动。app 逐时间层
# 打印 q/dpNetwork/dpBlock/R/残差/趟数/newton 与块收敛标志，外加 name-keyed 的节点压力与
# 管流序列；vtu 是最后一个时间层、最后一次子迭代的块场（步内解多次，只有最后一次可回放）。
# check_transient_coupled_1d3d.py 断言每步的动量恒等式
# dp = R·q + M·Δq/dt、三段流量衡算、首步对"突加速率按惯量分"闭式、由整条爬升回归出的
# τ 与 q∞、末端距工况点（且必须 >0），以及块侧守恒/已发展梯度/抛物剖面。
# 约 28 s（60 步共 95 次 1000 格块解）
./build/src/coupling/apps/chflow-coupled-1d3d-transient tests/cases/transientCoupled1d3d/case.dict
python3 tests/check_transient_coupled_1d3d.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient" \
    "$PWD/tests/cases/transientCoupled1d3d/case.dict"

# 多边分区耦合：一张管网同时挂两块三维详细区（X1 在并联支里、X2 在串联位置上）。
# app 按"初值集"逐个跑整趟外层回路，逐趟**对每条边**打印 R/q/dpBlock/残差/趟数与块收敛位，
# 再给 name-keyed 的逐边 dp/q、节点压力与管流、每边的块自测；每边各自的 vtu 只在最后一组
# 初值写。check_coupled_1d3d_multi.py 断言逐边接口恒等式、纯 1D 侧的衡算与压力链、
# "一趟一次网络解 + 每伙伴一次解"（数日志里 SIMPLE iter 1: 的条数）、边互不借力
# （X2 首趟压降在三组之间逐位相同）、三组初值逐边落同一点、串并联闭式工况点，
# 以及**两块各自**的守恒/已发展梯度/抛物剖面对各自的间隙。
# 约 65 s（1000 + 4800 单元、三组初值共 8 趟 = 16 次块解）
./build/src/coupling/apps/chflow-coupled-1d3d-multi tests/cases/coupled1d3dMulti/case.dict
python3 tests/check_coupled_1d3d_multi.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-multi" \
    "$PWD/tests/cases/coupled1d3dMulti/case.dict"

# 湍流伙伴块：同一套分区接口，X 换成 k-omega SST 的平面渠道——块的 dp/q 随速率变，
# 接口要追的不动点第一次不再是"一次就中"。同一个种子电阻跑两趟臂：未加速的
# block Gauss–Seidel 与准 Newton，差别只在系数怎么取。
# check_coupled_turbulent.py 断言接口恒等式（网络记的压降 == R·q、报的压降 == 伙伴自己
# 的、残差 == 两者的差）、把每趟的 (q, dpBlock) 拿来拟合 dp = c·q^m 得 m = 1.74（同一套
# 拟合放在层流块上是 1.06）、伙伴 apparent resistance 随速率变 2.12 倍、块内动量闭合
# |dp/dx|·H/2 == ρ·u_τ²（对 1.1%）、y+ 在壁函数带内、涡粘是分子粘的 153 倍、
# 剖面比抛物线平（u_中线/U_块 = 1.094 vs 1.5）、Pope 光壁律反推的来流速度对实测 −4.4%、
# 两臂落同一工况点（R 差 3.7e-4 相对）且加速臂 5 趟 ≤ 未加速臂 8 趟，
# 外加"块解次数 == 两臂趟数之和"（13 == 13）。约 42 s
./build/src/coupling/apps/chflow-coupled-1d3d tests/cases/coupledTurbulent/case.dict
python3 tests/check_coupled_turbulent.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d" \
    "$PWD/tests/cases/coupledTurbulent/case.dict"

# 瞬态湍流伙伴：同一个 RANS 渠道放进时间层内的接口，两条支路都随速率变——块是
# q^1.68（爬升两端的 (q, dpBlock) 直接拟合出来），旁路 1D 管是 q^1.73（Churchill 在它
# 自己 Re≈1.1e4 上的割线）。于是网络自己的时间常数不再是"两条割线之和"：
# τ = (M_X+M_By)/(m_X·R_X + m_By·R_By) 用斜率，层流算例里 m≡1 所以割线与斜率同值。
# check_transient_coupled_turbulent.py 断言每步动量行 dp == R·q + M·Δq/dt（到 7e-11）、
# 旁路那一行对**外面文献里的** Churchill 割线闭合（7e-14，等于把 1D 摩阻律本身也校了）、
# 拟合整条爬升得 B = 0.7842 对闭式 0.7835（−0.34%）而割线闭式读 0.3990 s（+62%，
# 偏差正是阻力加权的指数 1.70）、首步按惯量分（9.5e-4）且离阻力分工况点 2.3%、
# 末端落在拟合渐近线上（3.8e-5，与自身 τ 的指数尾巴 3.3e-5 同量级）、
# 每种子 44 趟对 27 步且块解次数 == 两种子趟数之和（88 == 88）、两组初值落同一点
# （速率 2.8e-5、系数 6.2e-5、τ 3e-5），外加块侧 y+/动量闭合/涡粘/剖面/Pope 与入口段。
# 约 10 s（320 单元、两组初值共 88 次块解）
./build/src/coupling/apps/chflow-coupled-1d3d-transient tests/cases/transientCoupledTurbulent/case.dict
python3 tests/check_transient_coupled_turbulent.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient" \
    "$PWD/tests/cases/transientCoupledTurbulent/case.dict"

# 瞬态多边同台：稳态多边算例的拓扑与两块网格，换成时间步进驱动（12 步 × dt 0.05）。
# X1 在并联支里 ⇒ 速率是动量平衡的结果，首步按惯量分、再向电阻分弛豫；X2 在串联位置上 ⇒
# 不可压网络每步把 Q 钉给它，它的界面只谈压降，但启动的账 M·Q/dt 仍在它自己的动量行里。
# 每趟除了逐边的 R/q/dpBlock/残差，还打印 **worst 的下标**，所以"停在最差边"可从日志重算。
# check_transient_coupled_multi.py 断言逐边逐步身份、两边各用自己的 M 闭合动量行
# （X1 3.0e-11 / X2 1.8e-13）、串联边速率每步 == Q、界面日程（每步恰 1 趟被采纳、
# 15 趟只有一条边超差、块解 == 趟数 × 2）、边互不借力（X2 首趟残差两组同为 0.6601）、
# 两组初值逐边落同一点（1.0e-7）、首步纯惯量分工（4.6e-12）、爬升对上**隐式欧拉自己的**
# `1/(1+dt/τ)` = 0.841045（实测 0.840671）而不是连续的 exp(−dt/τ)，外加两块各自的场。
# 约 235 s（1000 + 800 单元、两组初值共 41 趟 = 82 次块解）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMulti1d3d/case.dict
python3 tests/check_transient_coupled_multi.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMulti1d3d/case.dict"

# 多边 × 区域惯量归属同台：同一张管网、同样两条边，只把两臂的惯量来源换成
# `inertias ( network partner )`，种子逐位相同。两条边被"名义孔径"错估的幅度不同
# （X1 +27.3%、比值 1.273240；X2 +69.8%、比值 1.697653，都 == 实截面/内切圆到 1e-12），
# 所以一个全局修正因子做不出这两行。
# check_transient_coupled_multi_inertia.py 断言每边每臂装上的 M == 该来源（1e-12）而块自测的
# inertance 两臂逐位相同、界面系数两臂每步同一个数（X1 7.4e-4 / X2 1.3e-6）且 dpNetwork == R·q
# 到 0.0 Pa（归属不进界面系数）、串联边速率不动而首步节点压降差 == ΔM·Q/dt
# （446.497744 对 446.497745 Pa，启动账 1086.5 对 640.0 Pa）、并联边首步按各自 M 分工
# （4.6e-12 / 5.9e-13）、τ 0.264555 对 0.240207 s 且衰减各自对上 1/(1+dt/τ)、
# 目的地与归属无关（0.0012%）而 12 步末端仍差 +1.23%，外加块解 == 趟数 × 2。约 365 s
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiInertia/case.dict
python3 tests/check_transient_coupled_multi_inertia.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiInertia/case.dict"

# 湍流 + 多边 + 瞬态三者同台：§5.25 的两条边同时换成 §5.23 那类 RANS 伙伴，两臂只差
# `accelerates ( 0 1 )`（种子逐位相同 ⇒ 臂间差别只能归给接口算子）。代理管孔径取区域的
# 等效截面 `d = 2√(HW/π)` ⇒ §5.26 的归属自由度被中和，本算例只剩伙伴的本构可测。
# check_transient_coupled_multi_turbulent.py 断言两臂 24 个边×步上 dpNetwork == R·q 到 0.0 Pa、
# 每边动量行对发表的 Churchill 弦闭合（7.6e-14）、X1 的律 dp ∝ q^1.682 而 X2 读不出指数
# （速率被拓扑钉死）、斜坡的恢复项是**斜率和**（+0.23%）而不是弦（+64%）、`worst` 下标
# 53 趟全 0、块解 == 趟数 × 2、两臂末端同一状态与同一 τ，外加两块各自的壁面律
# （y+ 41.7 / 58.4、|dp/dx|·H/2 == ρu_τ²、Pope 光壁律）。约 90 s（106 次块解 / 53 趟）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulent1d3d/case.dict
python3 tests/check_transient_coupled_multi_turbulent.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulent1d3d/case.dict"

# 湍流 × 区域惯量归属同台：把上一算例中和掉的那格自由度放回去，并且放歪——两根代理管都写成
# 图纸上的名义圆且歪得不同（X1 1.273240×、X2 1.696946×），于是 `network` 与 `partner` 装进
# 动量行的是两个真的不同的数，而不存在一个全局修正因子能同时复现两臂。两臂只差
# `inertias ( network partner )`、种子逐位相同、都不加速。
# check_transient_coupled_multi_turbulent_inertia.py 断言归属只挪惯量账：两臂 24 个边×步上
# dpNetwork == R·q 仍到 0.0 Pa；每边行用本臂那份 M 闭合（4.0e-11 / 8.8e-11）；τ 是
# (M_X1+M_By)/斜率和，孔径臂 0.277231 s 对预测 −0.232%、区域臂 0.240829 s 对 +0.231%，弦式两臂
# 都读高 63.5%；τ 之比 1.15115 == 惯量之比 1.15649（−0.46%）；目的地不含 M（渐近速率差 0.0113%，
# 末端速率却差 1.492%）；串联边两臂逐位相同（1.4e-20）而首步踢出该臂的 M·Q/dt（差值 == ΔM·Q/dt）；
# 伙伴指数两臂读出同一条律（1.6762 / 1.6820）。约 45 s 串行（ctest -j8 下 142 s；114 次块解 / 57 趟）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentInertia/case.dict
python3 tests/check_transient_coupled_multi_turbulent_inertia.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentInertia/case.dict"

# 按边混合惯量归属：同一张管网、同样两块 RANS 伙伴、同样两个写歪的名义孔径，但把"这条边的
# 惯量从谁那儿取"从**臂**交到**边**——`couple{ inertias ( network partner network ) }` 之外再加
# `couple{ inertiaFrom{ mixed{ X1 partner; X2 network; } } }`，三臂种子逐位相同、都不加速。
# check_transient_coupled_multi_mixed_inertia.py 判的是"per-edge 表 = 两次整体选择的叠加，不是
# 第三种读法"：mixed 的 X1 全部 12 步与 region 臂逐位相同、X2 的首步节点压降与 bore 臂同一个数
# （17910.892818 Pa 对 region 臂 11638.381827 Pa，差 == ΔM·Q/dt）；界面仍是弦（三臂 36 个边×步
# dpNetwork == R·q 到 0.0 Pa）；每边动量行用本边那份 M 闭合（4.0e-11 / 8.8e-11 / 8.8e-11 并联、
# 2.8e-14 / 3.1e-14 / 3.1e-14 串联）；τ 跟着会爬那条边（mixed 0.240829 s == region 臂，逐位）；
# 目的地与来源无关（渐近速率差 0.0113%）；伙伴指数三臂同一条律（1.6762 / 1.6820 / 1.6820）；
# 日程 30 / 27 / 27 趟、块解 == 趟数 × 2、`worst` 下标 84 趟全 0。十族 132 条，约 147 s 串行
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiMixedInertia/case.dict
python3 tests/check_transient_coupled_multi_mixed_inertia.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiMixedInertia/case.dict"

# 湍流 × 归属 × 准 Newton 三者同台：把这条链上一直各自单独扳动的两条轴**交叉**成完整 2×2——
# 同一张管网、同样两块 k-ω SST 伙伴、同样两个写歪的名义孔径，四臂只差
# `inertias ( network network partner partner )` 与 `accelerates ( 0 1 0 1 )`，种子逐位相同。
# check_transient_coupled_multi_turbulent_inertia_newton.py 从两个方向判独立：算子进不了归属的表
# （同来源两臂装上逐位相同的 M：10185916.36 / 6363546.246 与 8000000 / 3750000 kg/m4）；界面在两种
# 算子下仍是弦（四臂 96 个边×步 dpNetwork == R·q 到 0.0 Pa，每边行按本边那份 M 闭合 4.0e-11 /
# 8.8e-11）；算子只买到日程（同对两臂末端速率差 3.2e-07 与 2.0e-06 相对、τ 差 8.9e-05 与 8.5e-04
# 相对、串联首步踢出差 0.044 与 0.053 Pa，加速臂各少 1 趟）；归属照旧只管惯量账（τ 预测 −0.232% /
# −0.247% / +0.231% / +0.158%，τ 之比 1.15115 == 惯量之比 1.15649 在两对里都成立，目的地 spread
# 0.0119%，伙伴指数 1.6762 / 1.6763 / 1.6820 / 1.6816）；准 Newton 作用在界面整体那条被钉死的边
# 上，残差地板 2.300e-07 → 1.670e-05 / 2.001e-05（73× / 87×）但仍在 tolerance 之内。十族 172 条，
# app 四臂 93 s、checker 串行约 99 s（224 次块解 / 112 趟界面 / 48 个谈拢的时间层）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentInertiaNewton/case.dict
python3 tests/check_transient_coupled_multi_turbulent_inertia_newton.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentInertiaNewton/case.dict"

# 湍流 × 归属 × 准 Newton × 按边：三条轴同台（8 臂，内核零改动）
# 同一张管网、同样两块 k-ω SST 伙伴、同样两个写歪的名义孔径，八臂只差 inertias、
# accelerates ( 0 0 0 0 1 1 1 1 ) 与四张 inertiaFrom{} 表（mixX1/mixX2 两张镜像表各跑两种算子）。
# check_transient_coupled_multi_three_axis.py 判的是「按边归属在准 Newton 下逐位存活」：装上 M 的表
# 八臂一致（X1 10185916.36 / 8000000、X2 6363546.246 / 3750000 kg/m4），192 个边×步 dpNetwork == R·q
# 到 0.0 Pa，每边行按本边那份 M 闭合（4.0e-11 / 8.8e-11、串联 3e-14）；同表两臂 X1 速率历史
# QN 四对与 GS 区域对逐位相同（0.0e+00），GS 孔径对 2.9e-10，异表 4.801e-02——分离度 1.6e8x..4.8e14x；
# 被拓扑钉死的 X2 用节点压降量（同表 1.1e-09/1.7e-09 GS、4.9e-07/7.5e-07 QN，异表 3.5e-01/5.4e-01）；
# 两张镜像表在两条边上都分开，混法臂不是折中；tau 由做斜坡那条边的来源决定（mixX2 == bore 0.277231 /
# 0.277255，mixX1 == region 0.240829 / 0.241034，差 4.0e-09..0.0 相对），串联首步踢差 6272.510991 /
# 6272.502249 Pa == deltaM*Q/dt；算子只买到日程（同对末端速率差 9.702e-07 / 4.428e-06，加速臂各少 1 趟）；
# 八臂弛豫到同一速率（渐近 spread 0.0119%）；被钉死那条边的残差地板 2.300e-07 → 1.670e-05 / 2.001e-05
# （73x / 87x / 87x / 72x）而谈判那条分毫不动。十族 373 条，app 八臂 275 s、checker 串行 267 s
# （448 次块解 / 224 趟界面 / 96 个谈拢的时间层）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiThreeAxis/case.dict
python3 tests/check_transient_coupled_multi_three_axis.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiThreeAxis/case.dict"

# 两条耦合边都在并联位置：八臂按边表换拓扑重跑（4 张表 × 2 个算子，内核零改动）
# 只改拓扑（Src -> A -> { X1 || X2 } -> B -> Snk，无旁管）与两个名义孔径，inertias /
# inertiaFrom{} / accelerates / starts / 容差 2e-4 / 6 趟预算整张照抄上一例。于是接口里有两个
# 都要谈的速率，三件事第一次可测：
#   1) worst 下标换手——279 趟逐趟等于两残差 argmax，边 1 当最差 4~13 次（上一例 224 趟全 0）；
#   2) 一个节点压降对两份动量行——同一个 J1-J2 各按本边那份 M 闭合（1.1e-10 / 7.6e-11），
#      从外面重算每层分流 q_e=(dp+M_e q_prev/dt)/(R_e+M_e/dt)、sum q_e=Q_IN 到 4.7e-12，
#      任一边换错账立刻偏 0.047~0.108（10 个数量级）；
#   3) 斜坡 τ 现在由 M_X1+M_X2 这一份分子排序——四张表 1.65495e7/1.43635e7/1.39359e7/1.175e7 kg/m4
#      对 τ 0.29115/0.25287/0.25075/0.21374 s，两算子内严格单调、τ 比 == 分子比到 3.40%，
#      两条混合表严格夹在两张整表之间（上一例混合臂 τ 与同源整表差 4e-9，因为只有一条边做斜坡）。
# 目的地仍不含惯量：八臂 0.409319~0.409364（spread 0.0206%）== 本臂两边幂律（X1 1.6601~1.6861、
# X2 1.6703~1.6842）的电阻分流 ≤0.070%；单一指数穿 12 层偏离 ≤0.258% 行程。「没有乘客」：八对
# 地板比 0.044~25.6 且方向不一致（5/8 平谈更差），最差地板 ≤ 容差 21.1%——对照上一例被钉死那条边
# 的单向 73x~87x。十族 332 条，含跑 app 的 ctest 单条 377.5 s（279 趟界面 / 558 次块解 / 96 层谈拢）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiBothParallel/case.dict
python3 tests/check_transient_coupled_multi_both_parallel.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiBothParallel/case.dict"

# 块内真实非定常 du/dt 进界面：§5.32 的拓扑与两个网格整张照抄，但这一例动了内核
# （SIMPLE 动量方程可选装配 fvm::ddt(u_old, dt)；界面的刷新与准 Newton 残差改用块报压降的
# 摩擦部分；加速度记进 1D 动量行并评分；归属新增第三本账 slab = rho*V_probe/A^2）。
# 六臂 = 准稳态控制（ddts 把块字典的 Euler 覆盖回 steady，于是它就是上一例的 bore_gs）
# + 三个来源 x 平谈 + 绘图/slab 各再跑一遍准 Newton。要点：
#   1) 界面只谈弦——R*q 落在 dpFriction 上最差 1.987e-04（容差 2e-4 的 1.0 倍），离整段压降
#      差 22x ~ 3039x 容差（最坏一次 61% 压降）：这就是"用整段压降刷新"那条 a-c/R 死路；
#   2) 一条不可迭代的评分 rowError = (R*q + M_row*dq/dt - dpBlock)/dpBlock——沉降层（4~12）它
#      恰等于 (M_row - M_probe)*(dq/dt)/dpBlock（绘图 0.968~1.015、区域 0.971~1.012），slab 臂
#      同一式预测 0、实测残留 <=0.23% 压降。归属第一次有「误差」而不是「路」；
#   3) 块给自己评分 chi = inertanceDynamic/inertanceProbe：起步层 92.8%~94.6%，第 4 层起
#      |chi-1| <= 2.8% 且五臂同 (边,层) 之比 <1.05，第 2~3 层 >5%（那两层集总块说不了话——
#      X2 第 2 层通量已翻负而场仍加速，是「块不是速率函数」的直接读数）；
#   4) 首层摩擦弦是准稳态弦的 2.46~3.05 倍（起始壁面剪），第 3 层起回到 ±3.8%——代价由下一层
#      用趟数付，maxSubIterations 因此从 6 提到 10（6 趟会让 slab_gs 第 2 层以 2.2e-3 谈不拢）；
#   5) 目的地仍由弦决定但位置随归属变：六臂末端 X1 分流 0.400006~0.406190 对电阻分流 0.409367
#      偏 0.78%~2.29%，三来源末端三条不同速率（两两差 >1e-9）。
# 写这一类判据时的坑：settled 层的挑选必须显式开窗——整段 12 层取中位会把绘图与区域在 X2 上压成
# 1.99% vs 1.96%（不可分辨）；chi 那种「两个小数之比」在 dq 趋零的层会失去意义，只能报区间不能定带。
# 十族 378 条，含跑 app 的 ctest 单条 ~36 s（259 趟界面 / 518 次块解 / 72 层谈拢）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledUnsteadyPartner/case.dict
python3 tests/check_transient_coupled_unsteady_partner.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledUnsteadyPartner/case.dict"

# profile 级界面映射（稳态多边）：§5.20 的拓扑、两块网格、relax、容差 1e-4 与 12 趟预算整张照抄，
# 但界面上交的不只是一个速率：每条边可选 profile ( flat | parabolic | power1over7 )，且 starts{}
# 的每一臂各自命名。这一例动了内核三处：volVectorField 的 fixedValue 支持逐面覆盖表、
# buildInletProfile 按面积分归一化使 ∫u dA = q 精确、HbyA.adoptFaceValuesFrom 让预测通量也看见
# 这份逐面 prescribe（缺它时抛物入口的 X2 连续性地板钉在 6.46e-2、400 扫不收敛）；外停机的 settled
# 参照由"历史最小值 1.01 倍"改为"两次扫之前"（冷启动前两个二循环会先跌破地板，热启动会从地板下方
# 起步，历史最小值被两件事同时污染）。六臂 = 塞状 x2 起点系数 + 抛物 + 1/7 + mixed（一次运行两条边
# 喂不同形状）+ 抛物走准 Newton。十族断言的要点：
#   1) 速率活过映射——fluxImposed == 网络要的数（1.2e-15 相对）、与 patch 的 fluxIn 逐位相同、
#      每个横截面通量 min/max 差 <=2e-5；
#   2) 形状就是声明的那个——块从自己的边界面上读回 betaImposed，对 checker 独立重算的中点求和
#      一致到 1e-9（抛物 20 面板 1.1970121617640874），flat 离 1.0 <=1.5e-15；
#   3) 还剩多少路——塞状 X1 从 1.0245788 单调爬到平台 1.1955544，1/7 从 1.0558 走同一条路但更短，
#      抛物入口在承载边界平面之外漂移 <=0.16%（平移不变本身就是形状被真喂进去的证据）；
#   4) 平台是区域的——三种入口的发展 betaOut X1 相对极差 2.9e-5、X2 1.4e-6，都贴连续 Poiseuille 1.2；
#   5) 阻力有自己的物理——对 12*mu*L_meas/(H^3*W)：抛物 -1.81%/-2.48%（壁面那一面的梯度只跨半个
#      单元，离散性质）、塞状 +22.97%/+14.08%、1/7 +19.03%/+11.27%，且 R_flat > R_17 > R_parab；
#      probe 窗内散射 X1 降 104x、X2 降 23x；
#   6) 工作点因此会动——只换形状，并联分流 +8.20%（4.9905859e-4 -> 5.3998308e-4），串联边恒等于
#      注入量（7.9e-15）；两域 Re 399~432 与 427，层流，所以臂间差异只可能是入口；
#   7) 不变量——起始系数相差三倍落同一个 R（4.4e-5）、mixed 的 X1 与 parab 的 X1 逐位相同
#      （dpPa 1961.8431819560024）、平谈与准 Newton 差 2.0e-7 只买到趟数（3 == 3）；
#   8) 形状进的是预测器——抛物入口的每边地板 1.997e-3 / 3.597e-3，比塞状同一边的 4.813e-2 /
#      5.430e-2 低 24x / 15x（§5.5 的"预测地板是 HbyA/压力通量固有失衡量"由此多一个形状依赖）。
# 写这类判据时的坑：per-edge 地板只能按每块 SIMPLE 序列的最后一条归位（块报告在所有 SIMPLE 行之后
# 才打印），用 len(floors) == 2*iterations 与 floorLens[slot]+1 == 块 sweeps 双向锁死；统一阈值是错
# 的——塞状臂坐在 4.8e-2 上本来就合法，只有形状臂该到 2e-3。
# 十族 296 条，含跑 app 的 ctest 单条 166.17 s（六臂、1500 + 3200 单元）
./build/src/coupling/apps/chflow-coupled-1d3d-multi tests/cases/coupled1d3dProfileMap/case.dict
python3 tests/check_coupled_1d3d_profile_map.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-multi" \
    "$PWD/tests/cases/coupled1d3dProfileMap/case.dict"

# profile 级界面映射（瞬态多边）：§5.25 的拓扑、驱动、12 步 × dt 0.05、容差 1e-4 与 8 趟预算整张
# 照抄，两条 duct 的几何也是上一例那两块（沿管向更粗：30x20x1 与 25x20x1），只是界面上每层都重交一次
# 形状。动的是耦合侧：瞬态 app 解析 profile（与稳态逐字同款，一个词=两边同形状、一个列表=逐边）、
# DetailBlock::profileLedger() 把形状台账交给 EdgeExchange，于是 transient step n: edge i 每层带
# inletShape/fluxImposed/betaImposed/betaIn/betaOut（fluxImposed 是界面写下的那个条件自己的面积分，
# 它等于该层 q 就是把"交了形状"和"交了另一份流体"分开的那条判据）。稳态那条路径直接读 measure()，没动。
# 五臂 = lo（两边塞状，对照）+ 抛物 + 1/7 + mixed（X1 抛物、X2 塞状）+ 抛物走准 Newton，全部用未加
# 形状那一版引用的同一对起点系数（X1 1.25e6 / X2 5.3e6，各低于自己 duct 值三分之一）。
# 十族 131 条断言的要点：
#   1) 形状是逐层的声明——5 臂 x 12 层 x 2 边 = 120 条边状态每条都报出自己那层的形状；
#   2) 速率活过映射，逐层都活——fluxImposed == 该层 q 最差 6.66e-16，betaImposed 每层回到该形状的
#      解析锚点（1 / 1.197012161764087 / 1.013430581098372）到 1.44e-15 且第 1 层 == 第 12 层；
#   3) 交进去的形状 != 区域内读到的形状——betaIn/betaOut 取**首末 cell 截面**而不是边界面，所以
#      塞状也读 1.09592（对 betaImposed 1）；抛物以 1.19701 到达、首截面 1.20274，是被往下修；
#   4) 两端不对称——末截面五臂极差 X1 6.86e-6 / X2 1.01e-8，首截面分开 0.107 / 0.077；
#   5) 去处与快慢同时随形状动——X1 分流末端 +0.4868% / +2.7647%，而 tau 里的 R_X1 就是这个落定的数，
#      于是 0.258135 s -> 0.266619 s（+3.29%）；每臂轨迹对上自己的离散惯量解 1/(1+dt/tau)；
#   6) 逐边到逐位——mixed 的 X1 与 parab 的 X1 逐位相同（X1/By 是电流分配器，下游动不了它），mixed
#      的 X2 跟着 lo 到 2.58e-6 而不是 0，其因由是界面在容差处停下、每层热启动取决于此前花了几趟；
#   7) 准 Newton 在此容差下什么都没买到（与平谈同样 17 趟），形状与算子两轴仍可分。
# 写这类判据时的坑：形状是否"还在发展"不能按**声明**分支（flat 之外还有 1/7：它到达 1.0134 却从
# 1.1213 起步，同样在爬），只能拿该形状的解析锚点和区域自己的平台比；"首截面极差 >=0.1"这种带子是
# 从 X1 抄来的，X2 短、只有 0.077，正确的陈述是"首截面极差 >> 末截面极差"（这里 15566x / 7.6e6x）；
# 抛物臂的每米压降散布在 X1 塌 65x 而在 X2 只塌 1.22x，所以承担结论的是"每块自己的窗宽对自己的压降
# <=0.053%"，不是散布之比。
# 十族 131 条，五臂共 101 趟 / 202 次块解 / 60 个谈拢的时间层；单独跑整个 case 35.5 s，
# 在 -j8 的 ctest 里两次分别记 41.65 s 与 63.64 s（同批负载决定，墙钟由最慢那条决定）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledProfileMap/case.dict
python3 tests/check_transient_coupled_profile_map.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledProfileMap/case.dict"

# profile 级界面映射 × 湍流伙伴（稳态多边）：§5.20 的拓扑 + §5.27 的两条 kOmegaSST 渠道（X2 由 0.12
# 加长到 0.24 m，好让抛物形的入口段落进块里），驱动 6.0e-5 m3/s 进水（Re_bulk 4879~10000、y+ 21.3~32.6），
# 块侧 simple maxIterations 300 / continuityTol 3e-3，界面 tolerance 1e-3（不是 1e-4，理由见下面第 9 条）。
# 动的是内核一处新工具：BlockMeasure 多一列 pressureStations——压强报在**与 beta 同一组**截面上
# （CoupledNetwork.h:209、coupledNetwork.cpp:798、multiCoupledApp.cpp:358）。五臂 = lo + 抛物 + 1/7 +
# mixed（X1 抛物、X2 塞状）+ 抛物走准 Newton。十族断言的要点：
#   1) 声明与伙伴无关——fluxImposed == 网络要的数（1.0e-12）、betaImposed 对 8 面板中点求和回到 1e-9
#      （抛物 1.181719848566793、1/7 1.0103203595676944）、离连续极限 <=2%；
#   2) 平台属于区域——十个「臂x边」的 betaOut 落在 1.020266~1.021045（极差 7.79e-4），离层流 duct 的
#      6/5 差 14.9%。这一条是整个增量的支点：形状的价格不能从层流那一版换算；
#   3) 判据是符号不是满度——十对里 5 过满 / 5 欠满（层流那版抛物是欠满的），过满者 beta 沿程单调降、
#      欠满者单调升，升的还包括 1/7 的第一个内部截面（交进去 1.01032、读回 1.01384）；
#   4) 账用压强付——等截面 duct 里壁面剪切只能往下拉，所以 parab X1 从 171.39 Pa **升**到 391.10 Pa
#      （x=0.045，+219.7）再降回 15.53 就是入口段在还账；lo/p17 四条曲线逐站不升；
#   5) 摩擦没动、割线动了——probeSlope 对模型自己的壁面平衡 2*rho*u_tau^2/H 落在 0.943~1.107（十对全测），
#      而 dpPa 对「同一平衡铺满整条 L」的外推：塞状 1.0118/1.0104、1/7 0.9800/0.9835、抛物 0.2835/0.5370。
#      抛物 X1 还回的 393.96 Pa 是 rho*(beta_in-beta_out)*U^2 = 225.7 Pa 的 1.75x（X2 2.19x），剩下那
#      一份是入口段壁面剪切低于发展值，本例只记差额不拆开；
#   6) 工作点跟着割线走——并联分流 +54.4%（抛物）/ +3.32%（1/7），串联边 X2 每臂恒等于注入量（1e-12），
#      J1 衡算 1.2e-12，旁路 1D 管以自己的 Churchill 弦闭合，J3 五臂同值（<1e-6 Pa）；
#   7) 逐边归属与算子各自仍然成立——mixed X1 与 parab X1 分流差 7.5e-13、压降差 1.3e-4（那一份是块的
#      扫数 29 对 32）；mixed X2 回到 lo X2（1.5e-7）；GS 与准 Newton 差 4.2e-4 / 3.0e-4，只买到趟数；
#   8) 壁面函数仍然合法——y+ 21.3~32.6（>11 未进尾迹区）、nut/nu 46~146x，否则「形状之间可比」不成立；
#   9) 容差有实测地板——抛物 X2 每趟花满 300 扫（blockConverged 0），在被拓扑钉死的速率下逐趟报的
#      dpBlock 在 719.35~720.30 Pa 之间摆（1.3e-3 相对），所以 1e-3 是伙伴分辨率的下界，不是愿望。
#      另有一堵墙：**1.2e-4 下塞状臂和 1/7 臂照样收敛（7 / 6 趟），只有抛物臂死在第 1 趟**——串联块 X2 把
#      整个 1.2e-4 推进 4 mm 缝隙（3.75 m/s），残差第 40 扫到 1.9e48、第 47 扫场非有限；同一次运行的并联块
#      X1 只是花满 300 扫停在 0.0163498。relaxP 0.3->0.1 前 32 扫逐位不变（aP = min(relaxP, max(pRelaxStart,
#      0.15*pScale/pStep))，simpleSolver.cpp:308）、此后 X1 提前在 101 扫停下但死因不变；model laminar 死在
#      并联块第 8 扫，所以不是旋钮问题；横向 8->16 加密后 13 次块解全部花满 300 扫、残差在 0.059~1.5 之间摆，
#      界面把 X1 推到 1.13e-4 时块报出**负的**压降 -0.4367 Pa 而停机。
# 写这类判据时的坑：betaImposed 与 betaIn 是两件事（前者取边界面、后者取第一个内部截面），过满/欠满的
# 分支只能按两者之差判、不能按声明的形状判——同一份 parabolic 在层流与湍流两版里落在符号的两边；每米
# 压降的独立锚点必须取模型自己报的 u_tau（2*rho*u_tau^2/H，只算那两面 y 向壁，z 向是 symmetry 不承载
# 剪切），拿 duct 律去套就把要检验的东西当成检验工具；外推压降要用 wslope*L 而不是 wslope*(末站-首站)，
# 否则块两端各半个单元的那段入口/出口会被抹进摩擦里。
# 十族 272 条，五臂共 27 趟 / 54 次块解；单独跑整个 case 29.7 s
./build/src/coupling/apps/chflow-coupled-1d3d-multi tests/cases/coupled1d3dTurbulentProfile/case.dict
python3 tests/check_coupled_1d3d_turbulent_profile.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-multi" \
    "$PWD/tests/cases/coupled1d3dTurbulentProfile/case.dict"

# 形状 × 湍流伙伴 × 瞬态三者同台（5.37）：§5.36 的两条 kOmegaSST 渠道原封不动搬进时间层内
# （X1 20x8x2 over 0.2x0.005x0.005、X2 24x8x2 over 0.24x0.004x0.008，intensity 0.05 / relaxation 0.7），
# 拓扑是 §5.34/5.35 的那张 Src -> A -> {X1 || By} -> X2 -> B -> Snk，驱动 6.0e-5 m3/s 恒定注入
# （斜坡 = 管网自己对并联分流的惯性弛豫，不是时间表源），time dt 0.05 / endTime 0.4 -> 8 层，
# 块侧 sdt steady + simple 300/3e-3，界面 tolerance 1e-3 / maxSubIterations 12 / accelerates (0 0 0 0 1)。
# 动的是内核一处：§5.36 那个「压强报在与 beta 同一组截面上」的工具多了一条**时间轴**——
# pressureStations 的归约（最高站压强减首站 = pRisePa，以及它在哪一站 = pRiseX）现在发生在
# profileLedger() 里（CoupledNetwork.h:296），随每层每边的 EdgeExchange 一起交出去
# （coupledNetwork.cpp:1224、transientCoupledMultiApp.cpp:486）。此前 app 是在整臂跑完后再读一次块
# 的测量，于是「入口段还压」只有臂末那一条曲线；现在层 n 报的是层 n 那次解自己的还压，80 个
# （层 x 边）状态各自的高峰可断言。五臂 = lo + 抛物 + 1/7 + mixed（X1 抛物、X2 塞状）+ 抛物走准 Newton。
# 十一族 652 条的要点：
#   1) 符号律逐层成立，不是端点巧合——80 个状态的 betaOut 全在 (1.02, 1.0213)（实测
#      1.02027438~1.02105813，比层流 duct 的 6/5 低 14.9%）；入口比平台满的恰好 40 个、差的恰好 40 个，
#      且与 pRisePa > 0 的那 40 个**是同一批**，没有一个状态骑在边界上；
#   2) 仪器不造上升量——pRisePa 为 0 的状态报的站就是区域第一站（0.005 m，单调下降曲线的最大值即起点），
#      非 0 的永远在同一站（X1 0.045 m、X2 0.035 m）；末层的 pRisePa/pRiseX 与该次求解自己的压强表逐位相符；
#   3) 还的就是动量那一笔且逐层——抛物 X1 上升是 rho*(beta_in-beta_out)*U^2 的 0.9888~0.9980
#      （136.19 预测对 134.66 实测），X2 落在 1.1748~1.1754；比值不随斜坡漂。上升量对速率的指数 2.1659
#      （101.7514 -> 134.6606 Pa，速率只涨 1.1381x），parab/mixed/parab_qn 三臂给同一个数；被钉死的
#      X2 坐在一个比值上到 6.3e-4（层 1 与层 8 都是 333.4297 Pa）；
#   4) 形状仍是逐边声明且每层重交——80 个状态的 inletShape 与该臂所名一致，betaImposed 对 8 面板重算到
#      1.1e-15、跨 8 层散布 <=9.4e-16，fluxImposed 对截面积分最差 6.66e-16（形状用形状买，不是别的流体）；
#   5) 形状不进 1D 动量行——rowDrop == dpNetwork + M*dq/dt 在 5 臂 x 8 层逐位成立，M 就是 rho*L/A
#      （8e+06 / 7.5e+06），归属这根轴本例固定在管网侧；dpInertia 恒 0、dpFriction 恒等于伙伴报的全量；
#   6) 首层是纯惯性分流——M_By/(M_X1+M_By) = 0.427275，实测 lo 0.424732 / parab 0.438354 /
#      p17 0.426449（±3% 内），次序仍是抛物 > 1/7 > 塞状；
#   7) tau 第一次被形状按**倍数**改掉——lo 0.3869 s / parab 0.7942 s（比 2.0531）/ p17 0.3897 s
#      （+0.72%，即抛物那份移动的 0.7%）；§5.35 对层流 duct 只移动 3.29%，因为那里的割线本来就贴着落定
#      斜率，这里 X1 落定割线被改到 3.1x。教材式 (M_X1+M_By)/(R_X1+R_By) 给 0.5823/1.2163/0.6160 s，
#      绝对值一律高估 1.5053~1.5807x，但抛物/塞状的**比值** 2.0887 对实测 2.0531 只差 1.2%；
#   8) 可分性在时间里保持——parab 与 mixed 的并联边 7 字段 x 8 层逐位相同（0.0）；mixed 的串联边就是 lo 的
#      （速率/压降相对差 <5e-5、pRisePa 恒 0）；串联边速率在 40 个（臂 x 层）状态上都是注入量（<=1e-11）；
#      parab 与 parab_qn 差 <1e-4，算子改路径不改物理；
#   9) 时间表与那堵墙——五臂 16/17/13/17/18 趟、SIMPLE 解次数恰为 2x 趟数、最坏残差
#      2.673e-06~6.040e-04 全在 1e-3 内；唯一不在自己门上停机的是**抛物形的串联块**（parab 8 层里 4 层、
#      parab_qn 6 层花满 300 扫，且花满与 blockConverged 0 同一批），同一钉死速率上塞状/1/7 只花 7 扫、
#      塑形并联边 24 扫；mixed 的串联边是平的所以即使邻居抛物也快——墙属于「喂给钉死速率的形状」，不属于
#      串联拓扑。网络握着速率，所以未收敛伙伴表现为摆动的压降而不是失败，这就是容差 1e-3 的理由；
#  10) 壁面函数仍然合法——y+ 21.4~32.6、nut/nu 47~134x、Re_Dh 4966~10000，发展窗 dp/dx 对模型自己的
#      壁面剪力 2*rho*u_tau^2/H 闭合到 -5.7%~+10.1%，每块入口进/出口出 <=5.2e-12；
#  11) 价格仍在割线里——末层 X2 塞状 1538.4003 / 1/7 1428.4563 / 抛物 719.9844 Pa（形状占 53.2% 的压降，
#      而这条边的速率从未变过）；抛物 X1 报 138.11 Pa 对发展外推 389.47 Pa = 0.3546，缺的 251.35 Pa 是
#      动量项 136.19 Pa 的 1.85x（X2 0.5368 / 2.19x）；塞状与 1/7 报自己 fully-developed 身的 0.9801~1.0124。
# 写这类判据时的坑：凡是「两臂之比」的断言都要先确认分母不会是 0——平入口的 pRisePa 按定义恒为 0，
# 拿它做相对差会 ZeroDivisionError，而恒 0 本身就是那条边该断言的东西（改用绝对差 + == 0）；
# 一个物理量的**区间**不要写死成单一实测值（同一次运行的 8 层给 0.9888~0.9980），断言要落在区间上；
# 逐层序列的衰减比只在斜坡中段是干净的几何级数，末端几层被界面容差截住，取前 6 层拟合、并断言散布 <5%。
# 十一族 652 条，五臂共 81 趟 / 162 次块解；单独跑整个 case 64.5 s
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentProfile/case.dict
python3 tests/check_transient_coupled_multi_turbulent_profile.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentProfile/case.dict"

# 形状 × 块内真 du/dt 同台（5.38）：上面那案的拓扑、两块 kOmegaSST 网格、壁面函数、操作点与 8 个时间层原封不动，
# 唯一变的是表述与时间层——time dt 0.2 / endTime 1.6（§5.37 是 0.05/0.4），两块 dict 都写 schemes{ ddt Euler; }，
# couple{ ddts ( steady steady Euler Euler Euler ) } 把前两臂**反向**覆盖回 steady，于是这个 override 两个方向都走一遍。
# 五臂 = {flat, parabolic} x {steady, Euler} 这张 2x2 加一格第三形状：lo_steady / parab_steady（= §5.37 那两臂）/
# lo_euler / parab_euler（交叉本身）/ p17_euler（让「加速度不认识形状」不只是一次巧合），
# accelerates ( 0 0 0 0 0 )（界面算子不是本例变量）、tolerance 1e-3 / maxSubIterations 12、惯量归属固定在管网侧。
# 内核一行没改：src/coupling/src/coupledNetwork.cpp:1182 的 chord = unsteady ? dpFrictionPa : dpPa 与 :1194 的
# residual = |R*q - chord|/chord 早就决定了「界面只迭代摩擦那半区」，本例做的是把这句话变成可断言的数。
# 要问的问题是一根轴各自问不出的：区域的压降现在是两笔之和 dpBlock = dpFriction + dpInertia，而界面只谈第一笔，
# 那么形状住在哪一项里。
# 十一族 490 条的要点：
#   1) 加速度对形状是盲的、摩擦半区不是——真正在动的那层（层 1）实测惯性对 rho*V_probe/A^2*dq/dt 的比 X1 三形状
#      0.9623/0.9658/0.9630（散布 0.34pp）、X2 0.9792..0.9794（0.07pp）；同三个形状把摩擦半区移动 X1
#      -14.8%（层1）/-51.2%（层8）、X2 -53.2%。层 3-8 三臂全在探针 slab 预测的 2.31% 内，且抛物在每个这种层最接近
#      1（L3 1.0084 对 1.0231 … L8 1.0004 对 1.0049）；末解 inertanceDynamic/inertanceProbe 1.00491/1.00036/1.00443。
#   2) 界面照常停机，因为它停的那半边不含惯量——80 个状态 residual 全 <1e-3，而 residual == |R*q - chord|/chord
#      在同一批状态到 0.0e+00；同处它们的惯量那半区最多 disagreement 62.7%。
#   3) 行的惯量在块也开始付费后被记了两次——rowError 层 1 = X1 53.9/62.7/54.7%、X2 44.1/56.0/45.3%（lo/parab/p17）；
#      末层 X1 仍差 12.3（抛物）/1.4（1/7）/1.0（塞状）%，被钉死的串联边降到 0.0051%（把它翻倍的那项停止运动了）。
#      准稳两臂 rowError 恒 0，且 rowDrop - dpBlock == M*dq/dt - (dpBlock - dpNetwork) 到 2.3e-14 Pa。
#   4) 行误差恰好是那笔归属算术、不含别的：rowError*dpBlock == (R*q - dpFriction) + (M*dq/dt - dpInertia) 在 48 个
#      步进状态到 2.3e-13 Pa。于是形状的还压不在行误差里——parab_euler 并联边末层 152.22 Pa、占自己 177.43 Pa 压降
#      的 85.8%，而行与块只差 -0.11 + 21.98 = 21.87 Pa：割线把还压吞进它谈判的摩擦半区，行从未看见这个形状。
#   5) 符号律与平台第一次是时间的函数——准稳 32 状态全在 (1.02, 1.0211) 且每臂每边跨 8 层散布 <=8.0e-5（一个数，
#      不是轨迹）；步进三臂层 1 在带之下（1.01427360..1.01750843 = 平台间隙的 68.4/71.0/68.6% X1、83.2% X2），而三臂
#      彼此差 <1e-4、强加 beta 相差 0.18：块回答 du/dt 时是「层」拥有发展，不是形状。末层落回与准稳孪生同一个平台
#      （1.02104568 对 1.02105374，不高于）。入口比平台满的恰好 32 个、差的 48 个，且与 pRisePa > 0 那批同一批。
#   6) 还压有时间签名——parab_euler X1 的 pRise/(rho*(beta_in-beta_out)*U^2) 从层 1 的 0.0704 爬到末层 0.9490，峰站随
#      区域发展下移 0.015 -> 0.035 -> 0.045 m；同一把尺子在准稳孪生读 0.9782..0.9928 且八层恒在 0.045 m（§5.36 那条
#      带），差的那一截是斜坡不是律。串联边从层 2 起坐在 1.1735..1.1800（散布 <0.01），层 1 的 0.4811 是族内唯一例外。
#   7) tau 被 du/dt 加长、而形状用同一个乘子——lo 0.4799 -> 0.5849 s（+21.9%）、parab 0.9786 -> 1.1683 s（+19.4%），
#      乘子步进 1.9973 / 准稳 2.0391，p17 只 +6.40%（与它离平台差 0.0071 的 beta 相称）；五臂尾段 |dq| 衰减都是几何
#      级数（0.65920..0.84267，层内散布 <5.6%）。教材式 (M_X1+M_By)/(R_X1+R_By) 是本例最干净的一格：对准稳两臂高估
#      +22.6%/+27.5%，对步进三臂到 -0.59%/-2.93%/+0.75%——被重复计的那笔惯量，正是它在准稳侧缺的那笔。
#   8) 分流——M_By/(M_X1+M_By) = 0.427275：准稳仍在附近（0.419889/0.464563，9% 内），步进系统性低于它
#      -18.92/-14.01/-18.40%（dt 0.2 s 短于这些臂自己的 tau，区域惯量不是第一脚唯一要推动的质量）。末层份额在两种
#      表述内部都保持形状次序（准稳 0.585000 > 0.407092；步进 0.537958 > 0.412185 > 0.400581），且每个步进臂都还在
#      自己准稳孪生之后（1.6 s 约三个 tau，不是斜坡尽头）。层 1 由表述移动 -17.5%/-20.9%，大于由形状移动
#      +10.6%/+6.1%：双重计费已经在掌舵，形状还没有。
#   9) 被钉死的边是干净实验室——串联速率 40 个（臂 x 层）状态 == 6.0e-5（<=1e-11）、层 2 起 |dq| <=1.36e-20，却照报
#      8.27/8.02/8.21 Pa 的 dpInertia（质量加权速度动过是因为剖面动了），它的动态惯量（q/0）因此只记录
#      5.5e13/6.4e16/1.4e14、不作判据。同一个开关在串联压降留 <=0.168%（层 3+），在并联速率留 -1.60%/-8.04%（5.0 倍）：
#      速率自由的地方双重计费不是可忽略项，它就是轨迹本身。X2 末层压降两种表述都保持抛物 < 1/7 < 塞状。
#  10) 时间表与成本——五臂 16/18/27/37/31 趟、SIMPLE 解次数恰为 2x 趟数、最坏残差 7.35e-4~9.58e-4 全在 1e-3 内、
#      管网线性代数每层到机器精度（newton <=9.1e-12）：迭代的是界面不是线性代数。du/dt 不免费，界面趟数 x1.69（lo）
#      /x2.06（parab），因为被谈判的那个系数现在随层一起漂；但步进不花块内扫数——平/1-7 的 Euler 并联边最差 17 扫。
#      300 扫花满与 blockConverged 0 仍是同一批（抛物串联边准稳 3/8、步进 1/8，其余 0/8），同一钉死速率塞状/1-7 7 扫。
#  11) 伙伴仍然合法、网络仍然闭合——y+ 21.0~32.6、nut/nu 46~134x、Re_Dh 4807~10000、每块进/出 <=5.2e-12；发展窗梯度
#      对模型自己的壁面剪力 2*rho*u_tau^2/H 在不带该层加速度时 -5.67%~+18.13%，把 rho*du/dt 记账后抛物两臂朝零走
#      （+10.53% -> +6.07%、+18.13% -> +9.63%，对应 +104.1/+169.1 Pa/m，因为它们并联边的速率到末层还在涨），塞状/1-7
#      的步进臂本来就 <4%——这项是真的、是剪力的几个百分点，既不用来救闭合也不需要藏。先导段 Churchill
#      11.264128 Pa 与 Src-J1、J3-Snk 到 1e-9，每层节点对逐边等于那一行动量式（含惯量项）。
# 另有一堵 dt 墙（记录，不作判据）：dt 0.05 下 X1 区域被要求比一层前更慢、于是真实地报出压强上升，撞上内核
# 「区域以正压降回答」的护栏；四种设法留下 -22.75 / -541.94 / -50.12 / -128.13 Pa 四个数，割线欠松到 0.4 的那档
# 活下来但 12 趟后停在 1.7035e-3 对 1e-3。本例把 dt 提到 0.2 s 让整案闭合，那几档仍未升格为判据。
# 写这类判据时的坑：平台族必须先按臂分家——准稳臂跨 8 层是「一个数」（散布 <=8e-5），步进臂层 1 是「一个没完成的
# 状态」（带下 68~83%），把两族塞进同一个区间断言会把层 1 判成违规；双重计费那种恒等式要断在**绝对 Pa** 而不是相对
# 差，因为塞状臂的摩擦半区差是 -0.11 Pa 这种量级；被钉死边的 inertanceDynamic 是 q/0，只能记录不能断言；app 打的
# transient step 行不都带 name/betaOut 字段，索引前先按字段存在过滤，否则 list.index 直接抛错。
# 十一族 490 条，五臂共 129 趟 / 258 次块解；单独跑整个 case 70.6 s
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentProfileUnsteady/case.dict
python3 tests/check_transient_coupled_multi_turbulent_profile_unsteady.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentProfileUnsteady/case.dict"

# 逐边 ddt 表（5.39）：上面那案的拓扑、两块 kOmegaSST 网格、壁面函数、操作点、dt 0.2 与 8 个时间层原封不动，
# 唯一新增的是表述可以逐边声明——couple{ ddtFrom { lo_x1free { X1 Euler; X2 steady; } … } }，
# 优先级 ddtFrom > ddts > details.<edge>.schemes{ ddt }，两个覆盖方向都走（三臂的逐臂词 Euler 被表把一条边收回
# steady；lo_steady 的逐臂词 steady 压过两块 dict 都写的 Euler）。七臂 = 塞状 {S/S, E/E, E/S, S/E} + 抛物
# {S/S, E/E, E/S}：后两格是新的且互为镜像；塑形族是为了让混合臂的每条边都有一个「同时固定该边表述与该边形状」
# 的对照，否则「这条边没变」只能拿去和一条平入口的边比。accelerates ( 0 0 0 0 0 0 0 )、tolerance 1e-3 /
# maxSubIterations 12、惯量归属固定在管网侧，所以臂间只剩「哪条边被表点了名」。
# 内核一行没改，改的是 app 声明面：src/coupling/apps/transientCoupledMultiApp.cpp:272（逐臂逐边解析，缺边或
# 非法词 cfFatal）、:361（三级覆盖）、:440（打印生效的逐边词而不是声明的那个）。
# 要问的问题是一根轴各自问不出的：chord 早就逐块取（coupledNetwork.cpp:1182），所以混合臂在一趟里谈判两种表述，
# 那么分流是否等于它各条边的和。
# 十族 202 条的要点：
#   1) 声明到达「边」——112 个（层 x 边）状态各带自己那条边的标志（56 步进 / 56 准稳，混合臂里两种坐在同一层），
#      块交回 dtLevel 0.2 s（步进）或 0（准稳）；每条边的行惯量是该边网格等效孔径的 rho*L/A（X1 8e6 / X2 7.5e6，
#      1e-14）且块自测同数，slab 也是它开场宣布的那一块（2.8e6 / 2.1875e6）。
#   2) 半区逐边切——dpFriction + dpInertia == dpBlock 到 2.2e-16；56 个准稳状态的惯性半区与 rowError 恰为 0、
#      56 个步进状态全带非零半区（含速率被钉死的）；混合臂每层恰一条边有惯性半区。行那笔 rho*L/A*dq/dt 在两条边上
#      都留在账上（112 状态到 0.0 Pa），不管该边的块是否自己付费。
#   3) 停机半边也逐边——按边重算 |R*q - chord|/chord（步进边取摩擦半区、准稳边取整压降）复现打印残差到 0.0e+00，
#      112 状态全 <1e-3；层 1 的 7 条步进边按摩擦读 0.0000~0.0137%、按整压降读 17.62~32.00%，准稳边两种读法同数。
#   4) 最强的一格——你留下的那条边就是它自己控制臂的算数：混合臂 X1 与「同表述同形状」的对照 15 字段 x 8 层逐位
#      0.0；X2 复现自己控制到 6.210e-4 Pa/1538.4（平）与 4.476e-1 Pa/719.8（塑）、速率到 6.8e-21；镜像臂的步进
#      钉死边距全 Euler 臂同边 8.6e-4 Pa 且逐层扫数相同——相等对该边自己的表述成立，不是对某个臂。尺度：同一条边
#      翻自己的表述动 662.0 / 2109.5 / 799.7 / 2148.9 Pa。
#   5) 两样不成加和——界面留在被留下那条边上的价移动 27.2 / 588.8 / 10138.7 Pa.s/m3（1.1e-6~8.5e-4 相对）而压降只
#      动 6.2e-4 相对：能相加的是边的物理不是谈判终点；塑形钉死块的预算差 290 扫（混合臂 6/8 花满 300 对控制 3/8），
#      那是 5.37 的容差地板透出来——速率被管网握着时塑形入口是「解到接近预算」的，够不够预算是邻居轨迹能翻的硬币。
#   6) 双重计费逐边出生也逐边停——层 1 X1 53.89%（平）/62.66%（塑）、X2 44.06%/56.01%，同臂被留准稳那 8 层恰 0.0
#      （不越过共享的节点对）；末层自由边仍 1.021%/12.328%，钉死边 0.000005~0.005102%；步进钉死边 49 状态速率恒
#      6.0e-5（|dq|<=1.36e-20）却照报 8.271/8.020 Pa 的 dpInertia。
#   7) 斜坡属于速率自由的那条边——混合臂的层 1 分流、末层去处、尾衰减比、时间表都与「自由边同表述同形状」的对照
#      逐位相同（0.419889/0.346431/0.367420；0.407092/0.400581/0.537958；0.64159/0.71416/0.84483），tau 也同数
#      （0.4506/0.5941/0.9897/1.1861 s）：du/dt 加在弛豫边上使斜坡 +31.8%（平）/+19.8%（塑），加在钉死边上加长 0；
#      混合臂内形状仍 x2.00。总趟数 16/27/27/16/18/37/37（178）→ x1.69（平，弛豫边）/x2.06（塑）/x1.00（钉死边），
#      SIMPLE 解恰 2x 趟数（356）。对 M_By/(M_X1+M_By) = 0.427275：准稳自由边 -1.73%/+8.73%，步进 -18.92/-14.01%。
#   8) 行被记的加速度仍形状盲——混合臂层 1 塑形 0.965755 / 塞状 0.962319（差 0.34pp，与 5.38 同），同一开关在同一条
#      边上把摩擦半区移 -14.8%（层1）/-51.2%（末层）；beta 声明值到 1.1e-15、截面积分到 6.7e-16。
#   9) 平台与还压也逐边——准稳边从层 1 就坐在平台（散布 3.0e-9~8.0e-5，一个数不是轨迹），步进边从带下
#      （1.0143/1.0173/1.0175）落到带上；还压状态恰 48/112（三个塑形臂两边八层）且与 beta_in > beta_out 同批；
#      同一个塑形臂两条边只差「自己的 du/dt 在不在产生压降的那条方程里」：准稳钉死边八层恒 333.4297 Pa @ 0.035 m
#      （动 1.2e-6），步进自由边是 4.6838 -> 152.2237 Pa 的曲线、峰站走 0.015 -> 0.035 -> 0.045 m；镜像
#      （parab_both 的钉死边一旦步进）峰走 [0.025, 0.035]、摆 200.9 Pa。分离性：混合 vs 控制在准稳钉死边差
#      1.0e-6 Pa，对照步进边上那 147.5 Pa 的曲线。
#  10) 墙第一次与失败不重合——准稳钉死边 7 扫（步进它的臂也 7 扫）、塑形自由边 [17,19,20,21,20,21,23,22] 两臂逐位
#      相同；墙属于形状（3/8、1/8、6/8 花满）而平臂与所有 X1 从不碰（其余十一列最差 28 扫），但 parab_x1free 那 6
#      层只有 5 层报 blockConverged 0；最坏残差 9.578e-4、newton <=9.09e-12。
#  11) 伙伴仍合法、网络仍闭合——y+ 23.0~32.6、nut/nu 46~134x、Re_Dh 4807~10000、每块进/出 <=3.8e-12；发展窗对
#      2*rho*u_tau^2/H 不带加速度时 -5.67%~+18.13%，记账后三种读法分别说出这笔账归谁：钉死边同一个数（速率停了就
#      没得加）、被留准稳的自由边仍移动（-5.51 -> -5.29%、+10.53 -> +6.07%，那项是行的、块的表述不拥有它）、步进的
#      自由边按形状分方向（平 -3.79 -> -5.41% 离零 +28.5 Pa/m；塑 +18.13 -> +9.63% 朝零 +169.1 Pa/m）。先导段
#      Churchill 11.264128 Pa 到 1e-9，每层节点对逐边等于那一行动量式（含惯量项）。
# 写这类判据时的坑：「逐位相同」要断成 == 0.0 而不是容差，而它的对照必须同时固定该边的表述与该边的形状——只固定
# 表述时，塑形混合边只能和平入口的臂比，任何差异都无法归给邻居；同一族的相对差要按被比较量自己的量级归一（钉死边
# 平入口差 6.2e-4 Pa/1538 Pa、塑形差 0.448 Pa/720 Pa，用同一个绝对界会把塑形那格判成违规、或用同一个相对界把平
# 那格放松到没有内容）；tau 的绝对值依赖估计窗口，本案取自由边层 3->8 的衰减比几何均值，不能与 5.38 的六层拟合
# 直接比大小，只能比同案内的两两之比；「花满预算」与「blockConverged 0」在 5.38 是同一批状态、在本案第一次不重合
# （6 对 5），所以这类一致性断言要按状态集合写、不能按计数写；case.dict 的注释行漏一个 # 会得到
# 「expected ';' near …」且不指行号，先查最近那一行。
# 十族 202 条，七臂共 178 趟 / 356 次块解；单独跑整个 case 129 s（ctest -j8 两次分别记 137.6 s 与 214.4 s，同批负载决定）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentEdgeDdt/case.dict
python3 tests/check_transient_coupled_multi_turbulent_edge_ddt.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentEdgeDdt/case.dict"

# 逐边 ddt × 逐边惯量归属同台（5.40）：上面那案的拓扑、两块 kOmegaSST 网格、probe、壁面函数、操作点、dt 0.2 与 8 个
# 时间层原封不动，唯一改动是把两根耦合管各自的图纸孔径故意写歪——X1 写 0.005（该截面等效直径 0.0056418958）、X2 写
# 0.006（等效 0.0063830765）。5.39 用的是精确等效孔径，于是 rho*L/A（图纸）与 rho*V/A^2（区域）本来就是同一个数，
# 归属那张表在任何结果里都看不见；写歪以后三本账成为三个数（X1 1.0186e7 / 8.0e6 / 2.8e6，X2 8.4883e6 / 7.5e6 /
# 2.1875e6 kg/m4，图纸把区域高估 1.273x 与 1.132x，slab 是区域的长度分数 0.350000 与 0.291667）。八臂 = 两边准稳
# {图纸, 区域} + 两边步进 {图纸, 区域, slab} + 5.39 的混合臂（es_nnn，处处图纸）+ 两个新格 es_mix / se_mix：同在
# X1 步进、X2 准稳，但步进那条边改吃 slab，两臂互为镜像——逐臂开关说不出「X1 是 Euler+slab 且 X2 是 steady+network」
# 这种话，所以这两张表必须能各自点名同一条边。accelerates 全 0、tolerance 1e-3，惯量归属在时间循环之前安装
# （coupledNetwork.cpp:1092），所以它对准稳边同样生效。内核一行没改，两张表都是既成机制（ddtFrom 解析在
# transientCoupledMultiApp.cpp:273、inertiaFrom 在 :303，非法词或缺边直接 cfFatal）。
# 要问的问题是内核那一行到底在数什么：rowDrop = R*q + M*dq/dt、rowError = (rowDrop - dpPa)/dpPa
# （coupledNetwork.cpp:1160/:1201，准稳边恰为 0，:1202），所以步进边上恒有
#   rowError * dpBlock = (R*q - dpFriction) + (M*dq/dt - dpInertia)
# 左边第一项是界面从不迭代的弦差（chord = unsteady ? dpFrictionPa : dpPa，:1182），第二项才是「这条边的惯量记在谁
# 账上」的账单——slab 只消掉后者。
# 十族 203 条的要点：
#   1) 两张表到达同一条边——128 个（层 x 边）状态各带自己那条边的两个词（72 步进 / 56 准稳，混合臂里两种坐在同一层），
#      es_mix 与 se_mix 的两份表在同一条边上互相不同意；装上行的 M 与该来源自己宣布的数到 1e-14 一致，块的回报同时带
#      inertance（区域）与 inertanceProbe（slab），所以接口装的是伙伴能核对的数。
#   2) 主判据：双重计费确实是那两项之和且逐边成立——72 个步进状态全部满足上面的恒等式，最差 6.4e-17（按该式四项绝对
#      值之和归一，es_mix 层 1 X1）；半区仍逐块切到 2.2e-16，56 个准稳状态的惯性半区与 rowError 恰为 0。层 1 两项不是
#      一个量级：X1 弦差 +0.0472 / +0.1015 / +0.6006 Pa 对账单 +718.87 / +551.40 / +11.79 Pa（15236x / 5430x / 20x），
#      X2 弦差 <=2.5e-4 Pa 对 +1903.86 / +1607.38 / +13.63 Pa；到平台顺序翻转（ss_both 层 8 X1：-0.3224 对 +0.0066）——
#      能被来源移掉的只有第二项。
#   3) 阶梯单调且第二项就是图纸多收的那笔——层 1 rowError X1 77.0836% -> 53.8941% -> 0.9364%、X2 52.1912% -> 44.0637%
#      -> 0.3737%，slab 消掉图纸账单的 98.36% / 99.28%；(M_图纸 - M_slab)*dq/dt 预测 708.19 / 1890.23 Pa 对实测 +718.87
#      / +1903.86，差额是块自己的 dpInertia。钉死边阶梯的绝对价更高（1903.9 对 718.9），因为它的 dq/dt 就是整步。
#   4) 归属只挪行侧那半——层 1 X2 的 dpInertia 三种归属下同为 642.6188 Pa（最差 1.1e-8 相对）；更强的形式是两次归属之差
#      恰为 (M_a - M_b)*dq/dt（296.4791 对 296.4791，4.3e-8 相对；1593.7496 对 1593.7500，2.7e-7；1890.2286 对 1890.2291，
#      2.4e-7），那点残差（1.3e-5 / 4.3e-4 / 4.5e-4 Pa）恰好被块自己两项之差吸收——所以 rowError 对装上的 M 线性，斜率
#      是该层的 dq/dt。反面：块的惯性半区是 M_slab*dq/dt 的 96.02 ~ 97.92%（9 个层 1 状态）——slab 是对的板、错的 du/dt。
#   5) 消掉之后剩下的正是界面本来在追的数——装 slab 的 24 个（层 3 起 x 步进）状态 |rowError| 最大 9.646e-4，低于停机
#      容差 1e-3，自由边还在同层残差 1.73x 内；图纸的自由边同层到 1.176e-1（容差 118x、自身残差 153.2x）——这项界面
#      既迭代不掉、容差也没要求它迭代。
#   6) 准稳边照付却一分不评——nn_steady 层 1 X1 装上 1165.3 Pa 的 M*dq/dt 对一块报 310.8 Pa 且无加速度的区域（评分会是
#      375%），rowError 打印 0.0；整个准稳族 56 状态都是几百到上千帕对恒零评分。
#   7) 两条局部性定律不是同一条——es_nnn/X1 == nn_both/X1、se_mix/X1 == nn_steady/X1 均 15 字段 x 8 层逐位 0.0（后者
#      邻居同时改了两张表）：表述跨边可交换。但 es_mix/X1 对 ss_both/X1（该边同表述同来源、只邻居不同）只 2/15 逐位，
#      速率差 7.41e-15（2.9e-10 相对）、压降差 4.54e-7 Pa：归属不可交换，因为装上的 M 是全局线性系统的一项。同一律在
#      准稳侧：三条准稳钉死边相对 nn_steady/X2 差 2.232e-4 ~ 5.900e-4 Pa（对 1538.4）。而评分为 0 不等于来源无所谓：
#      nn_steady 与 pp_steady 的自由边速率差 9.18%、压降差 50.5 Pa，rowError 是同一个 0.0；es_mix/X1 与 es_nnn/X1 同边
#      只差来源，速率差 6.735e-6 m3/s（28.3%）。钉死边压降沿阶梯不动（4.4e-8 相对）而 rowError 跨 51.8 个百分点。
#   8) 斜坡跟着自由边的两张表且单调——层 1 分流 slab 0.431868 > 区域 0.346431 > 图纸 0.319615（步进）、准稳侧 0.419889
#      > 0.381355；tau 0.2997 / 0.6014 / 0.6134 s（slab->图纸 2.0x），准稳 0.4329 / 0.4773 s 是管网自己的惯性松弛。三个
#      混合臂与「自由边同表述同来源」的对照逐位相同（0.319615/0.6134、0.431868/0.2997、0.381355/0.4773）——表没点名的
#      那条边把分流改变 0.0。去处仍聚（0.396201 ~ 0.407166）：归属在层 1 分开注入量的 11.2%、层 8 只剩 1.1%，平台散布
#      是容差不是来源选的物理。
#   9) 成本由自由边定且单调——趟数 16/16/26/27/30/26/30/16（187），SIMPLE 解恰 2x（374），最差扫数 15~17、无一碰 300
#      预算、blockConverged 全 1、最坏残差 1.861e-5 ~ 9.668e-4。准稳 16 -> 图纸 26 -> 区域 27 -> slab 30；se_mix（步进的
#      那条是钉死边）只花 16 趟，与完全准稳的臂同价。
#  10) 伙伴仍合法、网络仍闭合——y+ 22.8~32.6、nut/nu 45.2~74.9x、Re_Dh 5215~10938、beta_out 1.020947~1.021052、每块
#      进/出 <=2.6e-16，且八臂全部 flat 入口（下面任何差异都不是形状效应）；每层节点对逐边等于那一行动量式含惯量项
#      （1e-9）而不是块的压降，所以案件命名的那本账进了网络采纳的压力场；先导段 Churchill 11.264128 Pa 到 1e-9，层 1
#      行与区域最大 disagreement 2546 Pa、两项都有界。
# 写这类判据时的坑：恒等式两边是相消的量，归一标度必须取式子里四项的绝对值之和，不能取左右任一侧——pp_both 层 8 X2 的
# 两个括号相差三个数量级却几乎抵消（|lhs| ~1.5e-4 Pa、绝对误差 6e-14），按 max(|弦差|,|账单|) 归一会得到 3.99e-10 的
# 「假违规」；「线性于 M」这类定律要用括号差而不是 rowError 差来写，(rowError_a - rowError_b)*dpBlock_a 里混进了
# rowError_b*dpBlock 的变化，残差要由块自己那两项之差界定；「低于容差」不能对钉死边断言——ss_both 的钉死边残差本身是
# 5.66e-8，比值 1495x 毫无内容，只能对自由边断言比值、对全部状态断言绝对界；准稳边的评分是定义出来的 0（:1202），所以
# 任何「归属改不改结果」的准稳判据都必须读轨迹本身而不是读 rowError；本案的三本账要真的成为三个数，孔径必须先写歪，
# 5.39 那种精确等效孔径下归属是测不出来的。
# 十族 203 条，八臂共 187 趟 / 374 次块解；单独跑整个 case 约 28 s（ctest -j8 记 34.58 s）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentEdgeDdtInertia/case.dict
python3 tests/check_transient_coupled_multi_turbulent_edge_ddt_inertia.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentEdgeDdtInertia/case.dict"

# 塑形入口 × 逐边两表同台（5.41）：上面那案的拓扑、两块 kOmegaSST 网格、probe、壁面函数、操作点、故意写歪的两根 proxy
# 孔径（X1 0.005 对等效 0.0056419、X2 0.006 对 0.0063831；三本账因此是三个数：X1 1.018592e7 / 8.0e6 / 2.8e6，X2
# 8.488264e6 / 7.5e6 / 2.187500e6 kg/m4）、dt 0.2 与 8 个时间层原封不动，唯一加的东西是 starts.<arm>.profile——逐边的
# 入口形状列表（解析在 transientCoupledMultiApp.cpp:170-193，长度不等于耦合边数直接 cfFatal；每臂每条边生效的形状在 :433
# 打印；内核词表 flat/parabolic/power1over7 在 coupledNetwork.cpp:423-426、写入在 :491；ddtFrom :273、inertiaFrom :303）。
# 七臂 = 平口纯归属阶梯 {f_net, f_region, f_slab}（= 5.40 的 nn_both/pp_both/ss_both 换名，用来做跨例回归）+ 阶梯前两档
# 塑形 {p_net, p_region} + 两个混合臂 p_mix（X1 步进+图纸+塑形，X2 准稳+slab）与 p_pin（X1 准稳+图纸，X2 步进+slab）：
# p_pin 就是 5.39 记下的 parab_x2pin 那一格，p_mix/X2 就是 5.40 边界① 的准稳 x slab。accelerates 全 0，内核一行没改。
# 一堵墙写在 case.dict 的注释里：slab x 步进 x 自由边一旦同时遇上塑形入口，层 1 内发散（coupledNetwork.cpp:628
# "the detailed region reports a non-positive drop (-89.475517 Pa) at 0.000022 m3/s"），所以本案的 slab 只出现在准稳边
# 与钉死边上——不要为了"补全矩阵"把它放到塑形的自由步进边上。
# 要问的问题是 5.40 那两笔账单到底归谁定价：账单是 (M_装上 - M_块自己付)*dq/dt，在钉死边 dq 由拓扑给定、与轨迹无关，
# 所以如果形状能移动它，"记了谁的惯量"和"块做了什么"就不可分；如果它不动而块的 dpInertia/dpBlock 大动，归属阶梯就只是
# 所装 M 的性质，形状只给块那一半重新定价。
# 十一族 158 条的要点：
#   1) 两张表 + 形状逐边到达——112 个（层 x 边）状态各带自己那条边的三个词（96 步进 / 16 准稳），p_mix 与 p_pin 的两份
#      表在同一条边上互相不同意；装上的 M 与来源宣布的数逐字一致，块回报同时带 regionInertance 与 probeInertance。
#   2) 恒等式对形状免疫——rowError*dpBlock == (R*q - dpFriction) + (M*dq/dt - dpInertia) 在两种形状共 96 个步进状态全部
#      成立，最差 9.3e-17（按那四项绝对值之和归一，p_region 层 1 X1）；dpFriction + dpInertia == dpBlock 到 1.7e-16；
#      16 个准稳状态 rowError 恰为 0.0。
#   3) 主判据（钉死边）——X2 层 1 账单 1903.8603 Pa 平 对 1903.4177 Pa 塑（0.023%），而块的 dpBlock 掉 21.4%（沿 8 层最大
#      逐层相对变化 113.7%）、dpInertia 642.6188 -> 643.0614 Pa；形状能动弦差（2.47e-4 -> 2.38e-2 Pa，97x），但即使塑形
#      也只占账单的 1.3e-5；全案最大弦差 0.6423 Pa，最小的层 1 账单 11.79 Pa（f_slab/X1）仍是它的 18x。自由边换了主语：
#      账单跟着速率走（1.917688e-5 -> 2.014387e-5 m3/s，+5.04%；718.8665 -> 754.2383 Pa，+4.92%）。
#   4) 阶梯对形状盲（在它削掉的比例上）——图纸 -> 区域削掉 X1 层 1 行误差 30.08% 平 / 29.30% 塑（差 0.78 pp），X2 15.57% /
#      15.58%；平口三档逐字复现 5.40（X1 rowError 77.0836 / 53.8941 / 0.9364%，账单 718.8665 / 551.3981 / 11.7934 Pa；
#      X2 52.1912 / 44.0637 / 0.3737% 对 1903.8603 / 1607.3812 / 13.6312 Pa），塑形 region 档逐字复现 5.39 的 parab_both
#      （62.6559 / 56.0164%）——5.39 的 network 数就是本案的 partner 数，两案是同一臂的两种叫法。
#   5) 准稳边吃 slab 照付不评（5.40 边界① 填上）——p_mix/X2 付 656.25 Pa 的 M*dq/dt，块报 dpInertia 0.0、评分恒 0.0，
#      同 dq/dt 下图纸会收 2546.48 Pa，来源藏掉 74.2%；p_pin/X1 吃图纸，层 1 藏 1264.30 Pa、层 8 仍藏 35.06 Pa；准稳钉死
#      边的隐藏账单只在层 1 出现一次（656.25 然后 0.00），因为注入以后 dq 就是 0。块的回报逐边给 dtLevel：steady 那条恰
#      为 0 s、Euler 那条 0.2 s——形状仍然每层重写在这条看不见时间层的边上。
#   6) 斜率律只在拓扑钉住速率的地方精确——X2 上两档归属的账单差恰为 (M_a - M_b)*dq/dt（296.4791 对 296.4791 Pa，0.0000%；
#      1890.2276 对 1890.2291 Pa，0.0001%，0.0015 Pa 残差落在块自己那 0.0560 Pa 里）。反面：自由边同层两臂的 dq 本来就差
#      7.74%（平）/ 8.62%（塑），拿该式预测偏 22.6%——写这种判据时先看速率是不是给定的，否则"恒等式"只是巧合。归属也进
#      不到块里面：X2 的 dpInertia 沿 8 层在两种归属下不动到 8.1e-6 Pa（平）/ 1.8e-2 Pa（塑）。
#   7) 局部性多了一个触发——p_mix/X1 与 p_net/X1（自己三样全同、邻居两张表都改）在层 1..5 的 15 个字段逐字相同，只在层
#      6..8 分开（rowError 最差 3.1e-4、dpBlock 6.8e-3 Pa），而 p_net 的钉死块第一次填满 300 扫预算就在层 5：共享前缀
#      恰好在邻居的块放弃的那一层结束（这是实测巧合，不能写成因果）。反面：只改自己的表述（p_pin/X1 对 p_net/X1，M 逐字
#      相同）就把自由边速率动 4.681e-6 m3/s、层 1 分流 0.335731 -> 0.413742，而那本该报出 88.6217% 差额的评分在准稳边
#      恰是 0.0。
#   8) 墙是形状的、预算由归属定价——平口处处不碰 300（最差 17 扫），塑形只在钉死边碰（p_net/X2 层 5、p_mix/X2 层
#      1,2,4,5,6,8、p_pin/X2 层 4,7）；同一条塑形钉死边在图纸下要 300、在区域下峰值 162。这 9 个填满预算的状态与报
#      blockConverged 0 的 9 个完全重合（与 5.39 的 6 对 5 相反），9 个的界面残差最差 8.912e-4 仍在 1e-3 内、那一层照常
#      走完：一致性断言要按状态集合写，而且"重合"本身在这一案是新的观测，不是从 5.39 推出来的。
#   9) 成本由两个旋钮一起定，形状是贵的那个——趟数 26/27/30/35/37/36/17（208），SIMPLE 解恰 2x（416），每臂最差残差
#      7.279e-4 ~ 9.555e-4；同归属下塑形比平口贵（26 -> 35、27 -> 37）。
#  10) 斜坡由形状定尺度、由惯量定次序——tau 平口 0.6494 / 0.5932 / 0.4235 s（图纸/区域/slab），塑形 1.2824 / 1.2023 /
#      1.0907 s，同归属下形状把松弛时间翻倍（x1.97 / x2.03）；层 1 分流平口 0.319615 / 0.346431 / 0.431868、塑形
#      0.335731 / 0.367420 / 0.413742，终点平口 0.396201 / 0.400581 / 0.407166、塑形 0.518647 / 0.537958 / 0.518634 /
#      0.563791。
#  11) 伙伴仍合法、网络仍闭合——y+ 22.8~32.6、nut/nu 45.2~133.7x、Re_Dh 5365~11968、beta_out 1.020056~1.021052、每块
#      进/出 <=2.3e-12；beta_in 平口 1.0011 档、塑形 1.1162 / 1.1012（Poiseuille 的 1.18172 被壁面函数管道自己的发展段
#      松掉），所以任何一条臂间差异都能归到形状或归属其一。每层节点对逐边等于那一行的动量式含 M*dq/dt（1e-9），两段先导
#      管仍是管网自己的 Churchill 11.2641281367 Pa（1e-9），Snk 恒等于参考压力。
# 写这类判据时的坑：transient block 每臂每条边只有一行（块的最后回报），键是 name 不是层——把它当逐层表读会得到
# KeyError；"形状不动账单"只能对钉死边断言，自由边的 dq 随轨迹变，同一句话在那里必须换成"账单跟着速率走"，否则
# (M_a - M_b)*dq/dt 预测会偏 22.6% 而不是 0；min_bill 之类的界只能在层 1 的"步进"账单上取，全 8 层的最小值会抓到接近零
# 的平台尾（f_net 层 8 的 0.0066 Pa）；5.40 那种"墙与未收敛报告不重合"的结论不可移植——本案两者完全重合，所以判据要
# 写集合而不是写计数；app 头部那行 interface (...) 仍只描述第一臂，逐臂开关要读每臂自己的 transient start 行。
# 十一族 158 条，七臂共 208 趟 / 416 次块解；单独跑整个 case 143.37 s（ctest -j8 记 153.42 s，同批负载把最慢几条彼此推高）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentProfileInertia/case.dict
python3 tests/check_transient_coupled_multi_turbulent_profile_inertia.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentProfileInertia/case.dict"

# 对流动量账逐层报出（5.42）：**同一个 case、同一份 checker**，没有新增算例也没有新增 case 关键字。内核只做一件事——
# 在块的 profile 台账里用同一次求解的 betaIn / betaOut、该层速率与实测入口面积做一个减法
# rho*(betaOut-betaIn)*q|q|/A2（CoupledNetwork.h:314-:319，字段 :58 / :490），经 coupledNetwork.cpp:1226 交给边界，
# 打印成 transient step 行里 pRiseX 之后、unsteady 之前的 dpConv（transientCoupledMultiApp.cpp:527）。
# coupledNetwork.cpp:1200 的 rowDrop = R*q + M*dq/dt、弦的取法、rowError 的定义一处未改——5.41 那 158 条判据一条没改就
# 全绿，就是这件事的实测形式（第 12 族新增 11 条，共 169 条）。要看什么：
#   1) 恒等与符号——112 个状态 dpConv == rho*(betaOut-betaIn)*q|q|/A2 到 1.44e-15；48 个平口全正、64 个塑形全负，
#      没有例外。符号就是 betaOut-betaIn 的符号，所以这一族的第一条判据同时也是"形状"与"账"之间最硬的那道锁。
#   2) 它有多大——平口 f_net X1 7.632 -> 17.924 Pa、X2 57.615 -> 69.844 Pa（X2 那三个平口臂一致到 5.7e-8，X1 在层 1
#      是 7.632 / 9.157 / 15.044）；塑形 p_net X1 -55.239 -> -148.974、X2 -277.295 -> -283.272 Pa。占块自报压降：
#      平口最坏 5.73%（f_slab/X1 层 3 的 19.715 Pa 对 344.36 Pa）、p_net/X2 平台 39.3%、p_pin/X1 层 8 121.1%
#      （-179.391 Pa 对 148.15 Pa，那条是准稳边不是钉死边）。
#   3) 归属那张表进不来，但断言不能瞎选——钉死边三档平口的 dpConv/q2 八层一致到 5.7e-8（块的 dpInertia 也同为
#      642.6188 Pa）而那张表在该边装的惯量差 3.64x；自由边要断言的是"除掉 q2 之后"：dpConv 本身差 97.1%（层 1
#      7.632 / 9.157 / 15.044 Pa）而速率本身差 35.1%，系数差只有 7.96% 且逐层塌到 0.69%。
#   4) 形状一动它就是 4.06x（1.9401e10 -> 7.8687e10 Pa.s2/m6）——与 5.41"形状不动账单"正好互补：账单跟着 dq/dt，
#      这笔账跟着 beta 差与速率。
#   5) 与 pRise 同号同量级不等值——|dpConv|/pRise 在 56 个层 >=2、回升 >50 Pa 的塑形状态上是 0.8475~1.2989，X2 平台
#      （层 4..8 x 四个塑形臂）0.84760~0.85110；层 1 是 2.0786（回升才 133.40 Pa 而对账 -277.29 Pa）。平口处处
#      pRise 恰 0.0 而账不为 0（+69.844 Pa）。
#   6) 它不是 5.40 那笔双重计费——过滤条件必须是账单的量级而不是"账单非零"：dpConv/gap 在 |gap|>1 Pa 的集合上会炸到
#      +-3568%（slab 档把账单压到 0.005 Pa 而账仍 69.8 Pa）。按 |gap|>100 Pa 筛出 9 个层 1 状态，dpConv 至多占 17.3%，
#      且在塑形臂上符号与账单相反；另有 16 个步进状态账单已消失（<0.01 Pa）却仍各欠 69.84 Pa（平口 X2 平台 15 个）到
#      283.27 Pa（塑形 p_region/X2 层 7 那 1 个）——注意这个集合是按"步进"筛的，把准稳边算进来是 17 个。
#   7) 准稳边同样欠——p_mix/X2 每层 -283.670 ~ -283.783 Pa 而 unsteady 0、rowError 恰 0.0；p_pin/X1 -90.637 ->
#      -179.391 Pa。块的行不装惯性项的地方，这笔账照样在。
#   8) 行仍是那一行——96 个步进状态 rowDrop == R*q + M*dq/dt 到 1e-6；把 dpConv 加进 5.41 局部性那条的字段集（15 -> 16）
#      之后，p_mix/X1 对 p_net/X1 的层 1..5 仍逐位相同：新字段不引入新的非局部。
# 写这类判据时的坑：dpConv/q2 的散布与 (betaOut-betaIn) 的散布按构造恒等（它就是恒等式本身），所以"表不动系数"要断言
# 的是同一批状态跨臂的一致性与"系数散布 < 该层速率散布、且逐层塌"，不是那两个散布相等；pRise 那条必须把层 1 排除在带外
# 并单独报它的 2.0786，否则 min/max 带会被起始层拉宽到失去意义；"报出 != 收进"要落在判据里就是 rowDrop 那条 1e-6，
# 否则读者会以为行已经认识 beta；跨臂找"最坏那一格"时别写 `max((key, value) for ...)`——元组比较先看 key，它会返回字典序
# 最后的臂（本案是 f_slab/X2/层 8 的 4.54%）而不是最大比值（真值是 f_slab/X1 层 3 的 5.73%），必须写 `key=lambda kv: kv[1]`。
# 十二族 169 条（5.41 的 158 条未改），七臂 208 趟 / 416 次块解（96 步进 + 16 准稳状态）、每臂最差残差与 5.41 逐字相同；
# 单独跑整个 case 199.37 s（ctest -j8 记 284.17 s，全量 97 项 297.95 s。本机当时 load avg 7~11，而 5.41 记那对
# 143.37 / 153.42 s 时是 4.7，同批稳态 profile 映射也从 255.9 涨到 297.9 s——这两条墙钟要按负载打折读）

# 对流动量账收进 1D 那一行（5.43）：新算例 transientCoupledMultiTurbulentConvectiveRow，§5.42 的拓扑 / 两块湍流伙伴 / dt 0.2 /
# 8 层 / 1e-3 / GS / 三张逐边表全部照抄，唯一新增的声明面是第四张逐边表 convectiveToRow { <臂> { X1 charge; X2 report; } }，
# 四臂 c_off(report,report) / c_both(charge,charge) / c_pin(report,charge) / c_free(charge,report)。内核改的是**一对**：
# 行侧 transient.cpp:248-:252 加 C*|x|x 与 Jacobian 2C|x|、:257-:259 把 |C|x|x| 计入 Newton scale，系数由 setCoupledConvective
# （:429）装；界面侧 coupledNetwork.cpp:1232-:1235 令 chordPa = (unsteady ? dpFriction : dpPa) - convectiveRowPa、:1255 用它算
# residual、:1266 把它加回 rowDropPa、:1289/:1303 让 QN 与 GS 两处刷新都读 chordPa、:1317-:1318 每层把 C = convectiveRowPa/q2
# 装回网络。为什么要成对：界面对 R 迭代到 R*q == chord，只加行侧不改弦，就是 §5.40 那笔账单的 beta 版本（同一笔钱收两次）。
# 要看什么：
#   1) 惰性对照必须逐位——c_off 对 5.41 的 p_net：16 字段 x 16 状态差 0.0，且每层趟数 [3,6,5,5,4,4,4,4] 一字不差。字段全对但
#      趟数变了，说明开关动了轨迹而不是只动了记账。
#   2) 不收第三笔的判据是 rowDrop - (dpNetwork + M*dq/dt + rowBeta) == 0.0（64 状态，精确 0），不是某个 1e-6 的容差。
#   3) 弦真的让出来了：内核自己的 residual 字段逐字等于 |dpNetwork - (dpFriction - rowBeta)| / chord。把这条写成判据，
#      "重复计费"就变成一个可失败的断言而不是推理。
#   4) 改名的价格是算得出的：R_charged/R_off == 1 - dpConv/dpFriction_off（最差 1.722e-3）。X1 涨 9.5~90.5%、X2 涨 12.5~39.3%，
#      比例就是形状账占摩擦半区的比例——图纸那个系数不再是壁的数。
#   5) 钉死边 vs 自由边：前者 q 逐字仍是 6.0e-5、除 R/dpNetwork/rowBeta 外最大移动 3.521e-2 Pa（那是界面停机噪声）；后者才是
#      真轨迹（q 相对 1.423e-4、行总额 0.125 Pa 对搬走的 148.9 Pa）。By 自己吃掉剩余（层 8 差 4.427e-9）而 X1+By 恒等于注入。
#   6) 局部性按字段集断言：c_pin/X1 对 c_off/X1、c_both/X1 对 c_free/X1 都是 17/17 逐位 0.0——邻居的声明改不了我这行。
#   7) 新的第一堵墙在代码里：chord = dpFriction - dpConv 非正即 fatal（coupledNetwork.cpp:1235）。5.42 记下过 121.1% 那一格
#      （准稳边的弦是割线），本案最坏 90.5% 所以没触发——写文档时别把它说成"已验证的行为"。
# 写这类判据时的坑：①setCoupledConvective 必须和 setCoupledResistance 一样**只在不拢的那几趟**、且在逐边循环之后调用——
#   pendingPipeFlow() 一旦 pending_ 为假就 fatal，早一步整案崩在第一个时间层。②这项**带符号**，校验只能写 isfinite，
#   照抄阻力的 >0 会把所有塑形臂当场拒掉。③chordPa 要同时喂给 residual、QN 的 r[] 与 GS 的松弛三处；只改残差会让不动点
#   和判据分家，表面上"收敛"了但双侧相加已经不是摩擦半区。④跨臂比较一个有符号系数时 `max/min - 1` 在全体为负时给出**负数**
#   （比值反了），必须先取 abs——本案 parabolic 入口下 dpConv 恒负，第一次就打印出 -1.1e-06 而不是 2.4e-05。⑤逐层趟数要跟
#   出厂序列比顺序，别在打印时用 sorted()：[3,4,4,4,4,5,5,6] 和 [3,6,5,5,4,4,4,4] 是同一个多重集，排序会把"惰性"读成一样而
#   把真实差异藏掉。⑥脚本里别写 `re = None` 之类局部名，它会遮蔽模块 re 让后面每一次 value() 都炸。
# 十族 35 条；四臂 35/34/35/34 趟、70/68 次块解、最坏单级 7 趟（上限 12）；单独跑整个 case 154.31 s user / 155.6 s 墙钟
# （ctest -j8 记 153.92 s，全量 98 项 370.98 s。这批的负载比 5.41/5.42 那两批轻，所以同一条用例的墙钟从 284.17 s 回到
#  233.81 s——四臂的本案按同批比例读作 153.92 s，判据只跑一次，不与那两批的七臂互换）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentConvectiveRow/case.dict
python3 tests/check_transient_coupled_multi_turbulent_convective_row.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentConvectiveRow/case.dict"

# 形状账 × 归属账同台（5.44）：新算例 transientCoupledMultiTurbulentShapeAndInertiaSplit，**一行内核代码也没改**——
# 只是把 §5.40 逐边化的 inertiaFrom（network|partner|slab）与 §5.43 的 convectiveToRow（report|charge）同时指向
# **同一条边 X2**，做成 3 档 × 2 位置 = 六臂 n_rep / n_chr / g_rep / g_chr / s_rep / s_chr；X1 在六臂里恒"图纸 + 报出"，
# inertias ( network × 6 ) 保持全档图纸、让逐边表成为归属轴的全部陈述。n_rep ≡ 5.43 的 c_off ≡ 5.41 的 p_net，
# n_chr ≡ 5.43 的 c_pin——两颗已出厂的对照格正好是这张因式表的左上 2x2。
# 要看什么：
#   1) 两根杆子还是根一根：三档的层 1 账单在 report / charge 两臂差 **0.0 Pa**，而同一档内 M 逐层逐位相同——
#      形状的署名与惯量的署名是两个地址。反方向同样成立：dpConv 在六臂层 1..7 逐位相同，而它旁边的 M 差 3.8803×。
#   2) §5.40 那两只括号第一次同时收缩：1903.44 → 13.21 Pa（144.1×）、rowError 0.6635 → 0.0046；同档内装上与否只差
#      1.4e-7。这一步把 5.43 边界⑥那句"账单那半仍未消掉"从推理变成读数。
#   3) 阶梯的第一级必须**对上厂字**：n_rep/X2 层 1 的 1903.4177 Pa 由 §5.41 发表的四个字段重算（rowError×dpBlock +
#      dpFriction − dpNetwork）而不是硬写一个数——文档只允许引用测得的值，也允许引用由前一篇发表的字段算出的值。
#   4) 归属不会把弦谈松：六臂最坏相对弦差同一个 7.9585e-4。换档买到的是行里那个隐式项更小，不是界面更宽。
#   5) 两个声明都只到自己那条边：其余五臂的 X1 对 n_rep/X1 在 17 字段 × 8 层逐位 0.0，X1 账单恒 754.2383 Pa。
#   6) 仍不收第三笔：96 状态上 rowDrop − (dpNetwork + M·dq/dt + rowBeta) == 0.0，恒等式最差 6.821e-13。
# 写这类判据时的坑：①落位不是省——把 slab 同时落到两条边的第一版在 s_rep 臂 fatal 于
#   `non-positive drop (-89.475517 Pa) at 0.000022 m3/s`（148.47 s，与 5.41 记下的那两个数逐位相同）：slab × 步进 ×
#   **自由**边 × 塑形入口是那堵墙，本案只能落在钉死边上。②归属只在 dq ≠ 0 处可见，而钉死边除层 1 外 dq 恒 0——
#   所以"两只括号同时收缩"在本案是一条**首层**陈述，别把它写成全时间轴的性质。③逐边 inertiaFrom 会盖过位置式
#   inertias 列表，两者同时给出时以表为准；本案靠这个把 X1 恒钉在图纸上。④M 的闭合式按边不同（X1 是 d² 圆管、
#   X2 是 4×8 mm 矩形 + 0.07 m 探针），只对一条边写常量再去断言另一条会得出 10185916.3578813 对 8488263.63；
#   写成 M_of(edge, source)。⑤按 rung 过滤再做 max() 会遇到空生成器（X1 没有 partner/slab 臂），要么先补 0.0
#   要么改按 (arm, edge) 求最坏。⑥迭代一份判据不必重跑 250 s：importlib 载入 checker、把它的 subprocess.run
#   换成预先捕获的日志，所有断言在 1 s 内跑完——只在最后用真实运行确认一次。
# 十族 34 条；六臂各 35 趟界面 / 70 次块解（合计 420）、逐层 [3,6,5,5,4,4,4,4]、最坏单级 6 趟（上限 12）——
# 两张表都是写进行里的系数而不是迭代，所以交叉不花钱。单独跑整个 case（六臂）247.78 s user / 249.85 s 墙钟
# （ctest -j8 记 244.43 s（全量 99 项 524.88 s，load ~8）对 520.87 s（全量 100 项 525.22 s，load ~19）——本案是六臂，
#  同批耦合用例彼此推高，墙钟按负载读；单跑 249.85 s 是这份判据自己的地板）
./build/src/coupling/apps/chflow-coupled-1d3d-transient-multi tests/cases/transientCoupledMultiTurbulentShapeAndInertiaSplit/case.dict
python3 tests/check_transient_coupled_multi_turbulent_shape_and_inertia_split.py \
    "$PWD/build/src/coupling/apps/chflow-coupled-1d3d-transient-multi" \
    "$PWD/tests/cases/transientCoupledMultiTurbulentShapeAndInertiaSplit/case.dict"

# 非牛顿管流（5.45）：新算例 steadyPowerLawRheology。内核只做**一次替换**——src/1d/src/network.cpp:101-:117 的 generalizedReynolds
# （Metzner-Reed：rho u^(2-n) d^n / (K' 8^(n-1))，K' = K((3n+1)/(4n))^n，指数是 n，5.47 在 tau_0->0 极限上量出并改正）取代原先两处
# rho u d / mu：求解用的 pipeLossAt（:257）与打印用的 generalizedReynolds（调用点 equipmentFlow，:369）。
# 声明面 pipes.<name>{ flowIndex n; consistency K; }（Network.h:73-:74，0 = 未声明）。
# 为什么判据敢写成"精确"：Churchill 的层流半区是 f = 64/Re，换上 Re_MR 之后它就是幂律圆管的**解析解**
#   dp = (L/d) K (3n+1)/n (8u/d)^n   ->   Q = (pi d^3/32) [n dp d / ((3n+1) K L)]^(1/n)
# 而不是任何拟合；n=1、K=mu 时 K'=mu（4/4）、8^0=1、u^1=u，**逐项**塌回原式，所以既有 98 条测试结构上不可能被这一步移动。
# 要看什么：
#   1) declared-but-Newtonian 必须逐位——PowN1 对 Newt 五字段 × 两个文件全部 ==（不是 <1e-12，是相等）。这一步替错了代数，
#      这里立刻显形，而且它同时是新算例对旧算例的无害性证明。
#   2) f * Re == 64 在 18 个状态上最差 1.563e-10：这条把"闭式确实是内核用的那条律"从推理变成读数，别省掉。
#   3) 三个指数的方向各管一件事：孔径（16 对 32）、长度（÷2 对 ÷4，指数不动）、驱动差（2^(1/n)）。剪切稀化的**尺寸律**是
#      另一个幂——这是本案对工程最有用的那句话，写成比值而不是写成绝对数，这样它不依赖任何一次求解的巧合。
#   4) tau_w = dp d/4L 对材料无关（五根同几何的管逐位 2.5 Pa），而 K'(8u/d)^n 各自复现自己的 tau_w（5.413e-12，改正
#      常数后的实测；改前那个数是 8.040e-13，它同时说明旧常数在 n!=1 的管上并不是解自己含的那条律）。
#   5) 两条拒收要**被断言**（非零 rc + 消息），不要只在文档里描述：flowIndex 缺 consistency（量纲不同，默认即静默错数）、
#      耦合元带流变（它的损失就是伙伴的陈述——5.43 那条不重复计费）。
#   6) Re_MR 全程 <= 170：这是"敢用闭式"的前提，写成判据（worst_re < 1000），否则层流断言和湍流半区会被混着读。
# 写这类判据时的坑：①Churchill 返回的是 **Darcy** 因子（层流 64/Re），Fanning 文献里的 16/Re_MR 直接抄进来会把压降差 4 倍——
#   本案第一次推导就栽在这上面（闭式对不上，直到发现 n=1 时应退回 32 mu L u/d^2 而不是 8 mu L u/d^2）。②九根管的两端必须都
#   接**定压节点**，否则牛顿解把并联分配夹在断言和物理之间，闭式就不再加和。③打印精度是 12 位有效数字，所以 v == Q/A 的
#   残差地板是 1.8e-12 而不是 1e-16：判据要按"打印了什么"来定，不是按"算了什么"。④`worst_re` 的上界要按**翻倍**那份文件读
#   （Newtonian 的 Re 随 dp 线性走，170 出现在宽孔径管），只看 base 会写下一个假的安全界。⑤Python 里别把 max(...) 跨行写进
#   单个 f-string 表达式——语法错误，先在外面算好再插值。
# 八族 16 条；四个 case 文件、零求解成本（单次稳态毫秒级），整条 checker 墙钟 0.14 s（ctest -j8 记）
./build/src/1d/apps/chflow-1d-steady tests/cases/steadyPowerLawRheology/case.dict
python3 tests/check_steady_power_law_rheology.py \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/steadyPowerLawRheology/case.dict"

# 屈服应力（5.46）：新算例 steadyBinghamYield，五个 case 文件。内核改的是**一对**，而第一次容易只改前一半：
# src/1d/src/network.cpp:315-:327 在 Illinois 括号**之前**显式判掉死分支（`yieldStress > 0` 才生效，牛顿/幂律那两支
# 的 g(0) 本来就是送它们来这边的符号，一行未动），:131-:160 把 Bingham 层流解反解成网络要的 loss(q)——
#   F(tw) = tw - (4/3)tau_0 + tau_0^4/(3 tw^3) = B = 32 mu Q/(pi d^3)，括号 [tau_0, B + (4/3)tau_0] 是解析的
# （F(tau_0) = 0、F' = 1-(tau0/tw)^4 > 0、右端 F >= B），Newton + 二分落到舍入位；使用点 :252-:256。**同时** :272-:282
# 把报出的摩阻系数从 Churchill 换成解自己含的 f = 8 tw/(rho v^2)（调用点 :370），声明面 pipes.<name>{ yieldStress t0; }
# （Network.h:99，0 = 未声明），解析与拒收在 :1107-:1117。
# 为什么必须成对：留着 Churchill 的 f，日志同一行的系数就再也复现不了那一行的 dp，"报的是哪条律"重新变成推理。
# 要看什么：
#   1) 死分支报的是 **0** 不是"很小"：Stk 的 Q/v/Re/f 四个字段打印就是 0，而 dp 仍是**被问到的** 1.0e4 Pa 而不是阈值
#      1.04e4 Pa——那是静摩擦吃掉要它吃的那笔，不是被封顶的动摩擦。这条只能写成 == 0.0，写成 <1e-12 就读不到东西。
#   2) f * Re * g(xi) == 64（最差 1.843e-10，地板是 12 位打印）是"系数就是解里的系数"的可失败形式；配套再断言
#      dp == f (L/d) rho v^2/2（2.904e-12），否则同一行可以填两个各自成立但互相矛盾的数。
#   3) Q / Q_Hagen-Poiseuille **就是** g(xi)：0.7339 / 0.4752 / 0.2432 / 0.0699 / 0.0031 且全在 (0,1]——比值形式不依赖
#      任何一次求解的巧合，也比绝对数更容易看出"阈值只删流量、不弯抛物线"。
#   4) 壁应力八根逐位 2.500 Pa（tau_0 不在 dp d/4L 里），而 8 tw/(rho v^2) 复现报出的 f 到 2.904e-12。
#   5) 死分支要放上**网络 Newton** 再看：Sr1a/Sr1b 两条全同屈服管把自由中点落在恰好 2.5e4 Pa（iterations=1、残差恰 0）
#      且两根速率逐位相同；Sr2a 卡住时整条线清空、D 被下游定压拽到恰好 2.0e4。这一步不需要新 Jacobian，因为闭式给出
#      dQ/d(dp) = k[1-(dpStatic/dp)^4] <= k，solvePipeFlow 里那个 Hagen-Poiseuille 上限对屈服管仍是**上界**。
#   6) 反问题同台：case_driven 把 Q 钉成 5.0e-6 让 dp 成未知数，解出的 7692.34 Pa 对 checker 自己那套**二分**（不是
#      Newton）反演差 5.326e-13；反向那根逐位是正向的负数——死分支判的是被选中方向的符号，不是正向运气好。
#   7) 拒收要被断言（非零 rc + 消息）：屈服+耦合元（没有墙留给阈值 hold）、屈服+PVT/黑油（那是常物性单相层流律）、
#      负 tau_0（g(xi) > 1，这条律唯一能给出比牛顿更小压降的方式）。第四条墙——屈服+幂律＝Herschel-Bulkley——当时也
#      是拒收，5.47 落的就是它：这个文件已从 case_rejectRheology.dict 改名 case_rheologyNowSolves.dict，checker 现在断言
#      它**解得出**（rc=0、Q>0、dp 仍是被问到的 1.0e4），数值本身归 5.47 的 checker。
# 写这类判据时的坑：①死分支判据要写成 `forward ? fLo >= 0 : fLo <= 0` 这种**跟随已选方向**的形式，只判 dpTotal 与
#   阈值的大小会在反向流那根管上落错支。②哨兵是**代码级**的：只有 `yieldStress > 0` 才进括号，未声明与声明 0 那两根跑的还是原来那一行，所以 Newt/T0zero 五字段
#   逐位 == 是无害性证明而不是巧合；代数上 `tau_0 -> 0` 也确实让 F(tw) = tw（Newton 的导数在括号左端点为 0，那里只能二分，
#   代码靠 `next` 的括号包含判据避开），所以两支不会在最后一位上分岔。
#   ③内层解松紧会伪装：xi=0.96 那根管 Q ∝ (dp-dpStatic)² 是软着陆，括号左端扫过去也"差不多对"，所以要单独断言
#   它**仍带流量**（6.11731e-08）而不是只看闭式残差。④别把阈值那一格写成 dp == dpStatic 的精确相等：0.02 不是二进制
#   精确值，4*2.5*20/0.02 与 1.0e4 谁大由舍入决定，本案把 tau_0 放两侧各留 400 Pa。⑤Q=0 时 f 必须是 0（没有速率可归一），
#   照抄 Newtonian 那支会在 Re=0 上让 churchillFriction 的 1e-12 地板吐出 6.4e13（既有行为，本案没改它）。
#   ⑥checker 的 parse() 要按文件收名字：case_driven 的管名不在主 PIPES 表里，按 PIPES 过滤会静默收到 0 条然后
#   在 KeyError 上炸成一个假的成功。
# 九族 26 条；五个 case 文件、零求解成本（case.dict 12 根管 1 步收敛、case_driven 6 步），整条 checker 墙钟 0.28 s（ctest -j8 记）
./build/src/1d/apps/chflow-1d-steady tests/cases/steadyBinghamYield/case.dict
python3 tests/check_steady_bingham_yield.py \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/steadyBinghamYield/case.dict"

# Herschel-Bulkley（5.47）：新算例 steadyHerschelBulkley，三个 case 文件。这一格是 5.46 在解析期拒收的那一支——
# 阈值 + 指数同时声明。内核新增的是**一条律的两个方向**，不是一个分支：正定律 src/1d/src/network.cpp:173-:185 是剖面
# 积分的逐项闭式（tw = dp d/4L、D = tw - tau_0）
#   Q = pi R^3 tw^-3 K^(-1/n) n D^(1+1/n) [ tau0^2/(n+1) + 2 tau0 D/(2n+1) + D^2/(3n+1) ]
# 反解 :196-:221 用单调二分（左端 = tau_0 即阈值本身，右端从 tau_0=0 的幂律应力起，不够就加倍，速率按 tw^(1/n) 增长），
# 两支路的选择器 yieldedWallStress 在 :227-:232（有指数走 HB，只有阈值走 5.46 的 Buckingham 三次式 Newton）。
# **同时** generalizedReynolds（:101-:117，HB 那支 :104-:109）与 5.45 的层流常数一起改正：K' = K((3n+1)/(4n))^n，
# 指数是 n 而不是教科书塌缩约定的 1——差 (4n/(3n+1))^(1-1/n)，n=1 处恰好为 1（所以 PowN1 == Newt 那条逐位对照看不见它），
# n=0.5 处是 1.25。这一格是 HB 在 tau_0 -> 0 极限上量出来的，改正后 5.45 的 16 条仍全绿。
# 要看什么：
#   1) 两条极限必须在**打印位**上会合，而且要写成 ==：n=1 那根管（HbN1）与不写指数那根管（Bing）三个字段逐位相等——
#      三次式 Newton 与剖面积分二分是两条代码路径，会合才证明新支是 5.46 的**推广**而不是对手；tau_0=0 那根管对幂律
#      闭式 7.023e-13，这条是改常数之后才成立的，改常数前先跑它会红。
#   2) 同一个积分**别逐项**再算一遍：把 int tau^2 gdot(tau) d tau 数值积分（s^(1/n) 在屈服半径处导数无穷，[0,1] 必须按
#      2 的幂向 0 加密，本案用 21 个面板 × 12 节点高斯-勒让德）对内核仍是 1.511e-09；单个不加密面板差 3.296e-06，
#      加密一档就是 0.000e+00——这条和闭式是同一陈述的两次独立推导，(3n+1) 写错一处不会在另一处得到同样的错。
#   3) 软着陆那两格（HbThin2 阈值上方 0.5 Pa、HbNear 差 0.1 Pa）落在 1.511e-09 / 1.078e-10，比陡状态的 2.331e-12 差
#      三个量级。**不要**把这当成律放宽：把打印的速率喂回闭式反演，内核留下的压差是 1.112e-06 / 1.822e-08 Pa，乘上各自
#      的斜率 d ln Q/d ln tw = 13.593 / 59.151 恰好复现观测（第 4 位）。equipmentFlow 的速率括号门是 :364 的 1e-9 相对，
#      它先于 :347 的残差门触发；压差是算例钉住的，余量要乘着陆自己的斜率才变成速率差。判据因此分组（陡 1e-11 / 软 1e-8）
#      并且把余量本身写成一条断言（预测/观测 < 5e-3），这样下一次谁改了那两个门会在这里显形。
#   4) 阈值的几何性与 n 无关（dp_static = 4 tau_0 L/d = 2000 / 8000 / 10400 Pa），D 的主导幂次有关（1+1/n：Bingham 平方、
#      n=0.5 立方）；卡住那根在 1e4 下四个字段打印就是 0 而 dp 仍是被问到的 1e4，同一根管在 2e4 下流出 6.04532e-07
#      （闭式差 9.720e-14）——5.46 的"死分支由被问到的压差决定"在 n!=1 上重演。
#   5) 串接那一对要放在**自由节点**上：C 恰好 2.5e4（翻倍 3.0e4），两根 halves 五字段逐位相同；驱动差翻倍时它们从
#      1.17871e-08 涨到 3.58889e-07，也就是每半边的 D 从阈值上方 0.25 Pa 走到 1.25 Pa——软着陆在两个点上被读，不是在一个
#      点上被断言。
#   6) 带指数的屈服管没有材料黏度可报，Re 换成表观数 8 rho v^2/tw，于是 f * Re == 64 是**定义**（2.911e-10 = 打印精度）：
#      断言它的唯一理由是"两个数不同步说明日志在印另一个状态"，别把它写成物理证据。有内容的是 HbN1/Bing 之比就是
#      g(xi) = 0.475200（黏度之比即屈服挖掉的速率份额）和 dp == f(L/d)rho v^2/2 差 1.029e-10。
#   7) 拒收的账要一起更新：删掉一条（指数 + 阈值，现在由 steadyBinghamYield 的 case_rheologyNowSolves.dict 断言它解得出），
#      新写一条（声明指数与阈值而不给 consistency，K 的单位是 Pa.s^n，不能被 mu 静默默认）；其余五条照旧。
# 写这类判据时的坑：①高斯-勒让德的点权最容易写成 P_{m+1} 的根配 m 个权重（off-by-one），本案第一次就栽在这上面；核对办法是
#   验 sum(w) == 2 且节点对上公开表，别只看"跑通了"。②数值积分那个尺度因子 tw^-3 与 span^(1+1/n) 漏一个就会得到 1.25 / 7 / 81
#   这种荒谬比值——先在同一根管上把闭式和积分对齐到 1e-15，再去比内核。③tau_0=0 的管要判的是 5.45 的幂律闭式，如果 checker
#   自己还留着旧常数，改正内核后这条会红，而且红的是 checker。④parse() 要按 key/value 收字段而不是按下标：dp 之后新印的字段
#   不该把后面的数挪进前面的槽。⑤翻倍那份文件的管名要和主表同一套，否则 series 那一对被拿到并联的压差上比，误差直接是 1。
# 十族 22 条；三个 case 文件、零求解成本（毫秒级稳态求解），整条 checker 墙钟 0.21 s（ctest -j8 记）
./build/src/1d/apps/chflow-1d-steady tests/cases/steadyHerschelBulkley/case.dict
python3 tests/check_steady_herschel_bulkley.py \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/steadyHerschelBulkley/case.dict"

# 非牛顿壁面律进瞬态那一行（5.48）：新算例 transientPowerLawRheology，五个 case 文件（2 个求解 + 3 个拒收）。
# 这一格不是"再写一遍律"，而是**把稳态那条验过三次的律换成唯一的一份**：新增 pipeFrictionPa
# （src/1d/src/network.cpp:390，声明 Network.h:202，必须在 cf 作用域——瞬态内核是另一个翻译单元，放匿名作用域
# 就只有 pipeLossAt 看得见），稳态 pipeLossAt（:254）在它之上只加阀项，瞬态 frictionDp（src/1d/src/transient.cpp:72-:80）
# 除了调用与一次奇延拓什么都不算。之前那一行自己写 rho u d / mu，所以声明 flowIndex 0.5 **不报错也不警告**，只是
# 被忽略然后报 Newtonian 速率——这是最坏的一种静默数字，也是本格唯一真正的改动。
# 阻尼同步换：frictionSlope（:85-:105）在声明指数时取律自己的中心差分并在 :106-:114 早退，Newtonian 下限 0.2*128 mu L/(pi d^4)
# 只留在未声明的路径（:115-:116）。回读接口 pipeFlowIndex/pipeConsistencyOf（:556/:563）让 app 只在声明了的管上打印
# wall-law 那行（src/1d/apps/transientApp.cpp:42-:54），于是"未声明的案子输出逐字不变"是能在字符串上失败的断言。
# 要看什么：
#   1) 跨内核同一条律要打在**打印位**：t=4 时瞬态的 Q 与稳态 app 同一声明的 Q，Newt/Decl1 1.96349540849e-05、
#      Thin05 与 Thin05Rev 3.92699081699e-06、Thin07 9.75728733974e-07 是同一个字符串；只有 Thick13 差在第 10 位
#      （稳态 3.41904523976e-05 对瞬态 3.41904523728e-05，7.26e-10）。别把这一条写成"逐位相同"除非真的逐位相同。
#   2) 那个 7.26e-10 要**归因**而不是放宽：把闭式的速率喂回内核自己的 pipeFrictionPa，被问的 1.0e4 Pa 回到 3.638e-16，
#      说明层流半区与剖面积分是同一句话，余下的是每步 Newton 退出（transient.cpp:418 的 maxRel<1e-11，残差按
#      rho L/(A dt) 定标）在弛豫里累加。判据里因此有**两条**：终端差 <1e-8，根差 <1e-12。谁改了那个门，第二条会红。
#   3) 整条轨迹都要重解，不要只判终点：checker 自己写隐式欧拉递推逐步对照 6 根管 × 2001 状态，最差 7.246e-10。
#      这是"行读的是哪条律"与"某个形状碰巧落在对的终点"的区分。
#   4) 阻尼斜率从尾部拟合 s = (I/dt)(1/r-1)（r 是窗口两端 |delta|/Q_ss 的几何比），对律的 n dp/Q 差 9.6e-09~1.6e-03，
#      而弦 dp/Q 差 49.9%（n=0.5）/29.9%（0.7）/23.0%（1.3）、Newtonian 下限 1.018592e+08 差 >10%——三个候选一起写才
#      说明"用的不是另外两个"。注意窗口：n=0.7 的 tau 只有 0.0075 s，dt 必须显著小于它才有得拟合（本案 dt 0.002）。
#   5) 惰性对照打在字符串上：n=1, K=mu 那根管与未声明那根管 2001/2001 状态**打印逐字相同**。
#   6) 拒收即边界：负指数、指数无 consistency、耦合边带指数（三条 rc!=0 + 报文 needle，判据在 transient.cpp:212-:217）。
#      第四条原本也是拒收——"瞬态只有动壁律，阈值需要一个自己在时间里积分的静止支"（稳态 5.46 那次便宜是因为问到的
#      压差已知，这里 dp 正是行的未知数，approximate 是禁止的所以只能拒）；5.49 落了那条支之后它改名
#      case_yieldNowSolves.dict 并被断言**解得出**（Bing 1.4409438297e-05 对闭式 1.4409438304e-05 即 5.528e-10、
#      Hb 2.2250853568e-06 即 1.138e-12、Stuck 打印 Q 0），所以"边界"这一格从 needle 换成了数值断言。
# 写这类判据时的坑：①函数放错作用域会得到 "call to 'pipeFrictionPa' is ambiguous"（定义在匿名 namespace 里，声明在
#   Network.h）；放 cf 作用域后 pipeLossAt 仍能经外层作用域看见它。②case 注释里引闭式要引**改过常数之后**的那条
#   （5.47 的剖面积分），5.45 的 Metzner-Reed 常数在 n=0.5 处差 25%，先在稳态 app 上对一下再写。③面积算错一格
#   （pi d^2/4 写成 d^2/4 之类）会让手推的时间常数差 4 倍而对不上打印的 v——核对 v == Q/A，别核对自己的算术。
# 十族 45 条；五个 case 文件，checker 墙钟 6.64 s（ctest -j8 记，Test #77；5.48 那格当年同一条测试记的是 2.58 s）
./build/src/1d/apps/chflow-1d-transient tests/cases/transientPowerLawRheology/case.dict
python3 tests/check_transient_power_law.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/tests/cases/transientPowerLawRheology/case.dict"

# 屈服应力的静止支进瞬态那一行（5.49）：新算例 transientYieldStaticHold，五个 case 文件（1 个定压求解 +
# 1 个自由节点串联求解 + 1 个同拓扑稳态对照 + 2 个拒收）。5.48 留下的那一格现在有了物理：屈服管的 dp_wall(Q)
# 在 Q->0+ 时趋于 dp_static = 4 tau_0 L/d（transient.cpp:28-:31，只吃几何所以不认识 n），而 |drive| <= dp_static
# 时**那条奇延拓的律没有根**——稳态可以停在括号端点，瞬态的 dp 正是行在解的未知数，所以必须有一条静止支。
# 三件东西缺一不可：入口闩锁 holdsStatic（:40-:44，:310-:316 每步只决定一次；不闩就会被 hold 的管所见的 drive
# 随节点压力来回翻）、求解内的穿过钳制（:353-:362，x[k]*qPrev<=0 时写 x=0/r=0/J=1，这是有限时间停住那一半）、
# 静止处的值与斜率（frictionDp :72-:81 返回 copysign(dp_static, drive)；frictionSlope :89-:103 走"超额割线"
# (dp_wall-dp_static)/|q|，因为屈服律在静止点斜率是**无穷**，跨 q=0 中心差分会把阈值跳变除进商里得到 ~1e15，
# Newton 于是永远坐在静止点上——那是与 5.48 拒收方向相反的错）。drive 用 :52-:56 的 drivePa = pFrom-pTo-rho g dz，
# 与行自己的 dp-grav（:352）同一条约定，否则闩锁与行会对"哪端是上游"各有说法（本案 Hang/Roll 就是这条的判据）。
# 要看什么：
#   1) "恰好 0"要打在**字符串**上：四根管各自因不同理由打印单个 token `0`（壁高于所给压差 / 阈值不认识 n /
#      重力把 drive 拉到阈值下 / 静止在阈值之下没有根），而 2000/2000 步的 static hold 旗标是另一条独立打印
#      （transientApp.cpp:55-:64 报阈值、:84-:92 每样本报旗标）。用容差判"几乎不动"会让 1e-19 通过，也就会让
#      "静止支根本没生效"通过。
#   2) 停止必须**穿过**零并且与递推同步：StopBelow 从 2.0e-6 起步、第 17 步（t=0.034）达到恰好 0，停止前最后一个
#      速率 5.409e-08 m3/s，而"没有支时每步的位移" (dp_static-drive) dt/(rho L/A) = 7.392e-08 m3/s——两个数一起写
#      才说明差别是一条支而不是一个容差；checker 的递推自己也算出第 17 步（17 对 17）。
#   3) 节点账要量到**求解器自己的退出**：被 hold 的管另一端是自由节点时，体积账 V/(K dt)(p_i-p_{i-1})-(Q_in-Q_out)
#      关到 4.874e-13（C）/4.960e-13（D），对门 1e-11；注意 :414 的定标带一个 +1.0 m3/s，所以这里是**绝对**余量——
#      把它写成"相对 1e-11"是错的。同时给一个物理换算（4.960e-13 m3/s 是最大节点吞吐 1.657e-05 的 3.0e-08，作为漏
#      需要 1.008e+08 s 才排空 5000 Pa），否则读者无法判断这条余量是数值还是物理。
#   4) 静止支得能被**放开**：ReleaseLate 前 253 个状态打印 0，第 253 步（t=0.506）放开，判据打的是"入口 drive"：
#      节点 F 上一样本留 7996.7725 Pa、这一样本留 8006.8851 Pa，阈值 8000 Pa。只断言"后来不等于 0"不够——那只说明跑过。
#   5) 跨内核对照用**同一拓扑**的稳态文件：case_steadySeries.dict 与 case_series.dict 逐字同构，只差 solver 块；
#      放开那一对终点 8.8475721067e-07 对 8.8475719418e-07（1.864e-08）。同时必须断言两个被封节点**还没**到稳态
#      压力（22291.3 对 20000、25875.7 对 30000，4.0 s 只是 0.79 个弛豫时间 5.09 s）——不写这一条，读者会以为
#      速率对上就等于状态对上。
# 写这类判据时的坑：①drive 的重力符号写成 +rho g dz 时 Hang 会流动、Roll 会被 hold（与行的 dp-grav 相反）；
#   ②"每步位移"是 (dp_static-drive)*dt/(rho L/A)（单位 m3/s），误写成 /(I/dt) 得到的是 m3 体积，数字与单位一起错
#   三个量级，而这句话是要写进判据消息里的；③稳态 app 的 pipe 行比瞬态长（带 v/Re/f/dp）、节点行是 4 token，
#   两套 parse 不能共用一个 len() 条件。
# 十族 47 条；五个 case 文件，checker 墙钟 140.14 s（ctest -j8 记，Test #78；单独跑 122 s）。贵全部贵在 checker 自己那
# 10 x 2001 步的带支递推上：5 次 app 调用合计只有 0.31 s（case.dict 0.198 s、case_series.dict 0.107 s、
# case_steadySeries.dict 0.007 s），所以这一格的回归时间不是求解器的，剩下的都在 checker 自己的算术里
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldStaticHold/case.dict
python3 tests/check_transient_yield_static_hold.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldStaticHold/case.dict"

# 屈服壁与阀 K 第一次同一条行（5.50）：新算例 transientYieldValve，两个 case 文件（14 根管 + 3 个定压节点的
# 瞬态案，和同拓扑的稳态对照案）。零内核改动：5.49 的静止支与 M2 那一格就有的阀 K 损耗本来就在同一行里
# （transient.cpp:381 的 r[k] = I/dt (Q-Q_prev) - dp + grav + fr + vd + cc，vd 是 :370-:371 的
# K rho|Q|Q/(2A^2)，K 由 :293-:300 在 Newton 之前按 t+dt 定一次），这一格要证的是**两件账本没写错的事**：
# ①静止时阀不索要任何东西——dp_valve(0) = 0，所以闩锁/钳制看到的 drive 里没有 K（:52-:56 drivePa 只减重力）；
# ②凡是读壁律的地方，速率都是**两个损失合起来的根**（checker 在壁应力空间解 residual(tw) = I/dt(Q-Q_prev)
# + (4L/d) tw + dp_valve(Q(tw)) - drive，单调所以一次二分；residual(tau_0) >= 0 即无根即 Q=0）。
# 要看什么：
#   1) "阀不参与静止"要打在**同一拓扑不同 K 的字符串全等**上：Stuck 那组（tau_0 2.6、阈值 10400 Pa 对 10000 Pa
#      压差）K = 0 / 500 / 1e8 三根管 2001/2001 状态**打印逐字相同**、旗标 2000/2000 相同；Hang 那组（重力把
#      drive 拉到 2000 Pa 墙下）加 K=500 也逐字相同。注意**不能**写成"若把 K 折进闩锁的 drive，K=1e8 那根管
#      的数字会变"——dp_valve(0)=0 使任何这种变体都打印同一个 0，这条假的证伪性说法是错的。
#   2) 稳态侧的同一件事打在**缺失的 token** 上：8 根静止管里 app 一行都不印 dpValve
#      （steadyNetworkApp.cpp:65 只在 dpValve != 0 时追加），而行上仍带着 dp 10000 / dp 4000。
#      断言"dpValve 存在且为 0"会红，断言"没有 dpValve"才是这条打印约定的正确写法。
#   3) "同一个根"要**两笔都改**：Release 那组（阈值 4000 Pa、所给压差 10000 Pa）终端速率
#      9.33053017832e-06 > 9.25990098571e-06 > 8.99737384167e-06 > 7.20655693864e-06，阀分别占
#      0 / 36.9234 / 174.2973 / 1118.1880 Pa（后者是 1e4 Pa 的 11.18%，把速率砍掉 22.8%）。判据不是"阀越大速率越小"
#      而是"这个速率就是**两条律之和**在该 drive 下的根"：瞬态对二分根最差 5.419e-10（RollRef）、稳态 9.091e-13。
#      再补一条反向的：把壁应力从根上挪 tw = 2.4908/2.4564/2.2205/4.1751 Pa，速率相对差 <= 8.039e-13 是**根对**
#      的证据；只比大小关系会让"阀项写成线性"这种错通过。
#   4) 停住的那一步要由**两个数**夹住：Stop 组停在哪一步取决于 q* = (dp_static - drive) dt/(rho L/A)
#      = 7.391983e-08 m3/s——三根管停在第 17 / 17 / 5 步，停住那一步的 q_prev 是
#      5.408929e-08 / 4.647801e-08 / 2.607645e-08（< q*），前一样本是 1.390704e-07 / 1.306228e-07 /
#      1.184044e-07（> q*）。每步的减号量至少是 q*，因为 dp_wall(Q) >= 阈值对一切 Q>0 成立，所以"从上面穿过 q*
#      的那一步就是打印 0 的那一步"。停住瞬间阀的账单最多 292.8 Pa（是 2000 Pa 余量的 15%，K=5000 那根只有
#      4.651e-02 Pa）——这才是"K 进不了静止判据"的量级证据。
#   5) 幂次不能混：同一 K 在两根管上**花掉相同比例的 drive 是假的**（实测 8.92% vs 1.74%）。真不变量是阀的 Q^2：
#      速率比 3.0629 时阀降比是 9.3811 = 3.0629^2，而壁那笔的速率指数只有 0.474（屈服附近 Q ∝ (tw - tau_0)^2）。
#   6) 整条轨迹交给 checker 自己递推：14 x 2001 里动着的 6 根逐步对照，最差 6.890e-10（第 16 步的 StopK5000），
#      静止的 5 根逐位 0.000e+00；递推也算出同样的 [17,17,5]。
#   7) 跨内核只比**同拓扑**的稳态案：t=4 最差 5.428e-10（RollRef），且本案最大 Re 16.6492167083，闭式在读的是层流半区。
# 写这类判据时的坑：①checker 里 q_of_tw 必须带 5.47 剖面积分的 d^(1+1/n) 因子（写成 * d * d）。漏掉会得到**一半**
#   大小的速率、31 条断言同时红，而"壁应力空间二分"的残差看起来仍然单调——错在律不在解法。②稳态 app 的 pipe 行
#   在带阀时多出 dpValve，把 tok[-1] 当 dp 会把阀降读成总压降；按 key 定位。③瞬态 app 在第一个 t 头之前会打
#   每根管的 static-hold threshold，解析时要一个 started 标志，否则首状态被吃掉。④案例注释里阀降的量级要**实测**：
#   K=1e8、Q=2.0e-4 那笔是 1.7225e10 Pa（drive 的 4.3e6 倍），手推时把 A^2 少算一次会差四个数量级。
# 八族 81 条；两个 case 文件，ctest -j8 记 **8.51 s**（Test #79），单独跑 checker **5.0 s**，其中两次 app 调用只占
# 0.166 s + 0.043 s——这一格便宜在 checker 把隐式步解在**壁应力空间**（单调所以一次二分），而不是像 5.49 那样在速率空间
# 重解带支递推：同族那条 5.49 测试在**本轮同一负载**下记 338.38 s（上一轮轻负载记 140.14 s），差的从来不是求解器
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldValve/case.dict
python3 tests/check_transient_yield_valve.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldValve/case.dict"

# 屈服管的层流-湍流转捩 locus 第一次落在稳态那一行（5.51）：新算例 steadyYieldTransition，一个 case 文件
# （17 根管 + 6 个定压节点、5 个桌面压差 200/500/5000/3e4/1.5e5 Pa）。case.dict **没有新键**：内核新增
# yieldedTransitionOf（network.cpp:413）+ PipeFlowResult::transition，solvePipeFlow 逐管填一次（:501），
# app 打 He / ReB / ReBt / Qt / dpQt / Vrelt / band + past/belowTransition。这一格只报 locus、**不重定价摩阻**：
# Wilson & Thomas 2006 的湍流摩擦律（Eq 15/16，含 beta = tau_w/tau_0 的括号与 11.6/chi 钝化项）在无 poppler 的
# stdlib 抽取下上下标位置丢失，两篇原文抽出来的形式互相差一个 beta^0.5，且按抽出的写法 beta->无穷时 11.6/chi
# 发散、退不回牛顿极限——未确认的转写不进内核。要看什么：
#   1) locus 属于**流体和几何**而不属于状态：判据要打在"同几何孪生管在压差相差最多 750 倍时六个 locus 量逐位
#      相同（0.000e+00）"上，而不是打在"locus 等于某个数"上——后者对手写常数也通过，只有前者会在 locus 被从速率
#      反算时分开。被静态阈值顶住（Q 0）的管必须照样报出全套 locus，这是这条性质最锋利的一支。
#   2) dp_t 是**这根管自己的壁律**：把 Q_t 送回 5.46 的 Buckingham 括号重解 tw 再乘 4L/d（checker 自己二分 200 步），
#      不是相关式；同时断言 dp_t > 4 tau_0 L/d——只写"> 阈值"会让"转捩就是离开庭"这个错读法通过。
#   3) 判决要**三个入口交叉**：pastLaminarLimit == (Q > Qt) == (ReB > ReBt) == (dp > dpQt) 逐管，且案里必须有管真的
#      跨在两侧（实测 7 past / 8 below）；全侧一致的案证不了比较本身。
#   4) "读 locus 没碰任何一笔账"要**两条都断**：速率仍是所问压差在层流屈服律下的根（实测最差 3.2e-12），且
#      f Re_B g(xi) = 64 在**超出自己 locus 151 倍**的管上仍成立（3.4e-12）。第二条是留给下一格的：湍流摩擦闭包
#      落地那天它必须变红，所以这批断言要一起改掉，而不是到时候删掉。
#   5) 两个带的接缝是**量出来的向下台阶**（He=1e5 处 -1.17%、He=1700 处 -25.3%）。别写成"分段函数连续"，也别顺手
#      平滑——论文跨接缝用的 Eq (17) 常数没转写。注意方向：我第一版假设高带在中带上方，实测相反，跑一次再写判据。
#   6) 分组层面的反面：牛顿管与 HB 管**不打** He，所以判据断"token 缺失"而不是"值为 0"；HB 行打印的 Re 是表观数
#      8 rho v^2 / tau_w，与塑性数 rho v d / mu 不是一回事。
# 写这类判据时的坑：①app 用 setprecision(12) 打印，**所有"逐字"判据的下限就是打印舍入**（三族分别撞在
#   3.04e-12 / 2.14e-12 / 2.12e-12），把门写成 1e-13 会得到一屏假红；上限取 2e-11 并在消息里说明下限是打印。
#   ②checker 的 run() 用 cwd=CASE_DIR，参数必须传**绝对路径**，否则 FileNotFoundError。③不要用位置解析 pipe 行：
#   带 locus 的行多出八个 token，按 key 扫描。
# 七族 60 条；ctest -j8 记 **0.29 s**（Test #80）——全表最便宜的一格，因为 locus 是纯读出：一次 app 调用、零步进、
# 零重解网格；既有测试一字未改。
./build/src/1d/apps/chflow-1d-steady tests/cases/steadyYieldTransition/case.dict
python3 tests/check_steady_yield_transition.py \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/steadyYieldTransition/case.dict"

# 阀的时间表与屈服壁的静止支第一次同处一条行（5.52）：新算例 transientYieldValveSchedule，两个 case 文件
# （12 根管 + 3 个定压节点；4 根阈值以上带关阀/循环/关到底/写死关闭，4 根阈值以下从 2.0e-4 起步做停住对照，
# 2 根已经静止在阈值下，2 根去掉 yieldStress 的牛顿孪生；case_steady.dict 是同拓扑、每个 K 冻结在 t=4 值）。
# **内核一字未改**：阀项、穿过钳制、Newton 全是 5.38/5.49/5.50 装好的，这一格需要的是被量而不是被建。
# 两个开关读的不是同一个量——holdsStatic（transient.cpp:40-:43）比 drive 与 4 tau_0 L/d，drivePa（:52-:56）里没有 K；
# valveKAt（:125，按 t + dt 采样并在 Newton 期间保持，:293-:300）只给在动的那一行定价。要看什么：
#   1) 交叉两个开关时，判据要打在**哪个量碰得到哪个碰不到**上。最锋利的一支是"逐字不可见"（StopLate == StopPlain
#      在 2001 个 token 与 2001 面旗标上全同），但它必须与"同一形状的日程在别处把停住分成 54 / 72"同案存在才可
#      失败：只拿已经被 hold 的管做日程，恒等字符串是白拿的。
#   2) 日程的读数是**停住的步号**而不是某个时刻的速率：断 stop_step(app) == stop_step(递推) 且单调
#      （184 / 72 / 54 对 从不 / 0.6 / 0.15 收尾），比断速率更能说明"墙在停、日程只说什么时候"。附带一条物理：
#      被阀提前判停的管停在**更大**的速率上（7.08e-08 / 4.25e-08 对无阀的 3.28e-10），因为那次穿过是钳制给的。
#   3) 采样约定（t 对 t + dt）要**两个方向都证**：斜坡与动态重叠的五根管按 t 采样分别 miss 5.526e-02 / 3.055e-01 /
#      1.000e+00，而从不重叠的两根管两种采样差 **0.000e+00**。后者不是"失败的判据"而是"这个量在此处不可观测"的
#      陈述——把它写成 > 门会得到一条永远为假的断言。最干脆的是步号读数：app 停在 54 / 72，按 t 的递推停在 55 / 73。
#   4) **相对门只在有速率的地方成立**。内核 Newton 的退出定标（:389-:391）用的是 drive + 壁降（~1e4 Pa）而不是速率，
#      所以停住前那几个 Q ~ 1e-10 的状态，相对差是求解器精度而不是律：本案同时立两条——|Q| >= 1e-8 上相对最差
#      2.140e-09，与全部 24012 个状态绝对差 <= 2.666e-15。只写"全域 rel < 5e-10"会得到 2.192e-07 的假红（第一版就是这样）。
#   5) 阈值以上那一组要**同时**给"闭阀值多少钱"（K=1e6 在开阀速率上 98192 Pa = 9.8 倍整条 drive）与"held 0/2001"，
#      否则"阀没起作用"读起来像 bug 而不是判据；再补一条尾巴等值（关到底 vs 一开始就写死关闭 rel 0.000e+00），
#      说明日程在斜坡结束后不留痕迹。
#   6) 稳态交叉里"字段不出现"是被证的（阈值以下六根管打 Q 0 而**连 dpValve 都不出现**），别用 == 0 代替 token 缺失。
# 写这类判据时的坑：①run() 用 cwd=CASE_DIR，参数必须是**绝对路径**，第二个 case 文件用 f"{CASE_DIR}/case_steady.dict" 拼。
#   ②瞬态的 pipe 行有两种长度（`Q` 四个 token、`hold` 五个），按 (tok[0], tok[2]/tok[3]) 分派而不是按长度或位置。
#   ③斜坡值别做浮点等值：K(0.1) 手算与 valveKAt 差 2.3e-16（0.075/0.125 那两个二进制可精确表示所以恰好通过），
#      用 rel() 或改用二进制友好的常数。
# 七族 121 条；ctest -j8 记 **11.14 s**（Test #81；单独跑同量级 8.45 s，成本在 checker 自己那 24012 次壁应力空间的二分，
# 两次 app 调用可忽略）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldValveSchedule/case.dict
python3 tests/check_transient_yield_valve_schedule.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldValveSchedule/case.dict"

# 浮动节点 × 屈服壁的静止支 × 阀的时间表第一次同处一题（5.53）：新算例 transientYieldValveNodeRelease，三个 case
# 文件——主桌面 11 节点 / 18 根管 / 九个自由节点（每对：τ_0=2.0 的屈服管 A→F 串一根牛顿泄放管 F→B，阈值 8000 Pa、
# 入口 drive 5000 Pa、余量 3000 Pa，F 带 volume 0.1 / bulkModulus 1e8 即 Cc = 1e-9 m3/Pa）；case_long.dict 只跑 4 s
# 里不放开的那两对、endTime 30.0 / reportEvery 5；case_steady.dict 同拓扑、每个 K 冻结在 t=4 的值。
# **内核一字未改**：静止支、穿过钳制、阀项、节点行、采样约定全是 5.38/5.49/5.50/5.52 装好的。要看什么：
#   1) 桌面从"两端定压"换成"中间一个自由节点"之后，闩锁比的 drive 是**解出来的**，于是"日程在哪一行"成为一条
#      可判的差别：写在泄放管上的 K = 0/5e3/1e6 把放开推到第 253/299/1830 步（单调），写在被 hold 的屈服管上的
#      同一把阀连放开那一步都不移动（253，节点历史逐位相同到第 252 个状态）。放开步号是本案最好的观测量。
#   2) **同一个 closedK、同一段斜坡长度，两个结局**：0.2→0.3 s 到达的斜坡取消放开，1.0→1.1 s 到达的在 K0 那一步
#      之后才可见（第一个不同 token 在 501 = int(1.0/dt)+1）。别把这条写成"晚到的日程无效"——5.52 的恒等字符串
#      在这里不再成立，因为节点先动、墙后动。
#   3) **网络延迟要量，不要假设**：W 的泄放管 2.0→2.1 s 开回来，节点从第 1001 个状态分岔，而被 hold 的管仍打印
#      0 到第 1289 个状态，放开在 1290 步（比日程第一次动东西晚 0.578 s）。第一版按"斜坡结束即生效"写会红。
#   4) **稳态孪生带自由节点时必须有 sources 块**（还要 solver { type steadyNetwork; } 与自由节点写成 F {}）：
#      缺它残差不是"答案不对"而是冲到 3.8116e+279 并报 NOT converged，看起来像数值 bug 其实是定标缺项。
#   5) **跨内核终点的门要逐对写**：尾巴 rel 从 1.84e-08（放开后只剩 3.49 s 纯弛豫）到 1.24e-03（封到 2.1 s 才开），
#      且五格一律速率比压力离落点更远 5.9–53.8 倍（近壁 q ∝ (tw - tau_0)^2）。一个统一的 2e-9 门第一版就红了；
#      同理"还在路上"的两对要**承认在半路**（K6 rel 0.955 并给出按最后五采样 571.6 Pa/s 还欠 1.43 s），
#      而不是放宽门去假装收敛。
#   6) **"从不放开"要有第二个时钟 + 单调性**：30 s 那跑把 R1 在 24.25 s 交回、K8 到 30 s 仍差 304.35 Pa 且 drive
#      从 5000 Pa 起单调上升——这样 "never" 是一个读数而不是判决。再补一条算术（节点每步掉 0.2051 Pa = 102.55 Pa/s
#      = 泄放速率 / Cc）说明这个读数是可以外推的。
#   7) 节点密封残差（4.966e-14）要写成**求解器的数**：换算成占最大移动速率的比例（6.7e-09）与"当漏要 1.0e+08 s"，
#      别当作物理结论；跟 5.49⑦/5.52⑦ 同一格。
#   8) 断"走向稳态涓流"时**方向要先量**：K8 的泄放速率是从**上方**下降走向 6.10145e-08，写成 rising 会得到一条
#      永远为假的断言（本案第一版即如此）；同理"节点压回放开压力之上"结构上做不到，就只能断单向。
# 九族 104 条；ctest -j8 记 **8.51 s**（Test #82；单独跑 **4.61 s**，其中三次 app 调用合计只有 0.25 s，成本在 9 × 2001
# 加两对 2 × 15000 次壁应力空间二分与节点定点）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldValveNodeRelease/case.dict
./build/src/1d/apps/chflow-1d-steady tests/cases/transientYieldValveNodeRelease/case_steady.dict
python3 tests/check_transient_yield_valve_node_release.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldValveNodeRelease/case.dict"

# 分支节点 × 双向闩锁第一次同台（5.54）：新算例 transientYieldValveBranchRelease，三个 case 文件——主桌面
# 9 节点 / 21 根管 / 七个臂（每臂两根屈服壁 A→F 串一根带日程的泄放管 F→B；阈值 4τ_0L/d = 8000/6000/4000 Pa、
# 节点起点 2.7e4 对 pA 3.0e4 即起步 drive 3000 Pa，所以三档全从 held 出发；泄放全用同一段 close 0.8→1.0 /
# open 2.5→2.7 s，只有臂 e 提前到 0.2→0.4 s）；case_long.dict 只留臂 a、把重开推到 15.0→15.2 s、endTime 18.0；
# case_steady.dict 同拓扑、所有阀冻结在 K=0（仍要 solver { type steadyNetwork; } + 自由节点写成 F {} + sources 块）。
# **内核一字未改**：节点行对任意度数写（transient.cpp:395–:416）、闩锁是一条状态测试（:40–:43 的 Q_prev == 0 且
# drive 在阈值下）、穿过零的内步钳本就在（:353–:355）。要看什么：
#   1) 本案的地基是"第二根墙给节点一个天花板"：一根墙时节点从下方逼近**自己的**放开压力并停住（5.53 边界⑥），
#      两根墙时只要弱墙还在流，极限就是**弱墙的**放开压力，强墙在上行路上从它下面穿过、于是能被重新锁上。
#      第一版容易写成"每根墙各自闩锁"——臂 i（两根同为 8000）就是这条误判的反例：都在第 318 步放开、都不重锁。
#   2) 重锁与放开读**同一个数**，所以不要引入迟滞量去拟合：14 次放开事件（11 根壁）的起步 drive 超阈值 1.00–16.54 Pa、
#      3 次重锁差 3.10–3.55 Pa 不到同一个数，两侧都恰好是"一步的节点行程"。写成"回到阈值九成才锁"会红。
#   3) 严格不等式是**必要**而不是充分，两种落空各有各的读数，要分开写：St 从未被问（该节点最低 22593.53 Pa，
#      比它需要的 22000 高 593.53 Pa）；Se 放开后节点再没回去（放开之后最大 21979.46 Pa，差 20.54 Pa；
#      稳态落点 21379.73 差 620.27 Pa）。把"law ⟺ re-latch"当成全域等式在臂 t 与臂 e 上都立不住。
#   4) "到不了"要有第二个时钟加收缩因子才算读数：长跑 7001 步密封里弱墙缺口 1048.65 → 394.36 Pa，每千步收缩
#      0.55252 / 0.79305 / 0.91052 / 0.96289 / 0.98491 / 0.99392。只断"never"是判决，断因子趋 1 才是外推。
#   5) 静默的行按 **token** 断而不是按容差断：臂 t = 臂 s + 一根永不醒的壁，两者节点全程差 0.0 Pa、
#      Wt/Ws 与 Tt/Ts 与 Dt/Ds 在 2001 个状态里逐字相同。写 abs < 1e-12 会把"这两个桌面其实是同一条行"糊掉。
#   6) 长跑与短跑**不逐位相同**：共享的 1251 个状态里节点最大差 2.00e-07 Pa（相对 6.7e-12），因为 21 条边与
#      3 条边把同一张稠密行求和的次序不同。别照抄 5.53 的"逐位"措辞——这里对齐的是翻步 + 数值门。
#   7) hold 期间节点**不一定上升**：臂 e 在两根壁都还被 hold 时关上，节点随后整段向下（每步 −0.2201 ~ −0.2251 Pa
#      = dt·Q_泄放/Cc，那段穿过 1e8 的均值 1.11284e-07 m³/s）。"密封把节点抬到天花板"只在至少一根入流在动时成立。
#   8) 守恒门写成求解器的数：19658 个被 hold 的状态全打印单个 token 0，节点平衡只用打印 token 闭合到
#      4.951e-14 m³/s；稳态内核不带闩锁却恰好拒付 4 s 从未叫醒的那两根（Sp 5000.00 / St 7342.07 对 8000）。
# 九族 50 条；ctest -j8 记 **9.69 s**（Test #83；单独跑 **8.05 s**，其中三次 app 调用合计 0.76 s，成本在
# 7 × 2000 + 9000 次壁应力空间二分与节点定点）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldValveBranchRelease/case.dict
./build/src/1d/apps/chflow-1d-steady tests/cases/transientYieldValveBranchRelease/case_steady.dict
python3 tests/check_transient_yield_valve_branch_release.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldValveBranchRelease/case.dict"

# 屈服壁第一次挂在节点的流出侧（5.55）：新算例 transientYieldValveSplitRelease，五个 case 文件——主桌面
# 10 节点 / 19 根管 / 七个臂（`a` 保留 5.54 的两进一出作对照；`b`/`c`/`z` 一根进两根出；`d`/`f` 一进一出
# 且两侧都是壁、其中出侧那根才是本案的新格；`m` 一根带日程的行同时挂在两个自由节点上），14 根屈服壁 + 5 根牛顿行；泄放壁的 drive 是
# `pF - pB`，节点起点 2.7e4 对 pA 3.0e4 / pB 2.0e4，所以同一根 8000 Pa 的壁在入侧看到 3000（睡）、
# 在出侧看到 7000（也睡），而 6000/4000 的壁一开局就在流：第一次决策 9 睡 5 醒。case_long.dict 只留
# b/f/z、重开推到 15.0→15.2 s、endTime 18.0；case_steady.dict 同拓扑全部阀冻结在 K=0；
# case_steady_dead.dict 只留死带三臂；case_steady_slow.dict 是 7 臂把 maxIter 换成 50000。
# **内核一字未改**：holdsStatic 取 drive 绝对值（transient.cpp:40–:43）、drivePa 只是拓扑顺序
# （:52–:56）、Q=0 时 frictionDp 抄 drive 符号（:72–:81）、节点行照度数求和（:395–:416）。要看什么：
#   1) 镜像不是改名：入流壁的极限是天花板 `pA - 4τ_0L/d`（从下方逼近），泄放壁的是地板
#      `pB + 4τ_0L/d`（从上方逼近）。"有几根从 held 出发"因此本身就是一个读数（本案 9/14），别按 5.54 的
#      "全从 held 出发"照抄——`Tz`（4000）从第 2 个状态起就在流。
#   2) 重锁律的正确陈述是"**同一侧**、阈值**严格更小**的另一行"，载体可以是壁**也可以是一根牛顿行**：
#      本案预测 `Sa`←`Wa`、`Wz`←`Tz`、`Wb`←`Db` 三根，trace 命中三根，其余 11 根无载体、零重锁。写成
#      "同一节点有两根壁"会在臂 b 红；丢掉"同一侧"会在臂 c 红（节点度数 3，但 `Sc` 是唯一供给者，
#      第 173 步醒一次后再没睡回去）。
#   3) 单步穿越是一类翻、不是全部：`Sa`/`Wz` 重锁时 drive 在阈值下 3.10/1.73 Pa，`Wb` 却在 **99.04 Pa**
#      之下才被锁——牛顿行在 drive 已退回墙内之后还在减速液柱。把"|drive-阈值| < 一步"写成全域断言会红。
#   4) 地板"到不到"要有第二个时钟加收缩因子：18 s 里 `Ff` 475.39→82.54、`Fz` 585.51→60.44 Pa，臂 z 每千步
#      收缩 0.43708→0.73507（处处 <1）。只断 min() > 地板是"还没到"，断因子才是"到不了"。
#   5) 死带要按**整个解**判而不是按行：臂 d 两行都 hold 时节点行退化成 `0 = 0`，`Fd` 全程同一个 token
#      27000.0 Pa、旗标 2000 个全是 1。这类臂在稳态内核里是 `iterations=0 residual=0 (converged)` 且节点
#      打印求解器自己的猜测 25000——**那是秩亏不是失败**，所以 banner 不能当唯一判据。
#   6) 秩亏的**局部性**要用两个文件证：那 7 行 + 4 个节点在 7 臂与 3 臂文件里逐 token 相同，才敢说 7 臂
#      文件的不收敛来自**会流的两臂**（残差 0.0321364485 × 6e-5 = 1.928187e-06 m³/s 正是臂 a 的节点失配）。
#   7) 换 maxIter 时分开断"会动的行 / 不动的行 / 未定的行"：`Fa` 21898.4842→21379.7317、`Fz`→24000.0401
#      随预算走，`Fb`/`Fc` 逐 token 不变，`Fd`/`Ff`/`Fm*` 恒 25000；且 50000 步后残差 2.6336e-11 仍高于
#      tol 1e-12（banner 比的就是它，network.cpp:689）——这是**速率**不是拒绝，措辞要跟着数走。
#   8) march 与 steady 对同一臂可以给相反的真值：`Ef`/`Em` 在 march 末态 drive 6303.12/6520.43 对 6000
#      仍在流，steady 把它们归零（因为 steady 要落点、march 只逼近地板）。别把 steady 的 0 读成瞬态的 bug。
#   9) 守恒门写成求解器的一个末位：8 个节点行 × 2000 步只用打印 token 闭合到 4.961e-14 m³/s，即
#      1e-7 Pa × V/BC/dt = 5.0e-14 的一个单位；相对失配只长到 1.889e-03 且只在 `Fm1` 第 2 步（行两侧
#      都是 8e-12 量级），所以要按**绝对**门判、把相对门单列一条。
#  10) `t = 0` 那一列的旗标全打印 0（还没决定），"旗标开着却打印 0"的三行恰是三根壁各自的钳步
#      （`Sa` 506、`Wb` 82、`Wz` 192，都在自己重锁的前一个状态）。事件扫描从 i=2 起，否则初始列会被读成
#      一次翻，14 根壁会多出 14 次假翻转。
# 十一族 70 条；ctest -j8 记 **49.10 s**（Test #84；单独跑 **32.02 s**，其中 case_steady_slow.dict 那趟
# 50000 次迭代就占 22.1 s，成本在逐状态 38019 个速率的壁应力空间二分与节点定点）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldValveSplitRelease/case.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldValveSplitRelease/case_long.dict
./build/src/1d/apps/chflow-1d-steady tests/cases/transientYieldValveSplitRelease/case_steady.dict
./build/src/1d/apps/chflow-1d-steady tests/cases/transientYieldValveSplitRelease/case_steady_slow.dict
python3 tests/check_transient_yield_valve_split_release.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldValveSplitRelease/case.dict"

# 死带的停泊场：把自由节点自己的初值扫开（5.56）：新算例 transientYieldParkingBand，四个 case 文件——
# 主扫描 14 节点 / 24 根管 / 十二个臂（`a`–`k` 是节点初值 20000→30000 的十一个采样点，`m` 是 `f` 那份
# 再给两根行各挂一份 K→1e8 的关-开日程），每臂一副 `A → F → B` 串联、两根 τ_0=2.0 的壁（阈值都是
# 4τ_0L/d=8000）；case_long.dict 同十二臂推到 40.0 s（20001 个状态）；case_single.dict 只留 `c`/`a`/`d`
# 三臂，用来付"停泊场是节点自己两行的性质"那张账；case_steady.dict 是同一批串联壁交给稳态内核。
# **内核一字未改**：整个扫描只读两行代码——holdsStatic 比的是 `<=`（transient.cpp:40–:43）、frictionDp
# 在 Q=0 抄 drive 的符号（:72–:81），两者一相遇节点行就成了 `0 = 0`（:395–:416）。要看什么：
#   1) 死带的**宽度**只有扫出来才存在：5.55 那条臂只给"这一个初值会停"，本案 `[22000, 28000]`、宽
#      6000 Pa，且**两端各被夹到 1 Pa**（`Fc 21999` 动 / `Fd 22000` 睡，`Fh 28000` 睡 / `Fi 28001` 动）。
#      端点必须**当成臂来算**：醒着那根的 drive 恰等于阈值，本案第一次决策就 HELD，读的是 `<=` 不是 `<`。
#      只跑中间一个初值，`<` 与 `<=` 两种内核给出一样的输出。
#   2) "逼近但到不了"要三条腿一起断：不越界（4 s 里 max 21999.00196 / min 28000.99804）、40 s 后仍差
#      0.9807–51.5302 Pa、以及**预测**这条尾巴的系数。系数从内核自己的 Bingham 剖面来：`q_of_tw` 的方括号
#      `1 - 4x/3 + x⁴/3` 在 x=1 处**常数项与一次项同时相消**，所以近边速率对超出量是**二次**
#      （实测 d ln Q/d ln gap = 1.9647–1.9974），`d(gap)/dt = -k_Q gap²/Cc` 积分成 `1/gap` 随时间走直线、
#      斜率 K = k_Q/Cc = 4.9087e-4 Pa⁻¹s⁻¹。断斜率（最差 1.4468%）+ 断积分落点（最差 3.619%）+ 断每 2 s
#      收缩比单调升向 1 而不越界（0.22988→0.88896）才是"1/t 尾巴"；只断 min() > 边只是"还没到"。
#   3) 二次律不是想当然：先把方括号在 x=1 的展开核一遍再写预测式。本案第一版把系数写成层流阻力那支
#      `(THR/R_LAM)·2·(…)/Cc`，得到 0.0019635，整整大 4 倍——真正的来源是 `q_of_tw` 的小-d 极限。
#   4) 别把跨矩阵尺寸的对照写成"逐 token 相同"：睡着的臂（`d`）确实 2001/2001 token 相同，但会动的臂
#      两种矩阵把同一批行按不同顺序求和，实测差到一个打印量子（1.000e-07 Pa、速率 3.800e-16 m³/s、
#      `Fa` 只有 1306/2001 个 token 相同）。判据要分成"冻结臂按 token、滑动臂按数 + 一个量子兜底"。
#   5) 停泊是**只能选起点**的口袋：醒过的壁速率从来不是恰好的 0（本案最小打印值 4.78156e-13 m³/s），
#      所以闩锁再也不会武装；`m` 那份 K 涨到 1e8 的日程与 `f` 逐 token 相同（HOLD 的行里没有 K 项可读，
#      :353–:355）。这两条合起来才是"日历碰不到停泊的行"，缺一条就只是巧合。
#   6) checker 自己的打印精度会把结论吃掉：`%.4g` 把 21999 打成 2.2e+04、`%.7g` 把 21999.0000019 打成
#      21999。要断"它在动"，消息里的数至少给到 `%.6g`/`%.10g`。
#   7) 稳态没有闩锁，别把它的 banner 读成物理答案：**22** 根行全 `Q 0 dp 5000`、**11** 个节点全打印
#      25000、`iterations=0 residual=0 (converged)`——未定节点被留在定压的平均值，而
#      `(pA-4τ_0L/d + pB+4τ_0L/d)/2 = (pA+pB)/2` 正是被扫出来的带的**中心**，所以稳态的盲起点必然是
#      一个停泊态。秩亏不是失败；要说"稳态也能收敛在未定解上"，得像本案一样把那个值和瞬态扫出来的
#      带对上。
#   8) `t = 0` 那一列 **24** 个旗标全打印 0（决策还没发生），事件扫描从 i=1 起，否则 6 条停泊臂各报一次
#      假释放；同时 **48024** 个速率 token 里 **36024** 个恰为 `0`，节点行只用打印数字闭合到 1.189e-06
#      m³/s（每步 1.189e-15 Pa）。
# 十一族 49 条；ctest -j8 记 **26.09 s**（Test #85；单独跑 **14.56 s**，checker 直跑 14.24 s，成本在 20001 状态 × 24 行的
# 逐状态解析与壁应力空间二分复现），全量 **111/111** 绿（本轮实测 **354.01 s**）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldParkingBand/case.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldParkingBand/case_long.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldParkingBand/case_single.dict
./build/src/1d/apps/chflow-1d-steady tests/cases/transientYieldParkingBand/case_steady.dict
python3 tests/check_transient_yield_parking_band.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldParkingBand/case.dict"

# 死带 × Herschel-Bulkley：把同一张停泊桌面搬到认识指数的壁上（5.57，5.56 边界②）：新算例
# transientYieldIndexParkingBand，七份 case 文件——三副七臂矩阵（`a` 25000 / `b` 22000 / `c` 28000 在带里，
# `d`/`f` 在外面一格 1 Pa，`e`/`g` 在外面两千格），每臂 `A → F{tag}{x} → B` 串联两根 HB 壁
# （`length 20; diameter 0.02; yieldStress 2.0` + 逐管 `flowIndex`/`consistency`），三档只在
# (n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02) 上不同；`case_long_*` 把同样七臂推到 40.0 s（20001 个
# 状态）；`case_steady.dict` 把三档各一份串联壁交给稳态内核。**内核一字未改**。要看什么：
#   1) 先分清"哪一半认识指数"。阈值 `yieldThresholdPa = 4τ_0L/d`（transient.cpp:28–:31）里没有 n 也没有 K，
#      所以带、两个吸引子、`holdsStatic` 的 `<=`（:40–:43）在 HB 上与 Bingham 逐位相同（三档都停成单
#      token 25000/22000/28000，端点 `b`/`c` 第一次决策就 HELD、drive 恰为 8000）；认识指数的是闩锁松开
#      之后的那一行。把这两半混成一句"HB 改变了死带"会把本案读反。
#   2) 从打印 token 反解壁应力时**不要给下行加符号**：S 行的 drive 是 `pA - p_F`、E 行的 drive 是
#      `p_F - pB`，两支都是"从高侧往低侧流"，所以两支行打印的 Q 同为正。本案第一版按"下行反向"取了
#      `sgn = -1`，重构误差立刻涨到 5e1 相对；去掉那个翻转后括号**从第一次已决策那一步**就闭合
#      （1312 个采样速率最差 4.698e-07 / 3.612e-07 / 2.874e-07 = 打印量子）。
#   3) 第一步是斜坡不是稳态：醒着的行把 539.5 / 1635.8 / 1821.0 Pa（8000 的 6.7% / 20.4% / 22.8%）交给
#      流体惯性 `(inertia/dt)·ΔQ`（:381–:382），壁应力因此起步 2.365 / 2.091 / 2.045 而不是静止的 2.500。
#      这笔不是污染，是行的一部分——把它当残差丢掉，"括号只在中段闭合"的假结论就来了。
#   4) 指数有地理：`m = (n+1)/n` 只在近阈成立（括号在 δ=1 nanoPa 处读 3.000/2.000/1.500），远阈是纯幂律
#      `1/n`（δ=5000τ_0 处读 2.00007/1.00005/0.50003）。本案 gap 跨 0.4–2000 Pa、δ/τ_0 < 0.5，任何轨迹
#      都没走到远端，所以弦必须断"严格落在 (1/n, m) 内 + 与括号同两个 gap 的弦一致"（18/26/34 条，
#      0.027%/3.922%/2.023%），不能断"实测指数 = m"。
#   5) 采样要躲开噪声：峰值之前与 gap 比小于一个打印量子的配对都会造出 0.87 这种假弦。本案只取**速率峰值
#      之后**、且两点相对杠杆 > 1e-3 的配对；这样 `n=0.5` 的近边臂一条弦也不剩——那不是缺陷，是结论：
#      1 Pa 起点在 40 s 里动不了一个量子（0.000000 Pa），所以 5.56 的"1 Pa 夹边界"是 n=1 的便利。
#   6) 直线从 `1/gap` 挪到 `gap^(-1/n)`：预测式 `A_n/n`（`A_n = πR³(n/(n+1))(d/4L)^m/(τ_0 K^{1/n} Cc)`）
#      给 1.636246e-08 / 4.908739e-04 / 1.463503e-02，实测差 1.859% / 0.021% / 0.011%，且 n=1 那一档**就是**
#      5.56 的系数——所以 5.56 的"常数项与一次项同时相消"是 m=2 的特例，一支公式覆盖两节，不是两支。
#      注意本案 n 与 K 是**成对**改的，陡度差 2.9e4 / 30 倍要归给这一对。
#   7) 别把 HB 的 n=1 当成"另一个基线"去对照：`n=1, K=mu` 时通用括号与 5.56 的 Bingham 括号在 3999 个壁
#      应力上差 6.3e-16，4 s / 40 s 落点与 5.56 发布的数差 0.00e+00 – 3.00e-07 Pa。跨增量对**数**而不是对
#      token，是这种"同一支公式的第二半"最便宜的证法。
#   8) 稳态那侧的 banner 依旧不是物理答案：**6** 行全 `Q 0 dp 5000`、**3** 个未定节点三档都打印 25000、
#      `iterations=0`——`(pA+pB)/2` 是带中心的恒等式，与指数无关，因为阈值与指数无关。
# 十族 106 条；ctest -R 记 **26.61 s**、ctest -j8 记 **46 s**（Test #86；checker 直跑 26.45 s，成本在 3 档 × 2 个时钟 × 14 行的
# 逐状态壁应力二分复现，最尖的 `case_long_5` 单跑 10.31 s），全量 **112/112** 绿（本轮实测 **353.63 s**）；
# 既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexParkingBand/case_5.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexParkingBand/case_long_1.dict
./build/src/1d/apps/chflow-1d-steady tests/cases/transientYieldIndexParkingBand/case_steady.dict
python3 tests/check_transient_yield_index_parking_band.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/build/src/1d/apps/chflow-1d-steady" \
    "$PWD/tests/cases/transientYieldIndexParkingBand/case_5.dict"

# 阀日程 × HB 停泊带：日程只活在"闩锁松开之后的那一行"里（5.58，5.56 边界④ / 5.57 的"HB × 阀日程同台空白"）：
# 新算例 transientYieldIndexValveBand，三副十臂矩阵（`a`–`j`，20 根 HB 壁、12 节点、2001 个状态），三档只在
# (n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02) 上不同，`τ_0 = 2.0` 处处一样；日程 `0.5 → 1.0 s` 关到
# `closedK = 1e8`、`2.0 → 2.5 s` 重开（第一个受影响步 = 251）。**内核一字未改**。要看什么：
#   1) 先分清日程能不能碰到停泊那半边：`holdsStatic` 的判定里没有 `kv`，且持有行的钳位 `x[k]=0; r[k]=0; J[k][k]=1`
#      （transient.cpp:353–:355）在**读 `kv` 之前**执行，所以带日程的臂 `h` 与不带日程的双胞胎 `a` **逐 token 相同**
#      （2001 个节点 token + 2×2001 个速率 token），端点 `b`/`c` 在第一次决策就 HELD（drive 恰为 8000 = 阈值）、
#      K 爬到 1e8 也不醒。把这条读成"阀对 HB 没用"是过头的——它只对**停泊的行**没用。
#   2) 反解壁应力时那一笔 `K(t)·ρ|Q|Q/2A²` 不能省：带上它，括号**从第一次已决策那一步**闭合
#      （3.757e-10 / 1.001e-09 / 8.235e-10 相对）；省掉它，同一批状态读到 6.4e-01 / 1.4e+01 / 1.9e+01。
#      双胞胎（无日程行）在同一判据下是 3.658e-10 / 9.430e-10 / 6.909e-10——两读数一致才说明差值来自日程项而不是拟合。
#   3) 分割点就是日程的时钟，但**只有醒着的近边行按它分**：`a`/`h`、`b` 的速率在前 250 步逐位相同、第 251 步起分；
#      节点 token 却在 757 / 268 步才分（n=2 那档速率偏 4.90e-4，占自身 0.05%）。所以断言只能写"偏离 < 1% + 报告首个
#      不同的步"，写"逐 token 相同"会在最陡的档上假失败。
#   4) 成本要说成 Pa 而不是速率：节流比 **1.5926 / 7.2233 / 8.1974×**、40 s 落点少走的债
#      **24.1797 / 279.2209 / 253.1707 Pa**（占该起点到天花板 **167.8146 / 1540.2417 / 1822.4064 Pa** 的
#      14.4% / 18.1% / 13.9%）。用"峰值速率掉了多少倍"讲这件事会把 2.5 s 的开阀归还算进去，得出"几乎没损失"。
#   5) 日程给 drive 重新定价——这是本案唯一一处"日程改变了行里的账"：`K > 0` 那 999 步里阀价峰值
#      **323.0 / 2483.1 / 3701.8 Pa**，是阈值 8000 的 **4.0% / 31.0% / 46.3%**、该行自身 drive 的 **3.2% / 26.2% / 39.7%**
#      （分别出现在第 500 / 254 / 253 步——最陡那两档的峰值落在关闭爬坡**刚起步**处，那里 `K` 才 1.6e6 / 1.2e6，
#      但 `Q` 还没被压下来），壁应力因此比同一瞬间无日程的双胞胎错 **3.1822 / 10.2906 / 11.1758%**。
#      峰值不在 `closedK` 那一瞬而在爬坡途中——写检查时窗口要扫整个 `K > 0` 区间取最大，不能只取 t=1.0 s。
#   6) 开阀后的偿还单调且可数：重开后 **750/750** 步都比对照快（**1.03–1.04 / 1.59–3.57 / 1.93–5.31×**），
#      债从 26.3683 / 691.7019 / 848.1692 Pa 减到上面那三个数，**偿掉 8.3% / 59.6% / 70.2%**，且偿还是 n 的单调函数。
#   7) 带仍然进不去也穿不过：三档距边界 **1.00000 / 0.99804 / 0.89252 Pa**、关闭窗口内 20 根壁的最小速率
#      **8.1797e-18 / 4.8957e-13 / 2.6280e-11 m³/s**（非零、无一次 re-arm）。所以"日程 + 指数"没有造出第三个吸引子。
#   8) 镜像 `e`/`f`（同一个日程挂 S 行 vs 挂 E 行）依旧**精确**：节点和 50000 差 **0.0e+00**、醒着的速率
#      **2001/2001** 相同——因为损失项 `ρK|Q|Q/2A²` 是 Q 的奇函数，同一日程在两侧关掉同样的量。
#   9) 节点行的账要在**打印 token** 上闭合而不是在浮点上：`(Q_S − Q_E) − Cc·Δp/Δt` 三档最差
#      **4.9163e-14 / 4.9640e-14 / 4.9382e-14 m³/s**，单位要写"每步多少 Pa 的节点运动"而不是 Pa/s。
# 八族 115 条；ctest Test #87 单独跑记 **3.65 s**、`-j8` 全量里记 **6.91 s**（checker 直跑 3.5 s；0 内核改动），
# 全量 **113/113** 绿（本轮实测 **368.59 s**）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexValveBand/case_5.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexValveBand/case_1.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexValveBand/case_2.dict
python3 tests/check_transient_yield_index_valve_band.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/tests/cases/transientYieldIndexValveBand/case_5.dict"

# HB 壁 × 双向闩锁：把 5.54 的"重新锁上"搬到声明了 flowIndex 的壁上（5.59，5.54 边界③）：
# 新算例 transientYieldIndexValveRelatch，三副四臂文件（a/s/u/p，10 根管、6 节点、5 根 HB 壁、2001 个状态），
# (n, K) = (0.5, 1.0) / (1.0, 0.1) / (2.0, 0.02)，τ_0 = 2.0（阈值 8000、释放压力 22000）与 1.5（6000 → 24000），
# 日程关闭 0.8 → 1.0 s 到 closedK 1e8、重开 2.5 → 2.7 s（CLOSED = 500、OPENED = 1350），节点从 2.7e4 起步。
# **内核一字未改**（阈值 :28–:31、holdsStatic :40–:43、步内 clamp :340–:360、节点行 :395–:416）。要看什么：
#   1) 静止支有**两种**形状，别把它们读成一个：`held(i) = qPrev == 0 && |drive| <= 阈值`（:40–:43）与
#      `crossing(i) = qPrev != 0 && |drive| <= 阈值 && x[k]*qPrev <= 0`（:340–:360）。只有后者能让**正在动**的壁睡下。
#      本案 n=2 在第 465 步 drive 已落到 7998.4507 却还醒着打印 7.28766e-09 m³/s，旗标 7 步（1.4e-02 s）后才翻；
#      n=1 的 crossover 与 sleep 同步落在 507（滞后 0）。检查要**分别报** crossover 步与 sleep 步并断言 lag >= 0，
#      把"drive 下来了没睡"写成失败就是把内核的第二半删掉了。
#   2) "从未被问"不等于"被拒绝"：n=0.5 醒来之后入口 drive 的最小值仍是 8001.7883 Pa，桌面一次都没拿到判断机会。
#      判据写 `min(entry drive) > THR` 并另给原因——节点在整个密封窗 [500, 1350] 上最好只有 21400.49321 Pa，
#      距这根壁的释放压力 22000 差 599.507 Pa（n=1 / n=2 分别超出 1192.446 / 1544.178 Pa）。
#   3) 跨增量对账（n=1, K=mu 对 5.54 的 Bingham）要断言**事件身份 + 有界偏差**，不是"逐 token 相同"：旗标 2001 个状态
#      全等、Sa 阶梯 [1,333,507,1530] 与 Su 阶梯 [1,318] 两边同整数、节点 Fa/Fu 差 <= 1.0e-07 Pa（1920/1643 个 token
#      逐字相同）、速率相对差 <= 7.2e-09。只有壁从不流动的 p 臂才真逐位相同（Sp/Np/Dp 2001/2001）——那是"不碰任何律
#      的算术"该给的强度，把它推广到流动状态会假失败。
#   4) 睡眠窗口要**旗标与 token 双向一致**且是**一个连续块**：[507,1530) 与 [472,1568) 各 1023 / 1096 个恰好的 0。
#      另有一条只属于打印器的证据：唯一"醒着却打印 0"的状态是 506 / 471（= clamp 前一步），n=0.5 那 1682 个清醒状态
#      一个都没有；别处最小非零 9.243e-14 / 2.161e-13 / 1.651e-10 m³/s。判据写成 `zeros == [sleep_step - 1]`，
#      否则那一步会被读成第二个 held 事件。
#   5) 括号反解只在**有根**的步上成立：`tw = (drive − I_DT·Δq)/WALL_SCALE`，若 `tw <= τ_0` 截面积分没有根，必须跳过。
#      把跳过的步单独断言才是重点：整条轨迹里**恰好一个**未 held 步 implied 出低于屈服的壁应力，且它就是 clamp 前一步
#      （506 / 471；n=0.5 零个）——"held 不是残差，是行先失去根"。有根的那些闭合到 1.957e-08 / 4.561e-07 /
#      4.429e-08 相对（Sa）、3.463e-09 / 4.956e-09 / 2.818e-09（Wa/Ws）、1.539e-08 / 7.103e-09 / 7.153e-09（Su）。
#   6) 空断言要显式说出口：Sp 在 2000 个已决策状态全报 hold 1，括号一次都没被 offered。不要打印
#      "worst 0.000e+00 relative at step None"——改成"从未 offered"并断言 `all(flags[1:])`，或干脆不计入这条判据。
#   7) 供给/泄漏那笔账要在**第一个密封步**取（第 500 步：1.30083e-07 / 1.83744e-06 / 2.0945e-06 m³/s 对
#      5.43469e-08 / 6.73837e-08 / 7.0236e-08，比值 2.394 / 27.268 / 29.821），节点从 21292.23214 / 21982.03676 /
#      22152.29203 爬 108.2611 / 1210.4094 / 1391.8863 Pa。重锁的因在这个差里，不在速率里；弱壁独自峰值
#      1.27525e-07 对 1.83724e-06 / 2.0945e-06 就是"一个量级的供给 = 有重锁与没有"。
#   8) held 半边行的证据是"逐字相同直到强行第一次流动"（318 / 333 / 343 那一步分开），且 12 条节点行只用打印 token
#      闭合到 4.9451e-14 m³/s——单位写成"每步 9.89e-08 Pa 的节点运动"，别写 Pa/s。第二次醒来（1530 / 1568）距 K 归零
#      0.36 / 0.44 s，入口 drive 8001.0644 / 8002.4517 Pa：买回它的是重新打开的排放，不是流变。
#   9) Python 层的两个坑（本轮真踩过）：生成式最外层的 `range(ev[2], ev[3])` **在守卫为假时仍会求值**（先算可迭代对象
#      再过 if），要么整体走 if/else 分支；字面 `%%` 只有经过 `%` 运算的串才折叠成一个 `%`，硬编码文案里要直接写单 `%`
#      ——更稳的做法是把那个数也测出来插进去。
# 七族 87 条；ctest Test #88 单独跑记 4.92 s、`-j8` 全量里记 5.43 s（checker 直跑 4.96 s；0 内核改动），
# 全量 114/114 绿（本轮实测 372.08 s）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexValveRelatch/case_5.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexValveRelatch/case_1.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldIndexValveRelatch/case_2.dict
python3 tests/check_transient_yield_index_valve_relatch.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/tests/cases/transientYieldIndexValveRelatch/case_5.dict"

# τ_0 × flowIndex 交叉扫双向闩锁：把 5.59 的桌面再加一根 τ_0 轴（5.60，5.57/5.59 共同的空白格）：
# 新算例 transientYieldTauIndexRelatch，四臂 a/b/c/d 同拓扑只换强壁 τ_0 = 1.75/2.0/2.25/2.5（阈值
# 7000/8000/9000/10000 Pa，释放 23000/22000/21000/20000 Pa），弱壁固定 1.5，(n, K) 仍三对 ⇒ 12 格，
# 日程与 5.59 一致。另读 5.59 那三副文件做身份锚（臂 b 对 5.59 的臂 a）。**内核一字未改**。要看什么：
#   1) 重锁是**窗口**，不是半根轴：12 格普查 4 SLEPT / 2 AWAKE / 1 LATE / 5 NEVER。`n = 1 / 2` 的 a/b 睡下
#      （阶梯 [1,240,963,1403] / [1,333,507,1530] / [1,241,709,1439] / [1,343,472,1568]），τ_0 >= 2.25 那侧
#      **从不醒**（n=1/2 最大入口 drive 8674.0542 / 8577.5570 Pa，距 9000 差 325.9458 / 422.4430；n=0.5 的 d
#      最近 9802.5079，仍差 197.4921），n=0.5 那侧**从未被问**。断言要**分类 + 数数**（census），别写成布尔。
#   2) 桌面第三次盲，这次盲的是邻居的屈服应力：12 根弱壁全翻在 1 与 180 步，四臂前 239 / 240 / 241 个状态是
#      同一批字符串（弱壁速率 + 节点压力 alike），前缀结束的步号恰等于最早那次强壁醒来
#      （[239,318,1479,2001] / [240,333,2001,2001] / [241,343,2001,2001]）。held 行 = x=0, r=0, J=1
#      (:353-:355)，断言四臂在最早醒来之前**逐 token 相同**。
#   3) 分解是定量的：τ_0 买天花板（挪释放 3000 Pa，只挪封闭窗 [500,1350] 最好节点 51.36260 / 45.20672 /
#      44.18058 Pa），n 买补给（在固定释放压力下挪 2138.60589 / 2143.68512 / 2145.78792 / 2145.78792 Pa，
#      四臂跨度差 7.18202 Pa，比值 41.6）。最好节点单调下降于强壁 τ_0——更强的邻居供给更少，尽管它要的更低。
#      断言要分别报 τ_0 侧与 n 侧的跨度，别合并成一个"两轴都认识"的结论。
#   4) 正余量不够：n=0.5 臂 c 最好节点 21397.03982 Pa 比自家释放 21000 高 397.040 Pa，却**全部 851 个封闭
#      状态都 hold 1**——因为它醒在 1479，比密封破掉（1350）晚 129 步。桌面要有一根正在动的壁才有东西可 hold。
#      断言要**分两类**："margin > 0 but all flag 1"（LATE）与"margin < 0 and never asked"（AWAKE 那侧）
#      ——把它们读成一个就是删掉了桌面与打印机的区分。
#   5) 滑行同时认两条轴：n=1 滞后 0/0 步（crossing 那步 963 / 507 打印恰为 0），n=2 滞后 5/7 步
#      （704 → 709 打印 1.42139e-09、465 → 472 打印 7.28766e-09 m³/s）。同一段 x[k]*qPrev <= 0 里既无 n
#      也无 τ_0，四格四个滞后。断言要**分别报 cross 步与 sleep 步并断言 lag >= 0**，把"drive 下来了没睡"
#      写成失败就是把内核的第二半删掉了。睡眠窗 [963,1403)/[507,1530)/[709,1439)/[472,1568) 各 440/1023/730/1096
#      个恰好的 0、一个连续块；第二次放开距 K 归零 53/180/89/218 步。
#   6) 参数不醒时在打印里**不存在**：n=1/2 的 c 与 d（τ_0 2.25 对 2.5）**2001/2001** 个 token 逐字相同；
#      n=0.5 的 c/d 相同到 1479。断言要把"两档不同 τ_0 但同 2001 字符串"写成正面证据，不是当成没测到。
#   7) 身份锚（臂 b 对 5.59 的臂 a）要断言**事件身份 + 有界偏差**，不是"逐 token 相同"：旗标 2001 个状态全等、
#      阶梯 [1,318]/[1,333,507,1530]/[1,343,472,1568] 两边同整数、节点差 <= 1.000e-07 Pa（1989/1958/1663
#      个 token 逐字相同）、速率最差相对 1.780e-08。只有壁从不流动的 p 臂才真逐位相同；把它推广到流动状态会假失败。
#   8) 括号反解、sub-yield 判据、行诚实：与 5.59 同形（Sa 最差 8.500e-09 / 1.876e-06 / 8.512e-08，
#      Sb 1.968e-08 / 4.402e-07 / 8.066e-08；sub-yield 未 held 步 962/506/708/471 各一），5 格（n=0.5 d、
#      n=1/2 c/d）2000 个已决策状态全 hold 1、括号一次都没被 offered——空断言要显式说出口，别写
#      "worst 0.000e+00 at step None"。
#   9) 记账：12 条节点行只用打印 token 闭合到 4.9661e-14 / 4.8836e-14 / 4.9100e-14 m³/s（每步 9.9e-08 Pa）；
#      24 个 t=0 旗标全 0（事件扫描从 i=1 起，否则会凭空多出 24 次翻动）；任意两次决定至少隔 129 步（0.258 s，
#      无 dt 尺度抖动）。六次 app 调用合计 9.41 s（首次 4.06 s 含冷缓存）。
# 八族 146 条；ctest Test #89 单独跑记 10.61 s、`-j8` 全量里记 23.35 s（checker 直跑 10.93 s；0 内核改动），
# 全量 115/115 绿（本轮实测 933.34 s）；既有测试一字未改。
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldTauIndexRelatch/case_5.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldTauIndexRelatch/case_1.dict
./build/src/1d/apps/chflow-1d-transient tests/cases/transientYieldTauIndexRelatch/case_2.dict
python3 tests/check_transient_yield_tau_index_relatch.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/tests/cases/transientYieldTauIndexRelatch/case_5.dict"

# 一维瞬态：app 输出可解析的时间序列，check_transient.py 支持
# blowdown / rlc / slug 三种解析对照模式
./build/src/1d/apps/chflow-1d-transient tests/cases/transientBlowdown/case.dict
python3 tests/check_transient.py \
    "$PWD/build/src/1d/apps/chflow-1d-transient" \
    "$PWD/tests/cases/transientBlowdown/case.dict" blowdown

# 统一物性层：app 按 case.dict 的压力列表逐个闪蒸 + 泡点 + 制表，
# check_pvt.py 校验热力学恒等式（物料衡算、逸度相等、包络塌陷、表 vs 直接求值）
./build/src/pvt/apps/chflow-pvt tests/cases/pvtFlash/case.dict
python3 tests/check_pvt.py \
    "$PWD/build/src/pvt/apps/chflow-pvt" \
    "$PWD/tests/cases/pvtFlash/case.dict"
```

## 2. 代码规范

- 命名空间统一 `cf`；头文件放 `src/<mod>/include/cf/<sub>/<Name>.h`，
  实现放 `src/<mod>/src/...`，`#include` 一律用 `cf/...` 完整路径；
- 4 空格缩进，Allman 大括号（与现有代码一致）；
- `label`/`scalar` 只用来表达“这是 CFD 索引/量值”，不要裸用 `int/double`
  出现在公共接口上；
- 热路径避免无谓拷贝；`std::size_t` 下标转换在边界处集中处理
  （`static_cast<std::size_t>(label)`）；
- 注释只写“为什么”（物理/数值约定、踩过的坑），不复述代码。
  数值算法的关键不变量必须写明，例如 SIMPLE 的 rAU 缩放约定
  （见 `docs/architecture.md` §3.5）与 BiCGStab 的预处理约定——
  这些是回归历史换来的结论，回退它们代价极高。

## 3. 测试准则

- **每个新算子/模型必须带数值校验**，而不仅是“跑得通”：
  - 有解析解优先对比解析解（现有三个 verify 算例都是这种）；
  - 没有解析解的，用守恒量（通量守恒、散度定理 `Σ Sf = 0`）
    或 manufactured solution；
- 新增算例的目录约定：`tests/cases/<name>/case.dict` +
  `tests/check_<name>.py`（断言阈值、打印实测误差），并在
  `tests/CMakeLists.txt` 里 `add_test` 注册；
- **阈值分两类，不许混淆**：恒等式类（守恒量、壁面链闭合、离散精确解）用
  机器精度级的严格阈值；预测类（摩阻律、湍流模型给出的 `U_b`）按网格自身能
  支撑的精度放宽，并把**实测数字写进 case.dict 头注释**——`turbulentChannel`
  的 8 个对数层单元给出 −14.5%/+5.4% 的 `U_b` 偏差，这是网格的精度，不是
  模型的；把这种带宽写成 1% 只会让回归在下一次改网格时假失败；
- checker 需要**多份配置同一算例**时（例如两个湍流模型），把 case 目录复制
  到临时目录再改写，而不是往仓库里加第二份 case 或在仓库目录里留产物
  （求解器把 `.vtu` 写在 cwd）；
- 模型层要有 C++ 单元测试：断言闭合的**代数不变量**（如平衡壁面单元
  `k = u_τ²/√C_μ` 配上 `ε = u_τ³/(κy)` 必然给出 `nut = κu_τy`），而不只是
  "跑完没报错"；
- 求解器改动后必须全量 `ctest` 绿；**提交/推送的默认前提就是全量绿**，
  失败状态不入库；
- 临时调试脚手架（环境变量开关、`fprintf` 等）在合并前删除。

## 4. 算例字典格式

```
block { min (0 0 0); max (4 0.2 0.1); cells (100 20 1);
        xmin { name inlet; type inlet; } ... }
nu 0.01;
initial { U (8 0 0); }                     // 可选：起步速度场
patches { inlet { U (0.1 0 0); UType fixedValue; pType zeroGradient; } ... }
schemes { ddt steady; div upwind; laplacian corrected; grad corrected; }
turbulence { model kOmegaSST;              // laminar | kEpsilon | kOmegaSST
             intensity 0.05; lengthScale 0.0056; relaxation 0.7; }
simple { maxIterations 300; continuityTol 1e-6; pcorrTol 1e-4;
         momentumTol 1e-4; relaxU 0.9; relaxP 0.03; pRelaxStart 0.02;
         force (19.5 0 0); }
piso { dt 0.004; endTime 3; correctors 2; continuityTol 1e-8;
       reportEvery 250; force (19.5 0 0); }
output { vtk out_poiseuille.vtu; }
```

`block{}` 的 patch `type` 参与物理：`wall` 让湍流模型建立壁面层次与 `y+` 报告，
`inlet` 才是 `turbulence{}` 里 intensity/lengthScale 的注入面；缺 `turbulence{}`
即 laminar。体积力驱动（`piso{force}`）用 `symmetry` 封住流向两端做单列渠道，
`p` 只需一个 `fixedValue` 锚点。

`simple{momentumTol}` 是动量固定点变化 `max|U − Uold|/max|U|` 的门槛，两个
收敛分支都要过它：单列渠道这类流向均匀的算例所有 `phi` 恒为 0，只靠连续性
会在第 1 次外迭代就宣布收敛。`simple{force}` 与 `piso{force}` 同义。

VOF 算例（`chflow-vof`）复用同一套 `block{}`/`mesh{}`、`schemes{}`、`output{}`，
再加三个自己的块——速度是**给定的解析场**而不是求解出来的，验收的是输运与界面几何：

```
velocity { type uniform;  U (1 0 0); }   // 或 type rotation; omega 1.0; center (0 0 0);
                                         // 后者是 U = omega zhat x (x - center)
alpha { plane { point (0.2 0.05 0.05); normal (1 0 0); } }   // 单元填充的阶跃
//     或 alpha { circle { center (0.2 0 0); radius 0.16; smooth 4.0; } }
//        smooth > 0 把圆界面改成 tanh 抹平，厚度以 dx 计（曲率测量需要，见下）
patches { inlet  { alphaType fixedValue; alpha 1; }   // 入口必须喂相分数
          outlet { alphaType zeroGradient; } }
time    { dt 0.0125; endTime 0.2; reportEvery 4; interfaceThreshold 0.1; }
schemes { ddt Euler; div upwind; grad corrected; }    // ddt 只接受 Euler
solverTol 1e-13;                                     // solverMaxIter 缺省 2000
output  { vtk out.vtu; vtkInitial initial.vtu; }      // 后者写**初始**界面的 n̂/κ
```

`velocity{}` 的两种都是坐标的线性函数，所以逐单元 `sum_f U(x_f)·Sf = ∫div U dV`
精确为零——守恒与 α ∈ [0,1] 是隐式 Euler + 迎风作用在离散无散 phi 上的 M-矩阵
性质，不需要限制器，也不必把结果裁剪回去。`interfaceThreshold` 是界面带判据相对
`|∇alpha|` 全场最大值的比例（法向只在带内单元非零，曲率模板因此只读带内面）。
`vtkInitial` 存在的原因：对流会沿流向各向异性地抹开界面，**逐单元**法向随之倾斜，
所以逐单元几何断言只作用在解析的初始场上（带均值 κ 仍然守恒，见
`docs/architecture.md` §5.10）。

**动量模式**：用 `surfaceTension{}` 替换 `velocity{}`（两者必须恰好给一个），α
冻结、PISO 读出 CSF 隐含的压力-速度响应；`patches{}` 在同一块里加 `U (0 0 0)`
（无滑移）与一块 `p 0` 的 patch（闭合盒的压强常数靠它锚定）：

```
surfaceTension { sigma 0.08; rho 1; nu 0.01; correctors 2; continuityTol 1e-12; }
patches { left { U (0 0 0); p 0; } right { U (0 0 0); pType zeroGradient; } ... }
```

验收的是拉普拉斯跃变 `p_in − p_out = σ/(ρR)`（内核的 p 是运动学压强 p/ρ）、
斜坡只属于带、解场的镜像对称、以及伪速度幅值对 √(σ/(ρR)) 的比值。

**密度基可压缩 Euler**（`chflow-3d-compressible`）不读 `schemes{}`/`patches{}`：状态是每单元
的广延守恒量 `(ρV, ρu V, ρEV)`，边界一律是不可透滑移壁，所以字典只有四块——
`fluid{gamma}`、`initial{diaphragm; left/right{rho,u,p}}`、`time{dt,endTime,reportEvery}`、
`output{vtk}`（网格仍走下面的 `block{}`/`mesh{}` 二选一）：

```
fluid { gamma 1.4; }
initial {
    diaphragm 0.5;
    left  { rho 1;     u 0; p 1;   }
    right { rho 0.125; u 0; p 0.1; }
}
time { dt 0.00025; endTime 0.2; reportEvery 100; }
```

app 每 `reportEvery` 步报 `mass / mom / energy / impulse / pmin / pmax / cfl`：闭盒里质量与
总能量应当逐报告点不变，动量随 `impulse(t) = ∫dt Σ_b (p Sf)_x` 变化，所以
`check_euler.py` 既断言 `mom + impulse = mom(0)`，也用 Python 重放同一套 HLLC 递推与写出的场
逐格比对（1e-9）。CFL 越界（`dt·amax/dxMin > 1`）与非有限状态都直接报错，不做回退。

**三维区域分解**（`chflow-3d-parallel`）不读 `schemes{}`/`patches{}`：它把网格切成
1、2、3、5…个域，逐域跑显式 Gauss 算子，所以字典里只有网格、域数、两个测试多项式场
和验收阈值。`field{}` 的两个键都是十个系数
`T = c0 + c1x + c2y + c3z + c4x² + c5y² + c6z² + c7xy + c8yz + c9zx`——
`linear` 把 c4..c9 置零，是**精确性锚点**（有界线性插值对线性场精确，Gauss 和给出解析梯度）；
`quadratic` 是**不变性场**（它的离散梯度本就不是解析值，欠它的只是"分解 == 串行"）：

```
block { min (0 0 0); max (1 0.5 0.25); cells (12 8 5); }
decompose { domains (1 2 3 5 7 8 16); }        // 每个域数各切一遍、各写一份 vtu
field {
    linear    (1 2 -1 0.5 0 0 0 0 0 0);
    quadratic (1 2 -1 0.5 3 -0.25 0 1.5 -2 0.75);
}
verify {
    linearGradAnalyticTol    1e-10;   // 线性场对解析梯度（缺该键则只打印不断言）
    quadraticGradAnalyticTol 0.10;
    stencilDevTol            0;       // 分解算子 vs 串行：要求逐位相等
    domainCountDevTol        0;       // 1 域与 16 域写出的场：要求逐位相等
    divInteriorTol           1e-13;   // 无散通量在被切缝穿过的内区格上的净出流
}
output { vtk out_parallel.vtu; }      // 实际写 out_parallel_<域数>.vtu
```

app 报 `parallel stencil …` / `parallel flux …` / `parallel reduce …` 三族行（17 位有效
数字）；`verify{}` 只有 `check_parallel.py` 读，app 自己不判定。求解器把 `.vtu` 写在 **cwd**，
所以手动跑要在 case 目录里跑，或者接受根目录多出一份
`out_parallel_*.vtu`（已被 `.gitignore` 覆盖）；`check_parallel.py` 自己以 case 目录为 cwd。

**分区 1D ↔ 3D 耦合**（`chflow-coupled-1d3d`）把一个稳态管网字典与一个三维块字典叠在
同一个文件里，再加一块 `couple{}` 说明哪条边被详细化——1D 侧读 `fluid/nodes/sources/pipes/solve`，
3D 侧只在 `detail{}` 里读自己的网格、格式与 SIMPLE 控制：

```
fluid { rho 800; mu 0.1; }                    // 两侧共用；3D 侧的 nu = mu/rho
nodes { Src {}; J1 {}; J2 {}; Snk { pressure 1.0e6; } }
sources { feed { node Src; rate 2.4e-4; } }   // 速率驱动：分支分配是结果而不是输入
pipes {
    A  { from Src; to J1; length 5.0; diameter 0.1; }
    X  { from J1; to J2; length 0.5; diameter 0.02;   // length/diameter 只用于
         coupledResistance 3.75e6; }                  // 报告 v 与 Re，压降来自 3D
    By { from J1; to J2; length 1.5; diameter 0.03; }
    B  { from J2; to Snk; length 5.0; diameter 0.1; }
}
solve { maxIter 200; tol 1e-10; relax 1.0; }  // 网络节点残差的地板是注入速率的 ~1e-12 相对
couple {
    edge X;                                  // 被详细化的管名
    seeds ( 1.0e6 3.75e6 1.2e7 );            // 每个初值电阻各跑一整趟外层回路
    relax 1.0;                               // < 1 给非线性（湍流）伙伴阻尼
    tolerance 1e-4;                          // 接口相对残差；地板是伙伴自己解的误差
    maxOuterIterations 6;
    accelerate no;                           // yes ⇒ 界面上用准 Newton（Aitken/IQN）取系数
    accelerates ( 0 1 1 );                   // 可选：按 seeds 顺序逐组覆盖；长度 1 时摊平
}
detail {
    block { min (0 0 0); max (0.5 0.02 0.02); cells (50 20 1);
            xmin { name inlet;  type inlet;  }    // 进出口由 patch **类型**识别，
            xmax { name outlet; type outlet; }    // 不需要在外层 patches{} 再写一遍
            ymin { name bottom; type wall;   }    // 无滑移
            ymax { name top;    type wall;   }
            zmin { name front;  type symmetry; }  // 侧向：内部按 zeroGradient 处理
            zmax { name back;   type symmetry; } }
    probe { min 0.12; max 0.49; }            // 测已发展梯度的横截面窗口，须在入口段之外
    turbulence { model kOmegaSST; intensity 0.05; lengthScale 1.4e-3; }
                                             // 缺省 = laminar：块解的就是它一直解的那条方程
    schemes { ddt steady; div upwind; laplacian corrected; grad corrected; }
    simple { maxIterations 300; continuityTol 1e-6; pcorrTol 1e-5; momentumTol 1e-5;
             relaxU 0.9; relaxP 0.03; pRelaxStart 0.02; }
}
output { vtk out_coupled_block.vtu; }        // 最后一个种子的块场
```

三条容易踩的约定：① `coupledResistance > 0` 与黑油/PVT 管、`valve`、`pump` 互斥，
`loadNetworkFromDict` 直接报错——混合语义下"这条边的压降"不再是单一数。
② 块内**不能**给侧向 `symmetry`/`wall` 的动量条件：SIMPLE 把三个动量分量装配在同一
矩阵、共享对角，而这两个类型只固定法向分量，第一次装配就撞对角一致性断言；
`DetailBlock::applyBoundary` 因此对非 inlet/outlet/wall 面一律用 `zeroGradient`。
③ 三维核的 `p` 是运动压力（m²/s²），管网用 Pa，换算在 `lossAt` 里乘 ρ 完成——
`BlockMeasure::dpPa` 已经是 Pa，别在 checker 里再乘一次。
④ 挂 `turbulence{}` 时块的入口湍流强度由**界面第一趟送来的速率**定：`lossAt` 先
`applyBoundary(qMag)` 再 `RansModel::create(...)`，模型读 inlet patch 上的 `fixedValue` 速度
算出 `k0 = 1.5(U·I)²`，此后暖启动。代价是种子电阻给得离谱时 `k0` 也离谱，接口要靠几趟把
速率找回来（实测 `seeds 3.0e6` 时第一趟落在 2.28e-4，正确值的 0.44 倍）。另两条硬约束：
RANS 块的相对连续性残差停在 ~2e-3（修正后的压力残差到 1e-10，但通量带着不平衡），
所以 `couple.tolerance` 不能严于这个地板；网格第一格 `y+` 须在壁函数带（30–300）内，
否则块报的压降不再是壁律的压降。`tests/cases/coupledTurbulent` 就是按这两条配的。

**瞬态分区耦合**（`chflow-coupled-1d3d-transient`）把上面字典里的 `solve{}` 换成 1D 瞬态的
`time{}`，其余键（`fluid/nodes/sources/pipes/couple/detail/output`）与稳态侧逐字同形；变的
是这几处的语义：

```
nodes {
    Src { volume 1.0; liquidFraction 0.5; pressure 1.0e6; }  // 接住突加速率的节点要留持液余量
    Snk { fixedPressure 1.0e6; }                             // 瞬态侧的定压键名
}
pipes { X { ...; coupledResistance 3.75e6; } }               // 只替换本管**摩阻**
time  { dt 0.02; endTime 1.2; reportEvery 1; }              // dt/τ ≈ 0.077：一步跨 7.7% 时间常数
couple {
    edge X;  seeds ( 3.75e6 );  relax 1.0;
    tolerance 1e-4;               // 不得严于伙伴自己的分辨率：pcorrTol 1e-5 ⇒ 残差地板 ~1e-5…1e-4
    maxSubIterations 8;           // 取代稳态的 maxOuterIterations，用尽即 fatal
    accelerate no;                // yes ⇒ 每步的界面历史重建后用准 Newton（见 architecture.md §5.21）
    inertias ( network partner ); // 每个种子臂一个区域惯量来源（缺省全是 network，见下面 ⑥）
}
```

① `coupledResistance` 在瞬态侧**只**替换本管摩阻：重位和该区域的惯性 `ρL/A` 仍留在 1D 的
隐式动量行里，与网络其余部分一起被隐式解，伙伴每步只报"此刻该速率抵抗的压降"（准稳态
`lossAt`）。这不是省事：让块报含惯性的压降，接口映射就成了 `R → a − c/R`，其在本算例 stiff
区的线性化导数 `|c/R²| = 2.7e8 ≫ 1`，分区 Gauss–Seidel 直接发散；反过来把区域惯量算两遍，
时间常数 `τ = (M_X + M_By)/(R_X + R_By)` 会被抬高 `M_X/(M_X+M_By)`——本算例 43%。checker 因此
是拿整条爬升回归出的 τ（带 5%）来钉这个切分，而不是只验末端。
② 耦合管与 `valve`/`valveClose`/`valveOpen` 互斥，构造网络与 `setCoupledResistance` 刷新时
都会拒绝（`requireNoValve`）。接住注入速率的自由节点要同时给 `volume`（否则"free node needs a
volume"）与 `liquidFraction < 1`（缺省是 1.0，第一步突加的速率还没有任何管流带走它，提交时即
撞"liquid exceeds its volume"）。定压节点写成 `fixedPressure`。
③ 步内回路是"解而不采纳"：`net.solveStep(dt)` 把状态挂在 pending 上、`commitStep()` 才推进
时钟与持液输运，所以刷新系数走 `setCoupledResistance` 会**显式作废**挂起的步（那步是用旧
系数解的，已不是系统状态）。谈不拢就 fatal，不会提交一个两侧不一致的状态往下滚。
④ 瞬态接口每步都**热启动**块，而 SIMPLE 外收敛的 `settled` 判据若写回严格单调
（`contErr <= bestContErr`）会让首个 predictor 连续性低于该 case 最终落到的
HbyA/压力通量地板的热启动再也走不到 `|dp|` 分支（实测 0.00768788 对 0.00769587，首步烧满
300 迭代报未收敛）；`simpleSolver.cpp:343` 因此留 1% 余量，防发散的初衷不变。
⑤ `detail{}` 在瞬态侧同样能挂 `turbulence{}`（键与稳态逐字同形，稳态 ④ 的两条硬约束——
通量地板与 `y+` 带内——一字不改地带过来）。挂上之后变的不是接口协议而是**分工的动力学**：
伙伴的 `dp/q` 随速率变，分流的时间常数由两条支各自的**割线斜率** `m = dln(dp)/dln(q)` 定，
不是由电阻弦定：两条支相减消掉节点压降后 `M·dq/dt = f_By(Q−q) − f_X(q)`，线性化的恢复系数是
`S = m_X·R_X + m_By·R_By`（推导见 architecture.md §5.23）。稳态伙伴 `m ≡ 1`、弦即斜率，所以
`transientCoupled1d3d` 的 τ 用 `S_chord = R_X + R_By` 就能对上；两个伙伴都是湍流时
`m ≈ 1.7`，弦形式把 τ 读高 +62%。checker 因此拿两条路各算一遍（`tests/cases/
transientCoupledTurbulent`：斜率形式 −0.34%、弦形式 +61.7%），并且**不能**再拿 1D 侧自己的
`frictionDp` 反证旁路——那等于让网络跟自己核对，这里改用 Churchill 1977 的公开式子在
checker 里重写一份作外部锚（旁路动量行 7e-14）。逐步的界面残差实测落在 4.4e-5…2.0e-4
（最大在第 23 步，正压在 `tolerance 2e-4` 下面）——这个量级来自伙伴自己的通量不平衡，不是
调出来的：稳态 ④ 说的两条地板在时间步进里一条不减，反而每步只剩 `maxSubIterations` 趟去够它。
⑥ 区域的惯量有两个来源：`couple{ inertias ( network partner ) }` 按 `seeds` 的臂序给（长度 1
摊平到全部，缺省 `network`；多边算例同一个参数按 `starts` 的臂序给，两臂同台的实测见 §5.26）。`network` 用耦合管自己的 `rho*L/A`——它只和 case 里敲进去的名义
孔径一样好（方渠道配内切圆多报 **27.3%**）；`partner` 让每条耦合边取它的块在区域网格上量的
`rho*V/Area^2`。两条别改的约定：**惯量只装进 1D 的动量行，不写进界面系数**——伙伴若报
`R(q)*q + M*(q-qPrev)/(dt*q)`，映射就成了 `R -> a - c/R`，而"接口能分辨这条爬升"的条件正是
`dt < M/R`，于是任何减速步把系数推过零：实测 fatal 于 `Z = -2180.209163 Pa`，把 fatal 放宽后块
烧满 300 趟 SIMPLE 钉在连续性 0.00522332；**装载发生在全部边校验完之后**，且测量只取一次
（它是网格属性、与时间层无关），所以报不出区域的伙伴不会把网络留在半声明状态。截面在变的区域
直接 fatal（平面面积差过 1e-9 相对就没有单一惯量可报），而不是报一个错数。日志里每步的 `M` 是动量行
**实际带**的那个数（`EdgeExchange::inertance`，不在解时间层时为 0），块摘要的 `inertance` 是块
自测的那个——`tests/check_transient_coupled_inertia.py` 就是拿这两行分别对撞 `rho*L/A` 与
`rho*V/A^2` 的独立算式来验归属的（推导与实测表见 architecture.md §5.24）。
同一个决定还可以**逐边**给：`couple{ inertiaFrom{ <arm>{ <edge> network|partner; } } }` 是一张按
`starts`/`edges` 双重索引的表，缺席的边退回该臂的 `inertias` 那一项，于是**一次运行**里一条边的动量
行吃块测的 `ρV/A²`、另一条吃图纸的 `ρL/A`（表长不等于边数、或某条边点名了块却没区域可报，都在装载
之前 fatal；实测与逐边可分性见 architecture.md §5.29）。

**多边分区耦合**（`chflow-coupled-1d3d-multi`）把 `edge`/`seeds`/`detail` 三个单数键换成
复数：`couple{ edges (…) starts{} }` 列出被详细化的所有管名与"一整组初值"，`details{}` 按
**管名**给每块自己的网格、格式、SIMPLE 控制与 vtu。1D 侧的 `fluid/nodes/sources/pipes/solve`
与稳态单边逐字同形，一条边被详细化就把该管的摩阻换成它的界面电阻：

```
pipes {
    X1 { from J1; to J2; length 0.5; diameter 0.02; coupledResistance 3.75e6; }  // 并联支里
    By { from J1; to J2; length 1.5; diameter 0.03; }
    X2 { from J2; to J3; length 0.4; diameter 0.015; coupledResistance 1.6e7; }  // 串联位置上
}
couple {
    edges ( X1 X2 );            // 顺序即接口里的边顺序，日志与结果都按它排
    relax 1.0;  tolerance 1e-4;  maxOuterIterations 8;
    accelerate no;              // yes ⇒ 界面系数由准 Newton 给（architecture.md §5.21）
    accelerates ( 0 1 1 );      // 可选：按 starts 的键序逐组覆盖，长度 1 时摊平到全部
    starts {                    // 每组给**每条边**一个初值电阻，整组各跑一趟外层回路
        lo   { X1 1.25e6;  X2 5.3e6; }
        mix  { X1 3.75e6;  X2 5.3e6; }
        high { X1 1.125e7; X2 5.3e6; }
    }
}
details {
    X1 { block { … } probe { min 0.12; max 0.49; } schemes { … } simple { … }
         vtk out_multi_X1.vtu; }
    X2 { block { … } probe { min 0.15; max 0.38; } schemes { … } simple { … }
         vtk out_multi_X2.vtu; }
}
```

多边侧六条容易踩的约定：① 一层循环是**块 Gauss–Seidel**，不是每边各自迭代到收敛：一趟 =
一次网络解 + 每个伙伴各一次解，所有边的系数同一趟刷新，停止条件取**最差**的那条边
（`res.residual == max_i residual_i`）。逐边各自收敛会让先谈拢的边读到别的边的半状态，
"一趟对应一个系统状态"就没了；checker 因此数日志里 `SIMPLE iter 1:` 的条数来钉
"块解 == 趟数 × 边数"。② 接口结果里**没有** `dpNetwork` 字段——它恒等于 `resistance * q`，
两侧在 `CoupledNetwork.h` 里各自算同一式子；`EdgeExchange` 每条边一行，`OuterIteration.edges`
是同一趟所有边的快照。③ `starts{}` 用嵌套字典而不是列表：一组初值要给每条边各一个数，
而 `cf::common::Dictionary` 的列表只支持标量，读不了 `((a b) (c d))`；键名（`lo/mix/high`）
只用于日志与 checker 配对。④ 每边的 vtu **只在最后一组初值**写（每组的块解会覆盖同一文件），
`out_multi_*.vtu` 已被 `.gitignore` 覆盖。⑤ 伙伴非线性时接口映射会把外层的收缩率压到接近 −1：
并联边上的二次伙伴把电阻映射成 `R2 → R0 + k·q2`，其中 `q2 = Q·R1/(R1+R2)` 是分流结果，映射自身
的不动点处导数是 `−k·Q·R1/(R1+R2)²`。它**有界**，永远不到 −1，所以纯 Gauss–Seidel 不发散、只是
任意慢：`k=1e7` 是 −0.041（单元用例用的就是它，`tolerance 1e-8` 下 30 趟上限内落定），
`k=1.75e10` 是 −0.875，实测单边要 **183 趟**、两边一线性一二次要 **211 趟**才到 `1e-10`，
瞬态同一问题每步只给 8 趟时未加速臂直接谈不拢（fatal）。`accelerate yes` 走界面准 Newton
（Aitken/IQN，architecture.md §5.21）把这三笔账分别压到 **11 / 9 趟**与 20 步全过（67 趟对
122）；线性伙伴上它是无害的——两臂都在第 2 趟落定，首趟系数逐位相同。加速的收益只在
`coupling.unit` 的闭式伙伴上被证明；随仓库发的多边算例两块都是层流伙伴，加速臂只作回归
（同样的 8 趟 / 16 次块解，落到与未加速臂同一个工况点）。⑥ 开 `accelerate` 时有两条别改的约定：
割线基要**从最新的 pair 往回做 Gram–Schmidt** 重建，而不是"只追加、正交化时才决定留不留"——
单边接口所有 `dQ` 平行（秩 1），追加式会把窗口冻在第一个 pair 上（实测 14 趟 vs 11）；瞬态则
**每个时间层重建历史**（`qn->reset()`），因为接口问题含本步的惯性项，上一步的 pair 量的是另一个映射。

**多边 + 瞬态**（`chflow-coupled-1d3d-transient-multi`）把上面两套键原样搬进时间步进：`solve{}` 换成
`time{ dt endTime reportEvery }`，`maxOuterIterations` 换成 `maxSubIterations`（它现在是**每个时间层**
的趟数上限，不是整趟运行的），`accelerate`/`accelerates`/`inertias`/`inertiaFrom` 都照旧透传——`inertias`
因此**从一开始就是逐边生效的**（归属循环先量完所有边、全过校验才装），两臂只差来源的实测见 §5.26；
但它是每次运行（臂）一个整体选择，所以"这条边用孔径、那条边用区域"在同一次运行里要交给
`couple{ inertiaFrom{ <arm>{ <edge> network|partner; } } }`——按臂、按边的表，实测见 §5.29。
四条容易写错的约定：
① 一组 `starts{}` 必须给 `edges ()` 里**每条边**一个数，缺一个直接 fatal——半个网络被详细化等于
两侧描述的不是同一个系统；② 采纳是全网一起的：全部边谈拢才 `net.commitStep()`，任何一条没谈拢就
fatal 并点名最差那条边，所以**不能**把某条边的预算单独放宽；③ 每步的 `passes` 是 1 到
`maxSubIterations`，`subPasses` 与逐趟行的条数一致，checker 靠这些行重算调度，改日志格式会同时打断
"停在最差边"和"块解 == 趟数 × 边数"两条断言；④ 成本按**每层重谈一次**算，不是按边数——想压时间优先
动这三处：跨对称方向的格数（两侧都是 `symmetry` 的直渠道是 2D 问题，本算例 X2 从 6 格降到 1 格，
弦阻只移动 6.7e-5 相对）、步数（`dt` 与 τ 之比决定每常数几步，检查器若按离散闭式预测整条轨迹就不必
等爬升衰减完）、初值组数。**别**用"提高 `tolerance`"来省时间：接受趟的残差地板由伙伴自己的分辨率给
（本算例最差 8.531e-5 对 `tolerance 1e-4`），问得比块答得细只会多烧 SIMPLE 扫。

三维 app 的网格二选一（同时给出或都不给都会直接报错）：结构化用上面的
`block{}`；非结构用 Gmsh 导出的 **ASCII** `.msh`（2.2 或 4.1，一阶
tet/hex/prism，边界表面必须打好 physical group，否则只能按包围框平面命名）：

```
mesh {
    format gmsh;
    file "../../meshes/prismChannel_22.msh";   # 相对 cwd，找不到时按 .dict 所在目录解析
    patches { inlet { type inlet; } outlet { type outlet; } }   # 缺省类型 wall
}
```

`patches{}` 里的名字必须与 Gmsh physical group 同名；场边界条件仍写在外层
`patches{}`（两者是不同字典）。网格夹具放在 `tests/meshes/`，
`tests/unit/gmshMeshTests.cpp` 通过 `CF_MESH_DIR` 读取它们。

格式刻意兼容 OpenFOAM 字典语法（`key value;`、嵌套 block、`(a b c)` 列表），
解析器在 `src/common/src/Dictionary.cpp`。

## 5. 提交与分支

- 小步提交，commit message 说明动机（为什么改），而不是罗列文件；
- 只 `git add` 明确相关的文件路径，不要 `git add -A`（避免带入
  生成的 `out_*.vtu`、`build/`、临时脚本）；
- 推送前提：`ctest` 全绿；
- `.gitignore` 已覆盖构建目录、对象文件、Python 缓存、IDE 目录和生成
  的 VTU 结果；新增产物类型请补充 `.gitignore`，而不是提交产物。

## 6. 与 OpenFOAM-dev 的关系

参考 OpenFOAM-dev 只到**设计思想**层面（API 形状、模块划分、算法回路）。
禁止复制 OpenFOAM 源码或其 GPL 衍生代码；如借鉴论文算法（如 BiCGStab、
SIMPLE），在代码注释或文档里引用出处即可。
