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
ctest --test-dir build --output-on-failure -j8  # 全量回归（93 项）约 209 s ≈ 3.5 分钟，
                                             # 墙钟就是最慢那一条：稳态 profile 级映射 209 s；
                                             # 其余大头（瞬态惯量归属 89 s、VOF 振荡 73 s、
                                             # 形状 x 湍流 x 瞬态 70 s、稳态多边 63 s）
                                             # 在 -j8 下被并行掩盖
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
