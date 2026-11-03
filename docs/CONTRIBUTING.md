# 贡献指南（CONTRIBUTING）

## 1. 开发环境与构建

- C++20（Apple clang 14+ / GCC 12+ / Clang 14+ 实测通过），CMake ≥ 3.24；
- 内核层（`src/common`、`src/1d`、`src/3d`）**零第三方依赖**：只允许 C++ 标准库。
  Eigen/PETSc/spdlog 等只能出现在未来的可选后端模块，不进入内核；
- 测试需要 Python 3.10+（只用标准库）；
- 可选层：Python 绑定 `src/bindings` 经 FetchContent 拉取 pybind11
  （首次 configure 需要网络；离线或缺 Python 开发头文件时用
  `-DCHFLOW_WITH_PYTHON=OFF` 跳过，C++ 构建不受影响）；CLI `cli/` 为
  零依赖 Rust，未安装 cargo 时对应 CTest 用例自动不注册。

```bash
cmake -S . -B build
cmake --build build -j8
ctest --test-dir build --output-on-failure   # 全量回归约 1 分钟（poiseuille 占 37 s）
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
simple { maxIterations 1600; continuityTol 1e-6; pcorrTol 1.4e-4;
         momentumTol 1e-4; relaxU 0.9; relaxP 0.8; force (19.5 0 0); }
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
