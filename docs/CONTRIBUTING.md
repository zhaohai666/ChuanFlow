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
ctest --test-dir build --output-on-failure   # 全量回归（poiseuille 算例约 6 分钟）
```

单个验证算例可以直接手动跑，便于调试：

```bash
./build/src/3d/apps/chflow-simple-foam tests/cases/poiseuilleChannel/case.dict
python3 tests/check_poiseuille.py \
    "$PWD/build/src/3d/apps/chflow-simple-foam" \
    "$PWD/tests/cases/poiseuilleChannel/case.dict"

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
- 求解器改动后必须全量 `ctest` 绿；**提交/推送的默认前提就是全量绿**，
  失败状态不入库；
- 临时调试脚手架（环境变量开关、`fprintf` 等）在合并前删除。

## 4. 算例字典格式

```
block { min (0 0 0); max (4 0.2 0.1); cells (100 20 1);
        xmin { name inlet; type inlet; } ... }
nu 0.01;
patches { inlet { U (0.1 0 0); UType fixedValue; pType zeroGradient; } ... }
schemes { ddt steady; div upwind; laplacian corrected; grad corrected; }
simple { maxIterations 1600; continuityTol 1e-6; pcorrTol 1.4e-4;
         relaxU 0.9; relaxP 0.8; }
output { vtk out_poiseuille.vtu; }
```

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
