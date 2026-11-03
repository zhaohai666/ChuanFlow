# ChuanFlow
> 国产开源一体化多相流仿真平台
> 中文名称：川流
>
> 🎯 定位：集成**一维稳态管网生产系统（对标 Pipesim）**、**一维瞬态多相流动保障（对标 OLGA，段塞流、水合物、启停瞬态）**、**三维有限体积CFD（对标 Fluent）**，支持 1D ↔ 3D 耦合仿真，面向油气集输、井筒生产、能源装备多相流场景。
>
> > ⚠️ Important Notice
> > ChuanFlow **不是 OpenFOAM 的 Fork**。内核设计思路、有限体积框架、场/网格抽象参考 OpenFOAM-dev 的工程思想，但求解器、数据结构、离散格式、多相模型、一维管网模块全部独立自研，代码版权归属 ChuanFlow 项目自身，遵循独立开源协议。

[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)]()
[![Build Status](https://img.shields.io/badge/build-WIP-orange)]()

## ✨ 核心能力
- **1D Steady 稳态管网求解**：节点分析、水力热力计算、井筒/集输管网压力温度分布，对标 Pipesim
- **1D Transient 瞬态多相流动保障**：双流体模型、段塞追踪、持液率演化、水合物动力学、管线启停、泄压/充压瞬态工况，对标 OLGA
- **3D FVM 三维CFD求解器**：基于有限体积法，支持单相/多相流、湍流、换热；网格体系、场与离散框架参考 OpenFOAM-dev 思想，自研实现
- **统一 PVT 物性库**：一维、三维求解器共用同一套闪蒸、相平衡、流体状态方程，保证耦合仿真物性一致
- **1D ↔ 3D 耦合接口**：管网瞬态结果作为三维CFD进出口边界条件；三维局部装备计算回传给管网系统
- **多入口使用形态**：C++ 原生库 / Python API / CLI 命令行 / Tauri 桌面GUI，集群自动化、本地交互都支持
- **标准IO兼容**：支持 Gmsh msh、VTK/VTU、HDF5，可导入 OpenFOAM 网格，结果直接用 ParaView 可视化

## 📁 项目目录结构
```
ChuanFlow/
├── src/
│   ├── 1d/          # 一维管网、稳态/瞬态多相求解内核（C++20）
│   ├── 3d/          # 三维CFD有限体积求解内核（C++20，架构参考OpenFOAM-dev）
│   ├── pvt/         # 共用流体物性库：相平衡、闪蒸、PVT表、水合物模型
│   ├── bindings/    # pybind11 绑定层：暴露C++内核到Python
│   └── common/      # 公共基础组件：日志、单元、向量、稀疏工具、配置解析
├── cli/             # Rust 实现 chuanflow-cli：批量算例、集群调度、格式转换
├── gui/             # Tauri + egui 桌面前端：管网拖拽建模、曲线监控、算例管理
├── python/          # Python SDK、示例脚本、后处理工具、参数扫描工具
├── docs/            # 文档：理论手册、用户手册、开发手册、模型说明
└── tests/           # 测试算例集：稳态管网、段塞基准、CFD对流扩散、多相流校验算例
```

## 🧱 系统架构设计
ChuanFlow 采用 **内核与上层应用解耦、求解与UI进程隔离** 的分层架构，核心设计原则：
1. **求解内核纯库化、无UI依赖**：`libchuanflow-1d.so` / `libchuanflow-3d.so` 是无状态计算库，GUI 崩溃不影响正在运行的仿真任务；
2. **物性统一底座**：`pvt` 模块是全局唯一的流体计算层，一维、三维、耦合仿真共用，消除物性不一致带来的耦合误差；
3. **接口标准化**：所有求解器对外暴露统一的 case 输入（YAML）与结果输出（VTK/HDF5）；
4. **渐进式并行**：CPU 共享内存 + MPI 分布式并行为第一优先级；Kokkos 异构（NVIDIA GPU / AMD GPU）作为后续扩展。

### 分层说明
#### 1. 底层公共层 `src/common`
基础基础设施，所有模块依赖：
- 日志系统（spdlog）、格式化（fmt）
- 基础容器、几何向量、单位系统
- 配置解析、异常体系、时间步控制器
- 通用稀疏矩阵包装（Eigen / PETSc 适配层）

#### 2. 统一物性层 `src/pvt`
整个项目最重要的共享模块：
- 纯组分、黑油、组分闪蒸计算
- 流体粘度、密度、导热系数、表面张力
- 水合物生成/分解动力学模型
- PVT 数据表插值
> 一维管网持液率、压力损失、三维CFD相界面计算全部调用此模块，保证耦合仿真一致性。
>
> 首版已落地：纯组分与多组分的 Peng-Robinson 状态方程（密度、压缩因子、逸度系数）、
> 多组分等温 TP 闪蒸、混合物泡点、PVT 双线性表插值，以及 `chflow-pvt` 应用与
> `chflow.tp_flash` Python 接口；黑油分离器、黏度/导热/表面张力关联式、
> 水合物相平衡与生成条件仍待实现。

#### 3. 一维求解层 `src/1d`
独立子系统，分为稳态求解器与瞬态求解器：
- **SteadySolver**：管网节点法、环路迭代、水力热力耦合，对标 Pipesim；
- **TransientSolver**：一维双流体模型、段塞捕获、动量/质量/能量守恒方程组，适配 OLGA 典型流动保障工况；
- 管线单元库：水平管、倾斜管、井筒、阀门、分离器、节流装置；
- 边界条件：定压、定流量、井口生产曲线、瞬态压力边界。

#### 4. 三维CFD求解层 `src/3d`
> 框架思想参考 OpenFOAM-dev 的有限体积体系（网格拓扑、场存储、离散算子、时间推进），全部代码独立自研，不复用 OpenFOAM 源码。
- 网格模块：多面体网格、面/单元拓扑、域分解（MPI并行）；
- 场系统：标量场、矢量场、张量场；
- 离散算子：梯度、散度、拉普拉斯、插值；
- 时间推进：稳态SIMPLE、瞬态PISO；
- 多相模型：VOF、Euler-Euler双流体；
- 湍流模型：k-ε、k-ω SST；
- 线性求解后端：小规模 Eigen，大规模 PETSc。

> 首版已落地：多面体网格（blockMesh + Gmsh msh）、场与 BC、fvc/fvm 离散、
> 定常 SIMPLE（全压力回路 + Rhie-Chow 一致动量插值 + 自适应压力欠松弛，
> 连续性 + 动量残差 + 压力步长三判据）、瞬态 PISO、ASCII `.vtu` 输出，
> k-ε / k-ω SST 两方程 RANS 与 Spalding 壁面函数（含平衡壁面单元），以及
> VOF 首增量（相分数守恒输运 + 界面法向/曲率，速度给定的 `chflow-vof`）、
> VOF 的连续表面力（CSF 与 PISO 动量耦合，静止液滴算例验 Laplace 压差）与
> 两相耦合（解出通量输运 α，椭圆液滴振荡验密度匹配的 2D Lamb 周期）；
> VOF 的几何界面重构、Euler-Euler 双流体、张量场、
> MPI 域分解与 Eigen/PETSc 后端仍待实现。

#### 5. 绑定层 `src/bindings`
基于 pybind11：
- 将 1d / 3d / pvt 内核对象导出为 Python Class；
- 支持在 Python 脚本中构建管网、修改PVT、提交求解、读取场数据；
- 面向科研人员做批量DOE、参数反演、自动后处理。

#### 6. CLI 工具层 `cli/`（Rust）
作为自动化入口，适合 HPC、Slurm 集群、CI 回归测试：
```bash
# 示例
chuanflow run steady case.yaml
chuanflow run transient case.yaml
chuanflow run pvt case.dict
chuanflow mesh convert gmsh.msh chuan.cmesh
chuanflow post extract pressure along-pipe
```

#### 7. GUI 层 `gui/`（Tauri + egui）
桌面交互前端，定位是**建模与监控器，而不是求解器**：
- 一维管网可视化拖拽编辑器；
- 瞬态曲线实时绘图（压力、持液率、温度）；
- 算例树管理、日志流查看；
- 一键唤起 ParaView 浏览三维结果；
- 通过子进程调用 `chuanflow-cli`，和求解内核进程隔离。

#### 8. Python SDK & Examples `python/`
封装高层API，降低用户上手门槛：
```python
import chuanflow

# 快速示例示意
net = chuanflow.Network()
net.add_well("well-1", rate=100)
net.add_pipe(from_node="well-1", to_node="sep", length=2000)
solver = chuanflow.SteadySolver(net)
result = solver.run()
print(result.pressure_profile)
```

## 🛠️ 技术栈清单
- **C++20**：1d / 3d / pvt / common 内核，CMake 构建
- **Rust**：CLI、Tauri GUI
- **Python 3.10+**：上层脚本、后处理
- **pybind11**：C++ ↔ Python 绑定
- **Eigen**：稠密线性代数、小规模求解
- **PETSc**：大规模稀疏线性求解、MPI并行
- **Kokkos**（可选，后续）：异构CPU/GPU加速
- **VTK + HDF5**：网格与场数据IO
- **spdlog / fmt / nlohmann/json**：基础工具库
- **egui + Tauri**：跨平台桌面GUI

## 🚀 开发路线（Milestone）
- M1：一维稳态管网求解器 + Python绑定 + CLI，对标基础Pipesim算例，完成基准校验 ✅
- M2：一维瞬态多相求解，实现段塞基础模型、水合物计算，对标OLGA简单工况 🚧
  （瞬态水力 + 段塞持液率输运 + 水合物等温动力学 + 3 个解析回归已落地，
  真两相动量交换与热力耦合待做）
- M3：三维FVM基础求解器，单相对流扩散、k-ε，网格读写、MPI并行 🚧
  （多面体网格 + fvc/fvm 离散 + 对流扩散 + 不可压 SIMPLE（全压力回路 +
  Rhie-Chow 一致动量插值，Poiseuille 算例 1411 → 47 外迭代）+ 瞬态 PISO +
  blockMesh/Gmsh msh 读入 + vtu IO + k-ε / k-ω SST 湍流模型层与 Spalding 壁面
  函数（湍流渠道摩阻对 Pope 律回归，SIMPLE 与 PISO 两种耦合共用同一套断言）+
  VOF 首增量（相分数守恒输运的 M-矩阵性质、界面法向与曲率，平移与刚体旋转
  两个解析回归）已落地，MPI 并行与 Eigen/PETSc 后端待做）
- 统一物性层 `src/pvt`（上面两个里程碑共用的底座，见分层说明第 2 层）🚧
  （Peng-Robinson 状态方程 + 多组分等温 TP 闪蒸 + 混合物泡点 + PVT 双线性表 +
  `chflow-pvt` 应用 + Python 绑定 + CLI 已落地，EOS 内部不变量与解析回归全绿；
  水合物相平衡/生成条件、黏度与界面张力、非等温闪蒸、以及把密度实时供给
  1D/3D 求解器待做）
- M4：三维多相流模型，打通 1D ↔ 3D 耦合边界接口 🚧
  （VOF 连续表面力 CSF 与 PISO 动量耦合（静止液滴 Laplace 压差算例）与两相
  耦合——解出通量输运 α、椭圆液滴振荡验密度匹配 Lamb 周期——已落地；
  几何界面重构、Euler-Euler 双流体在此）
- M5：Tauri GUI 管网编辑器、实时曲线面板
- M6：Kokkos GPU 加速、更多工业模型扩充

已落地部分的实现细节与概念映射见 [docs/architecture.md](docs/architecture.md)。

## 📌 相关项目说明 & 免责
- ChuanFlow **借鉴 OpenFOAM-dev 的有限体积框架设计思路**，属于学习参考架构思想；**不包含 OpenFOAM 源码**，不衍生自 OpenFOAM；
- OpenFOAM® 是 OpenFOAM Foundation 的商标，本项目和 OpenFOAM Foundation、ESI Group 无关联；
- 项目目标是服务国内油气、能源行业多相流仿真需求，开源共建，欢迎高校、工业用户提交算例、模型与PR。

## 🤝 参与贡献
欢迎 Issue、Discussion、PR。
开发规范、代码风格、测试准则放在 `docs/CONTRIBUTING.md`。

## 📄 License
Apache-2.0，详见 LICENSE 文件。
