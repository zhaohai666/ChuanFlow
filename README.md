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
- **3D FVM 三维CFD求解器**：基于有限体积法，支持单相/多相流、湍流、换热（不可压能量方程：对流扩散 + 壁面热流 + Boussinesq 浮力 + 湍流释热 ε/cp）；网格体系、场与离散框架参考 OpenFOAM-dev 思想，自研实现
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
> 两相耦合（解出通量输运 α，椭圆液滴振荡验密度匹配的 2D Lamb 周期）与
> 几何界面重构 PLIC（逐单元平面切出精确截体积，α 带 60/60、平面读回 α 到 2e-15）与
> 几何输运（α 由平面带走：Zalesak 转盘整圈守恒到 9.4e-15、α 精确在 [0,1]、
> 锋面转一圈仍是 76 格而不是迎风摊出的 421 格）与邻域最小二乘平面法向 `fit`
> （法向对解析径向 1.16° 对 mixed 的 8.34°、随加密收敛到 0.30°，真平面场方向精确到
> 1e-6°；同一转盘带 65 格、形状误差 3.33% 对 3.85%）与截面几何曲率（带内截面点二次曲面，
> κ 带均值 −5.1 对 −1/R = −5 而散度模板停在 −3.5，平面锋至舍入级；`kappaPlic` 可喂给 CSF）；
> Euler-Euler 双流体、张量场、
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
  （黑油 PVT 相关式 + Beggs-Brill 两相管流已接入稳态管网：井筒/集输两相压力分布、
  流型与持液率，两相算例与干气对照回归全绿；**IPR–VLP 节点分析**已落地——
  linear/Vogel/Jones 闭式流入曲线 × Beggs-Brill 井筒 VLP，工作点取 f(q)=
  p_vlp(q)−p_ipr(q) 的最大根，活井/干油孪生/死井三算例回归；**设备单元**
  已落地——阀 kq² 损失系数 + 离心泵扬程曲线接入管 law，流量反解对流量
  括根（泵扬程定点迭代在扬程下降快于系统上升时不收缩），阀/闭环泵/
  泵+阀串联网三算例对闭式二次根回归；**表驱动单相物性**已落地——
  fluid{pvt{pressure;rho;mu}} 分段线性 ρ(p)/μ(p) 接入稳态管网，24 段局部
  能量平衡 march（质量流量守恒、气体沿程膨胀）+ 单调律 Illinois 割线，
  平表严格退化常物性、等温气柱静水闭式 exp(−gH/c) 回归）
- M2：一维瞬态多相求解，实现段塞基础模型、水合物计算，对标OLGA简单工况 🚧
  （集总瞬态水力 + 段塞持液率输运 + 水合物等温动力学，以及**分布式两流体
  内核**——逐管段 α 输运的 Godunov 漂移流模型（Wallis/Zuber 闭合）+ 混合
  动量，激波位置/稳态持液率/重力分异/层流启动 5 个解析回归——已落地；
  **井筒热力耦合**——两流体核的分布式混合温度方程（相加权 ρc 容量 +
  迎风焓输运 + 壁面/地层换热，恒壁温或地温线模型），弛豫闭式/换热器
  递推/关井井温对地温线 3 个回归——已落地；
  **动态滑移闭合**——逐段输运滑移态 s = v_G − v_L 的双相动量降阶模型
  （K_s 滑移惯性 + 界面阻力 σ = C·s|s| 半隐式，终端滑移 s* = √(drive/CD)
  精确不动点，混合动量带滞后滑移惯性项；重力分异/耦合稳态持液 2 个
  回归，另修正混合动量柱惯性系数为 ρ̄L）——已落地；
  气相 EOS 压缩性、接入管网待做）
- M3：三维FVM基础求解器，单相对流扩散、k-ε，网格读写、MPI并行 🚧
  （多面体网格 + fvc/fvm 离散 + 对流扩散 + 不可压 SIMPLE（全压力回路 +
  Rhie-Chow 一致动量插值，Poiseuille 算例 1411 → 47 外迭代）+ 瞬态 PISO +
  blockMesh/Gmsh msh 读入 + vtu IO + k-ε / k-ω SST 湍流模型层与 Spalding 壁面
  函数（湍流渠道摩阻对 Pope 律回归，SIMPLE 与 PISO 两种耦合共用同一套断言）+
  **能量方程**（对流扩散耦合连续性通量、`heatFlux` 壁面边界与离散焓平衡恒等式、
  Boussinesq 浮力；三道无流动传导 + 余弦模态衰减 + 闭盒静止 + 双壁热流渠道
  六个回归，渠道算例 `Nu_Dh` 8.2505 对 140/17）+
  **热壁面函数**（两支 Jayatilleke 温度律 `α_face = u_tau·y/T+` 接入 fixedValue
  壁面，副层极限退回分子导热；受热湍流渠道回归：报告 T+ 对 Python 重写律全同、
  场对重述离散解偏差 4e-5 K、壁函数通量为分子首格通量的 3.08 倍，双模型共用断言）+
  **湍流释热**（ε 热源接入能量方程：k-ε 直用输运变量、SST 按 β* k ω 换算；
  渠道离散焓平衡 SIMPLE 8.7e-16 / PISO 2.9e-9 W、source = ∫ε(vtu) 到 1e-8 相对、
  壁单元 ε = u_τ³/(κ y₁) = 475.23，source/pump 落在泵功可见性带 [0.5, 2]——
  壁函数壁单元非能量守恒，实测 0.76~1.18）+
  VOF 首增量（相分数守恒输运的 M-矩阵性质、界面法向与曲率，平移与刚体旋转
  两个解析回归）+ **`p_rgh` 表述**（均匀体力被 `p_rgh = p − force·x` 精确吸收，
  PISO/SIMPLE 共用 `pRgh` 开关；闭盒重力从 plain-p 的蠕升 uMax 1.83 m/s
  （持续增长）到精确静止 uMaxAll = 0，重建真压程 9.56475 = 9.81·0.975）
  已落地，热壁函数、MPI 并行与 Eigen/PETSc 后端待做）
- 统一物性层 `src/pvt`（上面两个里程碑共用的底座，见分层说明第 2 层）🚧
  （Peng-Robinson 状态方程 + 多组分等温 TP 闪蒸 + 混合物泡点 + PVT 双线性表 +
  黑油物性相关式（Standing 泡点、Vazquez-Beggs Rs/Bo/co、Beggs-Robinson 与
  Lee 黏度、DAK Z）+ `chflow-pvt` 应用 + Python 绑定 + CLI 已落地，
  EOS 内部不变量与解析回归全绿；
  水合物相平衡/生成条件、EOS 黏度与界面张力、非等温闪蒸、以及把密度实时供给
  1D/3D 求解器待做）
- M4：三维多相流模型，打通 1D ↔ 3D 耦合边界接口 🚧
  （VOF 连续表面力 CSF 与 PISO 动量耦合（静止液滴 Laplace 压差算例）与两相
  耦合——解出通量输运 α、椭圆液滴振荡验密度匹配 Lamb 周期——已落地；
  **几何界面重构 PLIC**（α 带逐单元平面、牛顿解偏移保证截体积、切面面积与质心，
  旋转圆盘算例 `vofPlic` 断言平面自洽与相体积守恒）与 **几何输运**（α 由平面自己带走：
  逐面扫过棱柱求交得通量、两侧共享同一数因此构造性守恒、出/入流双限迭代到不动点保 α ∈ [0,1]；
  Zalesak 转盘 `vofZalesak` 转完整一圈相体积守到 9.4e-15、锋面仍是 76 格而非迎风摊出的 421 格、
  形状误差 3.85% 圆盘面积 = 边界位移 1.92% R）已落地；
  **邻域最小二乘平面法向 `fit`**（把邻域各格 α 各自切准的平面间的倾斜当方向，阻尼 Gauss–Newton
  四轮重切：40/80/160 加密下均值 1.16°→0.488°→0.304° 而单元内取法停在 7.6° 不动、真平面场
  方向精确到 1e-6°；同一转盘 `vofFit` 带 65 格、Σ截面 0.1282 对周长 0.1257、形状误差 3.33%）
  与 **κ 由截面几何直接给出**（带内两环角邻域的截面点加权二次曲面，κ = div n̂，欠定如实报 0：
  圆盘带均值 −5.02…−5.12 对 −1/R = −5 而散度模板停在 −3.25…−3.90、平面锋 κ 至舍入级 2.77e-14；
  `surfaceTension{ kappaPlic true; }` 接入 CSF——锐利场上 κ 对到 0.7%，胶式压力对 delta 截面力
  的跃变位置由 `vofDropletPlic` 按实测断言）
  已落地；Euler-Euler 双流体在此）
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
