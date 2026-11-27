# Cloud GPU execution gate for the Kokkos backend (D13 remainder / F14)

本机无 NVIDIA/ROCm（唯一 GPU 是 RX 590 Polaris：无 CUDA、Windows ROCm 不支持该代、无 WSL、
Kokkos 4.x 已删 OpenCL 后端），所以 CUDA 执行验证在免费云 GPU 上跑。CI 侧的
`kokkos-cuda-compile` job 只编译不执行（ubuntu runner 没有 GPU 设备）；本目录的脚本负责实跑。

## 用法（Kaggle，推荐：T4/P100 免费额度大、开机即用）

1. 新建 Notebook，Settings → Accelerator 选 **GPU T4 x2**（或 P100），确认 internet 开启。
2. 第一个 cell 粘贴运行：

```bash
!git clone --depth 1 https://github.com/zhaohai666/ChuanFlow.git /kaggle/cf
!bash /kaggle/cf/scripts/gpu/kokkos_gpu_verify.sh
```

3. 期望输出（Colab 同理，去掉行首 `!`）：

```
GPU: Tesla T4 (sm_75)
...
probe backend: kokkos
probe execution space: cuda
PROBE PASS: parallelFor output == serialFor reference bitwise
GPU execution gate passed: ...
```

把完整 cell 输出回贴，即可按 §5.239 判据回填台账（F14 GPU 🟡 → ✅）。

## 判据

`cf::backend::parallelFor`（Kokkos CUDA 分支）写 `Kokkos::View`，`deep_copy` 回 host 后与
`serialFor` 参照 `memcmp` **逐位**对照（kernel 为 `i*1.5+0.25`，`n < 2^23` 内 IEEE 精确，
FMA 收缩不可能引入分歧）。这与 `verify.backendCpu` 的 active==serial 逐位口径一致。

## 失败排查

- `cmake >= 3.24 not found`：脚本会自动 `pip install cmake ninja`；若 pip 源被墙可先
  `pip config` 或改用 Modal/RunPod。
- `no suitable GPUs found`（编译能过、运行报此）：加速器没选 GPU 或驱动未就绪。
- `nvcc` 不在 `/usr/local/cuda`：设 `CUDA_HOME=<path>` 再跑。
- 编译过但 `PROBE FAIL`：这是真结果，原样回贴，不要重跑覆盖日志。
