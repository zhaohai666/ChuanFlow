#!/usr/bin/env bash
# Real-GPU execution gate for cf::backend::parallelFor (D13 remainder / F14).
#
# Designed for a free cloud GPU notebook (Kaggle T4/P100 "GPU" accelerator or
# Colab T4, runtime with internet enabled; also works on any CUDA box).
# It builds Kokkos 4.4.01 with the CUDA execution space, builds the device
# contract probe from scripts/ci/, and runs it: the CUDA-written array must
# match the serial reference index-for-index, bitwise.
#
# usage (Kaggle notebook cell):
#   !git clone https://github.com/zhaohai666/ChuanFlow.git /kaggle/cf
#   !bash /kaggle/cf/scripts/gpu/kokkos_gpu_verify.sh
# on Colab put the same two lines in one cell (drop the leading '!').
set -euo pipefail

KOKKOS_VERSION=4.4.01
N=${1:-1048576}

if ! command -v nvidia-smi >/dev/null || ! nvidia-smi -L >/dev/null 2>&1; then
    echo "FAIL: no NVIDIA GPU visible (nvidia-smi missing or empty)"; exit 1
fi
GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
GPU_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d .)
echo "GPU: ${GPU_NAME} (sm_${GPU_ARCH})"

REPO=$(cd "$(dirname "$0")/../.." && pwd)
WORK=${WORK_DIR:-$PWD/kokkos-gpu-verify}
mkdir -p "$WORK"

CMAKE=$(command -v cmake || true)
if [ -z "$CMAKE" ] || [ "$("$CMAKE" --version | head -1 | awk '{print $3}' | cut -d. -f2)" -lt 24 ]; then
    echo "cmake >= 3.24 not found; installing via pip"
    python3 -m pip install -q --upgrade cmake ninja
    hash -r
fi

if [ ! -d "$WORK/kokkos" ]; then
    git clone --branch "$KOKKOS_VERSION" --depth 1 https://github.com/kokkos/kokkos.git "$WORK/kokkos"
fi

CUDA_HOME=${CUDA_HOME:-/usr/local/cuda}

cmake -S "$WORK/kokkos" -B "$WORK/kokkos-build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DKokkos_ENABLE_CUDA=ON \
    -DKokkos_ENABLE_SERIAL=ON \
    -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" \
    -DCMAKE_CUDA_ARCHITECTURES="$GPU_ARCH" \
    -DCMAKE_INSTALL_PREFIX="$WORK/kokkos-install"
cmake --build "$WORK/kokkos-build" -j 4
cmake --install "$WORK/kokkos-build"

cmake -S "$REPO/scripts/ci" -B "$WORK/probe-build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DPROBE_WITH_KOKKOS=ON \
    -DKokkos_DIR="$WORK/kokkos-install/lib/cmake/Kokkos" \
    -DCMAKE_CUDA_ARCHITECTURES="$GPU_ARCH"
cmake --build "$WORK/probe-build" -j 4

echo "--- running probe (n = $N) ---"
"$WORK/probe-build/cf-kokkos-probe" "$N"
echo "GPU execution gate passed: cf::backend::parallelFor ran on the device"
echo "and matched the serial reference bitwise."
