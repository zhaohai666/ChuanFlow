// Kokkos CUDA compile/execution gate for the cf::backend::parallelFor
// contract (README M6 / D13 remainder).
//
// The five cf3d kernels wired through parallelFor (§5.126) are host lambdas
// over host memory; under Kokkos_ENABLE_CUDA their __host__-only bodies are
// rejected by the __host__ __device__ KOKKOS_LAMBDA wrapper. Full cf3d device
// offload is future work. This probe instead exercises the abstraction
// itself the way a GPU kernel will use it: a KOKKOS_FUNCTION functor over a
// raw pointer, driven through parallelFor, and cross-checked index-for-index
// against serialFor.
//
// Compile gate: builds clean with a Kokkos CUDA installation using nvcc even
// where no GPU exists (GitHub Actions ubuntu-latest). Execution: on a real
// device (Kaggle/Colab, see scripts/gpu/) the CUDA-written array must match
// the serial reference bitwise, mirroring verify.backendCpu's criterion.
//
// usage: cf-kokkos-probe [n]   (default n = 1<<20)

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "cf/parallel/Backend.h"

#if defined(CF_WITH_KOKKOS)
#include <Kokkos_Core.hpp>
#else
#define KOKKOS_FUNCTION
#endif

namespace
{

// Pure per-index map, device-callable, disjoint writes -- the same contract
// the domain-local kernels satisfy on the host backend.
//
// Both operations are exact in IEEE double for i < 2^23 (1.5*i needs one
// extra mantissa bit, +0.25 stays on the 2^-2 grid), so an FMA contraction
// on the device cannot diverge from the host serial reference bitwise.
struct ScaleKernel
{
    double* out;
    KOKKOS_FUNCTION void operator()(int i) const
    {
        out[i] = static_cast<double>(i) * 1.5 + 0.25;
    }
};

} // namespace

int main(int argc, char** argv)
{
#if defined(CF_WITH_KOKKOS)
    Kokkos::initialize(argc, argv);
    {
#endif
    const int n = (argc > 1) ? std::atoi(argv[1]) : (1 << 20);
    if (n <= 0 || n >= (1 << 23))
    {
        std::cerr << "usage: cf-kokkos-probe [n], 0 < n < 2^23\n";
        return 2;
    }

#if defined(CF_WITH_KOKKOS) && defined(KOKKOS_ENABLE_CUDA)
    Kokkos::View<double*> dev("probe", n);
    cf::backend::parallelFor(n, ScaleKernel{dev.data()});
    auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), dev);
    std::vector<double> par(host.data(), host.data() + static_cast<std::size_t>(n));
#else
    std::vector<double> par(static_cast<std::size_t>(n), 0.0);
    cf::backend::parallelFor(n, ScaleKernel{par.data()});
#endif

    std::vector<double> ser(static_cast<std::size_t>(n), 0.0);
    cf::backend::serialFor(n, ScaleKernel{ser.data()});

    std::size_t bad = 0;
    std::size_t firstBad = 0;
    for (int i = 0; i < n; ++i)
    {
        const std::size_t u = static_cast<std::size_t>(i);
        if (std::memcmp(&par[u], &ser[u], sizeof(double)) != 0)
        {
            if (bad == 0) firstBad = u;
            ++bad;
        }
    }

    std::cout << "probe backend: " << cf::backend::name() << ", n = " << n << "\n";
#if defined(CF_WITH_KOKKOS)
    std::cout << "probe execution space: "
              << Kokkos::DefaultExecutionSpace::name() << "\n";
#endif
    if (bad != 0)
    {
        std::cout << "PROBE FAIL: " << bad << " bitwise mismatches, first at i = "
                  << firstBad << " (" << par[firstBad] << " vs " << ser[firstBad] << ")\n";
        return 1;
    }
    std::cout << "PROBE PASS: parallelFor output == serialFor reference bitwise\n";
#if defined(CF_WITH_KOKKOS)
    }
    Kokkos::finalize();
#endif
    return 0;
}
