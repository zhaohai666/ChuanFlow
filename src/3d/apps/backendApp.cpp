// M6 GPU abstraction layer -- CPU-verifiable closed loop.
//
// Runs one domain-local, write-conflict-free kernel through cf::backend::parallelFor
// and compares its checksum, index by index and then as a serial reduction,
// against both the serial reference (cf::backend::serialFor) and the analytic
// closed form. With the default CPU build the two loops are the same machine
// code; with OpenMP the body is still per-index independent, so the trailing
// serial reduction is bitwise identical. The Kokkos branch is compiled and
// run-verified against a real Kokkos (Serial host execution space);
// CUDA/HIP GPU execution still requires GPU hardware absent on this host.
//
// usage: chflow-3d-backend <case.dict>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "cf/common/Dictionary.h"
#include "cf/common/Error.h"
#include "cf/common/Labels.h"
#include "cf/common/Log.h"
#include "cf/parallel/Backend.h"

#if defined(CF_WITH_KOKKOS)
#include <Kokkos_Core.hpp>
// Kokkos requires explicit runtime setup; the guard finalizes on every exit
// path so the checksum comparison below runs in a live default execution
// space (Serial on this host -- no CUDA/AMD GPU hardware is available).
namespace
{
struct KokkosRuntimeGuard
{
    ~KokkosRuntimeGuard() { Kokkos::finalize(); }
};
} // namespace
#endif


using namespace cf;

namespace
{
std::string num(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

struct Checksum
{
    double sum = 0.0;
    double max = 0.0;
};

Checksum reduce(const std::vector<double>& v)
{
    Checksum c;
    for (const double x : v)
    {
        c.sum += x;
        c.max = std::max(c.max, x);
    }
    return c;
}

double kernelValue(label i, const std::string& mode)
{
    if (mode == "square") return static_cast<double>(i) * static_cast<double>(i);
    // "divide": domainGradDivide style -- out[i] = sum[i] / volume[i].
    // sum[i] = i + 1, volume[i] = 2, so the result is a known closed form.
    return static_cast<double>(i + 1) / 2.0;
}

Checksum analytic(label n, const std::string& mode)
{
    const double nd = static_cast<double>(n);
    if (mode == "square")
    {
        const double m = nd - 1.0;
        return {m * nd * (2.0 * nd - 1.0) / 6.0, m * m};
    }
    return {nd * (nd + 1.0) / 4.0, nd / 2.0};
}
} // namespace

int main(int argc, char** argv)
{
#if defined(CF_WITH_KOKKOS)
    Kokkos::initialize(argc, argv);
    KokkosRuntimeGuard kokkosGuard;
#endif
    if (argc < 2)
    {
        std::cerr << "usage: chflow-3d-backend <case.dict>\n";
        return 2;
    }

    try
    {
        const Dictionary dict = Dictionary::parseFile(argv[1]);
        const label n = static_cast<label>(dict.getScalar("n"));
        const std::string mode = dict.getWordOr("mode", "square");

        if (!(n > 0)) cfFatal("backend: n must be positive (got " + std::to_string(n) + ")");
        if (mode != "square" && mode != "divide")
            cfFatal("backend: unknown mode '" + mode + "' (expected square|divide)");

        log::info("backend name ", backend::name(), " n ", n, " mode ", mode);

        std::vector<double> par(static_cast<std::size_t>(n), 0.0);
        std::vector<double> ser(static_cast<std::size_t>(n), 0.0);
        backend::parallelFor(n, [&](label i) { par[static_cast<std::size_t>(i)] = kernelValue(i, mode); });
        backend::serialFor(n, [&](label i) { ser[static_cast<std::size_t>(i)] = kernelValue(i, mode); });

        const Checksum cp = reduce(par);
        const Checksum cs = reduce(ser);
        const Checksum ca = analytic(n, mode);

        log::info("backend checksum active sum ", num(cp.sum), " max ", num(cp.max));
        log::info("backend checksum serial sum ", num(cs.sum), " max ", num(cs.max));
        log::info("backend checksum analytic sum ", num(ca.sum), " max ", num(ca.max));

        const bool bitwise = par == ser;
        const bool equalSum = (cp.sum == cs.sum) && (cs.sum == ca.sum);
        const bool equalMax = (cp.max == cs.max) && (cs.max == ca.max);
        log::info("backend elementwiseMatch ", bitwise ? 1 : 0,
                  " checksumMatch ", (equalSum && equalMax) ? 1 : 0);
        if (!bitwise || !equalSum || !equalMax)
            cfFatal("backend: parallel result differs from serial/analytic checksum");
        return 0;
    }
    catch (const std::exception& e)
    {
        log::error(e.what());
        return 1;
    }
}
