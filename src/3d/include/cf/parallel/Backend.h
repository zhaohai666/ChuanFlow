#pragma once

#include <cstddef>
#include <type_traits>

#if defined(CF_WITH_KOKKOS)
#include <Kokkos_Core.hpp>
#endif

namespace cf::backend
{

/// Compile-time selected parallel backend. The canonical CPU serial member is
/// the reference every other member must agree with, index by index.
enum class Kind
{
    CpuSerial,
    OpenMP,
    Kokkos,
};

/// Human-readable name of the backend compiled into this translation unit.
/// Defined in src/3d/src/parallel/backend.cpp (kept out of the header so the
/// build configuration has exactly one owner).
const char* name();
Kind kind();

/// Invoke f(i) exactly once for every i in [0, n), in ascending index order
/// on the serial backend. The loop is intentionally index-based (not a pair
/// of begin/end iterators) so the same call site maps onto a GPU flat index.
///
/// Contract for every backend:
///   * one call per index, no index skipped or repeated;
///   * f must be safe to call from multiple threads / device lanes, i.e. it
///     may only write to disjoint locations (the domain-local kernels in
///     domainOperators.cpp satisfy this: owner slots are unique per domain).
///
/// The default (no CF_WITH_* macro) is a plain serial for loop, which gives
/// byte-for-byte identical floating-point results for reduction-style callers
/// that accumulate afterwards in serial. CF_WITH_OPENMP keeps f(i) per-index
/// independent so that downstream serial reduction stays bitwise identical.
template<class Idx, class F>
void parallelFor(Idx n, F&& f)
{
    static_assert(std::is_integral_v<Idx>, "cf::backend::parallelFor needs an integral index");
    if (n <= Idx{0}) return;

#if defined(CF_WITH_KOKKOS)
    // Kokkos SDK is required: configure with
    //   cmake -DCFLOW_WITH_KOKKOS=ON -DKokkos_DIR=<kokkos-install>/lib/cmake/Kokkos
    // and call Kokkos::initialize(argc, argv) / Kokkos::finalize() in main().
    // The callable must be device-callable (KOKKOS_FUNCTION) and capture by
    // value; host-only lambdas will not compile here. This branch is not
    // built in this environment (no SDK, no CUDA) but is the full interface.
    Kokkos::parallel_for("cf::backend::parallelFor", static_cast<std::size_t>(n),
                         KOKKOS_LAMBDA(const std::size_t i)
                         {
                             f(static_cast<Idx>(i));
                         });
    Kokkos::fence();
#elif defined(CF_WITH_OPENMP)
    // -fopenmp is supplied by find_package(OpenMP) in the top-level CMake.
    // schedule(static) keeps each index's work together; the *body* f(i) is
    // what is parallelised, and any reduction is deliberately left to the
    // caller's serial pass so results stay bitwise equal to CpuSerial.
#pragma omp parallel for schedule(static)
    for (Idx i = 0; i < n; ++i) f(i);
#else
    for (Idx i = 0; i < n; ++i) f(i);
#endif
}

/// Serial reference: always the same ascending loop, independent of build
/// flags. Used by the verification app to assert the active backend matches
/// the reference index-for-index.
template<class Idx, class F>
void serialFor(Idx n, F&& f)
{
    static_assert(std::is_integral_v<Idx>, "cf::backend::serialFor needs an integral index");
    for (Idx i = 0; i < n; ++i) f(i);
}

} // namespace cf::backend
