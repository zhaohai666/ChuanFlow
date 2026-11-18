#include "cf/parallel/Backend.h"

namespace cf::backend
{

const char* name()
{
#if defined(CF_WITH_KOKKOS)
    return "kokkos";
#elif defined(CF_WITH_OPENMP)
    return "openmp";
#else
    return "cpu-serial";
#endif
}

Kind kind()
{
#if defined(CF_WITH_KOKKOS)
    return Kind::Kokkos;
#elif defined(CF_WITH_OPENMP)
    return Kind::OpenMP;
#else
    return Kind::CpuSerial;
#endif
}

} // namespace cf::backend
