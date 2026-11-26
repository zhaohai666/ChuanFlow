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

} // namespace cf::backend
