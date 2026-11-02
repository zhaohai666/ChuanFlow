#pragma once

#include <cstdint>
#include <cstddef>

namespace cf
{

using Label = std::int64_t;
using label = std::int64_t;

/// Physical quantity type: every scalar in the kernels is double.
using scalar = double;

constexpr label labelMax = 9223372036854775807LL;

} // namespace cf
