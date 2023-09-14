#ifndef _FLAMES_TYPE_HPP_
#define _FLAMES_TYPE_HPP_

#include <ap_fixed.h>
#include <type_traits>

namespace flames {

template <int I, int D>
using FxP_s = ap_fixed<I + D + 1, I>;

template <int I, int D>
using FxP = FxP_s<I, D>;

template <int I, int D>
using FxP_u = ap_ufixed<I + D, I>;

} // namespace flames

#endif
