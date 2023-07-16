/**
 * @file tensor.hpp
 * @author Wuqiong Zhao (me@wqzhao.org), et al.
 * @brief Tensor (3D Array) for FLAMES
 * @version 0.1.0
 * @date 2023-07-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _FLAMES_TENSOR_HPP_
#define _FLAMES_TENSOR_HPP_

#ifndef _FLAMES_CORE_HPP_
#    include "core.hpp"
#endif

#ifndef FLAMES_TENSOR_PARTITION_COMPLETE
#    ifdef FLAMES_MAT_PARTITION_COMPLETE
#        define FLAMES_TENSOR_PARTITION_COMPLETE
#    endif
#endif

namespace flames {
template <typename T, size_t n_rows, size_t n_cols, size_t n_slices, MatType type>
class Tensor {
  public:
    using element_type = T;
    using value_type   = T;
    using View         = MatView<T, n_rows, n_cols, type>;

    Tensor() {
#ifdef FLAMES_TENSOR_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
    }

    inline static constexpr size_t matSize() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    inline static constexpr size_t size() noexcept { return n_slices * size(); }

    inline View slice(size_t index) const {
        assert(index < n_slices && "Index should be within in range for MatView::slice(index).");
        return const_cast<T*>(_data + index * matSize());
    }

    inline View slice(size_t index) {
        assert(index < n_slices && "Index should be within in range for MatView::slice(index).");
        return const_cast<T*>(_data + index * matSize());
    }

    inline View operator[](size_t index) const { return slice(index); }

    inline View operator[](size_t index) { return slice(index); }

  private:
    T _data[type == MatType::NORMAL     ? n_slices * n_rows * n_cols
            : type == MatType::DIAGONAL ? n_slices * n_rows
            : type == MatType::SCALAR   ? n_slices * 1
            : type == MatType::SUPPER   ? n_slices * (n_rows - 1) * n_rows / 2
            : type == MatType::SLOWER   ? n_slices * (n_rows - 1) * n_rows / 2
            : type == MatType::ASYM     ? n_slices * (n_rows - 1) * n_rows / 2
                                        : n_slices * (1 + n_rows) * n_rows / 2];
};

} // namespace flames

#endif
