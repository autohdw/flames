/**
 * @file core.hpp
 * @author Wuqiong Zhao (me@wqzhao.org), et al.
 * @brief Core Utilities for FLAMES
 * @version 0.1.0
 * @date 2023-07-15
 *
 * @copyright Copyright (c) 2022-2023 Wuqiong Zhao
 *
 */

#ifndef _FLAMES_CORE_HPP_
#define _FLAMES_CORE_HPP_

#include <ap_fixed.h>
#include <ap_int.h>
#include <cassert>
#include <complex>
#include <cstddef>
#include <fstream>
#include <hls_vector.h>
#include <iostream>
#include <type_traits>
#include <vector>

#ifndef __VITIS_HLS__
#    error "FLAMES library can only be used for Vitis HLS."
#endif

/*
 * There are some unwanted warning messages, including:
 *   WARNING: [HLS 207-1462] template template parameter using 'typename' is a C++17 extension
 *   WARNING: [HLS 207-5292] unused parameter '...'
 * They are safe to ignore, so all warnings of the FLAMES library are suppressed.
 * You can restore them by defining `FLAMES_PRESERVE_WARNING`
 */
#ifndef FLAMES_KEEP_WARNING
#    ifdef __SYNTHESIS__
#        pragma GCC diagnostic push
#        pragma GCC diagnostic ignored "-Weverything"
#    endif
#endif

#ifndef PRAGMA_SUB
#    define PRAGMA_SUB(x) _Pragma(#    x)
#endif
#ifndef FLAMES_PRAGMA
// Alias for #pragma HLS with support for macro expansion.
#    define FLAMES_PRAGMA(x) PRAGMA_SUB(HLS x)
#endif

#ifndef FLAMES_MAT_PLUS_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_PLUS_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_PLUS_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_MINUS_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_MINUS_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_MINUS_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_SET_VALUE_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_SET_VALUE_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_SET_VALUE_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_POWER_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_POWER_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_POWER_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_COPY_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_COPY_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_COPY_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR
#    ifdef FLAMES_MAT_TIMES_UNROLL_FACTOR
#        define FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR FLAMES_MAT_TIMES_UNROLL_FACTOR
#    else
#        ifdef FLAMES_UNROLL_FACTOR
#            define FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#        else
#            define FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR 32
#        endif
#    endif
#endif
#ifndef FLAMES_MAT_TIMES_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_TIMES_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_TIMES_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_EMUL_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_EMUL_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_EMUL_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_ABS_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_ABS_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_ABS_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_INV_UNROLL_FACTOR
#    ifdef FLAMES_UNROLL_FACTOR
#        define FLAMES_MAT_INV_UNROLL_FACTOR FLAMES_UNROLL_FACTOR
#    else
#        define FLAMES_MAT_INV_UNROLL_FACTOR 32
#    endif
#endif
#ifndef FLAMES_MAT_PARTITION_COMPLETE
#    ifndef FLAMES_MAT_PARTITION_FACTOR
#        define FLAMES_MAT_PARTITION_FACTOR 8
#    endif
#endif
#ifndef FLAMES_SORT_PARTITION_COMPLETE
#    ifndef FLAMES_SORT_PARTITION_FACTOR
#        define FLAMES_SORT_PARTITION_FACTOR 32
#    endif
#endif

#if defined __SYNTHESIS__ && defined FLAMES_PRINT_PER_MAT_COPY
#    undef FLAMES_PRINT_PER_MAT_COPY
#endif

#ifdef INLINE
#    define DEFINED_INLINE
#    undef INLINE
#endif

/**
 * @brief Namespace for the FLAMES library.
 *
 * @details When you include 'flames.hpp', this namespace is automatically used.
 */
namespace flames {

/**
 * @brief Matrix type for storage.
 *
 */
enum MatType {
    NORMAL,   /**< Normal matrix */
    DIAGONAL, /**< Diagonal matrix */
    SCALAR,   /**< Scalar matrix */
    UPPER,    /**< Upper triangular matrix */
    LOWER,    /**< Lower triangular matrix */
    SUPPER,   /**< Strict Upper triangular matrix */
    SLOWER,   /**< Strict Lower triangular matrix */
    SYM,      /**< Symmetrical matrix */
    ASYM      /**< Antisymmetrical matrix */
};

/// Normal matrix as a class type.
using MATTYPE_NORMAL = std::integral_constant<int, MatType::NORMAL>;
/// Diagonal matrix as a class type.
using MATTYPE_DIAGONAL = std::integral_constant<int, MatType::DIAGONAL>;
/// Scalar matrix matrix as a class type.
using MATTYPE_SCALAR = std::integral_constant<int, MatType::SCALAR>;
/// Upper triangular matrix as a class type.
using MATTYPE_UPPER = std::integral_constant<int, MatType::UPPER>;
/// Lower triangular matrix as a class type.
using MATTYPE_LOWER = std::integral_constant<int, MatType::LOWER>;
/// Strict upper triangular matrix as a class type.
using MATTYPE_SUPPER = std::integral_constant<int, MatType::SUPPER>;
/// Strict lower triangular matrix as a class type.
using MATTYPE_SLOWER = std::integral_constant<int, MatType::SLOWER>;
/// Symmetrical matrix as a class type.
using MATTYPE_SYM = std::integral_constant<int, MatType::SYM>;
/// Antisymmetrical matrix as a class type.
using MATTYPE_ASYM = std::integral_constant<int, MatType::ASYM>;

/**
 * @brief Get the MatType value from the class form.
 *
 * @tparam T The MatType class form.
 * @return (constexpr MatType) The MatType value (as an enum item).
 */
template <typename T>
inline constexpr MatType matType() {
    return MatType(T::value);
}

/**
 * @brief Get the class form of MatType.
 *
 * @tparam type The MatType enum item value.
 */
template <int type>
using MType = std::integral_constant<int, type>;

/**
 * @brief Summation type of two matrices.
 *
 * @param type1 The MatType of the left matrix.
 * @param type2 The MatType of the right matrix.
 * @return (constexpr MatType) The summation matrix type.
 */
inline constexpr MatType sumType(MatType type1, MatType type2) noexcept {
    if (type1 == type2) return type1;
    else if ((type1 == MatType::DIAGONAL && type2 == MatType::SCALAR) ||
             (type1 == MatType::SCALAR && type2 == MatType::DIAGONAL))
        return MatType::DIAGONAL;
    else if ((type1 == MatType::DIAGONAL && type2 == MatType::UPPER) ||
             (type1 == MatType::UPPER && type2 == MatType::DIAGONAL))
        return MatType::UPPER;
    else if ((type1 == MatType::DIAGONAL && type2 == MatType::LOWER) ||
             (type1 == MatType::LOWER && type2 == MatType::DIAGONAL))
        return MatType::LOWER;
    else if ((type1 == MatType::DIAGONAL && type2 == MatType::SUPPER) ||
             (type1 == MatType::SUPPER && type2 == MatType::DIAGONAL))
        return MatType::UPPER;
    else if ((type1 == MatType::DIAGONAL && type2 == MatType::SLOWER) ||
             (type1 == MatType::SLOWER && type2 == MatType::DIAGONAL))
        return MatType::LOWER;
    else if ((type1 == MatType::SCALAR && type2 == MatType::SUPPER) ||
             (type1 == MatType::SUPPER && type2 == MatType::SCALAR))
        return MatType::UPPER;
    else if ((type1 == MatType::SCALAR && type2 == MatType::SLOWER) ||
             (type1 == MatType::SLOWER && type2 == MatType::SCALAR))
        return MatType::LOWER;
    else if ((type1 == MatType::SUPPER && type2 == MatType::UPPER) ||
             (type1 == MatType::UPPER && type2 == MatType::SUPPER))
        return MatType::UPPER;
    else if ((type1 == MatType::SLOWER && type2 == MatType::LOWER) ||
             (type1 == MatType::LOWER && type2 == MatType::SLOWER))
        return MatType::LOWER;
    else if ((type1 == MatType::DIAGONAL && type2 == MatType::SYM) ||
             (type1 == MatType::SYM && type2 == MatType::DIAGONAL))
        return MatType::SYM;
    else if ((type1 == MatType::SCALAR && type2 == MatType::SYM) || (type1 == MatType::SYM && type2 == MatType::SCALAR))
        return MatType::SYM;
    else return MatType::NORMAL;
}

/**
 * @brief Multiplication type of two matrices.
 *
 * @param type1 The MatType of the left matrix.
 * @param type2 The MatType of the right matrix.
 * @param n_rows The number of rows of the left matrix.
 * @param comm The number of columns of the left matrix and the number of rows of the right matrix.
 * @param n_cols The number of columns of the right matrix.
 * @return (constexpr MatType) The multiplication matrix type.
 */
inline constexpr MatType mulType(MatType type1, MatType type2, size_t n_rows, size_t comm, size_t n_cols) noexcept {
    if (n_rows == comm && comm == n_cols) {
        if (type1 == type2) {
            if (type1 != SYM && type1 != ASYM) return type1;
        }
        if (type1 == MatType::SCALAR) return type2;
        else if (type2 == MatType::SCALAR) return type1;
        else if (type1 == MatType::DIAGONAL && (type2 == MatType::SUPPER || type2 == MatType::UPPER ||
                                                type2 == MatType::SLOWER || type2 == MatType::LOWER))
            return type2;
        else if (type2 == MatType::DIAGONAL && (type1 == MatType::SUPPER || type1 == MatType::UPPER ||
                                                type1 == MatType::SLOWER || type1 == MatType::LOWER))
            return type1;
        else if (type1 == MatType::SUPPER && type2 == MatType::UPPER) return type1;
        else if (type1 == MatType::SLOWER && type2 == MatType::LOWER) return type1;
        else if (type1 == MatType::UPPER && type2 == MatType::SUPPER) return type2;
        else if (type1 == MatType::LOWER && type2 == MatType::SLOWER) return type2;
    }
    return MatType::NORMAL;
}

/**
 * @brief Transpose type of a matrix.
 *
 * @param type The MatType of the matrix.
 * @return (constexpr MatType) The transpose matrix type.
 */
inline constexpr MatType tType(MatType type) noexcept {
    if (type == MatType::SUPPER) return MatType::SLOWER;
    else if (type == MatType::SLOWER) return MatType::SUPPER;
    else if (type == MatType::UPPER) return MatType::LOWER;
    else if (type == MatType::LOWER) return MatType::UPPER;
    else return type;
}

/**
 * @brief Calculate the row index of a upper triangular matrix.
 *
 * @param index The data index.
 * @param N The matrix dimension.
 * @return (constexpr size_t) The row index.
 */
inline constexpr size_t upperRow(size_t index, size_t N) {
    size_t r = 0;
    while (index >= N - r) {
        index -= N - r;
        ++r;
    }
    return r;
}

/**
 * @brief Calculate the row index of a lower triangular matrix.
 *
 * @param index The data index.
 * @param N The matrix dimension.
 * @return (constexpr size_t) The row index.
 */
inline constexpr size_t lowerRow(size_t index, size_t N) {
    size_t r = 0;
    while (index >= r + 1) {
        index -= r + 1;
        ++r;
    }
    return r;
}

/**
 * @brief Calculate the row index of a strict upper triangular matrix.
 *
 * @param index The data index.
 * @param N The matrix dimension.
 * @return (constexpr size_t) The row index.
 */
inline constexpr size_t supperRow(size_t index, size_t N) {
    size_t r = 0;
    while (index >= N - 1 - r) {
        index -= N - 1 - r;
        ++r;
    }
    return r;
}

/**
 * @brief Calculate the row index of a strict lower triangular matrix.
 *
 * @param index The data index.
 * @param N The matrix dimension.
 * @return (constexpr size_t) The row index.
 */
inline constexpr size_t slowerRow(size_t index, size_t N) {
    size_t r = 0;
    while (index >= r) {
        index -= r;
        ++r;
    }
    return r;
}

/**
 * @brief Matrix.
 *
 * @details The matrix design is hardware optimized,
 *          with easy pragma configuration to be specified by the user.
 *          The use of this class follows the routine of C++,
 *          with some exceptions for the sake of better hardware performance.
 *          Please check out descriptions of class member functions for details.
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type matrix type.
 */
template <typename T, size_t n_rows, size_t n_cols, MatType type = MatType::NORMAL>
class Mat;

template <typename T, size_t n_rows, size_t n_cols, size_t n_slices, MatType type = MatType::NORMAL>
class Tensor;

/**
 * @brief Column vector.
 *
 * @details This is the type alias of Mat and the number of columns is set to 1.
 * @tparam T Element type.
 * @tparam N The vector dimension.
 */
template <typename T, size_t N>
using Vec = Mat<T, N, 1>;

/**
 * @brief Row vector.
 *
 * @details This is the type alias of Mat and the number of row is set to 1.
 * @tparam T Element type.
 * @tparam N The vector dimension.
 */
template <typename T, size_t N>
using RowVec = Mat<T, 1, N>;

/**
 * @brief Read only view version of a matrix.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type.
 */
template <typename T, size_t n_rows, size_t n_cols, MatType type = MatType::NORMAL>
class MatView;

/**
 * @brief Read only view version of the opposite of a matrix.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type.
 */
template <typename T, size_t n_rows, size_t n_cols, MatType type>
class MatViewOpp;

/**
 * @brief Read only view version of a transposed matrix.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type.
 */
template <typename T, size_t n_rows, size_t n_cols, MatType type>
class MatViewT;

/**
 * @brief Read only view version of a diagonal matrix.
 *
 * @tparam T Element type.
 * @tparam N Matrix dimension.
 * @tparam N_ Matrix dimension (unused, only to ensure it is a square matrix).
 * @tparam type Matrix type (surely DIAGONAL here).
 * @tparam type_parent Parent matrix (where it takes the diagonal) type.
 */
template <typename T, size_t N, size_t N_, MatType type, typename type_parent = MATTYPE_NORMAL>
class MatViewDiagMat;

/**
 * @brief Read only view version of a diagonal matrix as vector.
 *
 * @tparam T Element type.
 * @tparam N Matrix dimension.
 * @tparam N_ Matrix dimension (unused, only to ensure it is a square matrix).
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the diagonal) type.
 */
template <typename T, size_t N, size_t N_, MatType type, typename type_parent = MATTYPE_NORMAL>
class MatViewDiagVec;

/**
 * @brief Read only view version of a diagonal matrix as row vector.
 *
 * @tparam T Element type.
 * @tparam N Matrix dimension.
 * @tparam N_ Matrix dimension (unused, only to ensure it is a square matrix).
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the diagonal) type.
 */
template <typename T, size_t N, size_t N_, MatType type, typename type_parent = MATTYPE_NORMAL>
class MatViewDiagRowVec;

/**
 * @brief Read only view version of a off diagonal.
 *
 * @tparam T Element type.
 * @tparam N Matrix dimension.
 * @tparam N_ Matrix dimension (unused, only to ensure it is a square matrix).
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the off diagonal) type.
 */
template <typename T, size_t N, size_t N_, MatType type, typename type_parent = MATTYPE_NORMAL>
class MatViewOffDiag;

/**
 * @brief Read only view version of a certain column as colunm vector.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the colunm vector) type.
 */
template <size_t Index, typename T, size_t n_rows, size_t n_cols, MatType type, typename type_parent = MATTYPE_NORMAL>
class MatViewCol;

/**
 * @brief Read only view version of a certain row as row vector.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the row vector) type.
 */
template <size_t Index, typename T, size_t n_rows, size_t n_cols, MatType type, typename type_parent = MATTYPE_NORMAL>
class MatViewRow;

/**
 * @brief Read only view version of successive columns.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the colunm vectors) type.
 */
template <size_t first_col, size_t last_col, typename T, size_t n_rows, size_t n_cols, MatType type,
          typename type_parent = MATTYPE_NORMAL>
class MatViewCols;

/**
 * @brief Read only view version of successive rows.
 *
 * @tparam T Element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the row vectors) type.
 */
template <size_t first_row, size_t last_row, typename T, size_t n_rows, size_t n_cols, MatType type,
          typename type_parent = MATTYPE_NORMAL>
class MatViewRows;

/**
 * @brief Read only view version of a discrete columns.
 *
 * @tparam Cols The number of the discrete columns taken by the container.
 * @tparam T Element type.
 * @tparam M Container element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the colunm vectors) type.
 */
template <size_t Cols, typename T, typename M, size_t n_rows, size_t n_cols, MatType type,
          typename type_parent = MATTYPE_NORMAL>
class MatViewColsContainer;

/**
 * @brief Afterwards action with initialization.
 *
 * @note This should only be used for view classes.
 */
enum class InitAfterwards {
    NONE, /**< none */
    OPP,  /**< opposite */
    TR    /**< transpose */
};

/**
 * @brief Read only view version of a discrete rows.
 *
 * @tparam Rows The number of the discrete rows taken by the container.
 * @tparam T Element type.
 * @tparam M Container element type.
 * @tparam n_rows Number of rows.
 * @tparam n_cols Number of columns.
 * @tparam type Matrix type (surely NORMAL here).
 * @tparam type_parent Parent matrix (where it takes the row vectors) type.
 */
template <size_t Rows, typename T, typename M, size_t n_rows, size_t n_cols, MatType type,
          typename type_parent = MATTYPE_NORMAL>
class MatViewRowsContainer;

template <typename T, size_t n_rows, size_t n_cols, MatType type>
class Mat {
    friend class MatView<T, n_rows, n_cols, type>;
    friend class MatViewOpp<T, n_rows, n_cols, type>;
    friend class MatViewT<T, n_cols, n_rows, type>;
    friend class MatViewT<T, n_rows, n_cols, type>;
    template <typename View_T, size_t View_N, size_t View_N_, MatType View_type, typename type_parent>
    friend class MatViewDiagMat;
    template <typename View_T, size_t View_N, size_t View_N_, MatType View_type, typename type_parent>
    friend class MatViewDiagVec;
    template <typename View_T, size_t View_N, size_t View_N_, MatType View_type, typename type_parent>
    friend class MatViewDiagRowVec;
    template <typename View_T, size_t View_N, size_t View_N_, MatType View_type, typename type_parent>
    friend class MatViewOffDiag;
    template <size_t Index, typename View_T, size_t View_n_rows, size_t View_n_cols, MatType View_type,
              typename type_parent>
    friend class MatViewCol;
    template <size_t Index, typename View_T, size_t View_n_rows, size_t View_n_cols, MatType View_type,
              typename type_parent>
    friend class MatViewRow;
    template <size_t first_col, size_t last_col, typename View_T, size_t View_n_rows, size_t View_n_cols,
              MatType View_type, typename type_parent>
    friend class MatViewCols;
    template <size_t first_row, size_t last_row, typename View_T, size_t View_n_rows, size_t View_n_cols,
              MatType View_type, typename type_parent>
    friend class MatViewRows;
    template <size_t Cols, typename View_T, typename View_M, size_t View_n_rows, size_t View_n_cols, MatType View_type,
              typename type_parent>
    friend class MatViewColsContainer;
    template <size_t Rows, typename View_T, typename View_M, size_t View_n_rows, size_t View_n_cols, MatType View_type,
              typename type_parent>
    friend class MatViewRowsContainer;
    template <typename View_T_T, size_t T_n_rows, size_t T_n_cols, size_t T_n_slices, MatType T_type>
    friend class Tensor;

  public:
    using element_type = T;
    using value_type   = T;

    /**
     * @brief Construct a new Mat object.
     *
     * @details This is the default constructor.
     *          You can set the array partition using macro
     *          `FLAMES_MAT_PARTITION_COMPLETE` to set a complete array partition
     *          and `FLAMES_MAT_PARTITION_FACTOR` to set a block partition with the specific factor.
     * @note Data is stored as a row major sequence.
     */
    Mat() {
        static_assert(n_rows != 0, "'rows' should be no smaller than 1.");
        static_assert(n_cols != 0, "'n_cols' should be no smaller than 1.");
        static_assert(type == MatType::NORMAL || n_rows == n_cols, "Square matrix 'rows' should be equal to 'n_cols'.");
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
    }

    /**
     * @brief Construct a new Mat object with initial value.
     *
     * @details All values are set to the initial value you specify.
     *          You can configure macro `FLAMES_MAT_SET_VALUE_UNROLL_FACTOR`
     *          to set the initialization unroll factor.
     * @param val The initial value.
     */
    Mat(T val) {
        static_assert(n_rows != 0, "'rows' should be no smaller than 1.");
        static_assert(n_cols != 0, "'n_cols' should be no smaller than 1.");
        static_assert(type == MatType::NORMAL || n_rows == n_cols, "Square matrix 'rows' should be equal to 'n_cols'.");
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
        setValue(val);
    }

    /**
     * @brief Copy constructor from a Mat object.
     *
     * @details Macro `FLAMES_MAT_COPY_UNROLL_FACTOR` to configure the unrolling factor.
     *          You can set the array partition using macro
     *          `FLAMES_MAT_PARTITION_COMPLETE` to set a complete array partition
     *          and `FLAMES_MAT_PARTITION_FACTOR` to set a block partition with the specific factor.
     * @param mat The matrix to be copied.
     */
    Mat(const Mat& mat) {
    MAT_COPY:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = mat[i];
        }
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
#ifdef FLAMES_PRINT_PER_MAT_COPY
        std::cout << "Mat copy!" << std::endl;
#endif
    }

    template <typename T2, size_t _rows, size_t _cols, MatType _type,
              std::enable_if_t<!std::is_same<T, T2>::value && type == _type && n_rows == _rows && n_cols == _cols,
                               bool> = true>
    Mat(const Mat<T2, _rows, _cols, _type>& mat) {
    MAT_COPY:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = mat[i];
        }
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
#ifdef FLAMES_PRINT_PER_MAT_COPY
        std::cout << "Mat copy!" << std::endl;
#endif
    }

    template <typename T2, size_t _rows, size_t _cols, MatType _type,
              std::enable_if_t<type != _type && n_rows == _rows && n_cols == _cols, bool> = true>
    Mat(const Mat<T2, _rows, _cols, _type>& mat) {
    MAT_COPY:
        for (size_t r = 0; r != n_rows; ++r) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t c = 0; c != n_cols; ++c) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                _tryAssign(r, c, mat(r, c));
            }
        }
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
#ifdef FLAMES_PRINT_PER_MAT_COPY
        std::cout << "Mat copy!" << std::endl;
#endif
    }

    /**
     * @brief Construct a new Mat object from std::vector.
     *
     * @param vec The std::vector storing data in row major.
     * @details You can set the array partition using macro
     *          `FLAMES_MAT_PARTITION_COMPLETE` to set a complete array partition
     *          and `FLAMES_MAT_PARTITION_FACTOR` to set a block partition with the specific factor.
     */
    Mat(const std::vector<T>& vec) {
        assert(vec.size() == size() && "Initialization vector size disagrees.");
    MAT_COPY_FROM_STD_VEC:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = vec[i];
        }
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
#ifdef FLAMES_PRINT_PER_MAT_COPY
        std::cout << "Mat copy!" << std::endl;
#endif
    }

    //   public: // original private
  public:
    /**
     * @brief Construct a new Mat object from raw data pointer.
     *
     * @param ptr The raw data pointer.
     * @param opt Option if initialization.
     * @note This function should only be used by view classes and is therefore marked as private.
     */
    explicit Mat(const T* ptr, InitAfterwards opt = InitAfterwards::NONE) {
        if (opt == InitAfterwards::NONE) {
        MAT_COPY:
            for (size_t i = 0; i != size(); ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
                _data[i] = ptr[i];
            }
        } else if (opt == InitAfterwards::OPP) {
        MAT_COPY_OPP:
            for (size_t i = 0; i != size(); ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
                _data[i] = -ptr[i];
            }
        } else if (opt == InitAfterwards::TR) {
            Mat<T, n_cols, n_rows, type> tmp(ptr, InitAfterwards::NONE);
            this->t(tmp);
        }
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
#ifdef FLAMES_PRINT_PER_MAT_COPY
        std::cout << "Mat copy!" << std::endl;
#endif
    }

    explicit Mat(T* const ptr, InitAfterwards opt = InitAfterwards::NONE) {
        if (opt == InitAfterwards::NONE) {
        MAT_COPY:
            for (size_t i = 0; i != size(); ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
                _data[i] = ptr[i];
            }
        } else if (opt == InitAfterwards::OPP) {
        MAT_COPY_OPP:
            for (size_t i = 0; i != size(); ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
                _data[i] = -ptr[i];
            }
        } else if (opt == InitAfterwards::TR) {
            Mat<T, n_cols, n_rows, type> tmp(ptr, InitAfterwards::NONE);
            this->t(tmp);
        }
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = _data type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
#ifdef FLAMES_PRINT_PER_MAT_COPY
        std::cout << "Mat copy!" << std::endl;
#endif
    }

  public:
    /**
     * @brief Destroy the Mat object.
     *
     * @details So far there is nothing to do.
     */
    ~Mat() {
        // so far nothing to do
    }

    // template <typename... _unused, MatType _type = type,
    //           typename std::enable_if_t<_type == MatType::NORMAL, bool> = true>
    // inline constexpr size_t size_const() noexcept {
    //     static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::size()!");
    //     return n_rows * n_cols;
    // }

    /**
     * @brief Get the matrix size of storage.
     *
     * @details This is useful in that matrices of different MatType has the data array of different sizes.
     * @return (constexpr size_t) The data array size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Get read only element by row major index from the data array.
     *
     * @details If you are not sure which index it is
     *          operator() which takes tow index and column index is a better resort.
     * @param index The data array index in row major.
     * @return (T) The read only data element.
     */
    T operator[](size_t index) const {
        FLAMES_PRAGMA(INLINE)
        assert(index < size() && "Matrix index should be within range");
        return _data[index];
    }

    /**
     * @brief Get writeable element by row major index from the data array.
     *
     * @details If you are not sure which index it is
     *          operator() which takes tow index and column index is a better resort.
     * @param index The data array index in row major.
     * @return (T) The writeable data element.
     */
    T& operator[](size_t index) {
        FLAMES_PRAGMA(INLINE)
        assert(index < size() && "Matrix index should be within range");
        return _data[index];
    }

    /**
     * @brief Get read only data element by row index and column index.
     *
     * @param r The row index (starting from 0).
     * @param c The column index (staring from 0).
     * @return (T) The read only data element.
     */
    T operator()(size_t r, size_t c) const {
        FLAMES_PRAGMA(INLINE)
        assert(r < n_rows && "Matrix row index should be within range");
        assert(c < n_cols && "Matrix col index should be within range");
        if (type == MatType::NORMAL) {
            return _data[r * n_cols + c];
        } else if (type == MatType::DIAGONAL) {
            if (r == c) return _data[r];
            else return T(0);
        } else if (type == MatType::SCALAR) {
            if (r == c) return _data[0];
            else return T(0);
        } else if (type == MatType::UPPER) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return T(0);
        } else if (type == MatType::LOWER) {
            if (r >= c) return _data[(1 + r) * r / 2 + c];
            else return T(0);
        } else if (type == MatType::SUPPER) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 1 - 2 * r];
            else return T(0);
        } else if (type == MatType::SLOWER) {
            if (r > c) return _data[(1 + r) * r / 2 + c - r];
            else return T(0);
        } else if (type == MatType::SYM) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + r - c];
        } else if (type == MatType::ASYM) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 1 - 2 * r];
            else if (r > c) return -_data[(2 * n_cols + 1 - c) * c / 2 + r - 1 - 2 * c];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get writeable data element by row index and column index.
     *
     * @param r The row index (starting from 0).
     * @param c The column index (staring from 0).
     * @return (T) The writeable data element.
     */
    T& operator()(size_t r, size_t c) {
        FLAMES_PRAGMA(INLINE)
        assert(r < n_rows && "Matrix row index should be within range");
        assert(c < n_cols && "Matrix col index should be within range");
        if (type == MatType::NORMAL) {
            return _data[r * n_cols + c];
        } else if (type == MatType::DIAGONAL) {
            if (r == c) return _data[r];
            else assert(!"This element cannot be modified (DIAGONAL).");
        } else if (type == MatType::SCALAR) {
            assert(!"This element cannot be modified (SCALAR).");
        } else if (type == MatType::UPPER) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else assert(!"This element cannot be modified (UPPER).");
        } else if (type == MatType::LOWER) {
            if (r >= c) return _data[(1 + r) * r / 2 + c];
            else assert(!"This element cannot be modified (LOWER).");
        } else if (type == MatType::SUPPER) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 1 - 2 * r];
            else assert(!"This element cannot be modified (SUPPER).");
        } else if (type == MatType::SLOWER) {
            if (r > c) return _data[(1 + r) * r / 2 + c - r];
            else assert(!"This element cannot be modified (SLOWER).");
        } else if (type == MatType::SYM) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + r - c];
        } else if (type == MatType::ASYM) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 1 - 2 * r];
            else if (r > c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - 1 - 2 * c];

            // ATTENTION  This part needs to be perfected , missing a minus sign.
            // Because a minus sign will result in a error about reference.

            else assert(!"This element cannot be modified (ASYM).");
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
        // Just to make the compiler happy.
        return _data[0];
    }

    /**
     * @brief Set all elements of the matrix to a value.
     *
     * @param val The value.
     */
    void setValue(T val) {
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SET_VALUE_UNROLL_FACTOR)
            _data[i] = val;
        }
    }

    /**
     * @brief Set all elements of the matrix to zero.
     *
     */
    void setZero() { setValue(static_cast<T>(0)); }

    bool read(const std::string& file_name) {
#ifndef __SYNTHESIS__
        std::ifstream f(file_name);
        if (f.is_open()) {
            size_t in_rows, in_cols;
            f >> in_rows >> in_cols;
            if (n_rows != in_rows || n_cols != in_cols) return false; // dimension does not match
            std::string complex_real, mat_type;
            f >> complex_real >> mat_type;
            if (type == MatType::NORMAL && mat_type != "normal") return false;
            if (type == MatType::DIAGONAL && mat_type != "diagonal") return false;
            if (type == MatType::SCALAR && mat_type != "scalar") return false;
            if (type == MatType::UPPER && mat_type != "upper") return false;
            if (type == MatType::LOWER && mat_type != "lower") return false;
            if (type == MatType::SUPPER && mat_type != "supper") return false;
            if (type == MatType::SYM && mat_type != "sym") return false;
            if (type == MatType::ASYM && mat_type != "asym") return false;
            if (complex_real == "complex") {
                // std::string real_part, imag_part;
                // for (size_t i = 0; i != size(); ++i) {
                //     std::getline(f, real_part, ',');
                //     std::getline(f, imag_part, ',');
                //     _data[i].real(std::stod(real_part));
                //     _data[i].imag(std::stod(imag_part));
                // }
                assert(!"Read from a complex matrix is not currently supported.");
                return false;
            } else {
                std::string buf;
                for (size_t i = 0; i != size(); ++i) {
                    std::getline(f, buf, ',');
                    _data[i] = std::stod(buf);
                }
                return true;
            }
        } else {
            return false;
        }
#else
        return true;
#endif
    }

    /**
     * @brief Matrix plus matrix with same MatType.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2>
    Mat& add(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    MAT_PLUS_MAT_SAME:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            this->_data[i] = mat_L[i] + mat_R[i];
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into NORMAL.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::NORMAL && (type1 != MatType::NORMAL || type2 != MatType::NORMAL),
                               bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_NORMAL:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[i * n_cols + j] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into DIAGONAL.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::DIAGONAL && (type1 != MatType::DIAGONAL || type2 != MatType::DIAGONAL),
                               bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_LOWER:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            this->_data[i] = mat_L(i, i) + mat_R(i, i);
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into SCALAR.
     *
     * @details The result is stored to 'this'.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SCALAR && (type1 != MatType::SCALAR || type2 != MatType::SCALAR),
                               bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_SCALAR:
        this->_data[0] = mat_L(0, 0) + mat_R(0, 0);
        return *this;
    }

    /**
     * @brief Matrix plus matrix into UPPER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2,
        std::enable_if_t<type == MatType::UPPER && (type1 != MatType::UPPER || type2 != MatType::UPPER), bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_UPPER:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = i; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - i] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into LOWER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2,
        std::enable_if_t<type == MatType::LOWER && (type1 != MatType::LOWER || type2 != MatType::LOWER), bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_LOWER:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = i; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(1 + i) * i / 2 + j] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into SUPPER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SUPPER && (type1 != MatType::SUPPER || type2 != MatType::SUPPER),
                               bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_SUPPER:
        for (size_t i = 0; i != n_rows - 1; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - 2 * i - 1] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into SLOWER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SLOWER && (type1 != MatType::SLOWER || type2 != MatType::SLOWER),
                               bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_SLOWER:
        for (size_t i = 1; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols - 1; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(1 + i) * i / 2 + j - i] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into SYM.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SYM && (type1 != MatType::SYM || type2 != MatType::SYM), bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_SYM:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = i; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - i] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix plus matrix into ASYM.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_PLUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing addition in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2,
        std::enable_if_t<type == MatType::ASYM && (type1 != MatType::ASYM || type2 != MatType::ASYM), bool> = true>
    Mat& add(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_PLUS_MAT_ASYM:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_PLUS_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - i * 2 - 1] = mat_L(i, j) + mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix self plus a matrix.
     *
     * @tparam M The matrix type of the plus matrix.
     * @tparam _unused (unused)
     * @tparam T2 The plus matrix element type.
     * @tparam type2 The plus matrix MatType.
     * @param mat_R The plus matrix.
     * @return (Mat&) The addition result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& add(const M<T2, n_rows, n_cols, type2, _unused...>& mat_R) {
        FLAMES_PRAGMA(INLINE)
        return this->add(*this, mat_R);
    }

    /**
     * @brief Matrix minus matrix with same MatType.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2>
    Mat& sub(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    MAT_MINUS_MAT_SAME:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            this->_data[i] = mat_L[i] - mat_R[i];
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into NORMAL.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::NORMAL && (type1 != MatType::NORMAL || type2 != MatType::NORMAL),
                               bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_NORMAL:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[i * n_cols + j] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into DIAGONAL.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::DIAGONAL && (type1 != MatType::DIAGONAL || type2 != MatType::DIAGONAL),
                               bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_DIAGONAL:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            this->_data[i] = mat_L(i, i) - mat_R(i, i);
        }
        return *this;
    }
    /**
     * @brief Matrix minus matrix into SCALAR.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SCALAR && (type1 != MatType::SCALAR || type2 != MatType::SCALAR),
                               bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_SCALAR:
        this->_data[0] = mat_L(0, 0) - mat_R(0, 0);
        return *this;
    }

    /**
     * @brief Matrix minus matrix into UPPER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2,
        std::enable_if_t<type == MatType::UPPER && (type1 != MatType::UPPER || type2 != MatType::UPPER), bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_UPPER:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = i; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - i] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into LOWER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2,
        std::enable_if_t<type == MatType::LOWER && (type1 != MatType::LOWER || type2 != MatType::LOWER), bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_LOWER:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = i; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(1 + i) * i / 2 + j] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into SUPPER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SUPPER && (type1 != MatType::SUPPER || type2 != MatType::SUPPER),
                               bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_SUPPER:
        for (size_t i = 0; i != n_rows - 1; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - 2 * i - 1] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into SLOWER.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SLOWER && (type1 != MatType::SLOWER || type2 != MatType::SLOWER),
                               bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_SLOWER:
        for (size_t i = 1; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols - 1; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(1 + i) * i / 2 + j - i] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into SYM.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2,
              std::enable_if_t<type == MatType::SYM && (type1 != MatType::SYM || type2 != MatType::SYM), bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_SYM:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = i; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - i] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix minus matrix into ASYM.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix matrix MatType.
     * @tparam type2 The right matrix matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2,
        std::enable_if_t<type == MatType::ASYM && (type1 != MatType::ASYM || type2 != MatType::ASYM), bool> = true>
    Mat& sub(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
             const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    MAT_MINUS_MAT_ASYM:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                this->_data[(2 * n_cols + 1 - i) * i / 2 + j - i * 2 - 1] = mat_L(i, j) - mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Matrix self minus a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_MINUS_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing subtraction in parallel.
     * @tparam M The matrix type of the minus matrix.
     * @tparam _unused (unused)
     * @tparam T2 The minus matrix type.
     * @tparam type2 The minus matrix MatType.
     * @param mat_R The minus matrix.
     * @return (Mat&) The subtraction result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& sub(const M<T2, n_rows, n_cols, type2, _unused...>& mat_R) {
        FLAMES_PRAGMA(INLINE)
        // return this->sub(*this, mat_R);
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            _data[i] -= mat_R[i];
        }
        return *this;
    }

    /**
     * @brief Matrix times a scalar.
     *
     * @details This scalar is C++ arithmetic types, like double and int.
     *          The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M The matrix type.
     * @tparam _unused (unused)
     * @tparam ScalarT The scalar type.
     * @tparam T2 The matrix element type.
     * @param mat The matrix.
     * @param s The scalar.(The scalar here is a number. NOT the SCALAR which is a square matrix,despite their functions
     * are similar)
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename ScalarT,
        typename T2,
        std::enable_if_t<std::is_arithmetic<std::remove_cv_t<std::remove_reference_t<ScalarT>>>::value, bool> = true>
    Mat& mul(const M<T2, n_rows, n_cols, type, _unused...>& mat, ScalarT s) {
    MAT_SCALAR_TIMES:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] = mat[i] * s;
        }
        return *this;
    }

    /**
     * @brief Matrix times a complex number.
     *
     * @details This scalar is std::complex.
     *          The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M The matrix type.
     * @tparam _unused (unused)
     * @tparam ScalarT The scalar type.
     * @tparam T2 The matrix element type.
     * @param mat The matrix.
     * @param s The complex number.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename ScalarT,
        typename T2,
        std::enable_if_t<std::is_arithmetic<std::remove_cv_t<std::remove_reference_t<ScalarT>>>::value, bool> = true>
    Mat& mul(const M<T2, n_rows, n_cols, type, _unused...>& mat, std::complex<ScalarT> s) {
    MAT_SCALAR_TIMES:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] = mat[i] * s;
        }
        return *this;
    }

    /**
     * @brief Matrix times ap_int integer.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M The matrix type.
     * @tparam _unused (unused)
     * @tparam AP_W ap_int W param.
     * @tparam T2 The matrix element type.
     * @param mat The matrix.
     * @param s The integer.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, int AP_W,
              typename T2>
    Mat& mul(const M<T2, n_rows, n_cols, type, _unused...>& mat, ap_int<AP_W> s) {
    MAT_SCALAR_TIMES:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] = mat[i] * s;
        }
        return *this;
    }

    /**
     * @brief Matrix times ap_fixed float.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M The matrix type.
     * @tparam _unused (unused)
     * @tparam AP_W ap_fixed W param.
     * @tparam AP_I ap_fixed I param.
     * @tparam AP_Q ap_fixed Q param.
     * @tparam AP_O ap_fixed O param.
     * @tparam AP_N ap_fixed N param.
     * @tparam T2 The matrix element type.
     * @param mat The matrix.
     * @param s The float.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, int AP_W, int AP_I,
              ap_q_mode AP_Q, ap_o_mode AP_O, int AP_N, typename T2>
    Mat& mul(const M<T2, n_rows, n_cols, type, _unused...>& mat, ap_fixed<AP_W, AP_I, AP_Q, AP_O, AP_N> s) {
    MAT_SCALAR_TIMES:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] = mat[i] * s;
        }
        return *this;
    }

    /**
     * @brief Matrix self multiply a scalar.
     *
     * @details This scalar is C++ arithmetic types, like double and int.
     *          The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam ScalarT The scalar type.
     * @param s The scalar.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        typename ScalarT,
        std::enable_if_t<std::is_arithmetic<std::remove_cv_t<std::remove_reference_t<ScalarT>>>::value, bool> = true>
    Mat& mul(ScalarT s) {
    MAT_SCALAR_TIMES_SELF:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] *= s;
        }
        return *this;
    }

    /**
     * @brief Matrix self multiply a complex number.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam ScalarT The scalar type.
     * @param s The complex number.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        typename ScalarT,
        std::enable_if_t<std::is_arithmetic<std::remove_cv_t<std::remove_reference_t<ScalarT>>>::value, bool> = true>
    Mat& mul(std::complex<ScalarT> s) {
    MAT_SCALAR_TIMES_SELF:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] *= s;
        }
        return *this;
    }

    /**
     * @brief Matrix self multiply ap_int integer.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam AP_W ap_int W param.
     * @param s The integer.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <int AP_W>
    Mat& mul(ap_int<AP_W> s) {
    MAT_SCALAR_TIMES_SELF:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] *= s;
        }
        return *this;
    }

    /**
     * @brief Matrix self multiply ap_fixed float.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam AP_W ap_int W param.
     * @tparam AP_I ap_int I param.
     * @tparam AP_Q ap_int Q param.
     * @tparam AP_O ap_int O param.
     * @tparam AP_N ap_int N param.
     * @param s The float.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <int AP_W, int AP_I, ap_q_mode AP_Q, ap_o_mode AP_O, int AP_N>
    Mat& mul(ap_fixed<AP_W, AP_I, AP_Q, AP_O, AP_N> s) {
    MAT_SCALAR_TIMES_SELF:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            _data[i] *= s;
        }
        return *this;
    }

    /**
     * @brief General matrix matrix multiplication.(Including SYM or NORMAL matrix times a SYM or NORMAL matrix )
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     *          (This is implemented of using systolic array.)
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(!(std::is_same<T1, bool>::value) && !(std::is_same<T2, bool>::value)) &&
                                   ((type1 == MatType::NORMAL && type2 == MatType::NORMAL) ||
                                    (type1 == MatType::NORMAL && type2 == MatType::SYM) ||
                                    (type1 == MatType::SYM && type2 == MatType::NORMAL) ||
                                    (type1 == MatType::SYM && type2 == MatType::SYM)),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    GEMM:
        for (size_t i = 0; i != comm; ++i) {
        GEMM_r:
            for (size_t r = 0; r != n_rows; ++r) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            // FLAMES_PRAGMA(DEPENDENCE variable=mat_L inter false)
            // FLAMES_PRAGMA(DEPENDENCE variable=mat_L intra false)
            // FLAMES_PRAGMA(DEPENDENCE variable=mat_R inter false)
            // FLAMES_PRAGMA(DEPENDENCE variable=mat_R intra false)
            GEMM_c:
                for (size_t c = 0; c != n_cols; ++c) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    if (i == 0) (*this)(r, c) = T(0); // initialize
                    (*this)(r, c) += mat_L(r, i) * mat_R(i, c);
                }
            }
        }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @tparam comm The column number of the right matrix (should be the same as n_cols of 'this').
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_NORAML_OR_SYM_TIMES_MAT_ASYM:
        for (size_t i = 0; i != comm; ++i) {
        GEMM_r:
            for (size_t r = 0; r != n_rows; ++r) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            GEMM_c:
                for (size_t c = 0; c != n_cols; ++c) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    if (i == 0) (*this)(r, c) = T(0); // initialize
                    if (i != c) (*this)(r, c) += mat_L(r, i) * mat_R(i, c);
                }
            }
        }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix  times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::DIAGONAL),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_NORAML_OR_SYM_TIMES_DIAG:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat_L(i, j) * mat_R[j];
            }
        }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::SCALAR),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_NORAML_OR_MAT_SYM_TIMES_SCAL:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t j = 0; j != n_rows; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat_L(i, j) * mat_R[0];
            }
        }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix times a upper triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
        static const size_t i[288] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7
        };
        static const size_t c[288] = {
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7
        };
    MAT_NORAML_OR_MAT_SYM_TIMES_UPPER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix times a lower triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
        static const size_t i[288] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7
        };
        static const size_t c[288] = {
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0
        };
    MAT_NORAML_OR_MAT_SYM_TIMES_LOWER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix times a strict upper triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::SUPPER),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0,
            1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1,
            2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1,
            2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0,
            1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4,
            5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2,
            3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t c[224] = {
            1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3,
            3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4,
            4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5,
            5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
            6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7,
            7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_NORAML_OR_MAT_SYM_TIMES_SUPPER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t i = 0; i != n_rows; i++) { (*this)(i, 0) = 0; }
        return *this;
    }

    /**
     * @brief Normal matrix or symmetric matrix times a strict lower triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<((type1 == MatType::NORMAL || type1 == MatType::SYM) && type2 == MatType::SLOWER),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5,
            6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5,
            6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4,
            5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2,
            3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6,
            7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3,
            4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[224] = {
            6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4,
            4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3,
            3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2,
            2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1,
            1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0,
            0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_NORAML_OR_MAT_SYM_TIMES_SLOWER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n] + 1) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t i = 0; i != n_rows; i++) { (*this)(i, n_rows - 1) = 0; }
        return *this;
    }
    /**
     * @brief Diagonal matrix times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_DIAG_TIMES_MAT_DIAG:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(i, i) = mat_L[i] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_DIAG_TIMES_MAT_SCAL:
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[i] * mat_R[0];
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times a upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[36] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t i[36] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[36] = {
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_DIAG_TIMES_UPPER:
        for (size_t n = 0; n != n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times a lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[36] = {
            0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7,
        };
        static const size_t i[36] = {
            0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7,
        };
        static const size_t c[36] = {
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7,
        };
    MAT_DIAG_TIMES_LOWER:
        for (size_t n = 0; n != n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times a strict upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::SUPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[28] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6
        };
        static const size_t i[28] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6
        };
        static const size_t c[28] = {
            1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7
        };
    MAT_DIAG_TIMES_SUPPER:
        for (size_t n = 0; n != n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times a strict lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::SLOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[28] = {
            1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7
        };
        static const size_t i[28] = {
            1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7
        };
        static const size_t c[28] = {
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6
        };
    MAT_DIAG_TIMES_SLOWER:
        for (size_t n = 0; n != n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && (type2 == MatType::NORMAL || type2 == MatType::SYM)),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_DIAG_TIMES_MAT_NORMAL_or_SYM:
        for (size_t j = 0; j != n_cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t i = 0; i != n_rows; ++i) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat_L(i, i) * mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief Diagonal matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::DIAGONAL && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_DIAG_TIMES_MAT_ASYM:
        for (size_t j = 0; j != n_cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t i = 0; i != n_rows; ++i) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                if (i != j) (*this)(i, j) = mat_L(i, i) * mat_R(i, j);
                else (*this)(i, j) = T(0);
            }
        }
        return *this;
    }

    /**
     * @brief Scalar matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SCALAR && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_SCAL:
        _data[0] = mat_L[0] * mat_R[0];
        return *this;
    }

    /**
     * @brief  Scalar matrix times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SCALAR && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_DIAG:
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[0] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief  Scalar matrix times a upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type2 == MatType::UPPER && type1 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_UPPER:
        for (size_t i = 0; i != (1 + n_rows) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[0] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief  Scalar matrix times a lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type2 == MatType::LOWER && type1 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_LOWER:
        for (size_t i = 0; i != (1 + n_rows) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[0] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief  Scalar matrix times a strict upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type2 == MatType::SUPPER && type1 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_SUPPER:
        for (size_t i = 0; i != (n_rows - 1) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[0] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief  Scalar matrix times a strict lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type2 == MatType::SLOWER && type1 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_SLOWER:
        for (size_t i = 0; i != (n_rows - 1) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[0] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief  Scalar matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type2 == MatType::ASYM && type1 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_SLOWER:
        for (size_t i = 1; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                if (i != j) (*this)(i, j) = mat_L[0] * mat_R(i, j);
                else (*this)(i, j) = T(0);
            }
        }
        return *this;
    }

    /**
     * @brief  Scalar matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<((type2 == MatType::SYM || type2 == NORMAL) && type1 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SCAL_TIMES_MAT_NORMAL_OR_MAT_SYM:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat_L[0] * mat_R(i, j);
            }
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_UPPER_TIMES_MAT_SCALAR:
        for (size_t i = 0; i != (n_rows + 1) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[i] * mat_R[0];
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[36] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                                      2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 7 };
        static const size_t i[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4,
                                      5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7 };
        static const size_t c[36] = { 0, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 2, 3, 4,
                                      5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7 };
    MAT_UPPER_TIMES_DIAG:
        for (size_t n = 0; n != n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[120] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 0, 0,
                                       0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 1, 1, 1, 1,
                                       2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2,
                                       2, 2, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
                                       2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
        static const size_t i[120] = { 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 0, 1,
                                       2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6, 5, 6, 7, 0, 1, 2, 3, 1, 2, 3, 4,
                                       2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3,
                                       4, 5, 6, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 2, 3, 4, 5,
                                       6, 7, 0, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7 };
        static const size_t c[120] = { 0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 2, 2,
                                       2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 3, 3, 3, 3, 4, 4, 4, 4,
                                       5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6,
                                       6, 6, 6, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7,
                                       7, 7, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 };
    MAT_UPPER_TIMES_UPPER:
        for (size_t n = 0; n != n_rows * (n_rows + 1) * (n_rows + 2) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[120] = {
            0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 0, 0, 0, 1, 1, 1, 2, 2,
            2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
            0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[120] = {
            0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 0, 1, 2, 1, 2, 3, 2, 3,
            4, 3, 4, 5, 4, 5, 6, 5, 6, 7, 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7,
            0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4,
            5, 6, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[120] = {
            0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 2, 2, 2, 3, 3, 3, 4, 4,
            4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
            4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,
            6, 6, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_UPPER_TIMES_MAT_LOWER:
        for (size_t n = 0; n != n_rows * (n_rows + 1) * (n_rows + 2) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a strict upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::SUPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[84] = {
            0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2,
            3, 3, 3, 4, 4, 4, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[84] = {
            0, 1, 2, 3, 4, 5, 6, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 0, 1, 2, 1, 2, 3, 2, 3, 4,
            3, 4, 5, 4, 5, 6, 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 0, 1, 2, 3, 4, 1,
            2, 3, 4, 5, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[84] = {
            1, 2, 3, 4, 5, 6, 7, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 3, 3, 3, 4, 4, 4, 5, 5, 5,
            6, 6, 6, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_UPPER_TIMES_MAT_SUPPER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a strict lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::SLOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[168] = {
            7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2,
            1, 1, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[168] = {
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7,
            6, 7, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7,
            4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 3, 4,
            5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7,
            2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[168] = {
            0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0,
            0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_UPPER_TIMES_MAT_SLOWER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 3; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n] || i[n] == (c[n] + 1)) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::UPPER && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[288] = {
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t c[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[288] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
        };
    MAT_UPPER_TIMES_MAT_ASYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n]) (*this)(r[n], c[n]) = 0; // initialize
            if (i[n] != c[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Upper triangle matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<(type1 == MatType::UPPER && (type2 == MatType::SYM || type2 == MatType::NORMAL)), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[288] = {
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t c[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[288] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
        };
    MAT_LOWER_TIMES_NORAML_OR_MAT_SYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_LOWER_TIMES_MAT_SCALAR:
        for (size_t i = 0; i != (n_rows + 1) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[i] * mat_R[0];
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[36] = {
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[36] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[36] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
    MAT_LOWER_TIMES_DIAG:
        for (size_t n = 0; n != n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[204] = {
            7, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 7, 7, 7, 6, 6, 6, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[204] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
            1, 0, 1, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
            0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0,
            1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0,
            1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[204] = {
            0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4,
            4, 5, 5, 6, 6, 7, 7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5,
            6, 6, 6, 7, 7, 7, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
            6, 6, 6, 6, 7, 7, 7, 7, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7,
            7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_LOWER_TIMES_MAT_UPPER:
        for (size_t n = 0; n != (2 * n_rows * n_rows * n_rows + 3 * n_rows * n_rows + n_rows) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[120] = {
            0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 0, 0, 0, 1, 1, 1, 2, 2,
            2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4,
            0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[120] = {
            0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 0, 1, 2, 1, 2, 3, 2, 3,
            4, 3, 4, 5, 4, 5, 6, 5, 6, 7, 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7,
            0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4,
            5, 6, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[120] = {
            0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 2, 2, 2, 3, 3, 3, 4, 4,
            4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
            4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6,
            6, 6, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_LOWER_TIMES_MAT_LOWER:
        for (size_t n = 0; n != n_rows * (n_rows + 1) * (n_rows + 2) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a strict upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.bb
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::SUPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[168] = {
            0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
            6, 6, 7, 7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[168] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
            0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1,
            2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
            0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t c[168] = {
            7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 7, 7, 7, 6, 6, 6, 5, 5, 5, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7,
            7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_LOWER_TIMES_MAT_SUPPER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 3; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t r = 0; r != n_cols; r++) { (*this)(r, 0) = T(0); }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a strict lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */

    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::SLOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[84] = {
            0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2,
            3, 3, 3, 4, 4, 4, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[84] = {
            1, 2, 3, 4, 5, 6, 7, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 1, 2, 3, 2, 3, 4, 3, 4, 5,
            4, 5, 6, 5, 6, 7, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 1, 2, 3, 4, 5, 2,
            3, 4, 5, 6, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[84] = {
            1, 2, 3, 4, 5, 6, 7, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 3, 3, 3, 4, 4, 4, 5, 5, 5,
            6, 6, 6, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_LOWER_TIMES_MAT_SLOWER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n] + 1) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::LOWER && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[288] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[288] = {
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_LOWER_TIMES_MAT_ASYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            if (i[n] != c[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Lower triangle matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<(type1 == MatType::LOWER && (type2 == MatType::SYM || type2 == MatType::NORMAL)), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[288] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[288] = {
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_LOWER_TIMES_MAT_NORMAL_OR_MAT_SYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Strict upper triangle matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SUPPER_TIMES_MAT_SCALAR:
        for (size_t i = 0; i != (n_rows - 1) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[i] * mat_R[0];
        }
        return *this;
    }

    /**
     * @brief   Strict upper triangle matrix times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[28] = {
            1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7,
        };
        static const size_t i[28] = {
            1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 4, 5, 6, 7, 5, 6, 7, 6, 7, 7,
        };
        static const size_t r[28] = {
            0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6,
        };
    MAT_SUPPER_TIMES_MAT_DIAG:
        for (size_t n = 0; n != n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Strict upper triangle matrix times a upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[84] = {
            0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2,
            3, 3, 3, 4, 4, 4, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[84] = {
            1, 2, 3, 4, 5, 6, 7, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 1, 2, 3, 2, 3, 4, 3, 4, 5,
            4, 5, 6, 5, 6, 7, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 1, 2, 3, 4, 5, 2,
            3, 4, 5, 6, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[84] = {
            1, 2, 3, 4, 5, 6, 7, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 3, 3, 3, 4, 4, 4, 5, 5, 5,
            6, 6, 6, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SUPPER_TIMES_MAT_UPPER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n] + 1) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief   Strict upper triangle matrix times a lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[168] = {
            7, 7, 7, 7, 7, 7, 7, 6, 5, 4, 3, 2, 1, 0, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2,
            1, 1, 0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[168] = {
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7,
            6, 7, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7,
            4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 3, 4,
            5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7,
            2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[168] = {
            0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0,
            0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_UPPER_TIMES_MAT_SLOWER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 3; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n] || i[n] == (r[n] + 1)) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t c = 0; c != n_cols; c++) { (*this)(n_rows - 1, c) = T(0); }

        return *this;
    }

    /**
     * @brief   Strict upper triangle matrix times a strict upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && type2 == MatType::SUPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[56] = {
            0, 1, 2, 3, 4, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
            0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[56] = {
            1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6,
            1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
        };
        static const size_t c[56] = {
            2, 3, 4, 5, 6, 7, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7,
            5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SUPPER_TIMES_MAT_SUPPER:
        for (size_t n = 0; n != ((n_rows - 1) * (n_rows - 1) * (n_rows - 1) - (n_rows - 1)) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 1 + r[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t i = 0; i != n_rows - 1; i++) { (*this)(i, i + 1) = T(0); }
        return *this;
    }

    /**
     * @brief  Strict upper triangle times a strict lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<!(type1 == MatType::SUPPER && n_rows == 1) && !(type2 == MatType::SLOWER && cols_ == 1) &&
                             (type1 == MatType::SUPPER && type2 == MatType::SLOWER),
                         bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[140] = {
            6, 6, 6, 6, 6, 6, 6, 5, 4, 3, 2, 1, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 0, 0,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2, 2, 2, 1, 1, 1, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[140] = {
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7, 6, 7,
            5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7,
            4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 3, 4, 5, 6, 7,
            3, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[140] = {
            0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 0, 0, 1, 1, 1, 1,
            2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_SUPPER_TIMES_MAT_SLOWER:
        for (size_t n = 0;
             n != (2 * (n_rows - 1) * (n_rows - 1) * (n_rows - 1) + 3 * (n_rows - 1) * (n_rows - 1) + (n_rows - 1)) / 6;
             n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n] + 1 || i[n] == c[n] + 1) (*this)(r[n], c[n]) = T(0); // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t r = 0; r != n_rows; r++) { (*this)(r, n_cols - 1) = T(0); }
        for (size_t c = 0; c != n_cols - 1; c++) { (*this)(n_rows - 1, c) = T(0); }

        return *this;
    }

    /**
     * @brief  Strict upper triangle matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5,
            6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5,
            6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4,
            5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2,
            3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6,
            7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3,
            4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[224] = {
            6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4,
            4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3,
            3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2,
            2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1,
            1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0,
            0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_SUPPER_TIMES_MAT_ASYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n] + 1) (*this)(r[n], c[n]) = T(0); // initialize
            if (i[n] != c[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t c = 0; c != n_cols; c++) { (*this)(n_rows - 1, c) = T(0); }

        return *this;
    }

    /**
     * @brief  Strict upper triangle matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SUPPER && (type2 == MatType::SYM || type2 == MatType::NORMAL)),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5,
            6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5,
            6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4,
            5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2,
            3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6,
            7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3,
            4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t r[224] = {
            6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4,
            4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3,
            3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2,
            2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1,
            1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0,
            0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_SUPPER_TIMES_MAT_NORMAL_OR_SYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == r[n] + 1) (*this)(r[n], c[n]) = T(0); // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t c = 0; c != n_cols; c++) { (*this)(n_rows - 1, c) = T(0); }
        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_SLOWER_TIMES_MAT_SCALAR:
        for (size_t i = 0; i != (n_rows - 1) * n_rows / 2; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            _data[i] = mat_L[i] * mat_R[0];
        }
        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[28] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t i[28] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[28] = {
            1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_DIAG:
        for (size_t n = 0; n != n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            (*this)(r[n], c[n]) = mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[168] = {
            0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5,
            6, 6, 7, 7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[168] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
            0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1,
            2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
            0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[168] = {
            7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 7, 7, 7, 6, 6, 6, 5, 5, 5, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7,
            7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_UPPER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 3; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t c = 0; c != n_cols; c++) { (*this)(0, c) = T(0); }

        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[84] = {
            0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 0, 1, 1, 1, 2, 2, 2,
            3, 3, 3, 4, 4, 4, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 0, 1,
            1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[84] = {
            0, 1, 2, 3, 4, 5, 6, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 0, 1, 2, 1, 2, 3, 2, 3, 4,
            3, 4, 5, 4, 5, 6, 0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 0, 1, 2, 3, 4, 1,
            2, 3, 4, 5, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[84] = {
            1, 2, 3, 4, 5, 6, 7, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 3, 3, 3, 4, 4, 4, 5, 5, 5,
            6, 6, 6, 7, 7, 7, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_LOWER:
        for (size_t n = 0; n != (n_rows * n_rows * n_rows - n_rows) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a strict upper triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.bb
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::SUPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[140] = {
            1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[140] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
            0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3,
            0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4,
            0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[140] = {
            7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1, 1, 7, 7, 6, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            7, 7, 7, 6, 6, 6, 5, 5, 5, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 7, 7, 7, 6, 6, 6, 6,
            5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_SUPPER:
        for (size_t n = 0;
             n != (2 * (n_rows - 1) * (n_rows - 1) * (n_rows - 1) + 3 * (n_rows - 1) * (n_rows - 1) + (n_rows - 1)) / 6;
             n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t r = 0; r != n_rows; r++) { (*this)(r, 0) = T(0); }
        for (size_t c = 1; c != n_cols; c++) { (*this)(0, c) = T(0); }

        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a strict lower triangle matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */

    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::SLOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[56] = {
            0, 1, 2, 3, 4, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
            0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
        };
        static const size_t i[56] = {
            1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6,
            1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 1, 2, 3, 4, 5, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[56] = {
            2, 3, 4, 5, 6, 7, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7,
            5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_SLOWER:
        for (size_t n = 0; n != ((n_rows - 1) * (n_rows - 1) * (n_rows - 1) - (n_rows - 1)) / 6; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 1 + c[n]) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t i = 0; i != n_rows - 1; i++) { (*this)(i, i - 1) = T(0); }

        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0,
            1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1,
            2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1,
            2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0,
            1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4,
            5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2,
            3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[224] = {
            1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3,
            3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4,
            4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5,
            5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
            6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7,
            7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_ASYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            if (i[n] != c[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t c = 0; c != n_cols; ++c) { (*this)(0, c) = T(0); }
        return *this;
    }

    /**
     * @brief  Strict lower triangle matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::SLOWER && (type2 == MatType::SYM || type2 == MatType::NORMAL)),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t c[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0,
            1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1,
            2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1,
            2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0,
            1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4,
            5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2,
            3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t r[224] = {
            1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3,
            3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4,
            4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5,
            5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
            6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7,
            7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_SLOWER_TIMES_MAT_NORMAL_OR_MAT_SYM:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t c = 0; c != n_cols; ++c) { (*this)(0, c) = T(0); }
        return *this;
    }

    /**
     * @brief Anti-symmetric matrix  times a diagonal matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::DIAGONAL), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_ASYM_TIMES_MAT_DIAG:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t j = 0; j != n_rows; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                if (i != j) (*this)(i, j) = mat_L(i, j) * mat_R[j];
            }
        }
        for (size_t c = 0; c != n_cols; ++c) { (*this)(c, c) = T(0); }
        return *this;
    }

    /**
     * @brief Anti-symmetric matrix times a scalar matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::SCALAR), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_ASYM_TIMES_MAT_SCAL:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_rows; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                if (i != j) (*this)(i, j) = mat_L(i, j) * mat_R[0];
            }
        }
        return *this;
    }

    /**
     * @brief Anti-symmetric times a upper triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::UPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[288] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[288] = {
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
            0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_ASYM_TIMES_MAT_UPPER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            if (r[n] != i[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Anti-symmetric matrix times a lower triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::LOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[288] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[288] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[288] = {
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_ASYM_TIMES_MAT_LOWER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows + 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n]) (*this)(r[n], c[n]) = 0; // initialize
            if (r[n] != i[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        return *this;
    }

    /**
     * @brief Anti-symmetric matrix times a strict upper triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::SUPPER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0,
            1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1,
            2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1,
            2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0,
            1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4,
            5, 0, 1, 2, 3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2,
            3, 4, 5, 6, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
        };
        static const size_t c[224] = {
            1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3,
            3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4,
            4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5,
            5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
            6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6,
            6, 7, 7, 7, 7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7,
            7, 7, 7, 7, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7,
        };
    MAT_ASYM_TIMES_SUPPER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == 0) (*this)(r[n], c[n]) = 0; // initialize
            if (r[n] != i[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t r = 0; r != n_rows; r++) { (*this)(r, 0) = T(0); }

        return *this;
    }

    /**
     * @brief Anti-symmetric matrix times a strict lower triangular matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::SLOWER), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        static const size_t r[224] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
            6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
        };
        static const size_t i[224] = {
            7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5,
            6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5,
            6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4,
            5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2,
            3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6,
            7, 1, 2, 3, 4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3,
            4, 5, 6, 7, 7, 6, 7, 5, 6, 7, 4, 5, 6, 7, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
        };
        static const size_t c[224] = {
            6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4,
            4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3,
            3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2,
            2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1,
            1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
            1, 0, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0,
            0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
        };
    MAT_ASYM_TIMES_SLOWER:
        for (size_t n = 0; n != n_rows * n_rows * (n_rows - 1) / 2; n++) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
            if (i[n] == c[n] + 1) (*this)(r[n], c[n]) = 0; // initialize
            if (r[n] != i[n]) (*this)(r[n], c[n]) += mat_L(r[n], i[n]) * mat_R(i[n], c[n]);
        }
        for (size_t r = 0; r != n_rows; r++) { (*this)(r, n_cols - 1) = T(0); }

        return *this;
    }

    /**
     * @brief Anti-symmetric matrix times a normal matrix or a symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     *          (This is implemented of using systolic array.)
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<(type1 == MatType::ASYM && (type2 == MatType::NORMAL || type2 == MatType::SYM)), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_ASYM_TIMES_MAT_NORMAL:
        for (size_t i = 0; i != comm; ++i) {
        GEMM_r:
            for (size_t r = 0; r != n_rows; ++r) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            GEMM_c:
                for (size_t c = 0; c != n_cols; ++c) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    if (i == 0) (*this)(r, c) = T(0); // initialize
                    if (r != i) (*this)(r, c) += mat_L(r, i) * mat_R(i, c);
                }
            }
        }
        return *this;
    }

    /**
     * @brief Anti-symmetric matrix times an anti-symmetric matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @tparam comm The column number of the right matrix (should be the same as n_cols of 'this').
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(type1 == MatType::ASYM && type2 == MatType::ASYM), bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    MAT_ASYM_TIMES_MAT_ASYM:
    GEMM_r:
        for (size_t r = 0; r != n_rows; ++r) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
        GEMM_c:
            for (size_t c = 0; c != n_cols; ++c) {
                for (size_t i = 0; i != comm; ++i) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    if (i == 0) (*this)(r, c) = T(0); // initialize
                    if (i != c && i != r) (*this)(r, c) += mat_L(r, i) * mat_R(i, c);
                }
            }
        }
        return *this;
    }

    /**
     * @brief Bool matrix times a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     *          (This is implemented of using systolic array.)
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type (should be bool).
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <
        template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
        template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1, typename T2,
        MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
        std::enable_if_t<(std::is_same<T1, bool>::value) && ((type1 == MatType::NORMAL && type2 == MatType::NORMAL) ||
                                                             (type1 == MatType::NORMAL && type2 == MatType::SYM) ||
                                                             (type1 == MatType::SYM && type2 == MatType::NORMAL) ||
                                                             (type1 == MatType::SYM && type2 == MatType::SYM)),
                         bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
        this->setZero();
        for (size_t c = 0; c < n_cols; ++c) {
            for (size_t i = 0; i < comm; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TIMES_UNROLL_FACTOR)
                for (size_t r = 0; r < n_rows; ++r) {
                    FLAMES_PRAGMA(UNROLL)
                    (*this)(r, c) += mat_L(r, i) * mat_R(i, c);
                }
            }
        }

        return *this;
    }

    /**
     * @brief  Matrix times a bool matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure `FLAMES_MAT_TIMES_UNROLL_FACTOR`
     *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
     *          (This is implemented of using systolic array.)
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type .
     * @tparam T2 The right matrix element type (should be bool).
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam rows_ The row number of the left matrix (should be the same as rows of 'this').
     * @tparam cols_ The column number of the right matrix (should be the same as n_cols of 'this').
     * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t rows_, size_t cols_, size_t comm,
              std::enable_if_t<(std::is_same<T2, bool>::value && !std::is_same<T1, bool>::value) &&
                                   ((type1 == MatType::NORMAL && type2 == MatType::NORMAL) ||
                                    (type1 == MatType::NORMAL && type2 == MatType::SYM) ||
                                    (type1 == MatType::SYM && type2 == MatType::NORMAL) ||
                                    (type1 == MatType::SYM && type2 == MatType::SYM)),
                               bool> = true>
    Mat& mul(const M1<T1, rows_, comm, type1, _unused1...>& mat_L,
             const M2<T2, comm, cols_, type2, _unused2...>& mat_R) {
        FLAMES_PRAGMA(INLINE off)
        static_assert(n_rows == rows_, "Matrix dimension should meet.");
        static_assert(n_cols == cols_, "Matrix dimension should meet.");
    GEMM_r:
        for (size_t r = 0; r != n_rows; ++r) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
        GEMM_c:
            for (size_t c = 0; c != n_cols; ++c) {
            GEMM:
                for (size_t i = 0; i != comm; ++i) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    if (i == 0) (*this)(r, c) = T(0); // initialize
                    if (mat_R(i, c) != 0) (*this)(r, c) += mat_L(r, i);
                }
            }
        }
        return *this;
    }

    /**
     * @brief Element-wise product of two matrices.
     *
     * @details You can configure the macro `FLAMES_MAT_EMUL_UNROLL_FACTOR` to determine the parallelism.
     * @note It now only supports element-wise product of two matrices of the same dimension and MatType.
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam rows_ The number of rows.
     * @tparam cols_ The number of columns.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @return (Mat&) The element-wise product result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, size_t rows_, size_t cols_, MatType type1, MatType type2,
              std::enable_if_t<rows_ == n_rows && cols_ == n_cols && type1 == type && type2 == type, bool> = true>
    Mat& emul(const M1<T1, rows_, cols_, type1, _unused1...>& mat_L,
              const M2<T2, rows_, cols_, type2, _unused2...>& mat_R) {
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_EMUL_UNROLL_FACTOR)
            _data[i] = mat_L[i] * mat_R[i];
        }
        return *this;
    }

    /**
     * @brief Take a column of a matrix by index.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @param c The column index.
     * @return (Mat&) The certain column(a column vector) (a reference to 'this') .
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    Mat& col(const M<T2, rows_, cols_, type2, _unused...>& mat, size_t c) {
        assert(c < n_cols && "Take the specific col by index requires 'The index should be smaller than the number of "
                             "the matrix's columns.'.");
        static_assert(size() == rows_, "Element number should be rows_ in Mat::col(mat, index).");
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = mat(i, c);
        }
        return *this;
    }

    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    void col(size_t c, const M<T2, rows_, cols_, type2, _unused...>& mat) {
        assert(c < n_cols && "Take the specific col by index requires 'The index should be smaller than the number of "
                             "the matrix's columns.'.");
        assert(mat.size() == n_rows && "Element number should be n_rows in Mat::col(index, mat).");
        for (int i = 0; i != mat.size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            (*this)(i, c) = mat[i];
        }
    }

    /**
     * @brief Take a column of a matrix by index and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .col(Mat, int) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .col_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @param c The column index.
     * @return (Mat<T, n_rows, 1, MatType::NORMAL>)(column vector) The certain column vector copy.
     */
    Mat<T, n_rows, 1, MatType::NORMAL> col(size_t c) const {
        assert(c < n_cols && "Take the specific col by index requires 'The index should be smaller than the number of "
                             "the matrix's columns.'.");
        Mat<T, n_rows, 1, MatType::NORMAL> mat;
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat(i, 0) = (*this)(i, c);
        }
        return mat;
    }

    /**
     * @brief Take a column of a matrix by index as a read only view.
     *
     * @return (MatViewCol<Index, T, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only a column vector view.
     */
    template <size_t Index>
    MatViewCol<Index, T, n_rows, n_cols, MatType::NORMAL, MType<type>> Col_() const {
        static_assert(Index > int(0), "Take the specific col by index requires 'The index can't be smaller than 0'.");
        static_assert(Index < n_cols, "Take the specific col by index requires 'The index should be smaller than the "
                                      "number of the matrix's columns.'.");
        return *this;
    }

    /**
     * @brief Take a row of a matrix by index.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The matrix MatType.
     * @param mat The original matrix.
     * @param r The row index.
     * @return (Mat&) The certain column(a row vector) (a reference to 'this') .
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    Mat& row(const M<T2, rows_, cols_, type2, _unused...>& mat, size_t r) {
        assert(r < rows_ && "row should be smaller than rows_ of mat in Mat::row.");
        static_assert(size() == cols_, "The number of the Matrices' n_cols should meet.");
        for (size_t j = 0; j != size(); ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[j] = mat(r, j);
        }
        return *this;
    }

    /**
     * @brief Take a row of a matrix by index and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .row(Mat, int) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .row_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @param r The row index.
     * @return (Mat<T, 1, n_cols, MatType::NORMAL>)(row vector) The certain row vector copy.
     */
    Mat<T, 1, n_cols, MatType::NORMAL> row(size_t r) const {
        assert(r < n_rows && "Take the specific row by index requires 'The index should be smaller than the number of "
                             "the matrix's rows.'.");
        Mat<T, 1, n_cols, MatType::NORMAL> mat;
        for (size_t j = 0; j != n_cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat(0, j) = (*this)(r, j);
        }
        return mat;
    }

    /**
     * @brief Take a row of a matrix by index as a read only view.
     *
     * @return (MatViewRow<Index, T, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only a row vector view.
     */
    template <size_t Index>
    MatViewRow<Index, T, n_rows, n_cols, MatType::NORMAL, MType<type>> Row_() const {
        static_assert(Index > int(0), "Take the specific row by index requires 'The index can't be smaller than 0'.");
        static_assert(Index < n_rows, "Take the specific row by index requires 'The index should be smaller than the "
                                      "number of the matrix's rows.'.");
        return *this;
    }

    /**
     * @brief Take successive columns of a matrix by index.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The matrix MatType.
     * @param mat The original matrix.
     * @param first_col The first column index.
     * @return (Mat&) The successive columns(a normal matrix) (a reference to 'this') .
     */

    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    Mat& cols(const M<T2, rows_, cols_, type2, _unused...>& mat, size_t first_col) {
        assert(first_col < n_cols && "Take the successive n_cols by index requires 'The indexes should be smaller than "
                                     "the number of the matrix's columns.'.");
        assert(first_col > int(0) &&
               "Take the successive n_cols by index requires 'The indexes can't be smaller than 0'.");
        static_assert(n_rows == rows_, "The number of the Matrices' rows should meet.");
        for (size_t j = 0; j != n_cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t i = 0; i != n_rows; ++i) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat(i, j + first_col);
            }
        }
        return *this;
    }

    /**
     * @brief Take successive columns of a matrix by index and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .cols(Mat, int) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .cols_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam  _cols The number of columns taken.
     * @param first_col The first column index.
     * @return (Mat&) The successive columns(a normal matrix) (a reference to 'this') .
     */
    template <size_t _cols>
    Mat<T, n_rows, _cols, MatType::NORMAL> cols(size_t first_col) {
        assert(first_col < n_cols && "Take the successive n_cols by index requires 'The indexes should be smaller than "
                                     "the number of the matrix's columns.'.");
        assert(first_col > int(0) &&
               "Take the successive n_cols by index requires 'The indexes can't be smaller than 0'.");
        static_assert(_cols > int(1),
                      "Take the successive n_cols by index requires 'The total number can't be smaller than 1'.");

        Mat<T, n_rows, _cols, MatType::NORMAL> mat;
        for (size_t j = 0; j != _cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t i = 0; i != n_rows; ++i) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(i, j) = (*this)(i, j + first_col);
            }
        }
        return mat;
    }
    /**
     * @brief Take seccessive columns of a matrix by indexes as a read only view.
     *
     * @return (MatViewCols<first_row, last_row, T, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only columns
     * vector view.
     */
    template <size_t first_col, size_t last_col>
    MatViewCols<first_col, last_col, T, n_rows, n_cols, MatType::NORMAL, MType<type>> Cols_() const {
        static_assert(first_col > int(0),
                      "Take the successive cols by index requires 'The first index can't be smaller than 0'.");
        static_assert(
            first_col < last_col,
            "Take the successive cols by index requires 'The first index should be smaller than the last index'.");
        static_assert(last_col < n_cols, "Take the successive cols by index requires 'The last index should be smaller "
                                         "than the number of the matrix's cols.'.");
        return *this;
    }

    /**
     * @brief Take successive rows of a matrix by index.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The matrix MatType.
     * @param mat The original matrix.
     * @param first_row The first row index.
     * @return (Mat&) The successive rows(a normal matrix) (a reference to 'this') .
     */

    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    Mat& rows(const M<T2, rows_, cols_, type2, _unused...>& mat, size_t first_row) {
        assert(first_row < n_rows && "Take the successive n_rows by index requires 'The indexes should be smaller than "
                                     "the number of the matrix's rows.'.");
        assert(first_row > int(0) &&
               "Take the successive n_rows by index requires 'The indexes can't be smaller than 0'.");
        static_assert(n_cols == cols_, "The number of the Matrices' cols should meet.");
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat(i + first_row, j);
            }
        }
        return *this;
    }

    /**
     * @brief Take successive rows of a matrix by index and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .rows(Mat, int) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .rows_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam  _rows The number of columns taken.
     * @param first_row The first row index.
     * @return (Mat&) The successive rows(a normal matrix) (a reference to 'this') .
     */

    template <size_t _rows>
    Mat<T, _rows, n_cols, MatType::NORMAL> rows(size_t first_row) {
        assert(first_row < n_rows && "Take the successive n_rows by index requires 'The indexes should be smaller than "
                                     "the number of the matrix's rows.'.");
        assert(first_row > int(0) &&
               "Take the successive n_rows by index requires 'The indexes can't be smaller than 0'.");
        static_assert(_rows > int(1),
                      "Take the successive n_rows by index requires 'The total number can't be smaller than 1'.");

        Mat<T, _rows, n_cols, MatType::NORMAL> mat;
        for (size_t i = 0; i != _rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(i, j) = (*this)(i + first_row, j);
            }
        }
        return mat;
    }

    /**
     * @brief Take successive rows of a matrix by indexes as a read only view.
     *
     * @return (MatViewRows<first_row, last_row, T, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only rows
     * vector view.
     */
    template <size_t first_row, size_t last_row>
    MatViewRows<first_row, last_row, T, n_rows, n_cols, MatType::NORMAL, MType<type>> Rows_() const {
        static_assert(
            first_row < last_row,
            "Take the successive rows by index requires 'The first index should be smaller than the last index'.");
        static_assert(last_row < n_rows, "Take the successive rows by index requires 'The last index should be smaller "
                                         "than the number of the matrix's rows.'.");
        return *this;
    }

    /**
     * @brief Take discrete rows of a matrix by container.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M1 The original matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The vector type.
     * @tparam T1 The original matrix element type.
     * @tparam type1 The matrix MatType.
     * @tparam rows_ The row number of the matrix.
     * @tparam cols_ The column number of the matrix.
     * @param mat The original matrix.
     * @param vector The vector.
     * @return (Mat&) The discrete rows(a normal matrix) (a reference to 'this') .
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1, typename M2,
              typename T1, MatType type1, size_t rows_, size_t cols_>
    Mat& rows(const M1<T1, rows_, cols_, type1, _unused1...>& mat, M2 vector) {
        static_assert(n_cols == cols_, "The number of the Matrices' cols should meet.");
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            assert(vector[i] > int(0) &&
                   "Take the discrete n_rows by index requires 'The indexes can't be smaller than 0'.");
            assert(vector[i] < n_rows && "Take the discrete n_rows by index requires 'The indexes should be smaller "
                                         "than the number of the matrix's rows.'.");
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat(vector[i], j);
            }
        }
        return *this;
    }

    /**
     * @brief Take discrete rows of a matrix by container and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .rows(Mat, vector) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .rows_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M2 The vector type.
     * @param vector The vector.
     * @return (Mat&) The discrete rows(a normal matrix) (a reference to 'this') .
     */
    template <size_t _rows, typename M2>
    Mat<T, _rows, n_cols, MatType::NORMAL> rows(M2 vector) {
        Mat<T, _rows, n_cols, MatType::NORMAL> mat;
        for (size_t i = 0; i != _rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            assert(vector[i] > int(0) &&
                   "Take the discrete n_rows by index requires 'The indexes can't be smaller than 0'.");
            assert(vector[i] < n_rows && "Take the discrete n_rows by index requires 'The indexes should be smaller "
                                         "than the number of the matrix's rows.'.");
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(i, j) = (*this)(vector[i], j);
            }
        }
        return mat;
    }

    /**
     * @brief Take discrete rows of a matrix by container as a read only view.
     *
     * @return (MatViewRowsContainer<Rows, T, M, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only rows
     * vector view.
     */
    template <size_t Rows, typename M>
    MatViewRowsContainer<Rows, T, M, n_rows, n_cols, MatType::NORMAL, MType<type>> Rows_(const M container) const {
        MatViewRowsContainer<Rows, T, M, n_rows, n_cols, MatType::NORMAL, MType<type>> View(*this, container);
        return View;
    }

    /**
     * @brief Take discrete cols of a matrix by container.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M1 The original matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The vector type.
     * @tparam T1 The original matrix element type.
     * @tparam type1 The matrix MatType.
     * @tparam rows_ The row number of the matrix.
     * @tparam cols_ The column number of the matrix.
     * @param mat The original matrix.
     * @param vector The vector.
     * @return (Mat&) The discrete cols(a normal matrix) (a reference to 'this') .
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1, typename M2,
              typename T1, MatType type1, size_t rows_, size_t cols_>
    Mat& cols(const M1<T1, rows_, cols_, type1, _unused1...>& mat, M2 vector) {
        static_assert(n_rows == rows_, "The number of the Matrices' rows should meet.");
        for (size_t j = 0; j != n_cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            assert(vector[j] > int(0) &&
                   "Take the discrete n_cols by index requires 'The indexes can't be smaller than 0'.");
            assert(vector[j] < n_cols && "Take the discrete n_cols by index requires 'The indexes should be smaller "
                                         "than the number of the matrix's cols.'.");
            for (size_t i = 0; i != n_rows; ++i) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                (*this)(i, j) = mat(i, vector[j]);
            }
        }
        return *this;
    }

    /**
     * @brief Take discrete cols of a matrix by container and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .cols(Mat, vector) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .cols_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M2 The vector type.
     * @param vector The vector.
     * @return (Mat&) The discrete cols(a normal matrix) (a reference to 'this') .
     */
    template <size_t _cols, typename M2>
    Mat<T, n_rows, _cols, MatType::NORMAL> cols(M2 vector) {
        Mat<T, n_rows, _cols, MatType::NORMAL> mat;
        for (size_t j = 0; j != _cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            assert(vector[j] > int(0) &&
                   "Take the discrete n_cols by index requires 'The indexes can't be smaller than 0'.");
            assert(vector[j] < n_cols && "Take the discrete n_cols by index requires 'The indexes should be smaller "
                                         "than the number of the matrix's cols.'.");
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(i, j) = (*this)(i, vector[j]);
            }
        }
        return mat;
    }

    /**
     * @brief Take discrete cols of a matrix by container as a read only view.
     *
     * @return (MatViewColsContainer<Cols, T, M, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only columns
     * vector view.
     */
    template <size_t Cols, typename M>
    MatViewColsContainer<Cols, T, M, n_rows, n_cols, MatType::NORMAL, MType<type>> Cols_(const M container) const {
        MatViewColsContainer<Cols, T, M, n_rows, n_cols, MatType::NORMAL, MType<type>> View(*this, container);
        return View;
    }

    /**
     * @brief Apply 'abs' (absolute value) to itself.
     *
     * @note This does not support complex numbers yet.
     * @return (Mat&) The processed matrix (a reference to 'this').
     */
    Mat& abs_() {
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_ABS_UNROLL_FACTOR)
            if (_data[i] < 0) _data[i] = -_data[i];
        }
        return *this;
    }

    /**
     * @brief Transpose.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @tparam M The matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The original matrix MatType.
     * @param mat The original matrix.
     * @return (Mat&) The transpose matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& t(const M<T2, n_cols, n_rows, type2, _unused...>& mat) {
        if (type == MatType::DIAGONAL || type == MatType::SCALAR || type == MatType::SYM) {
            // do nothing
        } else if (type == MatType::NORMAL) {
        MAT_TRANSPOSE_NORMAL:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = 0; j != n_rows; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    _data[n_rows * i + j] = mat._data[n_rows * j + i];
                }
            }
        } else if (type == MatType::UPPER) {
        MAT_TRANSPOSE_UPPER:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = i; j != n_rows; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    (*this)(i, j) = mat(j, i);
                }
            }
        } else if (type == MatType::LOWER) {
        MAT_TRANSPOSE_LOWER:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = 0; j <= i; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    (*this)(i, j) = mat(j, i);
                }
            }
        } else if (type == MatType::SUPPER) {
        MAT_TRANSPOSE_SUPPER:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = i; j != n_rows; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    (*this)(i, j) = mat(j, i);
                }
            }
        } else if (type == MatType::SLOWER) {
        MAT_TRANSPOSE_SLOWER:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = i; j != n_rows; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    (*this)(i, j) = mat(j, i);
                }
            }
        } else if (type == MatType::ASYM) {
        MAT_TRANSPOSE_ASYM:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = i + 1; j != n_rows; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    (*this)(i, j) = mat(j, i);
                }
            }
        } else {
            assert(!"Impossible");
        }
        return *this;
    }

    /**
     * @brief Transpose NORMAL matrix as a copy.
     *
     * @details You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat<T, n_cols, n_rows, type>) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type,
              typename std::enable_if_t<_type == MatType::NORMAL, bool> = true>
    Mat<T, n_cols, n_rows, type> t() const {
        // Ref: https://stackoverflow.com/a/13401982/15080514
        static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::t()!");
        static_assert(n_cols > 0 && n_rows > 0, "The matrix should have size when transposing.");
        Mat<T, n_cols, n_rows, type> mat;
    MAT_TRANSPOSE_NORMAL:
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
            for (size_t j = 0; j != n_rows; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(i, j) = (*this)(j, i);
            }
        }
        return mat;
    }

    /**
     * @brief Transpose DIAGONAL/SCALAR/SYM matrix as a copy.
     *
     * @details This is actually *this (i.e. the trapsonse is the same as itself).
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type,
              typename std::enable_if_t<
                  (_type == MatType::DIAGONAL || _type == MatType::SCALAR || _type == MatType::SYM), bool> = true>
    Mat t() const {
        static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::t()!");
        return *this;
    }

    /**
     * @brief Transpose UPPER matrix as a copy.
     *
     * @details You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat<T, n_rows, n_cols, MatType::LOWER>) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type,
              typename std::enable_if_t<_type == MatType::UPPER, bool> = true>
    Mat<T, n_rows, n_cols, MatType::LOWER> t() const {
        static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::t()!");
        static_assert(n_cols > 0 && n_rows > 0, "The matrix should have size when transposing.");
        Mat<T, n_rows, n_cols, MatType::LOWER> mat;
    MAT_TRANSPOSE_NORMAL:
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
            for (size_t j = 0; j != n_rows; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(i, j) = (*this)(j, i);
            }
        }
        return mat;
    }

    /**
     * @brief Transpose LOWER matrix as a copy.
     *
     * @details You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat<T, n_rows, n_cols, MatType::UPPER>) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type,
              typename std::enable_if_t<_type == MatType::LOWER, bool> = true>
    Mat<T, n_rows, n_cols, MatType::UPPER> t() const {
        static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::t()!");
        static_assert(n_cols > 0 && n_rows > 0, "The matrix should have size when transposing.");
        Mat<T, n_rows, n_cols, MatType::UPPER> mat;
    MAT_TRANSPOSE_LOWER:
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
            for (size_t j = i; j <= i; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(j, i) = (*this)(i, j);
            }
        }
        return mat;
    }

    /**
     * @brief Transpose SUPPER matrix as a copy.
     *
     * @details You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat<T, n_rows, n_cols, MatType::SLOWER>) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type,
              typename std::enable_if_t<_type == MatType::SUPPER, bool> = true>
    Mat<T, n_rows, n_cols, MatType::SLOWER> t() const {
        static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::t()!");
        static_assert(n_cols > 0 && n_rows > 0, "The matrix should have size when transposing.");
        Mat<T, n_rows, n_cols, MatType::SLOWER> mat;
    MAT_TRANSPOSE_SUPPER:
        for (size_t i = 0; i != n_cols - 1; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_rows; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(j, i) = (*this)(i, j);
            }
        }
        return mat;
    }

    /**
     * @brief Transpose SLOWER matrix as a copy.
     *
     * @details You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat<T, n_rows, n_cols, MatType::SUPPER>) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type,
              typename std::enable_if_t<_type == MatType::SLOWER, bool> = true>
    Mat<T, n_rows, n_cols, MatType::SLOWER> t() const {
        static_assert(sizeof...(_unused) == 0, "Do not specify template arguments for Mat::t()!");
        static_assert(n_cols > 0 && n_rows > 0, "The matrix should have size when transposing.");
        Mat<T, n_rows, n_cols, MatType::SUPPER> mat;
    MAT_TRANSPOSE_SLOWER:
        for (size_t i = 1; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_rows - 1; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(j, i) = (*this)(i, j);
            }
        }
        return mat;
    }

    /**
     * @brief Transpose ASYM matrix as a copy.
     *
     * @details You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .t(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .t_() which creates a read only view,
     *       and this operation does not copy.
     * @return (Mat) The transpose result as a copy.
     */
    template <typename... _unused, MatType _type = type, typename std::enable_if_t<_type == MatType::ASYM, bool> = true>
    Mat t() const {
        static_assert(n_cols > 0 && n_rows > 0, "The matrix should have size when transposing.");
        Mat mat;
    MAT_TRANSPOSE_ASYM:
        for (size_t i = 1; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
            for (size_t j = i + 1; j != n_rows - 1; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                mat(j, i) = (*this)(i, j);
            }
        }
        return mat;
    }

    /**
     * @brief Transpose as a read only view.
     *
     * @return (MatViewT<T, n_cols, n_rows, tType(type)>) The read only view transpose.
     */
    MatViewT<T, n_cols, n_rows, tType(type)> t_() { return const_cast<T*>(_data); }

    MatViewT<T, n_cols, n_rows, tType(type)> t_() const { return const_cast<T*>(_data); }

    /**
     * @brief In-place transpose.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do transpose in parallel.
     * @return (Mat&) The transposed result (a reference to 'this').
     */
    Mat& tSelf() {
        static_assert(type != MatType::UPPER, "Upper matrix cannot perform in place transpose.");
        static_assert(type != MatType::LOWER, "Lower matrix cannot perform in place transpose.");
        static_assert(type != MatType::SUPPER, "Strict upper matrix cannot perform in place transpose.");
        static_assert(type != MatType::SLOWER, "Strict lower matrix cannot perform in place transpose.");
        static_assert(n_rows == n_cols, "In-place transpose requires a square matrix.");
        if (type == MatType::NORMAL) {
        MAT_NORMAL_TRANSPOSE_SELF:
            for (size_t i = 0; i != n_cols; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                for (size_t j = i; j != n_cols; ++j) {
                    FLAMES_PRAGMA(LOOP_FLATTEN)
                    std::swap(_data[n_cols * i + j], _data[n_cols * j + i]);
                }
            }
        } else if (type == MatType::ASYM) {
        MAT_AYSM_TRANSPOSE_SELF:
            for (size_t i = 0; i != n_cols * (n_cols - 1) / 2; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_TRANSPOSE_UNROLL_FACTOR)
                _data[i] = -_data[i];
            }
        }
        return *this;
    }

    /**
     * @brief Calculate the opposite of a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @return (Mat&) The opposite matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2>
    Mat& opp(const M<T2, n_rows, n_cols, type, _unused...>& mat) {
    MAT_OPP:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = -mat[i];
        }
        return *this;
    }

    /**
     * @brief Calculate the opposite of a matrix and make a copy.
     *
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .opp(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .opp_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @return (Mat) The opposite matrix copy.
     */
    Mat opp() const {
        Mat mat;
    MAT_OPP:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = -_data[i];
        }
        return mat;
    }

    /**
     * @brief Matrix opposite as a read only view.
     *
     * @return (MatViewOpp<T, n_rows, n_cols, type>) The read only view opposite.
     */
    MatViewOpp<T, n_rows, n_cols, type> opp_() const {
        FLAMES_PRAGMA(INLINE)
        return *this;
    }

    /**
     * @brief In-place matrix opposite.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @return (Mat&) The opposite matrix (a reference to 'this').
     */
    Mat& oppSelf() {
    MAT_OPP_SELF:
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = -_data[i];
        }
    }

    /**
     * @brief Take the diagonal of a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @return (Mat&) The diagonal matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& diagMat(const M<T2, n_rows, n_cols, type2, _unused...>& mat) {
        static_assert(n_rows == n_cols, "Take the diagonal requires 'n_rows == n_cols'.");
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            (*this)(i, i) = mat(i, i);
        }
        return *this;
    }

    /**
     * @brief Take the diagonal of a matrix and make a copy.
     *
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .diagMat(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .diagMat_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @return (Mat<T, n_rows, n_cols, MatType::DIAGONAL>) The diagonal matrix copy.
     */
    Mat<T, n_rows, n_cols, MatType::DIAGONAL> diagMat() const {
        Mat<T, n_rows, n_cols, MatType::DIAGONAL> mat;
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat(i, i) = (*this)(i, i);
        }
        return mat;
    }

    /**
     * @brief Take the diagonal of a matrix as a read only view.
     *
     * @return (MatViewDiagMat<T, n_rows, n_cols, MatType::DIAGONAL, MType<type>>) The read only diagonal matrix view.
     */
    MatViewDiagMat<T, n_rows, n_cols, MatType::DIAGONAL, MType<type>> diagMat_() const {
        static_assert(n_rows == n_cols, "Take the diagonal requires 'n_rows == n_cols'.");
        return *this;
    }

    /**
     * @brief Take the diagonal vector of a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @return (Mat&) The diagonal vector (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& diagVec(const M<T2, n_rows, n_rows, type2, _unused...>& mat) {
        static_assert(n_cols == 1, "Diagonal vector has column number as 1.");
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = mat(i, i);
        }
        return *this;
    }

    /**
     * @brief Take the diagonal vector of a matrix and make a copy.
     *
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .diagVec(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .diagVec_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @return (Vec<T, n_rows>) The diagonal vector copy.
     */
    Vec<T, n_rows> diagVec() const {
        static_assert(n_rows == n_cols, "Take the diagonal requires 'n_rows == n_cols'.");
        Vec<T, n_rows> mat;
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat._data[i] = (*this)(i, i);
        }
        return mat;
    }

    /**
     * @brief Take the diagonal vector of a matrix as a read only view.
     *
     * @return (MatViewDiagVec<T, n_rows, n_cols, MatType::DIAGONAL, MType<type>>) The read only diagonal vector view.
     */
    MatViewDiagVec<T, n_rows, 1, MatType::NORMAL, MType<type>> diagVec_() const {
        static_assert(n_rows == n_cols, "Take the diagonal requires 'n_rows == n_cols'.");
        return *this;
    }

    /**
     * @brief Take the diagonal row vector of a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @return (Mat&) The diagonal vector (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& diagRowVec(const M<T2, n_cols, n_cols, type2, _unused...>& mat) {
        static_assert(n_rows == 1, "Diagonal vector has column number as 1.");
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = mat(i, i);
        }
        return *this;
    }

    /**
     * @brief Take the diagonal row vector of a matrix and make a copy.
     *
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .diagRowVec(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .diagRowVec_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @return (RowVec<T, n_cols>) The diagonal row vector copy.
     */
    RowVec<T, n_cols> diagRowVec() const {
        static_assert(n_rows == n_cols, "Take the diagonal requires 'n_rows == n_cols'.");
        RowVec<T, n_cols> mat;
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat._data[i] = (*this)(i, i);
        }
        return mat;
    }

    /**
     * @brief Take the diagonal row vector of a matrix as a read only view.
     *
     * @return (MatViewDiagRowVec<T, 1, n_cols, MatType::NORMAL, MType<type>> ) The read only diagonal row vector view.
     */
    MatViewDiagRowVec<T, 1, n_cols, MatType::NORMAL, MType<type>> diagRowVec_() const {
        static_assert(n_rows == n_cols, "Take the diagonal requires 'n_rows == n_cols'.");
        return *this;
    }

    /**
     * @brief Take the off diagonal of a matrix.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @return (Mat&) The off diagonal matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& offDiag(const M<T2, n_cols, n_rows, type2, _unused...>& mat) {
        static_assert(n_rows == n_cols, "Take the off diagonal requires 'n_rows == n_cols'.");
        static_assert(type == MatType::NORMAL, "Take the off diagonal requires the matrix to be NORMAL.");
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                if (i != j) (*this)(i, j) = mat(i, j);
                else (*this)(i, j) = T(0);
            }
        }
        return *this;
    }

    /**
     * @brief Take the off diagonal of a matrix and make a copy.
     *
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .offDiag(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .offDiag_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @return (Mat<T, n_rows, n_cols, MatType::NORMAL>) The off diagonal matrix copy.
     */
    Mat<T, n_rows, n_cols, MatType::NORMAL> offDiag() const {
        static_assert(n_rows == n_cols, "Take the off diagonal requires 'n_rows == n_cols'.");
        static_assert(type == MatType::NORMAL, "Take the off diagonal requires the matrix to be NORMAL.");
        Mat<T, n_rows, n_cols, MatType::NORMAL> mat;
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            for (size_t j = 0; j != n_cols; ++j) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                if (i != j) mat(i, j) = (*this)(i, j);
                else mat(i, j) = T(0);
            }
        }
        return mat;
    }

    /**
     * @brief Take the diagonal of a matrix as a read only view.
     *
     * @return (MatViewOffDiag<T, n_rows, n_cols, MatType::NORMAL, MType<type>>) The read only off diagonal matrix view.
     */
    MatViewOffDiag<T, n_rows, n_cols, MatType::NORMAL, MType<type>> offDiag_() const {
        static_assert(n_rows == n_cols, "Take the off diagonal requires 'n_rows == n_cols'.");
        return *this;
    }

    /**
     * @brief Inverse the diagonal matrix.
     *
     * @details This can also be used for calculating the inverse of the diagonal part.
     * @note This should only be applied to matrix ('this') with MatType as DIAGONAL.
     *       If this rule is violated, a static assert will raise an error.
     *       This, however, allows the paramater matrix to be non DIAGONAL.
     *       But you need to make sure it is a diagonal matrix in mathematics
     *       otherwise the inverse only applies to its diagonal.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The original matrix MatType.
     * @param mat The original matrix.
     * @return (Mat&) The inverse matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& invDiag(const M<T2, n_cols, n_rows, type2, _unused...>& mat) {
        static_assert(type == MatType::DIAGONAL, "'invDiag' is only used for diagonal matrix.");
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = T(1.0) / mat(i, i);
        }
        return *this;
    }

    /**
     * @brief Inverse the diagonal matrix and makes a copy.
     *
     * @note This matrix should be DIAGONAL.
     * @return (Mat) The inverse matrix copy.
     */
    Mat invDiag() const {
        static_assert(type == MatType::DIAGONAL, "'invDiag' is only used for diagonal matrix.");
        Mat mat;
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_INV_UNROLL_FACTOR)
            mat._data[i] = T(1.0) / _data[i];
        }
    }

    /**
     * @brief Matrix inverse using Newton-Schulz iterative method (NSA).
     *
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The original matrix MatType.
     * @param mat The original matrix.
     * @param iter The number of iterations (default as 4).
     * @return (Mat&) The inverse matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    Mat& invNSA(const M<T2, n_rows, n_cols, type2, _unused...>& mat, size_t iter = 4) {
        static_assert(n_rows == n_cols, "Calculate inverse needs to be a square matrix.");
        assert(iter >= 1 && "At least one iteration is needed.");
        const auto D = mat.diagMat_(); // diagonal part
        const auto E = mat.offDiag_(); // off-diagonal part
        Mat<T, n_rows, n_cols, MatType::DIAGONAL> D_inv;
        D_inv.invDiag(D); // inverse of diagonal part
        // D_inv.print("D_inv:\n");
        // E.asMat().print("E:\n");
        auto D_inv_opp                                  = -D_inv;
        Mat<T, n_rows, n_cols, MatType::NORMAL> product = D_inv_opp * E;
        Mat<T, n_rows, n_cols, MatType::NORMAL> sum_tmp = *this = product; // the first iteration
        Mat<T, n_rows, n_cols, MatType::NORMAL> tmp;
    MAT_INV_NSA:
        for (size_t i = 1; i < iter; ++i) {
            tmp.mul(*this, product);
            *this = tmp;
            sum_tmp += tmp;
        }
        this->mul(sum_tmp, D_inv);
        return *this += D_inv;
    }

    /**
     * @brief Matrix inverse using Newton-Schulz iterative method (NSA) as a copy.
     *
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .invNSA(Mat) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     * @param iter The number of iterations (default as 4).
     * @return (Mat) The inverse matrix copy.
     */
    Mat invNSA(size_t iter = 4) const {
        static_assert(n_rows == n_cols, "Calculate inverse needs to be a square matrix.");
        Mat mat;
        mat.invNSA(*this, iter);
        return mat;
    }

    /**
     * @brief Matrix inverse using improved Newton-Schulz iterative method (INSA).
     *
     * @details This implements a coefficient to the original NSA method to accelerate convergence.
     * @note This function has not been implemented!
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @tparam type2 The original matrix MatType.
     * @tparam coeff_type The coefficient data type.
     * @param mat The original matrix.
     * @param iter The number of iterations (default as 3).
     * @param beta The coefficient applied to each iteration (default as 2).
     * @return (Mat&) The inverse matrix (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, typename coeff_type>
    Mat& invINSA(const M<T2, n_rows, n_cols, type2, _unused...>& mat, size_t iter = 3, coeff_type beta = 2) {
        static_assert(n_rows == n_cols, "Calculate inverse needs to be a square matrix.");
        assert(iter >= 1 && "At least one iteration is needed.");
        return *this;
    }

    /**
     * @brief Matrix inverse using improved Newton-Schulz iterative method (INSA) as a copy.
     *
     * @tparam coeff_type The coefficient data type.
     * @param iter The number of iterations (default as 3).
     * @param beta The coefficient applied to each iteration (default as 2).
     * @return (Mat) The inverse matrix copy.
     */
    template <typename coeff_type>
    Mat invINSA(size_t iter = 3, coeff_type beta = 1.0) const {
        static_assert(n_rows == n_cols, "Calculate inverse needs to be a square matrix.");
        Mat mat;
        mat.invINSA(*this, iter);
        return mat;
    }

    // Does not support complex number now.
    template <typename Tp = T>
    Tp power() const {
        Tp p = 0;
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_POWER_UNROLL_FACTOR)
            p += static_cast<Tp>(_data[i] * _data[i]);
        }
        return p;
    }

    /**
     * @brief Get the value from 1x1 matrix.
     *
     * @return (T) The element value.
     */
    T value() const {
        static_assert(n_rows == 1 && n_cols == 1, "This only applies to matrix of size 1x1.");
        return _data[0];
    }

    /**
     * @brief Get the writeable value from 1x1 matrix.
     *
     * @return (T&) The element reference.
     */
    T& value() {
        static_assert(n_rows == 1 && n_cols == 1, "This only applies to matrix of size 1x1.");
        return _data[0];
    }

    /**
     * @brief Print the matrix.
     *
     * @param str The string to print before the matrix.
     * @param os The out stream (default as std::cout).
     */
    void print(const std::string& str = "", std::ostream& os = std::cout) const {
#ifndef __SYNTHESIS__
        static_assert(n_rows > 0 && n_cols > 0, "Print requires this matrix be of valid size.");
        os << str << "[";
        for (size_t i = 0; i + 1 < n_rows; ++i) {
            os << "[";
            for (size_t j = 0; j + 1 < n_cols; ++j) os << static_cast<T>((*this)(i, j)) << ", ";
            os << static_cast<T>((*this)(i, n_cols - 1)) << "]," << std::endl;
        }
        os << "[";
        for (size_t j = 0; j + 1 < n_cols; ++j) os << static_cast<T>((*this)(n_rows - 1, j)) << ", ";
        os << static_cast<T>((*this)(n_rows - 1, n_cols - 1)) << "]]" << std::endl;
#endif
    }

    /**
     * @brief The unary plus operator.
     *
     * @return (MatView<T, n_rows, n_cols, type>) Return the MatView of the matrix.
     */
    MatView<T, n_rows, n_cols, type> operator+() const {
        FLAMES_PRAGMA(INLINE)
        return *this;
    }

    /**
     * @brief The unary minus operator.
     *
     * @details This is the same as .opp_().
     * @return (MatViewOpp<T, n_rows, n_cols, type>) The read only opposite matrix view.
     */
    MatViewOpp<T, n_rows, n_cols, type> operator-() const {
        FLAMES_PRAGMA(INLINE)
        return *this;
    }

    //   public: // original private
  public:
    /**
     * @brief Get the raw data array pointer.
     *
     * @details This operation is dangerous if you are not aware of what you are doing
     *          so this function is marked as private.
     * @return (const T*) Raw data pointer.
     */
    T* rawDataPtr() { return _data; }

    const T* rawDataPtr() const { return _data; }

    /**
     * @brief Try to assign a value to a specific position.
     *
     * @details This is useful when we do not want to check the MatType.
     * @note This will not check the boundary!
     *       If that is not sure, use ._tryAssignOutRange() instead.
     * @param r The row index.
     * @param c The column index.
     * @param value The value to be assigned.
     */
    inline void _tryAssign(size_t r, size_t c, T value) {
        FLAMES_PRAGMA(INLINE)
        if (type == MatType::DIAGONAL) {
            if (r == c) (*this)(r, c) = value;
        } else if (type == MatType::UPPER) {
            if (r <= c) (*this)(r, c) = value;
        } else if (type == MatType::LOWER) {
            if (r >= c) (*this)(r, c) = value;
        } else if (type == MatType::SUPPER) {
            if (r < c) (*this)(r, c) = value;
        } else if (type == MatType::SLOWER) {
            if (r > c) (*this)(r, c) = value;
        } else if (type == MatType::SYM) {
            if (r <= c) (*this)(r, c) = value;
        } else if (type == MatType::ASYM) {
            if (r < c) (*this)(r, c) = value;
        } else if (type == MatType::NORMAL) {
            (*this)(r, c) = value;
        } else {
            assert(!"Impossible! Unknown type!");
        }
    }

    /**
     * @brief Try to assign a value to a specific position that may have out-of-range index.
     *
     * @details This is useful when we do not want to check the MatType and the index range.
     * @param r The row index.
     * @param c The column index.
     * @param value The value to be assigned.
     */
    inline void _tryAssignOutRange(size_t r, size_t c, T value) {
        FLAMES_PRAGMA(INLINE)
        if (r >= n_rows || c >= n_cols) return;
        else _tryAssign(r, c, value);
    }

    /**
     * @brief Try to plus a value to a specific position.
     *
     * @details This is useful when we do not want to check the MatType.
     * @note This will not check the boundary!
     * @param r The row index.
     * @param c The column index.
     * @param value The value to add.
     */
    inline void _tryPlus(size_t r, size_t c, T value) {
        FLAMES_PRAGMA(INLINE)
        if (type == MatType::DIAGONAL) {
            if (r == c) (*this)(r, c) += value;
        } else if (type == MatType::UPPER) {
            if (r <= c) (*this)(r, c) += value;
        } else if (type == MatType::LOWER) {
            if (r >= c) (*this)(r, c) += value;
        } else if (type == MatType::SUPPER) {
            if (r < c) (*this)(r, c) += value;
        } else if (type == MatType::SLOWER) {
            if (r > c) (*this)(r, c) += value;
        } else if (type == MatType::ASYM) {
            if (r != c) (*this)(r, c) += value;
        } else if (type == MatType::SYM) {
            (*this)(r, c) += value;
        } else if (type == MatType::NORMAL) {
            (*this)(r, c) += value;
        } else {
            assert(!"Impossible! Unknown type!");
        }
    }

    /**
     * @brief Systolic array read the first column from the left matrix.
     *
     * @tparam Tmp The temp mat type.
     * @tparam M The left matrix type.
     * @tparam Comm The
     * @tparam Zero
     * @param i The iteration number.
     * @param tmp_L The temp matrix.
     * @param mat_L The left matrix.
     * @param begin_shift The begin shift value.
     * @param zero The zero value.
     */
    template <typename Tmp, typename M, typename Comm, typename Zero>
    inline void _sa_read_first_col_L(size_t i, Tmp& tmp_L, const M& mat_L, size_t begin_shift,
                                     Comm __attribute__((unused)) _comm_foo, Zero zero) const {
        FLAMES_PRAGMA(INLINE)
        constexpr size_t comm = Comm::value;
    read_first_col_L:
        for (size_t j = 0; j < n_rows; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            tmp_L[j][0] = ((j <= i + begin_shift) && (i + begin_shift - j < comm))
                              ? mat_L(j, comm - 1 - (i + begin_shift - j))
                              : zero;
        }
    }

    /**
     * @brief Systolic array read the first column from the right matrix.
     *
     * @tparam Tmp The temp mat type.
     * @tparam M The right matrix type.
     * @tparam Comm The common size (column number of the left matrix and the row number of the right matrix) as a type.
     * @tparam Zero The zero element type.
     * @param i The iteration number.
     * @param tmp_L The temp matrix.
     * @param mat_L The right matrix.
     * @param begin_shift The begin shift value.
     * @param zero The zero value.
     */
    template <typename Tmp, typename M, typename Comm, typename Zero>
    inline void _sa_read_first_row_R(size_t i, Tmp& tmp_R, const M& mat_R, size_t begin_shift,
                                     Comm __attribute__((unused)) _comm_foo, Zero zero) const {
        FLAMES_PRAGMA(INLINE)
        constexpr size_t comm = Comm::value;
    read_first_row_R:
        for (size_t j = 0; j < n_cols; ++j) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            tmp_R[0][j] = ((j <= i + begin_shift) && (i + begin_shift - j < comm))
                              ? mat_R(comm - 1 - (i + begin_shift - j), j)
                              : zero;
        }
    }

    /**
     * @brief Systolic array multiplication.
     *
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam HLSVec The HLS vector type for control whether assigning or adding.
     * @param tmp_L The left tmp matrix.
     * @param tmp_R The right tmp matrix.
     * @param use_assign The HLS vector for controlling whether assigning or adding.
     */
    template <typename T1, typename T2, typename HLSVec>
    inline void _sa_multiply(const T1& tmp_L, const T2& tmp_R, const HLSVec& use_assign) {
        FLAMES_PRAGMA(INLINE)
    multiply:
        for (size_t r = 0; r != n_rows; ++r) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
            for (size_t c = 0; c != n_cols; ++c) {
                FLAMES_PRAGMA(LOOP_FLATTEN)
                T result = tmp_L[r][c] * tmp_R[r][c];
                if (use_assign[r * n_cols + c]) {
                    _tryAssign(r, c, result);
                    // std::cout << "assign " << r << ' ' << c << " (" << tmp_L[r][c] << "x" << tmp_R[r][c] << "=" <<
                    // result << ")\n";
                } else {
                    _tryPlus(r, c, result);
                    // std::cout << "plus " << r << ' ' << c << " (" << tmp_L[r][c] << "x" << tmp_R[r][c] << "=" <<
                    // result << ")\n";
                }
            }
        }
    };

    /**
     * @brief The systolic array multiplication.
     *
     * @tparam M1 The left matrix type.
     * @tparam _unused1 (unused)
     * @tparam M2 The right matrix type.
     * @tparam _unused2 (unused)
     * @tparam T1 The left matrix element type.
     * @tparam T2 The right matrix element type.
     * @tparam type1 The left matrix MatType.
     * @tparam type2 The right matrix MatType.
     * @tparam comm The common size (column number of the left matrix and the row number of the right matrix) as a type.
     * @param mat_L The left matrix.
     * @param mat_R The right matrix.
     * @param begin_shift The begin shift value.
     * @param end_shift The end shift value.
     * @return (Mat&) The multiplication result (a reference to 'this').
     */
    template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
              template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
              typename T2, MatType type1, MatType type2, size_t comm>
    Mat& _systolicArrayMul(const M1<T1, n_rows, comm, type1, _unused1...>& mat_L,
                           const M2<T2, comm, n_cols, type2, _unused2...>& mat_R, size_t begin_shift = 0,
                           size_t end_shift = 0) {
        // constexpr size_t begin_shift = BeginShift::value;
        // constexpr size_t end_shift = EndShift::value;
        T1 tmp_L[n_rows][n_cols];
        T2 tmp_R[n_rows][n_cols];
        std::vector<bool> use_assign(n_rows * n_cols, true);
// hls::vector<bool, n_rows * n_cols> use_assign = true;
#ifdef FLAMES_MAT_PARTITION_COMPLETE
        FLAMES_PRAGMA(ARRAY_PARTITION variable = tmp_L type = complete)
        FLAMES_PRAGMA(ARRAY_PARTITION variable = tmp_R type = complete)
#else
        FLAMES_PRAGMA(ARRAY_PARTITION variable = tmp_L type = block factor = FLAMES_MAT_PARTITION_FACTOR)
        FLAMES_PRAGMA(ARRAY_PARTITION variable = tmp_R type = block factor = FLAMES_MAT_PARTITION_FACTOR)
#endif
        auto write_L = [&](size_t i) {
            if (i + 1 != n_rows + comm + n_cols - 2 - begin_shift - end_shift) {
            write_L:
                for (size_t c = n_cols; c > 1; --c) {
                    FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
                    for (size_t r = 0; r < n_rows; ++r) {
                        FLAMES_PRAGMA(LOOP_FLATTEN)
                        if (c <= i + 2 && r + c <= i + begin_shift + 2 && r <= i) { tmp_L[r][c - 1] = tmp_L[r][c - 2]; }
                    }
                }
            }
        };
        // systolic array passes data down (data from the right matrix)
        auto write_R = [&](size_t i) {
            if (i + 1 != n_rows + comm + n_cols - 2 - begin_shift - end_shift) {
            write_R:
                for (size_t r = n_rows; r > 1; --r) {
                    FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
                    for (size_t c = 0; c < n_cols; ++c) {
                        FLAMES_PRAGMA(LOOP_FLATTEN)
                        if (r <= i + 2 && c + r <= i + begin_shift + 2 && c <= i) { tmp_R[r - 1][c] = tmp_R[r - 2][c]; }
                    }
                }
            }
        };
        // systolic array passes data right (data from the left matrix)
        auto set_assign_ctl = [&](size_t i) {
            FLAMES_PRAGMA(LATENCY max = 0)
        set_assign_ctl:
            for (size_t r = 0; r < n_rows; ++r) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR)
                size_t c = i - r; // note this is unsigned so there can be overflow
                if (c < n_cols) { use_assign[r * n_cols + c] = false; }
            }
        };
        std::integral_constant<int, comm> _comm_foo = {};
    SYSTOLIC_ARRAY_MAT_MUL:
        for (size_t i = 0; i < n_rows + comm + n_cols - 2 - begin_shift - end_shift; ++i) {
            FLAMES_PRAGMA(PIPELINE)
            // std::cout << "ITER " << i << "\n";
            // read
            _sa_read_first_col_L(i, tmp_L, mat_L, begin_shift, _comm_foo, T1(0));
            _sa_read_first_row_R(i, tmp_R, mat_R, begin_shift, _comm_foo, T2(0));
            // multiply
            _sa_multiply(tmp_L, tmp_R, use_assign);
            set_assign_ctl(i);
            // write
            write_L(i);
            write_R(i);
        }
        return *this;
    }

  public: // original private
    /**
     * @brief The raw data array in row major.
     *
     * @details You can set the array partition using macro
     *          `FLAMES_MAT_PARTITION_COMPLETE` to set a complete array partition
     *          and `FLAMES_MAT_PARTITION_FACTOR` to set a block partition with the specific factor.
     *
     */
    T _data[type == MatType::NORMAL     ? n_rows * n_cols
            : type == MatType::DIAGONAL ? n_rows
            : type == MatType::SCALAR   ? 1
            : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
            : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
            : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                        : (1 + n_rows) * n_rows / 2];
};

template <typename T, size_t n_rows, size_t n_cols, MatType type>
class MatView {
  public:
    /**
     * @brief Construct a new MatView object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatView(const Mat<T, n_rows, n_cols, type>& m) : _data(m.rawDataPtr()) {}
    MatView(Mat<T, n_rows, n_cols, type>& m) : _data(m.rawDataPtr()) {}

    MatView(T* const ptr) : _data(ptr) {}

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatView object.
     */
    MatView(const MatView& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    MatView& operator=(const Mat<T, n_rows, n_cols, type>& m) {
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = m[i];
        }
        return *this;
    }

    template <typename M>
    void assign(M m) {
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = m[i];
        }
    }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        FLAMES_PRAGMA(INLINE)
        assert(r < n_rows && "Matrix row index should be within range");
        assert(c < n_cols && "Matrix col index should be within range");
        if (type == MatType::NORMAL) {
            return _data[r * n_cols + c];
        } else if (type == MatType::DIAGONAL) {
            if (r == c) return _data[r];
            else return T(0);
        } else if (type == MatType::SCALAR) {
            if (r == c) return _data[0];
            else return T(0);
        } else if (type == MatType::UPPER) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return T(0);
        } else if (type == MatType::LOWER) {
            if (r >= c) return _data[(1 + r) * r / 2 + c];
            else return T(0);
        } else if (type == MatType::SUPPER) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 2 * r - 1];
            else return T(0);
        } else if (type == MatType::SLOWER) {
            if (r >= c) return _data[(1 + r) * r / 2 + c - r];
            else return T(0);
        } else if (type == MatType::SYM) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + r - c];
        } else if (type == MatType::ASYM) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r * 2 - 1];
            else if (r > c) return -_data[(2 * n_cols + 1 - c) * c / 2 + r - c * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
            return T(0);
        }
    }

    /**
     * @brief Get writeable data element by row index and column index.
     *
     * @param r The row index (starting from 0).
     * @param c The column index (staring from 0).
     * @return (T) The writeable data element.
     */
    T& operator()(size_t r, size_t c) {
        FLAMES_PRAGMA(INLINE)
        assert(r < n_rows && "Matrix row index should be within range");
        assert(c < n_cols && "Matrix col index should be within range");
        if (type == MatType::NORMAL) {
            return _data[r * n_cols + c];
        } else if (type == MatType::DIAGONAL) {
            if (r == c) return _data[r];
            else assert(!"This element cannot be modified (DIAGONAL).");
        } else if (type == MatType::SCALAR) {
            assert(!"This element cannot be modified (SCALAR).");
        } else if (type == MatType::UPPER) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else assert(!"This element cannot be modified (UPPER).");
        } else if (type == MatType::LOWER) {
            if (r >= c) return _data[(1 + r) * r / 2 + c];
            else assert(!"This element cannot be modified (LOWER).");
        } else if (type == MatType::SUPPER) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 1 - 2 * r];
            else assert(!"This element cannot be modified (SUPPER).");
        } else if (type == MatType::SLOWER) {
            if (r > c) return _data[(1 + r) * r / 2 + c - r];
            else assert(!"This element cannot be modified (SLOWER).");
        } else if (type == MatType::SYM) {
            if (r <= c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + r - c];
        } else if (type == MatType::ASYM) {
            if (r < c) return _data[(2 * n_cols + 1 - r) * r / 2 + c - 1 - 2 * r];
            else if (r > c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - 1 - 2 * c];

            // ATTENTION  This part needs to be perfected , missing a minus sign.
            // Because a minus sign will result in a error about reference.

            else assert(!"This element cannot be modified (ASYM).");
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
        // Just to make the compiler happy.
        return _data[0];
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const {
        FLAMES_PRAGMA(INLINE)
        assert(index < size() && "[index] should be in range in MatView");
        return _data[index];
    }

    /**
     * @brief Get the writeable element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T& operator[](size_t index) {
        FLAMES_PRAGMA(INLINE)
        assert(index < size() && "[index] should be in range in MatView");
        return _data[index];
    }

    /**
     * @brief Set all elements of the matrix to a value.
     *
     * @param val The value.
     */
    void setValue(T val) {
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_SET_VALUE_UNROLL_FACTOR)
            _data[i] = val;
        }
    }

    /**
     * @brief Set all elements of the matrix to zero.
     *
     */
    void setZero() { setValue(static_cast<T>(0)); }

    /**
     * @brief Take a column of a matrix by index.
     *
     * @details The result is stored to 'this'.
     *          You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @tparam M The original matrix type.
     * @tparam _unused (unused)
     * @tparam T2 The original matrix element type.
     * @param mat The original matrix.
     * @param c The column index.
     * @return (Mat&) The certain column(a column vector) (a reference to 'this') .
     */
    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    MatView& col(const M<T2, rows_, cols_, type2, _unused...>& mat, size_t c) {
        assert(c < n_cols && "Take the specific col by index requires 'The index should be smaller than the number of "
                             "the matrix's columns.'.");
        static_assert(size() == rows_, "Element number should be rows_ in Mat::col(mat, index).");
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            _data[i] = mat(i, c);
        }
        return *this;
    }

    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2, size_t rows_, size_t cols_>
    void col(size_t c, const M<T2, rows_, cols_, type2, _unused...>& mat) {
        assert(c < n_cols && "Take the specific col by index requires 'The index should be smaller than the number of "
                             "the matrix's columns.'.");
        assert(mat.size() == n_rows && "Element number should be n_rows in Mat::col(index, mat).");
        for (int i = 0; i != mat.size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            (*this)(i, c) = mat[i];
        }
    }

    /**
     * @brief Take a column of a matrix by index and make a copy.
     * @note This function makes a copy,
     *       so if you want to optimize your design,
     *       always use .col(Mat, int) that takes an argument to avoid copy in place other than initialization.
     *       They are equivalent in initialization.
     *       Another choice is to use .col_() which creates a read only view,
     *       and this operation does not copy.
     * @details You may configure macro
     *          `FLAMES_MAT_COPY_UNROLL_FACTOR` or
     *          `FLAMES_UNROLL_FACTOR` to do the operation in parallel.
     * @param c The column index.
     * @return (Mat<T, n_rows, 1, MatType::NORMAL>)(column vector) The certain column vector copy.
     */
    Mat<T, n_rows, 1, MatType::NORMAL> col(size_t c) const {
        assert(c < n_cols && "Take the specific col by index requires 'The index should be smaller than the number of "
                             "the matrix's columns.'.");
        Mat<T, n_rows, 1, MatType::NORMAL> mat;
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat(i, 0) = (*this)(i, c);
        }
        return mat;
    }

    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    MatView& add(const M<T2, n_rows, n_cols, type2, _unused...>& mat_R) {
        FLAMES_PRAGMA(INLINE)
        // return this->sub(*this, mat_R);
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            _data[i] += mat_R[i];
        }
        return *this;
    }

    template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T2,
              MatType type2>
    MatView& sub(const M<T2, n_rows, n_cols, type2, _unused...>& mat_R) {
        FLAMES_PRAGMA(INLINE)
        // return this->sub(*this, mat_R);
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_MINUS_UNROLL_FACTOR)
            _data[i] -= mat_R[i];
        }
        return *this;
    }

    MatViewT<T, n_cols, n_rows, type> t_() const { return _data; }

    template <typename Tp = T>
    Tp power() const {
        Tp p = 0;
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_POWER_UNROLL_FACTOR)
            p += static_cast<Tp>(_data[i] * _data[i]);
        }
        return p;
    }

    template <typename Tp = T>
    Tp abssum() const {
        Tp p = 0;
        for (size_t i = 0; i != size(); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_POWER_UNROLL_FACTOR)
            auto d = _data[i];
            if (d > 0) p += static_cast<Tp>(d);
        }
        return p;
    }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (Mat<T, n_rows, n_cols, type>) The real Mat.
     */
    operator Mat<T, n_rows, n_cols, type>() const {
        FLAMES_PRAGMA(INLINE);
        return Mat<T, n_rows, n_cols, type>(const_cast<const T*>(_data), InitAfterwards::NONE);
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, n_rows, n_cols, type>) The real Mat.
     */
    Mat<T, n_rows, n_cols, type> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, n_rows, n_cols, type>>(*this);
    }

    //   public: // original private
  public:
    T* const _data;
};

template <typename T, size_t n_rows, size_t n_cols, MatType type>
class MatViewOpp {
  public:
    /**
     * @brief Construct a new MatViewOpp object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewOpp(const Mat<T, n_cols, n_rows, type>& m) : _data(m.rawDataPtr()) {}

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewOpp object.
     */
    MatViewOpp(const MatViewOpp& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        FLAMES_PRAGMA(INLINE)
        assert(r < n_rows && "Matrix row index should be within range");
        assert(c < n_cols && "Matrix col index should be within range");
        if (type == MatType::NORMAL) {
            return -_data[r * n_cols + c];
        } else if (type == MatType::DIAGONAL) {
            if (r == c) return -_data[r];
            else return T(0);
        } else if (type == MatType::SCALAR) {
            if (r == c) return -_data[0];
            else return T(0);
        } else if (type == MatType::UPPER) {
            if (r <= c) return -_data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return T(0);
        } else if (type == MatType::LOWER) {
            if (r >= c) return -_data[(1 + r) * r / 2 + c];
            else return T(0);
        } else if (type == MatType::SUPPER) {
            if (r < c) return -_data[(2 * n_cols + 1 - r) * r / 2 + c - 2 * r - 1];
            else return T(0);
        } else if (type == MatType::SLOWER) {
            if (r >= c) return -_data[(1 + r) * r / 2 + c - r];
            else return T(0);
        } else if (type == MatType::SYM) {
            if (r <= c) return -_data[(2 * n_cols + 1 - r) * r / 2 + c - r];
            else return -_data[(2 * n_cols + 1 - c) * c / 2 + r - c];
        } else if (type == MatType::ASYM) {
            if (r < c) return -_data[(2 * n_cols + 1 - r) * r / 2 + c - r * 2 - 1];
            else if (r > c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - c * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const {
        FLAMES_PRAGMA(INLINE)
        assert(index < size() && "[index] should be in range in MatViewOpp");
        return -_data[index];
    }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (Mat<T, n_rows, n_cols, type>) The real Mat.
     */
    operator Mat<T, n_rows, n_cols, type>() const {
        FLAMES_PRAGMA(INLINE);
        return Mat<T, n_rows, n_cols, type>(_data, InitAfterwards::OPP);
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, n_rows, n_cols, type>) The real Mat.
     */
    Mat<T, n_rows, n_cols, type> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, n_rows, n_cols, type>>(*this);
    }

  public: // original private
    const T* _data;
};

template <typename T, size_t n_rows, size_t n_cols, MatType type>
class MatViewT {
  public:
    /**
     * @brief Construct a new MatViewT object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewT(const Mat<T, n_cols, n_rows, tType(type)>& m) : _data(m.rawDataPtr()) {}
    MatViewT(Mat<T, n_cols, n_rows, tType(type)>& m) : _data(m.rawDataPtr()) {}

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewT object.
     */
    MatViewT(const MatViewT& m) : _data(m._data) {}

    // MatViewT(T* const ptr) : _data(ptr) {}

    MatViewT(T* ptr) : _data(ptr) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        FLAMES_PRAGMA(INLINE)
        assert(r < n_rows && "Matrix row index should be within range");
        assert(c < n_cols && "Matrix col index should be within range");
        if (type == MatType::NORMAL) {
            return _data[c * n_rows + r];
        } else if (type == MatType::DIAGONAL) {
            if (r == c) return _data[r];
            else return T(0);
        } else if (type == MatType::SCALAR) {
            if (r == c) return _data[0];
            else return T(0);
        } else if (type == MatType::UPPER) {
            if (r <= c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - c];
            else return T(0);
        } else if (type == MatType::LOWER) {
            if (r >= c) return _data[(1 + c) * c / 2 + r];
            else return T(0);
        } else if (type == MatType::SUPPER) {
            if (r < c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - 2 * c - 1];
            else return T(0);
        } else if (type == MatType::SLOWER) {
            if (r >= c) return _data[(1 + c) * c / 2 + r - c];
            else return T(0);
        } else if (type == MatType::SYM) {
            if (r <= c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - c];
            else return _data[(2 * n_cols + 1 - r) * r / 2 + c - r];
        } else if (type == MatType::ASYM) {
            if (r < c) return -_data[(2 * n_cols + 1 - r) * r / 2 + c - r * 2 - 1];
            else if (r > c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - c * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const {
        FLAMES_PRAGMA(INLINE)
        size_t r, c; // of the original matrix
        if (type == MatType::DIAGONAL) {
            r = c = index;
        } else if (type == MatType::SCALAR) {
            r = c = index;
        } else if (type == MatType::NORMAL) {
            r = index / n_cols;
            c = index % n_cols;
        } else if (type == MatType::UPPER) {
            r = upperRow(index, n_cols);
            c = r + index - (2 * n_cols + 1 - r) * r / 2;
        } else if (type == MatType::LOWER) {
            r = lowerRow(index, n_cols);
            c = index - (r + 1) * r / 2;
        } else if (type == MatType::SUPPER) {
            r = supperRow(index, n_cols);
            c = 2 * r + 1 + index - c - (2 * n_cols + 1 - r) * r / 2;
        } else if (type == MatType::SLOWER) {
            r = slowerRow(index, n_cols);
            c = index - (1 + r) * r / 2 + r;
        } else if (type == MatType::SYM) {
            r = upperRow(index, n_cols);
            c = r + index - (2 * n_cols + 1 - r) * r / 2;
        } else if (type == MatType::ASYM) {
            r = supperRow(index, n_cols);
            c = 2 * r + 1 + index - c - (2 * n_cols + 1 - r) * r / 2;
        }
        return (*this)(r, c);
    }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (Mat<T, n_rows, n_cols, type>) The real Mat.
     */
    operator Mat<T, n_rows, n_cols, type>() const {
        FLAMES_PRAGMA(INLINE);
        return Mat<T, n_rows, n_cols, type>(const_cast<const T*>(_data), InitAfterwards::TR);
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, n_rows, n_cols, type>) The real Mat.
     */
    Mat<T, n_rows, n_cols, type> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, n_rows, n_cols, type>>(*this);
    }

    //   public: // original private
  public:
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    T* const _data;
};

template <typename T, size_t N, size_t N_, MatType type, typename type_parent>
class MatViewDiagMat {
  public:
    /**
     * @brief Construct a new MatViewDiagMat object from a raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewDiagMat(const Mat<T, N, N, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(type == MatType::DIAGONAL && N == N_, "DiagMat is and comes from a square matrix.");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewDiagMat object.
     */
    MatViewDiagMat(const MatViewDiagMat& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? N * N
               : type == MatType::DIAGONAL ? N
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (N - 1) * N / 2
               : type == MatType::SLOWER   ? (N - 1) * N / 2
               : type == MatType::ASYM     ? (N - 1) * N / 2
                                           : (1 + N) * N / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        if (r != c) return T(0);
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[r * N + c];
        } else if (p_type == MatType::DIAGONAL) {
            return _data[r];
        } else if (p_type == MatType::SCALAR) {
            return _data[0];
        } else if (p_type == MatType::UPPER) {
            return _data[(2 * N + 1 - r) * r / 2];
        } else if (p_type == MatType::LOWER) {
            return _data[(1 + r) * r / 2 + r];
        } else if (p_type == MatType::SUPPER) {
            return T(0);
        } else if (p_type == MatType::SLOWER) {
            return T(0);
        } else if (p_type == MatType::SYM) {
            return _data[(2 * N + 1 - r) * r / 2];
        } else if (p_type == MatType::ASYM) {
            return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(index, index); }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (Mat<T, N, N, MatType::DIAGONAL>) The real Mat.
     */
    operator Mat<T, N, N, MatType::DIAGONAL>() const {
        if (pType() == MatType::DIAGONAL) {
            return this->_data;
        } else {
            Mat<T, N, N, MatType::DIAGONAL> mat;
        MAT_COPY_DIAG:
            for (size_t i = 0; i != N; ++i) {
                FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
                mat[i] = (*this)[i];
            }
            return mat;
        }
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, N, N, MatType::DIAGONAL>) The real Mat.
     */
    Mat<T, N, N, MatType::DIAGONAL> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, N, N, MatType::DIAGONAL>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <typename T, size_t N, size_t n_cols, MatType type, typename type_parent>
class MatViewDiagVec {
  public:
    /**
     * @brief Construct a new MatViewDiagVec object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewDiagVec(const Mat<T, N, N, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(n_cols == 1, "DiagVec is a column vector.");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewDiagVec object.
     */
    MatViewDiagVec(const MatViewDiagVec& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? N * N
               : type == MatType::DIAGONAL ? N
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (N - 1) * N / 2
               : type == MatType::SLOWER   ? (N - 1) * N / 2
               : type == MatType::ASYM     ? (N - 1) * N / 2
                                           : (1 + N) * N / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(c == 0 && "Column vector's column index should always be 0.");
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[r * N + r];
        } else if (p_type == MatType::DIAGONAL) {
            return _data[r];
        } else if (p_type == MatType::SCALAR) {
            return _data[0];
        } else if (p_type == MatType::UPPER) {
            return _data[(2 * N + 1 - r) * r / 2];
        } else if (p_type == MatType::LOWER) {
            return _data[(1 + r) * r / 2 + r];
        } else if (p_type == MatType::SUPPER) {
            return T(0);
        } else if (p_type == MatType::SLOWER) {
            return T(0);
        } else if (p_type == MatType::SYM) {
            return _data[(2 * N + 1 - r) * r / 2];
        } else if (p_type == MatType::ASYM) {
            return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(index, 0); }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (Vec<T, N>) The real Mat.
     */
    operator Vec<T, N>() const {
        Vec<T, N> mat;
    MAT_COPY_DIAG:
        for (size_t i = 0; i != N; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Vec<T, N>) The real Mat.
     */
    Vec<T, N> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Vec<T, N>>(*this);
    }

    /**
     * @brief Explicitly make a Mat (Vec) copy.
     *
     * @return (Vec<T, N>) The real Mat.
     */
    Vec<T, N> asVec() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Vec<T, N>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <typename T, size_t n_rows, size_t N, MatType type, typename type_parent>
class MatViewDiagRowVec {
  public:
    /**
     * @brief Construct a new MatViewDiagRowVec object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewDiagRowVec(const Mat<T, N, N, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(n_rows == 1, "DiagRowVec is a row vector.");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewDiagRowVec object.
     */
    MatViewDiagRowVec(const MatViewDiagRowVec& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? N * N
               : type == MatType::DIAGONAL ? N
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (N - 1) * N / 2
               : type == MatType::SLOWER   ? (N - 1) * N / 2
               : type == MatType::ASYM     ? (N - 1) * N / 2
                                           : (1 + N) * N / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(r == 0 && "Row vector's row index should always be 0.");
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[c * N + c];
        } else if (p_type == MatType::DIAGONAL) {
            return _data[c];
        } else if (p_type == MatType::SCALAR) {
            return _data[0];
        } else if (p_type == MatType::UPPER) {
            return _data[(2 * N + 1 - c) * c / 2];
        } else if (p_type == MatType::LOWER) {
            return _data[(1 + c) * c / 2 + c];
        } else if (p_type == MatType::SUPPER) {
            return T(0);
        } else if (p_type == MatType::SLOWER) {
            return T(0);
        } else if (p_type == MatType::SYM) {
            return _data[(2 * N + 1 - c) * c / 2];
        } else if (p_type == MatType::ASYM) {
            return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(0, index); }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (RowVec<T, N>) The real Mat.
     */
    operator RowVec<T, N>() const {
        RowVec<T, N> mat;
    MAT_COPY_DIAG:
        for (size_t i = 0; i != N; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (RowVec<T, N>) The real Mat.
     */
    RowVec<T, N> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<RowVec<T, N>>(*this);
    }

    /**
     * @brief Explicitly make a Mat (RowVec) copy.
     *
     * @return (RowVec<T, N>) The real Mat.
     */
    RowVec<T, N> asRowVec() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<RowVec<T, N>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <typename T, size_t N, size_t N_, MatType type, typename type_parent>
class MatViewOffDiag {
  public:
    /**
     * @brief Construct a new MatViewOffDiag object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewOffDiag(const Mat<T, N, N, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(type == MatType::NORMAL && N == N_, "OffDiag is and comes from a square and normal matrix.");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewOffDiag object.
     */
    MatViewOffDiag(const MatViewOffDiag& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? N * N
               : type == MatType::DIAGONAL ? N
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (N - 1) * N / 2
               : type == MatType::SLOWER   ? (N - 1) * N / 2
               : type == MatType::ASYM     ? (N - 1) * N / 2
                                           : (1 + N) * N / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        if (r == c) return T(0);
        constexpr MatType p_type = pType();
        return (*this)(r, c);
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const {
        if (index % (N + 1) == 0) return T(0);
        else return this->_data[index];
    }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (Mat<T, N, N, MatType::NORMAL>) The real Mat.
     */
    operator Mat<T, N, N, MatType::NORMAL>() const {
        Mat<T, N, N, MatType::NORMAL> mat;
    MAT_COPY_OFFDIAG:
        for (size_t i = 0; i != N * N; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            if (i % (N + 1) == 0) mat[i] = 0;
            else mat[i] = this->_m[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, N, N, MatType::NORMAL>) The real Mat.
     */
    Mat<T, N, N, MatType::NORMAL> asMat() const {
        FLAMES_PRAGMA(INLINE)
        return static_cast<Mat<T, N, N, MatType::NORMAL>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <size_t Index, typename T, size_t n_rows, size_t n_cols, MatType type, typename type_parent>
class MatViewCol {
  public:
    /**
     * @brief Construct a new MatViewCol object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewCol(const Mat<T, n_rows, n_cols, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(Index < n_cols, "Take the specific col by index requires 'The index should be smaller than the "
                                      "number of the matrix's cols.'");
        static_assert(Index > 0, "Take the specific col by index requires 'The index can't be smaller than 0.'");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewCol object.
     */
    MatViewCol(const MatViewCol& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(c == 0 && "Column vector's column index should always be 0.");
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[r * n_cols + Index];
        } else if (p_type == MatType::DIAGONAL) {
            if (r == Index) return _data[r];
            else return T(0);
        } else if (p_type == MatType::SCALAR) {
            if (r == Index) return _data[0];
            else return T(0);
        } else if (p_type == MatType::UPPER) {
            if (r <= Index) return _data[(2 * n_cols + 1 - r) * r / 2 + Index - r];
            else return T(0);
        } else if (p_type == MatType::LOWER) {
            if (r >= Index) return _data[(1 + r) * r / 2 + Index];
            else return T(0);
        } else if (p_type == MatType::SUPPER) {
            if (r < Index) return _data[(2 * n_cols + 1 - r) * r / 2 + Index - 2 * r - 1];
            else return T(0);
        } else if (p_type == MatType::SLOWER) {
            if (r >= Index) return _data[(1 + r) * r / 2 + Index - r];
            else return T(0);
        } else if (p_type == MatType::SYM) {
            if (r <= Index) return _data[(2 * n_cols + 1 - r) * r / 2 + Index - r];
            else return _data[(2 * n_cols + 1 - Index) * Index / 2 + r - Index];
        } else if (p_type == MatType::ASYM) {
            if (r < Index) return -_data[(2 * n_cols + 1 - r) * r / 2 + Index - r * 2 - 1];
            else if (r > Index) return _data[(2 * n_cols + 1 - Index) * Index / 2 + r - Index * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(index, 0); }

    /**
     * @brief Conversion from view to a real Mat.
     * @return (RowVec<T, n_rows>) The real Mat.
     */
    operator Vec<T, n_rows>() const {
        Vec<T, n_rows> mat;
    MAT_COPY_COL:
        for (size_t i = 0; i != n_rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Vec<T, n_rows>) The real Mat.
     */
    Vec<T, n_rows> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Vec<T, n_rows>>(*this);
    }

    /**
     * @brief Explicitly make a Mat (Vec) copy.
     *
     * @return (Vec<T, n_rows>) The real Mat.
     */
    Vec<T, n_rows> asVec() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Vec<T, n_rows>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <size_t Index, typename T, size_t n_rows, size_t n_cols, MatType type, typename type_parent>
class MatViewRow {
  public:
    /**
     * @brief Construct a new MatViewRow object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewRow(const Mat<T, n_cols, n_rows, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(Index < n_rows, "Take the specific row by index requires 'The index should be smaller than the "
                                      "number of the matrix's rows.'");
        static_assert(Index > 0, "Take the specific row by index requires 'The index can't be smaller than 0.'");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewRow object.
     */
    MatViewRow(const MatViewRow& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(r == 0 && "Row vector's row index should always be 0.");
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[Index * n_cols + c];
        } else if (p_type == MatType::DIAGONAL) {
            if (Index == c) return _data[Index];
            else return T(0);
        } else if (p_type == MatType::SCALAR) {
            if (Index == c) return _data[0];
            else return T(0);
        } else if (p_type == MatType::UPPER) {
            if (Index <= c) return _data[(2 * n_cols + 1 - Index) * Index / 2 + c - Index];
            else return T(0);
        } else if (p_type == MatType::LOWER) {
            if (Index >= c) return _data[(1 + Index) * Index / 2 + c];
            else return T(0);
        } else if (p_type == MatType::SUPPER) {
            if (Index < c) return _data[(2 * n_cols + 1 - Index) * Index / 2 + c - 2 * Index - 1];
            else return T(0);
        } else if (p_type == MatType::SLOWER) {
            if (Index >= c) return _data[(1 + Index) * Index / 2 + c - Index];
            else return T(0);
        } else if (p_type == MatType::SYM) {
            if (Index <= c) return _data[(2 * n_cols + 1 - Index) * Index / 2 + c - Index];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + Index - c];
        } else if (p_type == MatType::ASYM) {
            if (Index < c) return -_data[(2 * n_cols + 1 - Index) * Index / 2 + c - Index * 2 - 1];
            else if (Index > c) return _data[(2 * n_cols + 1 - c) * c / 2 + Index - c * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(0, index); }

    /**
     * @brief Conversion from view to a real Mat.
     *
     * @return (RowVec<T, n_cols>) The real Mat.
     */
    operator RowVec<T, n_cols>() const {
        RowVec<T, n_cols> mat;
    MAT_COPY_ROW:
        for (size_t i = 0; i != n_cols; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (RowVec<T, n_cols>) The real Mat.
     */
    RowVec<T, n_cols> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<RowVec<T, n_cols>>(*this);
    }

    /**
     * @brief Explicitly make a Mat (RowVec) copy.
     *
     * @return (RowVec<T, n_cols>) The real Mat.
     */
    RowVec<T, n_cols> asRowVec() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<RowVec<T, n_cols>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <size_t first_col, size_t last_col, typename T, size_t n_rows, size_t n_cols, MatType type,
          typename type_parent>
class MatViewCols {
  public:
    /**
     * @brief Construct a new MatViewCols object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewCols(const Mat<T, n_rows, n_cols, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(last_col < n_cols, "Take the successive cols by indexes requires 'The indexes should be smaller "
                                         "than the number of the matrix's cols.'");
        static_assert(first_col > 0, "Take the specific col by index requires 'The indexes can't be smaller than 0.'");
        static_assert(first_col <= last_col, "The first index should be smaller than the second one.");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewCols object.
     */
    MatViewCols(const MatViewCols& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(c < (last_col - first_col + 1) &&
               "The col index must be small than the number of the successive columns .");
        assert(r > 0 && "The row index can't be smaller than 0.");
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[r * n_cols + c + first_col];
        } else if (p_type == MatType::DIAGONAL) {
            if (r == (c + first_col)) return _data[r];
            else return T(0);
        } else if (p_type == MatType::SCALAR) {
            if (r == (c + first_col)) return _data[0];
            else return T(0);
        } else if (p_type == MatType::UPPER) {
            if (r <= (c + first_col)) return _data[(2 * n_cols + 1 - r) * r / 2 + c + first_col - r];
            else return T(0);
        } else if (p_type == MatType::LOWER) {
            if (r >= (c + first_col)) return _data[(1 + r) * r / 2 + c + first_col - r];
            else return T(0);
        } else if (p_type == MatType::SUPPER) {
            if (r < (c + first_col)) return _data[(2 * n_cols + 1 - r) * r / 2 + c + first_col - 2 * r - 1];
            else return T(0);
        } else if (p_type == MatType::SLOWER) {
            if (r >= (c + first_col)) return _data[(1 + r) * r / 2 + c + first_col - r];
            else return T(0);
        } else if (p_type == MatType::SYM) {
            if (r <= (c + first_col)) return _data[(2 * n_cols + 1 - r) * r / 2 + c + first_col - r];
            else return _data[(2 * n_cols + 1 - (c + first_col)) * (c + first_col) / 2 + r - (c + first_col)];
        } else if (p_type == MatType::ASYM) {
            if (r < (c + first_col)) return -_data[(2 * n_cols + 1 - r) * r / 2 + (c + first_col) - r * 2 - 1];
            else if (r > (c + first_col))
                return _data[(2 * n_cols + 1 - (c + first_col)) * (c + first_col) / 2 + r - (c + first_col) * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const {
        return (*this)(index / (last_col - first_col + 1), index % (last_col - first_col + 1));
    }

    /**
     * @brief Conversion from view to a real Mat.
     * @return (Mat<T, n_rows, last_col - first_col + 1, MatType::NORMAL>) The real Mat.
     */
    operator Mat<T, n_rows, last_col - first_col + 1, MatType::NORMAL>() const {
        Mat<T, n_rows, last_col - first_col + 1, MatType::NORMAL> mat;
    MAT_COPY_COLS:
        for (size_t i = 0; i != n_rows * (last_col - first_col + 1); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, n_rows, last_col - first_col + 1, MatType::NORMAL>) The real Mat.
     */
    Mat<T, n_rows, last_col - first_col + 1, MatType::NORMAL> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, n_rows, last_col - first_col + 1, MatType::NORMAL>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <size_t first_row, size_t last_row, typename T, size_t n_rows, size_t n_cols, MatType type,
          typename type_parent>
class MatViewRows {
  public:
    /**
     * @brief Construct a new MatViewRows object from raw data pointer.
     *
     * @param m The original matrix.
     */
    MatViewRows(const Mat<T, n_cols, n_rows, matType<type_parent>()>& m) : _data(m.rawDataPtr()) {
        static_assert(last_row < n_rows, "Take the successive rows by indexes requires 'The indexes should be smaller "
                                         "than the number of the matrix's rows.'");
        static_assert(first_row > 0, "Take the specific row by index requires 'The indexes can't be smaller than 0.'");
        static_assert(first_row <= last_row, "The first index should be smaller than the second one.");
    }

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewRows object.
     */
    MatViewRows(const MatViewRows& m) : _data(m._data) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Parent matrix as MatType.last_col
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(r < last_row - first_row && "The row index must be small than the number of the successive rows .");
        assert(c > 0 && "The row index can't be smaller than 0.");
        constexpr MatType p_type = pType();
        if (p_type == MatType::NORMAL) {
            return _data[(r + first_row) * n_cols + c];
        } else if (p_type == MatType::DIAGONAL) {
            if ((r + first_row) == c) return _data[r];
            else return T(0);
        } else if (p_type == MatType::SCALAR) {
            if ((r + first_row) == c) return _data[0];
            else return T(0);
        } else if (p_type == MatType::UPPER) {
            if ((r + first_row) <= c)
                return _data[(2 * n_cols + 1 - (r + first_row)) * (r + first_row) / 2 + c - (r + first_row)];
            else return T(0);
        } else if (p_type == MatType::LOWER) {
            if ((r + first_row) >= c) return _data[(1 + (r + first_row)) * (r + first_row) / 2 + c];
            else return T(0);
        } else if (p_type == MatType::SUPPER) {
            if ((r + first_row) < c)
                return _data[(2 * n_cols + 1 - (r + first_row)) * (r + first_row) / 2 + c - 2 * (r + first_row) - 1];
            else return T(0);
        } else if (p_type == MatType::SLOWER) {
            if ((r + first_row) >= c) return _data[(1 + (r + first_row)) * (r + first_row) / 2 + c - (r + first_row)];
            else return T(0);
        } else if (p_type == MatType::SYM) {
            if ((r + first_row) <= c)
                return _data[(2 * n_cols + 1 - (r + first_row)) * (r + first_row) / 2 + c - (r + first_row)];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + (r + first_row) - c];
        } else if (p_type == MatType::ASYM) {
            if ((r + first_row) < c)
                return -_data[(2 * n_cols + 1 - (r + first_row)) * (r + first_row) / 2 + c - (r + first_row) * 2 - 1];
            else if ((r + first_row) > c) return _data[(2 * n_cols + 1 - c) * c / 2 + r - c * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(index / n_cols, index % n_cols); }

    /**
     * @brief Conversion from view to a real Mat.
     * @return Mat<T, last_row - first_row, n_cols, MatType::NORMAL> The real Mat.
     */
    operator Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL>() const {
        Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL> mat;
    MAT_COPY_ROWS:
        for (size_t i = 0; i != n_rows * (last_row - first_row + 1); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL>)The real Mat.
     */
    Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL>>(*this);
    }

    /**
     * @brief Explicitly make a Mat (Vec) copy.
     *
     * @return (Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL>) The real Mat.
     */
    Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL> asVec() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, last_row - first_row + 1, n_cols, MatType::NORMAL>>(*this);
    }

  public: // original private
    /**
     * @brief Raw data pointer.
     *
     * @note This contents will not be modified.
     */
    const T* _data;
};

template <size_t Cols, typename T, typename M, size_t n_rows, size_t n_cols, MatType type, typename type_parent>
class MatViewColsContainer {
  public:
    /**
     * @brief Construct a new MatViewColsContainer object from a matrix raw data pointer and a container data pointer.
     *
     * @param m The original matrix.
     */
    MatViewColsContainer(const Mat<T, n_rows, n_cols, matType<type_parent>()>& m, const M Container)
        : _data(m.rawDataPtr()), container(Container) {}

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewColsContainer object.
     */
    MatViewColsContainer(const MatViewColsContainer& m) : _data(m._data), container(m.container) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(c < Cols && "The col index must be small than the number of the discrete columns .");
        assert(r < n_rows && "The row index should be smaller than the number of the matrix's rows.");
        assert(container[c] < n_cols && "The container index should be small than the number of the matrix's col.");
        constexpr MatType p_type = pType();
        if (type == MatType::NORMAL) {
            return _data[r * n_cols + container[c]];
        } else if (type == MatType::DIAGONAL) {
            if (r == container[c]) return _data[r];
            else return T(0);
        } else if (type == MatType::SCALAR) {
            if (r == container[c]) return _data[0];
            else return T(0);
        } else if (type == MatType::UPPER) {
            if (r <= container[c]) return _data[(2 * n_cols + 1 - r) * r / 2 + container[c] - r];
            else return T(0);
        } else if (type == MatType::LOWER) {
            if (r >= container[c]) return _data[(1 + r) * r / 2 + container[c]];
            else return T(0);
        } else if (type == MatType::SUPPER) {
            if (r < container[c]) return _data[(2 * n_cols + 1 - r) * r / 2 + container[c] - 2 * r - 1];
            else return T(0);
        } else if (type == MatType::SLOWER) {
            if (r >= container[c]) return _data[(1 + r) * r / 2 + container[c] - r];
            else return T(0);
        } else if (type == MatType::SYM) {
            if (r <= container[c]) return _data[(2 * n_cols + 1 - r) * r / 2 + container[c] - r];
            else return _data[(2 * n_cols + 1 - container[c]) * container[c] / 2 + r - container[c]];
        } else if (type == MatType::ASYM) {
            if (r < container[c]) return _data[(2 * n_cols + 1 - r) * r / 2 + container[c] - r * 2 - 1];
            else if (r > container[c])
                return -_data[(2 * n_cols + 1 - container[c]) * container[c] / 2 + r - container[c] * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(index / Cols, index % Cols); }

    /**
     * @brief Conversion from view to a real Mat.
     * @return (Mat<T, n_rows, Cols, MatType::NORMAL>) The real Mat.
     */
    operator Mat<T, n_rows, Cols, MatType::NORMAL>() const {
        Mat<T, n_rows, Cols, MatType::NORMAL> mat;
    MAT_COPY_COLS:
        for (size_t i = 0; i != n_rows * (Cols); ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, n_rows, Cols, MatType::NORMAL>) The real Mat.
     */
    Mat<T, n_rows, Cols, MatType::NORMAL> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, n_rows, Cols, MatType::NORMAL>>(*this);
    }

  public: // original private
    /**
     * @brief Matrix raw data pointer.
     * @brief Container.
     * @note This contents will not be modified.
     */
    const T* _data;
    const M container;
};

template <size_t Rows, typename T, typename M, size_t n_rows, size_t n_cols, MatType type, typename type_parent>
class MatViewRowsContainer {
  public:
    /**
     * @brief Construct a new MatViewRowsContainer object from a matrix raw data pointer and a container data pointer.
     *
     * @param m The original matrix.
     */
    MatViewRowsContainer(const Mat<T, n_rows, n_cols, matType<type_parent>()>& m, const M Container)
        : _data(m.rawDataPtr()), container(Container) {}

    /**
     * @brief Copy constructor.
     *
     * @param m Another MatViewRowsContainer object.
     */
    MatViewRowsContainer(const MatViewRowsContainer& m) : _data(m._data), container(m.container) {}

    /**
     * @brief The data element number.
     *
     * @return (constexpr size_t) The size.
     */
    inline static constexpr size_t size() noexcept {
        return type == MatType::NORMAL     ? n_rows * n_cols
               : type == MatType::DIAGONAL ? n_rows
               : type == MatType::SCALAR   ? 1
               : type == MatType::SUPPER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::SLOWER   ? (n_rows - 1) * n_rows / 2
               : type == MatType::ASYM     ? (n_rows - 1) * n_rows / 2
                                           : (1 + n_rows) * n_rows / 2;
    }

    /**
     * @brief Parent matrix as MatType.
     *
     * @return (constexpr MatType) The MatType.
     */
    inline static constexpr MatType pType() noexcept { return matType<type_parent>(); }

    /**
     * @brief Get the read only data element from row and column index.
     *
     * @param r The row index.
     * @param c The column index.
     * @return (T) The element value.
     */
    T operator()(size_t r, size_t c) const {
        assert(c < n_cols && "The col index must be small than the number of the matrix's column.");
        assert(r < Rows && "The row index should be smaller than the number of the discrete rows .");
        assert(container[r] < n_rows && "The container index should be small than the number of the matrix's row ");
        constexpr MatType p_type = pType();
        if (type == MatType::NORMAL) {
            return _data[container[r] * n_cols + c];
        } else if (type == MatType::DIAGONAL) {
            if (container[r] == c) return _data[c];
            else return T(0);
        } else if (type == MatType::SCALAR) {
            if (container[r] == c) return _data[0];
            else return T(0);
        } else if (type == MatType::UPPER) {
            if (container[r] <= c) return _data[(2 * n_cols + 1 - container[r]) * container[r] / 2 + c - container[r]];
            else return T(0);
        } else if (type == MatType::LOWER) {
            if (container[r] >= c) return _data[(1 + container[r]) * container[r] / 2 + c];
            else return T(0);
        } else if (type == MatType::SUPPER) {
            if (container[r] < c)
                return _data[(2 * n_cols + 1 - container[r]) * container[r] / 2 + c - 2 * container[r] - 1];
            else return T(0);
        } else if (type == MatType::SLOWER) {
            if (container[r] >= c) return _data[(1 + container[r]) * container[r] / 2 + c - r];
            else return T(0);
        } else if (type == MatType::SYM) {
            if (container[r] <= c) return _data[(2 * n_cols + 1 - container[r]) * container[r] / 2 + c - container[r]];
            else return _data[(2 * n_cols + 1 - c) * c / 2 + container[r] - c];
        } else if (type == MatType::ASYM) {
            if (container[r] < c)
                return _data[(2 * n_cols + 1 - container[r]) * container[r] / 2 + c - container[r] * 2 - 1];
            else if (container[r] > c) return -_data[(2 * n_cols + 1 - c) * c / 2 + container[r] - c * 2 - 1];
            else return T(0);
        } else {
            // Normally it is impossible to reach here.
            assert(!"Impossible! Unknown MatType!");
        }
    }

    /**
     * @brief Get the read only element by array row major index.
     *
     * @param index The index.
     * @return (T) The data.
     */
    T operator[](size_t index) const { return (*this)(index / n_cols, index % n_cols); }

    /**
     * @brief Conversion from view to a real Mat.
     * @return (Mat<T, Rows, n_cols, MatType::NORMAL>) The real Mat.
     */
    operator Mat<T, Rows, n_cols, MatType::NORMAL>() const {
        Mat<T, Rows, n_cols, MatType::NORMAL> mat;
    MAT_COPY_ROWS:
        for (size_t i = 0; i != n_cols * Rows; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            mat[i] = (*this)[i];
        }
        return mat;
    }

    /**
     * @brief Explicitly make a Mat copy.
     *
     * @return (Mat<T, Rows, n_cols, MatType::NORMAL>) The real Mat.
     */
    Mat<T, Rows, n_cols, MatType::NORMAL> asMat() const {
        FLAMES_PRAGMA(INLINE);
        return static_cast<Mat<T, Rows, n_cols, MatType::NORMAL>>(*this);
    }

  public: // original private
    /**
     * @brief Matrix raw data pointer.
     * @brief Container.
     * @note This contents will not be modified.
     */
    const T* _data;
    const M container;
};

/**
 * @brief Add two matrices and make a copy.
 *
 * @details This will call .add() function.
 * @note This function makes a copy so it should only by used for initialization.
 *       Otherwise use .add(Mat_L, mat_R) to avoid the copy operation.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam type1 The left matrix MatType.
 * @tparam type2 The right matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<T1, n_rows, n_cols, sumType(type1, type2)>) The addition result as a copy.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type1, MatType type2>
static inline Mat<T1, n_rows, n_cols, sumType(type1, type2)>
operator+(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
          const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    Mat<T1, n_rows, n_cols, type1> mat;
    return mat.add(mat_L, mat_R);
}

/**
 * @brief Matrix self plus a matrix.
 *
 * @tparam M The right side matrix type.
 * @tparam _unused (unused)
 * @tparam T1 The left side matrix element type.
 * @tparam T2 The right side matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type1 The left side matrix MatType.
 * @tparam type2 The right side matrix MatType.
 * @param mat_L The left side matrix.
 * @param mat_R The right side matrix.
 * @return (Mat<T1, n_rows, n_cols, type1>&) The addition result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T1, typename T2,
          size_t n_rows, size_t n_cols, MatType type1, MatType type2>
static inline Mat<T1, n_rows, n_cols, type1>& operator+=(Mat<T1, n_rows, n_cols, type1>& mat_L,
                                                         const M<T2, n_rows, n_cols, type2, _unused...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    return mat_L.add(mat_R);
}

/**
 * @brief Minus two matrices and make a copy.
 *
 * @details This will call .sub() function.
 * @note This function makes a copy so it should only by used for initialization.
 *       Otherwise use .sub(Mat_L, mat_R) to avoid the copy operation.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam type1 The left matrix MatType.
 * @tparam type2 The right matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<T1, n_rows, n_cols, sumType(type1, type2)>) The subtract result as a copy.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type1, MatType type2>
static inline Mat<T1, n_rows, n_cols, sumType(type1, type2)>
operator-(const M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
          const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    Mat<T1, n_rows, n_cols, type1> mat;
    return mat.sub(mat_L, mat_R);
}

/**
 * @brief Matrix self minus a matrix.
 *
 * @tparam M The right side matrix type.
 * @tparam _unused (unused)
 * @tparam T1 The left side matrix element type.
 * @tparam T2 The right side matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type1 The left side matrix MatType.
 * @tparam type2 The right side matrix MatType.
 * @param mat_L The left side matrix.
 * @param mat_R The right side matrix.
 * @return (Mat<T1, n_rows, n_cols, type1>&) The subtract result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type1, MatType type2>
static inline M1<T1, n_rows, n_cols, type1>& operator-=(M1<T1, n_rows, n_cols, type1, _unused1...>& mat_L,
                                                        const M2<T2, n_rows, n_cols, type2, _unused2...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    return mat_L.sub(mat_R);
}

/**
 * @brief Matrix self times a matrix.
 *
 * @details This operator calls .mul() function.
 *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
 *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
 * @note This is the correct way to do self multiplication.\n
 *       The wrong way: \code a.mul(a, b); \endcode
 *       because this will lead to writing into the read only matrix,
 *       leading to a access violation and wrong results.
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam ScalarT The scalar type.
 * @tparam T2 The matrix element type.
 * @param mat The matrix.
 * @param s The scalar.
 * @return (Mat&) The multiplication result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T1, typename T2,
          size_t n_rows, size_t n_cols, MatType type1, MatType type2,
          std::enable_if_t<(std::is_same<T1, bool>::value), bool> = true>
static inline Mat<T1, n_rows, n_cols, mulType(type1, type2, n_rows, n_cols, n_cols)>
operator*=(Mat<T1, n_rows, n_cols, type1>& mat, const M<T2, n_cols, n_cols, type2, _unused...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    Mat<T2, n_rows, n_cols, mulType(type1, type2, n_rows, n_cols, n_cols)> tmp;
    tmp.mul(mat, mat_R);
    return mat = tmp;
}

/**
 * @brief Matrix self times a matrix of a different type.
 *
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type1 The left side matrix MatType.
 * @tparam type2 The right side matrix MatType.
 * @param mat The matrix.
 * @param mat_R The right matrix.
 * @return Mat<T1, n_rows, n_cols, mulType(type1, type2, n_rows, n_cols, n_cols)>
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T1, typename T2,
          size_t n_rows, size_t n_cols, MatType type1, MatType type2,
          std::enable_if_t<!(std::is_same<T1, bool>::value), bool> = true>
static inline Mat<T1, n_rows, n_cols, mulType(type1, type2, n_rows, n_cols, n_cols)>
operator*=(Mat<T1, n_rows, n_cols, type1>& mat, const M<T2, n_cols, n_cols, type2, _unused...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    Mat<T1, n_rows, n_cols, mulType(type1, type2, n_rows, n_cols, n_cols)> tmp;
    tmp.mul(mat, mat_R);
    return mat = tmp;
}

/**
 * @brief Matrix self times a scalar.
 *
 * @details This operator calls .mul() function.
 *          This scalar is C++ arithmetic types, like double and int.
 *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
 *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam ScalarT The scalar type.
 * @tparam T The matrix element type.
 * @param mat The matrix.
 * @param s The scalar.
 * @return (Mat&) The multiplication result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T,
          typename ScalarT, size_t n_rows, size_t n_cols, MatType type,
          std::enable_if_t<std::is_arithmetic<std::remove_reference_t<ScalarT>>::value, bool> = true>
static inline Mat<T, n_rows, n_cols, type>& operator*=(const M<T, n_rows, n_cols, type, _unused...>& mat, ScalarT s) {
    FLAMES_PRAGMA(INLINE)
    return mat.mul(s);
}

/**
 * @brief Matrix times a scalar.
 *
 * @details This operator calls .mul() function.
 *          This scalar is C++ arithmetic types, like double and int.
 *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
 *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
 * @note This function makes a copy so it should only by used for initialization.
 *       Otherwise use .mul(Mat_L, s) to avoid the copy operation.
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam ScalarT The scalar type.
 * @tparam T The matrix element type.
 * @param mat_L The matrix.
 * @param s The scalar.
 * @return (Mat&) The multiplication result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T,
          typename ScalarT, size_t n_rows, size_t n_cols, MatType type,
          std::enable_if_t<std::is_arithmetic<std::remove_reference_t<ScalarT>>::value, bool> = true>
static inline Mat<T, n_rows, n_cols, type> operator*(const M<T, n_rows, n_cols, type, _unused...>& mat_L, ScalarT s) {
    FLAMES_PRAGMA(INLINE)
    Mat<T, n_rows, n_cols, type> mat;
    return mat.mul(mat_L, s);
}

/**
 * @brief Matrix times an ap_fixed scalar.
 *
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam T The matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @tparam AP_W ap_int W param.
 * @tparam AP_I ap_int I param.
 * @tparam AP_Q ap_int Q param.
 * @tparam AP_O ap_int O param.
 * @tparam AP_N ap_int N param.
 * @param mat_L The left matrix.
 * @param s The float.
 * @return (Mat&) The multiplication result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T,
          size_t n_rows, size_t n_cols, MatType type, int AP_W, int AP_I, ap_q_mode AP_Q, ap_o_mode AP_O, int AP_N>
static inline Mat<T, n_rows, n_cols, type> operator*(const M<T, n_rows, n_cols, type, _unused...>& mat_L,
                                                     ap_fixed<AP_W, AP_I, AP_Q, AP_O, AP_N> s) {
    FLAMES_PRAGMA(INLINE)
    Mat<T, n_rows, n_cols, type> mat;
    return mat.mul(mat_L, s);
}

/**
 * @brief Matrix times a scalar.
 *
 * @details This operator calls .mul() function.
 *          This scalar is C++ arithmetic types, like double and int.
 *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
 *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
 * @note This function makes a copy so it should only by used for initialization.
 *       Otherwise use .mul(Mat_R, mat_s) to avoid the copy operation.
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam ScalarT The scalar type.
 * @tparam T2 The matrix element type.
 * @param mat_R The matrix.
 * @param s The scalar.
 * @return (Mat&) The multiplication result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T,
          typename ScalarT, size_t n_rows, size_t n_cols, MatType type,
          std::enable_if_t<std::is_arithmetic<std::remove_reference_t<ScalarT>>::value, bool> = true>
static inline Mat<T, n_rows, n_cols, type> operator*(ScalarT s, const M<T, n_rows, n_cols, type, _unused...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    Mat<T, n_rows, n_cols, type> mat;
    return mat.mul(mat_R, s);
}

/**
 * @brief Matrix times an ap_fixed scalar.
 *
 * @tparam M The matrix type.
 * @tparam _unused (unused)
 * @tparam T The matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @tparam AP_W ap_int W param.
 * @tparam AP_I ap_int I param.
 * @tparam AP_Q ap_int Q param.
 * @tparam AP_O ap_int O param.
 * @tparam AP_N ap_int N param.
 * @param mat_R The right matrix.
 * @param s The float.
 * @return (Mat&) The multiplication result (a reference to 'this').
 */
template <template <class, size_t, size_t, MatType, class...> typename M, typename... _unused, typename T,
          size_t n_rows, size_t n_cols, MatType type, int AP_W, int AP_I, ap_q_mode AP_Q, ap_o_mode AP_O, int AP_N>
static inline Mat<T, n_rows, n_cols, type> operator*(ap_fixed<AP_W, AP_I, AP_Q, AP_O, AP_N> s,
                                                     const M<T, n_rows, n_cols, type, _unused...>& mat_R) {
    FLAMES_PRAGMA(INLINE)
    Mat<T, n_rows, n_cols, type> mat;
    return mat.mul(mat_R, s);
}

/**
 * @brief Matrix multiplication.
 *
 * @details This operator calls .mul() function.
 *          You may configure `FLAMES_MAT_SCALAR_TIMES_UNROLL_FACTOR`
 *          or `FLAMES_UNROLL_FACTOR` to doing multiplication in parallel.
 * @note This function makes a copy so it should only by used for initialization.
 *       Otherwise use .mul(Mat_L, mat_R) to avoid the copy operation.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The row number of the left matrix.
 * @tparam comm The common number (the column number of the left matrix and the row number of the right matrix).
 * @tparam n_cols The column number of the right matrix.
 * @tparam type1 The left matrix MatType.
 * @tparam type2 The right matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<T, n_rows, n_cols, mulType(type1, type2, n_rows, comm, n_cols)>) The multiplication result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t comm, size_t n_cols, MatType type1, MatType type2,
          std::enable_if_t<(std::is_same<T1, bool>::value), bool> = true>
static inline Mat<T2, n_rows, n_cols, mulType(type1, type2, n_rows, comm, n_cols)>
operator*(const M1<T1, n_rows, comm, type1, _unused1...>& mat_L,
          const M2<T2, comm, n_cols, type2, _unused2...>& mat_R) {
    Mat<T2, n_rows, n_cols, mulType(type1, type2, n_rows, comm, n_cols)> mat;
    mat.mul(mat_L, mat_R);
    return mat;
    // if "return mat.mul(mat_L, mat_R);" , then there will be a error. But I don't know why.
}

template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t comm, size_t n_cols, MatType type1, MatType type2,
          std::enable_if_t<!(std::is_same<T1, bool>::value), bool> = true>
static inline Mat<T1, n_rows, n_cols, mulType(type1, type2, n_rows, comm, n_cols)>
operator*(const M1<T1, n_rows, comm, type1, _unused1...>& mat_L,
          const M2<T2, comm, n_cols, type2, _unused2...>& mat_R) {
    Mat<T1, n_rows, n_cols, mulType(type1, type2, n_rows, comm, n_cols)> mat;
    mat.mul(mat_L, mat_R);
    return mat;
    // if "return mat.mul(mat_L, mat_R);" , then there will be a error. But I don't know why.
}

/**
 * @brief Element-wise product of two matrices.
 *
 * @details You can configure the macro `FLAMES_MAT_EMUL_UNROLL_FACTOR` to determine the parallelism.\n
 *          This internally calls Mat::emul(mat, mat).\n
 *          The return element type is that of the left matrix.
 * @note It now only supports element-wise product of two matrices of the same dimension and MatType.\n
 *       This is not the modulus operator.Ise .mod() for the element-wise modulus operation.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<T1, n_rows, n_cols, type>) The element-wise product result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<T1, n_rows, n_cols, type> operator%(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                        const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<T1, n_rows, n_cols, type> mat;
    return mat.emul(mat_L, mat_R);
}

/**
 * @brief Element-wise equal comparison
 *
 * @details This function returns a bool matrix.\n
 *          You can configure the macro `FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR` to determine the parallelism.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<bool, n_rows, n_cols, type>) The comparison result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<bool, n_rows, n_cols, type> operator==(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                           const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<bool, n_rows, n_cols, type> mat;
OPERATOR_EQUAL:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] == mat_R[i];
    }
    return mat;
}

/**
 * @brief Element-wise unequal comparison
 *
 * @details This function returns a bool matrix.\n
 *          You can configure the macro `FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR` to determine the parallelism.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<bool, n_rows, n_cols, type>) The comparison result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<bool, n_rows, n_cols, type> operator!=(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                           const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<bool, n_rows, n_cols, type> mat;
OPERATOR_UNEQUAL:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] != mat_R[i];
    }
    return mat;
}

/**
 * @brief Element-wise greater comparison
 *
 * @details This function returns a bool matrix.\n
 *          You can configure the macro `FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR` to determine the parallelism.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<bool, n_rows, n_cols, type>) The comparison result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<bool, n_rows, n_cols, type> operator>(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                          const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<bool, n_rows, n_cols, type> mat;
OPERATOR_GREATER:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] > mat_R[i];
    }
    return mat;
}

/**
 * @brief Element-wise less comparison
 *
 * @details This function returns a bool matrix.\n
 *          You can configure the macro `FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR` to determine the parallelism.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<bool, n_rows, n_cols, type>) The comparison result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<bool, n_rows, n_cols, type> operator<(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                          const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<bool, n_rows, n_cols, type> mat;
OPERATOR_LESS:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] < mat_R[i];
    }
    return mat;
}

/**
 * @brief Element-wise greater or equal comparison
 *
 * @details This function returns a bool matrix.\n
 *          You can configure the macro `FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR` to determine the parallelism.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<bool, n_rows, n_cols, type>) The comparison result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<bool, n_rows, n_cols, type> operator>=(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                           const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<bool, n_rows, n_cols, type> mat;
OPERATOR_GEQ:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] >= mat_R[i];
    }
    return mat;
}

/**
 * @brief Element-wise less or equal comparison
 *
 * @details This function returns a bool matrix.\n
 *          You can configure the macro `FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR` to determine the parallelism.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<bool, n_rows, n_cols, type>) The comparison result.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<bool, n_rows, n_cols, type> operator<=(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                           const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<bool, n_rows, n_cols, type> mat;
OPERATOR_LEQ:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] <= mat_R[i];
    }
    return mat;
}

/**
 * @brief Element-wise modulus.
 *
 * @details This requires the operator `%` is defined between both matrix elements.
 * @note This should not be misued with `operator%`, which is the element-wise multiplication.
 * @tparam M1 The left matrix type.
 * @tparam _unused1 (unused)
 * @tparam M2 The right matrix type.
 * @tparam _unused2 (unused)
 * @tparam T1 The left matrix element type.
 * @tparam T2 The right matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param mat_L The left matrix.
 * @param mat_R The right matrix.
 * @return (Mat<T1, n_rows, n_cols, type>) The modulus result matrix.
 */
template <template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t n_rows, size_t n_cols, MatType type>
Mat<T1, n_rows, n_cols, type> mod(const M1<T1, n_rows, n_cols, type, _unused1...>& mat_L,
                                  const M2<T2, n_rows, n_cols, type, _unused2...>& mat_R) {
    Mat<T1, n_rows, n_cols, type> mat;
OPERATOR_LEQ:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_BOOL_OPER_UNROLL_FACTOR)
        mat[i] = mat_L[i] % mat_R[i]; // Operator % Should be defined, otherwise compilation error.
    }
    return mat;
}

template <typename Tp = double, template <class, size_t, size_t, MatType, class...> typename M1, typename... _unused1,
          template <class, size_t, size_t, MatType, class...> typename M2, typename... _unused2, typename T1,
          typename T2, size_t L_rows, size_t L_cols, size_t R_rows, size_t R_cols, MatType type,
          std::enable_if_t<L_rows * L_cols == R_rows * R_cols, bool> = true>
Tp innerProd(const M1<T1, L_rows, L_cols, type, _unused1...>& mat_L,
             const M2<T2, R_rows, R_cols, type, _unused2...>& mat_R) {
    assert(mat_L.size() == mat_R.size() && "Dimension should meet for innerProd.");
    Tp result = 0;
INNER_PROD:
    for (size_t i = 0; i != mat_L.size(); ++i) {
        FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_EMUL_UNROLL_FACTOR)
        result += static_cast<Tp>(mat_L[i] * mat_R[i]);
    }
    return result;
}

/**
 * @brief Print matrix in a out stream.
 *
 * @details This function calls .print() function.
 * @tparam T The matrix element type.
 * @tparam n_rows The number of rows.
 * @tparam n_cols The number of columns.
 * @tparam type The matrix MatType.
 * @param os The out stream.
 * @param mat The matrix.
 * @return (std::ostream&) The out stream.
 */
template <typename T, size_t n_rows, size_t n_cols, MatType type>
static inline std::ostream& operator<<(std::ostream& os, const Mat<T, n_rows, n_cols, type>& mat) {
    mat.print("", os);
    return os;
}

} // namespace flames

#ifndef FLAMES_PRESERVE_WARNING
#    ifdef __SYNTHESIS__
#        pragma GCC diagnostic pop
#    endif
#endif

#ifdef DEFINED_INLINE
#    define INLINE inline
#    undef DEFINED_INLINE
#endif

#endif
