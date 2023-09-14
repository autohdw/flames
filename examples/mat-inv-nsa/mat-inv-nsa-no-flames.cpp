#include "ap_fixed.h"
#include <iostream>
#include <string>

using dtype = ap_fixed<17, 8>;

void print(dtype A[4][4], std::string s = "") {
#ifndef __SYNTHESIS__
    std::cout << s << "[";
    for (size_t i = 0; i + 1 < 4; ++i) {
        std::cout << "[";
        for (size_t j = 0; j + 1 < 4; ++j) std::cout << A[i][j] << ", ";
        std::cout << A[i][3] << "]," << std::endl;
    }
    std::cout << "[";
    for (size_t j = 0; j + 1 < 4; ++j) std::cout << A[3][j] << ", ";
    std::cout << A[3][3] << "]]" << std::endl;
#endif
}

void mat_copy(dtype from[4][4], dtype to[4][4]) {
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable = from complete
#pragma HLS ARRAY_PARTITION variable = to complete
    for (size_t i = 0; i != 4; ++i) {
#pragma HLS UNROLL
        for (size_t j = 0; j != 4; ++j) {
#pragma HLS LOOP_FLATTEN
            to[i][j] = from[i][j];
        }
    }
}

void top(dtype A[4][4], dtype A_inv[4][4]) {
#pragma HLS ARRAY_PARTITION variable = A complete
#pragma HLS ARRAY_PARTITION variable = A_inv complete
    dtype product[4][4];
    dtype D_inv[4];
    for (size_t i = 0; i != 4; ++i) {
#pragma HLS UNROLL
        D_inv[i] = static_cast<dtype>(1) / A[i][i];
    }
MAT_DIAG_TIMES_MAT_NORMAL:
    for (size_t j = 0; j != 4; ++j) {
#pragma HLS UNROLL
        for (size_t i = 0; i != 4; ++i) {
#pragma HLS LOOP_FLATTEN
            product[i][j] = i == j ? static_cast<dtype>(0) : static_cast<dtype>(-D_inv[i] * A[i][j]);
        }
    }
    dtype sum_tmp[4][4], tmp[4][4];
    mat_copy(product, sum_tmp);
    mat_copy(product, A_inv);
    const size_t iter = 4;
    for (size_t i = 1; i < iter; ++i) {
    // tmp = A_inv * product;
    GEMM:
        for (size_t _i = 0; _i != 4; ++_i) {
        GEMM_r:
            for (size_t r = 0; r != 4; ++r) {
#pragma HLS UNROLL
            GEMM_c:
                for (size_t c = 0; c != 4; ++c) {
#pragma HLS LOOP_FLATTEN
                    if (_i == 0) tmp[r][c] = static_cast<dtype>(0);
                    tmp[r][c] += A_inv[r][_i] * product[_i][c];
                }
            }
        }
        mat_copy(tmp, A_inv);
        // sum_tmp += tmp;
        for (size_t _i = 0; _i != 4; ++_i) {
#pragma HLS UNROLL
            for (size_t j = 0; j != 4; ++j) {
#pragma HLS LOOP_FLATTEN
                sum_tmp[_i][j] += tmp[_i][j];
            }
        }
    }
// A_inv = sum_tmp * D_inv;
MAT_NORMAL_TIMES_MAT_DIAG:
    for (size_t i = 0; i != 4; ++i) {
#pragma HLS UNROLL
        for (size_t j = 0; j != 4; ++j) {
#pragma HLS LOOP_FLATTEN
            A_inv[i][j] = sum_tmp[i][j] * D_inv[j];
        }
    }
    // A_inv += D_inv;
    for (size_t i = 0; i != 4; ++i) {
#pragma HLS UNROLL
        A_inv[i][i] += D_inv[i];
    }
}

int main() {
    dtype A[4][4] = { { 10, -2, 1, 0 }, { 1, -8, 2, 0 }, { 0, 0, 11, -1 }, { 0, 1, 2, 4 } };
    print(A, "A = ");
    dtype A_inv[4][4];
    top(A, A_inv);
    print(A_inv, "A_inv = ");
    return 0;
}
