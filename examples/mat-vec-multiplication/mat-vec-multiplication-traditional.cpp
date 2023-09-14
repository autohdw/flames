#include <ap_int.h>
#include <cstddef>

using dtype = ap_int<8>;

void top(dtype A[4][4], dtype b[4], dtype c[4]) {
#pragma HLS ARRAY_PARTITION variable = A complete
#pragma HLS ARRAY_PARTITION variable = b complete
    for (size_t i = 0; i != 4; ++i) {
        for (size_t j = 0; j != 4; ++j) {
#pragma HLS UNROLL
            c[i] = A[i][j] * b[j];
        }
    }
}

int main() {
    dtype A[4][4] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    dtype b[4]    = { 0, 1, 2, 3 }, c[4];
    top(A, b, c);
    return 0;
}
