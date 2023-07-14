#include "flames/flames.hpp"

using dtype = FxP<6, 2>;
using M     = Mat<dtype, 4, 4>;
using V     = Vec<dtype, 4>;

void main_task(const M& A, const V& b, V& c) {
#pragma HLS INLINE
    M tmp1;
    V tmp2;
    tmp1 = A * A;
    tmp2 = tmp1 * b;
    c    = tmp2 % tmp2;
}

void top(const M& A1, const M& A2, const M& A3, const V& b, V& c) {
    V c1, c2;
    M tmp;
#pragma HLS PIPELINE
    main_task(A1, b, c1);
    main_task(A2, c1, c2);
    main_task(A3, c2, c);
}

int main() {
    M A1, A2, A3;
    V b, c;
    top(A1, A2, A3, b, c);
    return 0;
}
