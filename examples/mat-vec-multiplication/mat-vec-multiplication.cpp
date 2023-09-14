#include "flames/flames.hpp"

using dtype = ap_int<8>;
using V     = Vec<dtype, 4>;
using M     = Mat<dtype, 4, 4>;

void top(const M& A, const V& b, V& c) { c = A * b; }

int main() {
    M A{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    V b{ 0, 1, 2, 3 }, c;
    top(A, b, c);
}
