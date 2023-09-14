#include "flames/flames.hpp"

using dtype = ap_int<8>;
using V     = Vec<dtype, 4>;
using M     = Mat<dtype, 4, 4>;

void top(M& A, const V& b) {
    for (size_t i = 0; i != 4; ++i) A.col_(i) = b;
    // // The following method is depracted and is not advised to use:
    // A.col(i, b);
}

int main() {
    M A{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    V b{ 0, 1, 2, 3 };
    A.print("A before = ");
    top(A, b);
    A.print("A after = ");
    return 0;
}
