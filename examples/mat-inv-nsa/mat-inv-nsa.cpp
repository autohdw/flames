#include "flames/flames.hpp"

using dtype = FxP<8, 8>;
using M     = Mat<dtype, 4, 4>;

M top(const M& A) { return A.invNSA(); }

int main() {
    M A{ 10, -2, 1, 0, 1, -8, 2, 0, 0, 0, 11, -1, 0, 1, 2, 4 };
    A.print("A = ");
    M A_inv = top(A);
    A_inv.print("A_inv = ");
    return 0;
}
