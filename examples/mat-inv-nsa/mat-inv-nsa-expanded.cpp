#include "flames/flames.hpp"

using dtype = FxP<8, 8>;
using M     = Mat<dtype, 4, 4>;

M top(const M& A) {
    M A_inv;
    const auto D = A.diagMat_(); // diagonal part
    const auto E = A.offDiag_(); // off-diagonal part
    Mat<dtype, 4, 4, MatType::DIAGONAL> D_inv;
    D_inv.invDiag(D); // inverse of diagonal part
    Mat<dtype, 4, 4, MatType::NORMAL> product = (-D_inv) * E;
    Mat<dtype, 4, 4, MatType::NORMAL> sum_tmp = A_inv = product; // the first iteration
    Mat<dtype, 4, 4, MatType::NORMAL> tmp;
    const size_t iter = 4;
MAT_INV_NSA:
    for (size_t i = 1; i < iter; ++i) {
        tmp.mul(A_inv, product);
        A_inv = tmp;
        sum_tmp += tmp;
    }
    A_inv.mul(sum_tmp, D_inv);
    return A_inv += D_inv;
}

int main() {
    M A{ 10, -2, 1, 0, 1, -8, 2, 0, 0, 0, 11, -1, 0, 1, 2, 4 };
    A.print("A = ");
    M A_inv = top(A);
    A_inv.print("A_inv = ");
    return 0;
}
