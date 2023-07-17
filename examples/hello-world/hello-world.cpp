#include "flames/flames.hpp"

using V = Vec<int, 2>;

/**
 * @brief A hello world top function returning the sum of two vectors.
 *
 * @param a Input vector 1.
 * @param b Input vector 2.
 * @return (V) The sum vector.
 */
V helloWorld(const V& a, const V& b) { return a + b; }

int main() {
#ifndef __SYNTHESIS__
    std::cout << "GCC Version: " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << std::endl;
    // GCC Version: 6.2.0
#endif
    V a = std::vector<int>{ 2, 3 }, b = std::vector<int>{ 5, 6 };
    auto sum = helloWorld(a, b);
    a.print("a = ");
    b.print("b = ");
    sum.print("a + b = ");
}
