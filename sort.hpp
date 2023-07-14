#ifndef _FLAMES_SORT_HPP_
#define _FLAMES_SORT_HPP_

#ifndef _FLAMES_CORE_HPP_
#    include "core.hpp"
#endif

namespace flames {

template <typename V>
static void mergeSort(V& vec) {
    constexpr size_t size = V::size();
    typename V::value_type temp[size];
MERGE_SORT_STAGE:
    for (int width = 1; width < size; width *= 2) {
        int f1 = 0;
        int f2 = width;
        int i2 = width;
        int i3 = 2 * width;
        if (i2 >= size) i2 = size;
        if (i3 >= size) i3 = size;
    MERGE_ARRAYS:
        for (int i = 0; i < size; ++i) {
            FLAMES_PRAGMA(PIPELINE II = 1)
            typename V::value_type t1 = vec[f1];
            typename V::value_type t2 = (f2 == i3) ? static_cast<typename V::value_type>(0) : vec[f2];
            if (f2 == i3 || (f1 < i2 && t1 <= t2)) {
                temp[i] = t1;
                ++f1;
            } else {
                assert(f2 < i3);
                temp[i] = t2;
                ++f2;
            }
            if (f1 == i2 && f2 == i3) {
                f1 = i3;
                i2 += 2 * width;
                i3 += 2 * width;
                if (i2 >= size) i2 = size;
                if (i3 >= size) i3 = size;
                f2 = i2;
            }
        }
    MERGE_SORT_COPY:
        for (int i = 0; i < size; ++i) {
            FLAMES_PRAGMA(UNROLL factor = FLAMES_MAT_COPY_UNROLL_FACTOR)
            vec[i] = temp[i];
        }
    }
}

template <typename V1, typename V2>
static void mergeSort(const V1& in, V2& out) {
    constexpr size_t size = V1::size();
    static_assert(size == V2::size(), "Sort in and out vectors/matrices should be of same size.");
    typename V1::value_type temp[size];
MERGE_SORT_STAGE:
    for (int width = 1; width < size; width *= 2) {
        int f1 = 0;
        int f2 = width;
        int i2 = width;
        int i3 = 2 * width;
        if (i2 >= size) i2 = size;
        if (i3 >= size) i3 = size;
    MERGE_ARRAYS:
        for (int i = 0; i < size; ++i) {
            FLAMES_PRAGMA(PIPELINE II = 1)
            typename V1::value_type t1 = in[f1];
            // TODO: here minimum is zero
            typename V1::value_type t2 = (f2 == i3) ? static_cast<typename V1::value_type>(0) : in[f2];
            if (f2 == i3 || (f1 < i2 && t1 <= t2)) {
                out[i] = t1;
                ++f1;
            } else {
                assert(f2 < i3);
                out[i] = t2;
                ++f2;
            }
            if (f1 == i2 && f2 == i3) {
                f1 = i3;
                i2 += 2 * width;
                i3 += 2 * width;
                if (i2 >= size) i2 = size;
                if (i3 >= size) i3 = size;
                f2 = i2;
            }
        }
    }
}

template <typename V>
static inline void sort(V& vec) {
    return mergeSort(vec);
}

template <typename V1, typename V2>
static inline void sort(const V1& in, V2& out) {
    return mergeSort(in, out);
}

template <typename T, typename I>
static inline void argmax_4_2(T in1, T in2, T in3, T in4, I i_in1, I i_in2, I i_in3, I i_in4, T& out1, T& out2,
                              I& i_out1, I& i_out2, bool sorted = false) {
    FLAMES_PRAGMA(INLINE)
    if (sorted) {
        assert(in1 >= in2 && in3 >= in4 && "Should be sorted in flames::argmax_4_2 with sorted=true.");
        if (in1 < in4) { // 3, 4
            out1   = in3;
            out2   = in4;
            i_out1 = i_in3;
            i_out2 = i_in4;
        } else if (in3 < in2) { // 1, 2
            out1   = in1;
            out2   = in2;
            i_out1 = i_in1;
            i_out2 = i_in2;
        } else if (in1 < in3) { // 3, 1
            out1   = in3;
            out2   = in1;
            i_out1 = i_in3;
            i_out2 = i_in1;
        } else { // 1, 3
            out1   = in1;
            out2   = in3;
            i_out1 = i_in1;
            i_out2 = i_in3;
        }
    } else {
        T _in1, _in2, _in3, _in4;
        I _i_in1, _i_in2, _i_in3, _i_in4;
        // pre sort 1 and 2
        if (in1 > in2) {
            _in1   = in1;
            _in2   = in2;
            _i_in1 = i_in1;
            _i_in2 = i_in2;
        } else {
            _in1   = in2;
            _in2   = in1;
            _i_in1 = i_in2;
            _i_in2 = i_in1;
        }
        if (in3 > in4) {
            _in3   = in3;
            _in4   = in4;
            _i_in3 = i_in3;
            _i_in4 = i_in4;
        } else {
            _in3   = in4;
            _in4   = in3;
            _i_in3 = i_in4;
            _i_in4 = i_in3;
        }
        // now same as sorted=true
        if (_in1 < _in4) { // 3, 4
            out1   = _in3;
            out2   = _in4;
            i_out1 = _i_in3;
            i_out2 = _i_in4;
        } else if (_in3 < _in2) { // 1, 2
            out1   = _in1;
            out2   = _in2;
            i_out1 = _i_in1;
            i_out2 = _i_in2;
        } else if (_in1 < _in3) { // 3, 1
            out1   = _in3;
            out2   = _in1;
            i_out1 = _i_in3;
            i_out2 = _i_in1;
        } else { // 1, 3
            out1   = _in1;
            out2   = _in3;
            i_out1 = _i_in1;
            i_out2 = _i_in3;
        }
    }
}

} // namespace flames

#endif
