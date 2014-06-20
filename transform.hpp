#pragma once

#include "matrix.hpp"

template <typename T, size_t N = 2>
class transform
{
    // matrix for left multiplication (direct and inverse)
    matrix<T, N+1, N+1> m, im;
    constexpr transform() {}
public:
    static transform identity()
    {
        transform tf;
        tf.m.fill(0);
        for (size_t i = 0; i < N+1; ++i) tf.m(i,i) = 1;
        tf.im = tf.m;
        return tf;
    }
    static transform translate(const array<T, N> &offset)
    {
        transform tf = identity();
        for (size_t i = 0; i < N; ++i)
        {
            tf.m(i, N) = offset[i];
            tf.im(i, N) = -offset[i];
        }
        return tf;
    }

    matrix<T, N, 1> apply(const matrix<T, N, 1>  &vec)
    {
        matrix<T, N+1, 1> avec;
        for (size_t i = 0; i < N; ++i) avec(i, 0) = vec(i, 0);
        avec(N, 0) = 1;
        auto ares = m * avec;
        matrix<T, N, 1> res;
        for (size_t i = 0; i < N; ++i) res(i, 0) = ares(i, 0);
        return res;
    }
};
