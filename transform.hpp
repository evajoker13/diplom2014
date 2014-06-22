#pragma once

#include "vec.hpp"
#include "matrix.hpp"
#include <array>
#include <algorithm>

template <typename T, size_t N = 2>
class tmatrix
{
    // matrix for left multiplication (direct and inverse)
    matrix<T, N+1, N+1> m, im;

public:
    constexpr const decltype(m) &as_matrix() { return m; }
    static tmatrix scale(const std::array<T, N> &factors)
    {
        tmatrix tf;
        tf.m.fill(0);
        for (size_t i = 0; i < N; ++i)
        {
            tf.m(i, i)  = factors[i];
            tf.im(i, i) = 1/factors[i];
        }
        tf.m(N,N) = 1; // augement
        return tf;
    }

    static tmatrix scale(T factor)
    {
        tmatrix tf;
        tf.m.fill(0);
        auto ifactor = 1/factor;
        for (size_t i = 0; i < N; ++i)
        {
            tf.m(i, i)  = factor;
            tf.im(i, i) = ifactor;
        }
        tf.m(N,N) = 1; // augement
        return tf;
    }

    static const tmatrix &identity()
    {
        static tmatrix tf = scale(1);
        return tf;
    }

    static tmatrix translate(const std::array<T, N> &offset)
    {
        tmatrix tf = identity();
        for (size_t i = 0; i < N; ++i)
        {
            tf.m(i, N)  = offset[i];
            tf.im(i, N) = -offset[i];
        }
        return tf;
    }

    static tmatrix rotate(T angle, size_t i, size_t j)
    {
        tmatrix tf = identity();
        tf.m(i, i) = std::cos(angle);
        tf.m(i, j) = std::sin(angle);
        tf.m(j, i) = -std::sin(angle);
        tf.m(j, j) = std::cos(angle);

        tf.im(i, i) = std::cos(-angle);
        tf.im(i, j) = std::sin(-angle);
        tf.im(j, i) = -std::sin(-angle);
        tf.im(j, j) = std::cos(-angle);
        return tf;
    }

    static tmatrix rotate(T angle)
    { return rotate(angle, 0, 1); }

    matrix<T, N, 1> apply(const matrix<T, N, 1>  &vec, bool direct = true)
    {
        using namespace std;
        auto &tm = direct ? m : im;
        matrix<T, N+1, 1> avec;
        copy(begin(vec), end(vec), begin(avec));
        avec(N, 0) = 1;
        avec = tm * avec;
        matrix<T, N, 1> res;
        copy(begin(avec), end(avec)-1, begin(res));
        return res;
    }

    vec<T, N> apply(const vec<T, N>  &v, bool direct = true) const
    {
        using namespace std;
        auto &tm = direct ? m : im;
        matrix<T, N+1, 1> avec;
        copy(begin(v), std::end(v), std::begin(avec));
        avec(N, 0) = 1;
        avec = tm * avec;
        vec<T, N> res;
        copy(begin(avec), end(avec)-1, begin(res));
        return res;
    }

    tmatrix apply(const tmatrix &tm, bool direct = true) const
    {
        tmatrix res;
        if (direct)
        {
            res.m = m * tm.m;
            res.im = im * tm.im;
        }
        else
        {
            res.m = m * tm.im;
            res.im = im * tm.m;
        }
        return res;
    }
};
