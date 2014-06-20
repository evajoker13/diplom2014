#pragma once

#include "matrix.hpp"
#include <cmath>
#include <ostream>

template <typename T, size_t N = 2>
struct vec : std::array<T, N>
{
    using std::array<T, N>::array;

    template <typename ... Args>
    vec(Args &&... args) :
        std::array<T,N>({ T(std::forward<Args>(args)) ... }
                        ) {}

    T length() const
    {
        T sum = 0;
        for (auto x : *this) sum += x*x;
        return std::sqrt(sum);
    }

    vec scale(T factor) const
    {
        vec res;
        for (size_t i = 0; i < N; ++i) res[i] = (*this)[i] * factor;
        return res;
    }

    vec operator-(const vec &v) const
    {
        vec res;
        for (size_t i = 0; i < N; ++i) res[i] = (*this)[i] - v[i];
        return res;
    }

    vec operator+(const vec &v) const
    {
        vec res;
        for (size_t i = 0; i < N; ++i) res[i] = (*this)[i] + v[i];
        return res;
    }

    explicit operator matrix<T, N, 1>() const { return *this; }
    explicit operator matrix<T, 1, N>() const { return *this; }
};

template <typename T, size_t N>
inline std::ostream &operator<<(std::ostream &os, const vec<T, N> &v)
{
    auto it = std::begin(v);
    os << '[' << *it;
    for (++it; it != std::end(v); ++it)
    {
        os << ", " << *it;
    }
    os << ']';
}
