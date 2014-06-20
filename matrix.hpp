#pragma once

#include <array>

template<typename T, size_t N, size_t M>
class matrix
{
    std::array<T, N*M> _elems;

public:
    constexpr matrix() {}

    constexpr size_t cols() const
    { return M; }

    constexpr size_t rows() const
    { return N; }

    T &operator() (size_t i, size_t j)
    {
        return _elems.at(i*M + j);
    }

    const T &operator() (size_t i, size_t j) const
    { return _elems.at(i*M + j); }

    void fill(const T &x)
    { _elems.fill(x); }

    class row
    {
        matrix &_matrix;
        size_t _m;

    public:
        constexpr row(matrix &origin, size_t m) :
            _matrix(origin), _m(m)
        {}

        typedef T* iterator;

        constexpr iterator begin()
        { return _matrix.data() + (_m*_matrix.cols()); }
        constexpr iterator end()
        { return _matrix.data() + ((_m+1)*_matrix.cols()); }

        row &operator++ () { ++_m; }
    };
};

template <typename T, size_t N, size_t M, size_t L>
inline matrix<T, N, M> operator*(const matrix<T, N, L> &a, matrix<T, L, M> &b)
{
    // a.cols() == b.rows()
    // c.rows() == a.rows()
    // c.cols() == b.cols()
    matrix<T, N, M> c;
    for (size_t i = 0; i < c.rows(); ++i)
    {
        for (size_t j = 0; j < c.cols(); ++j)
        {
            T product = 0;
            for (size_t k = 0; k < L; ++k)
            {
                product += a(i, k) * b(k, j);
            }
            c(i, j) = product;
        }
    }
    return c;
}

template <typename T, size_t N, size_t M>
inline matrix<T, N, M> operator+(const matrix<T, N, M> &a, matrix<T, N, M> &b)
{
    matrix<T, N, M> c;
    for (size_t i = 0; i < c.rows(); ++i)
    {
        for (size_t j = 0; j < c.cols(); ++j)
        {
            c(i, j) = a(i, j) + b(i, j);
        }
    }
    return c;
}

template <typename T, size_t N, size_t M>
inline matrix<T, N, M> operator-(const matrix<T, N, M> &a, matrix<T, N, M> &b)
{
    matrix<T, N, M> c;
    for (size_t i = 0; i < c.rows(); ++i)
    {
        for (size_t j = 0; j < c.cols(); ++j)
        {
            c(i, j) = a(i, j) - b(i, j);
        }
    }
    return c;
}
