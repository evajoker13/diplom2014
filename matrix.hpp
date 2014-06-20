#pragma once

#include <array>
#include <ostream>

template<typename T, size_t N, size_t M>
class matrix
{
    std::array<T, N*M> _elems;

public:
    constexpr matrix() {}

    static matrix from_elems(const std::array<T, N*M> &elems)
    {
        matrix m;
        m._elems = elems;
        return m;
    }

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

    constexpr const T *data() { return _elems.data(); }
    T *data() { return _elems.data(); }

    constexpr const T *begin() { return data(); }
    T *begin() { return data(); }

    constexpr const T *end() { return data() + N*M; }
    T *end() { return data() + N*M; }

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
        { return _matrix.data() + (_m*M); }
        constexpr iterator end()
        { return _matrix.data() + ((_m+1)*M); }

        row &operator++ () { ++_m; }
    };
};

template <typename T, size_t N, size_t M, size_t L>
inline matrix<T, N, M> operator*(const matrix<T, N, L> &a, const matrix<T, L, M> &b)
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
inline matrix<T, N, M> operator+(const matrix<T, N, M> &a, const matrix<T, N, M> &b)
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
inline matrix<T, N, M> operator-(const matrix<T, N, M> &a, const matrix<T, N, M> &b)
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

template <typename T, size_t N, size_t M>
inline std::ostream &operator<<(std::ostream &os, const matrix<T, N, M> &m)
{
    os << "[" << m(0,0);
    for (size_t j = 1; j < M; ++j)
        os << ", " << m(0,j);
    for (size_t i = 1; i < N; ++i)
    {
        os << ';' << std::endl << ' ' << m(i,0);
        for (size_t j = 1; j < M; ++j)
            os << ", " << m(i,j);
    }
    os << "]" << std::endl;
}
