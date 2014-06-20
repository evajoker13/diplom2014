#pragma once

#include "transform.hpp"
#include <array>
#include <cassert>

// find indexes of two most distant points
template <typename T = double, size_t N>
inline std::array<size_t, 2> diam(const std::vector<vec<T, N> > &points)
{
    std::array<size_t, 2> res = { 0, 0 };
    T longest                 = 0;
    for (size_t i = 0; i < points.size(); ++i)
    {
        const auto &a = points[i];
        for (size_t j = i + 1; j < points.size(); ++j)
        {
            const auto &b = points[j];

            auto dist = (b - a).length();
            if (dist > longest)
            {
                res     = { i, j };
                longest = dist;
            }
        }
    }
    return res;
}

template <typename T = double>
class ellipse
{
    // matrix that describes transformation into space where this ellipse looks
    // like a circle of radius 1.0
    tmatrix<T, 2> tm;

public:
    const decltype(tm) &as_tmatrix() const
    { return tm; }

    bool encloses(const vec<T,2> &point) const
    { return tm.apply(point).length() <= 1; }

    // approximate minimum volume enclosing ellipse
    // source: Journal of Mathematical Sciences, Vol. 97, No. 2, 1999
    //         "Pattern recognition with help of quadratic discriminant functions"
    static ellipse from_points(const std::vector<vec<T,2> > &points)
    {
        ellipse e;
        auto &tm = e.tm;

        // 1. find diameter
        auto d = diam(points);
        // order points by x coordinate
        if (points[d[0]][0] > points[d[1]][1])
            std::swap(d[0], d[1]);
        const auto &a = points[d[0]];
        const auto &b = points[d[1]];

        // 2. rotate to make diameter parallel to x
        // find angle for rotation
        auto ab    = b - a;
        auto angle = atan(ab[1] / ab[0]); // atan(y/x);
        tm = tm.rotate(angle);            // rotate
        // now we should work with points

        // 3. build a minimal volume enclosing rectangle
        T left, right;

        // points of diam forms left and right sides
        left  = tm.apply(a)[0];
        right = tm.apply(b)[0];

        // assume top and bottom for first point y' coordinate
        T top = tm.apply(points[0])[1], bottom = top;
        for (const auto &p : points)
        {
            auto y = tm.apply(p)[1];
            // adjust top and bottom to fit this point
            if (y > top)
                top = y;
            else if (y < bottom)
                bottom = y;
        }

        // 4. move everything to make center of rectangle be in (0,0)
        tm = tm.translate({ -(left + right)/2, -(bottom + top)/2 }).apply(tm);

        // 5. turn rectangle into square
        auto aspect_ratio = (right - left) / (top - bottom);
        tm = tm.scale({ 1, aspect_ratio }).apply(tm);

        // 6. find minimal radius of enclosing circle with center of square
        T radius = 0;
        for (const auto &p : points)
        {
            auto dist = tm.apply(p).length();
            if (dist > radius)
                radius = dist;
        }

        // 7. scale our space to make that circle of radius 1.0
        tm = tm.scale(1.0/radius).apply(tm);

        return e;
    }
};
