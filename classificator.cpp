#include "classificator.hpp"
#include "vec.hpp"
#include "ellipse.hpp"
//#include <iostream>
#include <vector>
#include <array>

using namespace std;

namespace {

constexpr size_t chars_count = std::tuple_size<properties_type>::value;

template <typename I1, typename I2>
std::vector<vec<double, 2>> zip(I1 first1, I1 last1, I2 first2)
{
    std::vector<vec<double, 2>> result;
    result.reserve(last1 - first1);
    auto i1 = first1;
    auto i2 = first2;
    for (; i1 != last1; ++i1, ++i2) result.emplace_back(*i1, *i2);
    return result;
}

// 2nd stage - distances with help of Penunin's statistics
struct distances
{
    array<vector<double>, chars_count>  aa, ab, ba, bb;

    distances(const psamples_type &a, const psamples_type &b)
    {
        const size_t n = a.size(), m = b.size();

        // build distance per characterestic
        for (size_t k = 0; k < chars_count; ++k)
        {
            aa[k].reserve(n);
            ab[k].reserve(n);
            for (auto i = a.begin(); i != a.end(); ++i)
            {
                double p = 0;
                for (auto j = a.begin(); j != a.end(); ++j)
                {
                    if (i == j)
                        continue;
                    p += (*i)[k].distance((*j)[k]);
                }
                aa[k].emplace_back(p / (n - 1));

                double q = 0;
                for (auto j = b.begin(); j != b.end(); ++j)
                {
                    q += (*i)[k].distance((*j)[k]);
                }
                ab[k].emplace_back(q / m);
                //cout << "k = " << k <<  ", aa = " << aa[k].back() << ", ab = " << ab[k].back() << endl;
            }

            ba[k].reserve(m);
            bb[k].reserve(m);
            for (auto i = a.begin(); i != a.end(); ++i)
            {
                double p = 0;
                for (auto j = a.begin(); j != a.end(); ++j)
                {
                    p += (*i)[k].distance((*j)[k]);
                }
                ba[k].emplace_back(p / n);

                double q = 0;
                for (auto j = b.begin(); j != b.end(); ++j)
                {
                    if (i == j)
                        continue;
                    q += (*i)[k].distance((*j)[k]);
                }
                bb[k].emplace_back(q / (m - 1));
                //cout << "k = " << k <<  ", ba = " << ba[k].back() << ", bb = " << bb[k].back() << endl;
            }
        }
    }
};

// 3rd stage - ellipses
struct ellipses
{
    matrix<ellipse<double>, chars_count, chars_count> aa, ab, ba, bb;
    ellipses(const distances &dists)
    {
        // couple together characteristics and build ellipsis out of those vectors
        for (size_t t = 0; t < chars_count; ++t)
        {
            for (size_t s = 0; s < chars_count; ++s)
            {
                aa(t,s) = ellipse<double>::from_points(zip(begin(dists.aa[t]), end(dists.aa[t]), begin(dists.aa[s])));
                ab(t,s) = ellipse<double>::from_points(zip(begin(dists.ab[t]), end(dists.ab[t]), begin(dists.ab[s])));
                ba(t,s) = ellipse<double>::from_points(zip(begin(dists.ba[t]), end(dists.ba[t]), begin(dists.ba[s])));
                bb(t,s) = ellipse<double>::from_points(zip(begin(dists.bb[t]), end(dists.bb[t]), begin(dists.bb[s])));
            }
        }
    }
};

}

void classificator::train(const psamples_type &a, const psamples_type &b)
{
    const size_t n = a.size(), m = b.size();

    distances stage2(a, b);

    ellipses stage3(stage2);
}
