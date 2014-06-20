#include "classificator.hpp"
#include "transform.hpp"
#include <iostream>
#include <vector>
#include <array>

using namespace std;

template <typename I1, typename I2>
std::vector<std::pair<double, double>> zip(I1 first1, I1 last1, I2 first2)
{
    std::vector<std::pair<double, double>> result;
    result.reserve(last1 - first1);
    auto i1 = first1;
    auto i2 = first2;
    for (; i1 != last1; ++i1, ++i2) result.emplace_back(*i1, *i2);
    return result;
}

void classificator::train(const psamples_type &a, const psamples_type &b)
{
    const size_t n = a.size(), m = b.size();

    array<vector<double>, chars_count>  aa, ab, ba, bb;

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

    // couple together characteristics
    for (size_t t = 0; t < chars_count; ++t)
    {
        for (size_t s = 0; s < chars_count; ++s)
        {
            // TODO: ellipsis aa, ab, ba, bb
            auto points = zip(begin(aa[t]), end(aa[t]), begin(aa[s]));
            auto tf = transform<double>::translate({1, 1});
        }
    }
}
