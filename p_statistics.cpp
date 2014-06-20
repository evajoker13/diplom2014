#include "p_statistics.hpp"
#include <algorithm>
#include <cassert>

double p_statistics(const vector<double> &x, const vector<double> &y)
{
    assert(std::is_sorted(x.begin(), x.end()));
    assert(std::is_sorted(y.begin(), y.end()));

    static constexpr double G  = 3;
    static constexpr double G2 = G * G;

    const size_t n = x.size();
    const size_t N = n * (n - 1) / 2;
    size_t L       = 0;

#ifdef PSTAT_HISTOGRAM
    // build histogram using x as a baseline
    size_t hist[n+1];
    {
        auto it = y.begin();

        size_t h = 0; // accumulator for amount of Y's on the left from X[i]
        size_t i = 0;
        for (; i < n && it != y.end(); )
        {
            if (*it < x[i])
            {
                auto it1 = std::upper_bound(it + 1, y.end(), x[i]);

                h      += it1 - it;
                hist[i] = h;

                it = it1;
                /*
                   hist[i] = ++h;
                   ++it;
                 */
            }
            else
            {
                // simply fill with same level
                hist[i++] = h;
            }
        }
        h += y.end() - it; // amount Ys left after cycle (ex. after last X[i])
        // fill rest of histogram with this level
        for (; i < (n+1); ++i) hist[i] = h;
    }
#endif

    for (size_t j = 1; j < n; ++j)
    {
        for (size_t i = 0; i < j; ++i)
        {
            // frequency h_{i,j}
#if defined(PSTAT_HISTOGRAM)
            double h = double(hist[j] - hist[i]) / n;

#elif defined(PSTAT_BSEARCH)
            auto lb  = std::lower_bound(y.begin(), y.end(), x[i]);
            auto ub  = std::upper_bound(lb, y.end(), x[j]);
            double h = double(ub - lb) / n;
#else

            size_t num = 0;
            for (size_t k = 0; k < n; ++k)
            {
                if (x[i] <= y[k] && y[k] <= x[j])
                    ++num;
            }
            double h = double(num) / n;
#endif

            // probability A_{i,j}
            double pA = double(j - i) / (n + 1);

            // conf. interval
            double middle = n * h + G2 / 2;
            double spread = G * sqrt(h * (1 - h) * n + G2 / 4);
            double norm   = n + G2;
            double p1     = (middle - spread) / norm,
                   p2     = (middle + spread) / norm;

            // check interval
            if (p1 <= pA && pA <= p2)
            {
                ++L;
            }
        }
    }

    return double(L) / N;
}

psamples_type toPSamples(const group_type &group)
{
    psamples_type psamples;
    for (auto person : group)
    {
        psamples.emplace_back();
        auto &props = psamples.back();
        for (size_t i = 0; i < props.size(); ++i)
        {
            vector<double> p;
            for (auto values : person)
            { p.push_back(values[i]); }
            sort(p.begin(), p.end());
            props[i] = std::move(p);
        }
    }
    return psamples;
}
