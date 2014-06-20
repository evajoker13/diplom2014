#pragma once

#include <vector>
#include "samples.hpp"

using namespace std;

#define PSTAT_HISTOGRAM
#define PSTAT_BSEARCH

double p_statistics(const vector<double> &x, const vector<double> &y);

template <typename C>
double p_statistics(const vector<double> &x, const C &ys)
{
    double sum = 0;
    for (const auto &y : ys) sum += p_statistics(x, y);
    return sum / ys.size();
}

class distribution
{
    std::vector<double> samples;
public:
    distribution() {}
    distribution(vector<double> &&source) :
        samples(std::forward<vector<double> >(source))
    {}

    double distance(const distribution &other) const
    { return ::p_statistics(samples, other.samples); }
};

typedef std::array<distribution, tuple_size<sample_type>::value> properties_type;
typedef std::deque<properties_type> psamples_type;

psamples_type toPSamples(const group_type &group);
