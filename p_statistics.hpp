#pragma once

#include "samples.hpp"
#include <vector>
#include <array>

#define PSTAT_HISTOGRAM
#define PSTAT_BSEARCH

std::array<double,2> p_interval(double h, size_t m);

double p_statistics(const std::vector<double> &x, const std::vector<double> &y);

template <typename C>
double p_statistics(const std::vector<double> &x, const C &ys)
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
    distribution(std::vector<double> &&source) :
        samples(std::forward<std::vector<double> >(source))
    {}

    double distance(const distribution &other) const
    { return ::p_statistics(samples, other.samples); }
};

typedef std::array<distribution, std::tuple_size<sample_type>::value> properties_type;
typedef std::deque<properties_type> psamples_type;

psamples_type toPSamples(const group_type &group);
