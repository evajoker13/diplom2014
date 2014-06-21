#pragma once

#include "p_statistics.hpp"
#include <functional>

struct classificator
{
    struct result {
        float cmg; // h_1
        float fam;
        float tcmg;
        float tfam; // h_4
    };
    void train(const psamples_type &a, const psamples_type &b);
    result classify(const properties_type &q) const
    { return classifier(q); }

    std::vector<result> classify(const psamples_type &qs) const
    {
        std::vector<result> res;
        res.reserve(qs.size());
        for (const auto &q : qs) res.emplace_back(classify(q));
        return res;
    }

private:
    std::function<result(const properties_type &)> classifier;
};
