#include "p_statistics.hpp"
#include "transform.hpp"
#include "samples.hpp"
#include "p_statistics.hpp"
#include "classificator.hpp"
#include <functional>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;

ostream &operator<<(ostream &os, const vector<classificator::result> &results)
{
    auto prec = os.precision(5);
    auto flags = os.flags();
    os.flags(flags | os.fixed);
    os << setw(12) << "CMG" << ' '
       << setw(12) << "FAM" << ' '
       << setw(12) << "TCMG" << ' '
       << setw(12) << "TFAM" << endl;

    for (const auto &result : results)
    {
        auto decision =
            (result.tcmg > result.tfam) ? "CMG" :
            (result.tcmg < result.tfam) ? "FAM" : "";
        os << setw(12) << result.cmg << ' '
           << setw(12) << result.fam << ' '
           << setw(12) << result.tcmg << ' '
           << setw(12) << result.tfam << "   "
           << decision
           << endl;
    }

    os.flags(flags);
    os.precision(prec);
    return os;
}

int main()
{
    ifstream a_file("a.pok"), b_file("b.pok");
    group_type a_group, b_group;
    load_group(a_group, a_file);
    load_group(b_group, b_file);

    // calibration

    for (auto group : {&a_group, &b_group})
    {
        bool is_cmg = group == &a_group;
        cout << "Calibration for " << (is_cmg ? "CMG" : "FAM") << " group" << endl;
        assert(group == &a_group || group == &b_group);
        vector<classificator::result> results;
        results.reserve(group->size());
        const size_t n = group->size();
        size_t k = 0;
        for (auto i = group->begin(); i != group->end(); ++i, ++k)
        {
            group_type sub_group;
            for (auto j = group->begin(); j != group->end(); ++j)
            {
                if (i == j) continue;
                sub_group.emplace_back(*j);
            }
            classificator cf;
            if (group == &a_group) cf.train(toPSamples(*group), toPSamples(b_group));
            else cf.train(toPSamples(a_group), toPSamples(*group));
            for (const auto &result : cf.classify(toPSamples({*i})))
            {
                results.emplace_back(result);
            }
            cout << k * 100 / n << "%\r    ";
            cout.flush();
        }
        cout << endl;
        cout << results << endl;
        size_t tcmg = 0, tfam = 0;
        for (auto &result : results)
        {
            if (result.tcmg > result.tfam) ++tcmg;
            if (result.tcmg < result.tfam) ++tfam;
        }
        size_t total = tcmg + tfam;
        if (is_cmg) cout << "correctness: " << tcmg * 100.0 / total << "%" << endl;
        else cout << "correctness: " << tfam * 100.0 / total << "%" << endl;
    }

    psamples_type a = toPSamples(a_group), b = toPSamples(b_group);
    classificator cf;
    cf.train(a, b);

    // check groups
    cout << "classification from c.pok" << endl;
    ifstream c_file("c.pok");
    group_type c_group;
    load_group(c_group, c_file);
    cout << cf.classify(toPSamples(c_group));
    return 0;
}
