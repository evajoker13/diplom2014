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
}

int main()
{
    ifstream a_file("a.pok"), b_file("b.pok");
    group_type a_group, b_group;
    load_group(a_group, a_file);
    load_group(b_group, b_file);

    psamples_type a = toPSamples(a_group), b = toPSamples(b_group);
    classificator cf;
    cf.train(a, b);

    // check groups
    ifstream c_file("c.pok");
    group_type c_group;
    load_group(c_group, c_file);
    cout << cf.classify(toPSamples(c_group));
    return 0;
}
