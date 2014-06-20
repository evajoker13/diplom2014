#include "p_statistics.hpp"
#include "transform.hpp"
#include "samples.hpp"
#include "p_statistics.hpp"
#include "classificator.hpp"
#include <functional>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert>

using namespace std;

int main()
{
    matrix<double, 2, 1> vec;
    vec(0,0) = 1;
    vec(1,0) = 2;
    cout << "[" << vec(0,0) << ", " << vec(1,0) << "]" << endl;
    tmatrix<double, 2> tf = tf.translate({ -1, -3 });
    //tmatrix<double, 2> tf = tf.identity();
    vec = tf.apply(vec);
    cout << "[" << vec(0,0) << ", " << vec(1,0) << "]" << endl;

    ifstream a_file("a.pok"), b_file("b.pok");
    group_type a_group, b_group;
    load_group(a_group, a_file);
    load_group(b_group, b_file);

    // check groups
    group_type c_group, d_group;
    c_group.emplace_back(std::move(a_group.back()));
    a_group.pop_back();
    d_group.emplace_back(std::move(b_group.back()));
    b_group.pop_back();

    psamples_type a = toPSamples(a_group), b = toPSamples(b_group);
    classificator cf;
    cf.train(a, b);
    return 0;

    for (size_t i = 0; i < tuple_size<properties_type>::value; ++i)
    {
        cout << "Property #" << i << endl;
        for (size_t j = 0; j < a.size(); ++j)
        {
            double p = 0;
            for (size_t k = 0; k < a.size(); ++k)
            {
                if (j == k )
                    continue; // skip same
                p += a[j][i].distance(a[k][i]);
            }
            p /= a.size() - 1;
            double np = 0;
            for (size_t k = 0; k < b.size(); ++k)
            {
                np += a[j][i].distance(b[k][i]);
            }
            np /= b.size();
            cout << j << ": p = " << p << ", np = " << np << endl;;
        }
    }

    return 0;
}

