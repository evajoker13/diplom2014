#include "classificator.hpp"
#include "vec.hpp"
#include "ellipse.hpp"
#include <vector>
#include <array>
#include <iostream>

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

struct patient_distance
{
    std::array<double, chars_count> qa, qb;
    patient_distance(const properties_type &q, const psamples_type &a, const psamples_type &b)
    {
        const size_t n = a.size(), m = b.size();
        for (size_t k = 0; k < chars_count; ++k)
        {
            double p = 0, d = 0;
            auto qk = q[k];
            for (auto j = a.begin(); j != a.end(); ++j)
            {
                p += qk.distance((*j)[k]);
            }
            for (auto j = b.begin(); j != b.end(); ++j)
            {
                d += qk.distance((*j)[k]);
            }
            qa[k] = p / n;
            qb[k] = d / m;
            // cout << qa[k] << ", " << qb[k] << endl;
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

    struct amatch {
        // qa
        bool a1 : 1; // aa
        bool a2 : 1; // ab
        bool a3 : 1; // aa - ab
        bool a4 : 1; // ab - aa
        // qb
        bool a1_ : 1; // bb
        bool a2_ : 1; // ba
        bool a3_ : 1; // bb - ba
        bool a4_ : 1; // ba - bb

        /*
        operator cmatch() const
        {
            return {
                a3 || a4,
                a4 || a3_,
                a1 || a2_,
                a2 || a1_
            };
        }
        */
    };

    struct cmatch {
        /*
        bool c1 : 1; // cmg = a3 || a4
        bool c2 : 1; // fam = a4 || a3*
        bool c3 : 1; // tcmg = a1 || a2*
        bool c4 : 1; // tfam = a2 || a1*
        */

        std::array<bool, 6> c;
        cmatch(const amatch &a)
        {
            c = {
                a.a3 || a.a4,
                a.a4 || a.a3_,
                a.a1 || a.a2_,
                a.a2 || a.a1_,
            };
        }
    };

    struct ematch {
        bool aa : 1; // qa
        bool ab : 1;
        bool ba : 1; // qb
        bool bb : 1;

        operator amatch() const
        {
            return {
                aa,
                ab,
                aa && !ab,
                ab && !aa,
                bb,
                bb && !ba,
                ba && !bb
            };
        }
    };

    ematch check(const patient_distance &pd, size_t t, size_t s) const
    {
        vec<double,2> va = { pd.qa[t], pd.qa[s] };
        vec<double,2> vb = { pd.qb[t], pd.qb[s] };
        return {
            aa(t,s).encloses(va),
            ab(t,s).encloses(va),
            ba(t,s).encloses(vb),
            bb(t,s).encloses(vb),
        };
    }

    classificator::result freq(const patient_distance &pd) const
    {
        classificator::result res = { 0, 0, 0, 0 };
        for (size_t t = 0; t < chars_count; ++t)
        {
            for (size_t s = t+1; s < chars_count; ++s)
            {
                cmatch c = amatch(check(pd, t, s));
                if (c.c[0]) ++res.cmg;
                if (c.c[1]) ++res.fam;
                if (c.c[2]) ++res.tcmg;
                if (c.c[3]) ++res.tfam;
            }
        }
        constexpr size_t m = chars_count * chars_count;
        res.cmg /= m;
        res.fam /= m;
        res.tcmg /= m;
        res.tfam /= m;
        return res;
    }

};

}

void classificator::train(const psamples_type &a, const psamples_type &b)
{
    const size_t n = a.size(), m = b.size();

    distances stage2(a, b);

    ellipses stage3(stage2);

    classifier = [=](const properties_type &q) -> result {
        return stage3.freq(patient_distance(q, a, b));
    };
}
