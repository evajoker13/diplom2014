#include "samples.hpp"
#include <cassert>

bool load_person(person_type &person, std::istream &f)
{
    assert(f.good());
    size_t n;
    f >> n;
    if (!f.good())
        return false;
    person.resize(n);
    for (auto &sample : person)
    {
        size_t i;
        f >> i;
        assert(f.good());
        assert(i <= n);
        for (auto &value : sample)
        {
            f >> value;
        }
    }
    return true;
}

void load_group(group_type &group, std::istream &f)
{
    assert(f.good());
    for(;; )
    {
        group.emplace_back();
        if (!load_person(group.back(), f))
        {
            group.pop_back();
            break;
        }
    }
}
