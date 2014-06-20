#pragma once

#include "p_statistics.hpp"

class classificator
{
    static constexpr size_t chars_count = std::tuple_size<properties_type>::value;
public:
    void train(const psamples_type &a, const psamples_type &b);
};
