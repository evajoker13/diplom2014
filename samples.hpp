#pragma once

#include <vector>
#include <deque>
#include <istream>
#include <array>

typedef std::array<double,15> sample_type;
typedef std::vector<sample_type> person_type;
typedef std::deque<person_type> group_type;

bool load_person(person_type &person, std::istream &f);

void load_group(group_type &group, std::istream &f);
