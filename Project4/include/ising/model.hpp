// include/ising/model.hpp
/*
For model parameters like E-field etc.
*/

#pragma once
#include <cmath>

namespace ising{
    struct Model{
        double J;
        bool double_count;

        constexpr Model(double J = 1.0, bool double_count = false)
        : J(J),
        double_count(double_count) {}
    };
}