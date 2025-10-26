// include/ising/model.hpp
/*
For model parameters like E-field etc.
*/

#pragma once
#include <cmath>

namespace ising{
    struct Model{
        double J;

        constexpr Model(double J = 1.0)
        : J(J) {}
    }
}