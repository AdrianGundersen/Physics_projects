// include/ising/model.hpp
/*
For model parameters of the Ising model. Coupling constant J and whether to double count interactions. (double count not implemented)
*/

#pragma once
#include <cmath>
#include <string>

namespace ising {
    struct Model {
        double J;
        bool double_count;
        std::string spin_config;

        Model(double J = 1.0, bool double_count = false, std::string spin_config = "random") // cant use constexpr with std::string
            : J(J),
              double_count(double_count),
              spin_config(std::move(spin_config)) // std::move to avoid copy
            {}
    };
}