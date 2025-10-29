// include/ising/metropolis.hpp
/*
For the metropolis algorithm
*/

#pragma once
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include <array>

namespace ising{

    struct Boltzmannfactors{
        std::array<double, 5> factors;

        void set(double beta, double J){
            factors[0] = std::exp(-beta * 8.0 * J); 
            factors[1] = std::exp(-beta * 4.0 * J);
            factors[2] = 1.0;
            factors[3] = std::exp(beta * 4.0 * J);
            factors[4] = std::exp(beta * 8.0 * J);
        }

    void Metropolis(Model& model, Lattice& lattice, int n_steps, double T, int seed);
    };
}