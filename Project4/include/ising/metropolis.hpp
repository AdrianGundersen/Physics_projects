// include/ising/metropolis.hpp
/*
For the metropolis algorithm
*/

#pragma once
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include <array>

namespace ising{

    struct simParams{
        int total_steps;
        double temperature;
        int seed;
        int burn_in_sweeps;
        int measure_sweeps;
        int total_sweeps;
        int cores;
        int walkers;

        // file writing
        bool write_enabled;
        std::string write_type;
    };

    struct Boltzmannfactors{
        std::array<double, 5> factors;

        void set(double beta, double J){
            factors[0] = 1.0;
            factors[1] = 1.0;
            factors[2] = 1.0;
            factors[3] = std::exp(-beta * 4.0 * J);
            factors[4] = std::exp(-beta * 8.0 * J);
        }
    };
    void Metropolis(Model& model, Lattice& lattice, simParams& parmas, std::mt19937& rng);
}