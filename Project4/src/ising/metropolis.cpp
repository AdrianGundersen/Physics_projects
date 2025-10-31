// src/ising/metropolis.cpp
/*
*/

#include "ising/metropolis.hpp"
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include <random>
#include <iostream>

namespace ising{
    void Metropolis(Model& model, Lattice& lattice, simParams& params, std::mt19937& generator) {
        const double T = params.temperature;
        const int n_steps = params.total_steps; // default N steps = N spins = 1 sweep
        const int seed = params.seed;
        
        const int L = lattice.size();
        const int N = lattice.num_spins();
        std::uniform_int_distribution<int> dist_pos(0, N-1); // select random spin
        std::uniform_real_distribution<double> dist_r(0.0, 1.0); // acceptance prob
        
        Boltzmannfactors boltz;
        const double beta = 1.0 / T;
        const double J = model.J;   
        boltz.set(beta, J); // precompute Boltzmann factors
        
        int i, j;
        int right, left, up, down;
        
        int s;
        int n;        

        for (int step = 0; step < n_steps; ++step) {
            int idx = dist_pos(generator);

            i = static_cast<int>(idx % L);
            j = static_cast<int>(idx / L);

            s = lattice(i, j);
            up = (i + 1) % L;
            down = (i + L - 1) % L; 
            right = (j + 1) % L;
            left = (j + L - 1) % L; // right/left in cols

            n = s * (lattice(up, j) + lattice(down, j) + lattice(i, right) +lattice(i, left));  //(-4, -2, 0, 2, 4)
            
            int factor_idx = (n+4) / 2; // n to idx in Boltzmann factors
            double r = dist_r(generator);
            if (r <= boltz.factors[factor_idx]) {
                lattice(i, j) = -s; // flip spin

                // std::cout << "Flipped spin at (" << i << ", " << j << ")\n";
            }

    }
    }
}