// src/ising/metropolis.cpp
/*
*/

#include "ising/metropolis.hpp"
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include <random>

namespace ising{
    void Metropolis(Model& model, Lattice& lattice, int n_steps, double T, int seed) {
        std::mt19937 generator(seed);
        const int L = lattice.size();
        const int N = lattice.num_spins();
        std::uniform_int_distribution<int> dist_pos(0, N-1); // select random spin
        std::uniform_real_distribution<double> dist_r(0.0, 1.0); // acceptance prob


        const double beta = 1.0 / T;
        const double J = model.J;
        
        int i, j;
        int right, left, up, down;
        double dE;
        for (int step = 0; step < n_steps; ++step) {
            int idx = dist_pos(generator);

            i = idx / L;
            j = idx % L;

            int s = lattice(i, j);
            right = (i + 1) % L;
            left  = (i - 1 + L) % L;
            down  = (j + 1) % L;
            up    = (j - 1 + L) % L;
        }

    }
}