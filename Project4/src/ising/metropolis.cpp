// src/ising/metropolis.cpp
/*
Implements the Metropolis algorithm for the Ising model.
Runs over a specified number of steps (standard is one sweep = N spins attempted flips).
Also updates energy and magnetization passed by reference so one does not need to recompute them each time.
*/

#include "ising/metropolis.hpp"
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include <random>
#include <iostream>

namespace ising{
    void Metropolis(Model& model, Lattice& lattice, simParams& params, std::mt19937& generator, double& E, int& M) {
        const double T = params.temperature;
        const int n_steps = params.total_steps; // default N steps = N spins = 1 sweep
        
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
        
        int s; // spin at (i,j)
        int n, sn, factor_idx;        

        double dE;
        int dM;
        for (int step = 0; step < n_steps; ++step) {
            int idx = dist_pos(generator);

            i = static_cast<int>(idx % L);
            j = static_cast<int>(idx / L);

            s = lattice(i, j);  
            up = (i + 1) % L;
            down = (i + L - 1) % L; 
            right = (j + 1) % L;
            left = (j + L - 1) % L; // right/left in cols

            n = lattice(up, j) + lattice(down, j) + lattice(i, right) + lattice(i, left);
            sn = s * n;
            factor_idx = (sn + 4) / 2; // n to idx in Boltzmann factors {-4,-2,0,2,4}

            dE = 2.0 * J * sn;
            dM = -2.0 * s;

        
            bool accept = (dE <= 0.0); // delta eps = 0, -4J, -8J -> always accept
            if (!accept) { // if dE > 0
                double r = dist_r(generator); // random number (0,1)
                accept = (r <= boltz.factors[factor_idx]); // delta eps = 4J, 8J -> accept if rng says so
            }
            if (accept){ // flip spin and update lattice, E, M
                lattice(i, j) = -s; // flip spin
                E += dE;
                M += dM;
            }
        }
    }
}