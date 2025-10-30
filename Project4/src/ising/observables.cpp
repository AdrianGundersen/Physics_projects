// src/ising/observables.cpp
/*

*/

#include "ising/observables.hpp"
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include <cmath>

namespace ising{
    // energy
    double total_energy(const Lattice& lat, const Model& model) {
        double J = model.J;
        int L = static_cast<int>(lat.size());
        double E = 0.0;
        bool double_count = model.double_count;
        int s;
        int right, down;
        int s_down, s_right;
        if (!double_count) {       
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                s = lat(i, j);
                right = (i+1) % L;
                down = (j+1) % L;
                s_right = lat(right, j);
                s_down = lat(i, down);
                E += s*(s_down + s_right);
            }
            }
        } 

        else {
        int s_left, s_up;
        int left, up;
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                s = lat(i, j);
                right = (i+1) % L;
                down = (j+1) % L;
                s_right = lat(right, j);
                s_down = lat(i, down);
                E += s*s_down;
                E += s*s_right;
                
                // also count left and up
                left = (i - 1 + L) % L;
                up = (j - 1 + L) % L;
                s_left = lat(left, j);
                s_up = lat(i, up);
                E += s*s_left;
                E += s*s_up;
            }
        }
        }
        E *= -J;
        return E;
    }

    double energy_per_spin(const Lattice& lat, const Model& model) {
        double N = lat.num_spins();
        double E_tot = total_energy(lat, model);
        return E_tot / N;
    }
    // magnetization
    double total_magnetization(const Lattice& lat) { 
        const int L = lat.size();
        double M = 0;
        for (int i = 0; i < L; ++i)
            for (int j = 0; j < L; ++j)
                M += lat(i, j);
        return M;
    }

    double magnetization_per_spin(const Lattice& lat) {
        double M_tot = total_magnetization(lat);
        double N = lat.num_spins();
        return M_tot / N;
    }
}
