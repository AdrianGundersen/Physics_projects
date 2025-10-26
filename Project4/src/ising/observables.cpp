// src/ising/observables.cpp
/*

*/

#include "ising/observables.hpp"
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include <cmath>

namespace ising{
    double total_energy(const Lattice& lat, const Model& model) {
        double J = model.J;
        int L = static_cast<int>(lat.size());
        double E = 0.0;
        return E;
    }

    double energy_per_spin(const Lattice& lat, const Model& model) {
        return total_energy(lat, model) / lat.num_spins();
    }
}
