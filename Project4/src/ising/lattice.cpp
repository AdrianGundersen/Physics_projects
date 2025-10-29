// src/ising/lattice.cpp
/*

*/

#include "ising/lattice.hpp"

namespace ising {
    void Lattice::init_spin_same(bool up) {
        if (up) {
            std::fill(spins.begin(), spins.end(), 1);     // all +1
        } else {
            std::fill(spins.begin(), spins.end(), -1);   // all -1
        }
    }
}