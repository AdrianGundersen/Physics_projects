// src/ising/lattice.cpp
/*

*/

#include "ising/lattice.hpp"

namespace ising {
    void Lattice::init_spin_same(bool up) {
        if (up) {
            spins.fill(1);     // all +1
        } else {
            spins.fill(-1);   // all -1
        }
    }
}