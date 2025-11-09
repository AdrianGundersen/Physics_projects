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
    void Lattice::init_spin_rand(std::mt19937& rng) {
        int spin;
        int size = spins.size();
        std::uniform_int_distribution<int> bit(0, 1);
        for (int i = 0; i < size; i++) {
            spin = bit(rng) * 2 - 1; // 0->-1, 1->+1
            spins[i] = spin;
        }
    }
}