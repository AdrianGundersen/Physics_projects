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
    void Lattice::init_spin_rand(int seed) {
        std::mt19937_64 rng(seed);
        for (int i = 0; i < spins.size(); i++) {
            std::uniform_int_distribution<int> bit(0, 1);
            int spin = bit(rng) * 2 - 1; // 0->-1, 1->+1
            spins[i] = spin;
        }
    }
}