// include/ising/lattice.hpp
#pragma once
#include <armadillo>
#include <vector>
#include <random>
#include <omp_rng.hpp>

namespace ising {
    class Lattice {
        private:
            int L; // Lattice size
            double N; // Number of spins


            arma::Mat<int> spins; // 2D matrix to hold spin values (+1 or -1)


            void init_spin_from_mat(const arma::Mat<int> spin_mat) {
                arma::Mat<int> spins = spin_mat;
            }

            void init_spin_rand(int seed);

        public:
            Lattice(int L) : L(L), N(L * L), spins(L, L) {} // constructor


            // helper functions
            int size() const { return L; } // return lattice size L
            double num_spins() const { return N; } // return number of spins N
    };
}