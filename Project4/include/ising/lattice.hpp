// include/ising/lattice.hpp

/*
For the Ising model. To represent the lattice and its spins.
*/
#pragma once
#include <vector>
#include <random>
#include <fstream>

namespace ising {
    class Lattice {
        private:
            int L; // Lattice size
            int N; // Number of spins


            std::vector<int> spins; // a flat vector to hold spin values (+1 or -1)

 

        public:
            Lattice() : L(0), N(0), spins() {}  // default constructor
            Lattice(int L) : L(L), N(L * L), spins(N) {} // constructor (assumes LxL lattice)
            void init_spin_from_mat(const std::vector<std::vector<int>>& spin_mat) {};

            void init_spin_rand(std::mt19937& rng); // initialize spins randomly

            void init_spin_same(bool up = true); // initialize all spins to same value (+1 by default)



            // operators

            // access spin value at (i, j)
            int& operator()(int i, int j) { return spins[i * L + j]; } // non-const version
            int  operator()(int i, int j) const { return spins[i * L + j]; } // const version

            // access spin value at flat vector idx
            int& operator[](int idx) { return spins[idx]; } // non-const version
            int  operator[](int idx) const { return spins[idx]; } // const version

            // helper functions
            int size() const { return L; } // return lattice size L
            double num_spins() const { return N; } // return number of spins N

            // QR-code to file
            inline void print(std::ofstream& os) const { // print lattice as "QR-code"
                for (int i = 0; i < L; ++i){
                    for (int j = 0; j < L; ++j){
                        int s = operator()(i, j);
                        char c = (s > 0 ? '1' : '0');
                        os << c;
                    }
                    os << "\n";
                }
            }
    };
}