// include/ising/metropolis.hpp
/*
For the metropolis algorithm
*/

#pragma once
#include "ising/lattice.hpp"
#include "ising/model.hpp"

namespace ising{
    void Metropolis(Model& model, Lattice& lattice, int n_steps, double T, int seed);
}