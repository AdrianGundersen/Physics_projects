// include/ising/mcmc_run.hpp
/*
Runs a Markov Chain Monte Carlo with Metropolis algorithm with different walkers in parallel.
*/
#pragma once
#include <vector>
#include <omp.h>
#include "ising/lattice.hpp
#include "ising/model.hpp"
#include "ising/io/json_util.hpp"
#include "ising/metropolis.hpp"
#include "ising/observables.hpp"

namespace ising {
    struct Walker { // results from each walker
        double sum_eps = 0.0;
        double sum_eps2 = 0.0;
        double sum_absM = 0.0;
        double sum_m = 0.0;
        double sum_m2 = 0.0;
        int n = 0;
    };

    struct Result { // overall result
        Walker avg_walker; // average over all walkers as walker structure
        std::vector<Walker> all_walkers; // individual walker results
    };
}
