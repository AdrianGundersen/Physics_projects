// include/ising/mcmc_run.hpp
/*
Runs a Markov Chain Monte Carlo with Metropolis algorithm with different walkers in parallel.
*/
#pragma once
#include <vector>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/io/json_util.hpp"
#include "ising/metropolis.hpp"
#include "ising/observables.hpp"


namespace ising {
    struct Walker { // results from each walker
        std::vector<double> eps_samples; // sum of energies
        std::vector<double> mabs_samples; // sum of absolute magnetizations
        std::vector<double> sum_m; // sum of magnetizations
        int n = 0; // number of measurements
        std::mt19937 rng;
    };

    struct Result { // overall result
        Walker avg_walker; // average over all walkers as walker structure
        std::vector<Walker> all_walkers; // individual walker results
    };
    Result mcmc_run(const Lattice& initial_lat, const Model& model, const nlohmann::json& j);
}
