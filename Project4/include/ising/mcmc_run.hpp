// include/ising/mcmc_run.hpp
/*
Runs a Markov Chain Monte Carlo with Metropolis algorithm with different walkers in parallel.
*/
#pragma once
#include <vector>
#include <random>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/io/json_util.hpp"
#include "ising/metropolis.hpp"
#include "ising/observables.hpp"


namespace ising {
    struct Walker { // results from each walker
        double E = 0.0; // energy temporary variable
        double Mabs = 0.0; // abs magnetization temporary variable
        std::vector<double> eps_samples; // energy samples
        std::vector<double> mabs_samples; // absolute magnetization samples
        std::vector<double> eps2_samples; // energy samples squared
        std::vector<double> mabs2_samples; // absolute magnetization samplessquared
        double heat_cap = 0.0; // heat capacity
        double susc = 0.0; // susceptibility
        int n = 0; // number of measurements
        Lattice lat; // lattice state
        Model model; // model parameters
        std::mt19937 rng;
    };

    struct Result { // overall result
        Walker avg_walker; // average over all walkers as walker structure
        std::vector<Walker> all_walkers; // individual walker results
    };
    Result mcmc_run(const Lattice& initial_lat, const Model& model, const nlohmann::json& j);
}
