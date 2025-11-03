// src/ising/mcmc_run.cpp
/*
Runs a Markov Chain Monte Carlo with Metropolis algorithm with different walkers in parallel.
*/
#include <vector>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/io/json_util.hpp"
#include "ising/metropolis.hpp"
#include "ising/observables.hpp"
#include "ising/mcmc_run.hpp"

namespace ising{
    Result mcmc_run(const Lattice& initial_lat, const Model& model, const nlohmann::json& j) {
        int n_walkers = omp_get_max_threads();
        std::vector<Walker> walkers(n_walkers);

        ising::simParams params;
        ising::io::simparams_from_json(j.at("simulation"),j.at("lattice"), params);
        int mother_seed = params.seed;
        std::mt19937 rng(mother_seed);
        double T = params.temperature;
        int n_steps = params.total_steps;
        int measure_sweeps = params.measure_sweeps;
        int total_sweeps = params.total_sweeps;
        int N = params.total_steps;
        int cores = params.cores;
        int walker_n = params.walkers;
        ising::Result result;
        #pragma omp parallel
        for (int w = 0; w < walker_n; ++w) {
            int thread_id = omp_get_thread_num();
            // std::mt19937 thread_rng(rng() + thread_id); // different seed for each thread
            Lattice lat = initial_lat; // copy initial lattice (maybe consider randomizing each walker?)        }

    }
    return result;
    }}