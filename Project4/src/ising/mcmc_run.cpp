// src/ising/mcmc_run.cpp
/*
Runs a Markov Chain Monte Carlo with Metropolis algorithm with different walkers in parallel.

What it does high level:
1. Parse simulation parameters from JSON.
2. Initialize multiple walkers with different random seeds and initial spin configurations.
3. For each walker, perform burn-in sweeps followed by measurement sweeps, collecting observables.
4. Use OpenMP to parallelize across multiple walkers.
5. Return results including all walkers and their average.
*/
#include <vector>
#include <omp.h>
#include <nlohmann/json.hpp>
#include <random>
#include <iostream>
#include "ising/lattice.hpp"
#include "ising/model.hpp"
#include "ising/io/json_util.hpp"
#include "ising/metropolis.hpp"
#include "ising/observables.hpp"
#include "ising/mcmc_run.hpp"
#include "omp_rng.hpp"


namespace ising{
    Result mcmc_run(const Lattice& initial_lat, const Model& model, const nlohmann::json& j) {
        nlohmann::json lattice_json = j.at("lattice");
        nlohmann::json sim_json = j.at("simulation");
        nlohmann::json model_json = j.at("model");
        std::string spin_config = model_json.at("spin_config");

       
        ising::simParams params;
        ising::io::simparams_from_json(sim_json, lattice_json, params);
        int mother_seed = params.seed;
        int measure_sweeps = params.measure_sweeps;
        int total_sweeps = params.total_sweeps;
        int cores = params.cores;
        int n_walkers = params.walkers;
        ising::Result result;

        std::vector<uint32_t> seeds = omp_rng::initialize_omp_rng_container(mother_seed, n_walkers);
        std::vector<Walker> walkers(n_walkers); // vector with n_walkers number of walkers 
        std::vector<std::mt19937> rng_vec;
        rng_vec.reserve(n_walkers); 
        for (int i = 0; i < n_walkers; i++) {
            rng_vec.push_back(std::mt19937(seeds[i]));
            walkers[i].rng = rng_vec[i];
            // initialize lattice for each walker
            ising::Lattice lat = ising::io::lattice_from_json(lattice_json);
            
            if (spin_config == "random")        lat.init_spin_rand(rng_vec[i]);
            else if (spin_config == "all_up")   lat.init_spin_same(true);
            else if (spin_config == "all_down") lat.init_spin_same(false);
            else {std::cerr << "Unknown spin configuration\n"; continue; }
            walkers[i].lat = std::move(lat); // move to avoid large copy
            walkers[i].model = model; 
        }
        omp_set_num_threads(cores);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n_walkers; ++i) {
            ising::Walker& walker = walkers[i];
            std::mt19937& rng_i = walker.rng;
            ising::Lattice& w_lat = walker.lat;
            ising::Model& w_model = walker.model;
            
            double E = ising::total_energy(w_lat, w_model);
            int M = total_magnetization(w_lat); 
            // burn-in sweeps
            for (int b = 0; b < params.burn_in_sweeps; ++b) {
                ising::Metropolis(w_model, w_lat, params, rng_i, E, M); 
            }
            
            walker.eps_samples.clear(); // just in case not clear
            walker.mabs_samples.clear(); 
            walker.n = 0;
            walker.eps_samples.reserve(total_sweeps / measure_sweeps + 1); // preallocate memory
            walker.mabs_samples.reserve(total_sweeps / measure_sweeps + 1); 

            double eps;
            double mabs;

            E = ising::total_energy(w_lat, w_model);
            M = total_magnetization(w_lat); 

           const int recompute_interval = 50000; // recompute E and M every 50000 sweeps to avoid accumulation errors (arbitrary choice) 

            int total_spins = w_lat.num_spins();
            for (int s = 0; s < total_sweeps; ++s) {
                ising::Metropolis(w_model, w_lat, params, rng_i, E, M);
                if (s % measure_sweeps == 0) {
                    eps  = E / static_cast<double>(total_spins);
                    mabs = std::abs(M / static_cast<double>(total_spins));
                    walker.eps_samples.push_back(eps);
                    walker.mabs_samples.push_back(mabs);
                    walker.n += 1;
                }
                if (s % recompute_interval == 0) { // recompute to avoid accumulation
                    E = ising::total_energy(w_lat, w_model);
                    M = total_magnetization(w_lat);
                }
            }
        }
        // write overall results
        result.all_walkers = walkers;
        // compute average walker
        Walker avg_walker;
        avg_walker.lat = walkers[0].lat; // just copy from first walker
        avg_walker.model = walkers[0].model;
        avg_walker.n = walkers[0].n; // all walkers should have same number
        for (int i = 0; i < avg_walker.n; ++i) {
            double eps_sum = 0.0;
            double mabs_sum = 0.0;
            for (const auto& walker : walkers) {
                eps_sum += walker.eps_samples[i];
                mabs_sum += walker.mabs_samples[i];
            }
            avg_walker.eps_samples.push_back(eps_sum / n_walkers);
            avg_walker.mabs_samples.push_back(mabs_sum / n_walkers);
    }
    result.avg_walker = avg_walker;
    // std::cout << result.all_walkers.size() << " walkers completed." << std::endl;
    return result;
}
}