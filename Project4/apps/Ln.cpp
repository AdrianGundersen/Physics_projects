#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <nlohmann/json.hpp>
#include <omp.h>
#include <cmath>
#include "ising/lattice.hpp"
#include "ising/io/json_util.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include "ising/metropolis.hpp"
#include "ising/mcmc_run.hpp"
#include "omp_rng.hpp"

using json = nlohmann::json;
using namespace ising;


void write_to_file(std::ofstream& ofile,
                const std::vector<double>& v1,
                const std::vector<double>& v2,
                const std::vector<double>& v3,
                const std::vector<double>& v4,
                int precision = 10){

    ofile << std::setprecision(precision);
    int size = std::size(v1);
    for (int i = 0; i < size; ++i) {
        ofile << v1[i] << ","
              << v2[i] << ","
              << v3[i] << ","
              << v4[i] << "\n";
    }
    ofile.close();
}

int main(int argc, char** argv) { // argc and argv to get JSON file path (argc is number of arguments, argv is array of arguments)
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configs/multiple_walkers.json>\n";
        return 1;
    }
    
    const std::string json_path = argv[1];

    std::ifstream f(json_path);
    if (!f) {
        std::cerr << "Could not open JSON file: " << json_path << "\n";
        return 1;
    }

    json j; f >> j; // parse JSON file into json object

    ising::Model model;
    ising::io::model_from_json(j.at("model"), model); // populate model


    std::string spin_config;
    if (j.at("model").contains("spin_config")) {
        spin_config = j.at("model").at("spin_config").get<std::string>();
    }
    else {
        std::cerr << "No spin config specified.\n";
        return 1;
    }


    ising::simParams params;
    ising::io::simparams_from_json(j.at("simulation"), j.at("lattice"), params);
    int seed = params.seed;
    std::mt19937 rng(seed);
    double T = params.temperature;
    int n_steps = params.total_steps;
    int measure_sweeps = params.measure_sweeps;
    int total_sweeps = params.total_sweeps;
    int burn_in_sweeps = params.burn_in_sweeps;
    int N = params.total_steps;
    int L = std::sqrt(N); // I aint pulling all that
    int cores = params.cores;
    int walker = params.walkers;

    std::vector<uint32_t> seeds = omp_rng::initialize_omp_rng_container(seed, walker);


    std::string filename = "data/outputs/L"+ std::to_string(L) + "_T=" + std::to_string(T) + "_spin=" + spin_config + ".txt";
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::setprecision(10);
    ofile.close();


  
    std::vector<ising::Walker> w(walker);
    std::vector<std::mt19937> rng_vec;
    for (int i = 0; i < walker; i++) {
        rng_vec.push_back(std::mt19937(seeds[i]));

    }
#pragma omp parallel for schedule(static)
for (int i = 0; i < walker; i++) {
    ising::Lattice lat = ising::io::lattice_from_json(j.at("lattice")); 
    auto rng_i = rng_vec[i];   // thread-local RNG

    if (spin_config == "random")        lat.init_spin_rand(rng_i);
    else if (spin_config == "all_up")   lat.init_spin_same(true);
    else if (spin_config == "all_down") lat.init_spin_same(false);
    else { std::cerr << "Unknown spin configuration\n"; continue; }

    double eps  = ising::energy_per_spin(lat, model);
    double mabs = std::abs(magnetization_per_spin(lat));
    w[i].eps_samples.push_back(eps);
    w[i].mabs_samples.push_back(mabs);

    for (int b = 0; b < burn_in_sweeps; ++b)
        ising::Metropolis(model, lat, params, rng_i);

    for (int s = 0; s < total_sweeps; ++s) {
        ising::Metropolis(model, lat, params, rng_i);
        if (s % measure_sweeps == 0) {
            eps  = ising::energy_per_spin(lat, model);
            mabs = std::abs(magnetization_per_spin(lat));
            w[i].eps_samples.push_back(eps);
            w[i].mabs_samples.push_back(mabs);
            w[i].n += 1;
        }
    }
}
    // ofile <<  << ","
    //       << std::abs(magnetization_per_spin(lattice)) << ","
    //       << ising::energy_per_spin(lattice, model) * ising::energy_per_spin(lattice, model) << ","
    //       << std::abs(magnetization_per_spin(lattice)) * std::abs(magnetization_per_spin(lattice)) << "\n";



    


    // double avg_eps = ising::average(eps_samples);
    // double avg_mabs = ising::average(mabs_samples);
    // double avg_eps2 = ising::average(eps2_samples);
    // double avg_mabs2 = ising::average(mabs2_samples);
    // double heat_cap = ising::heat_capacity(lattice, avg_eps2, avg_eps, T);
    // double susc = ising::susceptibility(lattice, avg_mabs2, avg_mabs, T);

    // std::cout << "After Metropolis sampling:\n";
    // std::cout << "Average absolute magnetization per spin <|m|>: " << avg_mabs << "\n"; 
    // std::cout << "Average energy per spin <ε>: " << avg_eps << "\n";
    // std::cout << "Heat capacity C_V: " << heat_cap << "\n";
    // std::cout << "Susceptibility χ: " << susc << "\n";

    // eps = ising::energy_per_spin(lattice, model);
    // M = total_magnetization(lattice);

    // std::cout << "\nAfter Metropolis:\n";
    // std::cout << "Total magnetization M: " << M << "\n";
    // std::cout << "Total energy per spin ε: " << eps << "\n";


   
    // write_to_file(ofile, eps_samples, mabs_samples, eps2_samples, mabs2_samples, 10);



    return 0;
}