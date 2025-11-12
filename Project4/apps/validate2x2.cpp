// apps/validate2x2.cpp
/*
Validate 2x2 against analytical results
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <nlohmann/json.hpp>
#include "ising/lattice.hpp"
#include "ising/io/json_util.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include "ising/metropolis.hpp"


using json = nlohmann::json;
using namespace ising;

int main(int argc, char** argv) { // argc and argv to get JSON file path (argc is number of arguments, argv is array of arguments)
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configs/2x2.json>\n";
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

    ising::Lattice lattice = ising::io::lattice_from_json(j.at("lattice")); // populate lattice

    if (j.at("model").contains("spin_config")) {
        std::string spin_config = j.at("model").at("spin_config").get<std::string>();
        if (spin_config == "all_up") {
            lattice.init_spin_same(true);
        } else if (spin_config == "all_down") {
            lattice.init_spin_same(false);
        } else {
            std::cerr << "Unknown spin configuration: " << spin_config << "\n";
            return 1;
        }
    }
    else {
        std::cerr << "No spin config specified.\n";
        return 1;
    }

    std::cout << "Model J: " << model.J << "\n";
    std::cout << "Lattice size L: " << lattice.size() << "\n";
    std::cout << "Spin (0,0): " << lattice(0, 0) << "\n";

    int M = total_magnetization(lattice);
    std::cout << "Total magnetization M: " << M << "\n";


    double eps = ising::energy_per_spin(lattice, model);

    std::cout << "Total energy per spin ε: " << eps << "\n";

    ising::simParams params;
    ising::io::simparams_from_json(j.at("simulation"), j.at("lattice"), params);
    int seed = params.seed;
    std::mt19937 rng(seed);
    double T = params.temperature;
    int burn_in = params.burn_in_sweeps;
    int measure_sweeps = params.measure_sweeps;
    int total_sweeps = params.total_sweeps;
    
    std::cout << "\nBurn-in sweeps: " << burn_in << "\n";
    double start_time = omp_get_wtime();
    // Running measure- and burn-in sweeps multiplied by N (number of spins)
    double E = ising::total_energy(lattice, model);
    for (int s = 0; s< burn_in; ++s) {
        ising::Metropolis(model, lattice, params, rng, E, M);
    }
    std::cout << "Burn-in complete.\n"; 
    std::cout << "Starting simulation with " << total_sweeps - burn_in << " sweeps and\n";
    std::cout << "Measure sweeps: " << measure_sweeps << "\n";
    std::vector<double> eps_samples, mabs_samples, eps2_samples, mabs2_samples;
    
    double mabs;
    const int N = lattice.num_spins();
    E = ising::total_energy(lattice, model);
    M = total_magnetization(lattice);
    
    for (int s = 0; s < total_sweeps; ++s) {
        ising::Metropolis(model, lattice, params, rng, E, M);
        if (s % measure_sweeps == 0) { // after each sweep
            //std::cout << "Sampling at sweep " << s << "...\n";
            eps = ising::energy_per_spin(lattice, model);
            mabs = std::abs(magnetization_per_spin(lattice));
            eps_samples.push_back(eps);
            mabs_samples.push_back(mabs);
            eps2_samples.push_back(eps * eps);
            mabs2_samples.push_back(mabs * mabs);
        }
    }
    double avg_eps = ising::average(eps_samples);
    double avg_mabs = ising::average(mabs_samples);
    double avg_eps2 = ising::average(eps2_samples);
    double avg_mabs2 = ising::average(mabs2_samples);
    double heat_cap = ising::heat_capacity(lattice, avg_eps2, avg_eps, T);
    double susc = ising::susceptibility(lattice, avg_mabs2, avg_mabs, T);

    double end_time = omp_get_wtime();
    double total_time = end_time - start_time;

    std::cout << "\nSimulation took in total " << total_time << " seconds" << "\n";
    std::cout << "After Metropolis sampling:\n";
    std::cout << "Average absolute magnetization per spin <|m|>: " << avg_mabs << "\n"; 
    std::cout << "Average energy per spin <ε>: " << avg_eps << "\n";
    std::cout << "Average energy squared per spin <ε²>: " << avg_eps2 << "\n";
    std::cout << "Average magnetization squared per spin <m²>: " << avg_mabs2 << "\n";
    std::cout << "Heat capacity per spin C_V/N: " << heat_cap << "\n";
    std::cout << "Susceptibility χ per spin/N: " << susc << "\n";

    eps = ising::energy_per_spin(lattice, model);
    M = total_magnetization(lattice);

    std::cout << "\nAfter Metropolis:\n";
    std::cout << "Total magnetization M: " << M << "\n";
    std::cout << "Total energy per spin ε: " << eps << "\n";

    const double J = model.J; 
    const double beta = 1.0 / T;
    double analytical_Z = 12.0 + 4.0 * std::cosh(8.0 * J * beta);
    double analytical_eps = -(32.0 *J) / N  * (std::sinh(8.0 * beta * J)) / analytical_Z;
    double analytical_eps2 = (256.0 * J * J) / N / N * (std::cosh(8.0 * J * beta)) / analytical_Z;
    double analytical_mabs = 8.0 / N *(std::exp(8.0 * J *beta) + 2.0) / analytical_Z;
    double analytical_m2 =  32.0 /( N * N ) * (std::exp(8.0 * J * beta) + 1.0) / analytical_Z;

    double analytical_Cv = N / (T * T) *(analytical_eps2 - analytical_eps * analytical_eps);
    double analytical_chi = N  / T * (analytical_m2 - analytical_mabs * analytical_mabs);

    std::cout << "\nAnalytical results:\n";
    std::cout << "Analytical average absolute magnetization per spin <|m|>: " << analytical_mabs << "\n";   
    std::cout << "Analytical average energy per spin <ε>: " << analytical_eps << "\n";
    std::cout << "Analytical average energy squared per spin <ε²>: " << analytical_eps2 << "\n";
    std::cout << "Analytical average magnetization squared per spin <m²>: " << analytical_m2 << "\n";
    std::cout << "Analytical heat capacity C_V/N: " << analytical_Cv << "\n";
    std::cout << "Analytical susceptibility χ/N: " << analytical_chi << "\n";

    return 0;
}


