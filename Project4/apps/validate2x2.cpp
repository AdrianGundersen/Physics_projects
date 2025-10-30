// apps/validate2x2.cpp
/*
Validate 2x2 against analytical results
*/

#include <iostream>
#include <fstream>
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

    double M = total_magnetization(lattice);
    std::cout << "Total magnetization M: " << M << "\n";


    double eps = ising::energy_per_spin(lattice, model);

    std::cout << "Total energy per spin ε: " << eps << "\n";

    ising::simParams params;
    ising::io::simparams_from_json(j.at("simulation"), params);
    int seed = params.seed;
    std::mt19937 rng(seed);
    double T = params.temperature;
    int n_steps = params.total_steps;
    int burn_in = params.burn_in_sweeps;
    int measure_steps = params.measure_sweeps;
    int N = lattice.num_spins();
    
    std::cout << "\nBurn-in sweeps: " << burn_in << ", Measure sweeps: " << measure_steps << "\n";

    // Running measure- and burn-in sweeps multiplied by N (number of spins)
    for (int s = 0; s< burn_in * N; ++s) {
        ising::Metropolis(model, lattice, params, rng);
    }
    std::vector<double> eps_samples, mabs_samples;

    for (int s = 0; s< measure_steps * N; ++s) {
        ising::Metropolis(model, lattice, params, rng);
        if (s % N == 0) { // after each sweep
            double eps = ising::energy_per_spin(lattice, model);
            double mabs = std::abs(magnetization_per_spin(lattice));
            eps_samples.push_back(eps);
            mabs_samples.push_back(mabs);
        }
    }
    double avg_eps = 0.0;
    double avg_mabs = 0.0;
    for (double e : eps_samples) {
        avg_eps += e;
    }
    avg_eps /= eps_samples.size();
    for (double m : mabs_samples) {
        avg_mabs += m;
    }
    avg_mabs /= mabs_samples.size();

    std::cout << "After Metropolis sampling:\n";
    std::cout << "Average absolute magnetization per spin <|m|>: " << avg_mabs << "\n"; 
    std::cout << "Average energy per spin <ε>: " << avg_eps << "\n";


    eps = ising::energy_per_spin(lattice, model);
    M = total_magnetization(lattice);

    std::cout << "\nAfter Metropolis:\n";
    std::cout << "Total magnetization M: " << M << "\n";
    std::cout << "Total energy per spin ε: " << eps << "\n";

    const double J = model.J; 
    const double beta = 1.0 / T;
    double analytical_Z = 12.0 + 4.0 * std::cosh(8.0 * J * beta);
    double analytical_eps = -(32.0 *J) / N  * (std::sinh(8 * beta * J)) / analytical_Z;
    double analytical_mabs = 8.0 / N *(std::exp(8.0 * J / T) + 2.0) / analytical_Z;

    std::cout << "\nAnalytical results:\n";
    std::cout << "Analytical average absolute magnetization per spin <|m|>: " << analytical_mabs << "\n";   
    std::cout << "Analytical average energy per spin <ε>: " << analytical_eps << "\n";

    return 0;
}


