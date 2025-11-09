#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <nlohmann/json.hpp>
#include "ising/lattice.hpp"
#include "ising/io/json_util.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include "ising/metropolis.hpp"


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
        std::cerr << "Usage: " << argv[0] << " <configs/L20.json>\n";
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
    ising::simParams params;
    ising::io::simparams_from_json(j.at("simulation"), j.at("lattice"), params);
    int seed = params.seed;
    std::mt19937 rng(seed);
    double T = params.temperature;
    int n_steps = params.total_steps;
    int measure_sweeps = params.measure_sweeps;
    int total_sweeps = params.total_sweeps;
    int N = lattice.num_spins();
    std::string spin_config;

    
    if (j.at("model").contains("spin_config")) {
        spin_config = j.at("model").at("spin_config").get<std::string>();
        if (spin_config == "all_up") {
            lattice.init_spin_same(true);
        } else if (spin_config == "all_down") {
            lattice.init_spin_same(false);
        } else if (spin_config == "random"){
            lattice.init_spin_rand(rng);
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

   


    std::string filename = "data/outputs/L20_T=" + std::to_string(T) + "_spin=" + spin_config + ".txt";
    std::ofstream ofile;
    ofile.open(filename);
    ofile << std::setprecision(10);
    ofile << ising::energy_per_spin(lattice, model) << ","
          << std::abs(magnetization_per_spin(lattice)) << ","
          << ising::energy_per_spin(lattice, model) * ising::energy_per_spin(lattice, model) << ","
          << std::abs(magnetization_per_spin(lattice)) * std::abs(magnetization_per_spin(lattice)) << "\n";



    std::vector<double> eps_samples, mabs_samples, eps2_samples, mabs2_samples;

    double mabs;
    double E = ising::total_energy(lattice, model);
    M = total_magnetization(lattice);
    for (int s = 0; s < total_sweeps; ++s) {    // her må metropolis kjøre N ganger for at info skal være riktig
        ising::Metropolis(model, lattice, params, rng, E, M);
        if (s % (measure_sweeps) == 0) { // after each sweep
            // std::cout << "Sampling at sweep " << s / N << "\n";
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

    std::cout << "After Metropolis sampling:\n";
    std::cout << "Average absolute magnetization per spin <|m|>: " << avg_mabs << "\n"; 
    std::cout << "Average energy per spin <ε>: " << avg_eps << "\n";
    std::cout << "Heat capacity C_V: " << heat_cap << "\n";
    std::cout << "Susceptibility χ: " << susc << "\n";

    eps = ising::energy_per_spin(lattice, model);
    M = total_magnetization(lattice);

    std::cout << "\nAfter Metropolis:\n";
    std::cout << "Total magnetization M: " << M << "\n";
    std::cout << "Total energy per spin ε: " << eps << "\n";


   
    write_to_file(ofile, eps_samples, mabs_samples, eps2_samples, mabs2_samples, 10);

    


    
}