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
#include "ising/io/write_files.hpp"
#include "omp_rng.hpp"

using json = nlohmann::json;
using namespace ising;


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

    ising::Lattice initial_lat = ising::io::lattice_from_json(j.at("lattice"));
    ising::Result result = ising::mcmc_run(initial_lat, model, j);

    // writing setup
    nlohmann::json write_json = j.at("write_to_file");

    // parameters for default filename
    int L = initial_lat.size();
    double T = j.at("simulation").at("temperature").get<double>();
    std::string spin_config = j.at("model").at("spin_config").get<std::string>();

    std::string default_file_name = "data/outputs/L" + std::to_string(L)
                                    + "_T=" + std::to_string(T)
                                    + "_spin=" + spin_config + ".txt";
    std::string output_file_name;
    if (write_json.value("output_filename", "default") == "default") {
        output_file_name = default_file_name;
    } else {
        output_file_name = write_json.at("output_filename").get<std::string>();
    }

    ising::io::write_results_to_file(write_json, result, output_file_name);

    std::vector<double> eps, eps2 , mabs, mabs2;


    const auto& walkers = result.all_walkers;
    for (const auto& w : walkers){
        const int n = std::min<int>(w.eps_samples.size(), w.mabs_samples.size());
        for (int i = 0; i < n; ++i) {
            eps.push_back(w.eps_samples[i]);
            mabs.push_back(w.mabs_samples[i]);
            eps2.push_back(w.eps_samples[i] * w.eps_samples[i]);
            mabs2.push_back(w.mabs_samples[i] * w.mabs_samples[i]);
        }
    }
    double avg_eps = ising::average(eps);
    double avg_mabs = ising::average(mabs);
    double avg_eps2 = ising::average(eps2);
    double avg_mabs2 = ising::average(mabs2);
    double heat_cap = ising::heat_capacity(initial_lat, avg_eps2, avg_eps, T);
    double susc = ising::susceptibility(initial_lat, avg_mabs2, avg_mabs, T);

    std::cout << "avgerage absolute magnetization per spin <|m|>: " << avg_mabs << "\n";
    std::cout << "average energy per spin <ε>: " << avg_eps << "\n";
    std::cout << "average energy squared per spin <ε²>: " << avg_eps << "\n";
    std::cout << "average magnetization squared per spin <m²>: " << avg_mabs2 << "\n";
    std::cout << "heat capacity: " << heat_cap << "\n";
    std::cout << "susceptibility: " << susc << "\n";

    //sdt::cout << "heat capasity: " << result.heat_cap << "\n";
    //sdt::cout << "susceptibility: " << result.susc << "\n";
    return 0;
}