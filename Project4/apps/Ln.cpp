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

// like numpy linspace
static std::vector<double> linspace(double start, double end,  int n){
    std::vector<double> v;
    v.reserve(std::max(1,n));
    if (n <= 1) {v.push_back(start); return v;}
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i){
        v.push_back(start + i * step);
    }
    return v;
}

int main(int argc, char** argv) { // argc and argv to get JSON file path (argc is number of arguments, argv is array of arguments)
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configs/*.json>\n";
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

    ising::simParams params;
    ising::io::simparams_from_json(j.at("simulation"), j.at("lattice"), params);
    const bool use_Trange = params.use_Trange;
    
    std::vector<double> T_values;
    if (use_Trange){
        T_values = linspace(params.Tmin, params.Tmax, params.Tsteps);

    // running the simulation for a single temperature if use_Trange is false
    } else {
        ising::Lattice initial_lat = ising::io::lattice_from_json(j.at("lattice"));
        double T = params.temperature;
        ising::Result result = ising::mcmc_run(initial_lat, model, j);

        // writing setup
        nlohmann::json write_json = j.at("write_to_file");
    
    
        // parameters for default filename
        //double T = j.at("simulation").at("temperature").get<double>();
        std::string spin_config = j.at("model").at("spin_config").get<std::string>();
        int L = initial_lat.size();
    
        std::string default_file_name = "data/outputs/L" + std::to_string(L)
                                        + "_T=" + std::to_string(T)
                                        + "_spin=" + spin_config + ".txt";
        std::string output_file_name;
        if (write_json.value("output_filename", "default") == "default") {
            output_file_name = default_file_name;
        } else {
            output_file_name = write_json.at("output_filename").get<std::string>();
        }
        ising::io::write_results_to_file(j, result, output_file_name); // takes the entire input json as argument

        double eps_sum = 0.0; double mabs_sum = 0.0;
        double eps2_sum = 0.0; double mabs2_sum = 0.0;
        int count = 0;
        for (const auto& w : result.all_walkers) {
                for (const auto& e : w.eps_samples) {
                    eps_sum += e;
                    eps2_sum += e * e;
                    count += 1;
                }
                for (const auto& m : w.mabs_samples) {
                    mabs_sum += m;
                    mabs2_sum += m * m;
                }
            }

            double avg_eps = eps_sum / count;
            double avg_mabs = mabs_sum / count;
            double avg_eps2 = eps2_sum / count;
            double avg_mabs2 = mabs2_sum / count;

            const double Cv   = ising::heat_capacity(result.all_walkers.front().lat, avg_eps2, avg_eps, T);
            const double chi = ising::susceptibility(result.all_walkers.front().lat, avg_mabs2, avg_mabs, T);
        std::cout << "Results for T = " << T << ":\n";
        std::cout << "c_V = " << Cv << "\nchi_s = " << chi << "\n";
        std::cout << "Results written to " << output_file_name << "\n";
        return 0;
    }

    std::cout << "\nRunning simulations for lattice size " << ising::io::lattice_from_json(j.at("lattice")).size() << "\n";
    std::cout << "Tmax = " << params.Tmax << ", Tmin = " << params.Tmin << ", Tsteps = " << params.Tsteps << "\n";

    for (double T : T_values){
        json jT = j;
        jT["simulation"]["temperature"] = T;

        ising::Lattice initial_lat = ising::io::lattice_from_json(jT.at("lattice"));
        ising::Result result = ising::mcmc_run(initial_lat, model, jT);
        ising::io::write_results_to_file(j, result, "data/outputs/tsweep_results.json", T);
    }
    return 0;
}