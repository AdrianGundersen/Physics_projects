// apps/runtime.cpp
/*
To compare time for parallell and serial execution.
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>
#include <nlohmann/json.hpp>
#include <omp.h>

#include "ising/lattice.hpp"
#include "ising/io/json_util.hpp"
#include "ising/model.hpp"
#include "ising/observables.hpp"
#include "ising/metropolis.hpp"
#include "ising/mcmc_run.hpp"
#include "ising/io/write_files.hpp"
#include "omp_rng.hpp"

using json = nlohmann::json;


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


    double start_time = omp_get_wtime();
    ising::Result result = ising::mcmc_run(initial_lat, model, j);

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    // print time and number of cores
    int cores = j.at("simulation").value("cores", 1);
    std::cout << "Elapsed time: " << elapsed_time << " seconds\n";
    std::cout << "Number of cores: " << cores << "\n";

    return 0;
}