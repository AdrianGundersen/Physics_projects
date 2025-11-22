//src/main.cpp

/* 

Main overview of program:

1. Read parameters from JSON file
2. Intialize data structures
3. Run simulation

Everything else is handled under the hood in other files.
*/

#include "io.hpp"
#include "constants.hpp"
#include "potential.hpp"
#include "solver.hpp"
#include <iostream>
#include <nlohmann/json.hpp>
#include <fstream>


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

    nlohmann::json j; f >> j; // parse JSON file into json object

    ds::simParams sim_params;
    ds::Grid grid(0,0); // temporary initialization
    ds::SolverParams solver_params;
    ds::PotentialParams potential_params;
    std::string filename;

    ds::params_from_json(j, sim_params, grid, solver_params, potential_params, filename);
    ds::simulation(sim_params, grid, solver_params, potential_params, filename);

    return 0;
}