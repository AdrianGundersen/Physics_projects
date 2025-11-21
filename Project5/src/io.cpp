// src/io.cpp

#include "io.hpp"
#include <fstream>
#include <iostream>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {

    void params_from_json(const nlohmann::json& j, ds::simParams& params, ds::Grid& grid, ds::SolverParams& solver_params, ds::PotentialParams& potential_params) {
        try {
            // Simulation parameters
            params.dt = j.at("dt").get<ds::Real>();
            params.N = j.at("N").get<ds::Index>();

            // Grid parameters
            ds::Index M = j.at("M").get<ds::Index>();
            ds::Real L = j.at("L").get<ds::Real>();
            grid = ds::Grid(M, L);

            // Solver parameters
            solver_params.max_iters = j.at("max_iterations").get<ds::Index>();
            solver_params.tol = j.at("tolerance").get<ds::Real>();

            // Potential params
            potential_params.type = j.at("potential").at("type").get<std::string>();
            potential_params.frequency = j.at("potential").at("frequency").get<ds::Real>();

            
        } catch (const nlohmann::json::exception& e) { // if mismatch or missing
            std::cerr << "Error reading JSON parameters: " << e.what() << std::endl;
            throw;
        }
    }
        
} // namespace ds 