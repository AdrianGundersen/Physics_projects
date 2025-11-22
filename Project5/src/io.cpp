// src/io.cpp

#include "io.hpp"
#include <fstream>
#include <iostream>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {

    void params_from_json(const nlohmann::json& j, ds::simParams& params, 
                        ds::Grid& grid, ds::SolverParams& solver_params, 
                        ds::PotentialParams& potential_params, std::string& filename) {
        try {
            // Simulation parameters
            nlohmann::json simulation_json = j.at("simulation");
            params.dt = simulation_json.at("dt").get<ds::Real>();
            params.N = simulation_json.at("N").get<ds::Index>();
            params.xC = simulation_json.at("xC").get<ds::Real>();
            params.sigma_x = simulation_json.at("sigma_x").get<ds::Real>();
            params.px = simulation_json.at("px").get<ds::Real>();
            params.yC = simulation_json.at("yC").get<ds::Real>();
            params.py = simulation_json.at("py").get<ds::Real>();
            params.sigma_y = simulation_json.at("sigma_y").get<ds::Real>();

            // Grid parameters
            nlohmann::json grid_json = j.at("grid");
            ds::Index M = grid_json.at("M").get<ds::Index>();
            ds::Real L = grid_json.at("L").get<ds::Real>();
            grid = ds::Grid(M, L);

            // Solver parameters
            nlohmann::json solver_json = j.at("solver");
            solver_params.max_iters = solver_json.at("max_iterations").get<ds::Index>();
            solver_params.tol = solver_json.at("tolerance").get<ds::Real>();

            // Potential params
            nlohmann::json potential_json = j.at("potential");
            potential_params.type = potential_json.at("type").get<std::string>();
            potential_params.frequency = potential_json.at("frequency").get<ds::Real>();
            potential_params.V0 = potential_json.at("V0").get<ds::Real>();
            // Slits params (part of potential)
            nlohmann::json slits_json = j.at("slits");
            potential_params.slits.enabled = slits_json.at("enabled").get<bool>();
            potential_params.slits.wall_center = slits_json.at("wall_center").get<ds::Real>();
            potential_params.slits.wall_thickness = slits_json.at("wall_thickness").get<ds::Real>();
            potential_params.slits.slit_aperture = slits_json.at("slit_aperture").get<ds::Real>();
            potential_params.slits.num_slits = slits_json.at("num_slits").get<ds::Index>();
            potential_params.slits.slit_starts = slits_json.at("slit_starts").get<std::vector<ds::Real>>();
            potential_params.slits.slit_ends = slits_json.at("slit_ends").get<std::vector<ds::Real>>();

            // Filename
            filename = j.at("output").at("file_name").get<std::string>(); 
        } catch (const nlohmann::json::exception& e) { // if mismatch or missing
            std::cerr << "Error reading JSON parameters: " << e.what() << std::endl;
            throw;
        }
    }

    void write_prob_to_file(const std::string& filename, const ds::rvec& prob_density, Index timestep) {
        std::ofstream file;
        file.open(filename, std::ios::app);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        file << "Timestep " << timestep << ":\n";
        for (const auto& val : prob_density) {
            file << val << "\n";
        }
        file << "\n"; // new time step
        file.close();
    }
        
} // namespace ds 