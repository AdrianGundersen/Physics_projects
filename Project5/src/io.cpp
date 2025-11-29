// src/io.cpp

#include "io.hpp"
#include <fstream>
#include <iostream>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {

    void params_from_json(const nlohmann::json& j, ds::simParams& params, 
                        ds::Grid& grid, ds::SolverParams& solver_params, 
                        ds::PotentialParams& potential_params, ds::OutputParams& output_params) {
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
            solver_params.method = solver_json.at("method").get<std::string>();
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
            potential_params.slits.slit_spacing = slits_json.at("slit_spacing").get<ds::Real>();

            // Output parameters
            output_params.filename_prob = j.at("output").at("file_name").get<std::string>();
            output_params.filename_wavefunction = j.at("output").at("file_name_wavefunction").get<std::string>();
            output_params.precision = j.at("output").at("precision").get<ds::Index>();
            output_params.write_interval = j.at("output").at("write_interval").get<ds::Index>();

            // Paralellization parameters
            nlohmann::json parallelization_json = j.at("parallelization");
            params.threads = parallelization_json.at("threads").get<ds::Index>();

        } catch (const nlohmann::json::exception& e) { // if mismatch or missing
            std::cerr << "Error reading JSON parameters: " << e.what() << std::endl;
            throw;
        }
    }

    void write_prob_to_file(const std::string& filename, const ds::rvec& prob_density, Index timestep, Index precision) {
        std::ofstream file;
        std::ios_base::openmode mode = std::ios::app;
        if (timestep == 0) {
            mode = std::ios::trunc; // overwrite file on first timestep
        }
        file.open(filename, mode);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        file << "Timestep " << timestep << ":\n";
        file << std::fixed << std::setprecision(precision);
        for (const auto& val : prob_density) {
            file << val << "\n";
        }
        file << "\n"; // new time step
        file.close();
    }
    void write_wavefunction_to_file(const std::string& filename, const ds::cvec& wavefunction, Index timestep, Index precision) {
        std::ofstream file;
        std::ios_base::openmode mode = std::ios::app; // append mode
        if (timestep == 0) { // if first timestep overwrite
            mode = std::ios::trunc; // overwrite file on first timestep
        }
        file.open(filename, mode); // append mode
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        file << "Timestep " << timestep << ":\n";
        file << std::fixed << std::setprecision(precision);
        for (const auto& val : wavefunction) {
            file << val.real() << "," << val.imag() << "\n";
        }
        file << "\n"; // new time step
        file.close();
    }
        
} // namespace ds 