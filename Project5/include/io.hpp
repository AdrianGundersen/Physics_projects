//include/io.hpp
/*
1. params_from_json read simulation parameters from a JSON object and populate the corresponding structures.
2. Function to write probability density to a file. To use uncomment in solver.cpp
3. Function to write wavefunction to a file. 
*/
#pragma once
#include <string>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {
   
    void params_from_json(const nlohmann::json& j, ds::simParams& params, ds::Grid& grid, ds::SolverParams& solver_params, ds::PotentialParams& potential_params, ds::OutputParams& output_params);
    
    void write_prob_to_file(const std::string& filename, const ds::rvec& prob_density, Index timestep, Index precision = 6);
    
    void write_wavefunction_to_file(const std::string& filename, const ds::cvec& wavefunction, Index timestep, Index precision = 6);
}