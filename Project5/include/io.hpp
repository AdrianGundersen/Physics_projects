//include/io.hpp

#pragma once
#include <string>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {
   /*
   Function to read simulation parameters from a JSON object and populate the corresponding structures.
   Same structure as from ising model project.
   */
    void params_from_json(const nlohmann::json& j, ds::simParams& params, ds::Grid& grid, ds::SolverParams& solver_params, ds::PotentialParams& potential_params, ds::OutputParams& output_params);
    /*
    Function to write probability density to a file.
    The function takes the filename, probability density vector, timestep, and precision as inputs.
    To use uncomment the call for the function in solver.cpp
    */
    void write_prob_to_file(const std::string& filename, const ds::rvec& prob_density, Index timestep, Index precision = 6);
    /*
    Same as above but for wavefunction includeing real and complex numbers.
    This is used in solver.cpp and the probabilities as calculated in python.
    */
    void write_wavefunction_to_file(const std::string& filename, const ds::cvec& wavefunction, Index timestep, Index precision = 6);
}