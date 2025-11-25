//include/io.hpp

#pragma once
#include <string>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {
    void params_from_json(const nlohmann::json& j, ds::simParams& params, ds::Grid& grid, ds::SolverParams& solver_params, ds::PotentialParams& potential_params, std::string& filename, std::string& filename_wavefunction);
    void write_prob_to_file(const std::string& filename, const ds::rvec& prob_density, Index timestep);
    void write_wavefunction_to_file(const std::string& filename, const ds::cvec& wavefunction, Index timestep);
}