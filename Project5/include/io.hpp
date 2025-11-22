//include/io.hpp

#pragma once
#include <string>
#include "constants.hpp"
#include <nlohmann/json.hpp>

namespace ds {
    void params_from_json(const nlohmann::json& j, ds::simParams& params, ds::Grid& grid, ds::SolverParams& solver_params, ds::PotentialParams& potential_params, std::string& filename);
    void write_prob_to_file(const std::string& filename, const ds::rvec& prob_density, Index timestep);
}