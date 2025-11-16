// include/ising/io/write_files.hpp
/*
Write results to files
*/

#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "ising/io/json_util.hpp"
#include "ising/observables.hpp"
#include "ising/mcmc_run.hpp"

namespace ising::io {


void write_results_to_file(const nlohmann::json& jin,
                                  const ising::Result& result,
                                  const std::string& filename,
                                  const double T = 1.0);
} // namespace ising::io