// include/ising/io/write_files.hpp
/*
Write results to files
*/

#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <nlohmann/json.hpp>
#include "ising/observables.hpp"
#include "ising/mcmc_run.hpp"

namespace ising::io {

inline void write_results_to_file(const nlohmann::json& jwrite,
                                  const ising::Result& result,
                                  const std::string& filename) {
    const std::string type   = jwrite.value("type", "txt");
    const std::string delim  = jwrite.value("delimiter", ",");
    const int precision      = jwrite.value("precision", 10);
    const bool header        = jwrite.value("header", true);
    const std::string mode   = jwrite.value("average_or_concatenate", "concatenate");
    const std::vector<std::string> observables =
        jwrite.value("observables", std::vector<std::string>{"energy_per_spin", "magnetization_per_spin"});

    if (type != "txt") {
        std::cerr << "Unknown file type: " << type << "\n";
        return;
    }

    std::ofstream ofile(filename);
    if (!ofile) {
        std::cerr << "Could not open file for writing: " << filename << "\n";
        return;
    }

    ofile << std::fixed << std::setprecision(precision);

    auto write_row = [&](const std::vector<double>& eps, const std::vector<double>& mabs, int i) {
        for (std::size_t k = 0; k < observables.size(); ++k) {
            const std::string& obs = observables[k];
            if (obs == "energy_per_spin" || obs == "energy" || obs == "e" || obs == "eps") {
                ofile << eps[i];
            } else if (obs == "magnetization_per_spin" || obs == "magnetization" || obs == "m" || obs == "mabs") {
                ofile << mabs[i];
            } else {
                // Unknown observable: write empty field to keep column count stable
                ofile << "";
            }
            if (k + 1 < observables.size()) ofile << delim;
        }
        ofile << "\n";
    };

    if (header) {
        for (std::size_t k = 0; k < observables.size(); ++k) {
            ofile << "# " << observables[k];
            if (k + 1 < observables.size()) ofile << delim;
        }
        ofile << "\n";
    }
    if (mode == "average") {
        const auto& w = result.avg_walker;
        const int n = std::min<int>(w.eps_samples.size(), w.mabs_samples.size());
        for (int i = 0; i < n; ++i) write_row(w.eps_samples, w.mabs_samples, i);
    } else if (mode == "concatenate") {
        const auto& walkers = result.all_walkers;
        for (const auto& w : walkers) {
            const int n = std::min<int>(w.eps_samples.size(), w.mabs_samples.size());
            for (int i = 0; i < n; ++i) write_row(w.eps_samples, w.mabs_samples, i);
        }
    } else {
        std::cerr << "Unknown average_or_concatenate mode: " << mode << "\n";
    }

    ofile.close();
}

} // namespace ising::io
