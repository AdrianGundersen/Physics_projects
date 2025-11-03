// include/ising/io/write_files.hpp
/*
Write results to files
*/

#pragma once
#include <string>
#include <nlohmann/json.hpp>
#include "ising/observables.hpp"
#include "ising/mcmc_run.hpp"
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>

namespace ising::io {
    void write_results_to_file(nlohmann::json& jwrite, const ising::Result& result, const std::string& filename) {
        std::string type = jwrite.value("type", "txt");
        std::string delimiter = jwrite.value("delimiter", ",");
        int precision = jwrite.value("precision", 10);
        std::vector<std::string> observables = jwrite.value("observables", std::vector<std::string>{"energy", "magnetization"});
        std::cout << "Writing results to file: " << filename << "\n";
        if (type == "txt") {
            std::ofstream ofile(filename);
            ofile << std::setprecision(precision);
            
            int n_samples = result.avg_walker.n;

            // write header
            for (const std::string& obs : observables) {
                ofile << obs << delimiter;
            }
            ofile << "\n";

            for (int i = 0; i < n_samples; ++i) {
                for (const std::string& obs : observables) {
                    // std::cout << "Writing observable: " << obs << "\n";
                    if (obs == "energy") {
                        ofile << result.avg_walker.eps_samples[i] << delimiter;
                    } else if (obs == "magnetization") {
                        ofile << result.avg_walker.mabs_samples[i] << delimiter;
                    }
                }
                ofile << "\n";
            }
            ofile.close();
        }
        // Additional file types can be implemented here
    }
}
