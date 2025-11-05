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
#include <filesystem>
#include <nlohmann/json.hpp>
#include "ising/io/json_util.hpp"
#include "ising/observables.hpp"
#include "ising/mcmc_run.hpp"

namespace ising::io {

inline void write_results_to_file(const nlohmann::json& jin,
                                  const ising::Result& result,
                                  const std::string& filename,
                                  const double T = 1.0) {
    const nlohmann::json& jwrite = jin.value("write_to_file", nlohmann::json::object()); // get "write" section

    const std::string type   = jwrite.value("type", "txt");
    const std::string delim  = jwrite.value("delimiter", ",");
    const int precision      = jwrite.value("precision", 10);
    const bool header        = jwrite.value("header", true);
    const std::string mode   = jwrite.value("average_or_concatenate", "concatenate");
    const std::vector<std::string> observables =
        jwrite.value("observables", std::vector<std::string>{"energy_per_spin", "magnetization_per_spin"});

    if (type != "txt" && type != "json") {
        std::cerr << "Unknown file type: " << type << "\n";
        return;
    }
    if (type == "txt") {
        std::ofstream ofile(filename);
        if (!ofile) {
            std::cerr << "Could not open file for writing: " << filename << "\n";
            return;
        }

        ofile << std::fixed << std::setprecision(precision);

        auto write_row = [&](const std::vector<double>& eps, const std::vector<double>& mabs, int i) {
            for (std::size_t k = 0; k < observables.size(); ++k) { // int did not work size_t needed as observables.size() is type size_t
                const std::string& obs = observables[k];
                if (obs == "Cv" || obs == "heat_capacity" || obs == "chi" || obs == "susc" || obs == "susceptibility") {
                    std::cerr << "Observable " << obs << " not supported in txt output. Use json output instead.\n";
                }
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
    } // end type txt

    if (type == "json") {
        nlohmann::json jout = nlohmann::json::object();
        if (std::ifstream in{filename}; in.good()) {
            try { in >> jout; } catch (...) { /* leave as empty object */ }
        } // copied from chatGPT (it never worked for me :( ))
        if (mode == "average") {
            std::cerr << "Warning! Not recommended.\n";
            const ising::Walker& w = result.avg_walker;
        // std::cout << "first epsilon of avg walker: " << w.eps_samples[0] << "\n";
        std::vector<double> eps = w.eps_samples;
        std::vector<double> mabs = w.mabs_samples;
        std::vector<double> eps2, mabs2;
        for (const auto& e : eps)  eps2.push_back(e * e);
        for (const auto& m : mabs) mabs2.push_back(m * m);

        double avg_eps = ising::average(eps);
        double avg_mabs = ising::average(mabs);
        double avg_eps2 = ising::average(eps2);
        double avg_mabs2 = ising::average(mabs2);
        // std::cout << "lattice size: " << w.lat.size() << ", avg energy per spin: " << avg_eps << ", avg abs magnetization per spin: " << avg_mabs << "\n";
        double heat_cap = ising::heat_capacity(w.lat, avg_eps2, avg_eps, jin.at("simulation").at("temperature").get<double>());
        double susc = ising::susceptibility(w.lat, avg_mabs2, avg_mabs, jin.at("simulation").at("temperature").get<double>());
        // std::cout << "Calculated heat capacity Cv: " << heat_cap << ", susceptibility chi: " << susc << "\n";
        ising::io::T_to_json(jout, jin, T, heat_cap, susc);

        std::ofstream out(filename);
        out << std::setw(jwrite.value("indent", 2)) << jout;
        out.close();
        }
        if (mode == "concatenate") {
            double eps_sum = 0.0; double mabs_sum = 0.0;
            double eps2_sum = 0.0; double mabs2_sum = 0.0;
            int count = 0;
            for (const auto& w : result.all_walkers) {
                for (const auto& e : w.eps_samples) {
                    eps_sum += e;
                    eps2_sum += e * e;
                    count += 1;
                }
                for (const auto& m : w.mabs_samples) {
                    mabs_sum += m;
                    mabs2_sum += m * m;
                }
            }

            double avg_eps = eps_sum / count;
            double avg_mabs = mabs_sum / count;
            double avg_eps2 = eps2_sum / count;
            double avg_mabs2 = mabs2_sum / count;

            const double Cv   = ising::heat_capacity(result.all_walkers.front().lat, avg_eps2, avg_eps, T);
            const double chi  = ising::susceptibility(result.all_walkers.front().lat, avg_mabs, avg_mabs2, T);
            ising::io::T_to_json(jout, jin, T, Cv, chi);

            std::ofstream out(filename);
            out << std::setw(jwrite.value("indent", 2)) << jout;
            out.close();
        }
    }
} // end write_results_to_file

    // inline void add_critical_temperature(const nlohmann::json& jwrite,
    //     const int L,
    //     const double& max_heat_cap,
    //     const double& max_susc,
    //     const double& T_max_heat_cap,
    //     const double& T_max_susc) {
    //     const std::string type   = jwrite.value("type", "txt");
    //     const std::string delim  = jwrite.value("delimiter", ",");
    //     const int precision      = jwrite.value("precision", 10);

    //     const std::string dir = "data/outputs/";
    //     std::string filename = dir + "critical_temperature.txt";
    //     std::filesystem::create_directories(dir);

    //     std::vector<std::string> lines;
    //     // looking for the line for the excisting lattice size L
    //     {
    //         std::ifstream in(filename);
    //         std::string line;
    //         while (std::getline(in, line)) {
    //             if (!line.empty()) lines.push_back(line);
    //         }
    //     }
    //     // new line to write
    //     std::ostringstream oss;
    //     oss << std::fixed << std::setprecision(precision)
    //         << L << delim 
    //         << T_max_heat_cap << delim 
    //         << max_heat_cap << delim 
    //         << T_max_susc << delim 
    //         << max_susc;
    //     const std::string new_line = oss.str();

    //     bool line_found = false;
    //     for (auto& line : lines) {
    //         size_t c = line.find(delim);
    //         if (c == std::string::npos) continue; // malformed line
    //         try {
    //             int existing_L = std::stoi(line.substr(0, c));
    //             if (existing_L == L) {
    //             line = new_line; // update line
    //             line_found = true;
    //             break;
    //             }
    //         }          
    //         catch (const std::invalid_argument& ia) {
    //             continue; // skip malformed line
    //         }
    //     }

    //     // If line not found, add new line
    //     if (!line_found) {
    //         lines.push_back(new_line);
    //     }
    //     // Write back to file
        
    //     std::ofstream out(filename, std::ios::trunc);
    //     for (const auto& line : lines) {
    //         out << line << "\n";
    //     }
        
    //     std::cout << "Critical temperatures written to " << filename << "\n";
    // }


} // namespace ising::io
