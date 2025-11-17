// include/ising/io/json_util.hpp

#pragma once
#include "ising/model.hpp"
#include "ising/lattice.hpp"
#include "ising/metropolis.hpp"
#include "ising/observables.hpp"
#include <string>
#include <nlohmann/json.hpp>
#include <cstdint>

namespace ising::io {
    // usage: model_from_json(root.at("model"), model);
    inline void model_from_json(const nlohmann::json& jm, ising::Model& m) {
        m.J = jm.value("J", 1.0); // default J=1.0
        m.double_count = jm.value("double_count", false); // default false
    }

    inline ising::Lattice lattice_from_json(const nlohmann::json& jl) {
        const int L = jl.at("L").get<int>();
        ising::Lattice lat(L);
        return lat;
    }

    inline void simparams_from_json(const nlohmann::json& js, const nlohmann::json& jl, ising::simParams& params) {
        const int L = jl.at("L").get<int>();
        const int N = L*L;
        if (js.value("total_steps", "N") == "N") {
            params.total_steps = N; // one sweep
        } else {
            params.total_steps = js.value("total_steps", N); // default N
        }

        params.temperature = js.value("temperature", 2.0); // default T=2.0
        params.use_Trange = js.value("use_Trange", false); // default false
        if (params.use_Trange && js.contains("Trange")) {
            const auto& jr = js.at("Trange");
            params.Tmin = jr.value("Tmin", 2.1); // default Tmin=1.0
            params.Tmax = jr.value("Tmax", 2.4); // default Tmax=3.0
            params.Tsteps = jr.value("Tsteps", 10); // default Tsteps=5
        }
        else {
            params.Tmin = params.temperature;
            params.Tmax = params.temperature;
            params.Tsteps = 1;
        }
        params.seed = js.value("seed", 67); // default seed=67
        params.burn_in_sweeps = js.value("burn_in_sweeps", 1000); // default 1000
        params.measure_sweeps = js.value("measure_sweeps", 5000);
        params.total_sweeps = js.value("total_sweeps", 10000);
        params.cores = js.value("cores", 1);
        params.walkers = js.value("walkers", 1);
        params.QR = js.value("qr", false); // default false
    }

    inline void write_to_file_from_json(const nlohmann::json& jwrite, ising::simParams& params) {
        params.write_enabled = jwrite.value("enabled", false); // default false
        params.write_type = jwrite.value("type", "txt"); // default txt
    }

    inline void observables_to_json(nlohmann::json& j, const ising::Observables& obs) {
        j["E"] = obs.E;
        j["M"] = obs.M;
    }

    inline void T_to_json(nlohmann::json& jout, const nlohmann::json& jin, 
            const double& T, const double& heat_cap, const double& chi, 
            const double& avg_eps = 0.0, const double& avg_mabs = 0.0, 
            const double& avg_eps2 = 0.0, const double& avg_mabs2 = 0.0, 
            const double& avg_eps3 = 0.0, const double& avg_mabs3 = 0.0,
            const double& avg_eps4 = 0.0, const double& avg_mabs4 = 0.0) {


        const int L = jin.at("lattice").at("L").get<int>();
        const int sweeps = jin.at("simulation").value("total_sweeps", 10000);
        const int walkers = jin.at("simulation").value("walkers", 1);

        double tol = 1e-5; // tolerance for overwriting existing T entry

        nlohmann::json& entry = jout[std::to_string(L)]; // order L : T : values
        
        // try to find existing T entry
        for (auto& [k, v] : entry.items()) { // k is T, v is values
            double Tk = std::stod(k); // string to double
            if (std::abs(Tk - T) < tol) {
                v["chi"]     = chi;
                v["Cv"]      = heat_cap;
                v["avg_eps"] = avg_eps;
                v["avg_mabs"] = avg_mabs;
                v["avg_eps2"] = avg_eps2;
                v["avg_mabs2"] = avg_mabs2;
                v["avg_eps3"] = avg_eps3;
                v["avg_mabs3"] = avg_mabs3;
                v["avg_eps4"] = avg_eps4;
                v["avg_mabs4"] = avg_mabs4;
                v["sweeps"]  = sweeps;
                v["walkers"] = walkers;
                return;
            }
        }

        // no match -> create a new exact T key
        nlohmann::json& T_entry = entry[std::to_string(T)];
        T_entry["chi"]     = chi;
        T_entry["Cv"]      = heat_cap;
        T_entry["avg_eps"] = avg_eps;
        T_entry["avg_mabs"] = avg_mabs;
        T_entry["avg_eps2"] = avg_eps2;
        T_entry["avg_mabs2"] = avg_mabs2;
        T_entry["avg_eps3"] = avg_eps3;
        T_entry["avg_mabs3"] = avg_mabs3;
        T_entry["avg_eps4"] = avg_eps4;
        T_entry["avg_mabs4"] = avg_mabs4;
        T_entry["sweeps"]  = sweeps;
        T_entry["walkers"] = walkers;
    }
}