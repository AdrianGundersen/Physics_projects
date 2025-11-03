// include/ising/io/json_util.hpp

#pragma once
#include "ising/model.hpp"
#include "ising/lattice.hpp"
#include "ising/observables.hpp"
#include "ising/metropolis.hpp"
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
        params.seed = js.value("seed", 67); // default seed=67
        params.burn_in_sweeps = js.value("burn_in_sweeps", 1000); // default 1000
        params.measure_sweeps = js.value("measure_sweeps", 5000);
        params.total_sweeps = js.value("total_sweeps", 10000);
        params.cores = js.value("cores", 1);
        params.walkers = js.value("walkers", 1);
    }

    inline void write_to_file_from_json(const nlohmann::json& jwrite, ising::simParams& params) {
        params.write_enabled = jwrite.value("enabled", false); // default false
        params.write_type = jwrite.value("type", "txt"); // default txt
    }

    inline void observables_to_json(nlohmann::json& j, const ising::Observables& obs) {
        j["E"] = obs.E;
        j["M"] = obs.M;
    }
}