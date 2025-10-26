// include/ising/io/json_util.hpp

#pragma once
#include "ising/model.hpp"
#include "ising/lattice.hpp"
#include "ising/observables.hpp"

#include <string>
#include <nlohmann/json.hpp>
#include <cstdint>

namespace ising::io {
    // usage: model_from_json(root.at("model"), model);
    inline void model_from_json(const nlohmann::json& jm, ising::Model& m) {
        m.J = jm.value("J", 1.0);
    }

    inline ising::Lattice lattice_from_json(const nlohmann::json& jl) {
        const int L = jl.at("L").get<int>();
        ising::Lattice lat(L);
        return lat;
    }

    inline void observables_to_json(nlohmann::json& j, const ising::Observables& obs) {
        j["E"] = obs.E;
        j["E2"] = obs.E2;
        j["M"] = obs.M;
        j["M2"] = obs.M2;
    }
}