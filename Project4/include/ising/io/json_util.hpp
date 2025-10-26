// include/ising/io/json_util.hpp

#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <nlohmann/json.hpp>

namespace ising::io {
    // usage: model_from_json(root.at("model"), model);
    inline void model_from_json(const nlohmann::json& jm, Model& m) {
        m.J = jm.value("J", 1.0);
}
}