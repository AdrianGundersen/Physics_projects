// include/potential.hpp
#pragma once
#include "constants.hpp"
#include <vector>

namespace ds {

    struct Potential {
        ds::rvec values;
    };

    void initialize_potential(ds::Potential& V, const ds::Grid& grid, const ds::PotentialParams& potential_params);
} // namespace ds