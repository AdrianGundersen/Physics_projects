// include/potential.hpp

/*
Functuin to initialize the potential based on the provided parameters.
Takes a Potential object to fill, the grid information, and potential parameters.
The potensial creates a box, a crosssection berrier with slits, or a harmonic oscillator potential.
*/

#pragma once
#include "constants.hpp"
#include <vector>

namespace ds {

    struct Potential {
        ds::rvec values;
    };

    void initialize_potential(ds::Potential& V, const ds::Grid& grid, 
                            const ds::PotentialParams& potential_params);
} // namespace ds