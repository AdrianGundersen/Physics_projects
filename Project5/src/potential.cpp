// src/potential.cpp

#include "potential.hpp"
#include "constants.hpp"
#include <cmath>


namespace ds {
    void initialize_potential(ds::Potential& V, const ds::Grid& grid, 
                            const ds::PotentialParams& potential_params) { 
        Index MM = grid.size();
        Index M = grid.M;
        Real mass = ds::mass;
        V.values.resize(MM);


        if (potential_params.type == "harmonic") {
            const Real omega = potential_params.frequency;

            for (Index i = 0; i < M; ++i) {
                const Real x = static_cast<Real>(i) * grid.h - 0.5 * grid.L;
                for (Index j = 0; j < M; ++j) {
                    const Real y = static_cast<Real>(j) * grid.h - 0.5 * grid.L;
                    const Real r2 = x * x + y * y;
                    V.values[grid.idx(i, j)] = 0.5 * mass * omega * omega * r2; // V = 0.5 * m * omega^2 * r^2
                }
            }
        }
        else {
            // Default: zero potential
            for (Index k = 0; k < MM; ++k) {
                V.values[k] = 0.0;
            }
        }
        const ds::SlitsParams Slit = potential_params.slits;
        
        if (Slit.enabled) {
            const Real V0 = potential_params.V0;
            const Real wall_x_start = Slit.wall_center - 0.5 * Slit.wall_thickness; // x-start of wall
            const Real wall_x_end = Slit.wall_center + 0.5 * Slit.wall_thickness; // x-end of wall
            const Real y_mid = 0.5 * grid.L; // middle of domain in y
            const Real total_span = (Slit.num_slits - 1) * Slit.slit_spacing;
            const Real first_center = y_mid - 0.5 * total_span; // center of first slit
            for (Index i = 0; i < M; ++i) {
                const Real x = static_cast<Real>(i) * grid.h; // x position
                if (x >= wall_x_start && x <= wall_x_end) { // within wall region
                    for (Index j = 0; j < M; ++j) {
                        const Real y = static_cast<Real>(j) * grid.h; 
                        bool in_slit = false; // check if in any slit
                        for (Index s = 0; s < Slit.num_slits; ++s) { // loop over slits
                            const Real slit_center = first_center + s * Slit.slit_spacing;
                            const Real half_apperture = 0.5 * Slit.slit_aperture;

                            if (std::abs(y - slit_center) <= half_apperture) { // y in [slit start, slit end]
                                in_slit = true;
                                break;
                            }
                        }
                        if (!in_slit) {
                            V.values[grid.idx(i, j)] = V0; // high potential for wall
                        }
                    }
                }
            }
    }
}
} // namespace ds