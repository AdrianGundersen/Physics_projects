// include/constants.hpp
#pragma once
#include <complex>
#include <vector>

namespace ds {
    // aliases for types we often use
    using Real = double;
    using Complex = std::complex<double>;
    using Index = size_t;

    using rvec = std::vector<Real>;
    using cvec = std::vector<Complex>;

    constexpr Real hbar = 1.0;
    constexpr Real mass = 1.0;
    struct Grid{
        Index M; // Number of grid points
        Real L; // Length of the domain
        Real h; // Grid spacing

        Grid() : M(0), L(0.0), h(0.0) {} // default constructor
        Grid(Index M, Real L) : M(M), L(L), h(L / (M - 1))  {}

        Index size() const { return M * M; } // total number of grid points

        inline Index idx(Index i, Index j) const {
            return i * M + j; // (i,j) to linear index
        }
    };

    struct simParams {
        Real dt; // time step
        Index N; // number of time steps
        Real xC; // initial x center
        Real sigma_x; // initial x width
        Real px; // initial x momentum
        Real yC; // initial y center
        Real py; // initial y momentum
        Real sigma_y; // initial y width
    };

    struct SolverParams { // only jacobi for now
        Index max_iters; // max number of iterations for jacobi
        Real tol; // tolerance for convergence
    };


    struct SlitsParams {
        bool enabled; // whether slits are enabled
        Real wall_center; // center position of the wall
        Real wall_thickness; // x-thickness of the wall
        Real slit_aperture; // height of each slit
        Index num_slits; // number of slits
        std::vector<Real> slit_starts; // starting y-positions of slits
        std::vector<Real> slit_ends; // ending y-positions of slits
    };

    struct PotentialParams {
        std::string type; // (harmonic, etc.)
        Real frequency; // for harmonic potential
        Real V0; // for constant potential
        SlitsParams slits; // parameters for slits potential
    };

}