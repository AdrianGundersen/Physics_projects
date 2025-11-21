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

        Grid(Index M, Real L) : M(M), L(L), h(L / (M - 1))  {}

        inline Index idx(Index i, Index j) const {
            return i * M + j;
        }
    };

    struct simParams {
        Real dt; // time step
        Index N; // number of time steps
    };

    struct SolverParams { // only jacobi for now
        Index max_iters; // max number of iterations for jacobi
        Real tol; // tolerance for convergence
    };

    struct PotentialParams {
        std::string type; // (harmonic, etc.)
        Real frequency; // for harmonic potential
    };
}