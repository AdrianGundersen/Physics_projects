// include/solver.hpp

/*

High level:

1. SolverData struct to hold current and next solution vector, RHS, and precomputed coeffs.
2. precompute_coeff function to calculate coefficients based on simulation params, grid, and potential.
3. build_RHS function to construct the right-hand side vector for the linear system.
4. jacobi_solve function to iteratively solve the linear system using the Jacobi method
5. crank_nicolson_step function to perform a single Crank-Nicolson time step.
6. simulation function to run the entire simulation loop. 

*/

#pragma once
#include "constants.hpp"
#include "potential.hpp"
#include <vector>

namespace ds::solver {

    struct SolverData {
        ds::cvec v_curr;
        ds::cvec v_next;
        ds::cvec rhs;
        ds::cvec invD;
        Complex beta;
    };
    

    void build_AB(const SolverData& data, const ds::Grid& grid, cvec& A, cvec& B); // build full A and B matrices as the task requires

    void precompute_coeff(SolverData& data, const ds::simParams& params, 
                        const ds::Grid& grid, const ds::Potential& V);

    void build_RHS(const SolverData& data, const ds::simParams& params, 
                   const ds::Grid& grid, const ds::Potential& V);

    void jacobi_solve(SolverData& data, const ds::Grid& grid, Index max_iters, Real tol);

    void crank_nicolson_step(SolverData& data,
                            const ds::simParams& params,
                            const ds::Grid& grid,
                            const ds::Potential& V,
                            const ds::SolverParams& solver_params);



} // namespace ds::solver

namespace ds {

    void initialize_wavefunction(solver::SolverData& data,
                                 const ds::Grid& grid,
                                 const ds::simParams& params);

    void probability_density(const ds::solver::SolverData& data, ds::rvec& prob_density, const ds::Grid& grid); 

    void simulation(const simParams& params,
                    const Grid& grid,
                    const SolverParams& solver_params,
                    PotentialParams& potential_params,
                    const std::string& filename);
} // namespace ds

