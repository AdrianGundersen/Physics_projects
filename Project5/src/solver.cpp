// src/solver.cpp

#include "solver.hpp"
#include "constants.hpp"
#include "potential.hpp"
#include "io.hpp"
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <omp.h>

namespace ds {
    void initialize_wavefunction(ds::solver::SolverData& data,
                                 const ds::Grid& grid, const ds::simParams& params)
    {
        const ds::Index M = grid.M;
        const ds::Real h = grid.h;

        const Real xC = params.xC;
        const Real sigma_x = params.sigma_x;
        const Real px = params.px;
        const Real yC = params.yC;
        const Real py = params.py;
        const Real sigma_y = params.sigma_y;

        ds::Real norm_factor = 0.0;

        // Fill v_curr with unormalized Gaussian
        //#pragma omp parallel for collapse(2) reduction(+:norm_factor)
        for (ds::Index i = 0; i < M; ++i) {
             const ds::Real x = static_cast<ds::Real>(i) * h;
             for (ds::Index j = 0; j < M; ++j) {
                const ds::Real y = static_cast<ds::Real>(j) * h;
                const ds::Index k = grid.idx(i, j);

                // Dirichlet BCs
                if (i == 0 || i == M - 1 || j == 0 || j == M - 1) {
                    data.v_curr[k] = ds::Complex(0.0, 0.0);
                    continue;
                }
                
                const ds::Real dx = x - xC;
                const ds::Real dy = y - yC;

                const ds::Real gauss = std::exp(-0.5 * ((dx * dx) / (sigma_x * sigma_x) + (dy * dy) / (sigma_y * sigma_y)));
                const ds::Complex phase = std::exp(ds::Complex(0.0, (px * dx + py*dy)));

                const ds::Complex psi_ij = gauss * phase;
                data.v_curr[k] = psi_ij;
                norm_factor += std::norm(psi_ij);
             }
        }
        norm_factor = std::sqrt(norm_factor * h * h); // normalization factor
        if (norm_factor > 0.0) {
            #pragma omp parallel for
            for (ds::Index k = 0; k < grid.size(); ++k) {
                data.v_curr[k] /= norm_factor;
            }
        }
    }

    void probability_density(const ds::solver::SolverData& data, ds::rvec& prob_density, const ds::Grid& grid) { // |psi|^2
        const ds::Index size = data.v_curr.size();
        #pragma omp parallel for
        for (ds::Index k = 0; k < size; ++k) {
            prob_density[k] = std::norm(data.v_curr[k]) * grid.h * grid.h; // multiply by area element
        }
    }
} // namespace

namespace ds::solver {
    void precompute_coeff(SolverData& data,
                          const ds::simParams& params,
                          const ds::Grid& grid,
                          const ds::Potential& V) {
        const Index M  = grid.M;
        const Index MM = M * M;

        data.v_curr.assign(MM, ds::Complex(0.0, 0.0));
        data.v_next.assign(MM, ds::Complex(0.0, 0.0));
        data.rhs.assign(MM, ds::Complex(0.0, 0.0));
        data.invD.assign(MM, ds::Complex(0.0, 0.0));

        const Real    h = grid.h;
        const Real    dt = params.dt;
        const Complex I(0.0, 1.0);
        const Complex lambda = I * (dt / (2.0 * h * h));

        data.beta = -lambda;

        for (Index i = 1; i < M-1; ++i) {
            for (Index j = 1; j < M-1; ++j) { // skip boundaries
                const Index k = grid.idx(i, j);

                // if (i == 0 || i == M - 1 || j == 0 || j == M - 1) {
                //     data.invD[k] = Complex(0.0, 0.0);
                //     continue; // Dirichlet BCs
                // }

                const Real Vij = V.values[k];
                const Complex alpha = Complex(1.0, 0.0) + 4.0 * lambda + I * (0.5 * dt * Vij);

                data.invD[k] = Complex(1.0, 0.0) / alpha;
            }
        }
    }
    void build_RHS(SolverData& data,
                   const ds::simParams& params,
                   const ds::Grid& grid,
                   const ds::Potential& V)
    {
        const Index M = grid.M;

        const cvec& u   = data.v_curr; // current solution as reference
        cvec&       rhs = data.rhs;       // right-hand side

        const Real    h  = grid.h; // grid spacing
        const Real    dt = params.dt; // time step
        const Complex I(0.0, 1.0);
        const Complex lambda = I * (dt / (2.0 * h * h));

        #pragma omp parallel for collapse(2)
        for (Index i = 0; i < M; ++i) {
            for (Index j = 0; j < M; ++j) { 
                const Index k = grid.idx(i, j); // idx in array

                if (i == 0 || i == M - 1 || j == 0 || j == M - 1) {
                    rhs[k] = u[k]; // Dirichlet BCs
                    continue;
                }
                const Real Vij = V.values[k];
                const Complex center = u[k]; // u_ij

                const Complex lap = // laplacian
                    u[grid.idx(i + 1, j)] +
                    u[grid.idx(i - 1, j)] +
                    u[grid.idx(i, j + 1)] +
                    u[grid.idx(i, j - 1)] -
                    Complex(4.0, 0.0) * center; 

                const Complex pot_term = I * (0.5 * dt * Vij) * center;

                rhs[k] = center + lambda * lap - pot_term; // RHS: u_ij + lambda*Laplacian - i*(0.5*dt*V_ij)*u_ij
            }
        }
    }


    void jacobi_solve(SolverData& data, const ds::Grid& grid, Index max_iters, Real tol) {
        const Index M = grid.M;
        
        const cvec& rhs = data.rhs; // right-hand side
        cvec& v = data.v_next; // next solution
        const cvec& invD = data.invD; // precomputed inverse diagonal
        
        if (v.size() != data.v_curr.size()) {
            v = data.v_curr; // initialize v_next if not same size (first call)
        }
        cvec v_old = v; // store old solution
        
        for (Index iter = 0; iter < max_iters; ++iter) {
            v_old = v; // store old solution
            Real max_diff = tol + 1.0; // maximum of all differences
            
            #pragma omp parallel for collapse(2) reduction(max:max_diff)
            for (Index i = 1; i < M - 1; ++i) { // skip boundaries
                for (Index j = 1; j < M - 1; ++j) { // skip boundaries
                    const Index k = grid.idx(i, j);
                    
                    
                    const Complex sumN = // sum of neighbors
                    v_old[grid.idx(i + 1, j)] +
                    v_old[grid.idx(i - 1, j)] +
                    v_old[grid.idx(i, j + 1)] +
                    v_old[grid.idx(i, j - 1)];
                    
                    const Complex new_val = invD[k] * (rhs[k] - data.beta * sumN);
                    const Real diff = std::abs(new_val - v_old[k]);
                    if (diff > max_diff) {
                        max_diff = diff;
                    }
                    v[k] = new_val;
                }
            }
            if (max_diff < tol) {
                break; // converged
            }
        }
    }

    void gauss_seidel_solve(SolverData& data, const ds::Grid& grid, Index max_iters, Real tol) {
    }

    void crank_nicolson_step(SolverData& data,
        const ds::simParams& params,
        const ds::Grid& grid,
        const ds::Potential& V,
        const ds::SolverParams& solver_params) 
        {
            build_RHS(data, params, grid, V);
            data.v_next = data.v_curr;
            jacobi_solve(data, grid, solver_params.max_iters, solver_params.tol);
            data.v_curr.swap(data.v_next);
        }
    } // namespace ds::solver
    
    namespace ds {
        
        void simulation(const simParams& sim_params,
            const Grid& grid,
            const SolverParams& solver_params,
            PotentialParams& potential_params,
            const std::string& filename,
            const std::string& filename_wavefunction)
            {
                Potential V;
                initialize_potential(V, grid, potential_params);
                
                solver::SolverData data;
                solver::precompute_coeff(data, sim_params, grid, V);
                ds::initialize_wavefunction(data, grid, sim_params);
                
                omp_set_num_threads(sim_params.threads);
                // ds::probability_density(data, prob_density, grid);
                ds::write_wavefunction_to_file(filename_wavefunction, data.v_curr, 0); // initial wavefunction
        ds::rvec prob_density(grid.size()); // allocates
        for (Index n = 1; n < sim_params.N; ++n) {
            solver::crank_nicolson_step(data, sim_params, grid, V, solver_params);
            ds::probability_density(data, prob_density, grid);
            ds::write_wavefunction_to_file(filename_wavefunction, data.v_curr, n);

            //uncomment to directryly calculate and write probability density to file
            //ds::write_prob_to_file(filename, prob_density, n);
            
        }
    }

} // namespace ds