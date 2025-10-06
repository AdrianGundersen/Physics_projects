// N_particles.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"
#include "parameters.hpp"

#include <armadillo>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
    arma::arma_rng::set_seed(parameters::seed);

    arma::vec f_list = {0.1, 0.4, 0.7};

    // omega_V list in MHz
    double w_min = 0.20, w_max = 2.50, w_step = 0.02;
    int n_omega = static_cast<int>((w_max - w_min) / w_step + 1.5);
    arma::vec omega_V_list = arma::linspace(w_min, w_max, n_omega);

    bool coulomb_on = false;          // off 
    int steps = static_cast<int>(parameters::total_time / parameters::dt);

    std::filesystem::create_directory("data");

    const int nw = static_cast<int>(omega_V_list.n_elem);
    const int nf = static_cast<int>(f_list.n_elem);


    // total number of runs
    const int total_runs = nw * nf;
    double avg_duration = 0.0;


    int done = 0; // number of completed runs
    double total_elapsed_sec = 0.0; 



    auto global_start = std::chrono::steady_clock::now();

    arma::mat frac(nw, nf, arma::fill::zeros); // store results

    // sweep omega_V
    std::cout << "Starting simulations. Ready your CPU...\n";
#pragma omp parallel for collapse(2) schedule(dynamic) // parallelize outer two loops and schedule dynamically for load balancing
    for (int iw = 0; iw < nw; ++iw) {
        for (int i = 0; i < nf; ++i) {
            double omega_MHz = omega_V_list(iw);
            double omega = omega_MHz * 2.0 * M_PI;   // rad/microsecond

            auto start = std::chrono::steady_clock::now();

            double f = f_list(i);
            PenningTrap trap(parameters::B0, parameters::V0, parameters::d, f, coulomb_on);
            trap.fill_random(parameters::N_particles,
                             constants::elementary_charge,
                             constants::atomic_mass_unit,
                             parameters::maxvel);

            double t = 0.0;
            for (int k = 0; k < steps; ++k) {
                t += parameters::dt;
                Integrator::RK4(trap, parameters::dt, t, omega);
            }

            double frac_i = static_cast<double>(trap.number_of_particles()) /
                            static_cast<double>(parameters::N_particles);
            frac(iw, i) = frac_i;

            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            #pragma omp atomic // ensure atomic update
            total_elapsed_sec += elapsed.count();

            int done_now; // 
            #pragma omp atomic capture // ensure atomic increment and capture
            done_now = ++done;

            double avg = total_elapsed_sec / done_now; // average duration per run
            double remaining = avg * (total_runs - done_now); // estimated remaining time
            #pragma omp critical // ensure only one thread prints at a time
            {
                std::cout << "Run " << done_now << "/" << total_runs
                          << " (f=" << f
                          << ", omega_V=" << omega_MHz << " MHz)"
                          << "  Fraction trapped: " << frac_i
                          << "  Duration: " << std::fixed << std::setprecision(2)
                          << elapsed.count() << " s"
                          << "  Estimated remaining: " << std::setprecision(2)
                          << (remaining / 60.0) << " min\n";
            }
        }
    }

    std::ofstream out("data/fraction_trapped_vs_omega.txt");
    out << std::fixed << std::setprecision(6);
    out << "# omega_V_MHz";
    for (int i = 0; i < nf; ++i) out << " frac_f" << f_list(i);
    out << "\n";
    for (int iw = 0; iw < nw; ++iw) {
        out << omega_V_list(iw);
        for (int i = 0; i < nf; ++i) out << " " << frac(iw, i);
        out << "\n";
    }
    out.close();

    auto global_end = std::chrono::steady_clock::now();
    std::chrono::duration<double> total_elapsed = global_end - global_start;

    std::cout << "All simulations complete.\n";
    std::cout << "Total runtime: " << total_elapsed.count() / 60.0 << " min\n";
    std::cout << "Saved: data/fraction_trapped_vs_omega.txt\n";
    return 0;
}
