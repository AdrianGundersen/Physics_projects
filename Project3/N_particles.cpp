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
    std::ofstream out("data/fraction_trapped_vs_omega.txt");
    out << std::fixed << std::setprecision(6);

    // header
    out << "# omega_V_MHz";
    for (int i = 0; i < f_list.n_elem; ++i) out << " frac_f" << f_list(i);
    out << "\n";

    // total number of runs
    const int total_runs = omega_V_list.n_elem * f_list.n_elem;
    int run_count = 0;
    double avg_duration = 0.0;

    auto global_start = std::chrono::steady_clock::now();

    // sweep omega_V
    for (int iw = 0; iw < omega_V_list.n_elem; ++iw) {
        double omega_MHz = omega_V_list(iw);
        double omega = omega_MHz * 2.0 * M_PI;   // rad/microsecond
        out << omega_MHz;

        

        // for each amplitude f run a new trap simulation
        for (arma::uword i = 0; i < f_list.n_elem; ++i) {
            run_count++;
            auto start = std::chrono::steady_clock::now();
            double f = f_list(i);

            PenningTrap trap(parameters::B0, parameters::V0, parameters::d, f, coulomb_on);
            trap.fill_random(parameters::N_particles,
                             constants::elementary_charge,
                             constants::atomic_mass_unit,
                             parameters::maxvel);
                             
            std::cout << "Run " << run_count << "/" << total_runs
                      << " (f = " << f
                      << ", omega_V = " << omega_MHz << " MHz)\n";

                
            double t = 0.0;
            for (int k = 0; k < steps; ++k) {
                t += parameters::dt;
                Integrator::RK4(trap, parameters::dt, t, omega);
            }

            double frac = static_cast<double>(trap.number_of_particles()) /
                          static_cast<double>(parameters::N_particles);
            out << " " << frac;

            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            // Update running average
            avg_duration += (elapsed.count() - avg_duration) / run_count;
            double remaining = (total_runs - run_count) * avg_duration;
            std::cout << "  Fraction trapped: N=" << frac << "\n";
            std::cout << std::fixed << std::setprecision(2)
                      << "  Duration: " << elapsed.count() << " s"
                      << ", Estimated remaining: " << remaining / 60.0 << " min\n\n";
        }

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
