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
#include <thread>
#include <atomic>
#include <cmath>
#include <algorithm>  

#ifdef _OPENMP
#include <omp.h>
#endif

int main() {
    arma::arma_rng::set_seed(parameters::seed);


    // Extract parameters
    auto sim_params = parameters::multi;
    double dt = sim_params.dt;
    double total_time = sim_params.total_time;
    int N = sim_params.N;
    bool coulomb_on = sim_params.coulomb_on;
    

    // particle info
    int N_particles = parameters::N_particles; // number of particles
    double pos_scaling = parameters::pos_scaling; // position scaling factor
    double vel_scaling = parameters::vel_scaling; // velocity scaling factor

    std::filesystem::create_directory("data");
    std::string filepath = "data/" + parameters::filename_multi; // output file path

    // Frequency and amplitude parameters for multi particle simulations
    arma::vec f_list = parameters::f_list; // Amplitude factors
    int nf = f_list.n_elem; // number of f values
    arma::vec omega_V_list = parameters::omega_V_list; // omega_V list in MHz
    int nw = omega_V_list.n_elem; // number of omega_V values



    // total number of runs
    const int total_runs = nw * nf;
    double avg_duration = 0.0;


    int done = 0; // number of completed runs
    double total_elapsed_sec = 0.0; 



    auto global_start = std::chrono::steady_clock::now();

    arma::mat frac(nw, nf, arma::fill::zeros); // store results


    std::cout << "Starting simulations on " << N_particles << " particles.\n";
    int total_threads = omp_get_max_threads();
    int N_threads = std::min(total_threads, parameters::max_threads); // use at most max_threads
    std::cout << "Using up to " << N_threads << "/" << total_threads << " threads.\n";
    std::cout << "Total runs: " << total_runs << "\n";
    std::cout << "Coulomb interaction: " << (coulomb_on ? "ON" : "OFF") << "\n";
    std::cout << "Ready your CPU...\n";
    
    char response;
    std::cout << "Proceed with simulation? (y/n): ";
    std::cin >> response;

    if (response != 'y' && response != 'Y') {
        std::cout << "Simulation aborted by user.\n";
        return 0;
    }
    
#pragma omp parallel for num_threads(N_threads) collapse(2) schedule(dynamic) // parallelize outer two loops and schedule dynamically for load balancing
    for (int iw = 0; iw < nw; ++iw){
          // rad/microsecond
        for (int i = 0; i < nf; ++i) {
            double omega_MHz = omega_V_list(iw);
            double omega = omega_MHz * 2.0 * M_PI; 
            auto start = std::chrono::steady_clock::now();

            double f = f_list(i);
            PenningTrap trap(parameters::B0, parameters::V0, parameters::d, f, coulomb_on, omega);
            trap.fill_random(N_particles,
                             constants::elementary_charge,
                             constants::atomic_mass_unit,
                             pos_scaling, vel_scaling);

            double t = 0.0;
            for (int k = 0; k < N; ++k) {
                t += dt;
                Integrator::RK4(trap, dt, t, omega);
            }

            double frac_i = static_cast<double>(trap.number_of_particles()) /
                            static_cast<double>(N_particles);
            frac(iw, i) = frac_i;

            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            #pragma omp atomic // ensure atomic update
            total_elapsed_sec += elapsed.count();
            int done_now; // 
            #pragma omp atomic capture // ensure atomic increment and capture
            done_now = ++done;

            double avg = total_elapsed_sec / done_now; // average duration per run
            double remaining = avg * (total_runs - done_now) / N_threads; // estimated remaining time assuming perfect load balancing
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

    std::ofstream out(filepath);
    out << std::fixed << std::setprecision(6);
    out << "# omega_V_MHz"; // makes header
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
    std::cout << "Saved: " << filepath << "\n";
    return 0;
}
