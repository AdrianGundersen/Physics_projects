// N_particles.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <iomanip>

int main() {
    using clock = std::chrono::steady_clock;
    const auto t_start = clock::now();

    std::cout << "Filling Penning trap with " << parameters::N_particles << " particles\n";

    arma::arma_rng::set_seed(parameters::seed);
    PenningTrap trap(parameters::B0, parameters::V0, parameters::d, parameters::frequency, parameters::omega_V, parameters::coulomb_on);

    trap.fill_random(parameters::N_particles, constants::elementary_charge, constants::atomic_mass_unit, parameters::maxvel);

    double time = 0;

    std::filesystem::create_directory("data");
    
    std::ofstream ofile1("data/particle_count.txt");
    ofile1 << std::fixed << std::setprecision(10);
    ofile1 << time << " " << trap.number_of_particles() << "\n";

    for (int step = 0; step < parameters::N; step++) {
        time += parameters::dt;
        Integrator::RK4(trap, parameters::dt, time);
        ofile1 << time << " " << trap.number_of_particles() << "\n";
    }
    ofile1.close();
    const auto t_end = clock::now();
    const std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << std::fixed << std::setprecision(6)
              << "Ran for " << elapsed.count() << " seconds\n";

    return 0; 
}
