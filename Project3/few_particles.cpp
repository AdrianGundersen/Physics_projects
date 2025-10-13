// few_particles.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>

int main() {
    arma::arma_rng::set_seed(parameters::seed);


    // Extract parameters
    auto sim_params = parameters::few;
    double total_time = sim_params.total_time;
    double dt = sim_params.dt;
    int N = sim_params.N;
    bool coulomb_on = sim_params.coulomb_on;

    PenningTrap trap(parameters::B0, parameters::V0, parameters::d, parameters::frequency, coulomb_on);

    // origin particle
    arma::vec pos = {20.0, 0.0, 20.0};
    arma::vec vel = {0.0, 25.0, 0.0};

    Particle p(constants::elementary_charge, constants::Ca_mass, pos, vel);
    trap.add_particle(p);

    // particle to the right
    arma::vec pos2 = {25.0, 25.0, 0.0};
    arma::vec vel2 = {0.0, 40.0, 5.0};

    Particle p2(constants::elementary_charge, constants::Ca_mass, pos2, vel2);


    trap.add_particle(p2);

    trap.print_particles();

    double time = 0;

    std::filesystem::create_directory("data");
    Particle& par1 = trap.particles[0];
    Particle& par2 = trap.particles[1];

    std::string filepath1 = "data/" + parameters::filename_few0;
    std::string filepath2 = "data/" + parameters::filename_few1;

    std::ofstream ofile1(filepath1);
    std::ofstream ofile2(filepath2);

    // makes headesr
    ofile1 << "# t x y z vx vy vz" << "\n";
    ofile2 << "# t x y z vx vy vz" << "\n";

    par1.write_to_file(ofile1, time, true);
    par2.write_to_file(ofile2, time, true);
    std::cout << "Starting simulation for N=" << N << " over simulation time t = " << total_time << " microseconds \n";
    for (int step = 0; step < N; step++) {
        time += dt;
        Integrator::RK4(trap, dt, time);
        par1.write_to_file(ofile1, time, true);
        par2.write_to_file(ofile2, time, true);
    }
    ofile1.close();
    ofile2.close();
    trap.print_particles();

    std::cout << "Wrote positions against time as: "<< filepath1 << "\n";
    std::cout << "Wrote positions against time as: "<< filepath2 << "\n";
    std::cout << dt << "\n";
    std::cout << time << "\n";
    return 0;
}
