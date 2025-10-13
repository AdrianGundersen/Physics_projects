
// single_particle.cpp
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
    PenningTrap trap_rk(parameters::B0, parameters::V0, parameters::d, parameters::frequency, parameters::coulomb_on);
    PenningTrap trap_eu = trap_rk;

    // Extract parameters
    auto sim_params = parameters::single;
    double total_time = sim_params.total_time;
    double dt = sim_params.dt;
    int N = sim_params.N;
    bool coulomb_on = sim_params.coulomb_on;

    // origin particle
    arma::vec pos = {20.0, 0.0, 20.0};
    arma::vec vel = {0.0, 25.0, 0.0};

    Particle p1(constants::elementary_charge, constants::Ca_mass, pos, vel);
    Particle p2 = p1;
    trap_rk.add_particle(p1);
    trap_eu.add_particle(p2);

    double time = 0;

    std::filesystem::create_directory("data");
    std::string filepath = "data/" + parameters::filename_single;
    std::ofstream ofile(filepath);


    Particle& rk = trap_rk.particles[0];
    Particle& eu = trap_eu.particles[0];

    // make header
    ofile << "# t rk.x rk.y rk.z rk.vx rk.vy rk.vz eu.x eu.y eu.z eu.vx eu.vy eu.vz" << "\n";

    rk.write_to_file(ofile, time, false);
    eu.write_to_file(ofile, time, true);


    for (int step = 0; step < N; step++){
        time += dt;
        Integrator::RK4(trap_rk, dt, time);
        Integrator::ForwardEuler(trap_eu, dt, time);
        rk.write_to_file(ofile, time, false);
        eu.write_to_file(ofile, time, true);
    }
    ofile.close();

    std::cout << "Wrote single particle agains time and analytical as:" << filepath << "\n";
    std::cout << dt << "\n";
    std::cout << time << "\n";
    return 0;   
}
