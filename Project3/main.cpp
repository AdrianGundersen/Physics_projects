// main.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>

int main() {

    PenningTrap trap(parameters::B0, parameters::V0, parameters::d, parameters::coulomb_on);


    // origin particle
    arma::vec pos = {20.0, 0.0, 20.0};
    arma::vec vel = {0.0, 25.0, 0.0};

    Particle p(constants::elementary_charge, constants::atomic_mass_unit, pos, vel);
    trap.add_particle(p);

    // particle to the right
    arma::vec pos2 = {25.0, 25.0, 0.0};
    arma::vec vel2 = {0.0, 40.0, 5.0};

    Particle p2(1.0, 1.0, pos2, vel2);


    trap.add_particle(p2);

    trap.print_particles();

    double time = 0;

    std::filesystem::create_directory("data");

    std::ofstream ofile1("data/pos_vel_1.txt");
    std::ofstream ofile2("data/pos_vel_2.txt");
    trap.write_file(ofile1, time, 0);
    trap.write_file(ofile2, time, 1);

    for (int step = 0; step < parameters::N; step++) {
        time += parameters::dt;
        Integrator::RK4(trap, parameters::dt);

        trap.write_file(ofile1, time, 0);
        trap.write_file(ofile2, time, 1);
    }
    ofile1.close();
    ofile2.close();
    trap.print_particles();

    std::cout << "Wrote z position agains time as: data/pos_vel.txt\n";
    std::cout << parameters::dt << "\n";
    std::cout << time << "\n";
    return 0;
}
