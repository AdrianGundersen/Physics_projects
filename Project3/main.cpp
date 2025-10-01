// main.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include <iostream>

int main() {

    PenningTrap trap(parameters::B0, parameters::V0, parameters::d);


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
    for (int step = 0; step < parameters::N; step++) {
        Integrator::RK4(trap, parameters::dt);
}
    trap.print_particles();
    return 0;
}
