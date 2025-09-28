// main.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"

int main() {

    PenningTrap trap(constants::B0, constants::V0, constants::d);


    // origin particle
    arma::vec pos = {0.0, 0.0, 0.0};
    arma::vec vel = {0.0, 0.0, 0.0};

    Particle p(1.0, 1.0, pos, vel);
    trap.add_particle(p);

    // particle to the right
    arma::vec pos2 = {1.0, 0.0, 0.0};
    arma::vec vel2 = {1.0, 1.0, 0.0};

    Particle p2(1.0, 1.0, pos2, vel2);


    trap.add_particle(p2);

    p.print();
    trap.print_particles();


    for (int step = 0; step < 100; step++) {
    Integrator::ForwardEuler(trap, constants::dt);
}
    trap.print_particles();
    return 0;
}
