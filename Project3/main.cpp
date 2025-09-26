// main.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"

int main() {

    double B0 = 1.0; 
    double V0 = 1.0; 
    double d = 1.0;  

    arma::vec pos = {0.0, 0.0, 0.0};
    arma::vec vel = {0.0, 0.0, 0.0};

    Particle p(1.0, 1.0, pos, vel);
    PenningTrap trap(B0, V0, d);
    trap.add_particle(p);

    p.print();
    trap.print_particles();
    return 0;
}
