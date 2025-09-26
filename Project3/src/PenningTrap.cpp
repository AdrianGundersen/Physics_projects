// PenningTrap.cpp
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include <cmath>


PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
// Constructor
    : B0(B0_in), V0(V0_in), d(d_in) {} 

// Add a particle
void PenningTrap::add_particle(const Particle& p) {
    particles.push_back(p); // add a copy of p to the particles vector
}

// External fields
arma::vec PenningTrap::external_E_field(const arma::vec& r) const {
    return arma::vec(); 
}

arma::vec PenningTrap::external_B_field(const arma::vec& r) const {
    return arma::vec(); 
}

// Forces
arma::vec PenningTrap::force_external(int i) const {
    return arma::vec(); 
}

arma::vec PenningTrap::force_particle(int i, int j) const {
    return arma::vec(); 
}

arma::vec PenningTrap::total_force(int i) const {
    return arma::vec(); 
}

// Debugging
void PenningTrap::print_particles() const {
    for (const Particle& p : particles) {
        p.print();
    }
}
