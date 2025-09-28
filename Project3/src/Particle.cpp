// Particle.cpp
#include "Particle.hpp"


Particle::Particle(double charge, double mass, 
                   const arma::vec& position, 
                   const arma::vec& velocity)
    : charge(charge), mass(mass), position(position), velocity(velocity) {}

void Particle::print() const {
    std::cout << "Charge: " << charge << ", Mass: " << mass 
                << ", Position: [" << position.t() << "]"
                << ", Velocity: [" << velocity.t() << "]" << std::endl;
}
