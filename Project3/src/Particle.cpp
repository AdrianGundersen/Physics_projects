// Particle.cpp
#include "Particle.hpp"
#include <fstream>
#include <iomanip>


Particle::Particle(double charge, double mass, 
                   const arma::vec& position, 
                   const arma::vec& velocity)
    : charge(charge), mass(mass), position(position), velocity(velocity) {}

void Particle::print() const {
    std::cout << "Charge: " << charge << ", Mass: " << mass 
                << ", Position: [" << position.t() << "]"
                << ", Velocity: [" << velocity.t() << "]" << std::endl;
}

 void Particle::write_to_file(std::ofstream& ofile, double time, bool new_line) const {
    ofile << std::fixed << std::setprecision(12)
    << time << "\t" 
    << position(0) << "\t" << position(1) << "\t" << position(2) << "\t" 
    << velocity(0) << "\t" << velocity(1) << "\t" << velocity(2) << "\t";
    if(new_line){
        ofile << "\n";}
}