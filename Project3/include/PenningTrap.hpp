// PenningTrap.hpp
#pragma once
#include <armadillo>
#include <vector>
#include "Particle.hpp"


class PenningTrap { 
    private:
        double B0; // Magnetic field strength
        double V0; // Electric potential
        double d;  // Characteristic dimension  
        bool coulomb_on; // Enable/disable Coulomb interaction

public:
    std::vector<Particle> particles;  // Store all particle objects

    PenningTrap(double B0, double V0, double d, bool coulomb_on); // constructor

    void add_particle(const Particle& p); // add a particle to trap

    void fill_random(int N, double q, double m, double max_vel); // fill trap with N random particles

    // External fields
    arma::vec external_E_field(const arma::vec& r) const;
    arma::vec external_B_field(const arma::vec& r) const;

    // Forces for particle i
    arma::vec force_external(int i) const;   // force from external fields
    arma::vec force_particle(int i, int j) const;  // force from j
    arma::vec total_force(int i) const;      // total force
    arma::mat acceleration_all(const arma::mat& R, const arma::mat& V) const; //temperary acceleration of partivcle

    // test functions
    int number_of_particles() const; // number of particles in trap

    // Debugging
    void print_particles() const; // print all particles
    void write_file(std::ofstream& ofile, double time, int particle_n) const; //write position, velocity and time to file
};


