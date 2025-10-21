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
        double f;  // Frequency of time-dependent potential
        bool coulomb_on; // Enable/disable Coulomb interaction
        double omega_V; // Angular frequency of time-dependent potential

public:
    std::vector<Particle> particles;  // Store all particle objects

    PenningTrap(double B0, double V0, double d, double f, bool coulomb_on, double omega_V = 0.0); // constructor

    void add_particle(const Particle& p); // add a particle to trap
    void delete_particle(int particle_i); // delete particle

    void fill_random(int N, double q, double m, double pos_scaling, double vel_scaling); // fill trap with N random particles

    // External fields
    arma::vec external_E_field(const arma::vec& r, double& r_norm, double t = 0.0) const;
    arma::vec external_B_field(const arma::vec& r, double& r_norm) const;

    // Forces for particle i
    arma::vec force_external(int i, double time) const;   // force from external fields
    arma::vec force_particle(int i, int j) const;  // force from j

    /* 
    - Calculates the acceleration caused by the electric- and magnetic field
    - Calculates the acceleration due to the Coulomb interactions
    - Returns a matrix A, acceleration vector for a particle p is the columns of A
    */
    arma::mat acceleration_all(const arma::mat& R, const arma::mat& V, double time) const;

    /*
    - Calculates the kinetick energy
    - Calculates the Electrical potential in the trap
    - Calculates the coulumb potential between the particles
    - Returning a vector with the Total energy of each particle
    */
    arma::vec total_energy(bool print_info = true) const;

    // test functions
    int number_of_particles(); // number of particles in trap (can add option to remove particles)

    // Debugging
    void print_particles() const; // print all particles
};


