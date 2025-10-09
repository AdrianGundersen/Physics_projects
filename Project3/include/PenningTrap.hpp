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
        double omega_V; // Angular frequency of time-dependent potential
        bool coulomb_on; // Enable/disable Coulomb interaction

public:
    std::vector<Particle> particles;  // Store all particle objects

    PenningTrap(double B0, double V0, double d, double f, bool coulomb_on); // constructor

    void add_particle(const Particle& p); // add a particle to trap
    void delete_particle(int particle_i); // delete particle

    void fill_random(int N, double q, double m, double pos_scaling, double vel_scaling); // fill trap with N random particles

    // External fields
    arma::vec external_E_field(const arma::vec& r, double& r_norm, double t = 0.0, double omega_V = 0.0) const;
    arma::vec external_B_field(const arma::vec& r, double& r_norm) const;

    // Forces for particle i
    arma::vec force_external(int i, double time, double omega_V = 0.0) const;   // force from external fields
    arma::vec force_particle(int i, int j) const;  // force from j
    arma::mat acceleration_all(const arma::mat& R, const arma::mat& V, double time, double omega_V = 0.0) const; // acceleration-matrix with column i being acceleration vector of particle i

    // test functions
    int number_of_particles(); // number of particles in trap (can add option to remove particles)

    // Debugging
    void print_particles() const; // print all particles
};


