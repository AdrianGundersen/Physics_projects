// PenningTrap.cpp
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"
#include "parameters.hpp"
#include <cmath>
#include <algorithm> 

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
// Constructor
    : B0(B0_in), V0(V0_in), d(d_in) {} 

// Add a particle
void PenningTrap::add_particle(const Particle& p) {
    particles.push_back(p); // add a copy of p to the particles vector
}

// External fields
arma::vec PenningTrap::external_E_field(const arma::vec& r) const {
    double V0_over_d2 = V0 / (d * d);
    arma::vec Efield = V0_over_d2 * arma::vec({r(0), r(1), -2.0 * r(2)});
    return Efield;
}

arma::vec PenningTrap::external_B_field(const arma::vec& r) const {
    return arma::vec({0.0, 0.0, B0});
}

// Forces
arma::vec PenningTrap::force_external(int i) const {
    const Particle& p = particles[i];

    arma::vec E = external_E_field(p.position);
    arma::vec B = external_B_field(p.position);

    // Lorentz force
    return p.charge * (E + arma::cross(p.velocity, B));
}

arma::vec PenningTrap::force_particle(int i, int j) const {
    const Particle& pi = particles[i];
    const Particle& pj = particles[j];

    arma::vec r_vec = pi.position - pj.position;
    double r_norm = arma::norm(r_vec);
    r_norm = std::max(r_norm, parameters::EPS); // avoid division by zero

    // Coulomb force
    return constants::ke * pi.charge * pj.charge * r_vec / std::pow(r_norm, 3);
}

arma::vec PenningTrap::total_force(int i) const {
    arma::vec F = force_external(i);

    for (int j = 0; j < (int)particles.size(); j++) {
        if (j != i) {
            F += force_particle(i, j);
        }
    }
    return F;

}

// returns the acceleration of the particles at a temperary psoition
arma::mat PenningTrap::acceleration_all(const arma::mat& R, const arma::mat& V)
    const {
    int N = particles.size();
    arma::mat A(3, N, arma::fill::zeros);

    for (int i = 0; i < N; i++){
        const Particle& p = particles[i];
        arma::vec r = R.col(i);
        arma::vec v = V.col(i);

        arma::vec E = external_E_field(r);
        arma::vec B = external_B_field(r);
        arma::vec F = p.charge * (E + arma::cross(v, B));

        for (int j = 0; j < N; j++) {
            if (j != i) {
                arma::vec r_vec = r - R.col(j);
                double r_norm = arma::norm(r_vec);
                r_norm = std::max(r_norm, parameters::EPS); // avoid division by zero
                F += constants::ke * p.charge * particles[j].charge * r_vec / std::pow(r_norm, 3);
            }
        }
        A.col(i) =  F/p.mass;
}   return A;
}

// Debugging
void PenningTrap::print_particles() const {
    for (const Particle& p : particles) {
        p.print();
    }
}

void PenningTrap::write_file(std::ofstream& ofile, double time, int particle_n) const {
    const Particle& p = particles[particle_n];
    p.write_to_file(ofile, time); 
}