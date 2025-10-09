// PenningTrap.cpp
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"
#include "parameters.hpp"
#include <cmath>
#include <algorithm> 

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f_in, bool coulomb_on_in)
// Constructor
    : B0(B0_in), V0(V0_in), d(d_in), f(f_in), coulomb_on(coulomb_on_in) {} 

// Add a particle
void PenningTrap::add_particle(const Particle& p) {
    particles.push_back(p); // add a copy of p to the particles vector
}

void PenningTrap::delete_particle(int particle_i) {
    particles.erase(particles.begin() + particle_i);
}

void PenningTrap::fill_random(int N, double q, double m, double pos_scaling, double vel_scaling) {
    for (int i = 0; i < N; i++) {
        arma::vec pos = arma::vec(3).randn() * pos_scaling * d; // typical pos is pos_scaling * d
        arma::vec vel = arma::vec(3).randn() * vel_scaling * d; // typical vel is vel_scaling * d / microsecond
        Particle p(q, m, pos, vel);
        add_particle(p);
    }
}


// External fields
arma::vec PenningTrap::external_E_field(const arma::vec& r, double& r_norm, double t, double omega_V) const {
    if (r_norm > d) {
        return arma::vec({0.0, 0.0, 0.0});  // no field outside the trap
    }
    if (f != 0.0) {
        double V_t = V0 * (1 + f * std::cos(omega_V * t));
        double V0_over_d2 = V_t / (d * d);
        arma::vec Efield = V0_over_d2 * arma::vec({r(0), r(1), -2.0 * r(2)});
        return Efield;
    }

    double V0_over_d2 = V0 / (d * d);
    arma::vec Efield = V0_over_d2 * arma::vec({r(0), r(1), -2.0 * r(2)});
    return Efield;
}

arma::vec PenningTrap::external_B_field(const arma::vec& r, double& r_norm) const {
    if (r_norm > d) {
        return arma::vec({0.0, 0.0, 0.0});  // no field outside the trap
    }
    return arma::vec({0.0, 0.0, B0});
}

// Forces
arma::vec PenningTrap::force_external(int i, double time, double omega_V) const {
    const Particle& p = particles[i];

    double r_norm = arma::vecnorm(p.position);

    arma::vec E = external_E_field(p.position, r_norm, time, omega_V);
    arma::vec B = external_B_field(p.position, r_norm);

    // Lorentz force
    return p.charge * (E + arma::cross(p.velocity, B));
}

arma::vec PenningTrap::force_particle(int i, int j) const {
    const Particle& pi = particles[i];
    const Particle& pj = particles[j];

    arma::vec r_vec = pi.position - pj.position;
    double r_norm = arma::vecnorm(r_vec);
    r_norm = std::max(r_norm, parameters::EPS); // avoid division by zero

    // Coulomb force
    return constants::ke * pi.charge * pj.charge * r_vec / (r_norm * r_norm * r_norm);
}

// returns the acceleration of the particles at a temporary position
arma::mat PenningTrap::acceleration_all(const arma::mat& R, const arma::mat& V, double time, double omega_V)
    const {
    int N = particles.size();
    arma::Mat<double> A(3, N);

    // external acceleration
    arma::rowvec r_norms = arma::vecnorm(R);
    arma::vec r(3);
    arma::vec E(3);
    arma::vec B(3);
    double r_norm;
    for (int i = 0; i < N; i++){
        const Particle& p = particles[i];
        r = R.col(i);
        arma::subview_col<double> v = V.col(i);

        double r_norm = r_norms(i);

        arma::vec E = external_E_field(r, r_norm, time, omega_V);
        arma::vec B = external_B_field(r, r_norm);

        arma::vec F = p.charge * (E + arma::cross(v, B));
        A.col(i) =  F/p.mass;
    }

    // Coulomb forces (symmetric interaction)
    if (coulomb_on && N > 1){
        double inv_r; // 1/r
        double inv_r3;
        arma::vec r_vec(3); // distance vector
        arma::vec a(3); // acceleration vector


        const double inv_mi = 1. / particles[0].mass; // assume constant mass
        const double qi = particles[0].charge; // assume constant
        const double ke_q2 = constants::ke * qi * qi; // useful prefactor

        for (int i = 0; i < N; i++) {
            arma::subview_col<double> col_i = R.col(i); // makes alias
            arma::subview_col<double> ai = A.col(i); // --||--
            for (int j = i + 1; j < N; j++) {
                //const double inv_mj = 1. / particles[j].mass; // removed due equal mass and charge for all particles
                //const double qj = particles[j].charge;

                arma::vec r_vec = col_i - R.col(j);
                
                double r_norm = arma::vecnorm(r_vec);
                r_norm = std::max(r_norm, parameters::EPS); // avoid division by 0
                const double inv_r = 1.0 / r_norm;
                const double inv_r3 = inv_r * inv_r * inv_r;
                
                arma::vec a = ke_q2 * inv_r3 * inv_mi * r_vec; // Coulomb-induced acceleration
                ai += a; // adds to external acceleration
                A.col(j) -= a; // symmetry
            }
        }
    }
return A;
}    

// test functions
int PenningTrap::number_of_particles() {
    int N = particles.size(); 
    int count = 0; // number of particles in trap
    for (int i = 0; i < N; i++) {
        const Particle& p = particles[i];
        double r_norm = arma::norm(p.position);
        if (r_norm < d) {
            count++;
            
        }
        // else {
        //     PenningTrap::delete_particle(i);
        // }
    }
    return count;
}

// Debugging
void PenningTrap::print_particles() const {
    for (const Particle& p : particles) {
        p.print();
    }
}
