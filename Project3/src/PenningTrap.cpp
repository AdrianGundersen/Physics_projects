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


void PenningTrap::fill_random(int N, double q, double m, double pos_scaling, double vel_scaling) {
    for (int i = 0; i < N; i++) {
        arma::vec pos = arma::vec(3).randn() * pos_scaling * d; // typical pos is pos_scaling * d
        arma::vec vel = arma::vec(3).randn() * vel_scaling * d; // typical vel is vel_scaling * d / microsecond
        Particle p(q, m, pos, vel);
        add_particle(p);
    }
}


// External fields
arma::vec PenningTrap::external_E_field(const arma::vec& r, double t, double omega_V) const {
    if (arma::norm(r) > d) {
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

arma::vec PenningTrap::external_B_field(const arma::vec& r) const {
    if (arma::norm(r) > d) {
        return arma::vec({0.0, 0.0, 0.0});  // no field outside the trap
    }
    return arma::vec({0.0, 0.0, B0});
}

// Forces
arma::vec PenningTrap::force_external(int i, double time, double omega_V) const {
    const Particle& p = particles[i];

    arma::vec E = external_E_field(p.position, time, omega_V);
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

arma::vec PenningTrap::total_force(int i, double time, double omega_V) const {
    arma::vec F = force_external(i, time, omega_V);

    for (int j = 0; j < (int)particles.size(); j++) {
        if (j != i) {
            F += force_particle(i, j);
        }
    }
    return F;

}

// returns the acceleration of the particles at a temporary position
arma::mat PenningTrap::acceleration_all(const arma::mat& R, const arma::mat& V, double time, double omega_V)
    const {
    int N = particles.size();
    arma::Mat<double> A(3, N);

    // external acceleration
    for (int i = 0; i < N; i++){
        const Particle& p = particles[i];
        arma::vec r = R.col(i);
        arma::vec v = V.col(i);

        arma::vec E = external_E_field(r, time, omega_V);
        arma::vec B = external_B_field(r);
        arma::vec F = p.charge * (E + arma::cross(v, B));
        A.col(i) =  F/p.mass;
    }

    // Coulomb forces (symmetric interaction)
    if (coulomb_on && N > 1) {
        const double inv_mi = 1. / particles[0].mass; // assume costante mass
        const double qi = particles[0].charge; // assume costante charge
        const double ke_q2 = constants::ke * qi*qi;

        for (int i = 0; i < N; i++) {
            const arma::vec col_i = R.col(i);
            arma::vec  ai = A.col(i);
            for (int j = i + 1; j < N; j++) {
                // const double inv_mj = 1. / particles[j].mass; // removed due to being equal to the others
                // const double qj = particles[j].charge;

                arma::vec r_vec = col_i - R.col(j);
            
                double r2    = arma::dot(r_vec, r_vec) + parameters::EPS2;    // to save FLOPs, but might cause some more error 
                double inv_r = 1.0 / std::sqrt(r2);
                double inv_r3 = inv_r * inv_r * inv_r;    
                
                arma::vec a = ke_q2 * inv_r3 * inv_mi * r_vec;
                ai += a;
                A.col(j) -= a;
            }
        }
}
return A;
}   

// test functions
int PenningTrap::number_of_particles() const {
    int N = particles.size(); 
    int count = 0; // number of particles in trap
    for (int i = 0; i < N; i++) {
        const Particle& p = particles[i];
        double r_norm = arma::norm(p.position);
        if (r_norm < d) {
            count++;
        }
    }
    return count;
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