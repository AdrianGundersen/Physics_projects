// PenningTrap.cpp
#include "PenningTrap.hpp"
#include "Particle.hpp"
#include "constants.hpp"
#include "parameters.hpp"
#include <cmath>
#include <algorithm> 
#include <vector>

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f_in, bool coulomb_on_in, double omega_V_in)
// Constructor
    : B0(B0_in), V0(V0_in), d(d_in), f(f_in), coulomb_on(coulomb_on_in), omega_V(omega_V_in) {} 

// Add a particle
void PenningTrap::add_particle(const Particle& p) {
    particles.push_back(p); // add a copy of p to the particles vector
}

void PenningTrap::delete_particle(int particle_i) {
    particles.erase(particles.begin() + particle_i);
}

// Fill the trap with random particles
void PenningTrap::fill_random(int N, double q, double m, double pos_scaling, double vel_scaling) {
    for (int i = 0; i < N; i++) {
        arma::vec pos = arma::vec(3).randn() * pos_scaling * d; // typical pos is pos_scaling * d
        arma::vec vel = arma::vec(3).randn() * vel_scaling * d; // typical vel is vel_scaling * d / microsecond
        Particle p(q, m, pos, vel);
        add_particle(p);
    }
}


// External fields
arma::vec PenningTrap::external_E_field(const arma::vec& r, double& r_norm, double t) const {
    if (r_norm > d) {
        return arma::vec({0.0, 0.0, 0.0});  // no field outside the trap
    }
    if (f != 0.0) {
        double V_t = V0 * (1 + f * std::cos(omega_V * t));
        double Vt_over_d2 = V_t / (d * d);
        arma::vec Efield = Vt_over_d2 * arma::vec({r(0), r(1), -2.0 * r(2)});
        return Efield;
    }

    double V0_over_d2 = V0 / (d * d);
    arma::vec Efield = V0_over_d2 * arma::vec({r(0), r(1), -2.0 * r(2)});
    return Efield;
}

// Magnetic field
arma::vec PenningTrap::external_B_field(const arma::vec& r, double& r_norm) const {
    if (r_norm > d) {
        return arma::vec({0.0, 0.0, 0.0});  // no field outside the trap
    }
    return arma::vec({0.0, 0.0, B0});
}

// Forces on the particle from the trtap
arma::vec PenningTrap::force_external(int i, double time) const {
    const Particle& p = particles[i];

    double r_norm = arma::vecnorm(p.position);

    arma::vec E = external_E_field(p.position, r_norm, time);
    arma::vec B = external_B_field(p.position, r_norm);

    // Lorentz force
    return p.charge * (E + arma::cross(p.velocity, B));
}

// Force on particle i from particle j, Not used since integrated into acceleration all
arma::vec PenningTrap::force_particle(int i, int j) const {
    const Particle& pi = particles[i];
    const Particle& pj = particles[j];

    arma::vec r_vec = pi.position - pj.position; // distance vector
    double r_norm = arma::vecnorm(r_vec);
    r_norm = std::max(r_norm, parameters::EPS); // avoid division by zero

    // Coulomb force
    return constants::ke * pi.charge * pj.charge * r_vec / (r_norm * r_norm * r_norm);
}

// returns the acceleration of the particles at a temporary position
arma::mat PenningTrap::acceleration_all(const arma::mat& R, const arma::mat& V, double time)
    const {
    int N = particles.size();
    arma::Mat<double> A(3, N);

    // external acceleration
    arma::rowvec r_norms = arma::vecnorm(R);
    arma::vec r(3);
    arma::vec E(3);
    arma::vec B(3);
    double r_norm;
    // loop over particles
    for (int i = 0; i < N; i++){
        const Particle& p = particles[i];
        arma::subview_col<double> r = R.col(i);
        arma::subview_col<double> v = V.col(i); // alias

        double r_norm = r_norms(i); 

        arma::vec E = external_E_field(r, r_norm, time);
        arma::vec B = external_B_field(r, r_norm);

        arma::vec F = p.charge * (E + arma::cross(v, B)); // Lorentz force
        A.col(i) =  F/p.mass; // accleration 
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

int PenningTrap::number_of_particles() {
    // also removes particles outside the trap (should maybe rename or be optional)
    std::erase_if(particles, [outside = d](const Particle& p){ // remove if outside the trap requires C++20
        return arma::norm(p.position) >= outside;  // all particles outside the trap
    });
    return static_cast<int>(particles.size());  // remaining inside 
}

// Calculatng total energy of the particles
arma::vec PenningTrap::total_energy() const {
    const int N = particles.size();

    //constants 
    const double V0_over_d2 = (V0) / (d * d); // v0 / d^2
    const double k_e = constants::ke; //8.9875517923e9;   // 1/4pi * epsilon_0

    arma::vec kinetic_energies(N);      // Kinetic energy
    arma::vec EPotential(N);            // Electical potential
    arma::vec coulumb_potential(N);     // Coulumb potential
    

    for (int i = 0; i<N; i++) {
        const Particle& p = particles[i];
        arma::vec r = p.position;       // pos in micrometer
        arma::vec v = p.velocity;   // vel in micrometer/microsecond
        double particle_charge = p.charge;  // charge in e
        double particle_mass = p.mass;      // mass in u

        kinetic_energies(i) = 0.5 * particle_mass * arma::dot(v, v);
        EPotential(i) =  0.5 * particle_charge * V0_over_d2 * (2 * r(2)*r(2) - r(0)*r(0) - r(1)*r(1)); // (1/2) * q * V0/d^2 * (2z^2 - x^2 - y^2)
    }

    //coulumb potential
    double coulumb;
    if (coulomb_on && N > 1){
        double r_norm;
        const double qi = particles[0].charge; // assume qj = qi constant
        const double ke_q2 = k_e * qi * qi; // 1/(4pi * epsilon_0) * q_i*q_j

        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                // calling two particles
                const Particle& pi = particles[i];
                const Particle& pj = particles[j];

                arma::vec r_vec = (pi.position - pj.position);   // distance vector
                double r_norm = arma::norm(r_vec, 2);     // distance norm
                r_norm = std::max(r_norm, parameters::EPS); // avoid division by 0

                double U = ke_q2 / r_norm;  //coulumb potential

                // giving each particle half the potential
                coulumb_potential(i) += 0.5 * U;
                coulumb_potential(j) += 0.5 * U;

        }}}
    std::cout   << "Kinetic:  " << arma::sum(kinetic_energies) << " u (m/s)^2" << "\n" 
                << "Electric: " << arma::sum(EPotential) << " u (m/s)^2" << "\n" 
                << "Coulumb:  " << arma::sum(coulumb_potential)<< " u (m/s)^2" << "\n";
    arma::vec total_energy = kinetic_energies + EPotential + coulumb_potential;

    return total_energy;
}

// Debugging
void PenningTrap::print_particles() const {
    for (const Particle& p : particles) {
        p.print();
    }
}
