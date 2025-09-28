#include "integrator.hpp"
#include <armadillo>

void Integrator::ForwardEuler(PenningTrap& trap, double dt) {
    int N = trap.particles.size();
    arma::vec new_position(3);
    arma::vec new_velocity(3);

    for (int i = 0; i < N; i++) {
        Particle& p = trap.particles[i];
        arma::vec F = trap.total_force(i);
        arma::vec acceleration = F / p.mass;

        // Update position and velocity using Forward Euler
        new_velocity = p.velocity + acceleration * dt;
        new_position = p.position + p.velocity * dt;

        p.velocity = new_velocity;
        p.position = new_position;
    }
};