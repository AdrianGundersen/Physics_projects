#include "integrator.hpp"
#include <armadillo>

void Integrator::ForwardEuler(PenningTrap& trap, double dt, double time) {
    int N = trap.particles.size(); // will return correct number even if particles are removed
    arma::vec new_position(3);
    arma::vec new_velocity(3);
    arma::mat r0(3, N);
    arma::mat v0(3, N);
    for (int i = 0; i < N; i++) {
        Particle& p = trap.particles[i];
        v0.col(i) = p.velocity;
        r0.col(i) = p.position;

    arma::mat a0 = trap.acceleration_all(r0, v0, time); // 

    arma::mat r1 = r0 + v0 * dt;
    arma::mat v1 = v0 + a0 * dt;

    for (int i = 0; i < N; i++) {
        trap.particles[i].position = r1.col(i);
        trap.particles[i].velocity = v1.col(i);
    }
}
}

void Integrator::RK4(PenningTrap& trap, double dt, double time) {
    int N = trap.particles.size(); // will return correct number even if particles are removed

    // temporary copy of current state
    arma::mat r0(3, N);
    arma::mat v0(3, N);
    for (int i = 0; i < N; i++){
        v0.col(i) = trap.particles[i].velocity;
        r0.col(i) = trap.particles[i].position;
    }

    // k-vectors containing velocity and position
    arma::mat k_r1(3, N), k_r2(3, N), k_r3(3, N), k_r4(3, N);
    arma::mat k_v1(3, N), k_v2(3, N), k_v3(3, N), k_v4(3, N);
    

    arma::mat a1 = trap.acceleration_all(r0, v0, time);
    k_r1 = v0 * dt;
    k_v1 = a1 * dt;

    arma::mat r2 = r0 + 0.5 * k_r1;
    arma::mat v2 = v0 + 0.5 * k_v1;
    arma::mat a2 = trap.acceleration_all(r2, v2, time + 0.5 * dt);
    k_r2 = v2 * dt;
    k_v2 = a2 * dt;

    arma::mat r3 = r0 + 0.5 * k_r2;
    arma::mat v3 = v0 + 0.5 * k_v2;
    arma::mat a3 = trap.acceleration_all(r3, v3, time + 0.5 * dt);
    k_r3 = v3 * dt;
    k_v3 = a3 * dt;

    arma::mat r4 = r0 + k_r3;
    arma::mat v4 = v0 + k_v3;
    arma::mat a4 = trap.acceleration_all(r4, v4, time + dt);
    k_r4 = v4 * dt;
    k_v4 = a4 * dt;

    arma::mat new_position = r0 + (k_r1 + 2*k_r2 + 2*k_r3 + k_r4)/6.0;
    arma::mat new_velocity = v0 + (k_v1 + 2*k_v2 + 2*k_v3 + k_v4)/6.0;


    // updates all particles
    for (int i = 0; i < N; i++){
        trap.particles[i].position = new_position.col(i);
        trap.particles[i].velocity = new_velocity.col(i);
    }
}
