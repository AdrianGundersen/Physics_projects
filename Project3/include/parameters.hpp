// parameters.hpp
#pragma once
#include "constants.hpp"
#include <armadillo>
#include <cmath>

namespace parameters {

struct SimulationParams {
    double total_time; // Total simulation time [µs]
    int N;             // Number of integration steps
    double dt;         // Time step [µs]
    bool coulomb_on;

    constexpr SimulationParams(double t_time, int n_steps, bool coulomb)
        : total_time(t_time),
          N(n_steps),
          dt(t_time / static_cast<double>(n_steps)),
          coulomb_on(coulomb) {}
};
 
// ---- GLOBAL SIMULATION PARAMETERS ----
// Seed
constexpr int seed = 67; // group number

// Allowed threads for OpenMP
constexpr int max_threads = 16; // Adjust based on your CPU



// Default Penning trap parameters
constexpr double B0 = constants::Tesla;   // Magnetic field strength [T]
constexpr double V0 = 2.41e6; // Electric potential [V]
constexpr double d  = 500.0;    // Characteristic dimension [µm]
constexpr bool coulomb_on = true; // Enable/disable Coulomb interaction
constexpr double frequency = 0.0; // Frequency of time-dependent potential [MHz]    
constexpr double omega_V = 0.0; // Angular frequency of time-dependent potential [MHz]


//Integration parameters
constexpr double EPS = 1e-12; // Avoid division by zero

// ----  SIMULATION PARAMETERS ----

// SINGLE PARTICLE PARAMETERS
inline constexpr SimulationParams single{
    500.0,  // total_time_single [µs]
    40000,  // N_single
    true    // coulomb_on_single
};

// FEW PARTICLE PARAMETERS
inline constexpr SimulationParams few{
    500.0,  // total_time_few [µs]
    40000,  // N_few
    true    // coulomb_on_few
};

// MULTI PARTICLE PARAMETERS
inline constexpr SimulationParams multi{
    500.0,   // total_time_multi [µs]
    40000,   // N_multi
    false   // coulomb_on_multi 
};

// Trap setup
constexpr int N_particles = 100; // Number of particles 
constexpr double pos_scaling = 0.1; // Position scaling factor so typical is pos_scaling * d 
constexpr double vel_scaling = 0.1; // Velocity scaling factor so typical is vel_scaling * d

// Frequency and amplitude parameters for multi-particle simulations
inline const arma::vec f_list = {0.1, 0.4, 0.7}; // Amplitude factors

// omega_V list in MHz
constexpr double w_min = 0.20, w_max = 2.50, w_step = 0.005;
constexpr int n_omega = static_cast<int>((w_max - w_min) / w_step + 1.5); // # of omega_points (+1.5 to avoid rounding issues)
inline const arma::vec omega_V_list = arma::linspace(w_min, w_max, n_omega);


} // namespace parameters
