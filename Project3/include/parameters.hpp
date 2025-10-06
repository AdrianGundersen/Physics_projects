// parameters.hpp
#pragma once
#include "constants.hpp"

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
// Integration time and step
inline constexpr SimulationParams single{
    500.0,  // total_time_single [µs]
    40000,  // N_single
    true    // coulomb_on_single
};

// FEW PARTICLE PARAMETERS
// Integration time and step
inline constexpr SimulationParams few{
    500.0,  // total_time_few [µs]
    40000,  // N_few
    true    // coulomb_on_few
};

// MULTI PARTICLE PARAMETERS
// Integration time and step
inline constexpr SimulationParams multi{
    500.0,   // total_time_multi [µs]
    5000,   // N_multi
    false   // coulomb_on_multi 
};

// Trap parameters
constexpr int N_particles = 100; // Number of particles
constexpr double maxvel = 40.0; // Maximum initial velocity [µm/µs]

} // namespace parameters
