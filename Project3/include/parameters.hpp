// parameters.hpp
#pragma once
#define PARAMETERS_HPP
#include "constants.hpp"

namespace parameters {
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
constexpr double total_time_single = 500.0; // [µs]
constexpr int N_single = 40000; // Number of integration steps minimum value: 10000
constexpr double dt_single = total_time_single / N_single; // Time step


// FEW PARTICLE PARAMETERS
// Integration time and step
constexpr double total_time_few = 500.0; // [µs]
constexpr int N_few = 40000; // Number of integration steps minimum value: 10000
constexpr double dt_few = total_time_few / N_few; // Time step


// MULTI PARTICLE PARAMETERS
// Integration time and step
constexpr double total_time_multi = 50.0; // [µs]
constexpr int N_multi = 5000; // Number of integration steps minimum value: 100
constexpr double dt_multi = total_time_multi / N_multi; // Time step

// Trap parameters
constexpr int N_particles = 100; // Number of particles
constexpr double maxvel = 40.0; // Maximum initial velocity [µm/µs]

}
