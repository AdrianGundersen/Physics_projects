// parameters.hpp
#pragma once
#define PARAMETERS_HPP
#include "constants.hpp"

namespace parameters {
// Seed

constexpr int seed = 67; // group number

// Default Penning trap parameters
constexpr double B0 = constants::Tesla;   // Magnetic field strength [T]
constexpr double V0 = 2.41e6; // Electric potential [V]
constexpr double d  = 500.0;    // Characteristic dimension [µm]
constexpr bool coulomb_on = true; // Enable/disable Coulomb interaction

//Integration parameters
constexpr double total_time = 50.0; // [µs]
constexpr int N = 100000; // Number of integration steps minimum value: 10000
constexpr double dt = total_time / N; // Time step

constexpr double EPS = 1e-12; // Avoid division by zero

// Multiparticle parameters
constexpr int N_particles = 100; // Number of particles
constexpr double maxvel = 40.0; // Maximum initial velocity [µm/µs]
}
