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

    constexpr SimulationParams(double t_time, int n_steps, bool coulomb) // calculates dt so not nescessary argument
        : total_time(t_time),
          N(n_steps),
          dt(t_time / static_cast<double>(n_steps)),
          coulomb_on(coulomb) {}
};
 
// ---- GLOBAL SIMULATION PARAMETERS ----
// Seed
constexpr int seed = 67; // group number

// Allowed threads for OpenMP
constexpr int max_threads = 10; // Adjust based on your CPU (will always cap on maximum number though)



// Default Penning trap parameters
constexpr double B0 = constants::Tesla;   // Magnetic field strength [T]
constexpr double V0 = 25.0e-3 * constants::Volt; // Electric potential [V]
constexpr double d  = 500.0;    // Characteristic dimension [µm]
constexpr bool coulomb_on = true; // Enable/disable Coulomb interaction
constexpr double frequency = 0.0; // Frequency of time-dependent potential [MHz]    
constexpr double omega_V = 0.0; // Angular frequency of time-dependent potential [MHz]

//Integration parameters
constexpr double EPS = 1e-12; // Avoid division by zero
constexpr double EPS2 = EPS*EPS; // Used when finding r2

// ----  SIMULATION PARAMETERS ----

// SINGLE PARTICLE PARAMETERS
inline constexpr SimulationParams single{
    50.0,  // total_time_single [µs]
    32000,  // N_single
    false    // coulomb_on_single (only one particle)
};


// filename for output data (if using another naming scheme, update plotting accordingly (not recommended))
const std::string filename_single = "single_particle_N" + std::to_string(single.N) + ".txt"; // particle 0



// FEW PARTICLE PARAMETERS
inline constexpr SimulationParams few{
    50.0,  // total_time_few [µs]
    10000,  // N_few
    false    // coulomb_on_few
};

// filename for output data (if using another naming scheme, update plotting accordingly (not recommended))
const std::string filename_few0 = "pos_vel_0_coulomb=" + std::to_string(few.coulomb_on) + "_N" + std::to_string(few.N) + ".txt"; // particle 0
const std::string filename_few1 = "pos_vel_1_coulomb=" + std::to_string(few.coulomb_on) + "_N" + std::to_string(few.N) + ".txt"; // particle 1

// MULTI PARTICLE PARAMETERS
inline constexpr SimulationParams multi{
    500.0,   // total_time_multi [µs]
    40000,   // N_multi
    false   // coulomb_on_multi coulumb forces
};

// Trap setup
constexpr int N_particles = 100; // Number of particles 
constexpr double pos_scaling = 0.1; // Position scaling factor so typical is pos_scaling * d 
constexpr double vel_scaling = 0.1; // Velocity scaling factor so typical is vel_scaling * d

// Frequency and amplitude parameters for multi-particle simulations
inline const arma::vec f_list = {0.1, 0.4, 0.7}; // Amplitude factors

// omega_V list in MHz
constexpr double w_min = 0.2, w_max = 2.5, w_step = 0.005; // omega_V range and step
constexpr int n_omega = static_cast<int>((w_max - w_min) / w_step + 1.5); // # of omega_points (+1.5 to avoid rounding issues)
inline const arma::vec omega_V_list = arma::linspace(w_min, w_max, n_omega);

// filename for output data
// recommended:

const std::string filename_multi = 
    std::string("trapped_w") + std::to_string(w_min) + "-" + std::to_string(w_max) +
    "_dw" + std::to_string(w_step) + "_N" + std::to_string(multi.N) + ".txt";


} // namespace parameters
