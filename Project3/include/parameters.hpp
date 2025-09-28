// parameters.hpp
#pragma once
#define PARAMETERS_HPP
#include "constants.hpp"

namespace parameters {
// Default Penning trap parameters
constexpr double B0 = constants::Tesla;   
constexpr double V0 = 2.41e6; 
constexpr double d  = 500.0;    

//Integration parameters
constexpr double total_time = 100.0; // [µs]
constexpr int N = 1000; // Number of integration steps
constexpr double dt = total_time / N; // Time step
}