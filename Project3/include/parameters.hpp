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
constexpr double total_time = 1050.0; // [µs]
constexpr int N = 100000; // Number of integration steps minimum value: 10000
constexpr double dt = total_time / N; // Time step

constexpr double EPS = 1e-12; // Avoid division by zero
}