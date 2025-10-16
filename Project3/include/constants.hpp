// constants.hpp
#pragma once

namespace constants {
// Base units
constexpr double micrometer = 1.0;  
constexpr double microsecond = 1.0;
constexpr double atomic_mass_unit = 1.0;
constexpr double elementary_charge = 1.0;

// Coulomb constant
constexpr double ke = 1.38935333e5;

// Derived units
constexpr double Tesla = 9.64852558e1;   
constexpr double Volt  = 9.64852558e7; 
constexpr double Joule = atomic_mass_unit * (micrometer * micrometer) / (microsecond * microsecond);


// Particle properties
constexpr double Ca_mass = 40.078 * atomic_mass_unit; // Calcium ion mass

}