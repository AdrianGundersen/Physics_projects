// integrator.hpp
#pragma once
#include "PenningTrap.hpp"
#include "constants.hpp"

class Integrator {
    public:
    static void ForwardEuler(PenningTrap& trap, double dt, double time = 0.0, double omega_V = 0.0);
    static void RK4(PenningTrap& trap, double dt, double time = 0.0, double omega_V = 0.0);
};