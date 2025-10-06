// integrator.hpp
#pragma once
#define INTEGRATOR_HPP
#include "PenningTrap.hpp"
#include "constants.hpp"

class Integrator {
    public:
    static void ForwardEuler(PenningTrap& trap, double dt, double time = 0.0);
    static void RK4(PenningTrap& trap, double dt, double time = 0.0);
};