// include/ising/observables.hpp
/*
For observables like <ε>, <|m|> etc.
*/

#pragma once
#include "ising/model.hpp"
#include "ising/lattice.hpp"
namespace ising{
    double total_energy(const Lattice& lat, const Model& model) {} // E
    double energy_per_spin(const Lattice& lat, const Model& model) {} // ε
    double energy2_per_spin(const Lattice& lat, const Model& model) {} // ε²

    double total_magnetization(const Lattice& lat, const Model& model) {} // M
    double magnetization_per_spin(const Lattice& lat, const Model& model) {} // m
    double magnetization2_per_spin(const Lattice& lat, const Model& model) {} // m²
    double heat_capacity_per_spin(const Lattice& lat, const Model& model) {} // C_V
    double susceptibility_per_spin(const Lattice& lat, const Model& model) {} // χ

    struct Observables {
        double E;
        double E2;
        double M;
        double M2;
        double absM;
        double absM2;
        double C_V;
        double chi;
    };
}