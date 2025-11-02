// include/ising/observables.hpp
/*
For observables like <ε>, <|m|> etc.
*/

#pragma once
#include "ising/lattice.hpp"
#include "ising/model.hpp"

namespace ising{
    double total_energy(const Lattice& lat, const Model& model);
    double energy_per_spin(const Lattice& lat, const Model& model);

    double total_magnetization(const Lattice& lat);
    double magnetization_per_spin(const Lattice& lat);

    double heat_capacity(const Lattice& lat, const Model& model, double T);
    double susceptibility(const Lattice& lat, const Model& model, double T);
    double average(const std::vector<double>& data);
    double heat_capacity(const Lattice& Lat, double eps2_avg, double eps_avg, double T);
    double susceptibility(const Lattice& Lat, double m2_avg, double m_avg, double T);
    struct Observables {
        double E{0.0};
        int   M{0};
        double absM{0.0};
    };
}