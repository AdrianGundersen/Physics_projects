// few_particles.cpp
#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "constants.hpp"
#include "integrator.hpp"
#include "parameters.hpp"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>

int main() {
    arma::arma_rng::set_seed(parameters::seed);
    double energy_unitless_to_eV = 1.66054e-27 / 1.602176634e-19; // u (m/s)^2 to eV

    // Extract parameters
    auto sim_params = parameters::few;
    double total_time = sim_params.total_time;
    double dt = sim_params.dt;
    int N = sim_params.N;
    bool coulomb_on = sim_params.coulomb_on;

    std::cout << "\nCoulumb forces: " << coulomb_on <<"\n";

    PenningTrap trap(parameters::B0, parameters::V0, parameters::d, parameters::frequency, coulomb_on);

    // origin particle
    arma::vec pos = {20.0, 0.0, 20.0};
    arma::vec vel = {0.0, 25.0, 0.0};

    Particle p(constants::elementary_charge, constants::Ca_mass, pos, vel);
    trap.add_particle(p);

    // particle to the right
    arma::vec pos2 = {25.0, 25.0, 0.0};
    arma::vec vel2 = {0.0, 40.0, 5.0};

    Particle p2(constants::elementary_charge, constants::Ca_mass, pos2, vel2);


    trap.add_particle(p2);

    //trap.print_particles();

    arma::vec ek = trap.total_energy();
    double tot_ek = arma::sum(ek);
    std::cout << "Total energy: " << tot_ek << " u (m/s)^2 or " << tot_ek * energy_unitless_to_eV << " eV\n";
    //std::cout << "Energies of the particles: "<< ek.t()<< " J" <<"\n";

    double time = 0;

    std::filesystem::create_directory("data");
    Particle& par1 = trap.particles[0];
    Particle& par2 = trap.particles[1];

    std::string filepath1 = "data/" + parameters::filename_few0;
    std::string filepath2 = "data/" + parameters::filename_few1;
    std::string energy_filepath = "data/" + parameters::filename_few_energy;

    std::ofstream ofile1(filepath1);
    std::ofstream ofile2(filepath2);
    std::ofstream energy_ofile(energy_filepath);

    // makes headesr
    ofile1 << "# t x y z vx vy vz" << "\n";
    ofile2 << "# t x y z vx vy vz" << "\n";
    energy_ofile << "# t E_total" << "\n";

    energy_ofile << std::setprecision(17);

    par1.write_to_file(ofile1, time, true);
    par2.write_to_file(ofile2, time, true);
    std::cout << "\nStarting simulation for N=" << N << " over simulation time t = " << total_time << " microseconds \n";
    for (int step = 0; step < N; step++) {
        Integrator::RK4(trap, dt, time);
        time += dt;
        par1.write_to_file(ofile1, time, true);
        par2.write_to_file(ofile2, time, true);
        ek = trap.total_energy(false);
        tot_ek = arma::sum(ek);
        energy_ofile << time << " " << tot_ek << "\n";
    }
    ofile1.close();
    ofile2.close();
    energy_ofile.close();
    //trap.print_particles();

    ek = trap.total_energy();
    tot_ek = arma::sum(ek);



    std::cout << "Total energy: " << tot_ek << " u (m/s)^2 or " << tot_ek * energy_unitless_to_eV << " eV\n";
    //std::cout << "Energies of the particles: "<< ek<< " J" <<"\n";

    std::cout << "Wrote positions against time as: "<< filepath1 << "\n";
    std::cout << "Wrote positions against time as: "<< filepath2 << "\n";
    std::cout << "Wrote total energy against time as: "<< energy_filepath << "\n";
    std::cout << dt << "\n";
    std::cout << time << "\n";
    return 0;
}
