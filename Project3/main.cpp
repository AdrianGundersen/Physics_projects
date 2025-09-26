#include "Particle.hpp"
#include "PenningTrap.hpp"

int main() {
    arma::vec pos = {0.0, 0.0, 0.0};
    arma::vec vel = {0.0, 0.0, 0.0};

    Particle p(1.0, 1.0, pos, vel);

    p.print();

    return 0;
}
