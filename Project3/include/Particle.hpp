#pragma once
#include <armadillo>

class PenningTrap; 

class Particle {
    friend class PenningTrap; // PenningTrap can access private

private:
    double charge;
    double mass;
    arma::vec position;  
    arma::vec velocity;  

    public:
    Particle(double charge, double mass, 
             const arma::vec& position, 
             const arma::vec& velocity); // constructor
    void print() const; // print particle info
};
