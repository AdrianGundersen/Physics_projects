#pragma once
#include <armadillo>

class Particle {
public:
    double charge;
    double mass;
    arma::vec position;  
    arma::vec velocity;  

    Particle(double charge, double mass, 
             const arma::vec& position, 
             const arma::vec& velocity); // constructor
    void print() const; // print particle info
};
