#ifndef ANALYTICAL_HPP
#define ANALYTICAL_HPP

#include <armadillo>
#include <cmath>

//calculating the analythical eigenvectors and eigenvalues
// N is size of matrix
// a is diagonal element [0,-1] and [1, 0]
// d is diagonal element [0,0]
// j is eigenvalue or eigenvector number j
double analytical_eigenvalue(int N, double a, double d, double j);
arma::vec analytical_eigenvector(int N, int j);
#endif