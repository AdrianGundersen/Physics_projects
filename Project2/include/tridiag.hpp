#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#include <armadillo>

arma::Mat<double> create_tridiagonal(int N, double a, double d);
double max_offdiag_symmetric(const arma::mat& A, int& k, int &l);
#endif