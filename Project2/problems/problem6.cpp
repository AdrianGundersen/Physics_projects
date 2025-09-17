#include "jacobi.hpp"
#include "tridiag.hpp"
#include <cmath>
#include <iostream>

// problems/problem6.cpp
#include "jacobi.hpp"
#include "tridiag.hpp"
#include <armadillo>
#include <iostream>

int main(){
    int N = 11;
    double a = 2.0;
    double d = -1.0;
    //maximum iterations
    int maxiter = 1e8;
    int iterations = 0;

    // lowest value vanted
    double eps = 1e-8;

    //argument to se in we converged
    bool converged = false;

    arma::mat R(N, N, arma::fill::eye);
    arma::Mat<double> A = create_tridiagonal(N, d, a);
    
    arma::vec eigenvalues;
    jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);

    for (int i = 0; i < 3 && i < eigenvalues.n_elem; ++i){
        std::cout << "eigenvalue:\n" << eigenvalues(i)
                  << "\n\neigenvector:\n" << R.col(i) << "\n";
    }
    return 0;
    
}