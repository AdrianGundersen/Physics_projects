
#include "jacobi.hpp"
#include "tridiag.hpp"
#include <armadillo>
#include <iostream>

int main(){
    int N = 6;
    double a = 2.0;
    double d = -1.0;
    //maximum iterations
    int maxiter = 1e8;
    int iterations = 0;

    // lowest value vanted
    double eps = 1e-8;

    //argument to se in we converged
    bool converged = false;

    arma::Mat<double> A = create_tridiagonal(N, d, a);

    for (int N = 5; N <= 100; ++N){
        iterations = 0;
        arma::mat R(N, N);
        arma::Mat<double> A = create_tridiagonal(N, d, a);
        arma::vec eigenvalues; 
        jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);
        std::cout << N  << "\t" << iterations << "\n";
        }
    return 0;
}