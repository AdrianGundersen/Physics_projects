#include "jacobi.hpp"
#include "tridiag.hpp"
#include <cmath>
#include <iostream>


int main(){
    int N = 6;
    double a = 2.0;
    double d = -1.0;
    int k = 0;
    int l = 0;
    arma::mat R(N, N);
    R.eye();
    //maximum iterations
    int maxiter = 1e8;
    int iterations = 0;

    // lowest value vanted
    double eps = 1e-8;

    //argument to se in we converged
    bool converged = false;

    arma::Mat<double> A = create_tridiagonal(N, d, a);
    
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, A);
    arma::vec eigenvalues_copy = eigenvalues;
    jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);
    
    // eigenvalues.print("Eigenvectors : ");
    // eigenvalues_copy.print("Armadillo sin egenvektor");


    // Problem 5
    // for(int i = 5; i <= 100; i += 1){
    //     iterations = 0;
    //     arma::mat R(i, i);
    //     R.eye();
    //     arma::Mat<double> A = create_tridiagonal(i, d, a);
    //     jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);
    //     std::cout << i << "\t" << iterations << "\n";
    // }

    //roblem 6
    int N_new = 11; // N = n+1
    arma::mat R_new(N_new, N_new);
    R_new.eye();
    arma::Mat<double> A_new = create_tridiagonal(N_new, d, a);
    jacobi_eigensolver(A_new, eps, eigenvalues, R_new, maxiter, iterations, converged);

    for(int i =0; i < 3; i++){
    std::cout << "eigenvalue: \n" << eigenvalues(i) << "\n \n Eigenvector: \n" <<  R_new.col(i) << "\n";
    }
    
}