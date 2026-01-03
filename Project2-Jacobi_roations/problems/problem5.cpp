
#include "jacobi.hpp"
#include "tridiag.hpp"
#include <armadillo>
#include <iostream>
#include <fstream>

int main(){
    int N = 6;
    double h = 1./(N+1); // step length
    double a = -1./(h*h); // sub- and superdiagonal elements
    double d = 2./(h*h);   // diagonal elements
    //maximum iterations
    int maxiter = 1e8;

    // lowest value vanted
    double eps = 1e-8;

    //argument to se in we converged
    bool converged = false;

    //Tridaiagonal matrix
    arma::Mat<double> A = create_tridiagonal(N, a, d);

    //write to file
    std::ofstream outfile("output/iterations.txt");

    //Run jacobi-eigensolver for different values of N to registrer iterations
    for (int i = 1; i <= 20; ++i){
        int N = i*5;

        //Tridiagonal matrix
        int iterations = 0;
        arma::mat R(N, N);
        arma::Mat<double> A = create_tridiagonal(N, a, d);
        arma::vec eigenvalues; 
        jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);

        //Dense random matrix
        int iterations_random = 0;
        arma::mat R_random(N,N, arma::fill::eye);
        arma::mat A_random = arma::mat(N,N).randn();
        A_random = arma::symmatu(A_random);
        arma::vec eigenvalues_random;
        jacobi_eigensolver(A_random, eps, eigenvalues_random, R_random, maxiter, iterations_random, converged);

        std::cout << N  << "\t" << iterations << "\t" << iterations_random << "\n";
        outfile   << N  << "\t" << iterations << "\t" << iterations_random << "\n";
        }
    outfile.close();
    return 0;
}