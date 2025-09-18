#include "jacobi.hpp"
#include "tridiag.hpp"
#include "analytical.hpp"
#include <cmath>
#include <iostream>

// problems/problem6.cpp
#include "jacobi.hpp"
#include "tridiag.hpp"
#include <armadillo>
#include <iostream>
#include <filesystem>
#include <fstream>

int main(){
    int N = 100 - 1;
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
    
    arma::mat R_analytical(N,N);

    for (int j = 1; j <= N; ++j){ // j = 1, 2, ..., N
        arma::vec v_analytical = analytical_eigenvector(N, j);
        R_analytical.col(j) = v_analytical;
    }
    
    //write to file
    namespace fs = std::filesystem;
    fs::create_directories("output");

    arma::mat R_save = R.cols(0,2);
    eigenvalues.save("output/problem6_eigenvalues100.csv", arma::csv_ascii);
    R_save.cols(0,2).eval().save("output/problem6_eigenvectors100.csv", arma::csv_ascii);

    
    return 0;
    
}