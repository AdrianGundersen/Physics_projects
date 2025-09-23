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
    for (int n : {1,10}){
        int N = 10*n - 1;
        double h = 1./(N+1); // step length
        double a = -1./(h*h); // sub- and superdiagonal elements
        double d = 2./(h*h);   // diagonal elements
        //maximum iterations
        int maxiter = 1e8;
        int iterations = 0;

        // lowest value vanted
        double eps = 1e-8;

        //argument to se in we converged
        bool converged = false;

        arma::mat R(N, N, arma::fill::eye);
        arma::Mat<double> A = create_tridiagonal(N, a, d);
        
        arma::vec eigenvalues;
        jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);

        
        for (int i = 0; i < 3 && i < eigenvalues.n_elem; ++i){
            std::cout << "eigenvalue:\n" << eigenvalues(i)
            << "\n\neigenvector:\n" << R.col(i) << "\n";
        }
        
        arma::mat R_analytical(N,N);

        for (int j = 1; j <= N; ++j){ // j = 1, 2, ..., N
            arma::vec v_analytical = analytical_eigenvector(N, j);
            v_analytical /= arma::norm(v_analytical, 2); // normalize
            R_analytical.col(j-1) = v_analytical;
        }

        //write to file
        namespace fs = std::filesystem;
        fs::create_directories("output");

        arma::mat R_save = R.cols(0,2);

        // We have n=N+1  in the filename to distinguish
        std::string filename_eigvals = "output/problem6_eigenvalues" + std::to_string(N+1) + ".csv";
        std::string filename_eigvec = "output/problem6_eigenvectors" + std::to_string(N+1) + ".csv";
        std::string filename_analytical = "output/problem6_eigenvectors_analytical" + std::to_string(N+1) + ".csv";

        eigenvalues.save(filename_eigvals, arma::csv_ascii);
        R_save.cols(0,2).eval().save(filename_eigvec, arma::csv_ascii);
        R_analytical.cols(0,2).eval().save(filename_analytical, arma::csv_ascii);

        
    }
    return 0;
}