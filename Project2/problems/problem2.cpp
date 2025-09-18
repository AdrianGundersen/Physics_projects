#include "tridiag.hpp"
#include "analytical.hpp"

#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>

int main(){
    int N = 6; // matrix size
    double h = 1./(N+1); // step length
    double a = -1./(h*h); // sub- and superdiagonal elements
    double d = 2./(h*h);   // diagonal elements
    arma::Mat<double> A = create_tridiagonal(N, a, d);

    // Numerical eigenvalues and eigenvectors
    arma::vec eigvals_num;
    arma::mat eigvecs_num;
    arma::eig_sym(eigvals_num, eigvecs_num, A); // calculates eigenvalues and eigenvectors

    A.print("Matrix A: "); 
    eigvals_num.print("Eigenvalues: ");

    for (int j = 1; j <= N; ++j){ // j = 1, 2, ..., N
        arma::vec v_analytical = analytical_eigenvector(N, j);
        v_analytical /= arma::norm(v_analytical, 2); // normalize

        arma::vec v_num = eigvecs_num.col(j-1); // j-th eigenvector 
        v_num /= arma::norm(v_num, 2); // normalize
        if (arma::dot(v_analytical, v_num) < 0){
            v_num = -v_num; // adjust sign if necessary
        };

        // Eigenvalues
        double lambda_analytical = analytical_eigenvalue(N, a, d, j);
        double lambda_num  = eigvals_num(j-1);

        // Differences
        double eigval_diff = std::abs(lambda_analytical - lambda_num);
        double eigvec_diff = arma::norm(v_analytical - v_num, 2);

        // Output
        std::cout << "j: " << j << "\n";
        std::cout << "Difference in eigenvalue: " << eigval_diff << "\n";
        std::cout << "Difference in eigenvector: " << eigvec_diff << "\n";
    }
    return 0;
}