#include "jacobi.hpp"
#include "tridiag.hpp"
#include "analytical.hpp"
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

    arma::Mat<double> A = create_tridiagonal(N, a, d);
    
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    arma::eig_sym(eigenvalues, eigenvectors, A);

    arma::mat R_analytical(N,N);
    arma::vec eigenvalue_anal(N);

    for (int j = 1; j <= N; ++j){ // j = 1, 2, ..., N
        arma::vec v_analytical = analytical_eigenvector(N, j);
        v_analytical /= arma::norm(v_analytical, 2); // normalize
        R_analytical.col(j-1) = v_analytical;
        eigenvalue_anal(j-1) = analytical_eigenvalue(N, a, d, j);
        }
    

    //sort analytical eivenvalues
    arma::uvec order_anal = arma::sort_index(eigenvalue_anal); //uvec hold unsigned ints
    eigenvalue_anal = eigenvalue_anal(order_anal);

    R_analytical = R_analytical.cols(order_anal);

    jacobi_eigensolver(A, eps, eigenvalues, R, maxiter, iterations, converged);
    
    std::cout << "Amount of iterations : " << iterations << "\n";
    
    eigenvalues.print("Eigenvalues : ");
    eigenvalue_anal.print("Analytical eigenvalues");

    eigenvectors.print("Eigenvectors : ");
    R_analytical.print("Analytical eigenvectors");

}