#include "jacobi.hpp"



// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    int N = A.n_rows;
    // trig varibles
    double tau = (A(l ,l) - A(k, k)) / (2*A(k, l));
    double tan_theta;
    if (tau > 0){
        tan_theta = 1/(tau + std::sqrt(1+std::pow(tau, 2)));
    }
    else{
        tan_theta = 1/(tau - std::sqrt(1+std::pow(tau, 2)));
    }
    double cos_theta = 1/std::sqrt(1+std::pow(tan_theta, 2));
    double sin_theta = cos_theta * tan_theta;
    

    // Flop and gather saves 
    double cos_square = std::pow(cos_theta, 2);
    double sin_square = std::pow(sin_theta, 2); 
    double cos_sin = cos_theta * sin_theta; 

    double A_kk = A(k, k);
    double A_kl = A(k, l);
    double A_ll = A(l, l);

    // Update A
    A(k, k) = A_kk * cos_square - 2*A_kl * cos_sin + A_ll * sin_square;
    A(l, l) = A_ll * cos_square + 2*A_kl * cos_sin + A_kk * sin_square;
    A(k, l) = A(l, k) = 0.0;
    for(int i = 0; i < N; i++){
        if (i != k && i != l){
            double A_ik = A(i, k);
            double A_il = A(i, l);
            A(i, k) = A(k, i) = A_ik * cos_theta - A_il * sin_theta;
            A(i, l) = A(l, i) = A_il * cos_theta + A_ik * sin_theta;
        };
    }
    // update R
    for (int i = 0; i < N; i++){
        double R_ik = R(i,k);  // in case l = k
        double R_il = R(i,l);
        R(i,k) = R_ik * cos_theta - R_il * sin_theta;
        R(i,l) = R_il * cos_theta - R_ik * sin_theta;
    }
    
}
    

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);

