#include "analytical.hpp"

double analytical_eigenvalue(int N, double a, double d, double j){
    // N = matrix size
    // j = index of eigenvalue
    return d + 2*a*cos((j)*M_PI/(N+1));
}

arma::vec analytical_eigenvector(int N, int j){
    // N = matrix size
    // j = index of eigenvector
    arma::vec vec(N);
    for (int i = 0; i < N; i++){
        vec(i) = sin((i+1)*(j)*M_PI/(N+1));
    }
    return vec;
}