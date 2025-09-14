#include "tridiag.hpp"

arma::Mat<double> create_tridiagonal(int N, double a, double d){
    arma::Mat<double> A(N, N, arma::fill::zeros);
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            if (i == j)            A(i, j) = d;
            else if (i - j == 1)   A(i, j) = a;
            else if (j - i == 1)   A(i, j) = a;
        }
    }
    return A;
}
