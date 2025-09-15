#include "tridiag.hpp"
#include <assert.h>

arma::Mat<double> create_tridiagonal(int N, double a, double d){
    arma::Mat<double> A(N, N, arma::fill::zeros);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i == j){
                A(i, j) = d;
            }
            else if (i-j == 1){
                A(i, j) = a;
            }
            else if (j-i == 1){
                A(i, j) = a;
            }
        }
    }
    return A;
}


double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    assert(A.is_square());
    assert(A.n_rows > 1);

    int N = A.n_rows;
    k = 0;
    l = 1;

    double max_val = std::abs(A(k,l)); // initial max value

    for (int i = 0; i < N; i++){
        for (int j = i+ 1 ; j < N; j++){ // only look at upper triangle
            double abs_val = std::abs(A(i, j));
            if (abs_val > max_val){
                max_val = abs_val;
                k = i;
                l = j;
            }
        }
    }
    return max_val;
}

