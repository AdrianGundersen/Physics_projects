#include "jacobi.hpp"
#include "tridiag.hpp"



int main(){
    int N = 6;
    double a = 2.0;
    double d = -1.0;
    int k = 0;
    int l = 0;
    arma::mat R(N, N);
    R.eye();
    arma::Mat<double> A = create_tridiagonal(N, a, d);
    arma::vec eigenval = arma::eig_sym(A);
    for(int i = 0; i < 1000; i++){
        max_offdiag_symmetric(A, k, l);
        jacobi_rotate(A, R, k, l);
    }
    A.print();
    eigenval.print();

}