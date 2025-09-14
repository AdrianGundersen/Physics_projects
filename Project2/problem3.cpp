#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>
#include "tridiag.hpp"

int main(){
    arma::Mat <double> A = arma::Mat<double>(4,4,arma::fill::zeros);
    for (int i = 0; i < 4; ++i){
        A(i,i) = 1.;
        }
    A(0,3) = A(3,0) = 0.5;
    A(1,2) = A(2,1) = -0.7;

    int k;
    int l;

    double max_val;
    max_val = max_offdiag_symmetric(A, k, l);

    std::cout << "Max off-diagonal value: " << max_val << "\n";
    std::cout << "(k,l) = " << k << "," << l << "\n";
}