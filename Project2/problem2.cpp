#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>


int main(){
    int N = 6;
    double a = -1.;
    double b = 1.;
    double c = -2.;
    arma::Mat<double> A(N, N, arma::fill::zeros);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            if (i ==j){
                A(i, j) = b;
            }
            else if (i-j == 1){
                A(i, j) = a;
            }
            else if (j-i == 1){
                A(i, j) = c;
            }
        }
    }
    A.print();
}