#include <iostream>
#include <vector>
#include <armadillo>
#include <cmath>


arma::Mat<double> create_tridiagonal(int N, double a, double b, double c){
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
    return A;
}

int main(){
    int N = 6;
    double h = 1./(N+1);
    double a = -1/(h*h);
    double b = 2/(h*h);
    double c = -1/(h*h);   
    arma::Mat<double> A = create_tridiagonal(N, a, b, c);
    arma::vec eigvals = arma::eig_sym(A);
    A.print("Matrix A: ");
    eigvals.print("Eigenvalues: ");
}