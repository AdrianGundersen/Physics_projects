// src/matrix_printer.cpp

/*
Simply builds a matrix and prints it because we where told to do it
*/

#include "constants.hpp"
#include <vector>
#include <complex>
#include <iostream>
#include <iomanip>
#include <string>
namespace ds {
    inline Index mat_idx(Index row, Index col, Index N) {
        return row * N + col;
    }

    inline Index grid_index(Index i, Index j, Index M) {
        return i * M + j;
    }

    void build_AB(Index M, Real h, Real dt, const rvec& V, cvec& A, cvec& B)
    {
        const Index m = M - 2;
        const Index N = m * m;

        A.assign(N * N, Complex(0.0, 0.0));
        B.assign(N * N, Complex(0.0, 0.0));

        const Complex I(0.0, 1.0);
        const Complex r = I * dt / (2.0 * h * h);

        auto pos = [m](Index i, Index j) {
            return i * m + j;
        };

        for (Index j = 0; j < m; ++j) {
            for (Index i = 0; i < m; ++i) {
                const Index k  = pos(i, j);
                const Index ip = i + 1;      
                const Index jp = j + 1;

                const Real vij = V[grid_index(jp, ip, M)];
                const Complex a = 1.0 + 4.0*r + I*0.5*dt*vij;
                const Complex b = 1.0 - 4.0*r - I*0.5*dt*vij;

                A[mat_idx(k,k, N)] = a;
                B[mat_idx(k,k, N)] = b;

                if (i < m - 1) {
                    const Index kR = pos(i + 1, j);
                    const Index idx1 = mat_idx(k, kR, N);
                    const Index idx2 = mat_idx(kR, k, N);

                    A[idx1] = A[idx2] = -r;
                    B[idx1] = B[idx2] =  r;
                }
                if (j < m - 1) {
                    const Index kU = pos(i, j + 1);
                    const Index idx1 = mat_idx(k, kU, N);
                    const Index idx2 = mat_idx(kU, k, N);

                    A[idx1] = A[idx2] = -r;
                    B[idx1] = B[idx2] =  r;
                }
            }
        }
    }
    void print_matrix(const cvec& mat, Index N, const std::string& name) {

        std::cout << "Matrix " << name << " imaginary part:" << std::endl;
        std::cout << std::fixed << std::setprecision(2);
        const Index width = 8;
        for (Index row = 0; row < N; ++row) {
            for (Index col = 0; col < N; ++col) {
                const Complex val = mat[mat_idx(row, col, N)];
                std::cout << std::setw(width) << std::right << std::imag(val);
            }
            std::cout << "\n";
        }
    }
} // namespace ds

int main()
{
    const ds::Index M = 5;
    ds::Real h  = 1.0 / (M - 1);
    h = 2.0;
    const ds::Real dt = 10;

    ds::rvec V;
    V.assign(M * M, 0.0);
    
    ds::cvec A, B;
    ds::build_AB(M, h, dt, V, A, B);
    const ds::Index m = M - 2;
    const ds::Index N = m * m;
    
    ds::print_matrix(A, N, "A");
    ds::print_matrix(B, N, "B");

    return 0;
}
