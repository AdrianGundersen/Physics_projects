#include <armadillo>
using namespace arma;

inline int pos(int i, int j, int m) { return i*m + j; }

void build_AB(int M, double h, double dt, const mat& V, cx_mat& A, cx_mat& B)
{
    int m = M - 2;
    int N = m * m;

    A.zeros(N, N);
    B.zeros(N, N);

    cx_double I(0.0, 1.0);
    cx_double r = I * dt / (2.0 * h * h);

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            int k  = pos(i, j, m);
            int ip = i + 1;      
            int jp = j + 1;

            double vij = V(jp, ip);
            cx_double a = 1.0 + 4.0*r + I*0.5*dt*vij;
            cx_double b = 1.0 - 4.0*r - I*0.5*dt*vij;

            A(k,k) = a;
            B(k,k) = b;

            if (i < m - 1) {
                int kR = pos(i + 1, j, m);
                A(k, kR) = A(kR, k) = -r;
                B(k, kR) = B(kR, k) =  r;
            }
            if (j < m - 1) {
                int kU = pos(i, j + 1, m);
                A(k, kU) = A(kU, k) = -r;
                B(k, kU) = B(kU, k) =  r;
            }
        }
    }
}

int main()
{
    int M = 5;
    double h  = 1.0 / (M - 1);
    h = 2.0;
    double dt = 10;

    mat V(M, M, fill::zeros);
    cx_mat A, B;
    mat A2, B2;

    build_AB(M, h, dt, V, A, B);
    
    // A2 = real(A); B2 = real(B);
    A2 = imag(A); B2 = imag(B);


    // complex vector was hard to comfirm
    A2.print("A:");
    B2.print("B:");

    return 0;
}
