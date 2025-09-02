#include <iostream>
#include <vector>
#include <cmath>

int main(){
    int n = 10000;
    std::vector<double> a(n, -1);
    std::vector<double> b(n, 2);
    std::vector<double> c = a;
    std::vector<double> b_tilde(n, 0);
    std::vector<double> g = b_tilde;
    std::vector<double> g_tilde = b_tilde;
    std::vector<double> v = b_tilde;

    // filling sourceing term
    double dx = (1.0/n);
    for(int i = 0; i < n; i++){
        double x = i*dx;
        g[i] = 100*std::exp(-10*x)*dx*dx;
}
    // initial states
    b_tilde[0] = b[0]; 
    g_tilde[0] = g[0];
    for(int i = 1; i < n; i++){
        // Forward sub
        b_tilde[i] = b[i] - a[i]*c[i-1]/b_tilde[i-1];
        g_tilde[i] = g[i] - a[i]*g_tilde[i-1]/b_tilde[i-1];
}   
    // last position of v
    v[n-1] = g_tilde[n-1]/b_tilde[n-1];
    // solving for v 
    for(int i = (n-2); i >= 0; i--){
        v[i] = (g_tilde[i] - c[i]*v[i+1])/b_tilde[i];
    }
    // initial conditions
    v.insert(v.begin(), 0.0);
    v.push_back(0.0);

    // output
    for(int i = 0; i < n+2; i++){
        std::cout << v[i] << " " << i*dx << "\n";
    }

    return 0;
}