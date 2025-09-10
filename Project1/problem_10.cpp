#include <iostream>
#include <vector>
#include <numbers> 
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <filesystem>
#include <algorithm>

std::ofstream ofile;

// Function for the optimized algorythm
void optimal(int n,
    std::vector<double>& v,
    const std::vector<double>& b,
    const std::vector<double>& g,
    std::vector<double>& btemp,
    std::vector<double>& gtemp){
        // Forward sub

        for (int i = 1; i < n; i++){
            btemp[i]   = 2.0 - 1.0 / btemp[i-1];
            gtemp[i] = gtemp[i] + gtemp[i-1] / btemp[i-1];
        }

        // backward sub
        v[n-1] = gtemp[n-1] / btemp[n-1];

        for (int i = n-2; i >= 0; i--) {
            v[i] = (gtemp[i] + v[i+1]) / btemp[i]; // back-substitute into v
        }
        // do not add boundary conditions as this increases time unnescessarily
}

// function for the original algorythm
void original(int n,
    const std::vector<double> & a,
    const std::vector<double> & b, 
    const std::vector<double> & c,
    std::vector<double> & v,
    const std::vector<double> & g,
    std::vector<double>& btemp,
    std::vector<double>& gtemp) 
    {    
        // Forward sub
        for (int i = 1; i < n; i++) {
            double k_i = a[i-1] / btemp[i-1];
            btemp[i]   = btemp[i] - c[i-1] * k_i;
            gtemp[i] = gtemp[i] - k_i * gtemp[i-1];
        }

        // backward sub
        v[n-1] = gtemp[n-1] / btemp[n-1];

        for (int i = n-2; i >= 0; i--) {
            v[i] = (gtemp[i] -c[i]* v[i+1]) / btemp[i]; // back-substitute into v
        }
        // do not add boundary conditions as this increases time unnescessarily
     }
void problem_10(){
    int power = 6; // number to power of 10

    // time taken for the diffrent algo
    std::vector<double> time_opt; 
    std::vector<double> time_org;

    //file 
    std::string folder_err = "output/";
    std::string filename_time_opt = "time_optimized_algo.txt";
    std::string filepath_time_opt = folder_err + filename_time_opt;

    ofile.open(filepath_time_opt);

    //start running 
    for(int j = 1; j < power+1; j++){

        //maaking vectors (matrix)

        int n = std::pow(10.0, j); // number of grid points 
        double h = 1.0/(n+1);

        
        std::vector<double> a(n-1, -1.0);  // superdiagonal a
        std::vector<double> b(n, 2.0);   // diagonal b
        std::vector<double> c(n-1, -1.0);  // subdiagonal c

        std::vector<double> v(n);        // approximate solution v
        std::vector<double> g(n, 0.0);   // right-hand side g

        std::vector<double> btemp = b;
        std::vector<double> gtemp = g;

        // builds g-vector (RHS)
        for (int i = 0; i < n; i++){
            double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
            g[i] = h * h * 100.0 * std::exp(-10.0 * x);
        }

        //Time optimized
        auto start = std::chrono::high_resolution_clock::now();

        for(int p = 0; p < 1000; p++){
            optimal(n, v, b, g, btemp, gtemp);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        time_opt.push_back(duration.count());

        //time original
        auto start2 = std::chrono::high_resolution_clock::now();
        for(int p = 0; p < 1000; p++){
            original(n, a, b, c, v, g, btemp, gtemp);
        }
        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration2 = end2 - start2;
        time_org.push_back(duration2.count());
    
    int width = 16;
    int prec = 6;

    // calculate the prosentage gain 
    double prosent;
    prosent = 100 * (time_org[j-1] - time_opt[j-1]) / time_org[j-1];

        //writing
        ofile 
        << std::setw(width) << std::setprecision(prec) << std::scientific << time_opt[j-1]
        << std::setw(width) << std::setprecision(prec) << std::scientific << time_org[j-1]
        << std::setw(width) << std::setprecision(2) << std::fixed << prosent
        << std::setw(width) << std::setprecision(0) << std::fixed << n
            << "\n";
}
    ofile.close();
}

int main(){
  
    problem_10();
    
    return 0;
}