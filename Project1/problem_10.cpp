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
    for(int j = 0; j < power; j++){

        //maaking vectors (matrix)

        int n = std::pow(10.0, j+1);  
        double h = 1.0/(n+1);
        
        std::vector<double> a(n-1, -1.0);  // subdiagonal a
        std::vector<double> b(n, 2.0);   // diagonal b
        std::vector<double> c(n-1, -1.0);  // superdiagonal c

        std::vector<double> v(n, 0.0);        // approximate solution v
        std::vector<double> g(n, 0.0);   // right-hand side g


        // builds g-vector (RHS)
        for (int i = 0; i < n; i++){
            double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
            g[i] = h * h * 100.0 * std::exp(-10.0 * x);
        }

        // start timing
        auto start = std::chrono::high_resolution_clock::now();

        // running optimalized algo 1000 times
        for(int p = 0; p < 1000; p++){
            std::vector<double> gtemp = g;  // temporary vector
            std::vector<double> btemp = b;   // diagonal b temp vector
            std::vector<double> v_10 = v;   // subdiagonal c temp vector

            // Forward sub

            for (int i = 1; i < n; i++){
                btemp[i]   = 2.0 - 1.0 / btemp[i-1];
                gtemp[i] = g[i] + gtemp[i-1] / btemp[i-1];
            }

            // backward sub
            v_10[n-1] = gtemp[n-1] / btemp[n-1];

            for (int i = n-2; i >= 0; i--) {
                v_10[i] = (gtemp[i] + v_10[i+1]) / btemp[i]; // back-substitute into v
            }
            // do not add boundary values here as we only look at time for algo
    }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        time_opt.push_back(duration.count());
        
        // Problem 8 vectors (standard algo)
        // start timing

        auto start2 = std::chrono::high_resolution_clock::now();
        for(int p = 0; p < 1000; p++){
            // builds g-vector (RHS)
            std::vector<double> v_8 = v;        // approximate solution v
            std::vector<double> gtemp_8 = g;  // temporary vector
            std::vector<double> btemp_8 = b;   // diagonal b temp vector
            
            // Forward sub

            for (int i = 1; i < n; i++) {
                double k_i = a[i] / btemp_8[i-1];
                btemp_8[i] = btemp_8[i] - c[i-1] * k_i;
                gtemp_8[i] = gtemp_8[i] - k_i * gtemp_8[i-1];
            }

            // backward sub
            v_8[n-1] = gtemp_8[n-1] / btemp_8[n-1];

            for (int i = n-2; i >= 0; i--) {
                v_8[i] = (gtemp_8[i] -c[i]* v_8[i+1]) / btemp_8[i]; // back-substitute into v
            }
            // do not add boundary values here as we only look at time for algo
    }
        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration2 = end2 - start2;
        time_org.push_back(duration2.count());
    
    int width = 16;
    int prec = 6;

    double prosent;
    prosent = 100 * (time_org[j] - time_opt[j]) / time_org[j];

        //writing
        ofile 
        << std::setw(width) << std::setprecision(prec) << std::scientific << time_opt[j]
        << std::setw(width) << std::setprecision(prec) << std::scientific << time_org[j]
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