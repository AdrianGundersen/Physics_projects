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

void problem_10(int k){
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

        //Running opt algo
        auto start = std::chrono::high_resolution_clock::now();

        int n = std::pow(10.0, j+1);  
        double h = 1.0/(n+1);
        
        std::vector<double> a(n, -1.0);  // superdiagonal a
        std::vector<double> b(n, 2.0);   // diagonal b
        std::vector<double> c(n, -1.0);  // subdiagonal c

        std::vector<double> v(n);        // approximate solution v
        std::vector<double> g(n, 0.0);   // right-hand side g

        std::vector<double> gtemp = g;  // temporary vector
        std::vector<double> btemp = b;   // diagonal b temp vector

        for(int p = 0; p < 1000; p++){

             // builds g-vector (RHS)
            for (int i = 0; i < n; i++){
                double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
                g[i] = h * h * 100.0 * std::exp(-10.0 * x);
        }

            // Forward sub

            for (int i = 1; i < n; i++){
                btemp[i]   = 2.0 - 1.0 / btemp[i-1];
                gtemp[i] = g[i] + gtemp[i-1] / btemp[i-1];
            }

            // backward sub
            v[n-1] = gtemp[n-1] / btemp[n-1];

            for (int i = n-2; i >= 0; i--) {
                v[i] = (gtemp[i] + v[i+1]) / btemp[i]; // back-substitute into v
            }
            v.insert(v.begin(), 0.0);
            v.push_back(0.0);  
    }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        time_opt.push_back(duration.count());
        
        // Problem 8 time 
        auto start2 = std::chrono::high_resolution_clock::now();

        std::vector<double> v_8(n);        // approximate solution v
        std::vector<double> g_8(n, 0.0);   // right-hand side g

        std::vector<double> gtemp_8 = g;  // temporary vector
        std::vector<double> btemp_8 = b;   // diagonal b temp vector

        for(int p = 0; p < 1000; p++){
            // builds g-vector (RHS)

            for (int i = 0; i < n; i++){
                double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
                g[i] = h * h * 100.0 * std::exp(-10.0 * x);
                }
            
            // Forward sub

            for (int i = 1; i < n; i++) {
                btemp_8[i]   = b[i] - c[i-1] * a[i] / btemp_8[i-1];
                gtemp_8[i] = g[i] - a[i] * gtemp_8[i-1] / btemp_8[i-1];
            }

            // backward sub
            v_8[n-1] = gtemp_8[n-1] / btemp_8[n-1];

            for (int i = n-2; i >= 0; i--) {
                v_8[i] = (gtemp_8[i] -c[i]* v_8[i+1]) / btemp_8[i]; // back-substitute into v
            }
            // adds boundary values
            v_8.insert(v_8.begin(), 0.0);
            v_8.push_back(0.0);
    }
        auto end2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration2 = end2 - start2;
        time_org.push_back(duration2.count());
    
    int width = 16;
    int prec = 6;

        //writing
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << time_opt[j]
            << std::setw(width) << std::setprecision(prec) << std::scientific << time_org[j]
            << "\n";
}
    ofile.close();
}

int main(){
    auto start = std::chrono::high_resolution_clock::now();

    // Writes duration of the program
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds\n";

    std::cout << "Run problem 10? Y/N \n";

    std::string input;
    std::cin >> input;

    if (input == "Y" || input == "y") {
        problem_10(1e6);
    }
    
    return 0;
}