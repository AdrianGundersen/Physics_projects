#include <iostream>
#include <vector>
#include <numbers> 
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <filesystem>

std::vector<double> u_x(std::vector<double> x_vals);

int main(){
    std::string filename = "diff_eq_sol.txt";

    // Starts timer
    auto start = std::chrono::high_resolution_clock::now(); 
    std::vector<double> x_values; // initialize x-array
    int n = 100000; // # of points
    double dx = 1.0/n; 

    // creates x-array
    for (int j = 0; j < n; j++){
        x_values.push_back(j*dx);
    };
    x_values.push_back(1.0); // adds boundary

    std::vector<double> results = u_x(x_values); // calculates y-values
    
    int width = 12; // aligns columns
    int prec = 4; // # of digits (precision)

    // Creating a folder to save the data inside
    std::ofstream ofile;
    std::string folder = "output/";
    namespace fs = std::filesystem;
    fs::create_directories(folder);

    // starts writing to file
    std::string filepath = folder + filename;
    ofile.open(filepath);

    for (int i = 0; i <= n; i++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x_values[i]
          << std::setw(width) << std::setprecision(prec) << std::scientific << results[i]
          << "\n";
        }
    ofile.close();
    // ends writing to file

    // Ends timer
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate duration
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds\n";
    
}

std::vector<double> u_x(std::vector<double> x_vals){
    std::vector<double> u_values;
    int size = x_vals.size();
    const double e = 2.71828182845904523536;

    for(int i = 0; i < size + 1; i++){
        double u = 1 - (1-std::exp(-10.0))*x_vals[i]-std::exp(-10.0*x_vals[i]); // exact solution
        u_values.push_back(u);
    }
    u_values.push_back(0.0); // boundary value

    return u_values;
};