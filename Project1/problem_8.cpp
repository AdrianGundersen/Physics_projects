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


// Insert argumentes and vectors
double problem_8ab(int k, bool check) // bool is to check if it is problem 10
{
    const int n = std::pow(10.0,k);  
    double h = 1.0/(n+1);

    if(!check){
        // Creating a directory and filename for the file to be saved
        std::string folder = "output/";
        std::string filename = "problem_8_";
        std::string txt = ".txt";
        namespace fs = std::filesystem;
        fs::create_directories(folder);
        std::string filepath = folder + filename + std::to_string(n) + txt;
        ofile.open(filepath);
}

    

    //step size

    // vector med x- verdier fra [0,1]
    std::vector<double> x_values;
    for (int j = 0; j < n +1; j++){
        x_values.push_back(j*h);
    };
    x_values.push_back(1.0);

    std::vector<double> a(n, -1.0);  // superdiagonal a
    std::vector<double> b(n, 2.0);   // diagonal b
    std::vector<double> c(n, -1.0);  // subdiagonal c

    std::vector<double> v(n);        // approximate solution v
    std::vector<double> g(n, 0.0);   // right-hand side g
    
    // builds g-vector (RHS)
    for (int i = 0; i < n; i++){
        double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
        g[i] = h * h * 100.0 * std::exp(-10.0 * x);
    }
    
    std::vector<double> gtemp = g;  // temporary vector
    std::vector<double> btemp = b;   // diagonal b temp vector

    // Forward sub

    for (int i = 1; i < n; i++) {
        btemp[i]   = b[i] - c[i-1] * a[i] / btemp[i-1];
        gtemp[i] = g[i] - a[i] * gtemp[i-1] / btemp[i-1];
    }


    // backward sub
    v[n-1] = gtemp[n-1] / btemp[n-1];

    for (int i = n-2; i >= 0; i--) {
        v[i] = (gtemp[i] -c[i]* v[i+1]) / btemp[i]; // back-substitute into v
    }

    

    // adds boundary values
    v.insert(v.begin(), 0.0);
    v.push_back(0.0);  

    // the theoretical vector u
    std::vector<double> u(n+2);

    // the absolute error between u and v
    std::vector<double> delta = u;
    std::vector<double> epsilon_log = delta;
    std::vector<double> epsilon = delta;

    // Calculatign the exact solution and the absolute error
    for(int i = 0; i < n +1; i++){
        u[i] = 1 - (1-std::exp(-10.0))*x_values[i]-std::exp(-10.0*x_values[i]); // exact solution
        //taiking their absolute values
        double diff = std::abs(u[i] - v[i]);
        double rel_err = std::abs((u[i]-v[i])/u[i]);

        // taking log_10
        delta[i] = log10(diff);   
        epsilon[i] = rel_err;
        epsilon_log[i] = log10(rel_err);
    }

    //Trekker ut maksimalveriden i epsilon vectoren
    double max_eps = *std::max_element(epsilon.begin() +1, epsilon.end() - 1);

    if(!check){ 
        // Removing the end ppoints, because we dont need the values where the function callaoses. 
        // sets precision and width of output
        int width = 16;
        int prec = 6;
        for (int i = 1; i < n+1; i++){
            ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x_values[i]
                << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
                << std::setw(width) << std::setprecision(prec) << std::scientific << u[i]
                << std::setw(width) << std::setprecision(prec) << std::scientific << delta[i]
                << std::setw(width) << std::setprecision(prec) << std::scientific << epsilon_log[i]
                << "\n";
        }
        ofile.close();
}

    return max_eps;
}

void problem_10(int k){
    std::vector<int> nums = {10, 100, 1000, 10000, 100000, 1000000};
    // time taken
    std::vector<double> time_opt;
    std::vector<double> time_org;

    
    std::string folder_err = "output/";
    std::string filename_time_opt = "time_optimized_algo.txt";
    std::string filepath_time_opt = folder_err + filename_time_opt;

    ofile.open(filepath_time_opt);

    for(int j = 0; j < nums.size(); j++){
        auto start = std::chrono::high_resolution_clock::now();     
        int n = std::pow(10.0, j+1);  
        double h = 1.0/(n+1);
        std::vector<double> a(n, -1.0);  // superdiagonal a
        std::vector<double> b(n, 2.0);   // diagonal b
        std::vector<double> c(n, -1.0);  // subdiagonal c

        std::vector<double> v(n);        // approximate solution v
        std::vector<double> g(n, 0.0);   // right-hand side g
        
        

        // builds g-vector (RHS)
        for (int i = 0; i < n; i++){
            double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
            g[i] = h * h * 100.0 * std::exp(-10.0 * x);
        }

        for(int p = 0; p < 100; p++){

            std::vector<double> gtemp = g;  // temporary vector
            std::vector<double> btemp = b;   // diagonal b temp vector

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
    }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        time_opt.push_back(duration.count());

        // Problem 8 time 
        auto start2 = std::chrono::high_resolution_clock::now();
        for(int p = 0; p < 100; p++){
            double dummy = problem_8ab(j+1, true);
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
    std::ofstream ofiles;

    // Creating a directory and filename for the file to be saved
    
    std::string folder_err = "output/";
    std::string filename_err = "max_relative_error.txt";
    std::string filepath_err = folder_err + filename_err;

    ofiles.open(filepath_err);

    int width = 16;
    int prec = 6;

    for (int n = 1; n < 6; n++){
        double maxeps = problem_8ab(n, false);
        double N = std::pow(10.0,n);
        double h_tabell = 1.0 / (N + 1);
        ofiles << std::setw(width) << std::setprecision(prec) << std::scientific << maxeps
              << std::setw(width) << std::setprecision(0) << std::scientific << N
              << std::setw(width) << std::setprecision(8) << std::scientific << h_tabell
              << "\n";
    }

    ofiles.close();

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