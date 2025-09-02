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

std::ofstream ofile;

void print_usage(const char* progname) {
    // made using chatgpt
    std::cout << "Usage: " << progname << " outputfile\n";
    std::cout << "\nDescription:\n"
              << "  Solves -u''(x) = 100 e^{-10x} on (0,1) with u(0)=u(1)=0\n"
              << "  using the Thomas algorithm for tridiagonal systems.\n"
              << "\nArguments:\n"
              << "  outputfile   Name of the file to write (x, u(x)) pairs\n"
              << "\nInput:\n"
              << "  The program will then ask you for the number of mesh points n.\n"
              << std::endl;
}

// Insert argumentes and vectors
int main(int argc, char* argv[])
{
    //  usage message
    if (argc > 1 &&
       (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
        print_usage(argv[0]);
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();
    char *outfilename;
    int i, j, n;
    double h;
    // Gir en feilmelding om du ikke har med rette argumenter
    if( argc <= 1)
    {
        std::cout << "Bad Usage: " << argv[0] << 
        " write the desired filename of output file" << std::endl;
        exit(1);
    }
    else{
        outfilename = argv[1];
    }

    // Creating a directory for the file to be saved
    std::string folder = "output/";
    namespace fs = std::filesystem;
    fs::create_directories(folder);
    std::string filepath = folder + outfilename;

    ofile.open(filepath);

    std::cout << "Read in number of mech points" << std::endl;

    //Reads # of mesh points
    std::cin >> n;

    //step size
    h = 1.0/( (double) n+1);

    // vector med x- verdier fra [0,1]
    std::vector<double> x_values;
    for (int j = 0; j < n; j++){
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

    // sets precision and width of output
    int width = 16;
    int prec = 6;

    // adds boundary values
    v.insert(v.begin(), 0.0);
    v.push_back(0.0);  

    // the theoretical vector u
    std::vector<double> u(n+2);

    // the absolute error between u and v
    std::vector<double> delta = u;

    // Calculatign the exact solution and the absolute error
    for(int i = 0; i < n; i++){
        u[i] = 1 - (1-std::exp(-10.0))*x_values[i]-std::exp(-10.0*x_values[i]); // exact solution
        double diff = std::abs(u[i] - v[i]);
        delta[i] = log10(diff);   // taking difference
    }
    u.push_back(0.0); // boundary value

    for (int i = 0; i <= n; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x_values[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << u[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << delta[i]
              << "\n";
    }
    ofile.close();

    // Writes duration of the program
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds\n";

    return 0;
}
