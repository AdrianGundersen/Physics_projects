#include <iostream>
#include <vector>
#include <numbers> 
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>

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
    ofile.open(outfilename);
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
    std::vector<double> temp(n, 0.0);  // temporary vector

    // builds g-vector (RHS)
    for (int i = 0; i < n; i++){
        double x = (i+1)*h;     // x-values from 0 to 1 (w/o boundaries)
        g[i] = h*h * 100.0 * std::exp(-10.0 * x);
    }

    // Forward sub
    double btemp = b[0];
    v[0] = g[0] / btemp; // initialize solution at first interior

    for (int i = 1; i < n; i++) {
        temp[i] = c[i-1] / btemp;
        btemp   = b[i] - a[i] * temp[i];
        v[i]    = (g[i] - a[i]*v[i-1]) / btemp;   // compute solution into v
    }

    // backward sub
    for (int i = n-2; i >= 0; i--) {
        v[i] -= temp[i+1] * v[i+1]; // back-substitute into v
    }

    // sets precision and width of output
    int width = 12;
    int prec = 4;

    v.push_back(0.0);  // adds boundary value

    for (int i = 0; i <= n; i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << x_values[i]
              << std::setw(width) << std::setprecision(prec) << std::scientific << v[i]
              << std::endl;
    }
    ofile.close();

    // Writes duration of the program
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "Time taken: " << duration.count() << " seconds\n";

    return 0;
}
