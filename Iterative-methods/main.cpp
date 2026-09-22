#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "Norms.h"
#include "SimpleIteration.h"
#include "Jacobi.h"

double EPS = 1e-4;
//double EPS = 1e-7;

int main() {
    //Reading from file
    std::string filename;
    std::cin >> filename;

    int n;
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return 0;
    }
    file >> n;

    double** A = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        A[i] = new double[n + 1];
    }
    double* b = new double[n + 1];
    double* x = new double[n + 1];

    double val;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            file >> val;
            if (j <= n) {
                A[i][j] = val;
            }
            if (j == n + 1) {
                b[i] = val;
            }
        }
    }
    file.close();

    //Checking accuracy for naming output file
    std::string accuracy;
    if (EPS == 1e-4) {
        accuracy = "1e-4";
    }
    if (EPS == 1e-7) {
        accuracy = "1e-7";
    }
    std::string file_number = filename.substr(4, filename.length() - 8);
    std::string output_filename = std::string("Simple_iteration") + file_number + '_' + accuracy + ".txt";

    //Experimental determination of tau for simple iteration method
    double tau_theor = 1 / matrix_norm_l1(A, n);
    double tau_limit = 2 * tau_theor;

    double step_tau = tau_limit / 20.0;
    double start_tau = step_tau;
    double end_tau = tau_limit * 1.2;

    double best_tau = -1;
    int min_steps = 999999;
    int converged_step = -1;
    double norm_C;

    for (double tau = start_tau; tau <= end_tau; tau += step_tau) {
        Simple_iteration(A, b, x, n, tau, &converged_step, &norm_C, output_filename);

        std::cout << "tau = " << tau << "\n\n";
        if (converged_step < min_steps and converged_step != -1 and norm_C < 1) {
            min_steps = converged_step;
            best_tau = tau;
        }
    }
    Simple_iteration(A, b, x, n, tau_theor, &converged_step, &norm_C, output_filename);
    std::cout << "\nTheoretical tau = " << tau_theor << " (converged in " << converged_step << ", norm_inf C = " << norm_C << ")\n\n";
    if (best_tau > 0) {
        Simple_iteration(A, b, x, n, best_tau, &converged_step, &norm_C, output_filename);
        std::cout << "\nExperimental tau: " << best_tau << " (converged in " << min_steps << ", norm_inf C = " << norm_C << ")\n\n";
    }

    //Jacobi method
    output_filename = std::string("Jacobi") + file_number + '_' + accuracy + ".txt";
    Jacobi(A, b, x, n, &converged_step, &norm_C, output_filename);

    //Tridiagonal matrix
    int N;
    std::cin >> N;
    int n_big = 200 + N;

    double* a = new double[n];
    double* b = new double[n];
    double* c = new double[n];
    double* d = new double[n];

    for (int i = 0; i < n_big; ++i) {
        b[i] = 4.0;
        a[i] = 1.0;
        c[i] = 1.0;
    }

    for (int idx = 0; idx < n_big; ++idx) {
        int i = idx + 1;

        if (i == 1) {
            d[idx] = 6.0;
        }
        else if (i == n_big) {
            d[idx] = 9.0 - 3.0 * (n_big % 2);
        }
        else {
            d[idx] = 10.0 - 2.0 * (i % 2);
        }
    }



    //Destructors
    for (int i = 1; i <= n; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    return 0;
}