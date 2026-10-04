#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "Norms.h"
#include "SimpleIteration.h"
#include "Jacobi.h"
#include "Seidel.h"
#include "Relaxation.h"

double EPS = 1e-4;
//double EPS = 1e-7;
double EPS0 = 1e-14;

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
    double tau_theor = 1 / matrix_norm_inf(A, n);
    double tau_limit = 2 * tau_theor;

    double step_tau = tau_limit / 50.0;
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
    std::cout << "\nInput N: ";
    int N;
    std::cin >> N;
    int n_big = 200 + N;

    double* a = new double[n_big + 1];
    double* b_big = new double[n_big + 1];
    double* c = new double[n_big + 1];
    double* d = new double[n_big + 1];

    for (int i = 1; i < n_big + 1; ++i) {
        b_big[i] = 4.0;
        if (i > 1) {
            a[i] = 1.0;
        }
        if (i < n_big) {
            c[i] = 1.0;
        }
    }

    for (int i = 1; i < n_big + 1; ++i) {
        if (i == 1) {
            d[i] = 6.0;
        }
        else if (i == n_big) {
            d[i] = 9.0 - 3.0 * (n_big % 2);
        }
        else {
            d[i] = 10.0 - 2.0 * (i % 2);
        }
    }

    //Seidel method
    double* x_big = new double[n_big + 1]();
    output_filename = std::string("Seidel") + file_number + '_' + accuracy + ".txt";
    Seidel(a, b_big, c, d, x_big, n_big, &converged_step, &norm_C, output_filename);

    //Relaxation method
    output_filename = std::string("Relaxation") + file_number + '_' + accuracy + ".txt";

    double start_omega = 1.0;
    double end_omega = 2;
    double step_omega = 0.05;

    double best_omega = -1;
    double min_steps_rel = 999999;
    int converged_step_rel = -1;
    double norm_C_rel;

    for (double omega = start_omega; omega < end_omega; omega += step_omega) {
        for (int i = 1; i <= n_big; i++) x_big[i] = 0.0;
        Relaxation_triag(a, b_big, c, d, x_big, n_big, omega, &converged_step_rel, &norm_C_rel, output_filename);
        if (converged_step_rel < min_steps_rel && converged_step_rel != -1) {
            min_steps_rel = converged_step_rel;
            best_omega = omega;
        }
    }

    std::cout << "\nBest experimental omega = " << best_omega << " (converged in " << min_steps_rel << " steps)\n";

    for (int i = 1; i <= n_big; i++) x_big[i] = 0.0;
    Relaxation_triag(a, b_big, c, d, x_big, n_big, best_omega, &converged_step_rel, &norm_C_rel, output_filename);

    //Destructors
    for (int i = 1; i <= n; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;
    delete[] a;
    delete[] b_big;
    delete[] c;
    delete[] d;
    return 0;
}