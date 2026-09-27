#include <cmath>
#include <iostream>
#include <fstream>
#include "Norms.h"

extern double EPS;

void Seidel(const double* a, const double* b, const double* c, const double* d, double* x, const int n, int* converged_step, \
    double* norm_C, const std::string output_filename) {

    double** C = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        C[i] = new double[n + 1]();
    }
    double* y = new double[n + 1]();

    for (int i = 1; i < n + 1; i++) {
        y[i] = d[i] / b[i];
        if (i < n) {
            C[i][i + 1] = -c[i] / b[i];
        }
        if (i > 1) {
            y[i] += (-a[i] / b[i]) * y[i - 1];
            for (int j = 1; j < n + 1; j++) {
                C[i][j] += (-a[i] / b[i]) * C[i - 1][j];
            }
        }
    }

    //Testing ||C|| < 1
    double norm_C_l1 = matrix_norm_l1(C, n);
    double norm_C_inf = matrix_norm_inf(C, n);

    std::cout << "\nSeidel method:\n";
    std::cout << "||C||_l1 = " << norm_C_l1 << "\n";
    std::cout << "||C||_inf = " << norm_C_inf << "\n";

    *norm_C = norm_C_inf;

    double** C_L = new double* [n + 1];
    double** C_U = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        C_L[i] = new double[n + 1]();
        C_U[i] = new double[n + 1]();
        for (int j = 1; j <= n; j++) {
            if (i > j) {
                C_L[i][j] = C[i][j]; //lower triangular
            }
            else if (i < j) {
                C_U[i][j] = C[i][j]; //upper triangular
            }
        }
    }

    double norm_C_L_inf = matrix_norm_inf(C_L, n);
    double norm_C_U_inf = matrix_norm_inf(C_U, n);
    double sum_norms = norm_C_L_inf + norm_C_U_inf;

    std::cout << "||C_L||_inf = " << norm_C_L_inf << "\n";
    std::cout << "||C_U||_inf = " << norm_C_U_inf << "\n";
    std::cout << "||C_L||_inf + ||C_U||_inf = " << sum_norms << "\n";

    if (sum_norms < 1.0) {
        double q = norm_C_U_inf/(1 - norm_C_L_inf);
        std::cout << "Convergence rate estimate (based on sum of norms) >= " << std::log(EPS)/std::log(q) << "\n";
    }
    else {
        std::cout << "Sum of norms >= 1. Convergence is not guaranteed by this criterion.\n";
    }

    for (int i = 1; i < n + 1; i++) {
        x[i] = y[i];
    }

    int max_iterations = 1000;
    double err;
    double* x_old = new double[n + 2]();

    for (int step = 0; step < max_iterations; step++) {
        err = residual_norm_three_diag(a, b, c, d, x, n);
        if (err < EPS * ((1 - norm_C_inf) / norm_C_U_inf)) {
            std::cout << "Method converged in " << step << " iterations\n";
            *converged_step = step;
            std::ofstream out(output_filename);
            for (int i = 1; i <= n; i++) {
                out << x[i] << ' ';
            }
            out << "\nsteps = " << step;
            out << "\nerr = " << err;
            out << "\nNorm_l1 C = " << norm_C_l1;
            out << "\nNorm_inf C = " << norm_C_inf;
            out.close();
            for (int i = 1; i <= n; i++) {
                delete[] C[i];
            }
            delete[] C;
            delete[] y;
            delete[] x_old;
            return;
        }
        for (int i = 1; i <= n; i++) {
            x_old[i] = x[i];
        }
        for (int i = 1; i < n + 1; i++) {
            double sum = d[i];
            if (i > 1) {
                sum -= a[i] * x[i - 1];
            }
            if (i < n) {
                sum -= c[i] * x_old[i + 1];
            }
            x[i] = sum / b[i];
        }
    }

    std::cout << "The maximum number of iterations has been reached\n";
    for (int i = 1; i <= n; i++) {
        delete[] C[i];
    }
    delete[] C;
    delete[] y;
    delete[] x_old;
    return;
}