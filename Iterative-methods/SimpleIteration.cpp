#include <cmath>
#include <iostream>
#include <fstream>
#include "Norms.h"

extern double EPS;

void Simple_iteration(const double* const* A_original, const double* b_original, double* x, const int n, const double tau, int* converged_step, \
    double* norm_C, const std::string output_filename) {
    double** A = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        A[i] = new double[n + 1];
    }
    double* b = new double[n + 1];

    for (int i = 1; i < n + 1; ++i) {
        for (int j = 1; j < n + 1; ++j) {
            A[i][j] = A_original[i][j];
            if (j == n) {
                b[i] = b_original[i];
            }
        }
    }

    //Making all diagonal elements positive
    for (int i = 1; i < n + 1; ++i) {
        if (A[i][i] < 0.0) {
            b[i] = -b[i];
            for (int j = 1; j < n + 1; ++j) {
                A[i][j] = -A[i][j];
            }
        }
    }
    
    double** C = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        C[i] = new double[n + 1];
    }
    double* y = new double[n + 1];

    for (int i = 1; i < n + 1; i++) {
        y[i] = tau * b[i];
        for (int j = 1; j < n + 1; j++) {
            if (i == j) {
                C[i][j] = 1.0 - tau * A[i][j]; // E - tau*A on diag
            }
            else {
                C[i][j] = -tau * A[i][j];      // -tau*A out of diag
            }
        }
    }

    //Testing ||C|| < 1
    double norm_C_l1 = matrix_norm_l1(C, n);
    double norm_C_inf = matrix_norm_inf(C, n);

    std::cout << "Simple iteration method:\n";
    std::cout << "||C||_l1 = " << norm_C_l1 << "\n";
    std::cout << "||C||_inf = " << norm_C_inf << "\n";

    *norm_C = norm_C_inf;

    double* x0 = new double[n + 1];
    for (int i = 1; i < n + 1; i++) {
        x[i] = y[i];
    }

    int max_iterations = 1000;
    double err;
    for (int step = 0; step < max_iterations; step++) {
        err = residual_norm(A_original, b_original, x, n);
        if (err < EPS) {
            std::cout << "Method converged in " << step << " iterations\n";
            *converged_step = step;
            std::ofstream out(output_filename);
            for (int i = 1; i <= n; i++) {
                out << x[i] << ' ';
            }
            out << "\nsteps = " << step;
            out << "\ntau = " << tau;
            out << "\nerr = " << err;
            out << "\nNorm_l1 C = " << norm_C_l1;
            out << "\nNorm_inf C = " << norm_C_inf;
            out.close();
            for (int i = 1; i <= n; i++) {
                delete[] A[i];
                delete[] C[i];
            }
            delete[] A;
            delete[] b;
            delete[] y;
            return;
        }

        // Вычисление x_next = C * x + y
        for (int i = 1; i < n + 1; i++) {
            double sum = 0.0;
            for (int j = 1; j < n + 1; j++) {
                sum += C[i][j] * x[j];
            }
            x[i] = sum + y[i];
        }
    }

    std::cout << "The maximum number of iterations has been reached\n";

    for (int i = 1; i <= n; i++) {
        delete[] A[i];
        delete[] C[i];
    }
    delete[] A;
    delete[] b;
    delete[] y;
    return;
}