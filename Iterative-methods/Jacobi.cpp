#include <cmath>
#include <iostream>
#include <fstream>
#include "Norms.h"

extern double EPS;
extern double EPS0;

void Jacobi(const double* const* A, const double* b, double* x, const int n, int* converged_step, \
    double* norm_C, const std::string output_filename) {

    double** C = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        C[i] = new double[n + 1];
    }
    double* y = new double[n + 1];

    for (int i = 1; i < n + 1; i++) {
        y[i] = b[i] / A[i][i];
        for (int j = 1; j < n + 1; j++) {
            if (i == j) {
                C[i][j] = 0.0;
            }
            else {
                C[i][j] = -A[i][j] / A[i][i];
            }
        }
    }

    //Testing ||C|| < 1
    double norm_C_l1 = matrix_norm_l1(C, n);
    double norm_C_inf = matrix_norm_inf(C, n);

    std::cout << "Jacobi method:\n";
    std::cout << "||C||_l1 = " << norm_C_l1 << "\n";
    std::cout << "||C||_inf = " << norm_C_inf << "\n";

    *norm_C = norm_C_inf;

    for (int i = 1; i < n + 1; i++) {
        x[i] = y[i];
    }

    double* x0 = new double[n + 1];
    for (int i = 1; i < n + 1; i++) {
        double sum = 0.0;
        for (int j = 1; j < n + 1; j++) {
            sum += A[i][j] * x[j];
        }
        x0[i] = sum - b[i];
    }
    double x0_norm = vector_norm_inf(x0, n);

    double k_est = (log(EPS) - log(x0_norm)) / log(norm_C_inf);
    std::cout << "Convergence rate estimate >= " << k_est << "\n";
    delete[] x0;

    int max_iterations = 1000;
    double* x_old = new double[n + 1];
    double* diff = new double[n + 1];
    double diff_norm;
    double criterion;
    for (int step = 0; step < max_iterations; step++) {
        for (int i = 1; i < n + 1; i++) {
            x_old[i] = x[i];
        }

        //x_next = C * x + y
        for (int i = 1; i < n + 1; i++) {
            double sum = 0.0;
            for (int j = 1; j < n + 1; j++) {
                sum += C[i][j] * x[j];
            }
            x[i] = sum + y[i];
        }

        //Different stopping criteria 
        //double x_old_norm = vector_norm_inf(x_old, n);
        for (int i = 1; i < n + 1; i++) {
            diff[i] = x[i] - x_old[i];
            //diff[i] = (x[i] - x_old[i]) / (x_old_norm + EPS0);
        }
        diff_norm = vector_norm_inf(diff, n);

        criterion = EPS * ((1 - norm_C_inf) / norm_C_inf);
        //criterion = EPS;
        //criterion = EPS * x_old_norm + EPS0;

        if (diff_norm <= criterion) {
            std::cout << "Method converged in " << step << " iterations\n";
            *converged_step = step;
            std::ofstream out(output_filename);
            for (int i = 1; i <= n; i++) {
                out << x[i] << ' ';
            }
            out << "\nk_est >= " << k_est;
            out << "\nsteps = " << step;
            out << "\ndiff_norm = " << diff_norm;
            out << "\nerr_norm = " << residual_norm(A, b, x, n);
            out << "\nNorm_l1 C = " << norm_C_l1;
            out << "\nNorm_inf C = " << norm_C_inf;
            out.close();
            for (int i = 1; i <= n; i++) {
                delete[] C[i];
            }
            delete[] y;
            delete[] x_old;
            return;
        }
    }

    std::cout << "The maximum number of iterations has been reached\n";

    for (int i = 1; i <= n; i++) {
        delete[] C[i];
    }
    delete[] y;
    delete[] x_old;
    return;
}