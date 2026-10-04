#include <cmath>
#include <iostream>
#include <fstream>
#include "Norms.h"

extern double EPS;
extern double EPS0;

void Simple_iteration(const double* const* A_original, const double* b_original, double* x, const int n, const double tau, int* converged_step, \
    double* norm_C, const std::string output_filename, bool flag = 0) {
    double** A = new double* [n + 1];
    for (int i = 1; i < n + 1; i++) {
        A[i] = new double[n + 1];
    }
    double* b = new double[n + 1];

    for (int i = 1; i < n + 1; ++i) {
        b[i] = b_original[i];
        for (int j = 1; j < n + 1; ++j) {
            A[i][j] = A_original[i][j];
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

    if (norm_C_inf >= 1.0) {
        std::cerr << "||C||_inf >= 1. Method may not converge\n";
    }

    *norm_C = norm_C_inf;

    for (int i = 1; i < n + 1; i++) {
        x[i] = y[i];
    }

    double* x1 = new double[n + 1];

    for (int i = 1; i <= n; i++) {
        double sum = 0.0;
        for (int j = 1; j <= n; j++) {
            sum += A[i][j] * x[j];
        }
        x1[i] = x[i] - tau * (sum - b[i]);
    }

    double* rho_vec = new double[n + 1]();
    for (int i = 1; i <= n; i++) {
        rho_vec[i] = x1[i] - x[i];
    }
    double rho_0 = vector_norm_inf(rho_vec, n);
 
    double k_est = (log(EPS) + log(1 - norm_C_inf) - log(rho_0)) / log(norm_C_inf);
    std::cout << "Convergence rate estimate >= " << k_est << "\n";
    delete[] x1;
    delete[] rho_vec;

    if (!flag) {
        const int max_iterations = std::max(1000, (int)(k_est * 2));
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
            double x_old_norm = vector_norm_inf(x_old, n);          //3
            for (int i = 1; i < n + 1; i++) {
                //diff[i] = x[i] - x_old[i];                        //1, 2
                diff[i] = (x[i] - x_old[i]) / (x_old_norm + EPS0);  //3
            }
            diff_norm = vector_norm_inf(diff, n);
            //criterion = EPS * ((1 - norm_C_inf) / norm_C_inf);    //1
            criterion = EPS;                                        //2, 3

            if (diff_norm <= criterion) {
                std::cout << "Method converged in " << step << " iterations\n";
                *converged_step = step;
                std::ofstream out(output_filename);
                for (int i = 1; i <= n; i++) {
                    out << x[i] << ' ';
                }
                out << "\nk_est >= " << k_est;
                out << "\nsteps = " << step;
                out << "\ntau = " << tau;
                out << "\ndiff_norm = " << diff_norm;
                out << "\nerr_norm = " << residual_norm(A_original, b_original, x, n);
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
                delete[] x_old;
                delete[] diff;
                return;
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
        delete[] x_old;
        delete[] diff;
        return;
    }
    else {
        const int max_iterations = k_est + 1;
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

        std::cout << "Method converged in " << floor(k_est + 1) << " iterations\n";
        std::ofstream out(output_filename);
        for (int i = 1; i <= n; i++) {
            out << x[i] << ' ';
        }
        out << "\nk_est >= " << k_est;
        out << "\nsteps = " << floor(k_est + 1);
        out << "\ntau = " << tau;
        out << "\ndiff_norm = " << diff_norm;
        out << "\nerr_norm = " << residual_norm(A_original, b_original, x, n);
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
        delete[] x_old;
        delete[] diff;
        return;
    }
}