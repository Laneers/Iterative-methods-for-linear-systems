#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Norms.h"

extern double EPS;
extern double EPS0;

void Relaxation_triag(const double* a, const double* b, const double* c, const double* d, double* x, const int n, const double omega, int* converged_step, \
    double* norm_C, const std::string output_filename, bool flag = 0) {

    double** C = new double* [n + 1];
    double** C_L = new double* [n + 1];     //lower triangular
    double** C_U = new double* [n + 1];     //upper triangular
    double** G2 = new double* [n + 1];      //C_U + C_D
    for (int i = 1; i < n + 1; i++) {
        C[i] = new double[n + 1]();
        C_L[i] = new double[n + 1]();
        C_U[i] = new double[n + 1]();
        G2[i] = new double[n + 1]();
    }
    double* y = new double[n + 1]();

    for (int i = 1; i < n + 1; i++) {
        y[i] = omega * d[i] / b[i];
        if (i < n) {
            C[i][i + 1] = -omega * c[i] / b[i];
        }
        if (i > 1) {
            y[i] += (-omega * a[i] / b[i]) * y[i - 1];
            for (int j = 1; j < n + 1; j++) {
                C[i][j] += (-omega * a[i] / b[i]) * C[i - 1][j];
            }
        }
        C[i][i] += (1.0 - omega);
    }

    //Testing ||C|| < 1
    double norm_C_l1 = matrix_norm_l1(C, n);
    double norm_C_inf = matrix_norm_inf(C, n);

    std::cout << "\nRelaxation method:\n";
    std::cout << "||C||_l1 = " << norm_C_l1 << "\n";
    std::cout << "||C||_inf = " << norm_C_inf << "\n";

    *norm_C = norm_C_inf;

    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            if (i > j) {
                C_L[i][j] = C[i][j];
            }
            else if (i < j) {
                C_U[i][j] = C[i][j];
            }
        }
    }

    double norm_G1_inf = matrix_norm_inf(C_L, n);
    for (int i = 1; i <= n; i++) {
        G2[i][i] = C[i][i];
        for (int j = 1; j <= n; j++) {
            G2[i][j] = C_U[i][j];
        }
    }
    double norm_G2_inf = matrix_norm_inf(G2, n);

    double norm_C_L_inf = matrix_norm_inf(C_L, n);
    double norm_C_U_inf = matrix_norm_inf(C_U, n);
    double sum_norms = norm_C_L_inf + norm_C_U_inf;

    std::cout << "||C_L||_inf = " << norm_C_L_inf << "\n";
    std::cout << "||C_U||_inf = " << norm_C_U_inf << "\n";
    std::cout << "||C_L||_inf + ||C_U||_inf = " << sum_norms << "\n";

    for (int i = 1; i < n + 1; i++) {
        x[i] = y[i];
    }

    double* x1 = new double[n + 1];
    for (int i = 1; i <= n; i++) {
        double sum = d[i];
        if (i > 1) {
            sum -= a[i] * x1[i - 1];
        }
        if (i < n) {
            sum -= c[i] * x[i + 1];
        }
        x1[i] = (1.0 - omega) * x[i] + (omega * sum) / b[i];
    }

    double* rho_vec = new double[n + 1]();
    for (int i = 1; i <= n; i++) {
        rho_vec[i] = x1[i] - x[i];
    }
    double rho_0 = vector_norm_inf(rho_vec, n);

    double k_est = 0.;
    if (sum_norms < 1.0) {
        double q = norm_G2_inf / (1.0 - norm_G1_inf);
        k_est = (log(EPS) + log(1 - q) - log(rho_0)) / log(q);
        std::cout << "Convergence rate estimate (based on sum of norms) >= " << k_est << "\n";
    }
    else {
        std::cout << "Sum of norms >= 1. Convergence is not guaranteed by this criterion.\n";
    }
    delete[] x1;
    delete[] rho_vec;

    if (!flag) {
        const int max_iterations = std::max(1000, (int)(k_est * 2));
        double* x_old = new double[n + 1];
        double* diff = new double[n + 1];
        double diff_norm;
        double criterion;
        for (int step = 0; step < max_iterations; step++) {
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
                x[i] = (1.0 - omega) * x_old[i] + (omega * sum) / b[i];
            }

            //Different stopping criteria 
            double x_old_norm = vector_norm_inf(x_old, n);            //3
            for (int i = 1; i < n + 1; i++) {
                //diff[i] = x[i] - x_old[i];                          //1
                diff[i] = (x[i] - x_old[i]) / (x_old_norm + EPS0);    //3
            }
            diff_norm = vector_norm_inf(diff, n);                 

            //criterion = EPS * ((1 - norm_C_inf) / norm_C_U_inf);    //1
            criterion = EPS;                                          //2, 3

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
                out << "\nerr_norm = " << residual_norm_three_diag(a, b, c, d, x, n);
                out << "\nNorm_l1 C = " << norm_C_l1;
                out << "\nNorm_inf C = " << norm_C_inf;
                out.close();
                for (int i = 1; i <= n; i++) {
                    delete[] C[i];
                    delete[] C_U[i];
                    delete[] C_L[i];
                    delete[] G2[i];
                }
                delete[] C;
                delete[] C_U;
                delete[] C_L;
                delete[] G2;
                delete[] y;
                delete[] x_old;
                delete[] diff;
                return;
            }
        }

        std::cout << "The maximum number of iterations has been reached\n";
        for (int i = 1; i <= n; i++) {
            delete[] C[i];
            delete[] C_U[i];
            delete[] C_L[i];
            delete[] G2[i];
        }
        delete[] C;
        delete[] C_U;
        delete[] C_L;
        delete[] G2;
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
                x[i] = (1.0 - omega) * x_old[i] + (omega * sum) / b[i];
            }
        }

        //Different stopping criteria 
        //double x_old_norm = vector_norm_inf(x_old, n);
        for (int i = 1; i < n + 1; i++) {
            diff[i] = x[i] - x_old[i];
            //diff[i] = (x[i] - x_old[i]) / (x_old_norm + EPS0);
        }
        diff_norm = vector_norm_inf(diff, n);

        criterion = EPS * ((1 - norm_C_inf) / norm_C_U_inf);
        //criterion = EPS * ((1 - norm_C_inf) / norm_C_inf);
        //criterion = EPS;
        //criterion = EPS * x_old_norm + EPS0;

        std::cout << "Method converged in " << floor(k_est + 1) << " iterations\n";
        std::ofstream out(output_filename);
        for (int i = 1; i <= n; i++) {
            out << x[i] << ' ';
        }
        out << "\nk_est >= " << k_est;
        out << "\nsteps = " << floor(k_est + 1);
        out << "\ndiff_norm = " << diff_norm;
        out << "\nerr_norm = " << residual_norm_three_diag(a, b, c, d, x, n);
        out << "\nNorm_l1 C = " << norm_C_l1;
        out << "\nNorm_inf C = " << norm_C_inf;
        out.close();
        for (int i = 1; i <= n; i++) {
            delete[] C[i];
            delete[] C_U[i];
            delete[] C_L[i];
            delete[] G2[i];
        }
        delete[] C;
        delete[] C_U;
        delete[] C_L;
        delete[] G2;
        delete[] y;
        delete[] x_old;
        delete[] diff;
        return;
    }
}