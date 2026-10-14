#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Norms.h"

double powerMethod(double** A, int n, double eps = 1e-6, int max_iter = 1000) {
    double* phi = new double[n + 1];
    double* next_phi = new double[n + 1];

    for (int i = 1; i <= n; i++) {
        phi[i] = 1.0;
    }

    double next_lambda = 0.0;
    double lambda = 0.0;

    for (int step = 0; step < max_iter; step++) {
        double norm_phi = vector_norm_2(phi, n);
        if (norm_phi < 1e-12) break;

        //next_phi = A * phi / ||phi||
        for (int i = 1; i <= n; i++) {
            double sum = 0.0;
            for (int j = 1; j <= n; j++) {
                sum += A[i][j] * phi[j];
            }
            next_phi[i] = sum / norm_phi;
        }

        next_lambda = vector_norm_2(next_phi, n);

        if (step > 0 && std::abs(next_lambda - lambda) < eps) {
            lambda = next_lambda;
            break;
        }

        for (int i = 1; i <= n; i++) {
            phi[i] = next_phi[i];
        }
    }
    delete[] phi;
    delete[] next_phi;

    return lambda;
}

void findEigenvalues(double** A, int n, double& lambda_max, double& lambda_min, double eps = 1e-6) {
    double shift = matrix_norm_inf(A, n) + 1.0;

    //A_shifted = A + shift * E
    double** A_shifted = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        A_shifted[i] = new double[n + 1];
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                A_shifted[i][j] = A[i][j] + shift;
            }
            else {
                A_shifted[i][j] = A[i][j];
            }
        }
    }

    double lambda_max_shifted = powerMethod(A_shifted, n);

    //B = lambda_max_shifted * E - A_shifted to find the minimum shifted eigenvalue 
    //lambda_min_shifted = lambda_max_shifted - power_method_max(B, n)
    double** B = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        B[i] = new double[n + 1];
        for (int j = 1; j <= n; j++) {
            if (i == j) {
                B[i][j] = lambda_max_shifted - A_shifted[i][j];
            }
            else {
                B[i][j] = -A_shifted[i][j];
            }
        }
    }

    double lambda_max_B = powerMethod(B, n, eps);

    double lambda_min_shifted = lambda_max_shifted - lambda_max_B;

	//Shift back to the original matrix
    lambda_max = lambda_max_shifted - shift;
    lambda_min = lambda_min_shifted - shift;

    for (int i = 1; i <= n; i++) {
        delete[] A_shifted[i];
        delete[] B[i];
    }
    delete[] A_shifted;
    delete[] B;
}