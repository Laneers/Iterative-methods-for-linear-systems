#pragma once

void Jacobi(const double* const* A, const double* b, double* x, const int n, int* converged_step, \
    double* norm_C, const std::string output_filename);