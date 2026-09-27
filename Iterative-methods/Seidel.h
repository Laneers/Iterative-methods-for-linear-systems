#pragma once

void Seidel(const double* a, const double* b, const double* c, const double* d, double* x, const int n, int* converged_step, \
    double* norm_C, const std::string output_filename);