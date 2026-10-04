#pragma once

void Relaxation_triag(const double* a, const double* b, const double* c, const double* d, double* x, const int n, const double omega, int* converged_step, \
    double* norm_C, const std::string output_filename, bool flag = 0);