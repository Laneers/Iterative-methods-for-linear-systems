#pragma once

void Simple_iteration(const double* const* A, const double* b, double* x, const int n, const double tau, int* converged_step, \
	double* norm_C, const std::string output_filename);