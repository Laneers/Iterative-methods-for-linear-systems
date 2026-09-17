#pragma once

double vector_norm_l1(const double* v, const int n);
double matrix_norm_l1(const double* const* A, const int n);
double vector_norm_inf(const double* v, const int n);
double matrix_norm_inf(const double* const* A, const int n);
double residual_norm(const double* const* A, const double* b, const double* x, const int n);