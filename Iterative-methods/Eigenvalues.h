#pragma once

double powerMethod(double** A, int n, double eps = 1e-6, int max_iter = 1000);
void findEigenvalues(double** A, int n, double& lambda_max, double& lambda_min, double eps = 1e-6);