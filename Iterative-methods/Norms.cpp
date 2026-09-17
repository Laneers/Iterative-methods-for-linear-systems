#include <cmath>

double vector_norm_l1(const double* v, const int n) {
    double sum = 0;
    for (int i = 1; i <= n; i++) {
        sum += std::fabs(v[i]);
    }
    return sum;
}

double matrix_norm_l1(const double* const* A, const int n) {
    double max_col_sum = 0;
    for (int j = 1; j <= n; j++) {      //идем по столбцам
        double col_sum = 0;
        for (int i = 1; i <= n; i++) {      //cуммируем строки в столбце
            col_sum += std::fabs(A[i][j]);
        }
        if (col_sum > max_col_sum) {
            max_col_sum = col_sum;
        }
    }
    return max_col_sum;
}

double vector_norm_inf(const double* v, const int n) {
    double max_val = 0;
    for (int i = 1; i <= n; i++) {
        if (std::fabs(v[i]) > max_val) {
            max_val = std::fabs(v[i]);
        }
    }
    return max_val;
}

double matrix_norm_inf(const double* const* A, const int n) {
    double max_row_sum = 0;
    for (int i = 1; i <= n; i++) {
        double row_sum = 0;
        for (int j = 1; j <= n; j++) {
            row_sum += std::fabs(A[i][j]);
        }
        if (row_sum > max_row_sum) {
            max_row_sum = row_sum;
        }
    }
    return max_row_sum;
}

double residual_norm(const double* const* A, const double* b, const double* x, const int n) {
    double max_res = 0.0;
    double r;
    for (int i = 1; i < n + 1; i++) {
        double current_row = 0.0;
        for (int j = 1; j < n + 1; j++) {
            current_row += A[i][j] * x[j];
        }
        r = std::fabs(current_row - b[i]);
        if (r > max_res) {
            max_res = r;
        }
    }
    return max_res;
}