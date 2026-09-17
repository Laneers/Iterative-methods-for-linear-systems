#include <iostream>
#include <fstream>
#include <cmath>

#include "Norms.h"
#include "SimpleIteration.h"

double EPS = 1e-4;
//double EPS = 1e-7;

int main() {

    std::string filename;
    std::cin >> filename;

    int n;
    std::ifstream file(filename + ".txt");
    if (!file.is_open()) {
        std::cerr << "Error: failed to open file " << filename + ".txt" << "\n";
        return 0;
    }
    file >> n;

    double** A = new double* [n + 1];
    for (int i = 1; i <= n; i++) {
        A[i] = new double[n + 1];
    }
    double* b = new double[n + 1];
    double* x = new double[n + 1];

    double val;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n + 1; j++) {
            file >> val;
            if (j <= n) {
                A[i][j] = val;
            }
            if (j == n + 1) {
                b[i] = val;
            }
        }
    }
    file.close();

    Simple_iteration(A, b, x, n);

    for (int i = 1; i <= n; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;
    return 0;
}