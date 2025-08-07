#include <iostream>
#include <Eigen/Sparse>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>

using namespace Eigen;
using namespace std;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

int main() {
    // Physical constants
    double mu0 = 4e-7 * M_PI;
    double B0 = 0.5;   // example magnetic field
    double kappa0 = 2.2;
    double q0 = 1.1;
    double R0 = 0.95;   // major radius

    // Constants from notes
    double c0 = B0 / (R0*R0 * kappa0 * q0);
    double c1 = B0 * (kappa0*kappa0 + 1) / (R0*R0 * kappa0 * q0);
    double c2 = 0.0;

    // Grid parameters
    int NR = 50;
    int NZ = 50;
    double Rmin = 0.2, Rmax = 1.4;
    double Zmin = -1.0, Zmax = 1.0;
    double dR = (Rmax - Rmin) / (NR - 1);
    double dZ = (Zmax - Zmin) / (NZ - 1);

    // Create R grid array
    vector<double> R(NR);
    for (int i = 0; i < NR; i++) {
        R[i] = Rmin + i * dR;
    }

    // Sparse matrix and RHS vector
    int N = NR * NZ;
    SpMat A(N, N);
    VectorXd b(N);
    vector<T> coefficients; // to build sparse matrix efficiently

    // Fill matrix A and vector b
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NZ; ++j) {
            int idx = i * NZ + j; // linear index

            double Ri = Rmin + i * dR;

            // Boundary conditions: for simplicity, Dirichlet Psi=0 at boundaries
            if (i == 0 || i == NR-1 || j == 0 || j == NZ-1) {
                coefficients.push_back(T(idx, idx, 1.0));
                b(idx) = 0.0; // boundary condition value
                continue;
            }

            // Coefficients
            double D1 = 1.0/(dR*dR) - 1.0/(2.0 * Ri * dR);
            double D2 = 1.0/(dR*dR) + 1.0/(2.0 * Ri * dR);
            double D3 = 1.0/(dZ*dZ);
            double D4 = 1.0/(dR*dR) + 1.0/(dZ*dZ);

            // Fill the matrix entries for Psi_{i+1,j}, Psi_{i-1,j}, Psi_{i,j+1}, Psi_{i,j-1}, Psi_{i,j}
            coefficients.push_back(T(idx, (i+1)*NZ + j, D1));
            coefficients.push_back(T(idx, (i-1)*NZ + j, D2));
            coefficients.push_back(T(idx, i*NZ + (j+1), D3));
            coefficients.push_back(T(idx, i*NZ + (j-1), D3));
            coefficients.push_back(T(idx, idx, -D4));  // Note the sign from rearranging

            // RHS from GS equation
            b(idx) = c1 * Ri * Ri + c2 * R0 * R0;
        }
    }

    // Build sparse matrix A
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Solve system Ax = b
    SparseLU<SpMat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    VectorXd x = solver.solve(b);

    Eigen::MatrixXd Psi(NR, NZ);
    for (int i = 0; i < NR; ++i){
        for (int j = 0; j < NZ; ++j){
            Psi(i, j) = x(i * NZ + j);
        }
    }

    std::ofstream output("solovev.dat");

    for (int j = 0; j < NZ; ++j) {
        for (int i = 0; i < NR; ++i) {
            output << R[i] << " " << Zmin + j*dZ << " " << Psi(i,j) << std::endl;
        }
        output<<std::endl;
    }

    return 0;
}
