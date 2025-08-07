#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>
#include <tuple>
#include <Eigen/Sparse>

//set up the grid to use for the finite difference scheme
const int Nr = 800; 
const int Nz = 800;
const int N = Nr * Nz;

//set up the constants c0,c1,c2, note that these are non-dimensionalised

const double c2 = 0.0;
const double kappa0 = 2.2;
const double B0 = 0.5;
const double R0 = 0.95;
const double q0 = 1.1;
const double c0 = B0 / (R0*R0* kappa0 * q0);
const double c1 = B0*(kappa0*kappa0 + 1) / (R0*R0*kappa0 * q0);

//set up ranges of r and z and find differences
double r0 = 0.2, r1 = 1.4;
double z0 = -1.0, z1 = 1.0;
double dr = (r1 - r0) / (Nr - 1);
double dz = (z1 - z0) / (Nz - 1);
double dr2 = dr * dr;
double dz2 = dz * dz;

//flatten the (i,j) into the index of a vector of size N (there will be Nz of each i with each j added on)
int vectorindx(int i, int j){
    return i * Nz + j;
}

//want a function that goes the other way as well for setting up my matrix
std::tuple<int,int> matrixindx(int k){
    int i = floor(k / Nz);
    int j = k - i * Nz;
    return {i,j};
}

//make the sparse matrix
Eigen::SparseMatrix<double> defineA(Eigen::SparseMatrix<double>& A){
    A.resize(N,N);
    A.reserve(N*5);
    //loop over the rows 
    for(int row=0; row<N; ++row){
        auto [i,j] = matrixindx(row);
        double r = r0 + i*dr;

        if(i==0 || j==0 || i==Nr-1 || j==Nz-1){
            A.insert(row,row) = 1.0;  // Dirichlet BC
        } else {
            double D1 = 1.0 / dr2 - 1.0 / (2.0 * r * dr);
            double D2 = 1.0 / dr2 + 1.0 / (2.0 * r * dr);
            double D3 = 1.0 / dz2;
            double D4 = 2.0 / dr2 + 2.0 / dz2;

            A.insert(row, row) = -D4;
            A.insert(row, vectorindx(i, j-1)) = D3;
            A.insert(row, vectorindx(i, j+1)) = D3;
            A.insert(row, vectorindx(i-1, j)) = D2;
            A.insert(row, vectorindx(i+1, j)) = D1;
        }
    }

    A.makeCompressed();

    return A;
}

Eigen::SparseMatrix<double> defineAquicker(Eigen::SparseMatrix<double>& A){
    A.resize(N,N);
    A.reserve(N*5);

    std::vector<Eigen::Triplet<double>> coefficients;
    int idx;
    double Ri,D1,D2,D3,D4;

    for (int i = 0; i < Nr; ++i) {
        for (int j = 0; j < Nz; ++j) {
            idx = i * Nz + j; 

            Ri = r0 + i*dr;

            //put in the boundary conditions
            if (i == 0 || i == Nr-1 || j == 0 || j == Nz-1) {
                coefficients.push_back(Eigen::Triplet<double>(idx, idx, 1.0));
                continue;
            }

            // coefficients
            D1 = 1.0/(dr*dr) - 1.0/(2.0 * Ri * dr);
            D2 = 1.0/(dr*dr) + 1.0/(2.0 * Ri * dr);
            D3 = 1.0/(dz*dz);
            D4 = 2.0/(dr*dr) + 2.0/(dz*dz);

            //matrix valus for Psii+1, Psii-1, Psij+1, Psij+1, Psiij
            coefficients.push_back(Eigen::Triplet<double>(idx, (i+1)*Nz + j, D1));
            coefficients.push_back(Eigen::Triplet<double>(idx, (i-1)*Nz + j, D2));
            coefficients.push_back(Eigen::Triplet<double>(idx, i*Nz + (j+1), D3));
            coefficients.push_back(Eigen::Triplet<double>(idx, i*Nz + (j-1), D3));
            coefficients.push_back(Eigen::Triplet<double>(idx, idx, -D4));  //sign from rearranging?

        }
    }

    // Build sparse matrix A
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

//make the b vector on the other side
Eigen::VectorXd defineb(Eigen::VectorXd b){
    b.resize(N);
    
    int idx;
    double Ri;

    for (int i = 0; i < Nr; ++i) {
        for (int j = 0; j < Nz; ++j) {
            idx = i * Nz + j; 

            Ri = r0 + i*dr;
            double z = z0 + j*dz;

            if (i == 0 || i == Nr-1 || j == 0 || j == Nz-1) {
                
                b(idx) = 0.5*(c2 + c0*Ri*Ri)*z*z + 0.125*(c1 - c0)*(Ri*Ri - 1)*(Ri*Ri - 1); // boundary condition value for dirichlet
                continue;
            }

            b(idx) = c1 * Ri * Ri + c2 * R0 * R0;
        }
    }
    return b;
}

//going to write a function for the exact solution
Eigen::VectorXd defineexact(Eigen::VectorXd b){
    b.resize(N);
    for(int row=0; row<N; ++row){
        auto[i,j] = matrixindx(row);
        double r = r0 + i*dr;
        double z = z0 + j*dz;

        b(row) = 0.5*(c2 + c0*r*r)*z*z + 0.125*(c1 - c0)*(r*r - 1)*(r*r - 1);
        
    }
    return b;
}


int main(){
    

    Eigen::SparseMatrix<double> A;
 
    A = defineA(A);
    Eigen::VectorXd b;
  
    b = defineb(b);
    Eigen::VectorXd exact;

    exact = defineexact(exact);

    // Solve system Ax = b
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Eigen::VectorXd x = solver.solve(b);

    Eigen::MatrixXd Psi(Nr, Nz);
    for (int i = 0; i < Nr; ++i){
        for (int j = 0; j < Nz; ++j){
            Psi(i, j) = x(i*Nz +j);
        }
    }

    std::ofstream output("solovev1.dat");

    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Nr; ++i) {
            output << r0 + i*dr << " " << z0 + j*dz << " " << Psi(i,j) << std::endl;
        }
        output<<std::endl;
    }
}
