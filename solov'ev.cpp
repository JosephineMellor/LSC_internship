#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>
#include <tuple>
#include <Eigen/Sparse>

//set up the grid to use for the finite difference scheme
const int Nr = 100; 
const int Nz = 100;
const int N = Nr * Nz;

//set up the constants c0,c1,c2

const double c0 = 1.0;
const double c1 = 0.5;
const double c2 = 0.0;

//set up ranges of r and z and find differences
double r0 = 0.1, r1 = 1.0;
double z0 = 0.0, z1 = 1.0;
double dr = (r1 - r0) / Nr;
double dz = (z1 - z0) / Nz;
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
    A.setZero();
    A.reserve(N*5);
    //loop over the rows 
    for(int row=0; row<N; ++row){
        auto [i,j] = matrixindx(row);
        double r = r0 + i*dr;

        // put in the identitity rows for boundary conditions (every i=0/j=0 and i=Nr/j=Nz make the row the identity), dirichlet bcs
        if(i==0 || j==0 || i== Nr-1 || j==Nz-1){
            A.insert(row,row) = 1;
        }

        // set up the rest of the matrix as per the linear system
        else{
            A.insert(row, row) = 1 / (dz2) + 1 / (dr2); //D4
            A.insert(row, i*Nz + j-1) = ( 1 / (dz2) ) * ( 1 / (dz2)); //1/dz2 D3
            A.insert(row, i*Nz + j+1) = 1 / (dz2); //D3
            A.insert(row, (i-1)*Nz + j) = 1 / (dz2) + 1 / (2*r*dr); //D2
            A.insert(row, (i+1)*Nz + j) = 1 / (dz2) - 1 / (2*r*dr); //D1
        }
    }

    A.makeCompressed();

    return A;
}

//make the b vector on the other side
Eigen::VectorXd defineb(Eigen::VectorXd b){
    b.resize(N);
    for(int row=0; row<N; ++row){
        auto[i,j] = matrixindx(row);
        double r = r0 + i*dr;
        double z = z0 + j*dz;

        //find out if were on the boundary or not
        if(i==0 || j==0 || i== Nr-1 || j==Nz-1){ 
            b(row) = 0.5*(c2*r0*r0 + c0*r*r)*z*z + 0.125*(c1 - c0)*(r*r - r0*r0)*(r*r - r0*r0);
        }
        else{
            b(row) = c1*r*r + c2*r0*r0;
        }
    }
    return b;
}

int main(){

    Eigen::SparseMatrix<double> A;
    A = defineA(A);
    Eigen::VectorXd b;
    b = defineb(b);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1;
    }

    Eigen::VectorXd x = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Solving failed!" << std::endl;
        return -1;
    }

    Eigen::MatrixXd Psi(Nr, Nz);
    for (int i = 0; i < Nr; ++i){
        for (int j = 0; j < Nz; ++j){
            Psi(i, j) = x(vectorindx(i, j));
        }
    }
}














