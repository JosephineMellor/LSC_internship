#include <iostream>
#include <Eigen/Dense>

int main() {
    // Define 2x2 matrices and vectors
    Eigen::Matrix2d A;
    Eigen::Vector2d b;

    // Initialize values
    A << 2, -1,
         -1, 3;
    b << 1, 2;

    // Solve Ax = b
    Eigen::Vector2d x = A.inverse() * b;

    std::cout << "Solution x:\n" << x << std::endl;

    return 0;
}