#include <iostream>

#include "Functions.h"
#include "Eigenvalues"

int main()
{
    Matrix m1;
    Matrix m2;

    Matrix A;
    A << 0.0, -1.0, -2.0, -3.0, 
         1.0,  0.0, -4.0, -5.0,
         2.0,  4.0,  0.0, -6.0,
         3.0,  5.0,  6.0,  0.0;


    Matrix K = Matrix::Identity() + A;
    Matrix Q = (Matrix::Identity() - A) * K.inverse();

    m1 << 1.0, 0.0, 0.0, 0.0,
          0.0, 2.0, 0.0, 0.0,
          0.0, 0.0, 3.0, 0.0,
          0.0, 0.0, 0.0, 4.0;

    m2 << 2.1, 0.0, 0.0, 0.0,
          0.0, 1.2, 0.0, 0.0,
          0.0, 0.0, 3.1, 0.0,
          0.0, 0.0, 0.0, 4.3;

    Matrix D = Q*m2*Q.transpose();
    Matrix R = su4::SetEigenvalues(D, m1);

    Eigen::ComplexEigenSolver<Matrix> cesD(D);
    Eigen::ComplexEigenSolver<Matrix> cesR(R);

    std::cout << "Dirty matrix\n" << D << std::endl << std::endl;
    for (int i=0; i<4; ++i)
    {
        std::cout << "val[" << i << "]: " << cesD.eigenvalues()[i] << "\n"
                  << "vec[" << i << "]: " << std::endl
                  << cesD.eigenvectors().col(i) << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Regularized matrix\n" << R << std::endl;
    for (int i=0; i<4; ++i)
    {
        std::cout << "val[" << i << "]: " << cesR.eigenvalues()[i] << "\n"
                  << "vec[" << i << "]: " << std::endl
                  << cesR.eigenvectors().col(i) << std::endl;
    }

    return 0;
}