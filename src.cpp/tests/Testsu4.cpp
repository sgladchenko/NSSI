#include "su4.h"

/*
    TESTED: 9th March
    Seems everything's OK
*/

int main()
{
    su4::Vector vec0({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
    Matrix mat1 = su4::Compose(vec0);
    su4::Vector vec1 = su4::Decompose(mat1);

    Matrix mat2 = su4::Normalise(mat1, 2);
    su4::Vector vec2 = su4::Decompose(mat2);

    Real norm1 = su4::Norm(mat1);
    Real norm2 = su4::Norm(mat2);

    std::vector<Real> v = {1, 2, 3, 4};

    Matrix m;
    m <<  1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12,
         13, 14, 15, 16;

    std::cout << "mat1 = " << std::endl << mat1 << std::endl;
    std::cout << "tr(mat1) = " << mat1.trace() << std::endl;
    std::cout << "vec1 = " << vec1 << std::endl;

    std::cout << "mat2 = " << mat2 << std::endl;
    std::cout << "vec2 = " << vec2 << std::endl;

    std::cout << "norm1 = " << norm1 << std::endl;
    std::cout << "norm2 = " << norm2 << std::endl;

    std::cout << "MakeDiagonal(v) = " << std::endl << su4::MakeDiagonal(v) << std::endl;
    std::cout << "MakeDiagonal(1, 2, 3, 4)" << std::endl << su4::MakeDiagonal(1, 2, 3, 4) << std::endl;

    std::cout << "Diag(m) = " << std::endl << su4::Diag(m) << std::endl; 
    std::cout << "Offdiag(m) = " << std::endl << su4::Offdiag(m) << std::endl;
    std::cout << "Round(m) = " << std::endl << su4::Round(m) << std::endl;

    return 0;
}