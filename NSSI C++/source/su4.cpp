#include "su4.h"
#include <iostream>

// This functions evaluates the decomposition of a given matrix
// into the su4Basis described in su4.h
su4::Vector su4::Decompose(const Matrix& m)
{
    Array15 tmpArr;
    Matrix tmp;
    Complex tr;
    for (int k = 0; k < 15; ++k)
    {
        tmp = m*su4::Basis[k];
        tr = tmp.trace();
        tmpArr[k] = tr.real(); // I presume that within the su(4) algebra those numbers must be pure real
        // Thus I can just omit the imaginary part; for an abritrary matrix this leads to unpredicatable values
    }
    return su4::Vector(tmpArr);
}

// The inverse operations; now I want to just evaluate the matrix from
// a given vector of coefficients
Matrix su4::Compose(const su4::Vector& v)
{
    Matrix res;
    res = Matrix::Zero();
    for (int k = 0; k < 15; ++k)
    {   
        res += v[k]*su4::Basis[k];
    }
    return res + 0.25*Matrix::Identity();
}

// Similar as above, but doesn't add an additional 0.25*Id
Matrix su4::ComposeTraceless(const su4::Vector& v)
{
    Matrix res;
    res = Matrix::Zero();
    for (int k = 0; k < 15; ++k)
    {   
        res += v[k]*su4::Basis[k];
    }
    return res;
}

// This operations provides a way to regularise the matrix and put it on
// the su(4) Lie algebra
Matrix su4::Normalise(const Matrix& m, Real initialNorm)
{
    // First let's evaluate the traceless hermitian part of the matrix
    Matrix hermitian; hermitian = 0.5*(m + m.adjoint());
    Matrix traceless; traceless = hermitian - 0.25*Matrix::Identity()*hermitian.trace();
    // Then finally  let's re-scale the norm
    Real norm = traceless.norm();
    if (norm == 0.0) // Formally, that may happen
    {
        return traceless + 0.25*Matrix::Identity();
    }
    else
    {
        return traceless*(initialNorm/norm) + 0.25*Matrix::Identity();
    }
}

// This function wipes the offdiagonal blocks from a given matrix
Matrix su4::Diag(const Matrix& m)
{
    Matrix res;
    res.topLeftCorner(2,2) = m.topLeftCorner(2,2);
    res.topRightCorner(2,2).setZero();
    res.bottomLeftCorner(2,2).setZero();
    res.bottomRightCorner(2,2) = m.bottomRightCorner(2,2);
    return res;
}

// This function wipes the diagonal blocks from a given matrix
Matrix su4::Offdiag(const Matrix& m)
{
    Matrix res;
    res.topLeftCorner(2,2).setZero();
    res.topRightCorner(2,2) = m.topRightCorner(2,2);
    res.bottomLeftCorner(2,2) = m.bottomLeftCorner(2,2);
    res.bottomRightCorner(2,2).setZero();
    return res;
}

// Some special operation needed for the construction of a NSSI Ham.
Matrix su4::Round(const Matrix& m)
{
    Matrix tmp;
    tmp.topLeftCorner(2,2)     = m.bottomRightCorner(2,2).transpose();
    tmp.topRightCorner(2,2)    = m.topRightCorner(2,2).transpose();
    tmp.bottomLeftCorner(2,2)  = m.bottomLeftCorner(2,2).transpose();
    tmp.bottomRightCorner(2,2) = m.topLeftCorner(2,2).transpose();
    return m - tmp;
}

// Construct a diagonal matrix
Matrix su4::MakeDiagonal(const std::vector<Real>& v)
{
    Matrix m;
    m << v[0], 0.0,  0.0,  0.0,
         0.0,  v[1], 0.0,  0.0,
         0.0,  0.0,  v[2], 0.0,
         0.0,  0.0,  0.0,  v[3];
    return m;
}
Matrix su4::MakeDiagonal(Real a, Real b, Real c, Real d)
{
    Matrix m;
    m <<   a,  0.0,  0.0,  0.0,
         0.0,    b,  0.0,  0.0,
         0.0,  0.0,    c,  0.0,
         0.0,  0.0,  0.0,    d;
    return m;
}

// Evaluates the Frobenius norm of the traceless part of a matrix
// whhich is nothing but the length of the corresponding 15D-vector
Real su4::Norm(const Matrix& m)
{
    Matrix traceless = m - 0.25*Matrix::Identity()*m.trace();
    return traceless.norm();
}

bool compareEigenpairs(const su4::Eigenpair& e1, const su4::Eigenpair& e2)
{
    return e1.val < e2.val;
}
