#include <algorithm>

#include "Functions.h"
#include "su4.h"

// Regularization of the eigenvalues
Matrix su4::SetEigenvalues(const Matrix& m, const Matrix& pre)
{
    // First let's solve the eigenvalue problem for the pre matrix
    // and save the solution by pairs value+vector into the std::vector
    Eigen::ComplexEigenSolver<Matrix> preEs(pre);
    std::vector<su4::Eigenpair> preEigs(4);
    for (int k=0; k<4; ++k)
    {
        preEigs[k].vec = preEs.eigenvectors().col(k);
        preEigs[k].val = preEs.eigenvalues()[k].real();
    }

    // Let's do the same but for the numerical approximation saved in m
    Eigen::ComplexEigenSolver<Matrix> mEs(m);
    std::vector<su4::Eigenpair> mEigs(4);
    for (int k=0; k<4; ++k)
    {
        mEigs[k].vec = mEs.eigenvectors().col(k);
        mEigs[k].val = mEs.eigenvalues()[k].real();
    }

    // The sort them in order to find the correspondence between approximate 
    // eigenvalues and true ones
    std::sort(preEigs.begin(), preEigs.end(), compareEigenpairs);
    std::sort(mEigs.begin(),   mEigs.end(), compareEigenpairs);

    // Then let's manually set the numerical approximate eigenvalues by those ones
    // that we see in pre
    for (int k=0; k<4; ++k)
    {
        mEigs[k].val = preEigs[k].val;
    }

    // Finally, let's re-assemble the matrix
    Matrix tmp = Matrix::Zero();
    for (int k=0; k<4; ++k)
    {
        tmp += mEigs[k].val * (mEigs[k].vec * mEigs[k].vec.adjoint());
    }

    return tmp;
}

// Regularization that sets the eigenvalues to ones that are set in the initial density matrices
void RegEigenvalues(const Constants& c, const Noise& n, TwoLines& tl, Real z)
{
    // The containers holding needed matrices
    #pragma omp parallel for
    for (int i=0; i<c.N_x; ++i)
    {
        Real x = i*c.dx;

        // Let's evaluate shifted initial conditions for the beams at this position in x
        Matrix lrho = c.MeanLeft  + FourierVector(originLeft(c,x,z),  c.X, n.sinCoeffsLeft,   n.cosCoeffsLeft);
        Matrix rrho = c.MeanRight + FourierVector(originRight(c,x,z), c.X, n.sinCoeffsRight,  n.cosCoeffsRight);

        tl.setleft(i)  = su4::SetEigenvalues(tl.left(i),  lrho);
        tl.setright(i) = su4::SetEigenvalues(tl.right(i), rrho);
    }
}