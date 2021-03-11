#include <omp.h>

#include "Containers.h"

// Allocate memory
TwoLines::TwoLines(int adim)
         : vleft(adim, Matrix::Zero()), vright(adim, Matrix::Zero()), vdim(adim) {}

// Allocate memory, initialise with copies of el
TwoLines::TwoLines(int adim, Matrix el)
         : vleft(adim, el), vright(adim, el), vdim(adim) {}

// Constructor taking parameters of the container
TwoLines::TwoLines(const std::vector<Matrix>& aleft, const std::vector<Matrix>& aright, int adim)
         : vleft(aleft), vright(aright), vdim(adim) {}

// Copy constructor
TwoLines::TwoLines(const TwoLines& twoLines) 
         : vleft(twoLines.vleft), vright(twoLines.vright), vdim(twoLines.vdim) {}

// A way to re-initialise data in the container
TwoLines& TwoLines::init(const std::vector<Matrix>& aleft, const std::vector<Matrix>& aright, int adim)
{
    // Let's first free the memory
    vleft.clear();
    vright.clear();
    // Then, let's make a copy assignments for the vectors
    vleft  = aleft;
    vright = aright;
    vdim   = adim;
    return *this;
}

// Copy assignment operator
TwoLines& TwoLines::operator=(const TwoLines& twoLines)
{
    // Let's first free the memory
    vleft.clear();
    vright.clear();
    // Then, let's make a copy assignments for the vectors
    vleft  = twoLines.vleft;
    vright = twoLines.vright;
    vdim   = twoLines.vdim;
    return *this;
}

// Elementwise summation of the density matrices
TwoLines TwoLines::operator+(const TwoLines& twoLines) const
{
    TwoLines tmp(vleft, vright, vdim);
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        tmp.vleft[i]  += twoLines.vleft[i];
        tmp.vright[i] += twoLines.vright[i];
    }
    return tmp;
}

// Elementwise subtraction of the density matrices
TwoLines TwoLines::operator-(const TwoLines& twoLines) const
{
    TwoLines tmp(vleft, vright, vdim);
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        tmp.vleft[i]  -= twoLines.vleft[i];
        tmp.vright[i] -= twoLines.vright[i];
    }
    return tmp;
}

// Elementwise multiplication on a complex number of the density matrices
TwoLines TwoLines::operator*(Complex z) const
{
    TwoLines tmp(vleft, vright, vdim);
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        tmp.vleft[i]  *= z;
        tmp.vright[i] *= z;
    }
    return tmp;
}

// Elementwise division by a complex number of the density matrices
TwoLines TwoLines::operator/(Complex z) const
{
    TwoLines tmp(vleft, vright, vdim);
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        tmp.vleft[i]  /= z;
        tmp.vright[i] /= z;
    }
    return tmp;
}

// Right-sided elementwise multiplication
TwoLines operator*(Complex z, const TwoLines& twoLines)
{
    return twoLines * z;
}

// Access a matrix of the left beam
const Matrix& TwoLines::left(int i) const
{
    return vleft[( (i % vdim) + vdim) % vdim];
}

// Access a matrix of the right beam
const Matrix& TwoLines::right(int i) const
{
    return vright[( (i % vdim) + vdim) % vdim];
}

// Edit an element in the left beam
Matrix& TwoLines::setleft(int i)
{
    return vleft[( (i % vdim) + vdim) % vdim];
}
// Edit an element in the right beam
Matrix& TwoLines::setright(int i)
{
    return vright[( (i % vdim) + vdim) % vdim];
}

TwoLines TwoLines::derivative2ndOrder(Real step) const
{
    TwoLines tmp(vdim);
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        tmp.vleft[i]  = difference2ndOrder(this->left(i+1),  this->left(i-1)) / step;
        tmp.vright[i] = difference2ndOrder(this->right(i+1), this->right(i-1)) / step;
    }
    return tmp;
}

TwoLines TwoLines::derivative4thOrder(Real step) const
{
    TwoLines tmp(vdim);
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        tmp.vleft[i]  = difference4thOrder(this->left(i+2),  this->left(i+1),  this->left(i-1),  this->left(i-2)) / step;
        tmp.vright[i] = difference4thOrder(this->right(i+2), this->right(i+1), this->right(i-1), this->right(i-2)) / step;
    }
    return tmp;
}

TwoLines TwoLines::derivative(Real step, int order) const
{
    if (order == 2)
    {
        return derivative2ndOrder(step);
    }
    else
    {
        return derivative4thOrder(step);
    }
}

// Elementwise multiplication on the matrices of twoLines
TwoLines dot(const TwoLines& twoLines1, const TwoLines& twoLines2)
{
    TwoLines tmp(twoLines1.vdim);
    #pragma omp parallel for
    for (int i = 0; i < twoLines1.vdim; ++i)
    {
        tmp.vleft[i]  = twoLines1.left(i)  * twoLines2.left(i);
        tmp.vright[i] = twoLines1.right(i) * twoLines2.right(i);
    }
    return tmp;
}

// Elementwise commutator
TwoLines com(const TwoLines& twoLines1, const TwoLines& twoLines2)
{
    return dot(twoLines1, twoLines2) - dot(twoLines2, twoLines1);
}

// Elementwise su4Normalise
void TwoLines::su4Normalise(const TwoNorms& norms)
{
    #pragma omp parallel for
    for (int i = 0; i < vdim; ++i)
    {
        vleft[i]  = su4::Normalise(vleft[i],  norms.lNorms[i]);
        vright[i] = su4::Normalise(vright[i], norms.rNorms[i]);
    }
}

// STL-like output
std::ostream& operator<<(std::ostream& os, const TwoLines& tl)
{
    // Output the left line
    for (int i=0; i < tl.dim(); ++i)
    {
        os << "left[" << i << "]: " << std::endl << tl.left(i) << std::endl;
    }
    // Output the right line
    for (int i=0; i < tl.dim(); ++i)
    {
        os << "right[" << i << "]: " << std::endl << tl.right(i) << std::endl;
    }
    return os;
}