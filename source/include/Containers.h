#pragma once

#include "Types.h"
#include "su4.h"
#include <iostream>

// A structure useful to contain and interchange pairs of
// the std::vector's containing the norms
struct TwoNorms
{
    std::vector<Real> lNorms;
    std::vector<Real> rNorms;
    TwoNorms()
            : lNorms(), rNorms() {}
    TwoNorms(int n)
            : lNorms(n, 0.0), rNorms(n, 0.0) {}
    TwoNorms(int n, Real l, Real r)
            : lNorms(n, l), rNorms(n, r) {}
    TwoNorms(const TwoNorms& tn)
            : lNorms(tn.lNorms), rNorms(tn.rNorms) {}
};

// The objects of this class are basically to contain
// the lines of a solution being calculated; that's much easier to work with such
// objects instead of direct usage of vector<T>'s, because here I aim to implement
// tons of similar operations which will be encountered in the scheme itself.

class TwoLines
{
    // Note: in all the binary operations (e.g. operator+) I admit that
    // the containers must be of the same size. I haven't done an exception for this because it seems
    // there'll be no situations when in my code I have two TwoLines of the different capacities
    public:
        // The default constructor
        TwoLines()  : vleft {}, vright {}, vdim {} {}
        // The constructor which just allocate the memory
        TwoLines(int adim);
        // Initialises with adim copies of el
        TwoLines(int adim, Matrix el);
        // Constructor that gains the initial data
        TwoLines(const std::vector<Matrix>& aleft, const std::vector<Matrix>& aright, int adim);
        // Copy constructor
        TwoLines(const TwoLines& twoLines);
        
        // A way to re-initialise data in the container
        TwoLines& init(const std::vector<Matrix>& aleft, const std::vector<Matrix>& aright, int adim);
        // Copy assignment operator
        TwoLines& operator=(const TwoLines& twoLines);

        // Basic arithmetic operations
        TwoLines operator+(const TwoLines& twoLines) const;
        TwoLines operator-(const TwoLines& twoLines) const;
        TwoLines operator*(Complex z) const;
        TwoLines operator/(Complex z) const;

        // Access methods; each of them has cycled index i
        const Matrix& left(int i) const;
        const Matrix& right(int i) const;

        // Edit methods; each of them has cycled index i
        Matrix& setleft(int i);
        Matrix& setright(int i);

        // Acces to the size
        int dim() const {return vdim;}

        // Finite differences of derivatives
        TwoLines derivative2ndOrder(Real step) const;
        TwoLines derivative4thOrder(Real step) const;
        TwoLines centralDerivative(Real step, int order) const;

        // Elementwise multiplication on the matrices of twoLines
        friend TwoLines dot(const TwoLines& twoLines1, const TwoLines& twoLines2);

        // Some specific ones; this does elementwise normalisation
        // of the matrices, required in the scheme
        void su4Normalise(const TwoNorms& norms);

        // Multiply by \pm
        TwoLines pm() const;

        // Some special directive needed for the enabling additional acceleration
        // through the vectorized operations
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
    private:
        std::vector<Matrix> vleft;
        std::vector<Matrix> vright;
        int vdim;
};

// Right-sided multiplication and division
TwoLines operator*(Complex z, const TwoLines& twoLines);
// Elementwise commutator
TwoLines com(const TwoLines& twoLines1, const TwoLines& twoLines2); // that will actually use a dot friend non-member function

// Short template function of the 2nd derivative of whatever
template <class T>
T difference2ndOrder(const T& next, const T& prev)
{
    return 0.5*next - 0.5*prev;
}

// Short template function of the 4th derivative of whatever
template <class T>
T difference4thOrder(const T& nnext, const T& next, const T& prev, const T& pprev)
{
    return -nnext/12.0 + 2.0*next/3.0 - 2.0*prev/3.0 + pprev/12.0;
}

// STL-like output
std::ostream& operator<<(std::ostream& os, const TwoLines& tl);