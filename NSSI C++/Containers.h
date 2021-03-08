#pragma once

#include "Types.h"
#include "su4.h"

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
        TwoLines(int adim) : vleft {adim}, vright {adim}, vdim {adim} {}
        // Constructor that gains the initial data
        TwoLines(std::vector<Matrix> aleft, std::vector<Matrix> aright, int adim)
                    : vleft  {std::move(aleft)},
                      vright {std::move(aright)},
                      vdim   {adim} {}
        // Copy constructor
        TwoLines(const TwoLines& twoLines) 
                    : vleft  {twoLines.vleft},
                      vright {twoLines.vright},
                      vdim   {twoLines.vdim} {}
        // Move constructor
        TwoLines(TwoLines&& twoLines)
                    : vleft  {std::move(twoLines.vleft)},
                      vright {std::move(twoLines.vright)},
                      vdim   {twoLines.vdim} {}
        
        // A way to re-initialise data in the container
        void init(std::vector<Matrix> aleft, std::vector<Matrix> aright, int adim)
        {
            vleft  = std::move(aleft);
            vright = std::move(aright);
            vdim   = adim;
        }
        // Copy assignment operator
        TwoLines(const TwoLines& twoLines)
        {
            vleft  = twoLines.vleft;
            vright = twoLines.vright;
            vdim   = twoLines.vdim;
        }
        // Move assignment operator
        TwoLines(TwoLines&& twoLines)
        {
            vleft  = std::move(twoLines.vleft);
            vright = std::move(twoLines.vright);
            vdim   = twoLines.vdim;
        }

        // Basic arithmetic operations
        TwoLines operator+(const TwoLines& twoLines) const;
        TwoLines operator-(const TwoLines& twoLines) const;
        TwoLines operator*(Complex z) const;
        TwoLines operator/(Complex z) const;

        // Access/edit methods; each of them has cycled index i
        const Matrix& left(int i) const;
        const Matrix& right(int i) const;

        // Finite differences for derivatives
        TwoLines derivative2ndOrder(Real step) const;
        TwoLines derivative4thOrder(Real step) const;
        TwoLines derivative(Real step, int order) const;

        // Elementwise multiplication on the matrices of twoLines
        friend TwoLines dot(const TwoLines& twoLines1, const TwoLines& twoLines2);

        // Some specific ones; this does elementwise normalisation
        // of the matrices, required in the scheme
        void su4Normalise(std::vector<Complex> normsLeft, std::vector<Complex> normsRight);
    
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