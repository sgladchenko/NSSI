#include <omp.h>

#include "Functions.h"
#include "Hamiltonians.h"

namespace ph = physicalConstants;
namespace hs = hamiltonians;

// Origin of the left beam going thorugh the point (x,z)
inline Real originLeft(const Constants& c, Real x, Real z)
{
    return x + z*c.tanChi;
}

// Origin of the left beam going thorugh the point (x,z)
inline Real originRight(const Constants& c, Real x, Real z)
{
    return x - z*c.tanChi;
}

// X-periodic sin function
inline Real Sin(Real x, Real X, int k)
{
    return std::sin(2.0*pi*k/X * x);
}
// X-periodic cos function
inline Real Cos(Real x, Real X, int k)
{
    return std::cos(2.0*pi*k/X * x);
}

// One specific function function which linearly interpolates
// a given discrete grid at some given point x
Real linearInterpolate(const Constants& c, Real x, const std::vector<Real>& vec)
{
	// First of all, let's measure the coordinates in the numbers of dx
	Real scaled = x / c.dx;
	// Then let's evaluate the neighbour points of the grid
	// and interpolate the value between them
	int lNode = std::floor(scaled);
	int rNode = lNode + 1;

	Real lWeight = (rNode - scaled); // / (rNode - lNode), but rNode - lNode = 1
	Real rWeight = (scaled - lNode);

	return lWeight*vec[((lNode % c.N_x) + c.N_x) % c.N_x] +
           rWeight*vec[((rNode % c.N_x) + c.N_x) % c.N_x];
}

// Evaluates a sum of Fourier harmonics of su4::Vector's
Matrix FourierVector(Real x, Real X, const std::vector<su4::Vector>& sinVecs,
                                     const std::vector<su4::Vector>& cosVecs)
{
    su4::Vector tmp;
    for (int k = 0; k < sinVecs.size(); ++k)
    {
        tmp += Sin(x, X, k+1)*sinVecs[k] + Cos(x, X, k+1)*cosVecs[k];
    }
    return su4::ComposeTraceless(tmp);
}

// Evaluates a sum of Fourier harmonics in the luminosity 
Real FourierLum(Real x, Real X, const std::vector<Real>& sinVecs,
                                const std::vector<Real>& cosVecs)
{
    Real tmp;
    for (int k = 0; k < sinVecs.size(); ++k)
    {
        tmp += Sin(x, X, k+1)*sinVecs[k] + Cos(x, X, k+1)*cosVecs[k];
    }
    return 1.0 + tmp;
}

// Evaluates and assembles TwoLines of the initial density matrices
TwoLines InitialConditions(const Constants& c, const Noise& n)
{
    TwoLines tmp(c.N_x);
    #pragma omp parallel for
    for (int i = 0; i < c.N_x; ++i)
    {
        Real x = i * c.dx;
        tmp.setleft(i)  = c.MeanLeft  + FourierVector(x, c.X, n.sinCoeffsLeft,  n.cosCoeffsLeft);
        tmp.setright(i) = c.MeanRight + FourierVector(x, c.X, n.sinCoeffsRight, n.cosCoeffsRight);

        if (c.PressFlag) // If it's needed that all the vectors lie on the same S^14
        {
            tmp.setleft(i)  = su4::Normalise(tmp.left(i),  c.NormMeanLeft);
            tmp.setright(i) = su4::Normalise(tmp.right(i), c.NormMeanRight);
        }
    }

    return tmp;
}

// The functions that generate the total Hamilatonians in a TwoLines instance
TwoLines Hamiltonians(const Constants& c, const TwoLines& rho, Real z)
{
    TwoLines tmp(c.N_x);
    #pragma omp parallel for
    for (int i=0; i<c.N_x; ++i)
    {
        Real x = i * c.dx;
        Real lumLeft  = FourierLum(originLeft(c, x, z),  c.X, c.lumSinLeft,  c.lumCosLeft);
        Real lumRight = FourierLum(originRight(c, x, z), c.X, c.lumSinRight, c.lumCosRight);
        // Note that opposite to the left is right
        tmp.setleft(i) = hs::VacuumHamiltonian(c) +
                         hs::MSWHamiltonian(c, z) +
                         hs::AMMHamiltonian(c, z) +
                         hs::VAHamiltonian(c, rho.right(i), z, lumRight) +
                         hs::NSSIHamiltonian(c, rho.right(i), z, lumRight);
        // Note that opposite to the right is left
        tmp.setright(i) = hs::VacuumHamiltonian(c) +
                          hs::MSWHamiltonian(c, z) +
                          hs::AMMHamiltonian(c, z) + 
                          hs::VAHamiltonian(c, rho.left(i), z, lumLeft) +
                          hs::NSSIHamiltonian(c, rho.left(i), z, lumLeft);
    }
    return tmp;
}

// Function that evaluates the norms of the initial conditions
TwoNorms InitialNorms(const Constants& c, const TwoLines& rho)
{
    // PressFlag means either the i.c. have same value of the su4::Norm
    // or not; if they have, TwoNorms has to be initialized with same values
    // given in c.NormMeanLeft and c.NormMeanRight
    if (c.PressFlag)
    {
        return TwoNorms(c.N_x, c.NormMeanLeft, c.NormMeanRight);
    }
    else
    {
        // Otherwise we have to calculate the directly ->
        TwoNorms tmp(c.N_x);
        #pragma omp parallel for
        for (int i=0; i<c.N_x; ++i)
        {
            tmp.lNorms[i] = su4::Norm(rho.left(i));
            tmp.rNorms[i] = su4::Norm(rho.right(i));
        }
        return tmp;
    }
}

// Function that approximately calculates the norms at a given position z
// (approximaely means that it uses linear interpolation of the initial norms
// in order to figure the norms between the nodes of the grid)
TwoNorms ApproximateNorms(const Constants& c, const TwoNorms& initNorms, Real z)
{
    // PressFlag means either the i.c. have same value of the su4::Norm
    // or not; if they have, TwoNorms has to be initialized with same values
    // given in c.NormMeanLeft and c.NormMeanRight, as well as in the InitialNorms
    if (c.PressFlag)
    {
        TwoNorms tmp(initNorms);
        return tmp;
    }
    else
    {
        // Here we want to obtain the norms at the position z, but we want to calculate them apprxoimately.
        // We already know the the initial norms at the points of the x grid,
        // and norms at the position z are same but shifted at \pm z * tahChi for left and right beams respectively.
        // Thus we need to know the values of norms between the points of the grids
        // and here we aim to obtain them just interpolating the values given on the main grid of x axis,
        // where the norms are known.
        TwoNorms tmp(c.N_x);
        #pragma omp parallel for
        for (int i=0; i<c.N_x; ++i)
        {
            Real x = i*c.dx;
            tmp.lNorms[i] = linearInterpolate(c, originLeft(c,x,z),  initNorms.lNorms);
            tmp.rNorms[i] = linearInterpolate(c, originRight(c,x,z), initNorms.rNorms);
        }
        return tmp;
    }
}

// Same as above, but it actually calculates them directly
TwoNorms TrueNorms(const Constants& c, const Noise& n, const TwoNorms& initNorms, Real z)
{
    // PressFlag means either the i.c. have same value of the su4::Norm
    // or not; if they have, TwoNorms has to be initialized with same values
    // given in c.NormMeanLeft and c.NormMeanRight, as well as in the InitialNorms
    if (c.PressFlag)
    {
        TwoNorms tmp(initNorms);
        return tmp;
    }
    else
    {
        // Else, we need to recalculate the values at the points of the grid;
        TwoNorms tmp(c.N_x);
        #pragma omp parallel for
        for (int i=0; i<c.N_x; ++i)
        {
            Real x = i*c.dx;
            // Let's evaluate frst the initial conditions
            Matrix lrho = c.MeanLeft  + FourierVector(originLeft(c,x,z),  c.X, n.sinCoeffsLeft,   n.cosCoeffsLeft);
            Matrix rrho = c.MeanRight + FourierVector(originRight(c,x,z), c.X, n.sinCoeffsRight,  n.cosCoeffsRight);
            // Then we need their su4::Norm's
            tmp.lNorms[i] = su4::Norm(lrho);
            tmp.rNorms[i] = su4::Norm(rrho);
        }
        return tmp;
    }
}