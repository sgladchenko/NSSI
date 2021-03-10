#include "Functions.h"
#include "Hamiltonians.h"

namespace ph = physicalConstants;

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

// Evaluates a sum of Fourier harmonics
Matrix FourierContract(Real x, Real X, const std::vector<su4::Vector>& sinVecs,
                                       const std::vector<su4::Vector>& cosVecs)
{
    su4::Vector tmp;
    for (int k = 0; k < sinVecs.size(); ++k)
    {
        tmp += Sin(x, X, k+1)*sinVecs[k] + Cos(x, X, k+1)*cosVecs[k];
    }
    return su4::ComposeTraceless(tmp);
}

// Evaluates and assembles TwoLines of the initial density matrices
TwoLines InitialConditions(const Constants& c, const Noise& n)
{
    std::vector<Matrix> leftRhos(c.N_x);
    std::vector<Matrix> rightRhos(c.N_x);
    Real x;

    for (int i = 0; i < c.N_x; ++i)
    {
        x = i * c.dx;
        leftRhos[i]  = c.MeanLeft  + FourierContract(x, c.X, n.sinCoeffsLeft,  n.cosCoeffsLeft);
        rightRhos[i] = c.MeanRight + FourierContract(x, c.X, n.sinCoeffsRight, n.cosCoeffsRight);

        if (c.PressFlag) // If it's needed that all the vectors lie on the same S^14
        {
            leftRhos[i]  = su4::Normalise(leftRhos[i],  c.NormMeanLeft);
            rightRhos[i] = su4::Normalise(rightRhos[i], c.NormMeanRight);
        }
    }

    return TwoLines(leftRhos, rightRhos, c.N_x);
}

// The functions that generate the total Hamilatonians in a TwoLines instance
TwoLines Hamiltonians(const Constants& c, const TwoLines& rho, Real z)
{
    
}
