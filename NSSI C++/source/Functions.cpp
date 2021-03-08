#include "Functions.h"

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


// Profiles, when they're needed
inline Real profileNu(Real z, Real z_0)
{
    return std::pow(z_0, 4) / std::pow(z + z_0, 4);
}
inline Real profileMSW(Real z, Real z_0)
{
    return std::pow(z_0, 2) / std::pow(z + z_0, 2);
}
inline Real profileAMM(Real z, Real z_0)
{
    return std::pow(z_0, 2) / std::pow(z + z_0, 2);
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
    return su4::Compose(tmp);
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

// Functions with particular formulae of the matrix flavour Hamiltonians
Matrix VacuumHamiltonian(const Constants& c)
{
    return (c.dm2/(4.0*c.E_0*ph::hbarc)) * (c.s2theta*su4::flavourSigma1 - c.c2theta*su4::flavourSigma3);
}

// MSW part of the flavour effective Hamiltonian
Matrix MSWHamiltonian(const Constants& c, Real z)
{
    Matrix tmp = su4::MakeDiagonal(c.V_e-0.5*c.V_n, -0.5*c.V_n, -c.V_e+0.5*c.V_n, 0.5*c.V_n);
    if (c.ProfileFlag)
    {
        return tmp*profileMSW(z, c.z_0);
    }
    else
    {
        return tmp;
    }
}

// AMM (interaction with magnetic field) part of the flavour effective Hamiltonian
Matrix AMMHamiltonian(const Constants& c, Real z)
{
    if (c.ProfileFlag)
    {
        return su4::AMMLike*c.V_AMM * profileAMM(z, c.z_0);
    }
    else
    {
        return su4::AMMLike*c.V_AMM;
    }
}

// Standard Model collective term
Matrix VAHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z)
{
    if (c.ProfileFlag)
    {
        return 
    }
}

Matrix NSSIHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z)
{
    
}
