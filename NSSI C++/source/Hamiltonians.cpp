#include "Hamiltonians.h"

namespace ph = physicalConstants;

// Profiles, when they're needed
inline Real profileNu(Real z, Real R)
{
    return std::pow(R, 4) / std::pow(z + R, 4);
}
inline Real profileMSW(Real z, Real R)
{
    return std::pow(R, 2) / std::pow(z + R, 2);
}
inline Real profileAMM(Real z, Real R)
{
    return std::pow(R, 2) / std::pow(z + R, 2);
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
        return tmp*profileMSW(z, c.R);
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
        return su4::AMMLike * (profileAMM(z, c.R)*c.V_AMM);
    }
    else
    {
        return su4::AMMLike * c.V_AMM;
    }
}

// Standard Model collective term
Matrix VAHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z)
{
    Matrix rhoG = rhoOpp * su4::G;
    Matrix tmp  = rhoG.trace()*su4::G + su4::Diag(su4::Round(rhoOpp));

    if (c.ProfileFlag)
    {
        return tmp * (profileNu(z, c.R)*c.V_Nu*(1 - c.cosOmega));
    }
    else
    {
        return tmp * (c.V_Nu*(1 - c.cosOmega));
    }
}

// The collective term that describes coupling through the scalar and 
// pseudoscalar fields
Matrix NSSIHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z)
{
    Matrix diagRound    = su4::Diag(su4::Round(rhoOpp));
    Matrix offdiagRound = su4::Offdiag(su4::Round(rhoOpp));
    Matrix tmp = c.gMinus*diagRound.transpose() + c.gPlus*offdiagRound.transpose();

    if (c.ProfileFlag)
    {
        return tmp * (profileNu(z, c.R)*c.V_Nu*(1 - c.cosOmega));
    }
    else
    {
        return tmp * (c.V_Nu*(1 - c.cosOmega));
    }
}