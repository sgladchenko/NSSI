#pragma once

#include <map>

#include "Types.h"

// My fancy version of pi

inline constexpr Real pi = 3.1415926535897932384626433832795;

// Definition of the CGS units and derived ones

namespace physicalConstants
{
    inline constexpr Real cm   = 1.0;
    inline constexpr Real g    = 1.0;
    inline constexpr Real sec  = 1.0;
    inline constexpr Real erg  = 1.0;
    inline constexpr Real G    = 1.0;

    inline constexpr Real eV   = 1.602e-12 * erg;
    inline constexpr Real MeV  = 1.0e6 * eV;
    inline constexpr Real km   = 1.0e5 * cm;

    inline constexpr Real muB  = 9.274009994e-21 * erg / G;

    inline constexpr Real G_F  = 1.166e-11  / (MeV*MeV);
    inline constexpr Real hbar = 1.0546e-27 * erg * sec;
    inline constexpr Real c    = 2.9979e10  * cm  / sec;
    inline constexpr Real deg  = pi / 180.0;

    inline constexpr Real hbarc = hbar*c;
}

// This class is needed to obtain all the parameters of the calculations
// from a given JSON file. It's necessary to make the parser same with one that
// I've already implemented in the python code
class Constants
{
    public:
        // The dafualt constructor; actually this is nontrivial. It initialises the inner
        // map of the physical units (this correspondance between the label and
        // final values in CGS is quite needed in the parser)
        // Added: gathering the data from JSON right here
        Constants(String filename);

        // Similarly, but vice versa; just dump the initial text file
        void dump(String filename);

        // Map that labels the units used in JSON file
        std::map<String, Real> units;
        // The text itself of the JSON file
        String parametersText;

        // "Basic"
        Real eta, theta, dm2_0, E_0, chi, dm2;

        // "V-A & NSSI & MSW & AMM"
        Real gMinus, gPlus, n_e, n_n, n_Nu, mu_Nu, Beff;

        // Initial Mixes"
        std::vector<Real> InitialProbsLeft;
        std::vector<Real> InitialProbsRight;
        Real MagneticEps;

        // "Scheme"
        int N_z, N_x, xDerivativeOrder, su4NormalisePeriod, RegEigenvaluesPeriod;
        bool PressFlag, su4NormaliseFlag, WhichNorms, RegEigenvaluesFlag;
        Real Z, X;

        // "Profile"
        bool ProfileFlag;
        Real R;

        // "Old-fashioned noise"
        bool OldNoiseFlag;
        int N_OldNoise;
        std::vector<Real> lumSinLeft;
        std::vector<Real> lumCosLeft;
        std::vector<Real> lumSinRight;
        std::vector<Real> lumCosRight;

        // Derived constants
        Real cosOmega, cosChi, sinChi, tanChi;
        Real c2theta, s2theta;
        Real dz, dx;
        Real V_Nu, V_e, V_n, V_AMM;

        Matrix MeanLeft;
        Matrix MeanRight;

        // su4::Norm of the MeanLeft and MeanRight
        Real NormMeanLeft; Real NormMeanRight;
};