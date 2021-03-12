#include "Constants.h"

#include <fstream>
#include <iostream>
#include "nlohmann/json.hpp"

using json = nlohmann::json;
namespace ph = physicalConstants;

#include "su4.h"

// Obtain a single entry of an arbitrary type
template <class T> 
T getEntry(json& j, String section, String key)
{
	return j[section][key].get<T>();
}

// Obtain an entry with floaing point, that doesn't have a unit
Real getDouble(json& j, String section, String key)
{
	// Cast the value to the Complex
	// (Thus I need to explcitly cast the value to Real at first)
	return (Real)j[section][key].get<double>();
}

// Obtain an entry with floating point, and multiply by the corresponding unit
Real getDouble(json& j, String section, String key, std::map<String, Real>& u)
{
	return (Real)(j[section][key][0].get<double>()) * u[j[section][key][1].get<String>()];
}

// Same as above but obtains an array of Real values and casts them to Complex
std::vector<Real> getList(json& j, String section, String key)
{
    std::vector<Real> res;
    for (auto it = j[section][key].begin(); it != j[section][key].end(); ++it)
    {
        res.push_back((Real)(it->get<double>()));
    }
    return res;
}

// Initialise the map of the units
Constants::Constants(String filename)
    : units { {"deg",   ph::deg},
              {"eV",    ph::eV},
              {"eV^2",  ph::eV*ph::eV},
              {"MeV",   ph::MeV},
              {"cm^-3", 1.0/(ph::cm*ph::cm*ph::cm)},
              {"km",    ph::km},
              {"G",     ph::G},
              {"muB",   ph::muB} }
{
    // Firstly, let's save the file into a buffer string
    std::ifstream jsonParameters(filename);
    std::stringstream buf;

    buf << jsonParameters.rdbuf();
    parametersText = buf.str();
    jsonParameters.close();

    // Secondly, let's parse it via nlohmann
    json j = json::parse(parametersText);

    // Then, the dull part of this; let's initialise manuallly all the needed values
    // Here: the basic ones 

	eta     = getDouble(j, "Basic", "eta");
	theta   = getDouble(j, "Basic", "theta", units);
	dm2_0   = getDouble(j, "Basic", "Delta m^2_0", units);
	dm2     = dm2_0 * eta;

	E_0     = getDouble(j, "Basic", "E_0", units);
	chi     = getDouble(j, "Basic", "chi", units);

	// Parameters of MSW potentials and NSSI/V-A couplings

	gMinus  = getDouble(j, "V-A & NSSI & MSW & AMM", "g_{-}");
	gPlus   = getDouble(j, "V-A & NSSI & MSW & AMM", "g_{+}");
	n_Nu    = getDouble(j, "V-A & NSSI & MSW & AMM", "n_Nu", units);
	n_e     = getDouble(j, "V-A & NSSI & MSW & AMM", "n_e", units);
	n_n     = getDouble(j, "V-A & NSSI & MSW & AMM", "n_n", units);
	Beff    = getDouble(j, "V-A & NSSI & MSW & AMM", "Beff", units);
	mu_Nu   = getDouble(j, "V-A & NSSI & MSW & AMM", "mu_Nu", units);

    // Initial parameters

  	InitialProbsLeft  = getList(j,  "Initial Mixes", "InitialProbsLeft");
	InitialProbsRight = getList(j,  "Initial Mixes", "InitialProbsRight");
	MagneticEps       = getDouble(j, "Initial Mixes", "MagneticEps");

	// Parameters of the numerical scheme

	N_z = getEntry<int>(j, "Scheme", "N_z");
	N_x = getEntry<int>(j, "Scheme", "N_x");
	Z   = getDouble(j,     "Scheme", "Z", units);
	X   = getDouble(j,     "Scheme", "X", units);

	WhichNorms         = getEntry<bool>(j, "Scheme", "WhichNorms");
	PressFlag          = getEntry<bool>(j, "Scheme", "PressFlag");
	su4NormaliseFlag   = getEntry<bool>(j, "Scheme", "su4NormaliseFlag");
	xDerivativeOrder   = getEntry<int>(j,  "Scheme", "xDerivativeOrder");
	su4NormalisePeriod = getEntry<int>(j,  "Scheme", "su4NormalisePeriod");

	// Profile

	ProfileFlag   = getEntry<bool>(j, "Profile", "toggle");
	R             = getDouble(j,      "Profile", "z_0", units);

	// Old Fashioned way of perturbing the lumonsity

	OldNoiseFlag = getEntry<bool>(j, "Old Fashioned Noise", "lumPerturbations");
	N_OldNoise   = getEntry<int>(j,  "Old Fashioned Noise", "N_lum");

	lumSinLeft  = getList(j, "Old Fashioned Noise", "lumSinLeft");
	lumCosLeft  = getList(j, "Old Fashioned Noise", "lumCosLeft");
	lumSinRight = getList(j, "Old Fashioned Noise", "lumSinRight");
	lumCosRight = getList(j, "Old Fashioned Noise", "lumCosRight");

	// Derived constants //

	// Cosines of the angles
	cosOmega = std::cos(2.0 * chi);
	cosChi   = std::cos(chi);
	tanChi   = std::tan(chi);

	// Functions of the mixing angles
	c2theta  = std::cos(2.0 * theta);
	s2theta  = std::sin(2.0 * theta);

	// Steps of the grid
	dz = Z / N_z;
	dx = X / N_x;

	// Constants in the collective terms of Hamiltonians
	// Note: all them are in [cm^{-1}]
	V_Nu  = ph::G_F * std::sqrt(2) * n_Nu * (ph::hbarc)*(ph::hbarc);
	V_e   = ph::G_F * std::sqrt(2) * n_e * (ph::hbarc)*(ph::hbarc);
	V_n   = ph::G_F * std::sqrt(2) * n_n * (ph::hbarc)*(ph::hbarc);
	V_AMM = mu_Nu * Beff / (ph::hbarc);

	// Mean initial density matrices
	MeanLeft  = su4::MakeDiagonal(InitialProbsLeft);
	MeanRight = su4::MakeDiagonal(InitialProbsRight);

	// Their norms
	NormMeanLeft  = su4::Norm(MeanLeft);
	NormMeanRight = su4::Norm(MeanRight);
}

// This will be used in the very end of the calculation. This method dumps the parameters of a setup
// alongside the results of the calculations (extremely useful, because it's highly important to see what the parameters
// led to the particular data set)
void Constants::dump(String filename)
{
	std::ofstream out(filename);
	out << parametersText;
	out.close();
}