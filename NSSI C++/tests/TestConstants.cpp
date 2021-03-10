#include "Constants.h"

/*
	TESTED: 9th March. Seems everything's OK
*/

#include <iostream>

namespace ph = physicalConstants;

// Some output of the initialised values
void BaseConstants(Constants& c)
{
	std::cout << "> Basic:"    << std::endl;
	std::cout << "\t eta   = " << c.eta << std::endl
			  << "\t theta = " << c.theta/ph::deg << " deg" << std::endl
			  << "\t dm2_0 = " << c.dm2_0/ph::eV/ph::eV << " eV^2" << std::endl
			  << "\t E_0 = "   << c.E_0/ph::MeV << " MeV" << std::endl
			  << "\t chi = "   << c.chi/ph::deg << " deg" << std::endl;

	std::cout << "> Potentials:" << std::endl;
	std::cout << "\t gMinus = "  << c.gMinus << std::endl
			  << "\t gPlus  = "  << c.gPlus  << std::endl
			  << "\t n_Nu   = "  << c.n_Nu/(ph::cm*ph::cm*ph::cm) << " cm^{-3}" << std::endl
			  << "\t n_e   = "   << c.n_e/(ph::cm*ph::cm*ph::cm) << " cm^{-3}" << std::endl
			  << "\t n_n   = "   << c.n_n/(ph::cm*ph::cm*ph::cm) << " cm^{-3}" << std::endl
			  << "\t Beff = "    << c.Beff/ph::G << " G" << std::endl
			  << "\t mu_Nu = "   << c.mu_Nu/ph::muB << " muB" << std::endl;

	std::cout << "> Initial parameters" << std::endl;
	std::cout << "\t InitialProbsLeft = [";
	for (auto it = c.InitialProbsLeft.begin(); it != c.InitialProbsLeft.end(); ++it)
		std::cout << *it << " ";
	std::cout << "]" << std::endl;

	std::cout << "\t InitialProbsRight = [";
	for (auto it = c.InitialProbsRight.begin(); it != c.InitialProbsRight.end(); ++it)
		std::cout << *it << " ";
	std::cout << "]" << std::endl;

	std::cout << "\t MagenticEps = " << c.MagneticEps << std::endl;

	std::cout << "> Scheme parameters:" << std::endl;
	std::cout << "\t N_x = "                << c.N_x                << std::endl
			  << "\t N_z  = "               << c.N_z                << std::endl
			  << "\t X   = "                << c.X/ph::km << " km"  << std::endl
			  << "\t Z   = "                << c.Z/ph::km << " km"  << std::endl
			  << "\t su4NormalisePeriod = " << c.su4NormalisePeriod << std::endl
			  << "\t su4NormaliseFlag = "   << c.su4NormaliseFlag   << std::endl
			  << "\t PressFlag = "          << c.PressFlag          << std::endl
			  << "\t xDerivativeOrder = "   << c.xDerivativeOrder   << std::endl;

	std::cout << "> Profile:" << std::endl;
	std::cout << "\t z_0   = "       << c.z_0/ph::km << " km"  << std::endl
			  << "\t ProfileFlag = " << c.ProfileFlag << std::endl;

	std::cout << "> Old fashioned noise" << std::endl;
	std::cout << "\t OldNoiseFlag = " << c.OldNoiseFlag << std::endl;
	std::cout << "\t N_OldNoise = "   << c.N_OldNoise << std::endl;
	std::cout << "\t lumSinLeft = [";
	for (auto it = c.lumSinLeft.begin(); it != c.lumSinLeft.end(); ++it)
		std::cout << *it << " ";
	std::cout << "]" << std::endl;

	std::cout << "\t lumCosLeft = [";
	for (auto it = c.lumCosLeft.begin(); it != c.lumCosLeft.end(); ++it)
		std::cout << *it << " ";
	std::cout << "]" << std::endl;

	std::cout << "\t lumSinRight = [";
	for (auto it = c.lumSinRight.begin(); it != c.lumSinRight.end(); ++it)
		std::cout << *it << " ";
	std::cout << "]" << std::endl;

	std::cout << "\t lumCosRight = [";
	for (auto it = c.lumCosRight.begin(); it != c.lumCosRight.end(); ++it)
		std::cout << *it << " ";
	std::cout << "]" << std::endl;
}

// Same but fot the derived ones
void DerivedConstants(Constants& c)
{
    std::cout << "> Derived constants, angles:"    << std::endl;
    std::cout << "\t cosOmega = " << c.cosOmega  << std::endl
              << "\t cosChi = "   << c.cosChi    << std::endl
              << "\t tanChi = "   << c.tanChi    << std::endl
              << "\t s2theta = "  << c.s2theta   << std::endl
              << "\t c2theta = "  << c.c2theta   << std::endl;

    std::cout << "> Scheme steps:" << std::endl;
    std::cout << "\t dz = " << c.dz / ph::km << " km" << std::endl
              << "\t dx = " << c.dx / ph::km << " km" << std::endl;

	std::cout << "> Matrices and norms:" << std::endl;
	std::cout << "\t MeanLeft = "  << std::endl
			  << c.MeanLeft        << std::endl;
	std::cout << "\t MeanRight = " << std::endl
			  << c.MeanRight       << std::endl;
	std::cout << "\t NormLeft = "  << c.NormMeanLeft << std::endl
			  << "\t NormRight = " << c.NormMeanRight << std::endl;
}

int main()
{
    Constants c;
    c.gather("./Parameters copy.json");
    BaseConstants(c);
    DerivedConstants(c);

    return 0;
}