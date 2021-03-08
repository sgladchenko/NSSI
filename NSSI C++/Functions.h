#pragma once

// The module contains the functions needed throughout the main 
// calculational process

#include "Noise.h"
#include "Constants.h"
#include "Containers.h"

// The function that initializes the TwoLines of the initial density matrices
TwoLines InitialConditions(const Constants& c, const Noise& n);

// Functions with particular formulae of the matrix flavour Hamiltonians
Matrix VacuumHamiltonian(const Constants& c);
Matrix MSWHamiltonian(const Constants& c, Real z);
Matrix AMMHamiltonian(const Constants& c, Real z);

// and the collective ones
Matrix VAHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z);
Matrix NSSIHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z);

// The functions that generate the total Hamilatonians in a TwoLines instance
TwoLines Hamiltonians(const Constants& c, const TwoLines& rho, Real z);

// Some smaller functions needed
Real originLeft(Real x, Real z);
Real originRight(Real x, Real z);

// Function that evaluates the norms (and saves them in the vector put in the arguements
// by references)
void Norms(const TwoLines& rho, std::vector<Real>& leftNorms, std::vector<Real>& rightNorms);