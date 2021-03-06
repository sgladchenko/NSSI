#pragma once

// The module contains the functions needed throughout the main 
// calculational process

#include "Noise.h"
#include "Constants.h"
#include "Containers.h"

// Origin of the left beam going thorugh the point (x,z)
Real originLeft(const Constants& c, Real x, Real z);
// Origin of the right beam going thorugh the point (x,z)
Real originRight(const Constants& c, Real x, Real z);

// Evaluates a sum of Fourier harmonics of su4::Vector's
Matrix FourierVector(Real x, Real X, const std::vector<su4::Vector>& sinVecs,
                                     const std::vector<su4::Vector>& cosVecs);

// The function that initializes the TwoLines of the initial density matrices
TwoLines InitialConditions(const Constants& c, const Noise& n);

// The functions that generate the total Hamilatonians in a TwoLines instance
TwoLines Hamiltonians(const Constants& c, const TwoLines& rho, Real z);

// Function that evaluates the norms of the initial conditions
TwoNorms InitialNorms(const Constants& c, const TwoLines& rho);

// Function that approximately calculates the norms at a given position z
// (approximaely means that it uses linear interpolation of the initial norms
// in order to figure the norms between the nodes of the grid)
TwoNorms ApproximateNorms(const Constants& c, const TwoNorms& initialNorms, Real z);

// Same as above, but it actually calculates them directly
TwoNorms TrueNorms(const Constants& c, const Noise& n, const TwoNorms& initNorms, Real z);

// Regularization that sets the eigenvalues to ones that are set in the initial density matrices
void RegEigenvalues(const Constants& c, const Noise& n, TwoLines& tl, Real z);