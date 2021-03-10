#pragma once

// The module contains the functions needed throughout the main 
// calculational process

#include "Noise.h"
#include "Constants.h"
#include "Containers.h"

// The function that initializes the TwoLines of the initial density matrices
TwoLines InitialConditions(const Constants& c, const Noise& n);

// The functions that generate the total Hamilatonians in a TwoLines instance
TwoLines Hamiltonians(const Constants& c, const TwoLines& rho, Real z);

// Function that evaluates the norms of the initial conditions
void InitialNorms(const TwoLines& rho,
                  std::vector<Real>& leftNorms,
                  std::vector<Real>& rightNorms);

// Function that approximately calculates the norms at a given position z
// (approximaely means that it uses linear interpolation of the initial norms
// in order to figure the norms between the nodes of the grid)
void ApproximateNorms(const std::vector<Real>& initLeftNorms,
                      const std::vector<Real>& initRightNorms,
                      std::vector<Real>& leftNorms,
                      std::vector<Real>& rightNorms,
                      Real z);

// Same as above, but it actually calculates them directly
void TrueNorms(const Constants& c,
               const Noise& n,
               std::vector<Real>& leftNorms,
               std::vector<Real>& rightNorms,
               Real z);