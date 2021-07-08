#pragma once

#include "Types.h"
#include "Constants.h"
#include "Noise.h"
#include "Containers.h"

namespace hamiltonians
{
    // Vacuum Hamiltonian of the oscillations
    Matrix VacuumHamiltonian(const Constants& c);
    // MSW Term of the oscillations
    Matrix MSWHamiltonian(const Constants& c, Real z);
    // The term that describes possible interaction with magnetic field
    // through the anomalous magnetic moment of neutrino
    Matrix AMMHamiltonian(const Constants& c, Real z);

    // Standard Nodel collective term
    Matrix VAHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z, Real factor);
    // The collective term that describes coupling through the scalar and 
    // pesudoscalar fields
    Matrix NSSIHamiltonian(const Constants& c, const Matrix& rhoOpp, Real z, Real factor);
}