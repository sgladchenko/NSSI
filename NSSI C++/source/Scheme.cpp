#include "Scheme.h"

// Evaluates next lines in the scheme RK4 (in z direction) + central approximation
// of the x derivative
TwoLines RK4AndCentral(const Constants& c, const TwoLines& prev)
{

}

Scheme::Scheme(String fParams, String fNoise, String ddumpDir)
       : c(fParams), n(fNoise), dumpDir(ddumpDir)
{
    // First. let's generate the initial conditions
    RhoPrev = InitialConditions(c, n);
    // Let's evaluate the initial norms over here
    InitNorms = InitialNorms(c, RhoPrev);
}

void Scheme::Solve(bool whichNormsFlag)
{
    
}