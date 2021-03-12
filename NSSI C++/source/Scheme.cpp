#include "Scheme.h"

// RHS in the equation that derfines the z-derivatives.
// In this case we use the central 2th or 4th approximation of the
// x-derivative of the density matrices
TwoLines RHSCentral(const Constants& c, const TwoLines& rho, Real z)
{
    TwoLines xDer = rho.centralDerivative(c.dx, c.xDerivativeOrder);
    TwoLines ham  = Hamiltonians(c, rho, z);
    return c.tanChi*(xDer.pm()) + com(ham, rho)*(-1.0i / c.cosOmega);
}

// Evaluates next line in the scheme RK4 (in z direction) + central approximation
// of the x derivative
TwoLines RK4AndCentral(const Constants& c, const TwoLines& rho, Real z)
{
    // RK4 components
    TwoLines K1 = RHSCentral(c, rho, z);
    TwoLines K2 = RHSCentral(c, rho + K1*(c.dz/2.0), z + c.dz/2.0);
    TwoLines K3 = RHSCentral(c, rho + K2*(c.dz/2.0), z + c.dz/2.0);
    TwoLines K4 = RHSCentral(c, rho + K3*c.dz, z + c.dz);
    // Final approximation at z+c.dz
    return rho + (K1 + 2.0*K2 + 2.0*K3 + K4)*(c.dz/6.0);
}

Scheme::Scheme(String fParams, String fNoise, String root, int px, int pz)
       : c(fParams), n(fNoise), periodN_x(px), periodN_z(pz)
{
    // Set the directories
    dir = setDir(root);
    // Dump the JSON's
    c.dump(dir + "Parameters.json");
    n.dump(dir + "Noise.json");
    // Dump XGrid
    dumpXGrid(c, periodN_x, dir);
    // First. let's generate the initial conditions
    RhoPrev = InitialConditions(c, n);
    // Let's evaluate the initial norms over here
    InitNorms = InitialNorms(c, RhoPrev);
    // And dump initial conditions as a first line
    dumpTwoLines(c, 0.0, periodN_x, RhoPrev, dir);
}

// Subsequently find all the lines, from z=c.dz to z=c.Z
void Scheme::Solve()
{
    for (int j=0; j<c.N_z; ++j)
    {
        Real z = j*c.dz;
        // Find an approximation of the next line, at z+c.dz
        RhoNext = RK4AndCentral(c, RhoPrev, z);
        // If needed, su4::Normalise them to the initial norms
        if (c.su4NormaliseFlag)
        {
            if (c.WhichNorms)
            { RhoNext.su4Normalise(TrueNorms(c, n, InitNorms, z+c.dz)); }
            else
            { RhoNext.su4Normalise(ApproximateNorms(c, InitNorms, z+c.dz)); }
            // Note that the position is actully z+c.dz, because RhoNext
            // is placed a the next line
        }
        // Dump the new lines (periodically)
        if ((j+1) % periodN_z == 0)
        {
            dumpTwoLines(c, z+c.dz, periodN_x, RhoNext, dir);
        }
    }
}