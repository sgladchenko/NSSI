#pragma once

#include "Containers.h"
#include "Functions.h"
#include "Constants.h"
#include "Noise.h"

class Scheme
{
    public:
        // This constructor initializes the inner variables
        Scheme(String fParams, String fNoise, String dumpDir);

        // The main function that calculates each new line
        void Solve(bool whichNormsFlag);


    private:
        // The TwoLines object that will be calculating...
        TwoLines Rho;
        // ... taking known values  of the previous lines from here
        TwoLines RhoPrev;
        // su4::Norm's at the initial conditions
        TwoNorms InitNorms;

        // Objects containing data obtained from JSON's
        Constants c;
        Noise n;

        // The directory in which the resultimg data will be dumped
        String dumpDir;
};
