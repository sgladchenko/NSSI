#pragma once

#include "Containers.h"
#include "Functions.h"
#include "Constants.h"
#include "Noise.h"
#include "Files.h"

class Scheme
{
    public:
        // This constructor initializes the inner variables
        Scheme(String fParams, String fNoise, String root, int px, int pz);

        // The main function that calculates each new line
        void Solve();

    private:
        // The TwoLines object that will be calculated...
        TwoLines RhoNext;
        // ... taking known values of the previous two lines from here
        TwoLines RhoPrev;
        // su4::Norm's at the initial conditions
        TwoNorms InitNorms;

        // Objects containing data obtained from JSON's
        Constants c; Noise n;

        // The directory in which the resultimg data will be dumped
        String dir;

        // Periods of the displayed grids
        int periodN_x;
        int periodN_z;
};
