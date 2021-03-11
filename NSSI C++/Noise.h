#pragma once

#include "Types.h"
#include "su4.h"

// The class inilar to the class Constants, but it captures and dumps
// the parameters of the noise in the initial conditions
class Noise
{
    public:
        // Meta
        int  N_Noise;
        bool NoiseFlag;
        // Particular Fourier coefficients
        std::vector<su4::Vector> sinCoeffsLeft;
        std::vector<su4::Vector> cosCoeffsLeft;
        std::vector<su4::Vector> sinCoeffsRight;
        std::vector<su4::Vector> cosCoeffsRight;
        // The text in the JSON file
        String noiseText;

        // Added: gathering the data from JSON right here
        Noise(String filename);
        // Method that particularly gathers the data from JSON
        void gather(String filename);
        // And finally dump the initial JSON file into the final one
        // alongside the results of the calcuations
        void dump(String filename);
};