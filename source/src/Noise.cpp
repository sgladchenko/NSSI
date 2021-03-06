#include "Noise.h"

#include <fstream>
#include <iostream>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

void getVectors(json& j, String section, String key, std::vector<su4::Vector>& stdvec)
{
    for (auto it  = j[section][key].begin(); it != j[section][key].end(); ++it)
    {
        Array15 tmp;
        for (int i = 0; i < 15; ++i)
        {
            tmp[i] = (Real)((*it)[i].get<double>());
        }
        stdvec.push_back(su4::Vector(tmp));
    }
}

// The default constructor
Noise::Noise(String filename) : NoiseFlag {false}
{
    // Similarly to the Constants class
    std::ifstream jsonParameters(filename);
    std::stringstream buf;

    buf << jsonParameters.rdbuf();
    noiseText = buf.str();
    jsonParameters.close();

    // Secondly, let's parse it via nlohmann
    json j = json::parse(noiseText);

    // Finally, let's save the data
    N_Noise = j["Meta"]["N_Noise"].get<int>();

    // And the Fourier components themselves
    getVectors(j, "Harmonics", "sinCoeffsLeft",  sinCoeffsLeft);
    getVectors(j, "Harmonics", "cosCoeffsLeft",  cosCoeffsLeft);
    getVectors(j, "Harmonics", "sinCoeffsRight", sinCoeffsRight);
    getVectors(j, "Harmonics", "cosCoeffsRight", cosCoeffsRight);
}

// Similarly to Constants
void Noise::dump(String filename)
{
	std::ofstream out(filename);
	out << noiseText;
	out.close();
}
