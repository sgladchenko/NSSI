#!/usr/bin/env python3

from Lambdas import flavourSigma1,flavourSigma3
import json
import numpy as np

# Definition of the CGS units and derived ones

cm   = 1.
g    = 1.
sec  = 1.
erg  = 1.
G    = 1.

eV   = 1.602e-12 * erg
MeV  = 1.0e6 * eV
km   = 1.0e5 * cm

muB  = 9.274009994e-21 * erg / G

G_F  = 1.166e-11  / MeV**2
hbar = 1.0546e-27 * erg * sec
c    = 2.9979e10  * cm  / sec
deg  =  np.pi / 180.

hbarc = hbar*c

# A dictionary that connects string-labels of units in a JSON-file 
# with the variables above; only actual units are listed

dictUnits = {"deg": deg, "eV": eV, "eV^2": eV*eV, "MeV": MeV, "cm^-3": 1.0/cm**3, "km": km, "G": G, "muB": muB}
dictParameters = {}
dictNoise = {}

# A small parser, that converts strings to the floats and ints
# multiplying on their dimensional constants

def getValue(dictData, section, valKey, valType=float, hasUnit=True):
    if hasUnit:
        return valType(dictData[section][valKey][0]) * dictUnits[dictData[section][valKey][1]]
    else:
        return valType(dictData[section][valKey])

def getArray(dictData, section, valKey, valType=float, hasUnit=True):
    if hasUnit:
        unit = dictUnits[dictData[section][valKey][1]]
        return [valType(x)*unit for x in dictData[section][valKey][0]]
    else:
        return [valType(x) for x in dictData[section][valKey]]

def getDoubleArray(dictData, section, valKey, valType=float, hasUnit=True):
    if hasUnit:
        unit = dictUnits[dictData[section][valKey][1]]
        return [np.array([valType(y)*unit for y in x]) for x in dictData[section][valKey][0]]
    else:
        return [np.array([valType(y) for y in x]) for x in dictData[section][valKey]]


""" Parameters needed in the linear stability analysis """

class Constants:
    def __init__(self, filename="./Parameters.json"):
        with open(filename, "r") as fileJSON:
            dictParameters = json.load(fileJSON)

            # Basic physical parameters

            self.eta     = getValue(dictParameters, "Basic", "eta", hasUnit=False)

            self.theta   = getValue(dictParameters, "Basic", "theta")
            self.dm2_0   = getValue(dictParameters, "Basic", "Delta m^2_0")
            self.dm2     = self.dm2_0 * self.eta

            self.chi     = getValue(dictParameters, "Basic", "chi")
            self.E_0     = getValue(dictParameters, "Basic", "E_0")

            # Parameters of MSW potentials and NSSI/V-A couplings

            self.gMinus  = getValue(dictParameters, "V-A & NSSI & MSW & AMM", "g_{-}", hasUnit=False)
            self.gPlus   = getValue(dictParameters, "V-A & NSSI & MSW & AMM", "g_{+}", hasUnit=False)
            self.n_Nu    = getValue(dictParameters, "V-A & NSSI & MSW & AMM", "n_Nu")
            self.n_e     = getValue(dictParameters, "V-A & NSSI & MSW & AMM", "n_e")
            self.n_n     = getValue(dictParameters, "V-A & NSSI & MSW & AMM", "n_n")

            # Scheme

            self.X = getValue(dictParameters, "Scheme", "X")

            # Initial parameters

            self.InitialProbsLeft  = getArray(dictParameters, "Initial Mixes", "InitialProbsLeft", hasUnit=False)
            self.InitialProbsRight = getArray(dictParameters, "Initial Mixes", "InitialProbsRight", hasUnit=False)

            # Constants in the collective terms of Hamiltonians
            self.V_Nu  = G_F * np.sqrt(2) * self.n_Nu * (hbarc)**2
            self.V_e   = G_F * np.sqrt(2) * self.n_e * (hbarc)**2
            self.V_n   = G_F * np.sqrt(2) * self.n_n * (hbarc)**2 

            # Construct the initial (mean) density matrix
            self.InitialRho = np.matrix(np.diag(self.InitialProbsLeft))

            # Functions of the angle chi
            self.cosOmega = np.cos(2.0 * self.chi)
            self.cosChi   = np.cos(self.chi)
            self.tanChi   = np.tan(self.chi)

            self.vacOmega = self.dm2_0/(2.0*self.E_0*hbarc)

            # Flag of hierarchy
            if self.eta == 1.0:
                self.Hierarchy = "NH"
            else:
                self.Hierarchy = "IH"

        with open(filename, "r") as fileJSON:
            self.text = fileJSON.read()

    # Write the JSON file with these constants
    def dump(self, filename):
        with open(filename, "w") as f:
            f.write(self.text)

# # # Test section

if __name__ == "__main__":
    c = Constants()
    omega = c.dm2_0 / (2.0 * c.E_0 * hbarc)
    print(c.V_Nu * (1.0 - c.cosOmega) / omega)
    print(omega)