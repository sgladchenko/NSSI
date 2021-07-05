#!/usr/bin/env python3

import os, json, itertools
from copy import deepcopy

# Very important constant; it means what's the value of
# mu/omega corresponding to the n_Nu=1e29 cm^-3
muOf1e29 = 13.976129602265964

# A simple converter
mu_to_n = lambda mu: 1.0e29 * (mu / muOf1e29)

# It uses a JSON file and just changes in it only
# entries of mu and gPlus, then dumping into the target file
def makeJSONfromDict(eta, mu, N, gPlus, initialdict, targetJSON):
    pdict = deepcopy(initialdict)
    # Then let's just substitute the information in the entries
    pdict["Basic"]["eta"] = eta
    pdict["V-A & NSSI & MSW & AMM"]["g_{+}"] = gPlus
    pdict["V-A & NSSI & MSW & AMM"]["n_Nu"][0] = mu_to_n(mu)
    # Dimensions of the grid
    pdict["Scheme"]["N_z"] = N*2
    pdict["Scheme"]["N_x"] = N
    # Coordinates on the grid
    pdict["Scheme"]["Z"][0] = 20.0
    pdict["Scheme"]["X"][0] = 10.0
    # And finally let's dump it
    with open(targetJSON, "w") as f:
        json.dump(pdict, f, indent=4)

# This function assembles the JSON-files for our purposes
# and it wll make a two JSON files (for both values of hierarchy)
# for each pair of mu and gPlus from the lists below mentioned in the
# arguements of the function; note that the specified value of N_x/N_z
# will be gathered from the list of Ns below, and the particular value of
# N will be chosen according to the value of mu, so len(mus) must be equal to
# len(Ns).
def makeTask(muNPairs, gPluses, initialJSON, foldername):
    # First of all, let's make a tree of folders
    try:
        os.mkdir(foldername)
    except FileExistsError:
        print(f"{foldername} already exists; everything will have been dumped inside.")

    try:
        os.mkdir(f"{foldername}/NH")
    except FileExistsError:
        print(f"{foldername}/NH already exists; the files will be dumped here.")

    try:
        os.mkdir(f"{foldername}/IH")
    except FileExistsError:
        print(f"{foldername}/IH already exists; the files will be dumped here.")

    # Then let's setup the initial JSON file in which the entries of n_Nu, gPlus and eta will be changed
    with open(initialJSON) as f:
        initialdict = json.load(f)

    # And then let's make a structure of JSON files
    for muNPair, gPlus in itertools.product(muNPairs, gPluses):
        mu,N = muNPair
        makeJSONfromDict( 1.0, mu, N, gPlus, initialdict, f"{foldername}/NH/Parameters_mu={mu:.1f}_gPlus={gPlus:.2f}_NH.json")
        makeJSONfromDict(-1.0, mu, N, gPlus, initialdict, f"{foldername}/IH/Parameters_mu={mu:.1f}_gPlus={gPlus:.2f}_IH.json")


if __name__ == "__main__":
    gPluses = [0.00,  0.25,  0.50,  0.75,  1.00]
    muNPairs = [(10.0, 10000), (20.0, 14000), (30.0, 18000), (40.0, 22000), (50.0, 26000)]

    makeTask(muNPairs, gPluses, "./Parameters.json", "./Scans/Scan2D_20km_Parameters")