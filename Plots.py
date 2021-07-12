#!/usr/bin/env python3

import numpy as np
import os, sys, json
from Modules import FourPlots1D, OnePlot2D_PNG, OnePlot2D_EPS, Grids, Data, Rare, ObtainParameters, Symm, Averages

# Additional constants
labels = [r"$\nu_e$", r"$\nu_x$", r"$\bar{\nu}_e$", r"$\bar{\nu}_x$"]

def min2D(arr):
    return [min(each) for each in arr]

# Draws the basic plots that show the probabilities of the
# x,z-plane
def BasicPlots(directory, pX=1, pZ=1, fmt="png"):
    # Location of the plots
    plotDir = directory + "plots/"
    # Obtain grids
    XGrid, ZGrid, D_x, D_z = Grids(directory)
    # Obtain data
    print("Obtaining the data...")
    data = Data(directory, D_x, D_z, pX, pZ)
    # Make rare meshes
    XGridRare, ZGridRare = Rare(XGrid, ZGrid, pX, pZ)

    # Then let's draw the plots
    print("Drawing compressed plots located at: {}".format(plotDir))
    filelabels  = ["eNu", "xNu", "eANu", "xANu"]
    titlelabels = ["Electron neutrino probability",
                   "X neutrino probability",
                   "Electron antineutrino probability",
                   "X antineutrino probability"]

    filenames = [f"{each}L.{fmt}" for each in filelabels] + [f"{each}R.{fmt}" for each in filelabels]
    titles    = [f"{each} (Left)"for each in titlelabels] + [f"{each} (Right)"for each in titlelabels]

    filenamesSymm = [f"{each}Symm.{fmt}" for each in filelabels]
    titlesSymm = [f"{each} (Symmetric)" for each in titlelabels]
    
    # Left/right dependencies
    for k,fn in enumerate(filenames):
        if fmt == "png":
            OnePlot2D_PNG(XGridRare, ZGridRare, data[k], plotDir  + fn, titles[k], dpi=512)
        else:
            OnePlot2D_EPS(XGridRare, ZGridRare, data[k], plotDir  + fn)
    
    # Symmetric functions
    symm = Symm(data)
    for k,fn in enumerate(filenamesSymm):
        if fmt == "png":
            OnePlot2D_PNG(XGridRare, ZGridRare, symm[k], plotDir  + fn, titlesSymm[k], dpi=512)
        else:
            OnePlot2D_EPS(XGridRare, ZGridRare, symm[k], plotDir  + fn)

    # And then let's draw the average values
    avsL, avsR = Averages(data)
    if fmt == "png":
        FourPlots1D(ZGridRare, avsL, plotDir + f"avsL.{fmt}", labels, "Average probabilities (Left)")
        FourPlots1D(ZGridRare, avsR, plotDir + f"avsR.{fmt}", labels, "Average probabilities (Right)")
    else:
        FourPlots1D(ZGridRare, avsL, plotDir + f"avsL.{fmt}", labels)
        FourPlots1D(ZGridRare, avsR, plotDir + f"avsR.{fmt}", labels)

    # And one more nice functionaly: let's save the average
    # probabilities in a special JSON file
    averagesDir = directory + "averages/"
    if "averages" not in [f.name for f in os.scandir(directory)]:
        os.mkdir(averagesDir)
    
    with open(averagesDir + "avsL.json", "w") as fout:
        json.dump({"ZGridRare": ZGridRare, "avs": avsL}, fout, indent=4)
    with open(averagesDir + "avsR.json", "w") as fout:
        json.dump({"ZGridRare": ZGridRare, "avs": avsR}, fout, indent=4)

    return data, avsL, avsR, XGridRare, ZGridRare

if __name__ == "__main__":
    directory, pX, pZ = ObtainParameters()

    BasicPlots(directory, pX, pZ, fmt="eps")