#!/usr/bin/env python3

import numpy as np
import os
import sys
from Modules import OnePlot2D_PNG, Grids, Data, Rare, ObtainParameters

def min2D(arr):
    return [min(each) for each in arr]

if __name__ == "__main__":
    directory, pX, pZ = ObtainParameters()        

    # Location of the plots
    plotDir = directory + "plots/"

    # Obtain grids
    XGrid, ZGrid, D_z, D_x = Grids(directory)

    # Obtain data
    print("Obtaining the data...")
    data = Data(directory, D_x, D_z, pX, pZ)

    # Make rare meshes
    XGridRare, ZGridRare = Rare(XGrid, ZGrid, pX, pZ)

    # Then let's draw the plots
    print("Drawing compressed plots located at: {}".format(plotDir))
    OnePlot2D_PNG(XGridRare, ZGridRare, data[0], plotDir  + "eNuL.png", "Electron neutrino probability (Left)")
    OnePlot2D_PNG(XGridRare, ZGridRare, data[1], plotDir  + "xNuL.png", "X neutrino probability (Left)")
    OnePlot2D_PNG(XGridRare, ZGridRare, data[2], plotDir + "eANuL.png", "Electron antineutrino probability (Left)")
    OnePlot2D_PNG(XGridRare, ZGridRare, data[3], plotDir + "xANuL.png", "X antineutrino probability (Left)")

    OnePlot2D_PNG(XGridRare, ZGridRare, data[4], plotDir  + "eNuR.png", "Electron neutrino probability (Right)")
    OnePlot2D_PNG(XGridRare, ZGridRare, data[5], plotDir  + "xNuR.png", "X neutrino probability (Right)")
    OnePlot2D_PNG(XGridRare, ZGridRare, data[6], plotDir + "eANuR.png", "Electron antineutrino probability (Right)")
    OnePlot2D_PNG(XGridRare, ZGridRare, data[7], plotDir + "xANuR.png", "X antineutrino probability (Right)")

    #for i, each in enumerate(min2D(data[1])): print(i, each)