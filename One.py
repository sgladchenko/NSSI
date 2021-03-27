#!/usr/bin/env python3

import numpy as np
import os
import sys
from Modules import OnePlot2D_PNG, Grids, Data, Rare, ObtainParameters

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
    print("Drawing high res plot of eNuL at: {}".format(plotDir))
    OnePlot2D_PNG(XGridRare, ZGridRare, data[0], plotDir  + "eNuL_Precise.png",  "Electron neutrino probability (Left)", dpi=512)