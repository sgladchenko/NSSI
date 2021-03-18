#!/usr/bin/env python3

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker, rcParams, cm

import matplotlib as mpl
import numpy as np
import os
import sys

# Colours and other constants
colour_e = "#010003"
colour_x = "#2D3F88"
colours  = colour_e, colour_x

viridis  = cm.get_cmap('viridis')
twilight = cm.get_cmap('twilight')
bone     = cm.get_cmap('bone')

MainFontSize = 12

# Some customisation of matplotlib
ticker.rcParams['xtick.direction'] = 'in'
ticker.rcParams['ytick.direction'] = 'in'
ticker.rcParams['xtick.labelsize'] = MainFontSize
ticker.rcParams['ytick.labelsize'] = MainFontSize

mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif']  = ['Times New Roman'] # Computer Modern Roman
mpl.rcParams['legend.handlelength'] = 2

# Standard 2D-plot
def OnePlot2D(xCoordinates, zCoordinates, Mesh, filePlot, title, cs=bone, dpi=None):
    MeshArray = np.array([[np.real(x) for x in line] for line in Mesh])

    xArray = np.array(xCoordinates)
    zArray = np.array(zCoordinates)
    xArray, zArray = np.meshgrid(xArray, zArray)

    fig = Figure(figsize=(8, 6))
    FigureCanvas(fig)

    axs  = fig.add_subplot(111)
    plot = axs.pcolor(xArray, zArray, MeshArray, cmap=cs, vmin=0.0, vmax=1.0)

    axs.set_title(title, fontsize=MainFontSize)
    axs.set_xlabel(r"$x, \mathrm{km}$")
    axs.set_ylabel(r"$z, \mathrm{km}$")

    fig.colorbar(plot)

    if dpi:
        fig.savefig(filePlot, fmt="png", bbox_inches='tight', dpi=dpi)
    else:
        fig.savefig(filePlot, fmt="png", bbox_inches='tight')

# Obtain the grids
def Grids(dir):
    XGrid = []; D_x = 0 
    ZGrid = []; D_Z = 0

    # Grid of the x-axis that's been dumped right after
    # the calculations by c++ code
    with open(dir + "XGrid.txt") as f:
        listed = f.read().split(" ")
        if listed[-1] == "":
            XGrid = [float(x) for x in listed[:-1]]
        else:
            XGrid = [float(x) for x in listed]

        D_x = len(XGrid)

    # Grid of the z-axis that's been dumped right after
    # the calculations by c++ code
    with open(dir + "ZGrid.txt") as f:
        listed = f.read().split(" ")
        if listed[-1] == "":
            ZGrid = [float(x) for x in listed[:-1]]
        else:
            ZGrid = [float(x) for x in listed]
        
        D_z = len(ZGrid)

    return XGrid, ZGrid, D_x, D_z

def Compressed():
    pass

if __name__ == "__main__":
    dir = ""    # Obtain the name of the directoy where the binaries are dumped
    pX = 1      # These are optional periods, needed if we didn't se such periods
    pZ = 1      # in th beginning of calculations
    onePrecise = None
    plots = True

    for each in sys.argv[1:]:
        if '=' in each:
            key, val = each.split("=")[0].replace("--",""), each.split("=")[1]

            if key == "dir":
                dir = val
                if dir[-1] != "/": dir += "/" # Add slash in th end if needed
                print("Location used: " + dir)
                
            elif key == "pX":
                pX = int(val)
                print("Additional period in x-axis: {}".format(pX))

            elif key == "pZ":
                pZ = int(val)
                print("Additional period in z-axis: {}".format(pZ))

            elif key == "onePrecise":
                onePrecise = int(val)
                print("Additional plot of eNuL with periods: {}".format(onePrecise))

            else:
                print("Undefined key")
            
        else:
            if each == "--noPlots":
                plots = False
    
    # Location of the plots
    plotDir = dir + "plots/"

    # Then let's read the data; obtain the binaries
    LData  = np.fromfile(dir + "bin/left.bin",  dtype=np.float64)
    RData  = np.fromfile(dir + "bin/right.bin", dtype=np.float64)


    XGridR = []
    ZGridR = []

    if plots:
        XGridRare = [x for ix, x in enumerate(XGrid) if ix % pX == 0]
        ZGridRare = [z for iz, z in enumerate(ZGrid) if iz % pZ == 0]

        # Then let's deserialize the data obtained in binaries
        eNuL  = [[LData[4*D_x*j + 4*i + 0] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
        xNuL  = [[LData[4*D_x*j + 4*i + 1] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
        eANuL = [[LData[4*D_x*j + 4*i + 2] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
        xANuL = [[LData[4*D_x*j + 4*i + 3] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]    

        eNuR  = [[RData[4*D_x*j + 4*i + 0] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
        xNuR  = [[RData[4*D_x*j + 4*i + 1] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
        eANuR = [[RData[4*D_x*j + 4*i + 2] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
        xANuR = [[RData[4*D_x*j + 4*i + 3] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]

        # Then let's draw the plots
        print("Main plots...")
        OnePlot2D(XGridRare, ZGridRare, eNuL, plotDir  + "eNuL.png",  "Electron neutrino probability (Left)")
        OnePlot2D(XGridRare, ZGridRare, xNuL, plotDir  + "xNuL.png",  "X neutrino probability (Left)")
        OnePlot2D(XGridRare, ZGridRare, eANuL, plotDir + "eANuL.png", "Electron antineutrino probability (Left)")
        OnePlot2D(XGridRare, ZGridRare, xANuL, plotDir + "xANuL.png", "X antineutrino probability (Left)")

        OnePlot2D(XGridRare, ZGridRare, eNuR, plotDir  + "eNuR.png",  "Electron neutrino probability (Right)")
        OnePlot2D(XGridRare, ZGridRare, xNuR, plotDir  + "xNuR.png",  "X neutrino probability (Right)")
        OnePlot2D(XGridRare, ZGridRare, eANuR, plotDir + "eANuR.png", "Electron antineutrino probability (Right)")
        OnePlot2D(XGridRare, ZGridRare, xANuR, plotDir + "xANuR.png", "X antineutrino probability (Right)")

    if onePrecise:
        p = onePrecise
        XGridPrec = [x for ix, x in enumerate(XGrid) if ix % p == 0]
        ZGridPrec = [x for ix, x in enumerate(ZGrid) if ix % p == 0]

        # Then let's deserialize the data obtained in binaries
        eNuLPrec  = [[LData[4*D_x*j + 4*i + 0] for i in range(D_x) if i % p == 0] for j in range(D_z) if j % p == 0]

        print("Drawing one precise plot...")
        OnePlot2D(XGridPrec,
                  ZGridPrec,
                  eNuLPrec,
                  plotDir  + "eNuL_p.png",
                  "Electron neutrino probability (Left)", dpi=512)    