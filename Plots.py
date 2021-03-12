#!/usr/bin/env python3

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
colours = colour_e, colour_x

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
def OnePlot2D(xCoordinates, zCoordinates, Mesh, filePlot, title, cs=bone):
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
	fig.savefig(filePlot, fmt="png", bbox_inches='tight') #, dpi=512

# Standard 1D-plot
def OnePlot1D(zCoordinates, Mesh, filename, title):
	pass

if __name__ == "__main__":
    # Obtain the name of the directoy where the binaries are dumped
    dir = sys.argv[1]
    if dir[-1] != "/": dir += "/"

    # Then let's read the data; obtain the binaries
    LData  = np.fromfile(dir + "bin/left.bin",  dtype=np.float64)
    RData  = np.fromfile(dir + "bin/right.bin", dtype=np.float64)

    # Then let's read the grids
    XGrid = []
    ZGrid = []
    with open(dir + "XGrid.txt") as f:
        listed = f.read().split(" ")
        if listed[-1] == "":
            XGrid = [float(x) for x in listed[:-1]]
        else:
            XGrid = [float(x) for x in listed]

    with open(dir + "ZGrid.txt") as f:
        listed = f.read().split(" ")
        if listed[-1] == "":
            ZGrid = [float(x) for x in listed[:-1]]
        else:
            ZGrid = [float(x) for x in listed]

    # Then let's deserialize the data obtained in binaries
    D_x = len(XGrid) # number of the displayed points
    D_z = len(ZGrid) # number of the displayed points

    eNuL  = [[LData[4*D_x*j + 4*i + 0] for i in range(D_x)] for j in range(D_z)]
    xNuL  = [[LData[4*D_x*j + 4*i + 1] for i in range(D_x)] for j in range(D_z)]
    eANuL = [[LData[4*D_x*j + 4*i + 2] for i in range(D_x)] for j in range(D_z)]
    xANuL = [[LData[4*D_x*j + 4*i + 3] for i in range(D_x)] for j in range(D_z)]    

    eNuR  = [[RData[4*D_x*j + 4*i + 0] for i in range(D_x)] for j in range(D_z)]
    xNuR  = [[RData[4*D_x*j + 4*i + 1] for i in range(D_x)] for j in range(D_z)]
    eANuR = [[RData[4*D_x*j + 4*i + 2] for i in range(D_x)] for j in range(D_z)]
    xANuR = [[RData[4*D_x*j + 4*i + 3] for i in range(D_x)] for j in range(D_z)]

    # Then let's draw the plots
    plotDir = dir + "plots/"

    OnePlot2D(XGrid, ZGrid, eNuL, plotDir  + "eNuL.png",  "Electron neutrino probability (Left)")
    OnePlot2D(XGrid, ZGrid, xNuL, plotDir  + "xNuL.png",  "X neutrino probability (Left)")
    OnePlot2D(XGrid, ZGrid, eANuL, plotDir + "eANuL.png", "Electron antineutrino probability (Left)")
    OnePlot2D(XGrid, ZGrid, xANuL, plotDir + "xANuL.png", "X antineutrino probability (Left)")

    OnePlot2D(XGrid, ZGrid, eNuR, plotDir  + "eNuR.png",  "Electron neutrino probability (Right)")
    OnePlot2D(XGrid, ZGrid, xNuR, plotDir  + "xNuR.png",  "X neutrino probability (Right)")
    OnePlot2D(XGrid, ZGrid, eANuR, plotDir + "eANuR.png", "Electron antineutrino probability (Right)")
    OnePlot2D(XGrid, ZGrid, xANuR, plotDir + "xANuR.png", "X antineutrino probability (Right)")