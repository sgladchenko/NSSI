from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker, rcParams, cm

import matplotlib as mpl
import numpy as np

import sys

# Colours and other constants
colour_e = "#324851"
colour_x = "#86AC41"
colour_ae = "#34675C"
colour_ax = "#7DA3A1"
colours  = colour_e, colour_x, colour_ae, colour_ax

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

# The functions that makes the rare grids
# (so the points are actually displayed with some periods)
def Rare(XGrid, ZGrid, pX, pZ):
    XGridRare = [x for ix, x in enumerate(XGrid) if ix % pX == 0]
    ZGridRare = [z for iz, z in enumerate(ZGrid) if iz % pZ == 0]

    return XGridRare, ZGridRare

# Obtain dir and periods from the console
def ObtainParameters():
    directory = ""    # Obtain the name of the directoy where the binaries are dumped
    pX = 1            # These are optional periods, needed if we didn't set such periods
    pZ = 1            # in the beginning of calculations

    for each in sys.argv[1:]:
        if '=' in each:
            key, val = each.split("=")[0].replace("--",""), each.split("=")[1]

            if key == "dir":
                directory = val
                if directory[-1] != "/": directory += "/" # Add slash in th end if needed
                print("Location used: " + directory)

            elif key == "periodN_x":
                pX = int(val)
                print("Additional period in x-axis: {}".format(pX))

            elif key == "periodN_z":
                pZ = int(val)
                print("Additional period in z-axis: {}".format(pZ))

            else:
                print("Undefined key")

    return directory, pX, pZ

# Standard 2D-plot
def OnePlot2D_PNG(xCoordinates, zCoordinates, Mesh, filePlot, title, cs=viridis, dpi=None):
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
    XGrid = []
    ZGrid = []

    # Grid of the x-axis that's been dumped right after
    # the calculations by c++ code
    with open(dir + "XGrid.txt") as f:
        listed = f.read().split(" ")
        if listed[-1] == "":
            XGrid = [float(x) for x in listed[:-1]]
        else:
            XGrid = [float(x) for x in listed]

    # Grid of the z-axis that's been dumped right after
    # the calculations by c++ code
    with open(dir + "ZGrid.txt") as f:
        listed = f.read().split(" ")
        if listed[-1] == "":
            ZGrid = [float(x) for x in listed[:-1]]
        else:
            ZGrid = [float(x) for x in listed]

    return XGrid, ZGrid, len(XGrid), len(ZGrid)

# Draw four one dimensional plots
def FourPlots1D_PNG(zCoordinates, values, filePlot, title, labels, maxV=True):
    fig1 = Figure(figsize=(7, 5))
    FigureCanvas(fig1)

    axs1 = fig1.add_subplot(1, 1, 1)
    axs1.set_xlabel(r"$x, \ \mathrm{km}$", fontsize=MainFontSize)
    axs1.set_ylabel(r"$\nu$ Probabilities", fontsize=MainFontSize)
    axs1.grid(False)
    axs1.set_xlim([zCoordinates[0], zCoordinates[-1]])
    if maxV:
        axs1.set_ylim([.0, 1.])
    else:
        axs1.set_ylim([.0, max([max(each) for each in values])])
    axs1.set_title(title, fontsize=MainFontSize)

    curves = [None for i in range(4)]
    for i, each in enumerate(values):
        curves[i] = axs1.plot(
            zCoordinates,
            each,
            '-',
            markersize=1.0,
            color=colours[i]
        )[0]

    axs1.legend(curves, labels, fontsize=MainFontSize, loc="upper right")
    fig1.savefig(filePlot, fmt="png", bbox_inches='tight')

# Obtain the full data
def Data(dir, D_x, D_z, pX, pZ):
    # Then let's read the data; obtain the binaries
    LData  = np.fromfile(dir + "bin/left.bin",  dtype=np.float64)
    RData  = np.fromfile(dir + "bin/right.bin", dtype=np.float64)
    
    # Deserialize
    eNuL  = [[LData[4*D_x*j + 4*i + 0] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
    xNuL  = [[LData[4*D_x*j + 4*i + 1] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
    eANuL = [[LData[4*D_x*j + 4*i + 2] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
    xANuL = [[LData[4*D_x*j + 4*i + 3] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]    

    eNuR  = [[RData[4*D_x*j + 4*i + 0] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
    xNuR  = [[RData[4*D_x*j + 4*i + 1] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
    eANuR = [[RData[4*D_x*j + 4*i + 2] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]
    xANuR = [[RData[4*D_x*j + 4*i + 3] for i in range(D_x) if i % pX == 0] for j in range(D_z) if j % pZ == 0]

    return eNuL, xNuL, eANuL, xANuL, eNuR, xNuR, eANuR, xANuR