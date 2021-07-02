from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker, rcParams, cm

import matplotlib as mpl
import numpy as np
from scipy.interpolate import griddata

import sys

# Default values of the dimensions of the plots
defaultDims = (6, 4.5)

# Colours and other constants
colour_e = "#FF0000"
colour_x = "#45BA45"
colour_ae = "#1034A6"
colour_ax = "#7F00FF"
colour_NB = "#050505"
colours  = colour_e, colour_x, colour_ae, colour_ax, colour_NB
viridis  = cm.get_cmap('viridis')
MainFontSize = 12

# Some customisation of matplotlib
ticker.rcParams['xtick.direction'] = 'in'
ticker.rcParams['ytick.direction'] = 'in'
ticker.rcParams['xtick.labelsize'] = MainFontSize
ticker.rcParams['ytick.labelsize'] = MainFontSize
mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif']  = ['Times New Roman'] # Computer Modern Roman
mpl.rcParams['legend.handlelength'] = 2

rcParams['mathtext.fontset'] = 'cm'#'dejavuserif'

# The functions that makes the rare grids
# (so the points are actually displayed with some periods)
def Rare(XGrid, ZGrid, pX, pZ):
    XGridRare = [x for ix, x in enumerate(XGrid) if ix % pX == 0]
    ZGridRare = [z for iz, z in enumerate(ZGrid) if iz % pZ == 0]
    return XGridRare, ZGridRare

# Make a symmetric function
def Symm(data):
    symm = []
    nz = len(data[0]); nx = len(data[0][0])
    for f in range(4):
        symm.append([[0.5*(data[f][iz][ix] + data[f+4][iz][ix]) for ix in range(nx)] for iz in range(nz)])

    return symm

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

# Aveage probabilities channel by channel
def Averages(data):
    # Average values of the left beam
    avProbsL = [[sum(line)/len(line) for line in mesh] for mesh in data[:4]]
    avProbsR = [[sum(line)/len(line) for line in mesh] for mesh in data[4:]]

    return avProbsL, avProbsR

# Standard 2D-plot
def OnePlot2D_PNG(xCoordinates, zCoordinates, Mesh, filePlot, title=None, cs=viridis, dpi=None, dims=defaultDims):
    MeshArray = np.array([[np.real(x) for x in line] for line in Mesh])

    xArray = np.array(xCoordinates)
    zArray = np.array(zCoordinates)
    xArray, zArray = np.meshgrid(xArray, zArray)

    fig = Figure(figsize=dims)
    FigureCanvas(fig)

    axs  = fig.add_subplot(111)
    plot = axs.pcolor(xArray, zArray, MeshArray, cmap=cs, vmin=0.0, vmax=1.0)

    if title: axs.set_title(title, fontsize=MainFontSize)
    axs.set_xlabel(r"$x$, km", fontsize=MainFontSize)
    axs.set_ylabel(r"$z$, km", fontsize=MainFontSize)
    axs.set_aspect('equal')
    
    fig.colorbar(plot)
    if dpi: fig.savefig(filePlot, fmt="png", bbox_inches='tight', dpi=dpi)
    else:   fig.savefig(filePlot, fmt="png", bbox_inches='tight')

# Same standard 2D-plot, but in EPS vector
def OnePlot2D_EPS(xCoordinates, zCoordinates, Mesh, filePlot, dpi=200, dims=defaultDims):
    MeshArray = np.array([[np.real(x) for x in line] for line in Mesh])

    xArray = np.array(xCoordinates)
    zArray = np.array(zCoordinates)
    xArray, zArray = np.meshgrid(xArray, zArray)

    fig = Figure(figsize=dims)
    FigureCanvas(fig)

    axs  = fig.add_subplot(111)
    plot = axs.pcolor(xArray, zArray, MeshArray, cmap=viridis, vmin=0.0, vmax=1.0, rasterized=True)

    axs.set_xlabel(r"$x$, km", fontsize=MainFontSize)
    axs.set_ylabel(r"$z$, km", fontsize=MainFontSize)
    axs.set_aspect('equal')

    fig.colorbar(plot)
    fig.savefig(filePlot, fmt="eps", bbox_inches='tight', dpi=dpi)

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
def FourPlots1D(zCoordinates, values, filePlot, labels, title=None, fmt="png", maxV=True, dims=defaultDims):
    fig1 = Figure(figsize=dims)
    FigureCanvas(fig1)

    axs1 = fig1.add_subplot(1, 1, 1)
    axs1.set_xlabel(r"$z$, km", fontsize=MainFontSize)
    axs1.set_ylabel(r"$\bar{P}_{f}(z)$", fontsize=MainFontSize)
    axs1.grid(False)
    axs1.set_xlim([zCoordinates[0], zCoordinates[-1]])

    if maxV: axs1.set_ylim([0.0, 1.0])
    else: axs1.set_ylim([0.0, max([max(each) for each in values])])

    if title: axs1.set_title(title, fontdict={'fontsize':MainFontSize,})

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
    fig1.savefig(filePlot, fmt=fmt, bbox_inches='tight', dpi=192)

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

# Functions that makes a plot
def PlotAverageFlavours(xs, xlims, xlabel, flavours, ylims, labels, fileplot, fmt="eps", dims=defaultDims):
    fig1 = Figure(figsize=dims)
    FigureCanvas(fig1)

    axs1 = fig1.add_subplot(1, 1, 1)
    axs1.set_xlabel(xlabel, fontsize=MainFontSize)
    #axs1.set_ylabel(r"$\langle \bar{P}_f \rangle$", fontsize=MainFontSize)
    axs1.grid(False)
    axs1.set_xlim(xlims)
    axs1.set_ylim(ylims)

    curves = [None for i in range(4)]
    for i, each in enumerate(flavours):
        curves[i] = axs1.plot(
            xs,
            each,
            'D-',
            markersize=3.0,
            color=colours[i]
        )[0]

    axs1.legend(curves, labels, fontsize=MainFontSize, loc="upper right")
    fig1.savefig(fileplot, fmt=fmt, bbox_inches='tight')

# Functions that makes a plot
def PlotAverageNuNubars(xs, xlims, xlabel, nunubars, ylims, fileplot, fmt="eps", dims=defaultDims):
    fig1 = Figure(figsize=dims)
    FigureCanvas(fig1)

    axs1 = fig1.add_subplot(1, 1, 1)
    axs1.set_xlabel(xlabel, fontsize=MainFontSize)
    #axs1.set_ylabel(r"$\langle \bar{P}_{\nu} \rangle / \langle \bar{P}_{\bar{\nu}} \rangle$", fontsize=MainFontSize)
    axs1.grid(False)
    axs1.set_xlim(xlims)
    axs1.set_ylim(ylims)

    axs1.plot(
        xs,
        nunubars,
        'D--',
        markersize=3.0,
        color="black"
    )

    fig1.savefig(fileplot, fmt=fmt, bbox_inches='tight')

#cmaps = [cm.get_cmap('viridis')]

def PlotSurface(xs, xlims, ys, ylims, zs, zlims, flav, fileplot, fmt="eps", dims=defaultDims):
    fig = Figure(figsize=dims)
    FigureCanvas(fig)

    axs = fig.add_subplot(111, projection='3d')
    axs.set_xlabel(r"$g_{+}$", fontsize=MainFontSize)
    axs.set_ylabel(r"$\mu / \omega$", fontsize=MainFontSize)
    axs.set_xlim(xlims)
    axs.set_ylim(ylims)
    axs.set_zlim(zlims)

    points = []; values = []
    for ix,x in enumerate(xs):
        for iy,y in enumerate(ys):
            points.append([x,y])
            values.append(zs[iy][ix])

    xArray = np.linspace(xlims[0], xlims[1], 100)
    yArray = np.linspace(ylims[0], ylims[1], 100)
    xArray, yArray = np.meshgrid(xArray, yArray)
    zArray = griddata(points, values, (xArray,yArray), method="linear")
    axs.plot_surface(xArray, yArray, zArray, cmap=viridis)

    #X, Y = np.meshgrid(xs, ys)
    #Z    = np.array(zs)
    #axs.plot_surface(X, Y, Z, cmap=viridis)

    fig.savefig(fileplot, fmt=fmt, bbox_inches='tight')