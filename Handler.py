#!/usr/bin/env python3

import numpy as np
import os
import sys
from Modules import FourPlots1D_PNG, Grids, Data, Rare, ObtainParameters

from scipy.optimize import curve_fit

# Aveage probabilities channel by channel
def Averages(data):
    # Average values of the left beam
    avProbsL = [[sum(line)/len(line) for line in mesh] for mesh in data[:4]]
    avProbsR = [[sum(line)/len(line) for line in mesh] for mesh in data[4:]]

    return avProbsL, avProbsR

# Covariances at some given height
def Covariances(data, iLine):
    # Pick the needed lines at some height z, possition at iLine index
    lines = [mesh[iLine] for mesh in data]
    # Evaluate the means probabilities
    means = [sum(each)/len(each) for each in lines]

    # Cycled index
    N = len(lines[0])
    cycled = lambda i: ((i % N) + N) % N

    covs = []
    for l,m in zip(lines, means):
        covs.append([ np.abs(sum([(l[j] - m)*(l[cycled(j+i)] - m) for j in range(N)]) / N) for i in range(N)])

    return covs

if __name__ == "__main__":
    directory, pX, pZ = ObtainParameters()

    # Let's make a directory fot the plots
    avDir = directory + "averages/"
    if not(os.path.exists(avDir)):
        os.mkdir(avDir)    

    # Obtain grids
    XGrid, ZGrid, D_z, D_x = Grids(directory)

    # Obtain data
    print("Obtaining the data...")
    data = Data(directory, D_x, D_z, pX, pZ)

    # Make rare meshes
    XGridRare, ZGridRare = Rare(XGrid, ZGrid, pX, pZ)

    print("Drawing the average probabilities...")
    avProbsL, avProbsR = Averages(data) # Evalaue them and draw

    names  = (r"$\nu_e$", r"$\nu_x$", r"$\bar{\nu}_e$", r"$\bar{\nu}_x$")
    labelsL = [n + r" (L)" for n in names]
    labelsR = [n + r" (R)" for n in names]

    FourPlots1D_PNG(ZGridRare, avProbsL, avDir + "avsL.png", "Average probabilities (Left)",  labelsL)
    FourPlots1D_PNG(ZGridRare, avProbsR, avDir + "avsR.png", "Average probabilities (Right)", labelsR)

    print("Drawing the covariances...")

    # Evaluate the covariances at two different positions
    covs10kmRaw = Covariances(data, int(len(ZGridRare) * 10.0/50.0))
    covs40kmRaw = Covariances(data, int(len(ZGridRare) * 40.0/50.0))

    # Change the order of the halves of the meshes of covariances
    # in order to set the central point of the covariances in the middle
    # of the plots
    covs10km = [mesh[int(len(XGridRare)/2):] + mesh[:int(len(XGridRare)/2)] for mesh in covs10kmRaw]
    covs40km = [mesh[int(len(XGridRare)/2):] + mesh[:int(len(XGridRare)/2)] for mesh in covs40kmRaw]
    # And make a suitable x grid 
    XGridMiddle = [x - XGrid[-1]/2.0 for x in XGridRare]

    FourPlots1D_PNG(XGridMiddle, covs10km[:4], avDir + "covs10kmL.png", "Covariances at z = 10 km (Left)",  labelsL, False)
    FourPlots1D_PNG(XGridMiddle, covs10km[4:], avDir + "covs10kmR.png", "Covariances at z = 10 km (Right)", labelsR, False)

    FourPlots1D_PNG(XGridMiddle, covs40km[:4], avDir + "covs40kmL.png", "Covariances at z = 40 km (Left)",  labelsL, False)
    FourPlots1D_PNG(XGridMiddle, covs40km[4:], avDir + "covs40kmR.png", "Covariances at z = 40 km (Right)", labelsR, False)

    #func = lambda x,a,s: a * np.exp(-x**2 / (2.0*s))
    #kekw  = curve_fit(func, XGridMiddle, covs10km[1])
    #print(kekw)