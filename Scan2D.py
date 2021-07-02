#!/usr/bin/env python3

from Scan import Scan, extractMu, extractgPlus, extractProbs, labels, filelabels
from Modules import PlotAverageFlavours, PlotAverageNuNubars, PlotSurface
import os

class WrongNumberPoints(Exception):
    def __init__(self, dir, eta, NFound, NDesired):
        self.NFound = NFound
        self.NDesired = NDesired
        hierarchy = "NH" if eta == 1.0 else "IH"
        self.message = f"At {dir}/{hierarchy} {NFound} points found ({NDesired} expected)"
    
    def __str__(self):
        return self.message

class Setup:
    """
        A little class that contains the information
        about a calculated setup (it saves mu and gPlus).
    """
    def __init__(self, setupfolder, roundMu=2):
        self.mu     = round(extractMu(setupfolder), roundMu)
        self.gPlus  = extractgPlus(setupfolder)
        self.probs  = extractProbs(setupfolder)
        self.folder = setupfolder
        # Note: roundMu is needed when mu values are a bit ugly
        # and we don't want to see all the nonzero digits in its
        # decimal form

    def __lt__(self, stp):
        if self.mu < stp.mu:
            return True
        if self.mu > stp.mu:
            return False
        # And if we're here, mus are equal,
        # so we're to just compare the gPlus
        return self.gPlus < stp.gPlus

def makeMatrices(setups, mus, gPluses, dir, eta):
    # This function merely makes the matrices of probabilities
    # that can be found setups
    # If there's not enough points to make the whole matrix
    # or for a given pair of parameters there're more than one
    # setups, it'll tell there's an error in the data
    probsMatrices = [[[0.0 for g in gPluses] for mu in mus] for f in range(5)]

    # Check whether we're indeed able to assemble the whole matrix
    # (so for each pair we have only one evaluated setup)
    Nmus = len(mus)
    NgPluses = len(gPluses)
    if len(setups) != Nmus * NgPluses:
        raise WrongNumberPoints(dir, eta, len(setups), Nmus * NgPluses)

    # Then let's sort them and recover from this sorted serialised form
    # the entire matrices
    sortedSetups = sorted(setups)
    for i in range(len(mus)):
        for j in range(len(gPluses)):
            for f in range(5):
                probsMatrices[f][i][j] = sortedSetups[NgPluses*i + j].probs[f]

    return probsMatrices

def fixFirst(matrices, i):
    # Just returns a slice of the matrix
    # with elements laying at [i][whatever] (i.e. at the ith row)
    return [matrices[f][i] for f in range(5)]

def fixSecond(matrices, j):
    # Just returns a slice of the matrix
    # with elements laying at [whatever][j] (i.e. at the jth column)
    return [[line[j] for line in matrices[f]] for f in range(5)]


# Plot a slice of the two-dimensional scan with a fixed mu value
def PlotMuFixed(matrices, imu, mus,
                          gPluses, gPlusLims,
                          probLims, nunubarLims,
                          dir, filetitle, hierarchy):
    curves = fixFirst(matrices, imu)
    print(f"The chosen value of mu/omega is {mus[imu]} ({hierarchy} chosen)")

    flavoursfn = f"{dir}/{filetitle}_Flavours_at_mu={mus[imu]}_({hierarchy}).eps"
    nunubarsfn = f"{dir}/{filetitle}_NuNubars_at_mu={mus[imu]}_({hierarchy}).eps"
    print(f"The following files will be generated:\n[*] {flavoursfn}\n[*] {nunubarsfn}")

    PlotAverageFlavours(gPluses, gPlusLims, r"$g_{+}$",
                        curves[:4], probLims,
                        labels, flavoursfn)

    PlotAverageNuNubars(gPluses, gPlusLims, r"$g_{+}$",
                        curves[4], nunubarLims,
                        nunubarsfn)

    print(f"Success.")

# Plot a slice of the two-dimensional scan with a fixed gPlus value
def PlotgPlusFixed(matrices, mus, muLims,
                             gPluses, igPlus,
                             probLims, nunubarLims,
                             dir, filetitle, hierarchy):
    curves = fixSecond(matrices, igPlus)
    print(f"The chosen value of gPlus is {gPluses[igPlus]} ({hierarchy} chosen)")

    flavoursfn = f"{dir}/{filetitle}_Flavours_({hierarchy})_at_gPlus={gPluses[igPlus]}.eps"
    nunubarsfn = f"{dir}/{filetitle}_NuNubars_({hierarchy})_at_gPlus={gPluses[igPlus]}.eps"
    print(f"The following files will be generated:\n[*] {flavoursfn}\n[*] {nunubarsfn}")

    PlotAverageFlavours(mus, muLims, r"$\mu / \omega$",
                        curves[:4], probLims,
                        labels, flavoursfn)

    PlotAverageNuNubars(mus, muLims, r"$\mu / \omega$",
                        curves[4], nunubarLims, nunubarsfn)

    print(f"Success.")

# Plots a a surface of values lying at flav inside each setup.probs
def PlotSurfaces(matrices, mus, muLims,
                          gPluses, gPlusLims,
                          probLims, nunubarLims,
                          dir, filetitle, hierarchy, fmt="eps"):
    print(f"Plotting surface of the {hierarchy}")

    fns = [f"{dir}/{filetitle}_{filelabels[f]}_{hierarchy}.{fmt}" for f in range(5)]
    print(f"The following files will be generated:")
    print("[*] " + "\n[*] ".join(fns))

    for f in range(5):
        PlotSurface(xs=gPluses,     xlims=gPlusLims,
                    ys=mus,         ylims=muLims,
                    zs=matrices[f], zlims=probLims if f!=4 else nunubarLims, flav=f,
                    fileplot=fns[f], fmt=fmt)
    
    print("Success.")

# Draws all the plots
def PlotAll(matrices, mus, muLims,
                      gPluses, gPlusLims,
                      probLims, nunubarLims,
                      dir, filetitle, hierarchy):
    # All the fixed-mu slices
    for imu in range(len(mus)):
        PlotMuFixed(matrices, imu, mus, gPluses, gPlusLims, probLims, nunubarLims, dir, filetitle, hierarchy)
    # All the fixed-gPlus slices
    for igPlus in range(len(gPluses)):
        PlotgPlusFixed(matrices, mus, muLims, gPluses, igPlus, probLims, nunubarLims, dir, filetitle, hierarchy)

    # All the surfaces
    PlotSurfaces(matrices, mus, muLims, gPluses, gPlusLims, probLims, nunubarLims, dir, filetitle, hierarchy)

class Scan2D(Scan):
    """
        Object that reads the two-dimensional scan 
        (in mu/g_+ axes) and draws the plots.
    """

    def __init__(self, dir):
        Scan.__init__(self, dir)
        # In further we are presuming that in the target dir we've
        # come up with two-dimensional scan in mu and gPlus parameters
        # and each of the pair (mu,gPlus) is enountered just once.
        self.NHSetups = [Setup(f) for f in self.NHFolders]
        self.IHSetups = [Setup(f) for f in self.IHFolders]

        # All the possible mus, in increasing order
        self.NHmus = sorted({each.mu for each in self.NHSetups})
        self.IHmus = sorted({each.mu for each in self.IHSetups})
        # All the possible gPluses, in increasing order
        self.NHgPluses = sorted({each.gPlus for each in self.NHSetups})
        self.IHgPluses = sorted({each.gPlus for each in self.IHSetups})

        # The matrices of probabilities
        self.NHMatrices = makeMatrices(self.NHSetups, self.NHmus, self.NHgPluses, self.dir, 1.0)
        self.IHMatrices = makeMatrices(self.IHSetups, self.IHmus, self.IHgPluses, self.dir, -1.0)

    def PlotData(self, muLims, gPlusLims,
                       probLims, nunubarLims, filetitle,
                       NHFlag=True, IHFlag=True):
        
        # First of all, let's make a directory if it doesn't exist
        plotsdir = f"{self.dir}/Scan2D_Plots"
        if "Scan2D_Plots" not in os.listdir(self.dir): os.mkdir(plotsdir)
        if "NH" not in os.listdir(plotsdir): os.mkdir(f"{plotsdir}/NH")
        if "IH" not in os.listdir(plotsdir): os.mkdir(f"{plotsdir}/IH")

        # Plot the normal hierarchy data
        if NHFlag:
            PlotAll(self.NHMatrices, self.NHmus, muLims,
                                     self.NHgPluses, gPlusLims,
                                     probLims, nunubarLims,
                                     f"{plotsdir}/NH", filetitle, "NH")
        # Plot the inverted hierarchy data
        if IHFlag:
            PlotAll(self.IHMatrices, self.IHmus, muLims,
                                     self.IHgPluses, gPlusLims,
                                     probLims, nunubarLims,
                                     f"{plotsdir}/IH", filetitle, "IH")


if __name__ == "__main__":
    scan2d = Scan2D("./Scans/Scan2D")
    scan2d.PlotData(muLims=[10.0, 50.0], gPlusLims=[0.0, 1.0],
                    probLims=[0.0, 0.6], nunubarLims=[1.0, 1.6], filetitle="10km_Scan2D")