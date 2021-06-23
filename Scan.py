#!/usr/bin/env python3

import os, json, re
from Constants import Constants
from Modules import PlotAverageFlavours, PlotAverageNuNubars

class NoAveragesFound(Exception):
    def __init__(self, setupfolder):
        self.message = f"Make sure to run ./Plots.py at {setupfolder}"
        self.folder  = setupfolder
    def __str__(self):
        return self.message

class NoFolderStructureFound(Exception):
    def __init__(self, dir):
        self.message = f"No folder structure found at {dir}"
        self.dir = dir
    def __str__(self):
        return self.message

# Says whether folder is in dir or not
def findFolder(dir, folder):
    for each in os.scandir(dir):
        if each.is_dir() and each.name == folder:
            return True
    return False

# Scan for the folders with setups
def findSetups(dir):
    fullfolders = []
    for each in os.scandir(dir):
        if each.is_dir():
            if re.findall(r"NSSI NLM .*", each.name):
                fullfolders.append(f"{dir}/{each.name}")
    return fullfolders

# Extract the asymptotic probabilities
# from a specified setup folder
def extractProbs(setupfolder, percentage=0.1):
    # First let's chech whether these averages have already been
    # evaluated and saved by ./Plots.py at a special directory
    if not(findFolder(setupfolder, "averages")):
        raise NoAveragesFound(setupfolder)
    averagesdir = setupfolder + "/averages"
    # Then let's extract the average probabilities;
    # in order to obtain clear asymptotic probabilities without possible
    # fluctuations they should be averaged as well; next we evaluate exactly
    # these mean probabilities at the termalisation stage, averaging them over 
    # last percentage of points
    probs = []
    for filename in ["avsL.json", "avsR.json"]:
        with open(f"{averagesdir}/{filename}") as f:
            d = json.load(f)
            numpoints = int(percentage * len(d["ZGridRare"]))
            probs.append([sum(flav[-numpoints:])/numpoints for flav in d["avs"]])
    # Average left-right probabilitie
    meanLR = [0.5*(probs[0][f] + probs[1][f]) for f in range(4)]
    # And nu-nubar ratio
    nunubar = (meanLR[0] + meanLR[1]) / (meanLR[2] + meanLR[3])
    return meanLR + [nunubar,]

# Extract the level of mu
def extractMu(setupfolder):
    # First, let's initialise the Constants object
    c = Constants(f"{setupfolder}/Parameters.json")
    return c.muOverOmega

# Extract the coupling constant gPlus
def extractgPlus(setupfolder):
    # First, let's initialise the Constants object
    c = Constants(f"{setupfolder}/Parameters.json")
    return c.gPlus

class Scan:
    """
        Base class that scans the output directories
        and save the asymptotic values of probabilities.
    """
    def __init__(self, dir):
        # Let's scan firstly the specified directory
        self.dir = dir
        if not(findFolder(dir, "NH") and findFolder(dir, "IH")):
            raise NoFolderStructureFound(dir)
        # Then let's search for the setup directories;
        # they all start with "NSSI NLM <...>",
        # where <...> is a time stamp
        self.NHFolders = findSetups(f"{dir}/NH")
        self.IHFolders = findSetups(f"{dir}/IH")
        # And then finally let's save the thermalied probabilities
        # for each setup
        self.NHProbs = [extractProbs(ff) for ff in self.NHFolders]
        self.IHProbs = [extractProbs(ff) for ff in self.IHFolders]

# Additional constants
labels = [r"$\nu_e$", r"$\nu_x$", r"$\bar{\nu}_e$", r"$\bar{\nu}_x$"]

class MuScan(Scan):
    """
        Collective potential mu scans.
    """
    def __init__(self, dir):
        Scan.__init__(self, dir)
        # Then here let's obtain the values
        # of gPlus in the setups
        self.NHmus = [extractMu(ff) for ff in self.NHFolders]
        self.IHmus = [extractMu(ff) for ff in self.IHFolders]
        # Then let's save everything pairwise
        self.NHs = list(zip(self.NHmus, self.NHProbs))
        self.IHs = list(zip(self.IHmus, self.IHProbs))
        # And let's sort them by the value of gPlus
        self.NHs = sorted(self.NHs, key=lambda each: each[0])
        self.IHs = sorted(self.IHs, key=lambda each: each[0])

    def plot(self, muLims=[10.0, 150.0], probLims=[0.0, 0.6], nunubarLims=[1.0, 1.6]):
        # For each hierarchy draw the plot
        flavFiles    = f"{self.dir}/NH_Flav_MuScan.eps",\
                       f"{self.dir}/IH_Flav_MuScan.eps"

        nunubarFiles = f"{self.dir}/NH_NuNubar_MuScan.eps",\
                       f"{self.dir}/IH_NuNubar_MuScan.eps"
        
        for eta, data in enumerate([self.NHs, self.IHs]):
            flavours = [[each[1][f] for each in data] for f in range(4)]
            nunubars = [each[1][4] for each in data]
            mus      = [each[0] for each in data]
            # And plot them
            PlotAverageFlavours(mus,
                                muLims,
                                r"$\mu / \omega$",
                                flavours,
                                probLims,
                                labels,
                                flavFiles[eta])

            PlotAverageNuNubars(mus,
                                muLims,
                                r"$\mu / \omega$",
                                nunubars,
                                nunubarLims,
                                nunubarFiles[eta])

class gPlusScan(Scan):
    """
        gPlus scans.
    """
    def __init__(self, dir):
        Scan.__init__(self, dir)
        # Then here let's obtain the values
        # of gPlus in the setups
        self.NHgPlus = [extractgPlus(ff) for ff in self.NHFolders]
        self.IHgPlus = [extractgPlus(ff) for ff in self.IHFolders]
        # Then let's save everything pairwise
        self.NHs = list(zip(self.NHgPlus, self.NHProbs))
        self.IHs = list(zip(self.IHgPlus, self.IHProbs))
        # And let's sort them by the value of gPlus
        self.NHs = sorted(self.NHs, key=lambda each: each[0])
        self.IHs = sorted(self.IHs, key=lambda each: each[0])

    def plot(self, gPlusLims=[0.0, 1.0], probLims=[0.0, 0.6], nunubarLims=[1.0, 1.6]):
        # For each hierarchy draw the plot
        flavFiles    = f"{self.dir}/NH_Flav_gPlusScan.eps",\
                       f"{self.dir}/IH_Flav_gPlusScan.eps"

        nunubarFiles = f"{self.dir}/NH_NuNubar_gPlusScan.eps",\
                       f"{self.dir}/IH_NuNubar_gPlusScan.eps"
        
        for eta, data in enumerate([self.NHs, self.IHs]):
            flavours = [[each[1][f] for each in data] for f in range(4)]
            nunubars = [each[1][4] for each in data]
            gpluses  = [each[0] for each in data]
            # And plot them
            PlotAverageFlavours(gpluses,
                                gPlusLims,
                                r"$g_{+}$",
                                flavours,
                                probLims,
                                labels,
                                flavFiles[eta])

            PlotAverageNuNubars(gpluses,
                                gPlusLims,
                                r"$g_{+}$",
                                nunubars,
                                nunubarLims,
                                nunubarFiles[eta])

if __name__ == "__main__":
    #scan = gPlusScan("./Scans/gPlusscan_1st")
    #scan.plot()
    scan = MuScan("./Scans/muscan_dummy")
    scan.plot()
