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
    with os.scandir(dir) as entries:
        for each in entries:
            if each.is_dir() and each.name == folder:
                return True
        return False

# Scan for the folders with setups
def findSetups(dir):
    fullfolders = []
    with os.scandir(dir) as entries:
        for each in entries:
            if each.is_dir():
                if re.search(r"NSSI NLM .*", each.name):
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
filelabels = ["eNu", "xNu", "eANu", "xANu", "NuNuBar"]