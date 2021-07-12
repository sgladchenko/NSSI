#!/usr/bin/env python3

import numpy as np

class Data:
    """ Reading the binaries from setup folders with
        possible bufferization """
    def __init__(self, setupfolder, D_x, D_z, pX, pZ):
        # The folder in which the subfolder /bin
        # must be found, with the files /bin/left.bin
        # and bin/right.bin
        self.setupfolder = setupfolder if setupfolder[-1] != "/" else setupfolder[:-1]

        # These lists are for the probabilities
        # in the left beam
        self.eNuL  = []
        self.xNuL  = []
        self.eANuL = []
        self.xANuL = []

        # These lists are for the probabilities
        # in the right beam
        self.eNuR  = []
        self.xNuR  = []
        self.eANuR = []
        self.xANuR = []

        # Dimensions of the grid saved in the binaries
        self.D_x = D_x
        self.D_z = D_z

        # The steps of the grids with which the probabilities
        # will be shown on the plots
        self.pX = pX
        self.pZ = pZ

        # The number of z-lines that have already been loaded in memory
        # in this instance
        self.savedZ = 0

    def read(self, divisor=10):
        """
        LData = np.fromfile(f"{self.setupfolder}/left.bin",
                            offset=self.savedZ*self.D_x*4,
                            count=numLines*self.D_x*4)
        
        RData = np.fromfile(f"{self.setupfolder}/right.bin",
                            offset=self.savedZ*self.D_x*4,
                            count=numLines*self.D_x*4)
        """
        pass
        
if __name__ == "__main__":
    pass