#!/usr/bin/env python3

import numpy as np
import os, json, struct
import warnings

# The size of the types used for operating over the
# floating point numbers (in bytes)
SIZE_FLOAT = 8

class Data:
    """ Reading dumps of the calculated setups """

    def __init__(self, folder, periodN_x=1, periodN_z=1):
        # The folder in which the subfolder bin
        # must be found, containing the files bin/left.bin and bin/right.bin
        self.folder = folder.rstrip('/')

        if not os.path.isdir(f"{folder}/bin"):
            raise FileNotFoundError(f"No 'bin' directory found at {self.folder}")
        
        if not os.path.isfile(f"{folder}/bin/left.bin"):
            raise FileNotFoundError(f"No 'left.bin' file found at {self.folder}/bin")
        
        if not os.path.isfile(f"{folder}/bin/right.bin"):
            raise FileNotFoundError(f"No 'right.bin' file found at {self.folder}/bin")

        # Then check whether the configuration JSONs are found
        if not os.path.isfile(f"{folder}/XGrid.txt"):
            raise FileNotFoundError(f"No 'XGrid.txt' file found at {self.folder}")

        if not os.path.isfile(f"{folder}/ZGrid.txt"):
            raise FileNotFoundError(f"No 'ZGrid.txt' file found at {self.folder}")

        # Main buffer containing the probabilities from the binary files
        self._Probabilities = [[] for i in range(8)]
        self.loadedBin = False

        # Buffer containing (evaluated) average probabilities
        self._Averages = [None for i in range(8)]
        self.loadedAverages = False

        # First step is to obtain the dimensions of the saved grids
        # and the particular poiints at the axes
        with open(f"{folder}/XGrid.txt") as f:
            grid = f.read().rstrip().split(" ")
            self.XGrid = [float(each) for each in grid]
            self.Length_x_bin = len(self.XGrid)
            self.N_x_bin = self.Length_x_bin - 1

        with open(f"{folder}/ZGrid.txt") as f:
            grid = f.read().rstrip().split(" ")
            self.ZGrid = [float(each) for each in grid]
            self.Length_z_bin = len(self.ZGrid)
            self.N_z_bin = self.Length_z_bin - 1

        # Note: Of course, N_{axes}_bin mean the number of unique points
        # (as it is done in the numerical scheme), thus N_{axes}_bin = len(...) - 1,
        # and not just len(...) (i.e. excluding the copy of the zeroth point at the very end)
        # These N's can be interpreted as numbers of the steps between the points

        # Save periods and check whether they are divisors of the lengths of the lines
        self.periodN_x = periodN_x
        self.periodN_z = periodN_z

        if self.N_x_bin % periodN_x != 0:
            raise ValueError(f"Period of the x axis ({periodN_x}) is not a divisor of the horizontal dimension ({self.N_x_bin})")

        if self.N_z_bin % periodN_z != 0:
            raise ValueError(f"Period of the z axis ({periodN_z}) is not a divisor of the vertical dimension ({self.N_z_bin})")

        # Dimensions of the final displayed grid (also excluding the last point that is retained in the grids)
        self.N_x_displayed = self.N_x_bin // self.periodN_x
        self.N_z_displayed = self.N_z_bin // self.periodN_z

        # And including the last point -- these are the real dimensions of the matrices with grids
        # to be displayed via mpl
        self.Length_x_displayed = self.N_x_displayed + 1
        self.Length_z_displayed = self.N_z_displayed + 1

        # Finally, the displayed grids are
        self.XGrid_displayed = self.XGrid[::self.periodN_x]
        self.ZGrid_displayed = self.ZGrid[::self.periodN_z]

    # Access to the solution grids of probabilities
    @property
    def Probabilities(self):
        if not self.loadedBin:
            raise RuntimeError(f"The binary data haven't been loaded yet at {str(self)}.")
        else:
            return self._Probabilities

    @property
    def Averages(self):
        if not self.loadedAverages:
            raise RuntimeError(f"The averages haven't been loaded yet at {str(self)}.")
        else:
            return self._Averages

    # Same special function that actually is a generator
    # yielding the lines of solution one-by-one
    def lines(self):
        # Here we are going to directly save only the needed entries from the binaries
        fl = open(f"{self.folder}/bin/left.bin", "rb")
        fr = open(f"{self.folder}/bin/right.bin", "rb")

        # Offset from the beginning of each file where the desired entries are saved
        def offset(iz_displayed, ix_displayed):
            return ((iz_displayed*self.periodN_z)*self.Length_x_bin + ix_displayed*self.periodN_x)*4*SIZE_FLOAT

        # And then let's pick out exactly the desired values
        for iz in range(self.Length_z_displayed):
            line = [np.empty(self.Length_x_displayed, np.float64) for i in range(8)]

            # Evaluate the elements of the line
            for ix in range(self.Length_x_displayed):
                # Move to the needed entry
                fl.seek(offset(iz, ix))
                fr.seek(offset(iz, ix))
                # Read at each file 4 floats containing the probabilities of each of the
                # 4 floavours (in total: 8 floats) and let's unpack then simultaneously
                probsl = struct.unpack("dddd", fl.read(4*SIZE_FLOAT))
                probsr = struct.unpack("dddd", fr.read(4*SIZE_FLOAT))
                probs = probsl + probsr
                # And finally, let's save them in the matrices
                for i in range(8):
                    line[i][ix] = probs[i]

            # Finally, return this part of the solution
            yield line
        
        fl.close()
        fr.close()

    # Load the data contained in the binaries left.bin and right.bin
    def LoadBin(self):
        # Firsly, let's initialize the empty matrices
        for i in range(8):
            self._Probabilities[i] = np.empty((self.Length_z_displayed, self.Length_x_displayed), np.float64)

        # Let's use the generator for reading the binaries line by line
        for iline, line in enumerate(self.lines()):
            for i in range(8):
                self._Probabilities[i][iline] = line[i]

        # Toggle the flag that the binaries are loaded at this instance
        self.loadedBin = True

    # Evaluate the average probabilities (dependent on the z coordinate)
    def EvaluateAverages(self):
        # Evaluate the averages
        self._Averages = [[sum(line)/len(line) for line in sol] for sol in self.Probabilities]
        # And toggle the flag that the average probabilities have been evaluated
        self.loadedAverages = True

    # Dump the average probabilities into a folder
    def DumpAverages(self):
        # If the folder doesn't yet exists, make it
        averagesfolder = f"{self.folder}/averages"
        if not os.path.isdir(averagesfolder):
            os.mkdir(averagesfolder)
        
        # Then search for the files with the average probabilities
        # and warn that the following will rewrite the average data
        if os.path.isfile(f"{averagesfolder}/avsL.json"):
            warnings.warn(f"File {averagesfolder}/avsL.json has been found; it's going to be re-written", RuntimeWarning)
        
        if os.path.isfile(f"{averagesfolder}/avsR.json"):
            warnings.warn(f"File {averagesfolder}/avsR.json has been found; it's going to be re-written", RuntimeWarning)

        # Save them in the JSON files
        with open(f"{averagesfolder}/avsL.json", "w") as f:
            json.dump({"ZGridRare": self.ZGrid_displayed, "avs": [list(arr) for arr in self.Averages[:4]]}, f, indent=4)

        with open(f"{averagesfolder}/avsR.json", "w") as f:
            json.dump({"ZGridRare": self.ZGrid_displayed, "avs": [list(arr) for arr in self.Averages[4:]]}, f, indent=4)

    # The method that reads in a lazy way the binary data  line by line an evaluates
    # for each line the average probabilities
    def LazyEvaluateAverages(self):
        # Allocate the empty matrices
        self._Averages = [np.empty(self.Length_z_displayed, np.float64) for i in range(8)]

        # Read the lines and save the average probabilities
        for iline, line in enumerate(self.lines()):
            for i in range(8):
                self._Averages[i][iline] = sum(line[i]) / len(line[i])

        # And toggle flag that now we possess the averages probabilities
        self.loadedAverages = True

if __name__ == "__main__":
    data = Data("build/Data/NSSI NLM Mon Aug  2 14:18:26 2021/", periodN_x=1)
    #data.EvaluateAverages()
    #data.LoadBin()
    #print(data.Probabilities[0])
    #data.DumpAverages()
    #print(data[0])
    data.LazyEvaluateAverages()
    print(data.Averages[0])