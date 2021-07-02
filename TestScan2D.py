#!/usr/bin/env python3

import unittest, os, random, string, json, shutil
from Scan2D import Setup, makeMatrices, WrongNumberPoints
from Scan2D import Scan2D

randomString = lambda N: "".join(random.choices(string.ascii_lowercase + string.digits, k=N))
func = lambda mu, gPlus, k: (mu + gPlus)*(k+1)

def makeTestAverages(folder, mu, gPlus):
    with open(f"{folder}/averages/avsL.json", "w") as fout:
        json.dump({"ZGridRare": [float(k) for k in range(10)], "avs": [[func(mu, gPlus, f) for k in range(10)] for f in range(4)]}, fout)
    with open(f"{folder}/averages/avsR.json", "w") as fout:
        json.dump({"ZGridRare": [float(k) for k in range(10)], "avs": [[func(mu, gPlus, f) for k in range(10)] for f in range(4)]}, fout)


class TestScan2D(unittest.TestCase):

    def test_setup(self):
        setup = Setup("./TestData/NSSI NLM Test")
        self.assertEqual(setup.folder, "./TestData/NSSI NLM Test")
        self.assertEqual(setup.gPlus, 0.0)
        self.assertEqual(setup.mu, 40.0)
        self.assertEqual(setup.probs, [func(40.0, 0.0, f) for f in range(4)] + [3.0/7.0, ])

    def test_makeMatrices(self):
        setups = []
        mus = [10.0, 20.0, 30.0, 40.0, 50.0]
        gPluses = [0.0, 0.25, 0.50, 0.75, 1.0]

        for mu in mus:
            for gPlus in gPluses:
                s = Setup("./TestData/NSSI NLM Test")
                s.probs = [func(mu, gPlus, k) for k in range(5)]
                s.folder = randomString(10)
                s.mu = mu; s.gPlus = gPlus
                setups.append(s)
        random.shuffle(setups)

        self.assertEqual(makeMatrices(setups, mus, gPluses, "testentry", 1.0),
                         [[[func(mu, gPlus, k) for gPlus in gPluses] for mu in mus] for k in range(5)]) 
        self.assertRaises(WrongNumberPoints, makeMatrices, setups[:-2], mus, gPluses, "testentry", 1.0)

    def test_Scan2D(self):
        # Initialisation part
        testfolder = "./TestData/ScanTest"
        os.mkdir(testfolder)
        os.mkdir(f"{testfolder}/NH")
        os.mkdir(f"{testfolder}/IH")

        try:
            NHFolders = []
            IHFolders = []
            
            mus = [10.0, 20.0, 30.0, 40.0, 50.0]
            gPluses = [0.0, 0.25, 0.50, 0.75, 1.0]

            for mu in mus:
                for gPlus in gPluses:
                    folderNH = f"{testfolder}/NH/NSSI NLM {randomString(15)}"; os.mkdir(folderNH); NHFolders.append(folderNH)
                    folderIH = f"{testfolder}/IH/NSSI NLM {randomString(15)}"; os.mkdir(folderIH); IHFolders.append(folderIH)

                    # Generate Parameters.json files
                    with open("./Parameters.json") as fin:
                        d = json.load(fin)
                        d["V-A & NSSI & MSW & AMM"]["n_Nu"][0] = mu * (1.0e29/13.976129602265964)
                        d["V-A & NSSI & MSW & AMM"]["g_{+}"] = gPlus
                        dNH = d; dNH["eta"] =  1.0
                        dIH = d; dIH["eta"] = -1.0
                    
                    with open(f"{folderNH}/Parameters.json", "w") as fout: json.dump(dNH, fout)
                    with open(f"{folderIH}/Parameters.json", "w") as fout: json.dump(dIH, fout)

                    # Generate then averages
                    os.mkdir(f"{folderNH}/averages")
                    os.mkdir(f"{folderIH}/averages")

                    makeTestAverages(folderNH, mu, gPlus)
                    makeTestAverages(folderIH, mu, gPlus)

            # Finally, test section of the Scan2D class
            scan2d = Scan2D(testfolder)
            self.assertEqual(sorted(scan2d.NHFolders), sorted(NHFolders))
            self.assertEqual(sorted(scan2d.IHFolders), sorted(IHFolders))
            self.assertEqual(scan2d.NHmus, mus)
            self.assertEqual(scan2d.IHmus, mus)
            self.assertEqual(scan2d.NHgPluses, gPluses)
            self.assertEqual(scan2d.IHgPluses, gPluses)

            self.assertEqual(scan2d.NHMatrices, [[[func(mu, gPlus, f) for gPlus in gPluses] for mu in mus] for f in range(4)] + [[[3.0/7.0 for gPlus in gPluses] for mu in mus],])
            self.assertEqual(scan2d.IHMatrices, [[[func(mu, gPlus, f) for gPlus in gPluses] for mu in mus] for f in range(4)] + [[[3.0/7.0 for gPlus in gPluses] for mu in mus],])

        finally:        
            # Remove test folders
            shutil.rmtree(testfolder)

if __name__ == "__main__":
    if "TestData" not in os.listdir():
        os.mkdir("./TestData/")
        # And make dummy setup folders and files
        os.mkdir("./TestData/NSSI NLM Test/")
        os.mkdir("./TestData/NSSI NLM Test/averages")

        with open("Parameters.json") as fin:
            d = json.load(fin)
            d["V-A & NSSI & MSW & AMM"]["n_Nu"][0] = 40.0 * (1.0e29/13.976129602265964)
            d["V-A & NSSI & MSW & AMM"]["g_{+}"] = 0.0
        with open("./TestData/NSSI NLM Test/Parameters.json", "w") as fout: json.dump(d, fout)
        makeTestAverages("./TestData/NSSI NLM Test", 40.0, 0.0)

    unittest.main()