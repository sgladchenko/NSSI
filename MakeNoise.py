#!/usr/bin/env python3

import json
import random
import numpy as np

import sys

# Convert format from old one to the new one
def MakeFromOld(filename="./Noise.json"):
	dictNew = {"Meta":{}, "Harmonics":{}}

	with open("./NoiseOld.json") as f:
		dictOld = json.load(f)

		dictNew["Meta"]["N_Noise"] = int(dictOld["Meta"]["N_perturb"])
		dictNew["Meta"]["sigma"] = float(dictOld["Meta"]["sigma"])

		labels = "sinCoefficientsLeft", "cosCoefficientsLeft", "sinCoefficientsRight", "cosCoefficientsRight"
		h = "Harmonics"

		for lab in labels:
			dictNew[h][lab] = [[float(y) for y in x] for x in dictOld[h][lab]]

	with open(filename, "w") as f:
		json.dump(dictNew, f, indent=4)

def MakeNoise(N_Noise, sigma, filename="./Noise.json"):
	# Note: We should keep in mind that we originally aim to generate noise with
	# some kind of fixed norm in L_2([0, X]) space. And we aim that in average this noise should not exceed
	# the mean norm in L_2([0,X]) more than L \times sigma^2, which should not depend on the number of harmonics,
	# i.e. N_perturb

	# This means that we have to rescale the standard deviations in the generation of random Fourier coefficients
	# below as \sigma_R = \sigma / \sqrt{N_perturb \times 15}, in order to achieve (in average)

	#  || \delta \vec{\rho}_{\pm}(x) ||^2 ~ L/2 \times (\sigma_R^2 + 15 \times N_perturb \times \sigma_R^2) ~ 
	#                                     ~ L/2 \times \sigma^2

	# Desired rescaling
	sigma_R = sigma / (np.sqrt(N_Noise) * 15)

	# Generate 2*N harmonics for the left beam
	sinCoeffsLeft = [[random.gauss(0.0, sigma_R) for f in range(15)] for k in range(N_Noise)]
	cosCoeffsLeft = [[random.gauss(0.0, sigma_R) for f in range(15)] for k in range(N_Noise)]

	# Generate 2*N harmonics for the right beam
	sinCoeffsRight = [[random.gauss(0.0, sigma_R) for f in range(15)] for k in range(N_Noise)]
	cosCoeffsRight = [[random.gauss(0.0, sigma_R) for f in range(15)] for k in range(N_Noise)]

	# Dump into the file
	with open(filename, "w") as out:
		dictData = {"Meta": {}, "Harmonics": {}}

		# Save meta
		dictData["Meta"]["N_Noise"] = N_Noise
		dictData["Meta"]["sigma"]   = sigma

		# Save harmonics
		dictData["Harmonics"]["sinCoeffsLeft"] = sinCoeffsLeft
		dictData["Harmonics"]["cosCoeffsLeft"] = cosCoeffsLeft

		dictData["Harmonics"]["sinCoeffsRight"] = sinCoeffsRight
		dictData["Harmonics"]["cosCoeffsRight"] = cosCoeffsRight

		json.dump(dictData, out, indent=4)

if __name__ == "__main__":
	# Obtaining them from a command line

	# Default ones
	N_Noise = 0
	sigma = 0.05

	if len(sys.argv) > 1:
		for each in sys.argv[1:]:
			# Here I am parisng only nong options
			if '=' in each:
				key, val = each.split("=")[0].replace("--",""), each.split("=")[1]

				if key == "N_Noise":
					N_Noise = int(val)
				elif key == "sigma":
					sigma = float(val)
				else:
					print("Undefined key")
				
				MakeNoise(N_Noise, sigma)
			else:
				if each == "--old":
					MakeFromOld()

	MakeNoise(N_Noise, sigma)