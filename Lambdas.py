import numpy as np

# Super Ultra Turbo basis in su(4) Lie algebra
# It's actually even an orthonormal basis

lambda1 = np.matrix([[0.0, 1.0, 0.0, 0.0],
					 [1.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])

lambda2 = np.matrix([[0.0, -1.0j, 0.0, 0.0],
					 [1.0j, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])


lambda3 = np.matrix([[1.0, 0.0, 0.0, 0.0],
					 [0.0, -1.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])

lambda4 = np.matrix([[0.0, 0.0, 1.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [1.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])

lambda5 = np.matrix([[0.0, 0.0, -1.0j, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [1.0j, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])

lambda6 = np.matrix([[0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 1.0, 0.0],
					 [0.0, 1.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])

lambda7 = np.matrix([[0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, -1.0j, 0.0],
					 [0.0, 1.0j, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]])

lambda8 = np.matrix([[1.0, 0.0, 0.0, 0.0],
					 [0.0, 1.0, 0.0, 0.0],
					 [0.0, 0.0, -2.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0]]) / np.sqrt(3)

lambda9 = np.matrix([[0.0, 0.0, 0.0, 1.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [1.0, 0.0, 0.0, 0.0]])

lambda10 = np.matrix([[0.0, 0.0, 0.0, -1.0j],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0],
					  [1.0j, 0.0, 0.0, 0.0]])

lambda11 = np.matrix([[0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 1.0],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 1.0, 0.0, 0.0]])

lambda12 = np.matrix([[0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, -1.0j],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 1.0j, 0.0, 0.0]])

lambda13 = np.matrix([[0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 1.0],
					  [0.0, 0.0, 1.0, 0.0]])

lambda14 = np.matrix([[0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, 0.0],
					  [0.0, 0.0, 0.0, -1.0j],
					  [0.0, 0.0, 1.0j, 0.0]])

lambda15 = np.matrix([[1.0, 0.0, 0.0, 0.0],
					  [0.0, 1.0, 0.0, 0.0],
					  [0.0, 0.0, 1.0, 0.0],
					  [0.0, 0.0, 0.0, -3.0]]) / np.sqrt(6)

lambdas = [lambda1, lambda2,  lambda3,  lambda4,  lambda5,  lambda6,  lambda7,  lambda8,
		   lambda9, lambda10, lambda11, lambda12, lambda13, lambda14, lambda15]

# Orthogonal; it means 

su4Basis = [x / np.sqrt(2) for x in lambdas]

# Pauli matrices
Id2x2   = np.matrix([[1.0,  0.0],  [0.0,  1.0]])
sigma_1 = np.matrix([[0.0,  1.0],  [1.0,  0.0]])
sigma_2 = np.matrix([[0.0, -1.0j], [1.0j, 0.0]])
sigma_3 = np.matrix([[1.0,  0.0],  [0.0, -1.0]])
sigmas  = (sigma_1, sigma_2, sigma_3)

# A big 4x4 identity matrix
Id4x4   = np.matrix(np.diag([1.0, 1.0, 1.0, 1.0]))

# Some useful matrices neded  in the construction of the total effective Hamiltonians
AMMLike = np.matrix([[ 0.0,   0.0,  0.0, -1.0j],
					 [ 0.0,   0.0, 1.0j,   0.0],
					 [ 0.0, -1.0j,  0.0,   0.0],
					 [1.0j,   0.0,  0.0,   0.0]])

flavourSigma1 = np.matrix([[0.0, 1.0, 0.0, 0.0],
						   [1.0, 0.0, 0.0, 0.0],
						   [0.0, 0.0, 0.0, 1.0],
						   [0.0, 0.0, 1.0, 0.0]])

flavourSigma3 = np.matrix(np.diag([1.0, -1.0, 1.0, -1.0]))

G = np.matrix(np.diag([1.0, 1.0, -1.0, -1.0]))