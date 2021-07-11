import numpy as np

# Super Ultra Turbo basis in su(4) Lie algebra
# It's actually even an orthonormal basis

lambda1 = np.array([[0.0, 1.0, 0.0, 0.0],
					[1.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])

lambda2 = np.array([[0.0, -1.0j, 0.0, 0.0],
					[1.0j, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])


lambda3 = np.array([[1.0, 0.0, 0.0, 0.0],
					[0.0, -1.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])

lambda4 = np.array([[0.0, 0.0, 1.0, 0.0],
					[0.0, 0.0, 0.0, 0.0],
					[1.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])

lambda5 = np.array([[0.0, 0.0, -1.0j, 0.0],
					[0.0, 0.0, 0.0, 0.0],
					[1.0j, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])

lambda6 = np.array([[0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 1.0, 0.0],
					[0.0, 1.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])

lambda7 = np.array([[0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, -1.0j, 0.0],
					[0.0, 1.0j, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]])

lambda8 = np.array([[1.0, 0.0, 0.0, 0.0],
					[0.0, 1.0, 0.0, 0.0],
					[0.0, 0.0, -2.0, 0.0],
					[0.0, 0.0, 0.0, 0.0]]) / np.sqrt(3)

lambda9 = np.array([[0.0, 0.0, 0.0, 1.0],
					[0.0, 0.0, 0.0, 0.0],
					[0.0, 0.0, 0.0, 0.0],
					[1.0, 0.0, 0.0, 0.0]])

lambda10 = np.array([[0.0, 0.0, 0.0, -1.0j],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [1.0j, 0.0, 0.0, 0.0]])

lambda11 = np.array([[0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 1.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 1.0, 0.0, 0.0]])

lambda12 = np.array([[0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, -1.0j],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 1.0j, 0.0, 0.0]])

lambda13 = np.array([[0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 1.0],
					 [0.0, 0.0, 1.0, 0.0]])

lambda14 = np.array([[0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, 0.0],
					 [0.0, 0.0, 0.0, -1.0j],
					 [0.0, 0.0, 1.0j, 0.0]])

lambda15 = np.array([[1.0, 0.0, 0.0, 0.0],
					 [0.0, 1.0, 0.0, 0.0],
					 [0.0, 0.0, 1.0, 0.0],
					 [0.0, 0.0, 0.0, -3.0]]) / np.sqrt(6)

lambdas = [lambda1, lambda2,  lambda3,  lambda4,  lambda5,  lambda6,  lambda7,  lambda8,
		   lambda9, lambda10, lambda11, lambda12, lambda13, lambda14, lambda15]

# Orthogonal; it means 

su4Basis = [x / np.sqrt(2) for x in lambdas]

# Some useful matrices neded  in the construction of the total effective Hamiltonians
AMMLike = np.array([[ 0.0,   0.0,  0.0, -1.0j],
					[ 0.0,   0.0, 1.0j,   0.0],
					[ 0.0, -1.0j,  0.0,   0.0],
					[1.0j,   0.0,  0.0,   0.0]])

flavourSigma1 = np.array([[0.0, 1.0, 0.0, 0.0],
						  [1.0, 0.0, 0.0, 0.0],
						  [0.0, 0.0, 0.0, 1.0],
						  [0.0, 0.0, 1.0, 0.0]])

flavourSigma3 = np.diag([1.0, -1.0, 1.0, -1.0])

G = np.diag([1.0, 1.0, -1.0, -1.0])


# Returns the diagonal blocks
def su4Diag(m):
	return np.array([[m[0,0], m[0,1], 0.0,    0.0   ],
					 [m[1,0], m[1,1], 0.0,    0.0   ],
					 [0.0,    0.0,    m[2,2], m[2,3]],
					 [0.0,    0.0,    m[3,2], m[3,3]]])

# Returns the offdiagonal blocks
def su4Offdiag(m):
	return np.array([[0.0,     0.0,   m[0,2], m[0,3]],
					 [0.0,     0.0,   m[1,2], m[1,3]],
					 [m[2,0], m[2,1], 0.0,    0.0   ],
					 [m[3,0], m[3,1], 0.0,    0.0   ]])

# Evaluate the round matrix
def su4Round(m):
	# The matrix that corresponds to (rho^C)^T
	mCT = np.array([[m[2,2], m[3,2], m[0,2], m[1,2]],
		            [m[2,3], m[3,3], m[0,3], m[1,3]],
		            [m[2,0], m[3,0], m[0,0], m[1,0]],
		            [m[2,1], m[3,1], m[0,1], m[1,1]]])

	return m - mCT

# Simple commutator
def com(a,b):
    return a @ b - b @ a

# Find the biggest real part among the eigenvalues of a matrix
def BiggestRealEigPart(m):
    # Eigenvalues
    eigs = np.linalg.eig(m)[0]
    return np.max([np.real(e) for e in eigs])