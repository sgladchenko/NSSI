#!/usr/bin/env python3

import numpy as np
from copy import deepcopy
from Lambdas import su4Diag, su4Offdiag, su4Round, com, BiggestRealEigPart, flavourSigma3, G

# In all the further fomulae zeta=+1 corresponds to the 'left' beam
# and zeta=-1 to the 'right' one; the matrices will be packed as (left,right)
# tuples.
#
# The Hamiltonians, and q,mu values are all the units of omega_vac.

# Send index to the position in the matrices
indexMapper_Offdiag = {0:(0,3), 1:(1,2), 2:(2,1), 3:(3,0), 4:(0,2), 5:(1,3), 6:(2,0), 7:(3,1)}
indexMapper_Diag    = {0:(0,0), 1:(1,1), 2:(2,2), 3:(3,3), 4:(0,1), 5:(1,0), 6:(2,3), 7:(3,2)}
indexMapper = (indexMapper_Offdiag, indexMapper_Diag)

# Send position in the matrices to an index
maskLeft = np.array([[16, 20, 4,  0 ],
                     [21, 17, 1,  5 ],
                     [6,  2,  18, 22],
                     [3,  7,  23, 19]])
maskRight = np.array([[24, 28, 12, 8 ],
                      [29, 25, 9,  13],
                      [14, 10, 26, 30],
                      [11, 15, 31, 27]])

class Pair:
    def __init__(self, left, right):
        self.pair = [deepcopy(left), deepcopy(right)]

    @property
    def left(self):
        return self.pair[0]
    
    @property
    def right(self):
        return self.pair[1]

    # Note: the indexing works as follows:
    #   0..7:   left beam, offdiagonal components
    #   8..15:  right beam, offdiagonal components
    #   16..23: left beam, diagonal components
    #   24..31: right beam, diagonal components

    # Returns a matrix element in the pair saved in this instance
    def __getitem__(self, index):
        # Left/right index
        leftrightIndex = (index // 8) % 2
        # Diag/offdiag index
        diagoffdiagIndex = index // 16
        # Matrix indices
        i,j = indexMapper[diagoffdiagIndex][index % 8]

        return self.pair[leftrightIndex][i,j]

    # Sets a matrix element at index position
    def __setitem__(self, index, value):
        # Left/right index
        leftrightIndex = (index // 8) % 2
        # Diag/offdiag index
        diagoffdiagIndex = index // 16
        # Matrix indices
        i,j = indexMapper[diagoffdiagIndex][index % 8]

        self.pair[leftrightIndex][i,j] = deepcopy(value)

    # Alternative constructor that sets only one nonzero matrix element
    # at the index position
    @classmethod
    def singleEntry(cls, index):
        tmp = cls(np.zeros((4,4)), np.zeros((4,4)))
        tmp[index] = 1.0
        return tmp

class Setup:
    def __init__(self, eta, q, mu, gPlus, chi, rhoInit):
        self.eta = eta
        self.q = q
        self.mu  = mu
        self.gPlus = gPlus
        self.tan = np.tan(chi)
        self.cos = np.cos(chi)
        self.rhoInit = rhoInit

# Returns a vacuum Hamiltonian in the units of omega_{vac}
def VacHam(setup):
    return -0.5 * setup.eta * flavourSigma3

# A function that states in the collective part of the total Hamiiltonian
def CollHam(setup, rho):
    SMPart   = np.trace(rho @ G)*G + su4Round(rho)
    NSSIPart = setup.gPlus * su4Offdiag(su4Round(rho)).T
    return setup.mu * (SMPart + NSSIPart)

# Returns a corresponding pair of z-derivatives in the linearized equations.
def L(setup, pair):
    # Let's subsequently define all the parts of the L operator
    # first: i*\zeta*q*tan(\chi) \delta \rho_{\zeta}
    lDerivativeX = ( 1.0j*setup.q*setup.tan)*pair.left
    rDerivativeX = (-1.0j*setup.q*setup.tan)*pair.right

    # 1st commutator
    hvac = VacHam(setup) 
    lCommutator1 = (-1.0j/setup.cos)*com(hvac + CollHam(setup, setup.rhoInit), pair.left)
    rCommutator1 = (-1.0j/setup.cos)*com(hvac + CollHam(setup, setup.rhoInit), pair.right)

    # 2nd commutator
    lCommutator2 = (-1.0j/setup.cos)*com(CollHam(setup, pair.right), setup.rhoInit)
    rCommutator2 = (-1.0j/setup.cos)*com(CollHam(setup, pair.left),  setup.rhoInit)

    # And finally, the Pair instance containg the matrices
    return Pair(lDerivativeX + lCommutator1 + lCommutator2,
                rDerivativeX + rCommutator1 + rCommutator2)