#!/usr/bin/env python3

import numpy as np
from copy import deepcopy
from Lambdas import su4Diag, su4Offdiag, su4Round, com, BiggestRealEigPart, flavourSigma3, G
from Modules import PlotStability
from Elapsed import elapsed

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

# Permutator of the indices
inds = [0, 1, 10, 11, 2, 3, 8, 9, 4, 5, 6, 7, 12, 13, 14, 15]
pmt = lambda i: inds[i]

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

    def __str__(self):
        return str(self.pair[0]) + "\n" + str(self.pair[1])

class Setup:
    def __init__(self, eta, q, mu, gPlus, chi, rhoInit):
        self.eta = eta
        self.q = q
        self.mu  = mu
        self.gPlus = gPlus
        self.tan = np.tan(chi * np.pi / 180.0) # chi is in degrees
        self.cos = np.cos(chi * np.pi / 180.0)
        self.rhoInit = rhoInit

# Returns a vacuum Hamiltonian in the units of omega_{vac}
def VacHam(setup):
    return -0.5 * setup.eta * flavourSigma3

# A function that states in the collective part of the total Hamiiltonian
def CollHam(setup, rho):
    SMPart   = np.trace(rho @ G)*G + su4Diag(su4Round(rho))
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

def LMatrix(setup):
    # A matrix of the linear map of L, acting on the entire
    # space of the \delta \rho matrices
    LM = np.zeros((32,32), dtype=np.complex128)
    for i in range(32):
        tmp = L(setup, Pair.singleEntry(i))
        for j in range(32):
            LM[j,i] = tmp[j]
    return LM

def LMatrixOffdiagonal(setup):
    #A matrix of the linear map of L, acting on the subspace
    # of the offdiagonal matrices \delta \rho
    LOff = np.zeros((16,16), dtype=np.complex128)
    for i in range(16):
        tmp = L(setup, Pair.singleEntry(i))
        for j in range(16):
            LOff[j,i] = tmp[j]
    return LOff

# Returns max Re(k/\omega_{vac})
def KMax(setup):
    return BiggestRealEigPart(LMatrixOffdiagonal(setup))

class Stability:
    """ And finally the class containing different stability diagrams """
    def __init__(self, dir):
        self.rhoInit = np.diag([0.5, 0.1, 0.3, 0.1])
        self.dir = dir

    @elapsed
    def MuQ(self, mulims, qlims, eta, gPlus, Nmu=100, Nq=100, filetitle="MuQ"):
        Grid = np.zeros((Nq,Nmu))
        mus = list(np.linspace(mulims[0], mulims[1], Nmu))
        qs  = list(np.linspace(qlims[0], qlims[1], Nq))
        
        for iq,q in enumerate(qs):
            print(f"Evaluating {q=:.2f} in the range ({qlims[0]},{qlims[1]}), {iq}/{Nq} points")
            for imu,mu in enumerate(mus):
                setup = Setup(eta=eta,q=q,mu=mu,gPlus=gPlus,chi=15.0,rhoInit=self.rhoInit)
                Grid[iq,imu] = KMax(setup)

        hrchy = "NH" if eta == 1.0 else "IH"

        PlotStability(mus, mulims, r"$\mu/\omega$",
                      qs, qlims, r"$q/\omega$",
                      Grid, f"{self.dir}/{filetitle}.eps",
                      hrchy)

    @elapsed
    def gPlusQ(self, gPluslims, qlims, eta, mu, NgPlus=100, Nq=100, filetitle="gPlusQ"):
        Grid = np.zeros((Nq,NgPlus))
        gPluses = list(np.linspace(gPluslims[0], gPluslims[1], NgPlus))
        qs  = list(np.linspace(qlims[0], qlims[1], Nq))
        
        for iq,q in enumerate(qs):
            print(f"Evaluating {q=:.2f} in the range ({qlims[0]},{qlims[1]}), {iq}/{Nq} points")
            for igPlus,gPlus in enumerate(gPluses):
                setup = Setup(eta=eta,q=q,mu=mu,gPlus=gPlus,chi=15.0,rhoInit=self.rhoInit)
                Grid[iq,igPlus] = KMax(setup)

        hrchy = "NH" if eta == 1.0 else "IH"

        PlotStability(gPluses, gPluslims, r"$g_{+}$",
                      qs, qlims, r"$q/\omega$",
                      Grid, f"{self.dir}/{filetitle}.eps",
                      hrchy)

# Draw the basic stability diagrams used in the article
def StabilityDiagrams(N, dir):
    stab = Stability(dir=dir)

    stab.MuQ(mulims=(0,50), qlims=(0,100), eta=1.0,  gPlus=1.0, Nmu=N, Nq=N, filetitle="MuQ_NH")
    stab.MuQ(mulims=(0,50), qlims=(0,100), eta=-1.0, gPlus=1.0, Nmu=N, Nq=N, filetitle="MuQ_IH")

    stab.gPlusQ(gPluslims=(0,1.0), qlims=(0,100), eta=1.0,  mu=20.0, NgPlus=N, Nq=N, filetitle="gPlusQ_NH")
    stab.gPlusQ(gPluslims=(0,1.0), qlims=(0,100), eta=-1.0, mu=20.0, NgPlus=N, Nq=N, filetitle="gPlusQ_IH")

# A short function that casts 16x16 matrices to text changing the order of the indices
# with some formatting
def MyFancyArrayToString(nparray, permutate=True):
    indices = [0, 1, 10, 11, 2, 3, 8, 9, 4, 5, 6, 7, 12, 13, 14, 15]
    if permutate:
        textarray = [[f"{nparray[i,j]:6.2f}" for j in indices] for i in indices]
    else:
        textarray = [[f"{nparray[i,j]:6.2f}" for j in range(16)] for i in range(16)]
    
    finalstring = ""

    # A simple function constructing a line with a separator between the blocks
    def line(arr):
        return " ".join(arr[:8]) + "  " + " ".join(arr[8:16]) + "\n"
    
    # Write the matrix by its blocks 8x8
    # Top blocks
    for i in range(8): finalstring += line(textarray[i])
    # Separator beetween the top and the bottom blocks
    finalstring += "\n"
    # Bottom blocks
    for i in range(8,16): finalstring += line(textarray[i])

    return finalstring

# The direct formula of the matrix of the linear map standing in the linear stability analysis
def LMatrixOffdiagonal_Direct(setup):
    # First, let's evaluate some additional constants
    delta = -setup.q*setup.tan
    alphaPlus  = delta + setup.eta/setup.cos
    alphaMinus = delta - setup.eta/setup.cos
    mub = setup.mu / setup.cos
    gPlus = setup.gPlus

    # Some functions of the probabilities
    se = setup.rhoInit[0,0]; sx = setup.rhoInit[1,1]; seb = setup.rhoInit[2,2]; sxb = setup.rhoInit[3,3]
    De = se - seb
    Dx = sx - sxb
    Offex = mub*gPlus*(se - sxb)
    Offxe = mub*gPlus*(sx - seb)
    Diag  = mub*(De + Dx)

    # Then, let's define the blocks 4x4
    L_A = np.array([[ alphaMinus+3.0*Diag, 0.0,                   Offex,              -Offex              ],
                    [ 0.0,                 alphaPlus+3.0*Diag,   -Offxe,               Offxe              ],
                    [-Offxe,               Offxe,                -alphaPlus-3.0*Diag, 0.0                 ],
                    [ Offex,              -Offex,                 0.0,                -alphaMinus-3.0*Diag]])

    L_B = np.array([[ alphaMinus-3.0*Diag, 0.0,                  -Offxe,               Offxe              ],
                    [ 0.0,                 alphaPlus-3.0*Diag,   Offex,               -Offex              ],
                    [ Offex,              -Offex,                -alphaPlus+3.0*Diag,  0.0                ],
                    [-Offxe,               Offxe,                 0.0,                -alphaMinus+3.0*Diag]])

    L_C = np.diag([ delta + 2.0*mub*(2.0*De + Dx),
                    delta + 2.0*mub*(De + 2.0*Dx),
                    delta - 2.0*mub*(2.0*De + Dx),
                    delta - 2.0*mub*(De + 2.0*Dx)])

    L_D = np.diag([-delta + 2.0*mub*(2.0*De + Dx),
                   -delta + 2.0*mub*(De + 2.0*Dx),
                   -delta - 2.0*mub*(2.0*De + Dx),
                   -delta - 2.0*mub*(De + 2.0*Dx)])

    # And the final matrix
    zeros = np.zeros((4,4), dtype=np.float64)

    return np.block([[L_A,   zeros, zeros, zeros],
                     [zeros, L_B,   zeros, zeros],
                     [zeros, zeros, L_C,   zeros],
                     [zeros, zeros, zeros, L_D]])

# Save particular matrix elements of the operator used to conduct the stability analysis
def SaveL(eta, q, mu, gPlus, fileout):
    setup = Setup(eta=eta, q=q, mu=mu, gPlus=gPlus, chi=15.0, rhoInit=np.diag([0.5, 0.1, 0.3, 0.1]))

    # Cook this matrix
    Matrix = np.real(LMatrixOffdiagonal(setup)*1.0j)
    # Update: now we have a direct formula!
    Direct_Matrix = LMatrixOffdiagonal_Direct(setup)
    
    # Save to the text file
    with open(fileout, "w") as f:
        f.write("Numerically evaluated matrix elements:\n")
        f.write(MyFancyArrayToString(Matrix))
        f.write("\nAnd via the direct formula:\n")
        f.write(MyFancyArrayToString(Direct_Matrix, permutate=False))

# A function that writes parts of the matrix of the linear map L
def PartsL(filename):
    # Zeroth part: the derivative part
    DerivativeMatrix = np.diag([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0])

    def CommutatorMap(matrix):
        result = np.empty((16, 16), dtype=np.float64)
        for i in range(16):
            b = Pair.singleEntry(i)
            tmp = Pair(com(matrix, b.left), com(matrix, b.right))
            for j in range(16):
                result[j,i] = tmp[j]
        return result

    # First part: the vacuum matrix
    VacuumMatrix = CommutatorMap(flavourSigma3)

    # Second part: first collective terms
    # IIIA term:
    IIIAMatrix = CommutatorMap(G)
    # IIIB term:
    IIIBMatrix = CommutatorMap(np.diag([1.0, 0.0, -1.0, 0.0]))
    # IIIC term:
    IIICMatrix = CommutatorMap(np.diag([0.0, 1.0, 0.0, -1.0]))

    # The third part: second collective terms:
    IVA = [np.empty((16,16), dtype=np.float64) for i in range(4)]
    IVB = [np.empty((16,16), dtype=np.float64) for i in range(4)]

    purefloavours = [np.diag([1.0, 0.0, 0.0, 0.0]),
                     np.diag([0.0, 1.0, 0.0, 0.0]),
                     np.diag([0.0, 0.0, 1.0, 0.0]),
                     np.diag([0.0, 0.0, 0.0, 1.0])]

    for i in range(16):
        b = Pair.singleEntry(i)
        # For each flavour
        for f in range(4):
            # A section
            tmpA = Pair(com(np.trace(b.right@G)*G + su4Diag(su4Round(b.right)), purefloavours[f]),
                        com(np.trace(b.left@G)*G + su4Diag(su4Round(b.left)), purefloavours[f]))
            # B section
            tmpB = Pair(com(su4Offdiag(su4Round(b.right)).T,purefloavours[f]),
                        com(su4Offdiag(su4Round(b.left)).T, purefloavours[f]))
            
            # Save the components
            for j in range(16):
                IVA[f][j,i] = tmpA[j]
                IVB[f][j,i] = tmpB[j]

    # Finally, let's dump all of these matrices
    with open(filename, "w") as f:
        f.write("Derivative Part:\n")
        f.write(MyFancyArrayToString(DerivativeMatrix))
        f.write("\n")
        f.write("Vacuum Part:\n")
        f.write(MyFancyArrayToString(VacuumMatrix))

        f.write("\n")
        f.write("IIIA:\n")
        f.write(MyFancyArrayToString(IIIAMatrix))
        f.write("\n")
        f.write("IIIB:\n")
        f.write(MyFancyArrayToString(IIIBMatrix))
        f.write("\n")
        f.write("IIIC:\n")
        f.write(MyFancyArrayToString(IIICMatrix))

        for fl in range(4):
            f.write("\n")
            f.write(f"IVA[{fl}]:\n")
            f.write(MyFancyArrayToString(IVA[fl]))
            f.write(f"IVB[{fl}]:\n")
            f.write(MyFancyArrayToString(IVB[fl]))

if __name__ == "__main__":
    #StabilityDiagrams(N=500, dir="build/StabilityPlots")
    SaveL(eta=1.0, q=20.0, mu=25.0, gPlus=1.0, fileout="build/StabilityPlots/LOffdiagonal.txt")
    #PartsL("build/StabilityPlots/PartsL.txt")