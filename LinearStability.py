#!/usr/bin/env python3

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker, rcParams, cm

import matplotlib as mpl
import numpy as np

import sys

# Colours and other constants
colour_e = "#324851"
colour_x = "#86AC41"
colour_ae = "#34675C"
colour_ax = "#7DA3A1"
colours  = colour_e, colour_x, colour_ae, colour_ax

viridis  = cm.get_cmap('viridis')
twilight = cm.get_cmap('twilight')
bone     = cm.get_cmap('bone')

MainFontSize = 12

# Some customisation of matplotlib
ticker.rcParams['xtick.direction'] = 'in'
ticker.rcParams['ytick.direction'] = 'in'
ticker.rcParams['xtick.labelsize'] = MainFontSize
ticker.rcParams['ytick.labelsize'] = MainFontSize

mpl.rcParams['font.family'] = ['serif']
mpl.rcParams['font.serif']  = ['Times New Roman'] # Computer Modern Roman
mpl.rcParams['legend.handlelength'] = 2


import numpy as np
np.set_printoptions(precision=2)

import Constants as ph
from Lambdas import flavourSigma1, flavourSigma3, G

# Index maps
directOffdiagonal = {0:(0,3), 1:(1,2), 2:(2,1), 3:(3,0), 4:(0,2), 5:(1,3), 6:(2,0), 7:(3,1)}
inverseOffdiagonal = [[None, None, 4, 0],
                      [None, None, 1, 5],
                      [6, 2, None, None],
                      [3, 7, None, None]]

# Offdiagonal indecies
offdiagonalInds = directOffdiagonal.values()
# Diagonal indecies
diagonalInds = [(0,0), (0,1), (1,0), (1,1), (2,2), (2,3), (3,2), (3,3)]




""" Some operations from su4 modules """



# Returns the diagonal blocks
def su4Diag(m):
	M = np.array(m)
	return np.matrix([[M[0][0], M[0][1], 0.0,     0.0],
					  [M[1][0], M[1][1], 0.0,     0.0],
					  [0.0,     0.0,     M[2][2], M[2][3]],
					  [0.0,     0.0,     M[3][2], M[3][3]]])

# Returns the offdiagonal blocks
def su4Offdiag(m):
	M = np.array(m)
	return np.matrix([[0.0,     0.0,     M[0][2], M[0][3]],
					  [0.0,     0.0,     M[1][2], M[1][3]],
					  [M[2][0], M[2][1], 0.0,     0.0],
					  [M[3][0], M[3][1], 0.0,     0.0]])

# Evaluate the round matrix
def su4Round(m):
	M = np.array(m)
	# The matrix that corresponds to (rho^C)^T
	M_C_T = np.matrix([[M[2][2], M[3][2], M[0][2], M[1][2]],
		               [M[2][3], M[3][3], M[0][3], M[1][3]],
		               [M[2][0], M[3][0], M[0][0], M[1][0]],
		               [M[2][1], M[3][1], M[0][1], M[1][1]]])

	return M - M_C_T

# Simple commutator
def com(a,b):
    return a*b - b*a

# Find the biggest real part among the eigenvalues of a matrix
def findBiggest(m):
    # Eigenvalues
    eigs = np.linalg.eig(m)[0]
    return np.max([np.real(e) for e in eigs])


""" Then Hamiltonians' defintiions """


# Vacuum Hamiltonian --- in the linear stability anaysis theta = 0!!!
def VacuumHamiltonian(c):
    return (-c.vacOmega*c.eta/2.0) * flavourSigma3

# MSW part of the flavour effective Hamiltonian
def MSWHamiltonian(c):
    return np.diag([c.V_e-0.5*c.V_n, -0.5*c.V_n, -c.V_e+0.5*c.V_n, 0.5*c.V_n])

# Standard Model collective term;
def VAHamiltonian(c, rhoOpp, relMu):
    tmp      = np.trace(rhoOpp*G)*G + su4Diag(su4Round(rhoOpp))
    return tmp * (relMu * c.vacOmega)

# The collective term that describes coupling through the scalar and 
# pseudoscalar fields
def NSSIHamiltonian(c, rhoOpp, gPlus, relMu):
    diagRound    = su4Diag(su4Round(rhoOpp))
    offdiagRound = su4Offdiag(su4Round(rhoOpp))
    tmp          = c.gMinus*diagRound.getT() + gPlus*offdiagRound.getT()
    return tmp * (relMu * c.vacOmega)

# And the total collective Hamiltonian
def CollectiveHam(c, rho, gPlus, relMu):
    return (VAHamiltonian(c,rho,relMu) + NSSIHamiltonian(c,rho,gPlus,relMu))



""" The function that evaluates 16x16 matrix of the coefficients """



def AssembleOne(dzeta, i, j):
    # Assemble delta with only one "1" at dzeta,i,j
    delta = [np.zeros((4,4)), np.zeros((4,4))]
    delta[dzeta][i][j] = 1.0
    return np.matrix(delta[0]), np.matrix(delta[1])


def Derivatives(c, q, deltaLeft, deltaRight, gPlus, relMu): # delta means \delta \rho
    # x-derivative contribution
    derLeft  = (1.0j*q*c.tanChi)*deltaLeft
    derRight = (-1.0j*q*c.tanChi)*deltaRight  # Note: different signs are caused by different directions

    # Noncollective contribution
    hVac = VacuumHamiltonian(c) + MSWHamiltonian(c)
    derLeft  += (-1.0j/c.cosChi)*com(hVac, deltaLeft)
    derRight += (-1.0j/c.cosChi)*com(hVac, deltaRight)

    # Collective contribution
    derLeft  += (-1.0j/c.cosChi)*com(CollectiveHam(c,c.InitialRho,gPlus,relMu), deltaLeft)
    derRight += (-1.0j/c.cosChi)*com(CollectiveHam(c,c.InitialRho,gPlus,relMu), deltaRight)

    derLeft  += (-1.0j/c.cosChi)*com(CollectiveHam(c,deltaRight,gPlus,relMu), c.InitialRho)
    derRight += (-1.0j/c.cosChi)*com(CollectiveHam(c,deltaLeft,gPlus,relMu), c.InitialRho)

    return derLeft, derRight

# Assmeble the matrix
def Coeffs(c, q, gPlus,relMu):
    L = np.zeros((16,16), dtype=np.complex)

    for dzeta in (0,1): # 0 stands for left, 1 stands for right
        for i,j in offdiagonalInds:
            I = inverseOffdiagonal[i][j] + 8*dzeta
            # Assemble the delta matrices with only one nonzero component
            deltaLeft, deltaRight = AssembleOne(dzeta, i, j)
            # Evaluate the deriavtives
            ders = Derivatives(c, q, deltaLeft, deltaRight,gPlus,relMu)

            for psi in (0,1):
                for n,m in offdiagonalInds:
                    J = inverseOffdiagonal[n][m] + 8*psi
                    L[J][I] = np.array(ders[psi])[n][m]

    return L




""" Plot the rates """



def PlotRates(Rates, xs, ys, filePlot):
    xGrid, yGrid = np.meshgrid(xs, ys)

    fig = Figure(figsize=(8, 6))
    FigureCanvas(fig)

    axs  = fig.add_subplot(111)
    plot = axs.pcolor(xGrid, yGrid, Rates, cmap=viridis, rasterized=True)

    axs.set_xlabel(r"$g_{+}$")
    axs.set_ylabel(r"$q / \omega$")

    fig.colorbar(plot)
    fig.savefig(filePlot, fmt="eps", bbox_inches='tight', dpi=192)



if __name__ == "__main__":
    # Read the Parameters.json
    c = ph.Constants()   

    # Main calculational part

    print(c.V_Nu * (1.0 - c.cosOmega) / c.vacOmega)

    N = 201 # Number of points in each axis
    M = 201
    # Grids of axes
    #relMus  = np.linspace(0.0, 50.0, N)
    gPluses = np.linspace(0.0, 1.0, N)
    ks = np.linspace(0.0, 100, M)

    # The mesh of rates
    Rates = np.zeros((M,N))

    # Calculate the mesh of rates
    for i, gPlus in enumerate(gPluses):
        print(f"Evaluating: i={i} out of {N}")
        for j, k in enumerate(ks):
            q = k * c.vacOmega
            Rates[j][i] = findBiggest(Coeffs(c, q, gPlus, 20.0)) / c.vacOmega
            #Rates[j][i] = i

    # Plot the rates
    PlotRates(Rates, gPluses, ks, "./TalkPlots/LinearStability_relMu_20_IH.eps")

    # Assemble the delta matrices with only one nonzero component
    #deltaLeft, deltaRight = AssembleOne(0, 0, 3)
    #print(deltaLeft)
    #print(deltaRight)

    #print(NSSIHamiltonian(c,deltaLeft))

    # Evaluate the deriavtives
    #derLeft, derRight = Derivatives(c, 0.0, deltaLeft, deltaRight)

    """

    # Little check -- equations on diagonal and offdiagonal components are separate
    print("Diagonal coefficients \n\n")
    for dzeta in (0,1): # 0 stands for left, 1 stands for right
        for i,j in diagonalInds:
            # Assemble the delta matrices with only one nonzero component
            deltaLeft, deltaRight = AssembleOne(dzeta, i, j)
            # Evaluate the deriavtives
            derLeft, derRight = Derivatives(c, 0.0, deltaLeft, deltaRight)

            print(f"dzeta={dzeta}, i={i}, j={j}")
            print(">derLeft:")
            print(derLeft)
            print(">derRight:")
            print(derRight, "\n")

    print("Offdiagonal coefficients \n\n")
    for dzeta in (0,1): # 0 stands for left, 1 stands for right
        for i,j in offdiagonalInds:
            # Assemble the delta matrices with only one nonzero component
            deltaLeft, deltaRight = AssembleOne(dzeta, i, j)
            # Evaluate the deriavtives
            derLeft, derRight = Derivatives(c, 0.0, deltaLeft, deltaRight)

            print(f"dzeta={dzeta}, i={i}, j={j}")
            print(">derLeft:")
            print(derLeft)
            print(">derRight")
            print(derRight, "\n")

    """