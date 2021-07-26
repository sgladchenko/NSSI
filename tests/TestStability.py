#!/usr/bin/env python3

from copy import deepcopy
import unittest, sys
import numpy as np

sys.path.append("..")
from Stability import Pair, indexMapper, maskLeft, maskRight, LMatrix, Setup
from Stability import LMatrixOffdiagonal, LMatrixOffdiagonal_Direct, MyFancyArrayToString, pmt

class TestStability(unittest.TestCase):

    def test_indexMapper(self):
        mleft  = np.zeros((4,4))
        mright = np.zeros((4,4))

        for key in indexMapper[0]:
            i,j = indexMapper[0][key]
            mleft[i,j]  = key
            mright[i,j] = key+8

        for key in indexMapper[1]:
            i,j = indexMapper[1][key]
            mleft[i,j]  = key+16
            mright[i,j] = key+24

        dleft  = maskLeft
        dright = maskRight

        for i in range(4):
            for j in range(4):
                self.assertEqual(dleft[i,j], mleft[i,j])
                self.assertEqual(dright[i,j], mright[i,j])

    def test_Pair(self):
        # Just basic construction of an instance from two matrices
        mleft  = np.array([[1, 2, 3, 4],[5, 6, 7, 8],[9, 10, 11, 12],[13, 14, 15, 16]])
        mright = np.array([[17, 18, 19, 20],[21, 22, 23, 24],[25, 26, 27, 28],[29, 30, 31, 32]])
        pair   = Pair(mleft, mright)

        for i in range(4):
            for j in range(4):
                self.assertEqual(pair.left[i,j], mleft[i,j])
                self.assertEqual(pair.right[i,j], mright[i,j])

        ml = deepcopy(mleft)
        mr = deepcopy(mright)
        mleft[0,0] = -2.0  # # Test of the deepcopy inside the Pair's constructor:
        mright[0,0] = -1.0 # what's gonna happen if I change an element in the matrix from which
        # I have previously intialized the Pair instance? I have to make sure that it does nothing
        # to the Pair instance and its 

        for i in range(4):
            for j in range(4):
                self.assertEqual(pair[maskLeft[i,j]],  ml[i,j])
                self.assertEqual(pair[maskRight[i,j]], mr[i,j])

        # And then finally let's test the setitem operator
        for i in range(4):
            for j in range(4):
                pair[maskLeft[i,j]]  = mr[i,j]
                pair[maskRight[i,j]] = ml[i,j]

        for i in range(4):
            for j in range(4):
                self.assertEqual(pair[maskLeft[i,j]],  mr[i,j])
                self.assertEqual(pair[maskRight[i,j]], ml[i,j])

    def test_singleEntry(self):
        # Test of the alternative constructor that sets only one matrix element
        # in the pair of matrices to 1.
        for i in range(4):
            for j in range(4):
                m = np.zeros((4,4)); m[i,j] = 1.0

                # Set to the [i,j] in the left matrix
                pairLeft  = Pair.singleEntry(maskLeft[i,j])
                # Set to the [i,j] in the right matrix
                pairRight = Pair.singleEntry(maskRight[i,j])

                for a in range(4):
                    for b in range(4):
                        self.assertEqual(pairLeft[maskLeft[a,b]],   m[a,b])
                        self.assertEqual(pairRight[maskRight[a,b]], m[a,b])
                        self.assertEqual(pairLeft[maskRight[a,b]], 0.0)
                        self.assertEqual(pairRight[maskLeft[a,b]], 0.0)
    
    def test_LDirect(self):
        rhoInit = np.diag([0.5, 0.1, 0.3, 0.1])

        # Normal hierarchy  
        for q in np.linspace(0.0, 70.0, 10):
            for mu in np.linspace(0.0, 50.0, 10):
                for gPlus in np.linspace(0.0, 1.0, 10):
                    setup = Setup(eta=1.0,q=q,mu=mu,gPlus=gPlus,chi=15.0,rhoInit=rhoInit)
                    N = np.real(LMatrixOffdiagonal(setup)*1.0j)
                    D = LMatrixOffdiagonal_Direct(setup)
                    for i in range(16):
                        for j in range(16):
                            self.assertAlmostEqual(N[pmt(i),pmt(j)], D[i,j], 13, msg=f"{q=},{mu=},{gPlus=},{i=},{j=}")

        # Inverted hierarchy
        for q in np.linspace(0.0, 70.0, 10):
            for mu in np.linspace(0.0, 50.0, 10):
                for gPlus in np.linspace(0.0, 1.0, 10):
                    setup = Setup(eta=-1.0,q=q,mu=mu,gPlus=gPlus,chi=15.0,rhoInit=rhoInit)
                    N = np.real(LMatrixOffdiagonal(setup)*1.0j)
                    D = LMatrixOffdiagonal_Direct(setup)
                    for i in range(16):
                        for j in range(16):
                            self.assertAlmostEqual(N[pmt(i),pmt(j)], D[i,j], 13, msg=f"{q=},{mu=},{gPlus=},{i=},{j=}")

        # And for another choice of the initial unperturbed density matrcies
        rhoInit = np.diag([0.5, 0.0, 0.5, 0.0])

        # Normal hierarchy
        for q in np.linspace(0.0, 70.0, 10):
            for mu in np.linspace(0.0, 50.0, 10):
                for gPlus in np.linspace(0.0, 1.0, 10):
                    setup = Setup(eta=1.0,q=q,mu=mu,gPlus=gPlus,chi=15.0,rhoInit=rhoInit)
                    N = np.real(LMatrixOffdiagonal(setup)*1.0j)
                    D = LMatrixOffdiagonal_Direct(setup)
                    for i in range(16):
                        for j in range(16):
                            self.assertAlmostEqual(N[pmt(i),pmt(j)], D[i,j], 13, msg=f"{q=},{mu=},{gPlus=},{i=},{j=}")

        # Inverted hierarchy
        for q in np.linspace(0.0, 70.0, 10):
            for mu in np.linspace(0.0, 50.0, 10):
                for gPlus in np.linspace(0.0, 1.0, 10):
                    setup = Setup(eta=-1.0,q=q,mu=mu,gPlus=gPlus,chi=15.0,rhoInit=rhoInit)
                    N = np.real(LMatrixOffdiagonal(setup)*1.0j)
                    D = LMatrixOffdiagonal_Direct(setup)
                    for i in range(16):
                        for j in range(16):
                            self.assertAlmostEqual(N[pmt(i),pmt(j)], D[i,j], 13, msg=f"{q=},{mu=},{gPlus=},{i=},{j=}")

if __name__ == "__main__":
    unittest.main()