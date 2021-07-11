#!/usr/bin/env python3

import unittest, sys
import numpy as np

sys.path.append("..")
from Lambdas import su4Diag, su4Offdiag, su4Round, BiggestRealEigPart

class TestLambdas(unittest.TestCase):

    def test_su4Diag(self):
        m = np.array([[1, 2, 3, 4],[5, 6, 7, 8],[9, 10, 11, 12],[13, 14, 15, 16]])
        d = np.array([[1, 2, 0, 0],[5, 6, 0, 0],[0, 0,  11, 12],[0,  0,  15, 16]])
        r = su4Diag(m)

        for i in range(4):
            for j in range(4):
                self.assertEqual(d[i,j], r[i,j])

    def test_su4Offdiag(self):
        m = np.array([[1, 2, 3, 4],[5, 6, 7, 8],[9, 10, 11, 12],[13, 14, 15, 16]])
        d = np.array([[0, 0, 3, 4],[0, 0, 7, 8],[9, 10, 0,  0 ],[13, 14, 0,  0 ]])
        r = su4Offdiag(m)

        for i in range(4):
            for j in range(4):
                self.assertEqual(d[i,j], r[i,j])
        
    def test_su4Round(self):
        m = np.array([[1, 2, 3, 4],[5, 6, 7, 8],[9, 10, 11, 12],[13, 14, 15, 16]])
        C = np.array([[0, 0, 1, 0],[0, 0, 0, 1],[1, 0, 0, 0],[0, 1, 0, 0]])
        d = m - ((C @ m) @ C).T
        r = su4Round(m)

        for i in range(4):
            for j in range(4):
                self.assertEqual(d[i,j], r[i,j])

    def test_BiggestRealEigPart(self):
        m = np.diag([1.0+1.0j, 2.0+4.0j, -3.0+5j])
        self.assertEqual(BiggestRealEigPart(m), 2.0)

if __name__ == "__main__":
    unittest.main()