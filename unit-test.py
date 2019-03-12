#!/usr/bin/env python3

import os
import unittest
import math

from linepy.matrix import *
from molecule import *

class TestPBC(unittest.TestCase):
    def setUp(self):
        self.Bound2 = PeriodicBox()
        self.Bound2.setWalls(Vector(0,0), Vector(1,2))
        self.Bound3 = PeriodicBox()
        self.Bound3.setWalls(Vector3D(0,0,0), Vector3D(1,2,3))
        self.Bound10 = PeriodicBox()
        self.Bound10.setWalls(Vector(*([0,]*10)), Vector(*(range(1, 11))))
        self.Bound3p = PeriodicBox()
        self.Bound3p.setWalls(Vector3D(-0.4, -0.5, -0.5), Vector3D(0.6,0.5,0.5))

    def test_add_walls(self):
        self.assertEqual(self.Bound2.Dims, 2)
        self.assertEqual(self.Bound2.spanOfDim(0), 1)
        self.assertEqual(self.Bound2.spanOfDim(1), 2)

        self.assertEqual(self.Bound3.Dims, 3)
        self.assertEqual(self.Bound3.spanOfDim(0), 1)
        self.assertEqual(self.Bound3.spanOfDim(1), 2)
        self.assertEqual(self.Bound3.spanOfDim(2), 3)
        self.assertEqual(self.Bound3p.spanOfDim(0), 1)
        self.assertEqual(self.Bound3p.spanOfDim(1), 1)
        self.assertEqual(self.Bound3p.spanOfDim(2), 1)

        self.assertEqual(self.Bound10.Dims, 10)
        self.assertEqual(self.Bound10.spanOfDim(0), 1)
        self.assertEqual(self.Bound10.spanOfDim(9), 10)

    def test_offset1D(self):
        self.assertEqual(self.Bound3.offset1D(0, 0, 0), 0)
        self.assertEqual(self.Bound3.offset1D(0, 0, 1), 0)
        self.assertEqual(self.Bound3.offset1D(0, 0, 0.5), 0.5)
        self.assertAlmostEqual(self.Bound3.offset1D(0, 0, 0.50001), 0.49999)
        self.assertAlmostEqual(self.Bound3.offset1D(0, 0, 0.9), 0.1)
        self.assertEqual(self.Bound3p.offset1D(0, -0.4, 0.6), 0)
        self.assertAlmostEqual(self.Bound3p.offset1D(0, -0.41, 0.6), 0.01)
        self.assertEqual(self.Bound3p.offset1D(0, 0.1, 0.6), 0.5)
        self.assertAlmostEqual(self.Bound3p.offset1D(0, 0.10001, 0.6), 0.49999)

    def test_dist(self):
        self.assertEqual(self.Bound2.dist(Vector(0, 0), Vector(0, 1)), 1)
        self.assertEqual(self.Bound2.dist(Vector(0, 0), Vector(0, 2)), 0)
        self.assertEqual(self.Bound2.dist(Vector(0, 0), Vector(0.5, 2)), 0.5)
        self.assertEqual(self.Bound2.dist(Vector(0, 0), Vector(0.5, 2.5)),
                         0.5 * math.sqrt(2))

class TestMolecule(unittest.TestCase):
    def test_atom_import(self):
        Lines = ["ATOM   6128  O4  HPO    18      -1.362  18.749  19.735  1.00  0.00      HPO1",
                 "ATOM     33  H17 HEX X   1     -30.828 -28.792   4.852  1.00  0.00      HEX1          "]

        a1 = Atom.loadFromPDBLine(Lines[0])
        a2 = Atom.loadFromPDBLine(Lines[1])

        self.assertEqual(a1.x, float("-1.362"))
        self.assertEqual(a2.x, float("-30.828"))
        self.assertEqual(a2.z, float("4.852"))
        self.assertEqual(a1.ID, "O4")
        self.assertEqual(a2.ID, "H17")
        self.assertEqual(a1.Residue, "HPO")
        self.assertEqual(a1.ResidueNum, 18)
        self.assertEqual(a2.Residue, "HEX")
        self.assertEqual(a1.Occupancy, 1.0)
        self.assertEqual(a2.Temp, 0.0)

    def test_molecule_import(self):
        Lines = ["ATOM   6128  O4  HPO    18      -1       2       3      1.00  0.00      HPO1",
                 "ATOM     33  H17 HEX X   1       0       8       4      1.00  0.00      HEX1          ",
                 "END"]

        m = Molecule.loadFromPDB(Lines)
        self.assertEqual(m.Size, 2)
        self.assertEqual(m[0].x, -1.0)
        self.assertEqual(m[1].z, 4.0)

    def test_molecule_export(self):
        Lines = ["ATOM   6128  O4  HPO    18      -1       2       3      1.00  0.00      HPO1",
                 "ATOM     33  H17 HEX X   1       0       8       4      1.00  0.00      HEX1          ",
                 "END"]

        m = Molecule.loadFromPDB(Lines)
        self.assertEqual(m.Size, 2)

        with open("test.lalala", 'w') as File:
            m.saveAsPDB(File)

        with open("test.lalala", 'r') as File:
            m1 = Molecule.loadFromPDB(File)

        self.assertEqual(m[0], m1[0])
        self.assertEqual(m[1], m1[1])

        os.remove("test.lalala")

    def test_molecule_math(self):
        Lines = ["ATOM   6128  O4  HPO    18      -1       2       3      1.00  0.00      HPO1",
                 "ATOM     33  H17 HEX X   1       0       8       4      1.00  0.00      HEX1          ",
                 "END"]

        m = Molecule.loadFromPDB(Lines)
        self.assertEqual(m.GeoCenter, Vector3D(-0.5, 5.0, 3.5))

    def test_molecule_translate(self):
        Lines = ["ATOM   6128  O4  HPO    18      -1       2       3      1.00  0.00      HPO1",
                 "ATOM     33  H17 HEX X   1       0       8       4      1.00  0.00      HEX1          ",
                 "END"]

        m = Molecule.loadFromPDB(Lines)
        m.translate(Vector3D(2, 3, 4))
        self.assertEqual(m[0].x, 1)
        self.assertEqual(m[0].y, 5)
        self.assertEqual(m[0].z, 7)
        self.assertEqual(m[1].x, 2)
        self.assertEqual(m[1].y, 11)
        self.assertEqual(m[1].z, 8)

if __name__ == '__main__':
    unittest.main()
