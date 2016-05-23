#!/usr/bin/env python3

import os
import unittest
from linepy.matrix import *
from molecule import *

class TestMolecule(unittest.TestCase):
    def test_atom_import(self):
        Lines = ["ATOM   6128  O4  HPO    18      -1.362  18.749  19.735  1.00  0.00      HPO1",
                 "ATOM     33  H17 HEX X   1     -30.828 -28.792   4.852  1.00  0.00      HEX1          "]

        a1 = Atom.loadFromPDBLine(Lines[0])
        a2 = Atom.loadFromPDBLine(Lines[1])

        self.assertEqual(a1.x, float("-1.362"))
        self.assertEqual(a2.x, float("-30.828"))
        self.assertEqual(a2.z, float("4.852"))
        self.assertEqual(a1.ID, " O4")
        self.assertEqual(a2.ID, " H17")
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
