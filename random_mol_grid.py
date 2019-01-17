#!/usr/bin/env python3

import sys, os

import linepy.matrix as Matrix
import molecule
import random_mol
import random

def placeByGrid(ns, big_box_size, mol):
    Result = molecule.Molecule()
    for ix in range(ns):
        for iy in range(ns):
            for iz in range(ns):
                Pos =  Matrix.Vector3D(
                    big_box_size / ns * ix - big_box_size * 0.5,
                    big_box_size / ns * iy - big_box_size * 0.5,
                    big_box_size / ns * iz - big_box_size * 0.5)

                Single = random_mol.randPlaceAround(mol, 1, None, big_box_size / ns, True)
                Single.translate(Pos)
                Result.addMolecule(Single)

    return Result

if __name__ == "__main__":
    random.seed()
    # Test
    with open(sys.argv[1], 'r') as PdbFile:
        Mol = molecule.Molecule.loadFromPDB(PdbFile)
    Result = placeByGrid(2, 200, Mol)
    Result.saveAsPDB(sys.stdout)
