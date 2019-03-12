# Read PDB file from stdin, add the residue number by 1 for each atom, and write
# the result to stdout as a PDB file.

import sys, os
import molecule

Mol = molecule.Molecule.loadFromPDB(sys.stdin)
for Atom in Mol:
    Atom.ResidueNum += 1

Mol.saveAsPDB(sys.stdout)
