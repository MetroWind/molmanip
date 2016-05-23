#!/usr/bin/env python3

# MolManip, PDB processing library in Python.
# Copyright (C) 2016 Mingyang Sun (aka. MetroWind)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# This module is pdoc-ready.  Generate document by
# "PYTHONPATH=.:$PYTHONPATH pdoc --html --html-no-source --overwrite THIS_FILE".

"""Utility API mainly for manipulating atoms and molecules in PDB
files.  Use the `Atom` and `Molecule` class.
"""

from __future__ import print_function
import logging

import linepy.matrix as matrix

def getLogger(name="Main", level=logging.DEBUG):
    import sys
    Logger = logging.getLogger(name)
    Logger.setLevel(level)
    Handler = logging.StreamHandler(sys.stderr)
    Format = logging.Formatter("[%(levelname)s %(asctime)s] %(message)s",
                                  "%Y-%m-%d %H:%M:%S")
    Handler.setFormatter(Format)
    Logger.addHandler(Handler)
    return Logger

Logger = getLogger("Molecule", logging.INFO)

class Atom(object):
    """Information of an atom in a PDB file."""
    Warned = {"Residue": False,}
    def __init__(self):
        self.Loc = matrix.Vector3D()
        """Location of the atom."""
        self.ID = ""
        """ID of the atom, 3rd column in PDB files.  If name of the element is
        only 1 character, this should have 1 leading space.
        """
        self.Residue = ""
        """Name of residue."""
        self.ResidueNum = 0
        """Residue number."""
        self.ResidueID = ""
        """This field is not in the PDB specification, but charmm uses it."""

    def __eq__(self, rhs):
        return self.Loc == rhs.Loc and self.ID == rhs.ID and \
            self.Residue == rhs.Residue and \
            self.ResidueNum == rhs.ResidueNum and \
            self.ResidueID == rhs.ResidueID

    def __str__(self):
        return "{} @({})".format(self.ID, self.Loc)

    @classmethod
    def loadFromPDBLine(cls, line):
        """Load an atom from one line in the PDB file."""
        NewStuff = cls()
        # See specification
        # http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
        x, y, z = [float(x) for x in line[31:54].split()]
        NewStuff.Loc = matrix.Vector3D(x, y, z)
        NewStuff.ID = line[12:16].rstrip() # Could have one leading space
        if line[20] == ' ':
            NewStuff.Residue = line[17:20]
        else:
            NewStuff.Residue = line[17:21]
            if not cls.Warned["Residue"]:
                Logger.warn("Reading from non-standard PDB.")
                cls.Warned["Residue"] = True
        NewStuff.ResidueNum = int(line[22:26])
        NewStuff.Occupancy = float(line[54:60])
        NewStuff.Temp = float(line[60:66])
        return NewStuff

    def translate(self, vec):
        """Translate the atom by a vector `vec`."""
        self.Loc += vec

    def rotate(self, rot_mat):
        """Rotate the atom around the origin with rotation matrix `rot_mat`."""
        NewLoc = rot_mat * self.Loc
        self.Loc = matrix.Vector3D(NewLoc[0], NewLoc[1], NewLoc[2])

    @property
    def x(self):
        """The x corrdinate of the atom."""
        return self.Loc.x

    @x.setter
    def x(self, value):
        self.Loc.x = value

    @property
    def y(self):
        """The y corrdinate of the atom."""
        return self.Loc.y

    @y.setter
    def y(self, value):
        self.Loc.y = value

    @property
    def z(self):
        """The z corrdinate of the atom."""
        return self.Loc.z

    @z.setter
    def z(self, value):
        self.Loc.z = value

class Molecule(object):
    def __init__(self):
        self._Atoms = []
        self._GeoCenter = None

    def addAtom(self, a):
        """Add an atom `a` to the molecule."""
        self._Atoms.append(a)
        self._GeoCenter = None
        return self

    def duplicate(self):
        """Return a deep copy of self."""
        import copy
        return copy.deepcopy(self)

    def addMolecule(self, m):
        """Add all atoms in `m` to this molecule."""
        for a in m:
            self.addAtom(a)

    @property
    def Size(self):
        """Number of atoms in the molecule."""
        return len(self._Atoms)

    def __len__(self):
        return self.Size

    @classmethod
    def loadFromPDB(cls, f):
        """Load all atoms from file object `f` as a whole molecule."""
        NewMol = cls()
        for Line in f:
            if Line.startswith("ATOM"):
                NewMol.addAtom(Atom.loadFromPDBLine(Line))
        return NewMol

    def saveAsPDB(self, f):
        """Save the current molecule to file object `f` in PDB format."""
        Warned = {"resSeq": False, "Residue": False}
        for i in range(self.Size):
            a = self[i]
            ResNumStr = str(a.ResidueNum)
            if len(ResNumStr) > 4:
                ResNum = ResNumStr[:4]
                ResNumOverflow = ResNumStr[4:]
                if not Warned["resSeq"]:
                    Logger.warn("resSeq too long, generating non-standard PDB.")
                    Warned["resSeq"] = True
            else:
                ResNum = str(a.ResidueNum)
                ResNumOverflow = ""
            if len(a.Residue) > 3 and not Warned["Residue"]:
                Logger.warn("Residue name too long, generating non-standard PDB.")
                Warned["Residue"] = True
            f.write("ATOM  {: >5d} {: <4s} {: <4s} {: >4s}{: <4s}"
                    "{:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}"
                    "      {: <4s}\n".format(
                        i+1, a.ID, a.Residue, ResNum, ResNumOverflow,
                        a.x, a.y, a.z, a.Occupancy, a.Temp, a.ResidueID))
        f.write("END\n")

    def __iter__(self):
        return self._Atoms.__iter__()

    @property
    def GeoCenter(self):
        """The geometric center of the molecule as a `Vector3D` object.  This
        is cached, so using it repeatedly will not recalculate.  Any
        change on the molecule (except `translate()`) will empty the
        cache.
        """
        if self._GeoCenter is None:
            x = sum(a.x for a in self) / float(self.Size)
            y = sum(a.y for a in self) / float(self.Size)
            z = sum(a.z for a in self) / float(self.Size)

            self._GeoCenter = matrix.Vector3D(x, y, z)

        return self._GeoCenter

    def radius(self):
        """Treat the molecule as a sphere centered at its geometric center,
        with radius calculated by this function, so that all the atoms
        of the molecule are fully contained by this sphere.  Return
        the minimum radius of such a sphere.
        """
        return max(a.Loc.distanceFrom(self.GeoCenter) for a in self)

    def boundingBox(self):
        """Return a `Vector3D` object `v`, so that the molecule is fully
        contained in a minimum box centered at its geometric center,
        with length, width, height of `v.x/2`, `v.y/2`, and `v.z/2`.
        """
        x = max(a.x for a in self) - min(a.x for a in self)
        y = max(a.y for a in self) - min(a.y for a in self)
        z = max(a.z for a in self) - min(a.z for a in self)
        return matrix.Vector3D(x, y, z)

    def translate(self, vec):
        """Translate the molecule in space by some 3D vector."""
        for Atom in self:
            Atom.translate(vec)
        if self._GeoCenter is not None:
            self._GeoCenter += vec
        return self

    def rotate(self, rot_mat):
        """Rotation the molecule with rotation matrix `rot_mat`."""
        for Atom in self:
            Atom.rotate(rot_mat)
        return self

    def __getitem__(self, i):
        return self._Atoms[i]

    def __setitem__(self, i, x):
        self._Atoms[i] = x
        self._GeoCenter = None
