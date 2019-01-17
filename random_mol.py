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

import sys
import math
import random
import copy
import logging

import linepy.matrix as matrix
import molecule

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

Logger = getLogger(level=logging.INFO)

def randRotation():             # The most naive way...
    """Return a random rotation matrix in 3D."""
    Theta = random.uniform(0, 2*math.pi)
    MatX = matrix.Matrix(3, 3)
    Cos = math.cos(Theta)
    Sin = math.sin(Theta)
    MatX[0, 0] = 1.0
    MatX[1, 1] = Cos
    MatX[1, 2] = -Sin
    MatX[2, 1] = Sin
    MatX[2, 2] = Cos

    Theta = random.uniform(0, 2*math.pi)
    MatY = matrix.Matrix(3, 3)
    Cos = math.cos(Theta)
    Sin = math.sin(Theta)
    MatY[1, 1] = 1.0
    MatY[0, 0] = Cos
    MatY[0, 2] = Sin
    MatY[2, 0] = -Sin
    MatY[2, 2] = Cos

    Theta = random.uniform(0, 2*math.pi)
    MatZ = matrix.Matrix(3, 3)
    Cos = math.cos(Theta)
    Sin = math.sin(Theta)
    MatZ[2, 2] = 1.0
    MatZ[0, 0] = Cos
    MatZ[0, 1] = -Sin
    MatZ[1, 0] = Sin
    MatZ[1, 1] = Cos

    return MatX * MatY * MatZ

def randPlaceAround(small, count, center, box_size, rand_rot=False):
    """Randomly place `count` number of `small` molecules around molecule
    `center` inside a cubic box with `box_size`, so that all instances
    of `small` will not be inside `center`.

    `Center` can be `None`, which means no center molecule will be
    considered.
    """

    if center is None:
        Radius = 0.0
    else:
        Radius = center.radius()
    Box = small.boundingBox() * 0.5
    Logger.debug("Bounding box of small molecule: " + str(Box))
    SmallZeroed = copy.deepcopy(small)
    SmallZeroed.translate(-(small.GeoCenter))
    # SmallZeroed is now center at (0,0,0).
    if center is not None:
        Logger.info("Center molecule is at ({}) with radius {:.2f}.".format(
            center.GeoCenter, Radius))

    if center is None:
        CenterCenter = matrix.Vector3D(0, 0, 0) # This doesn't actually do anything.
    else:
        CenterCenter = center.GeoCenter
    Result = molecule.Molecule()
    iRaw = 0
    i = 0
    while i < count:
        TransVec = matrix.Vector3D(
            random.uniform(-box_size / 2.0 + Box.x, box_size / 2.0 - Box.x),
            random.uniform(-box_size / 2.0 + Box.y, box_size / 2.0 - Box.y),
            random.uniform(-box_size / 2.0 + Box.z, box_size / 2.0 - Box.z))
        if TransVec.distanceFrom(CenterCenter) >= Radius:
            Logger.debug("Translation vector: " + str(TransVec))
            NewMol = SmallZeroed.duplicate()
            if rand_rot:
                RotMat = randRotation()
                NewMol.rotate(RotMat)
            for a in NewMol:
                a.ResidueNum = i+1
                a.ResidueID = "HEX1"
            NewMol.translate(TransVec)
            Logger.debug("Translated center: " + str(NewMol.GeoCenter))
            for a in NewMol:
                Logger.debug(str(a))
            Result.addMolecule(NewMol)
            i += 1
        iRaw += 1

    Logger.info("Dropped {:.1f}% random samples.".format((iRaw - i) / iRaw * 100.0))
    return Result

def main():
    import argparse

    Parser = argparse.ArgumentParser(
        description="Put some copies of MOL1 randomly around MOL2, "
        "inside some cubic boundary box, without intruding the space occupied "
        "by MOL2.")
    Parser.add_argument('Small', metavar='MOL1', type=str,
                        help='The PDB file name of molecule to randomly put.')
    Parser.add_argument('Center', metavar='MOL2', type=str, nargs='?',
                        default="",
                        help="The PDB file name of the molecule around which "
                        "copies of MOL1 are put.  If this is not given, "
                        "put MOL1 all over the place.")
    Parser.add_argument('-n', "--count", dest='Count', type=int, default=497,
                        help='Number of copies.  Default: %(default)s')
    Parser.add_argument('-o', "--output", dest='Dest', type=str,
                        default="output.pdb",
                        help="Result PDB will be written to this file.  "
                        "If it's '-', write to stdandard output.  "
                        "Default: %(default)s")
    Parser.add_argument('-s', "--box-size", dest='Size', type=float, default=90,
                        help='Size of the boundary box.  Default: %(default)s')

    Parser.add_argument('-r', "--random-rotate", dest="Rotate", default=False,
                        action="store_true",
                        help="Apply random rotation to small mulecules.")

    Args = Parser.parse_args()

    random.seed()
    Logger.info("Loading molecules...")
    with open(Args.Small, 'r') as File:
        Small = molecule.Molecule.loadFromPDB(File)
    if Args.Center:
        with open(Args.Center, 'r') as File:
            Center = molecule.Molecule.loadFromPDB(File)
    else:
        Center = None

    Logger.info("Placing molecules...")
    Result = randPlaceAround(Small, Args.Count, Center, Args.Size, Args.Rotate)

    if Args.Dest == '-':
        Result.saveAsPDB(sys.stdout)
    else:
        with open(Args.Dest, 'w') as File:
            Result.saveAsPDB(File)

    return 0

if __name__ == "__main__":
    sys.exit(main())
