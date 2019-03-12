# MolManip

PDB processing library in Python.

## Why

Manipulating PDB files by hand is tedious, and nobody should ever need
to do so.

## Installation

The easiest way to install this library is to not do that, but instead
clone the git repository to somewhere in your project, by running

    git clone --recursive https://github.com/MetroWind/molmanip.git

You can also clone this as a git submodule.

## Usage

The main file of this library is `molecule.py`.  It models atoms and
molecules into classes, and provides functions to e.g. shift an atom
and rotate a molecule.  The `Molecule` class takes care of loading and
saving to PDB files.

For mole detail, generate documentation by running

    PYTHONPATH=.:$PYTHONPATH pdoc --html --html-no-source --overwrite molecule.py

given that [pdoc](http://pdoc.burntsushi.net/pdoc) is installed.

A simple example is provided in `resnum-add1.py`. It simply reads PDB
file from stdin, adds the residue number by 1 for each atom, and
writes the result to stdout as a PDB file.

Another example program using this library is given in
`random_mol.py`. This is a more involved example concerning all
aspects of the library, including linear algebra. It takes two
molecules, and randomly place copied of the 1st molecule around the
second in a cubic box. Random rotation can be applied to the copies of
1st molecule. For detailed usage, run `random_mol.py -h`.
