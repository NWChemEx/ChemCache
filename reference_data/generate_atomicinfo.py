#!/usr/bin/env python3

"""This script is used to create the experimental data look up tables for the
atom class.

Original author: Ben Pritchard
Modified by: Zachery Crandall

In order to run, this script needs to know the location of the data directory
to read from and the `src` code directory to output the result into. One can 
also optionally provide a nondefault value for the electron mass to Dalton ratio.

For readability and convenience we use a few abbreviations throughout this
script:

- Z: the atomic number of an atom
- Sym: the atomic symbol of an atom (e.g. H for hydrogen, He for helium)

This script looks for the following file(s):

+---data_dir
|       ElementNames.txt
|       CIAAW-ISOTOPEMASSES.txt
|       CIAAW-MASSES.txt

This script creates the following file(s):

+---src_dir
|       load_elements.hpp
"""

import argparse
import os
import re

import helper_fxns as helpers


class AtomicData:
    def __init__(self):
        self.sym = ""
        self.name = ""
        self.Z = 0
        self.mass = 0.0
        self.isotopes = []       # List of isotope mass numbers
        self.isotope_masses = {}  # Isotope mass values, indexed by mass number

    def add_isotope(self, num: int, mass: float) -> None:
        """Add a new isotope mass to the list of isotopes for this element.

        :param num: Isotope mass number (Z + number of neutrons, N)
        :type num: int

        :param mass: Isotope mass value
        :type mass: float
        """

        # Set a new isotope mass for the given isotope
        self.isotope_masses[num] = mass

        # If the isotope does not already exist in the isotope atomic number
        # list, add it
        if not num in self.isotopes:
            self.isotopes.append(num)

    def __repr__(self):
        """Return a formatted text representation of the atomic data.

        :return: Formatted representation of the atomic data
        :rtype: str
        """

        rv = "{} {} {}[".format(self.Z, self.name, self.mass)

        for x in self.isotopes:
            rv += "{}: {}".format(x, self.isotope_masses[x])
            if x != self.isotopes[-1]:
                rv += ", "

        rv += "]"

        return rv


def parse_symbols(name_file: str, atoms: dict) -> None:
    """Parse the given symbols file and add them to the existing atom 
    collection.

    :param name_file: File with atomic numbers, symbols, and names for atoms.
    :type name_file: str

    :param atoms: Current collection of atoms. Loaded atoms will be added here.
    :type atoms: dict
    """

    with open(name_file, 'r') as fin:
        for line in fin:
            z, sym, name = line.strip().split()

            if not z in atoms:
                atoms[z] = AtomicData()

            atoms[z].sym = sym
            atoms[z].name = name
            atoms[z].Z = z


def parse_ciaaw_isotopes(iso_file: str, atoms: dict) -> None:
    """Parses an isotope mass file from the Commission on Isotopic Abundances
    and Atomic Weights (CIAAW) and adds it to a given atomic collection.

    :param iso_file: CIAAW isotope mass file
    :type iso_file: str

    :param atoms: Collection of atoms. Loaded isotopes will be added here.
    :type atoms: dict of class:`AtomicData`
    """

    new_atom = r"(\d+)\s+[a-zA-Z]{1,2}\s+[a-zA-Z]+\s+"
    new_iso = r"(\d+)\**\s+(\d+\.\d+)\s((\d+\s?)+)+"
    new_atom += new_iso

    def par_fxn(match: tuple, atom: AtomicData) -> None:
        """Parse the isotope match and add it to the given atom.

        :param match: Regex match groups for an isotope
        :type match: tuple

        :param atom: Atom to add the isotope to
        :type atom: class:`AtomicData`
        """
        as_str = match[1].replace(" ", "") + match[2].replace(" ", "")
        atom.add_isotope(match[0], float(as_str))

    with open(iso_file, 'r') as fin:
        Z = 0
        for line in fin:
            if re.match(new_atom, line):
                match = re.match(new_atom, line).groups()
                Z = match[0]
                par_fxn(match[1:], atoms[Z])
            elif re.match(new_iso, line):
                match = re.match(new_iso, line).groups()
                par_fxn(match, atoms[Z])
            # else: this is not a line containing isotope information


def parse_ciaww_mass(mass_file: str, atoms: dict) -> None:
    """Parses a mass file from the Commission on Isotopic Abundances
    and Atomic Weights (CIAAW) and adds it to a given atomic collection.

    :param iso_file: CIAAW mass file
    :type iso_file: str

    :param atoms: Collection of atoms. Loaded masses will be added here.
    :type atoms: dict of class:`AtomicData`
    """

    mass = r"((?:\d+\.\d+(?:\s\d+)*,*\s?)+)"

    with open(mass_file, 'r') as fin:
        for line in fin:
            if re.search(mass, line):
                Z = line.split()[0]

                masses = re.search(mass, line).groups()[0].replace(" ", "")

                # Split the mass string and convert to floats
                masses = [float(mass) for mass in masses.split(',')]

                # Average masses since it can be single mass or min/max from
                # error bar
                atoms[Z].mass = sum(masses) / len(masses)


def write_load_elements(out_dir: str, amu2me: float, atoms: dict) -> None:
    """Generate a file which will load the atomic info into a 
    chemist::PeriodicTable instance.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of :class:`AtomicData`
    """

    out_file = os.path.join(out_dir, "load_elements.cpp")

    sorted_Z = sorted([int(x) for x in atoms.keys()])

    tab = "    "
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file
        fout.write(
            """#include "chemcache/chemcache.hpp"

namespace chemcache {

void load_elements(chemist::PeriodicTable& pt) {
"""
        )

        # Add atoms and isotopes to the PeriodicTable
        for Z in sorted_Z:
            ai = atoms[str(Z)]

            # Comment atomic number being handled
            fout.write("{}// Z = {}\n".format(tab, Z))

            # Add abundance-weighted atom
            fout.write("{}pt.insert({}, ".format(tab, Z))
            fout.write("chemist::Atom({}ul, {}, \"{}\"));\n".format(
                Z, ai.mass * amu2me, ai.sym
            ))

            # Add isotope atoms
            sorted_mass_numbers = sorted([int(x) for x in ai.isotopes])
            for mn in sorted_mass_numbers:
                mi = ai.isotope_masses[str(mn)] * amu2me

                fout.write("{}pt.add_isotope({}, {}, ".format(tab, Z, mn, ))
                fout.write("chemist::Atom({}ul, {}, \"{}\"));\n"
                           .format(
                               Z, mi, ai.sym
                           ))

            # Extra whitespace between atomic numbers being handled
            fout.write("\n")

        # End of the file
        fout.write(
            """} // function load_elements

} // namespace chemcache
"""
        )


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic info files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    # Get and set some paths
    data_dir = os.path.abspath(args.data_dir)
    out_dir = os.path.abspath(args.src_dir)
    name_file = os.path.join(data_dir, "ElementNames.txt")
    iso_file = os.path.join(data_dir, "CIAAW-ISOTOPEMASSES.txt")
    mass_file = os.path.join(data_dir, "CIAAW-MASSES.txt")
    #cov_file = os.path.join(data_dir, "CovRadii.txt")
    #vdw_file = os.path.join(data_dir, "VanDerWaalRadius.txt")
    #mult_file = os.path.join(data_dir, "NIST-ATOMICION.formatted.txt")

    # Parse atomic data
    atoms = dict()
    parse_symbols(name_file, atoms)
    parse_ciaaw_isotopes(iso_file, atoms)
    parse_ciaww_mass(mass_file, atoms)

    write_load_elements(out_dir, args.amu2me, atoms)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('data_dir', type=str,
                        help="Data directory for atomic information files.")
    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")
    parser.add_argument('--amu2me', type=float,
                        default=1822.888486192,
                        help="""Ratio of mass of electron to one Dalton. 
                             (Default: 1822.888486192)""")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    main(args)
