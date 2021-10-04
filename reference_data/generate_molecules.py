#!/usr/bin/env python3
"""This script will read each molecule file in the provided directory
and generate a C++ source file with commands to make each molecule.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

+---src
|       load_molecules.cpp
"""

import argparse
import io
import os
import re
import sys

from generate_atomicinfo import parse_symbols
import helper_fxns as helpers


class Molecule:
    """Representation of a molecule.
    """
    def __init__(self) -> None:
        self.carts = []
        self.atoms = []

    def add_atom(self, Z: int, carts: list) -> None:
        """Adds an atom at the specified cartesian coordinates.

        :param Z: Atomic number
        :type Z: int
        :param carts: Cartesian coordinates
        :type carts: list of float
        """

        self.atoms.append(Z)
        for q in carts:
            self.carts.append(q)

    def __repr__(self) -> str:
        """Text representation of the molecule.

        :return: Text representation of the molecule
        :rtype: str
        """

        rv = ""
        for i, ai in enumerate(self.atoms):
            rv += ai + " "
            for j in range(3):
                rv += str(self.carts[i * 3 + j]) + " "
            rv += "\n"
        return rv

    def cxxify(self,
               fout: io.TextIOWrapper,
               mol: str,
               pt: str,
               indent: str,
               tab: str = "    ") -> None:
        """C++ representation of the molecule.

        :param fout: Text stream to output the C++ representation to
        :type fout: :class:io.TextIOWrapper

        :param mol: Molecule instance to fill with this molecule's data.
        :type mol: str

        :param pt: PeriodicTable instance to use to get atomic data.
        :type pt: str

        :param indent: The current indentation
        :type indent: str

        :param tab: The current tab character
        :type tab: str

        :return: C++ string representing the molecule
        :rtype: str
        """

        fout.write("{}Molecule {} = Molecule();\n".format(indent, mol))

        for i, ai in enumerate(self.atoms):
            # Start writing atom details
            fout.write("{}{}.push_back(".format(indent, mol))
            fout.write("Atom{{{0}.get_atom({1}).mass(), {1}ul,\n".format(
                pt, ai))

            # Write cartesian coordinates and end atom
            fout.write("{}Atom::coord_type{{{}, {}, {}}}}});\n".format(
                indent + tab, self.carts[i * 3], self.carts[i * 3 + 1],
                self.carts[i * 3 + 2]))


def parse_molecules_xyz(filepaths: list, sym2Z: dict,
                        ang2au: float) -> Molecule:
    """Parses an XYZ formatted molecule file.

    :param file_name: Full paths to molecule files.
    :type file_name: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param ang2au: Ratio of angstroms to atomic units
    :type ang2au: float

    :return: Molecule parsed from file.
    :rtype: :class:Molecule
    """

    an_atom = r"^\s*(\S{1,2})((?:\s+-?\d*.\d+)+)"

    molecules = {}

    # Parse each molecule file
    for filepath in filepaths:

        # Let the user know which molecule is being parsed
        print("Parsing {}...".format(os.path.basename(filepath)), end='')
        sys.stdout.flush()

        molecule = Molecule()
        with open(filepath, 'r') as fin:
            for line in fin:
                is_match = re.match(an_atom, line)

                if is_match:
                    (sym, str_carts) = is_match.groups()
                    Z = sym2Z[sym.lower()]
                    molecule.add_atom(
                        Z, [float(x) * ang2au for x in str_carts.split()])

        # Add the molecule to the dictionary of molecules
        molecule_name = os.path.splitext(os.path.basename(filepath))[0]
        molecules[molecule_name] = molecule

        # Let the user know the molecule is done being parsed
        print("complete")

    return molecules


def parse_molecules(filepaths: list,
                    sym2Z: dict,
                    ang2au: float,
                    extension: str = ".xyz") -> dict:
    """Parse molecule files of the specified format.

    :param filepaths: Full paths to molecule files.
    :type filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param ang2au: Ratio of angstroms to atomic units
    :type ang2au: float

    :param extension: File format extension to parse, defaults to ".xyz"
    :type extension: str, optional

    :raises RuntimeError: Unsupported atomic density file format.

    :return: Collection of molecules
    :rtype: dict of :class:Molecule
    """

    if (extension == ".xyz"):
        return parse_molecules_xyz(filepaths, sym2Z, ang2au)
    else:
        raise RuntimeError(
            "Unsupported molecule file format: {}".format(extension))


def write_load_molecules(src_dir: str, mols: dict, tab: str = "    ") -> None:
    """Write the load_molecules.cpp file with all parsed molecules.

    :param src_dir: ``src`` directory to write load_molecules.cpp to
    :type src_dir: str

    :param mols: Dictionary of names to :class:Molecule
    :type mols: dict of str to :class:Molecule

    :param tab: String representing a tab character
    :type tab: str
    """

    with open(os.path.join(src_dir, "load_molecules.cpp"), 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file
        fout.write("""#include "chemcache/chemcache.hpp"
#include <libchemist/libchemist.hpp>

using namespace libchemist;

namespace chemcache {

void load_molecules(MoleculeManager& mm, const PeriodicTable& pt) {
""")

        # Write molecules
        for name, molecule in sorted(mols.items()):
            molecule.cxxify(fout, name, "pt", tab, tab)
            fout.write("{0}mm.insert(\"{1}\", {1});\n\n".format(tab, name))

        # End of the file
        fout.write("""} // function load_molecules

} // namespace chemcache
""")


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic density files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    extensions = [".xyz"]

    # Create some paths
    my_dir = os.path.dirname(os.path.realpath(__file__))
    name_file = os.path.join(my_dir, "physical_data", "ElementNames.txt")
    src_dir = os.path.abspath(args.src_dir)

    # Discover molecule files
    molecule_dir = os.path.abspath(args.molecule_dir)
    molecule_filepaths = helpers.find_files(molecule_dir, extensions,
                                            args.recursive)

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = {ai.sym.lower(): ai.Z for ai in atoms.values()}

    molecules = {}
    for extension in extensions:
        # NOTE: Extension order CAN matter!
        #       If the same molecule exists in `molecules`
        #       and the new dict returned from `parse_molecules()`, the
        #       `molecules` version will be replaced by the
        #       `parse_molecules()` version.
        molecules.update(
            parse_molecules(molecule_filepaths[extension], sym2Z, args.ang2au,
                            extension))

    print("Writing molecules to file...", end='')
    sys.stdout.flush()

    write_load_molecules(src_dir, molecules)

    print("complete")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('molecule_dir',
                        type=str,
                        help="""Data directory for molecule files. If
                             combined with the \"-r\" flag, this directory
                             will be recursively searched.""")
    parser.add_argument(
        'src_dir',
        type=str,
        help="Destination directory for generated source files.")

    parser.add_argument('--ang2au',
                        type=float,
                        default=1.8897161646320724,
                        help="""Ratio of angstroms to atomic units.
                             (Default: 1.8897161646320724)""")
    parser.add_argument('-r',
                        '--recursive',
                        action="store_true",
                        help="""Toggle on recursive search through molecule_dir
                             directory. Default OFF.""")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    main(args)
