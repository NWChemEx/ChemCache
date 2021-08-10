#!/usr/bin/env python3

"""This script will read each molecule file in the provided directory
and generate a C++ source file with commands to make each molecule.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

+---src
|       nwx_molecule_manager_pimpl.cpp
"""

import argparse
import os
import re

from generate_atomicinfo import parse_symbols
import helper_fxns as helpers

class Molecule:
    """Representation of a molecule.
    """

    def __init__(self):
        self.carts = []
        self.atoms = []

    def add_atom(self, Z, carts):
        """Adds an atom at the specified cartesian coordinates.

        :param Z: Atomic number
        :type Z: int
        :param carts: Cartesian coordinates
        :type carts: list of float
        """

        self.atoms.append(Z)
        for q in carts:
            self.carts.append(q)

    def __repr__(self):
        """Text representation of the molecule.

        :return: Text representation of the molecule
        :rtype: str
        """

        rv = ""
        for i, ai in enumerate(self.atoms):
            rv += ai + " "
            for j in range(3):
                rv += str(self.carts[i*3 + j]) + " "
            rv += "\n"
        return rv

    def cxxify(self, indent, tab, f):
        """C++ representation of the molecule.

        :param indent: The current indentation
        :type indent: str

        :param tab: The current tab character
        :type tab: str

        :param f: Text IO object to output text.
        :type f: :class:io.TextIOBase
        """

        for i, ai in enumerate(self.atoms):
            f.write(
"""{}mol.push_back(Atom{{mass_t(ptable_.get_atom({}).mass()), {}ul,
{}cart_t{{""".format(indent, ai, ai, indent + tab*4 + "   "))
            line = ""
            line = "{}cart_t{{".format(indent + tab*4 + "   ")
            line+="{}, {},".format(self.carts[i*3], self.carts[i*3+1])
            line +=" {}}},".format(self.carts[i*3+2])            
            f.write("{}, {},".format(self.carts[i*3], self.carts[i*3+1]))
            #Write third coordinate on newline to stay under 80 character column limit
            if(len(line)>80):
                f.write("\n{}{}".format(indent + tab*6 + "  ", self.carts[i*3+2]))
            else:
                f.write(" {}".format(self.carts[i*3+2]))
            f.write("}});")
            if i != len(self.atoms) -1:
                f.write("\n")


def parse_molecules_xyz(filepaths, sym2Z, ang2au):
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

    for filepath in filepaths:
        # Parse each molecule file
        molecule = Molecule()
        with open(filepath, 'r') as fin:
            for line in fin:
                is_match = re.match(an_atom, line)

                if is_match:
                    (sym, str_carts) = is_match.groups()
                    Z = sym2Z[sym.lower()]
                    molecule.add_atom(
                        Z, [float(x)*ang2au for x in str_carts.split()]
                    )

        # Add the molecule to the dictionary of molecules
        molecule_name = os.path.splitext(os.path.basename(filepath))[0]
        molecules[molecule_name] = molecule

    return molecules

def parse_molecules(filepaths, sym2Z, ang2au, extension=".xml"):
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
            "Unsupported molecule file format: {}".format(extension)
        )

def print_source(src_dir, mols):
    tab = "    "
    with open(os.path.join(src_dir, "nwx_molecule_manager_pimpl.cpp"), 'w') as f:
        f.write(
"""/*
 * This file was autogenerated by generate_mols.py. Any edits to it will be
 * lost the next time it is generated.
 */

#include <libchemist/managers/detail_/molecule_manager_pimpl.hpp>
#include <stdexcept>

namespace chemcache::detail_ {

class HardCodedMolsPIMPL : public MoleculeManagerPIMPL {
public:
    using MoleculeManagerPIMPL::MoleculeManagerPIMPL;

protected:
    HardCodedMolsPIMPL(const HardCodedMolsPIMPL& rhs) = default;

private:
    unique_pimpl clone_() const override {
        using unique_me = std::unique_ptr<HardCodedMolsPIMPL>;
        return unique_me(new HardCodedMolsPIMPL(*this));
    }

    value_type at_(const key_type& name) const override {
        using cart_t = typename Atom::coord_type;
        using mass_t = typename Atom::mass_type;
        """)
        for mname, m in sorted(mols.items()):
            f.write("if(name == \"{}\") {{\n{}auto mol = Molecule();\n".format(mname, tab*3))
            m.cxxify(tab*3, tab, f)
            f.write("\n{}return mol;\n{}}} else ".format(tab*3,tab*2))
        f.write("\nthrow std::out_of_range(\"Unknown molecule name\");\n")
        f.write(
"""    } // end at_
};    // end HardCodedMolsPIMPL

std::unique_ptr<MoleculeManagerPIMPL> nwx_default_mols() {
    return std::make_unique<HardCodedMolsPIMPL>();
}

} // namespace chemcache::detail_
""")


def main(args):
    """Entry point function to generate atomic density files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    extensions = [ ".xyz" ]

    # Create some paths
    my_dir    = os.path.dirname(os.path.realpath(__file__))
    name_file = os.path.join(my_dir, "physical_data", "ElementNames.txt")
    src_dir   = os.path.abspath(args.src_dir)
    test_dir  = os.path.abspath(args.test_dir)

    # Discover molecule files
    molecule_dir = os.path.abspath(args.molecule_dir)
    molecule_filepaths = helpers.find_files(
        molecule_dir, extensions, args.recursive
    )

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = {ai.sym.lower() : ai.Z for ai in atoms.values()}

    molecules = {}
    for extension in extensions:
        # NOTE: Extension order CAN matter!
        #       If the same molecule exists in `molecules`
        #       and the new dict returned from `parse_molecules()`, the 
        #       `molecules` version will be replaced by the 
        #       `parse_molecules()` version.
        molecules.update(parse_molecules(
            molecule_filepaths[extension], sym2Z, args.ang2au, extension
        ))
    
    print_source(src_dir, molecules)

def parse_args():
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """
    
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('molecule_dir', type=str,
                        help="""Source directory for molecule files. If 
                             combined with the \"-r\" flag, this directory 
                             will be recursively searched for basis sets.""")
    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")
    parser.add_argument('test_dir', type=str,
                        help="Destination directory for generated unit tests.")

    parser.add_argument('--ang2au', type=float,
                        default=1.8897161646320724,
                        help="""Ratio of angstroms to atomic units. 
                             (Default: 1.8897161646320724)""")
    parser.add_argument('-r', '--recursive', action="store_true",
                        help="""Toggle on recursive search through the basis 
                             set source directory. Default OFF.""")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    main(args)