#!/usr/bin/env python3
# Copyright 2022 NWChemEx-Project
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""This script will read each molecule file in the provided directory
and generate a C++ source file with commands to make each molecule.

Usage
-----

::

   usage: generate_molecules.py [-h] [--ang2au ANG2AU] [-r] [-a ATOMS_DIR] molecule_dir src_dir

   positional arguments:
     molecule_dir          Data directory for molecule files. If combined with the "-r" flag, this directory will be recursively searched.
     src_dir               Destination directory for generated source files.

   options:
     -h, --help            show this help message and exit
     --ang2au ANG2AU       Ratio of angstroms to atomic units. (Default: 1.8897161646320724)
     -r, --recursive       Toggle on recursive search through molecule_dir directory. Default OFF.
     -a ATOMS_DIR, --atoms_dir ATOMS_DIR
                           The path to where ElementNames.txt can be found.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

::

   +---src
   |       load_molecules.cpp
"""

import argparse
import io
import os
import re
import sys

import data_management.helper_fxns as helpers
from data_management.generate_atomicinfo import parse_symbols


class Molecule:
    """Representation of a molecule."""

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

    def cxxify(self, tab: str = "    ") -> str:
        """C++ representation of the molecule.

        :param indent: The current indentation
        :type indent: str

        :param tab: The current tab character, defaults to "    "
        :type tab: str, optional

        :return: C++ string representing the molecule
        :rtype: str
        """

        atom = "{t}{t}atm = atoms_mod.run_as<atom_pt>(z_t{{{z}}});"
        xyz = "{t}{t}translate(atm, {x}, {y}, {z});"
        push = "{t}{t}mol.push_back(atm);\n"

        lines = []
        for i, ai in enumerate(self.atoms):
            # Start writing atom details
            x = self.carts[i * 3]
            y = self.carts[i * 3 + 1]
            z = self.carts[i * 3 + 2]
            lines.append(atom.format(t=tab, z=ai))
            lines.append(xyz.format(t=tab, x=x, y=y, z=z))
            lines.append(push.format(t=tab))

        return "\n".join(lines)


def _parse_molecules_xyz(
    filepaths: list, sym2Z: dict, ang2au: float
) -> Molecule:
    """Parses an XYZ formatted molecule file.

    :param file_name: Full paths to molecule files.
    :type file_name: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param ang2au: Ratio of angstroms to atomic units
    :type ang2au: float

    :return: Molecule parsed from file.
    :rtype: Molecule
    """

    an_atom = r"^\s*(\S{1,2})((?:\s+-?\d*.\d+)+)"

    molecules = {}

    # Parse each molecule file
    for filepath in filepaths:
        # Let the user know which molecule is being parsed
        print("Parsing {}...".format(os.path.basename(filepath)), end="")
        sys.stdout.flush()

        molecule = Molecule()
        with open(filepath, "r") as fin:
            for line in fin:
                is_match = re.match(an_atom, line)

                if is_match:
                    (sym, str_carts) = is_match.groups()
                    Z = sym2Z[sym.lower()]
                    molecule.add_atom(
                        Z, [float(x) * ang2au for x in str_carts.split()]
                    )

        # Add the molecule to the dictionary of molecules
        molecule_name = os.path.splitext(os.path.basename(filepath))[0]
        molecules[molecule_name] = molecule

        # Let the user know the molecule is done being parsed
        print("complete")

    return molecules


def _parse_molecules(
    filepaths: list, sym2Z: dict, ang2au: float, extension: str = ".xyz"
) -> dict:
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
    :rtype: dict of Molecule
    """

    if extension == ".xyz":
        return _parse_molecules_xyz(filepaths, sym2Z, ang2au)
    else:
        raise RuntimeError(
            "Unsupported molecule file format: {}".format(extension)
        )


def _write_load_molecules(src_dir: str, mols: dict, tab: str = "    ") -> None:
    """Write the load_molecules.cpp file with all parsed molecules.

    :param src_dir: ``src`` directory to write load_molecules.cpp to
    :type src_dir: str

    :param mols: Dictionary of names to :class:Molecule
    :type mols: dict of str to :class:Molecule

    :param tab: String representing a tab character
    :type tab: str
    """

    # The template we'll be filling values into
    src_template = """
#include "molecules.hpp"
#include <simde/chemical_system/atom.hpp>
#include <simde/chemical_system/molecule_from_string.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using molecule_pt = simde::MoleculeFromString;
using atom_pt     = simde::AtomFromZ;
using molecule_t  = simde::type::molecule;
using atom_t      = simde::type::atom;
using z_t         = simde::type::atomic_number;

static constexpr auto module_desc = R"(
NWChemEx Molecules
---------------------------------

This module returns an instance of a particular module.
This module was autogenerated.
)";

MODULE_CTOR(default_molecules) {{
    description(module_desc);
    satisfies_property_type<molecule_pt>();
    add_submodule<atom_pt>("Atoms");
}}

MODULE_RUN(default_molecules) {{
    const auto& [name] = molecule_pt::unwrap_inputs(inputs);
    auto& atoms_mod    = submods.at("Atoms");

    auto translate = [](auto& a, double x, double y, double z) {{
        a.x() = x;
        a.y() = y;
        a.z() = z;
    }};

    atom_t atm;
    molecule_t mol;
    {entries} else {{
        throw std::out_of_range("No molecule found for name: " + name);
    }}

    auto rv = results();
    return molecule_pt::wrap_results(rv, mol);
}}

}} // namespace chemcache
"""

    entry_template = 'if(name == "{n}") {{\n{c}\n{t}}}'

    # Write molecules
    entries = []
    for name, molecule in sorted(mols.items()):
        contents = molecule.cxxify(tab)
        entries.append(entry_template.format(n=name, c=contents, t=tab))

    with open(os.path.join(src_dir, "molecules.cpp"), "w") as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(src_template.format(entries=" else ".join(entries)))


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic density files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    extensions = [".xyz"]

    # Create some paths
    my_dir = os.path.dirname(os.path.realpath(__file__))
    src_dir = os.path.abspath(args.src_dir)
    name_file = os.path.abspath(
        os.path.join(args.atoms_dir, "ElementNames.txt")
    )

    # Discover molecule files
    molecule_dir = os.path.abspath(args.molecule_dir)
    molecule_filepaths = helpers.find_files(
        molecule_dir, extensions, args.recursive
    )

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
            _parse_molecules(
                molecule_filepaths[extension], sym2Z, args.ang2au, extension
            )
        )

    print("Writing molecules to file...", end="")
    sys.stdout.flush()

    _write_load_molecules(src_dir, molecules)

    print("complete")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """

    parser = argparse.ArgumentParser(
        description="This script will read each molecule file in the provided directory"
        "and generate a C++ source file with commands to make each molecule."
    )

    parser.add_argument(
        "molecule_dir",
        type=str,
        help="""Data directory for molecule files. If
                             combined with the \"-r\" flag, this directory
                             will be recursively searched.""",
    )

    parser.add_argument(
        "src_dir",
        type=str,
        help="Destination directory for generated source files.",
    )

    parser.add_argument(
        "--ang2au",
        type=float,
        default=1.8897161646320724,
        help="""Ratio of angstroms to atomic units.
                             (Default: 1.8897161646320724)""",
    )

    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="""Toggle on recursive search through molecule_dir
                             directory. Default OFF.""",
    )

    parser.add_argument(
        "-a",
        "--atoms_dir",
        action="store",
        type=str,
        default="reference_data/physical_data",
        help="The path to where ElementNames.txt can be found.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
