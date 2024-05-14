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
"""This script is used to create the experimental data look up tables for the
atom class.

| Original author: Ben Pritchard
| Modified by: Zachery Crandall

In order to run, this script needs to know the location of the data directory
to read from and the `src` code directory to output the result into. One can
also optionally provide a nondefault value for the electron mass to Dalton ratio.

For readability and convenience we use a few abbreviations throughout this
script:

- Z: the atomic number of an atom
- Sym: the atomic symbol of an atom (e.g. H for hydrogen, He for helium)

Usage
-----

::

   usage: generate_atomicinfo.py [-h] [--amu2me AMU2ME] data_dir src_dir

   positional arguments:
     data_dir         Data directory for atomic information files.
     src_dir          Destination directory for generated source files.

   optional arguments:
     -h, --help       show this help message and exit
     --amu2me AMU2ME  Ratio of mass of electron to one Dalton. (Default: 1822.888486192)

This script looks for the following file(s)::

   +---data_dir
   |       ElementNames.txt
   |       CIAAW-ISOTOPEMASSES.txt
   |       CIAAW-MASSES.txt

This script creates the following file(s)::

   +---src_dir
   |       load_elements.hpp
"""

import argparse
import os
import re

import data_management.helper_fxns as helpers


class AtomicData:

    def __init__(self):
        self.sym = ""
        self.name = ""
        self.Z = 0
        self.mass = 0.0
        self.isotopes = []  # List of isotope mass numbers
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


def _parse_ciaaw_isotopes(iso_file: str, atoms: dict) -> None:
    """Parses an isotope mass file from the Commission on Isotopic Abundances
    and Atomic Weights (CIAAW) and adds it to a given atomic collection.

    :param iso_file: CIAAW isotope mass file
    :type iso_file: str

    :param atoms: Collection of atoms. Loaded isotopes will be added here.
    :type atoms: dict of AtomicData
    """

    new_atom = r"(\d+)\s+[a-zA-Z]{1,2}\s+[a-zA-Z]+\s+"
    new_iso = r"(\d+)\**\s+(\d+\.\d+)\s((\d+\s?)+)+"
    new_atom += new_iso

    def par_fxn(match: tuple, atom: AtomicData) -> None:
        """Parse the isotope match and add it to the given atom.

        :param match: Regex match groups for an isotope
        :type match: tuple

        :param atom: Atom to add the isotope to
        :type atom: AtomicData
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


def _parse_ciaww_mass(mass_file: str, atoms: dict) -> None:
    """Parses a mass file from the Commission on Isotopic Abundances
    and Atomic Weights (CIAAW) and adds it to a given atomic collection.

    :param iso_file: CIAAW mass file
    :type iso_file: str

    :param atoms: Collection of atoms. Loaded masses will be added here.
    :type atoms: dict of AtomicData
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


def _write_z_from_sym(out_dir: str, amu2me: float, atoms: dict) -> None:
    """Generate the Z_from_sym.cpp source file.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData
    """

    # The template we'll be filling values into
    src_template = '''
#include "atoms.hpp"
#include <simde/chemical_system/Z_from_symbol.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using z_pt = simde::ZFromSymbol;
using z_t  = simde::type::atomic_number;

static constexpr auto module_desc = R"(
Atomic Number from Atomic Symbol
--------------------------------

This module returns atomic number associated with the atomic symbol.
This module was autogenerated.
)";

MODULE_CTOR(Z_from_sym) {{
    description(module_desc);
    satisfies_property_type<z_pt>();
}}

MODULE_RUN(Z_from_sym) {{
    const auto& [sym] = z_pt::unwrap_inputs(inputs);

    z_t Z;
    {entries} else {{
        throw std::out_of_range("Z not available for Symbol");
    }}

    auto rv = results();
    return z_pt::wrap_results(rv, Z);
}}

}} // namespace chemcache
'''

    # Gather code entries for Z values
    entries = []
    tab = "    "
    sorted_Z = sorted([int(x) for x in atoms.keys()])
    for Z in sorted_Z:
        ai = atoms[str(Z)]
        entries.append('if(sym == "{}") {{\n{}Z = {};\n{}}}'.format(
            ai.sym, tab + tab, Z, tab))

    # Print the filled out src to the output
    out_file = os.path.join(out_dir, "Z_from_sym.cpp")
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(src_template.format(entries=" else ".join(entries)))


def _write_sym_from_z(out_dir: str, amu2me: float, atoms: dict) -> None:
    """Generate the sym_from_Z.cpp source file.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData
    """

    # The template we'll be filling values into
    src_template = '''
#include "atoms.hpp"
#include <simde/chemical_system/symbol_from_Z.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using sym_pt = simde::SymbolFromZ;
using sym_t  = simde::type::atomic_symbol;

static constexpr auto module_desc = R"(
Atomic Symbol from Atomic Number
---------------------------------

This module returns atomic symbol associated with the atomic number.
This module was autogenerated.
)";

MODULE_CTOR(sym_from_Z) {{
    description(module_desc);
    satisfies_property_type<sym_pt>();
}}

MODULE_RUN(sym_from_Z) {{
    const auto& [Z] = sym_pt::unwrap_inputs(inputs);

    sym_t sym;
    {entries} else {{
        throw std::out_of_range("Symbol not available for Z");
    }}

    auto rv = results();
    return sym_pt::wrap_results(rv, sym);
}}

}} // namespace chemcache
'''

    # Gather code entries for Z values
    entries = []
    tab = "    "
    sorted_Z = sorted([int(x) for x in atoms.keys()])
    for Z in sorted_Z:
        ai = atoms[str(Z)]
        entries.append('if(Z == {}) {{\n{}sym = "{}";\n{}}}'.format(
            Z, tab + tab, ai.sym, tab))

    # Print the filled out src to the output
    out_file = os.path.join(out_dir, "sym_from_Z.cpp")
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(src_template.format(entries=" else ".join(entries)))


def _write_atoms_average(out_dir: str, amu2me: float, atoms: dict) -> None:
    """Generate the atoms_average.cpp source file.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData
    """

    # The template we'll be filling values into
    src_template = '''
#include "atoms.hpp"
#include <simde/chemical_system/atom.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using atom_pt = simde::AtomFromZ;
using atom_t  = simde::type::atom;

static constexpr auto module_desc = R"(
Atoms with Abundance-Weighted Mass
----------------------------------

This module returns an atom with abundance-weighted mass.
This module was autogenerated.
)";

MODULE_CTOR(atoms_average) {{
    description(module_desc);
    satisfies_property_type<atom_pt>();
}}

MODULE_RUN(atoms_average) {{
    const auto& [Z] = atom_pt::unwrap_inputs(inputs);
    auto rv         = results();

    switch(Z) {{
{cases}
        default: {{
            throw std::out_of_range("Atom not available for Z");
        }}
    }}
}}

}} // namespace chemcache
'''

    case_template = """{two_tabs}case({Z}): {{
{three_tabs}atom_t atom{{"{sym}", Z, {mass}, 0.0, 0.0, 0.0}};
{three_tabs}return atom_pt::wrap_results(rv, atom);
{two_tabs}}}"""

    # Gather code entries for Z values
    entries = []
    tab = "    "
    sorted_Z = sorted([int(x) for x in atoms.keys()])
    for Z in sorted_Z:
        ai = atoms[str(Z)]
        entries.append(
            case_template.format(Z=Z,
                                 sym=ai.sym,
                                 mass=ai.mass * amu2me,
                                 two_tabs=tab + tab,
                                 three_tabs=tab + tab + tab))

    # Print the filled out src to the output
    out_file = os.path.join(out_dir, "atoms_average.cpp")
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(src_template.format(cases="\n".join(entries)))


def _write_atoms_isotope(out_dir: str, amu2me: float, atoms: dict) -> None:
    """Generate the atoms_isotopes.cpp source file.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData
    """

    # The template we'll be filling values into
    src_template = '''
#include "atoms.hpp"
#include <simde/chemical_system/atom.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using z_t        = simde::type::atomic_number;
using atom_t     = simde::type::atom;
using isotope_pt = simde::Atom<std::pair<z_t, z_t>>;

static constexpr auto module_desc = R"(
Atoms with Isotope Mass
---------------------------------

This module returns an instance of an atomic isotope.
This module was autogenerated.
)";

MODULE_CTOR(atoms_isotope) {{
    description(module_desc);
    satisfies_property_type<isotope_pt>();
}}

MODULE_RUN(atoms_isotope) {{
    const auto& [Z_and_N] = isotope_pt::unwrap_inputs(inputs);
    auto rv               = results();

    auto [Z, N]   = Z_and_N;
    auto message1 = "Isotopes not available for Z";
    auto message2 = "Isotope not available for Z and mass number";
    atom_t atom;
    {entries} else {{
        throw std::out_of_range(message1);
    }}

    return isotope_pt::wrap_results(rv, atom);
}}

}} // namespace chemcache
'''

    entry_template = '''if({n} == {v}) {{
{t1}{c}
{t2}}}'''

    atom_ctor = 'atom = atom_t{{"{}", Z, {}, 0.0, 0.0, 0.0}};'

    error_msg = "{{\n{}throw std::out_of_range(message2);\n{}}}"

    # Gather code entries for Z values
    entries = []
    tab = "    "
    sorted_Z = sorted([int(x) for x in atoms.keys()])
    for Z in sorted_Z:
        # Gather code entries for N values
        ai = atoms[str(Z)]
        internal_entries = []
        sorted_mass_numbers = sorted([int(x) for x in ai.isotopes])
        if not sorted_mass_numbers:
            continue
        for mn in sorted_mass_numbers:
            mi = ai.isotope_masses[str(mn)] * amu2me
            atom = atom_ctor.format(ai.sym, mi)
            internal_entries.append(
                entry_template.format(n="N",
                                      v=mn,
                                      c=atom,
                                      t1=tab * 3,
                                      t2=tab * 2))
        internal_entries.append(error_msg.format(tab * 3, tab * 2))
        entries.append(
            entry_template.format(n="Z",
                                  v=Z,
                                  c=' else '.join(internal_entries),
                                  t1=tab * 2,
                                  t2=tab))

    # Print the filled out src to the output
    out_file = os.path.join(out_dir, "atoms_isotope.cpp")
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(src_template.format(entries=" else ".join(entries)))


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic info files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    # Get and set some paths
    data_dir = os.path.abspath(args.data_dir)
    out_dir = os.path.abspath(args.src_dir)
    name_file = os.path.join(data_dir, "ElementNames.txt")
    iso_file = os.path.join(data_dir, "CIAAW-ISOTOPEMASSES.txt")
    mass_file = os.path.join(data_dir, "CIAAW-MASSES.txt")
    # cov_file = os.path.join(data_dir, "CovRadii.txt")
    # vdw_file = os.path.join(data_dir, "VanDerWaalRadius.txt")
    # mult_file = os.path.join(data_dir, "NIST-ATOMICION.formatted.txt")

    # Parse atomic data
    atoms = dict()
    parse_symbols(name_file, atoms)
    _parse_ciaaw_isotopes(iso_file, atoms)
    _parse_ciaww_mass(mass_file, atoms)

    _write_z_from_sym(out_dir, args.amu2me, atoms)
    _write_sym_from_z(out_dir, args.amu2me, atoms)
    _write_atoms_average(out_dir, args.amu2me, atoms)
    _write_atoms_isotope(out_dir, args.amu2me, atoms)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description=
        "This script is used to create the experimental data look up tables "
        "for the atom class.")

    parser.add_argument('data_dir',
                        type=str,
                        help="Data directory for atomic information files.")
    parser.add_argument(
        'src_dir',
        type=str,
        help="Destination directory for generated source files.")
    parser.add_argument('--amu2me',
                        type=float,
                        default=1822.888486192,
                        help="""Ratio of mass of electron to one Dalton. 
                             (Default: 1822.888486192)""")

    return parser.parse_args()


if __name__ == '__main__':
    main(parse_args())
