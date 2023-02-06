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


"""This script is used to create the data array containing ground state atomic
electronic configurations used for the SAD guess module.

For readability and convenience we use a few abbreviations throughout this
script:

- Z: the atomic number of an atom
- Sym: the atomic symbol of an atom (e.g. H for hydrogen, He for helium)

Usage
-----

::

   usage: generate_atomconfigs.py [-h] data_dir src_dir

   positional arguments:
     data_dir         Data directory for atomic information files.
     src_dir          Destination directory for generated source files.

   optional arguments:
     -h, --help       show this help message and exit

This script looks for the following file(s)::

   +---data_dir  
   |       ElementNames.txt  
   |       NIST-ATOMICION.txt  

This script creates the following file(s)::

   +---src_dir/atomic_configurations  
   |       atomconfigs.cpp  
"""

import argparse
import os
import re
from collections import defaultdict

import helper_fxns as helpers

LMAX = 3
NMAX = 7
i_to_lchar = 'spdfgh'[:LMAX+1]
lchar_to_i = {c: i for i, c in enumerate(i_to_lchar)}

class AtomicData:
    def __init__(self):
        self.sym = ""
        self.name = ""
        self.Z = 0
        self.confstr = ""
        self.confdict = defaultdict(int)

    @property
    def config_full(self):
        confs = []
        for li in range(LMAX+1):
            tmp = []
            for ni in range(NMAX):
                tmp.append(self.confdict[ni+li+1,i_to_lchar[li]])
            confs.append(tmp)
        return confs

    @property
    def config(self):
        confs = [0]*(LMAX+1)
        for li in range(LMAX+1):
            for ni in range(NMAX):
                confs[li] += self.confdict[ni+li+1,i_to_lchar[li]]
        return tuple(confs)

    def __repr__(self):
        return f"{self.Z} {self.name} {self.sym}\n{self.confstr}\n{self.config}\n{self.confdict}\n{self.config_full}"


def parse_symbols(name_file: str, atoms: dict) -> None:
    """Parse the given symbols file and add them to the existing atom 
    collection.
    Atoms will be added twice (same data with two keys) to allow access via
    either Z or Sym

    :param name_file: File with atomic numbers, symbols, and names for atoms.
    :type name_file: str

    :param atoms: Current collection of atoms. Loaded atoms will be added here.
    :type atoms: dict
    """

    with open(name_file, 'r') as fin:
        for line in fin:
            z, sym, name = line.strip().split()
            z = int(z)

            if not z in atoms:
                atoms[z] = AtomicData()
            if not sym in atoms:
                atoms[sym] = AtomicData()

            for k in (z, sym):
                atoms[k].sym = sym
                atoms[k].name = name
                atoms[k].Z = z


def parse_nist_configs(ip_file: str, atoms: dict) -> None:
    """Parses file containing ionization energy data from the NIST Atomic
    Spectra Database
    each line should contain "Z element_name atomic_config term_symbol"

    :param ip_file: NIST Ionization Energies Data file
    :type iso_file: str

    :param atoms: Collection of atoms. Loaded configs will be added here.
    :type atoms: dict of AtomicData
    """

    new_atom = r"(\d+)\s+[a-zA-Z]{1,2}\s+[a-zA-Z]+\s+"
    new_iso = r"(\d+)\**\s+(\d+\.\d+)\s((\d+\s?)+)+"
    new_atom += new_iso

    def parse_cfg_line(line: str, atoms: dict) -> None:
        """Parse the isotope match and add it to the given atom.

        :param match: Regex match groups for an isotope
        :type match: tuple

        :param atom: Atom to add the isotope to
        :type atom: AtomicData
        """
        try:
            Z, name, conf_s0, _ = line.split()
            Z = int(Z)
            conf_s = conf_s0.strip('[').split(']')  # separate core if present
            # get core config
            if len(conf_s) == 2:
                conf_dict = atoms[conf_s[0]].confdict.copy()
            else:
                conf_dict = defaultdict(int)

            # get string representing remaining elec config
            cval_s = conf_s[-1]
            revcval_s = cval_s[::-1]  # reverse config (easier to parse)
            shells = re.findall(r'(\d*[a-z]\d)', revcval_s)
            for shellrev in shells:
                shell = shellrev[::-1]  # put back in correct order
                matches = re.match(r'(?P<n>\d)(?P<l>\w)(?P<e>\d*)', shell)
                # 1-elec shells have implicit 'e'
                nelec = int(matches['e']) if matches['e'] else 1
                conf_dict[int(matches['n']),matches['l']] += nelec
            atoms[Z].confdict = conf_dict.copy()
            atoms[atoms[Z].sym].confdict = atoms[Z].confdict
            atoms[Z].confstr = conf_s0 
            atoms[atoms[Z].sym].confstr = atoms[Z].confstr
        except Exception as e:
            print(f"skipping header line: {line.strip()}")

    with open(ip_file, 'r') as fin:
        for line in fin:
            parse_cfg_line(line, atoms)


def _write_pt_configs(out_dir: str, atoms: dict) -> None:
    """Generate a file containing a function which returns the atomic
    configurations as a std::array<std::array<size_t,n_l>,n_elements>

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData {Z: (config, Sym, name)
    """

    out_file = os.path.join(out_dir, "load_elec_configs.cpp")

    sorted_Z = list(sorted([x for x in atoms.keys() if isinstance(x, int)]))
    n_elements = len(sorted_Z)
    n_l = len(atoms[sorted_Z[0]].config)

    tab = "    "
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file
        fout.write(
            """#include "chemcache/chemcache.hpp"
#include <array>

namespace chemcache {

void load_elec_configs(chemist::PeriodicTable& pt) {
"""
#    return std::array<std::array<size_t, N_L>, N_ELEMENTS> {{
#""".replace("N_L",str(n_l)).replace("N_ELEMENTS",str(n_elements))
        )

        # Add atoms and isotopes to the PeriodicTable
        for Z in sorted_Z:
            #conf_i, sym_i, name_i = atoms[Z]
            ai = atoms[Z]
            # Comment atomic number, symbol, name
            comment_str = f"// Z = {Z:>3d}, {ai.sym:<2s}, {ai.name:<14s}"

            # Print config and comment
            fout.write(tab + comment_str + "\n")

            fout.write(
                tab + f"pt.add_elec_config({Z:>3d}, {{"
            + ",".join(f"{ni:>2d}" for ni in ai.config)
            + "});\n\n") 


        # End of the file
        fout.write(
            """} // function load_elec_configs

} // namespace chemcache
"""
        )

def _write_array_configs(out_dir: str, atoms: dict) -> None:
    """Generate a file containing a function which returns the atomic
    configurations as a std::array<std::array<size_t,n_l>,n_elements>

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData {Z: (config, Sym, name)
    """

    out_file = os.path.join(out_dir, "atomic_configurations", "atomconfigs.cpp")

    sorted_Z = list(sorted([x for x in atoms.keys() if isinstance(x, int)]))
    n_elements = len(sorted_Z)
    n_l = len(atoms[sorted_Z[0]].config)

    tab = "    "
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file
        fout.write(
            """#include "chemcache/chemcache.hpp"
#include <array>

namespace chemcache {

inline auto atomconfigs() {
    return std::array<std::array<size_t, N_L>, N_ELEMENTS> {{
""".replace("N_L",str(n_l)).replace("N_ELEMENTS",str(n_elements))
        )

        # Add atoms and isotopes to the PeriodicTable
        for Z in sorted_Z:
            #conf_i, sym_i, name_i = atoms[Z]
            ai = atoms[Z]
            # Comment atomic number, symbol, name
            comment_str = f" // Z = {Z:>3d}, {ai.sym:<2s}, {ai.name:<14s}"

            # Print config and comment
            fout.write(tab + "{" 
            + ", ".join(f"{ni:>2d}" for ni in ai.config)
            + "},"
            + comment_str + "\n") 


        # End of the file
        fout.write(
            tab + """}};
} // function atomconfigs

} // namespace chemcache
"""
        )

def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic info files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    # Get and set some paths
    data_dir = os.path.abspath(args.data_dir)
    out_dir = os.path.abspath(args.src_dir)
    name_file = os.path.join(data_dir, "ElementNames.txt")
    ip_file = os.path.join(data_dir, "NIST-ATOMICION.txt")
    #iso_file = os.path.join(data_dir, "CIAAW-ISOTOPEMASSES.txt")
    #mass_file = os.path.join(data_dir, "CIAAW-MASSES.txt")
    #cov_file = os.path.join(data_dir, "CovRadii.txt")
    #vdw_file = os.path.join(data_dir, "VanDerWaalRadius.txt")
    #mult_file = os.path.join(data_dir, "NIST-ATOMICION.formatted.txt")

    # Parse atomic data
    atoms = dict()
    parse_symbols(name_file, atoms)
    parse_nist_configs(ip_file, atoms)
    # add dummy element zero for nicer indexing
    #atoms[0].config = (0,0,0,0)
    # remove elements without configs
    # TODO: check that remaining elements are contiguous from 0?
    atoms = {z:a for z,a in atoms.items() if sum(a.config) == z}
    for z,atom in atoms.items():
        print(atom)

    _write_pt_configs(out_dir, atoms)
    _write_array_configs(out_dir, atoms)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description=
        "This script is used to create the experimental data look up tables "
        "for the atom class."
    )

    parser.add_argument('data_dir', type=str,
                        help="Data directory for atomic information files.")
    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    main(args)
