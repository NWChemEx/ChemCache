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
"""This script parses a file containing ground state atomic electronic
configurations and generates a function to load that data into a PeriodicTable
object.

For readability and convenience we use a few abbreviations throughout this
script:

- Z: the atomic number of an atom
- Sym: the atomic symbol of an atom (e.g. H for hydrogen, He for helium)

Usage
-----

::

   usage: generate_ptable_configs.py [OPTIONS]... data_dir src_dir

   positional arguments:
     data_dir         Data directory for atomic information files.
     src_dir          Destination directory for generated source files.

   optional arguments:
     -h, --help       show this help message and exit

This script looks for the following file(s)::

   <data_dir>
   ├── ElementNames.txt
   └── NIST-ATOMICION.txt

This script creates the following file(s)::

   <src_dir>
   └── atomic_configurations
       └── atomconfigs.cpp
"""

import argparse
import os
import re
from collections import defaultdict

import data_management.helper_fxns as helpers

LMAX = 3
NMAX = 7
i_to_lchar = "spdfgh"[: LMAX + 1]
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
        """Full electronic configuration (number of electrons for each (l, n))

        :return: number of electrons for each (l, n); indexed by [l][n-(l+1)]
        :rtype: list[list[int]]
        """

        confs = []
        for li in range(LMAX + 1):
            tmp = []
            for ni in range(NMAX):
                tmp.append(self.confdict[ni + li + 1, i_to_lchar[li]])
            confs.append(tmp)
        return confs

    @property
    def config(self):
        """Reduced config (number of electrons per l)

        :return: (Ns, Np, Nd, ...)
        :rtype: tuple[int]

        """

        confs = [0] * (LMAX + 1)
        for li in range(LMAX + 1):
            for ni in range(NMAX):
                confs[li] += self.confdict[ni + li + 1, i_to_lchar[li]]
        return tuple(confs)

    def __repr__(self):
        repr_lines = [f"{self.Z} {self.name} {self.sym}"]
        repr_lines.append(f"{self.confstr}")
        repr_lines.append(f"{self.config}")
        repr_lines.append(f"{self.confdict}")
        repr_lines.append(f"{self.config_full}")

        return repr_lines.join("\n")


def parse_symbols(name_file: str, atoms: dict) -> None:
    """Parse the given symbols file and add them to the existing atom
    collection.
    Atoms will be added twice (same data with two keys) to allow access via
    either Z or Sym

    :param name_file: File with atomic numbers, symbols, and names for atoms.
    :type name_file: str

    :param atoms: Current collection of atoms. Loaded atoms will be added here.
    :type atoms: dict

    :raises RuntimeError: Z value was not able to be parsed properly.
    """

    with open(name_file, "r") as fin:
        line_number: int = 0
        for line in fin:
            line_number += 1
            z, sym, name = line.strip().split()

            try:
                z = int(z)
            except ValueError:
                raise RuntimeError(
                    (
                        f"Invalid Z value: Failed to convert Z value of {z} "
                        f"to an integer in {name_file}:{line_number}."
                    )
                )

            if z not in atoms:
                atoms[z] = AtomicData()
            if sym not in atoms:
                atoms[sym] = AtomicData()

            for k in (z, sym):
                atoms[k].sym = sym
                atoms[k].name = name
                atoms[k].Z = z


def parse_config_str(sconf: str) -> dict:
    """Parse an electronic configuration.

    :param sconf: Electronic configuration string
                  (TODO: specify exact format further)
    :type sconf: str

    :raises RuntimeError: Shell was unable to be parsed properly.

    :return: dictionary with number of electrons in each shell, indexed by
             (n, l) where n is an int and l is a str with length 1
    :rtype: dict[tuple[int,str],int]
    """

    conf_dict = defaultdict(int)
    rev_sconf = sconf[::-1]  # reverse config (easier to parse)
    shells = re.findall(r"(\d*[a-z]\d)", rev_sconf)

    # TODO: Do we want the case of no shells found (len(shells) == 0) to be
    #       an exception here?

    for shellrev in shells:
        shell = shellrev[::-1]  # put back in correct order
        matches = re.match(r"(?P<n>\d)(?P<l>\w)(?P<e>\d*)", shell)

        if matches is None:
            # TODO: Make this a more specific exception
            raise RuntimeError(
                (
                    "Invalid shell: unable to parse 'n', 'l', and 'e'."
                    f"Shell: {shell}"
                )
            )

        # 1-elec shells have implicit 'e'
        nelec = int(matches.group("e")) if matches.group("e") else 1
        conf_dict[int(matches.group("n")), matches.group("l")] += nelec

    return dict(conf_dict)


def parse_nist_configs(ip_file: str, atoms: dict) -> None:
    """Parses file containing ionization energy data from the NIST Atomic
    Spectra Database or other file with similar format
    each line should contain "Z element_name atomic_config [...]"

    If an element symbol is used to represent a closed-shell core in a config,
    the config for that element must be parsed and added to `atoms` before
    it is used to represent a core configuration.
    e.g. Ne config must be added before [Ne]3s2

    :param ip_file: path to file containing configs
    :type iso_file: str

    :param atoms: Collection of atoms. Loaded configs will be added here.
    :type atoms: dict of AtomicData

    """

    def parse_cfg_line(line: str, atoms: dict) -> None:
        """Parse a line from the data file and add the config to the
        corresponding atom in the dict atoms.

        :param line: line from NIST file: Z, name, config, term symbol
        :type line: str

        :param atoms: dict of atoms
        :type atoms: dict[tuple[int,str],int]

        :raises RuntimeError: Z value was not able to be parsed properly.
        """

        try:
            # ATOMCONFIGS-HF[S].txt have "Z, name, conf"
            # NIST-ATOMICION.txt has "Z, name, conf, IP sym"
            Z, name, conf_s0 = line.split()[:3]

            try:
                Z = int(Z)
            except ValueError:
                raise RuntimeError(
                    (
                        f"Invalid Z value: Failed to convert Z value of {Z} "
                        f"to an integer in {line}."
                    )
                )

            # get symbol for core if present
            conf_s = conf_s0.strip("[").split("]")

            # get core config as defaultdict
            if len(conf_s) == 2:
                conf_dict = atoms[conf_s[0]].confdict.copy()
            else:
                conf_dict = defaultdict(int)

            # get valence config as dict
            cval_dict = parse_config_str(conf_s[-1])

            # add valence elecs to existing core
            for (n, l), nelec in cval_dict.items():
                conf_dict[n, l] += nelec

            # insert config dict and string into atoms dict
            atoms[Z].confdict = conf_dict.copy()
            atoms[atoms[Z].sym].confdict = atoms[Z].confdict
            atoms[Z].confstr = conf_s0
            atoms[atoms[Z].sym].confstr = atoms[Z].confstr
        # TODO: This should be moved out of this function and specifically name
        #       exceptions that should not be fatal. Using the generic
        #       Exception to catch everything can hide errors and is generally
        #       considered bad practice.
        except Exception:
            print(f"skipping line: {line.strip()}")

    with open(ip_file, "r") as fin:
        # parse config lines
        for line in fin:
            parse_cfg_line(line, atoms)


def _write_configs(out_dir: str, atoms: dict) -> None:
    """Generate the electronic_configurations.cpp source file.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param atoms: Collection of atoms.
    :type atoms: dict of AtomicData
    """

    src_template = """
#include "electronic_configurations.hpp"
#include <simde/simde.hpp>

namespace chemcache {{

using elec_config_pt = simde::ElecConfigFromZ;
using elec_config_t  = std::vector<simde::type::size>;

static constexpr auto module_desc = R"(
Electronic Configurations
---------------------------------

This module returns the electron configuration associated with the atomic
number. This module was autogenerated.
)";

MODULE_CTOR(elec_configs) {{
    description(module_desc);
    satisfies_property_type<elec_config_pt>();
}}

MODULE_RUN(elec_configs) {{
    const auto& [Z] = elec_config_pt::unwrap_inputs(inputs);
    auto rv         = results();

    switch(Z) {{
{cases}
        default: {{
            throw std::out_of_range("Atomic Density not available for Z");
        }}
    }}
}}

}} // namespace chemcache"""

    case_template = """{t}{t}case({Z}): {{
{t}{t}{t}elec_config_t elec_config{{{values}}};
{t}{t}{t}return elec_config_pt::wrap_results(rv, elec_config);
{t}{t}}}"""

    tab = "    "
    sorted_Z = list(sorted([x for x in atoms.keys() if isinstance(x, int)]))
    entries = []
    for Z in sorted_Z:
        ai = atoms[Z]
        values = ", ".join(str(ni) for ni in ai.config)
        entries.append(case_template.format(Z=Z, t=tab, values=values))

    out_file = os.path.join(out_dir, "electronic_configurations.cpp")
    with open(out_file, "w") as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(src_template.format(cases="\n".join(entries)))


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic info files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    # Get and set some paths
    data_dir = os.path.abspath(args.data_dir)
    out_dir = os.path.abspath(args.src_dir)
    name_file = os.path.join(data_dir, "ElementNames.txt")
    ip_file = os.path.join(data_dir, "ATOMCONFIGS-HF.txt")
    # ip_file = os.path.join(data_dir, "ATOMCONFIGS-HFS.txt")
    # ip_file = os.path.join(data_dir, "NIST-ATOMICION.txt")

    # Parse atomic data
    atoms = {}
    parse_symbols(name_file, atoms)
    parse_nist_configs(ip_file, atoms)
    # remove elements without configs
    atoms = {z: a for z, a in atoms.items() if sum(a.config) == z}

    _write_configs(out_dir, atoms)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(
        description=(
            "This script is used to create the experimental data look up "
            "tables for the atom class."
        )
    )

    parser.add_argument(
        "data_dir",
        type=str,
        help="Data directory for atomic information files.",
    )
    parser.add_argument(
        "src_dir",
        type=str,
        help="Destination directory for generated source files.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())
