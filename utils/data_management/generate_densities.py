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
"""Reads files with atomic densities and write to cpp files.

Usage
-----

::

   usage: generate_densities.py [-h] [-r] [-a ATOMS_DIR] atomic_density_dir src_dir

   positional arguments:
     atomic_density_dir    Source directory for basis set files. If combined with the "-r" flag, this directory will be recursively searched for basis sets.
     src_dir               Destination directory for generated source files.
   
   options:
     -h, --help            show this help message and exit
     -r, --recursive       Toggle on recursive search through the basis set source directory. Default OFF.
     -a ATOMS_DIR, --atoms_dir ATOMS_DIR
                           The path to where ElementNames.txt can be found.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

::

   +---include
   |       nwx_atomic_densities.hpp
   |
   +---src
   |   \---atomic_densities (this directory must be made manually before run this python script)
   |           add_density.cmake
   |           <all_basis_set_files>
   |       nwx_atomic_densities.cpp
"""

import argparse
import os
import xml.etree.ElementTree as ET

from data_management.generate_atomicinfo import parse_symbols
import data_management.helper_fxns as helpers


def make_square_arr(a: list, spacer='\n'):
    n2 = len(a)
    n = int(n2**0.5)
    if n == 1:
        return a[0]
    assert (n * n == n2)
    return f"{spacer}{{" + f"}},{spacer}{{".join(
        [", ".join(a[i * n:(i + 1) * n]) for i in range(n)]) + "}"


def _write_den_files(out_file: str,
                     bs_name: str,
                     basis_set: dict,
                     tab: str = "    ") -> None:

    source_template = '''
#include "atomic_densities.hpp"
#include <simde/simde.hpp>

namespace chemcache {{

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using atomic_den_pt   = simde::AtomDenFromZ;
using atomic_den_t    = simde::type::el_density;
using tensor_t        = simde::type::tensor;
using atomic_bs_t     = simde::type::atomic_basis_set;
using ao_bs_t         = simde::type::ao_basis_set;
using aos_t           = simde::type::ao_space;

static constexpr auto module_desc = R"(
{bs_name} atomic densities
---------------------------------

This module returns precomputed atomic densities in this basis set.
This module was autogenerated.
)";

MODULE_CTOR({s_name}_atom_dm) {{
    description(module_desc);
    satisfies_property_type<atomic_den_pt>();
    add_submodule<atomic_basis_pt>("AO Basis");

}}

MODULE_RUN({s_name}_atom_dm) {{
    const auto& [Z] = atomic_den_pt::unwrap_inputs(inputs);
    auto& ao_mod    = submods.at("AO Basis");
    auto atom_bs    = ao_mod.run_as<atomic_basis_pt>(Z);
    // convert atomic_basis_set -> ao_basis_set -> ao_space
    ao_bs_t ao_bs;
    ao_bs.add_center(atom_bs);
    aos_t aos(ao_bs);
    auto rv         = results();

    switch(Z) {{
{cases}
        default: {{
            throw std::out_of_range("Atomic Density not available for Z");
        }}
    }}
}}

}} // namespace chemcache
'''

    cases_template = '''{t}{t}case({Z}): {{
{t}{t}{t}tensor_t atom_dm{{{values}}};
{t}{t}{t}atomic_den_t rho(atom_dm, aos);
{t}{t}{t}return atomic_den_pt::wrap_results(rv, rho);
{t}{t}}}'''

    s_name = helpers.sanitize_basis_name(bs_name)

    cases = []
    for z in sorted([int(x) for x in basis_set.keys()]):
        # values = ", ".join(basis_set[str(z)].split())
        values = make_square_arr(basis_set[str(z)].split(),
                                 spacer="\n" + tab * 4)
        cases.append(
            cases_template.format(t=tab, Z=z, bs_name=bs_name, values=values))

    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(
            source_template.format(bs_name=bs_name,
                                   s_name=s_name,
                                   cases="\n".join(cases)))


def _write_densities(src_dir: str, bases: dict, tab="    ") -> None:
    """Writes basis set data to C++ files.

    :param src_dir: Source directory for source files.
    :type src_dir: str

    :param bases: Collection of basis sets parsed from files
    :type bases: dict

    :param tab: String representing a tab, defaults to "    "
        :type tab: str, optional
    """

    bases_template = """
#pragma once
#include <pluginplay/pluginplay.hpp>

namespace chemcache {{

// Module declarations will go here
{d}

namespace atom_dm_mods {{

inline void set_defaults(pluginplay::ModuleManager& mm) {{
    // Default submodules within this subcollection will be set here
    {ao}
}}

inline void load_modules(pluginplay::ModuleManager& mm) {{
    // Modules will be added to the ModuleManager here
    {m}

    set_defaults(mm);
}}

}} // namespace atom_dm_mods

}} // namespace chemcache
"""

    d_template = "DECLARE_MODULE({}_atom_dm);"
    m_template = 'mm.add_module<{}_atom_dm>("{} atomic dm");'
    ao_template = 'mm.change_submod("{} atomic dm", "AO Basis", "{} atomic basis");'

    ntab = "\n" + tab

    ds = []
    ms = []
    aos = []

    for bs_name, basis_set in sorted(bases.items()):
        d_name = helpers.desanitize_basis_name(bs_name)
        den_file = os.path.join(src_dir, bs_name + ".cpp")
        _write_den_files(den_file, bs_name, basis_set)

        s_name = helpers.sanitize_basis_name(bs_name)
        ds.append(d_template.format(s_name))
        ms.append(m_template.format(s_name, d_name))
        aos.append(ao_template.format(d_name, d_name))

    bases_file = os.path.join(src_dir, "atomic_densities.hpp")
    with open(bases_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(
            bases_template.format(d="\n".join(ds),
                                  m=ntab.join(ms),
                                  ao=ntab.join(aos)))


def _parse_densities_xml(filepaths, sym2Z) -> dict:
    """Parse atomic density files in XML format. 
    Deprecated as new atomic density format is adopted. Y. Z. 02/23/23

    :param filepaths: Full paths to atomic density files.
    :type filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :return: Collection of atomic densities sorted by basis set and element
    :rtype: dict
    """

    basis_sets = {}
    for filepath in filepaths:
        basis_set = os.path.splitext(os.path.basename(filepath))[0]

        basis_sets[basis_set] = {}

        root = ET.parse(filepath)

        for atom in root.findall('atomicguess'):
            atom_z = sym2Z[atom.get('symbol').lower()]

            basis_sets[basis_set][atom_z] = atom.find(
                'guessdensitymatrix').text

    return basis_sets


def _parse_densities_dat(filepaths, sym2Z) -> dict:
    """Parse atomic density files in .dat format. 

    :param filepaths: Full paths to atomic density files.
    :type filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :return: Collection of atomic densities sorted by basis set and element
    :rtype: dict
    """

    basis_sets = {}
    for filepath in filepaths:
        basis_set = os.path.splitext(os.path.basename(filepath))[0]

        basis_sets[basis_set] = {}

        with open(filepath, "r") as f:
            for line in f:
                line1 = line.strip().split()
                if ((len(line1) == 1) and (line1[0].isalpha())):
                    atom_z = sym2Z[line1[0].lower()]
                    guessDM = ''
                    line2 = f.readline().strip()
                    while (len(line2) > 0):
                        guessDM += (line2 + '\n')
                        line2 = f.readline().strip()
                basis_sets[basis_set][atom_z] = guessDM

    return basis_sets


def _parse_densities(filepaths, sym2Z, extension=".dat") -> dict:
    """Parse atomic density files of the specified format.

    :param filepaths: Full paths to atomic density files.
    :type filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param extension: File format extension to parse, defaults to".dat"
    :type extension: str, optional

    :raises RuntimeError: Unsupported atomic density file format.

    :return: Collection of atomic densities sorted by basis set and element
    :rtype: dict
    """

    if (extension == ".dat"):
        return _parse_densities_dat(filepaths, sym2Z)
    elif (extension == ".xml"):
        return _parse_densities_xml(filepaths, sym2Z)
    else:
        raise RuntimeError(
            "Unsupported atomic density file format: {}".format(extension))


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic density files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    extensions = [".dat"]

    # Create some paths
    my_dir = os.path.dirname(os.path.realpath(__file__))
    src_dir = os.path.abspath(args.src_dir)
    name_file = os.path.abspath(
        os.path.join(args.atoms_dir, "ElementNames.txt"))

    # Discover atomic density files
    atomic_density_dir = os.path.abspath(args.atomic_density_dir)
    atomic_density_filepaths = helpers.find_files(atomic_density_dir,
                                                  extensions, args.recursive)

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = {ai.sym.lower(): ai.Z for ai in atoms.values()}

    # Gather atomic densities
    basis_sets = {}
    for extension in extensions:
        # NOTE: Extension order CAN matter!
        #       If the same basis set exists in atomic_densities
        #       and the new dict returned from parse_densities(), the
        #       atomic_densities version will be replaced by the
        #       parse_densities() version.
        basis_sets.update(
            _parse_densities(atomic_density_filepaths[extension], sym2Z,
                             extension))

    _write_densities(src_dir, basis_sets)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser(
        description="Reads files with atomic densities and write to cpp files."
    )

    parser.add_argument(
        'atomic_density_dir',
        type=str,
        help="""Source directory for basis set files. If combined
                             with the \"-r\" flag, this directory will be
                             recursively searched for basis sets.""")

    parser.add_argument(
        'src_dir',
        type=str,
        help="Destination directory for generated source files.")

    parser.add_argument('-r',
                        '--recursive',
                        action="store_true",
                        help="""Toggle on recursive search through the basis 
                             set source directory. Default OFF.""")

    parser.add_argument(
        "-a",
        "--atoms_dir",
        action="store",
        type=str,
        default="reference_data/physical_data",
        help="The path to where ElementNames.txt can be found.")

    return parser.parse_args()


if __name__ == "__main__":
    main(parse_args())