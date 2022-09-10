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

   usage: generate_densities.py [-h] [-i INC] [-r] atomic_density_dir src_dir test_dir

   positional arguments:
     atomic_density_dir  Source directory for basis set files. If combined with the "-r" flag, this directory will be recursively searched for basis sets.
     src_dir             Destination directory for generated source files.
     test_dir            Destination directory for generated unit tests.

   optional arguments:
     -h, --help          show this help message and exit
     -i INC, --inc INC   Destination include directory, if different than the required "destination" argument.
     -r, --recursive     Toggle on recursive search through the basis set source directory. Default OFF.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

::

   +---include
   |       nwx_atomic_densities.hpp
   |
   +---src
   |   \---atomic_densities
   |           add_density.cmake
   |           <all_basis_set_files>
   |       nwx_atomic_densities.cpp
"""

import argparse
import os
import xml.etree.ElementTree as ET

from generate_atomicinfo import parse_symbols
import helper_fxns as helpers

def _print_pimpl_header(f):
    helpers.write_warning(f, os.path.basename(__file__))

    f.write(
"""#include <string>
#include <stdexcept>

#include \"chemcache/nwx_atomic_densities.hpp\"

namespace chemcache::detail_ {
            
    std::vector<std::vector<double>> get_atomic_density_(const std::string& name, std::size_t Z) {         
""")

def _print_pimpl_footer(f):
    f.write(
"""throw std::out_of_range(\"Basis not available for SAD guess\");
    }//end get_atomic_density_
     
} // namespace chemcache::detail_
""")

def _print_basis_header(f, bs_name):
    helpers.write_warning(f, os.path.basename(__file__))

    f.write(
"""#include <stdexcept>
 
#include \"chemcache/nwx_atomic_densities.hpp\"

namespace chemcache::detail_ {{
 
std::vector<std::vector<double>> {}_density(std::size_t Z) {{
    switch(Z) {{         
""".format(bs_name))

def _print_basis_list(f):
    helpers.write_warning(f, os.path.basename(__file__))

    f.write(
"""#pragma once
#include <vector>
         
namespace chemcache::detail_ {
""")

def _print_basis_footer(f):
    tab = "    "
    f.write(
"""{}default : {{ 
{}throw std::out_of_range(\"Atomic density not available for Z\");
{}}}\n{}}} // end switch\n
}} //end function
}} //end chemcache::detail_""".format(tab*2, tab*3, tab*2, tab))

def _print_atom_basis(f, z, density):
    tab = "    "
    f.write("{}case({}) : {{\n{}return std::vector<std::vector<double>>{{\n".format(tab*2, z, tab*3))
    for line in density.lstrip().rstrip().splitlines():
        f.write("{}{{".format(tab*4))
        for x in line.split():
            f.write(" {},".format(x))
        f.write("},\n")
    f.write("{}}}; //End atomic density\n{}}} //End case\n".format(tab*3, tab*2))

def _write_bases(inc_dir, src_dir, bases):
    tab = "    "
    with open(os.path.join(src_dir,"nwx_atomic_densities.cpp"),'w') as f:
        _print_pimpl_header(f)
        with open(os.path.join(inc_dir, "nwx_atomic_densities.hpp"), 'w') as g:
            _print_basis_list(g)
            f.write("{}".format(tab*2))
            for bs_name, bs in sorted(bases.items()):
                s_name = helpers.sanitize_basis_name(bs_name)
                d_name = helpers.desanitize_basis_name(bs_name)
                f.write("if(name == \"{}\") {{ ".format(d_name))
                f.write("return {}_density(Z); ".format(s_name))
                g.write("std::vector<std::vector<double>> {}_density(std::size_t Z);\n".format(s_name))
                bs_file_name = "{}.cpp".format(bs_name)
                bs_path = os.path.join(src_dir,"atomic_densities", bs_file_name)
                with open(bs_path, 'w') as h:
                    _print_basis_header(h, s_name)
                    for z in sorted([int(x) for x in bs.keys()]):
                        _print_atom_basis(h, z, bs[str(z)])
                    _print_basis_footer(h)
                f.write("}}\n{}else ".format(tab*2))
            g.write("} //end namespace\n")
        _print_pimpl_footer(f)
    
    with open(os.path.join(src_dir,"atomic_densities", "add_density.cmake"), "w") as f:
        helpers.write_warning(f, os.path.basename(__file__), prefix = "# ")

        f.write("set(CHEMCACHE_DENSITY_SOURCE\n")
        for bs_name, bs in sorted(bases.items()):
            f.write("    defaults/atomic_densities/{}.cpp\n".format(bs_name))
        f.write(")")

def _parse_densities_xml(filepaths, sym2Z) -> dict:
    """Parse atomic density files in XML format.

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
                'guessdensitymatrix'
            ).text

    return basis_sets

def _parse_densities(filepaths, sym2Z, extension=".xml") -> dict:
    """Parse atomic density files of the specified format.

    :param filepaths: Full paths to atomic density files.
    :type filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param extension: File format extension to parse, defaults to".xml"
    :type extension: str, optional

    :raises RuntimeError: Unsupported atomic density file format.

    :return: Collection of atomic densities sorted by basis set and element
    :rtype: dict
    """

    if (extension == ".xml"):
        return _parse_densities_xml(filepaths, sym2Z)
    else:
        raise RuntimeError(
            "Unsupported atomic density file format: {}".format(extension)
        )

def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic density files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    extensions = [ ".xml" ]

    # Set include directory to default src_dir path if no path is specified
    if (args.inc == ""):
        args.inc = args.src_dir

    # Create some paths
    my_dir    = os.path.dirname(os.path.realpath(__file__))
    name_file = os.path.join(my_dir, "physical_data", "ElementNames.txt")
    inc_dir   = os.path.abspath(args.inc)
    src_dir   = os.path.abspath(args.src_dir)
    test_dir  = os.path.abspath(args.test_dir)

    # Discover atomic density files
    atomic_density_dir = os.path.abspath(args.atomic_density_dir)
    atomic_density_filepaths = helpers.find_files(
        atomic_density_dir, extensions, args.recursive
    )

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = {ai.sym.lower() : ai.Z for ai in atoms.values()}

    # Gather atomic densities
    basis_sets = {}
    for extension in extensions:
        # NOTE: Extension order CAN matter!
        #       If the same basis set exists in atomic_densities
        #       and the new dict returned from parse_densities(), the 
        #       atomic_densities version will be replaced by the 
        #       parse_densities() version.
        basis_sets.update(_parse_densities(
            atomic_density_filepaths[extension], sym2Z, extension
        ))

    _write_bases(inc_dir, src_dir, basis_sets)

def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """
    
    parser = argparse.ArgumentParser(description=
        "Reads files with atomic densities and write to cpp files."
    )
    
    parser.add_argument('atomic_density_dir', type=str,
                        help="""Source directory for basis set files. If combined
                             with the \"-r\" flag, this directory will be
                             recursively searched for basis sets.""")
    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")
    parser.add_argument('test_dir', type=str,
                        help="Destination directory for generated unit tests.")

    parser.add_argument('-i', '--inc', type=str,
                        default="",
                        help="""Destination include directory, if different 
                             than the required \"destination\" argument.""")
    parser.add_argument('-r', '--recursive', action="store_true",
                        help="""Toggle on recursive search through the basis 
                             set source directory. Default OFF.""")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    main(args)