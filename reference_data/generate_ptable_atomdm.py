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


"""This script generates a function to load pre-calculated atomic density 
matrices into a PeriodicTable object.
Usage
-----
::
   usage: generate_ptable_atomdm.py [-h] src_dir
   positional arguments:
     src_dir          Destination directory for generated source files.
   optional arguments:
     -h, --help       show this help message and exit
This script creates the following file:
   +---src_dir/load_atom_dm.cpp  
"""

import argparse
import os

import helper_fxns as helpers

def _write_pt_atomdm(out_dir: str) -> None:
    """Generate a file containing a function which loads pre-calculated
    atomic density matrices to an existing PeriodicTable object
    :param out_dir: Output directory for the generated header file.
    :type out_dir: str
    """

    out_file = os.path.join(out_dir, "load_atom_dm.cpp")

    tab = "    "
    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file
        fout.write(
            """#include "chemcache/chemcache.hpp"
#include <chemist/chemist.hpp>
#include "nwx_atomic_densities.hpp"

namespace chemcache {

    void load_atom_dm(chemist::PeriodicTable::size_type Z, std::string basis_name, chemist::PeriodicTable& pt) {
        pt.add_atom_dm(Z, basis_name, chemcache::detail_::get_atomic_density_(basis_name, Z));
    } // function load_atom_dm for the non-default case

    void load_atom_dm(chemist::PeriodicTable& pt) {
""")

        # Default case: loading the STO-3G densities to the PeriodicTable
        for Z in range(1,54):
            fout.write(tab + tab +"pt.add_atom_dm({}, \"STO-3G\", chemcache::detail_::get_atomic_density_(\"STO-3G\", {}));\n".format(Z,Z))
        # End of the file
        fout.write(
            tab+"""} // function load_atom_dm for the default case

} // namespace chemcache
""")

def main(args: argparse.Namespace) -> None:
    """Entry point function to generate atomic info files.
    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """
    # Get and set some paths
    out_dir = os.path.abspath(args.src_dir)
    
    _write_pt_atomdm(out_dir)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description=
        "This script is used to create the experimental data look up tables "
        "for the atom class."
    )

    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    main(args)
