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


"""This script will loop over a series of basis sets and write out a file
that will fill them in. The format of the resulting basis sets is suitable
for use with the BasisSetExchange class.

Usage
-----

::

   usage: generate_basis.py [-h] [-r] basis_set_source src_dir

   positional arguments:
     basis_set_source  Source directory for basis set files. If combined with the "-r" flag, this directory will be recursively searched for basis sets.
     src_dir           Destination directory for generated source files.

   optional arguments:
     -h, --help        show this help message and exit
     -r, --recursive   Toggle on recursive search through the basis set source directory. Default OFF.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

::

   +---src_dir
   |   \---bases
   |           <all_basis_set_files>
   |       basis_set_list.hpp
   |       load_basis_sets.cpp
"""

import argparse
import io
import os
import re
import sys

from generate_atomicinfo import AtomicData, parse_symbols
import helper_fxns as helpers


class Shell:
    """Class representing a shell for an element.
    """

    def __init__(self, ls: list, num_format: str = ".10e") -> None:
        """Initialization function.

        :param ls: List of angular momenta
        :type ls: list of int

        :param num_format: Format string for numbers being output, 
            defaults to ".10e"
        :type num_format: str, optional
        """

        self.ls = ls
        self.exp = []
        self.coefs = []
        self.gen = 0
        self.number_format = num_format

    def add_prim(self, exp: str, coefs: list) -> None:
        """Add a primitive for the shell.

        :param exp: Primitive exponent
        :type exp: str

        :param coefs: Primitive contraction coefficients
        :type coefs: list of str
        """

        self.exp.append(exp)
        self.coefs.append(coefs)
        self.gen = max(len(coefs), self.gen)

    def cxxify(self, fout: io.TextIOWrapper, center: str,
               tab: str = "    ") -> None:
        """Create a C++ source representation of the shell.

        :param fout: C++ source file opened for writing
        :type fout: class: _io.TextIOWrapper

        :param center: chemist::Center to add the shell to
        :type center: str

        :param tab: String representing a tab, defaults to "    "
        :type tab: str, optional
        """

        for i in range(self.gen):
            l = self.ls[i]
            fout.write(
                "{}{}.add_shell(chemist::ShellType::pure, {},\n".format(tab, center, l))
            cs = "std::vector<double>{"
            es = "std::vector<double>{"
            for j, ai in enumerate(self.exp):
                ci = format(float(self.coefs[j][i].replace('D', 'E')
                                  .replace('E', 'e')), self.number_format)
                ai_f = format(float(ai.replace('D', 'E')
                                    .replace('E', 'e')), self.number_format)
                cs += ci
                es += ai_f
                if j < len(self.exp) - 1:
                    cs += ','
                    es += ','
                else:
                    cs += '}'
                    es += '}'
            fout.write("{}{},\n".format(tab*2, cs))
            fout.write("{}{});\n".format(tab*2, es))


def _print_atom_basis(fout: io.TextIOWrapper, bs_name: str, z: int,
                     basis: list, tab: str = "    ") -> None:
    """Print basis data for the given atom.

    :param fout: C++ source file opened for writing
    :type fout: io.TextIOWrapper

    :param bs_name: Name of the basis set that this basis belongs to
    :type bs_name: str

    :param z: Atomic number of the element that the basis is assigned to
    :type z: int

    :param basis: List of shells representing the basis of an atom
    :type basis: list

    :param tab: String representing a tab, defaults to "    "
    :type tab: str, optional
    """

    d_name = helpers.desanitize_basis_name(bs_name)

    fout.write("{}basis_map.emplace({}, ".format(tab, z))
    fout.write("chemist::AtomicBasisSet<double>(\"{}\", {}, 0.0, 0.0, 0.0));\n\n".format(d_name, z))

    center = "basis_map.at({})".format(z)

    for shell in basis:
        shell.cxxify(fout, center, tab)


def _write_basis_files(out_file: str, bs_name: str, basis_set: dict,
                      tab: str = "    ") -> None:

    with open(out_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        s_name = helpers.sanitize_basis_name(bs_name)

        # Start of the file
        fout.write(
            """#include "../basis_set_list.hpp"
# include <chemist/chemist.hpp>

namespace chemcache::basis_sets {{

void load_{}(chemist::BasisSetManager& bsm) {{
{}chemist::BasisSetManager::ao_basis_map basis_map;

""".format(s_name, tab)
        )

        # Write each atomic basis in a basis set
        for z in sorted([int(x) for x in basis_set.keys()]):
            # Comment basis set name and atomic number being handled
            fout.write("{}// Z = {}\n".format(tab, z))

            _print_atom_basis(fout, bs_name, z, basis_set[str(z)], tab)

            # White space between bases
            fout.write("\n")

        # White space between basis sets
        fout.write("\n")

        # End the file
        fout.write("""{}bsm.insert(\"{}\", basis_map);
}} // function load_{}

}} // namespace chemcache::basis_sets
""".format(tab, bs_name, s_name)
        )


def _write_basis_list(src_dir: str, bases: dict, tab="    ") -> None:
    # Write out the basis list file
    basis_list_file = os.path.join(src_dir, "basis_set_list.hpp")
    with open(basis_list_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file the load_basis_file
        fout.write(
            """#include <chemist/chemist.hpp>

namespace chemcache::basis_sets {

"""
        )

        # Create function declarations for load_<basis_name> functions
        for bs_name, _ in sorted(bases.items()):
            # Make call to load basis set into bsm
            fout.write(
                "void load_{}(chemist::BasisSetManager& bsm);\n"
                .format(helpers.sanitize_basis_name(bs_name)))

        # End of the file
        fout.write("\n} // namespace chemcache")


def _write_bases(src_dir: str, bases: dict, tab="    ") -> None:
    """Writes basis set data to C++ files.

    :param src_dir: Source directory for source files.
    :type src_dir: str

    :param bases: Collection of basis sets parsed from files
    :type bases: dict

    :param tab: String representing a tab, defaults to "    "
        :type tab: str, optional
    """

    load_basis_file = os.path.join(src_dir, "load_basis_sets.cpp")
    with open(load_basis_file, 'w') as fout:
        helpers.write_warning(fout, os.path.basename(__file__))

        # Start of the file the load_basis_file
        fout.write(
            """#include "chemcache/chemcache.hpp"
#include "./basis_set_list.hpp"
#include <chemist/chemist.hpp>

namespace chemcache {

void load_basis_sets(chemist::BasisSetManager& bsm) {
"""
        )

        # Add basis sets
        for bs_name, basis_set in sorted(bases.items()):
            # Make call to load basis set into bsm
            fout.write("{}basis_sets::load_{}(bsm);"
                       .format(tab, helpers.sanitize_basis_name(bs_name)))

            basis_file = os.path.join(src_dir, "bases", bs_name + ".cpp")
            _write_basis_files(basis_file, bs_name, basis_set)

            # White space between basis sets
            fout.write("\n")

        # End of the file
        fout.write(
            """} // function load_basis_sets

} // namespace chemcache
"""
        )

    _write_basis_list(src_dir, bases, tab)


def _parse_bases_gbs(basis_set_filenames: list, sym2Z: dict,
                    l2num: "function") -> dict:
    """Parses basis set files from the filepaths given.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list

    :param sym2Z: Dictionary associating atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Function associating orbital letters with a number
    :type l2num: function

    :return: Collection of basis sets and the supported elements of each.
    :rtype: dict
    """

    new_atom = re.compile("^\s*\D{1,2}\s*0\s*$")
    new_shell = re.compile("^\s*[a-zA-Z]+\s*\d+\s*1.00\s*$")
    same_shell = re.compile("^\s*(?:-?\d+.\d+(?:(E|e)(\+|-)\d\d)*\s*)+")
    bases = {}
    for filepath in basis_set_filenames:
        # Extract file name without extension
        bs = os.path.splitext(os.path.basename(filepath))[0]

        # Let the user know which basis set is being parsed
        print("Parsing {}...".format(os.path.basename(filepath)),
              end='')
        # Print immediately
        sys.stdout.flush()

        bases[bs] = {}
        with open(filepath, 'r') as f:
            atom_z = 0
            for line in f:
                if re.search(new_atom, line):
                    atom_z = sym2Z[line.split()[0].lower()]
                    bases[bs][atom_z] = []
                elif re.search(new_shell, line):
                    ls = [l2num(l.lower()) for l in line.split()[0]]
                    bases[bs][atom_z].append(Shell(ls))
                elif re.search(same_shell, line):
                    prim = line.split()
                    bases[bs][atom_z][-1].add_prim(prim[0], prim[1:])

        # Let the user know which basis set is being parsed
        print("complete")
        # Print immediately
        sys.stdout.flush()

    return bases


def _parse_bases_nw(basis_set_filepaths: list, sym2Z: dict, l2num: "function") -> dict:
    """Parses basis set files from the filepaths given.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list

    :param sym2Z: Dictionary associating atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Function associating orbital letters with a number
    :type l2num: function

    :return: Collection of basis sets and the supported elements of each.
    :rtype: dict
    """

    atom_shell_line = re.compile("^[a-zA-Z]{1,3}\s+[a-zA-Z]+\s*$")
    shell_data = re.compile(
        "^\s*(?:-?\d*\.\d+(?:(?:E|e|D|d)(?:\+|-)\d+)*\s*)+$")
    basis_start = re.compile("^BASIS")
    block_end = re.compile("^END$")

    bases = {}
    for filepath in basis_set_filepaths:
        # Extract file name without extension
        basis_set = os.path.splitext(os.path.basename(filepath))[0]

        # Let the user know which basis set is being parsed
        print("Parsing {}...".format(os.path.basename(filepath)),
              end='')
        # Print immediately
        sys.stdout.flush()

        # Read data from the file into memory
        tmp_basis = {}
        with open(filepath, 'r') as fin:
            atom_z = 0     # Atomic number of current element
            ls = ""        # Angular momenta of current shells

            basis_block = False

            for line in fin:
                if basis_block:
                    if re.search(atom_shell_line, line):
                        atom_sym, ls = line.split()

                        # Update current atomic number and angular momentum
                        atom_z = sym2Z[atom_sym.lower()]

                        # Add the new element to basis set list
                        if (not atom_z in tmp_basis.keys()):
                            tmp_basis[atom_z] = {}
                        if (not ls in tmp_basis[atom_z].keys()):
                            tmp_basis[atom_z][ls] = []

                        tmp_basis[atom_z][ls].append([])

                    elif re.search(shell_data, line):
                        # Build up shell data to generate primitives
                        tmp_basis[atom_z][ls][-1].append(line.split())

                    elif re.search(block_end, line):
                        # Done reading basis; ignore the rest
                        break

                # Check for basis block start
                elif re.search(basis_start, line):
                    # Indicate that we are in a basis block and should parse
                    basis_block = True

        # Process the data read from the file
        bases[basis_set] = {}

        for element in tmp_basis.keys():
            bases[basis_set][element] = []

            for ls in tmp_basis[element].keys():
                for exp_coefs in tmp_basis[element][ls]:
                    # Transpose the array
                    exp_coefs_t = list(map(list, zip(*exp_coefs)))

                    rowspan = len(ls)

                    exp = exp_coefs_t[0]

                    coefs = [row for row in exp_coefs_t[1:]]

                    for j in range(len(coefs) - rowspan + 1):
                        shell = Shell([l2num(l.lower()) for l in ls])

                        for i in range(len(exp)):
                            coef = [row[i] for row in coefs[j:j + rowspan:]]

                            # No coefficients should be zero
                            if (all([float(c) != 0.0 for c in coef])):
                                shell.add_prim(exp[i], coef)

                        bases[basis_set][element].append(shell)

        # Let the user know the basis set is done being parsed
        print("complete")
        # Print immediately
        sys.stdout.flush()

    return bases


def _parse_bases(basis_set_filepaths: list, sym2Z: dict, l2num: "function", format: str = "nwchem") -> dict:
    """Parse basis set files of the specified format.

    This function redirects to the correct parsing function based on the file
    format given. The data structure returned is a dict with basis set
    names as keys, effectively making a mapping from basis set names to
    basis sets. The basis set values are another dict with atomic
    numbers as keys and a list of Shell as the values. This makes a structure,
    `basis_sets`, where the basis for `Z` in `basis_name` can be accessed as:

    .. code-block:: Python

       list_of_shells = basis_sets[basis_name][Z]

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Conversion function from shell symbol (s, p, d, f, etc)
        to the corresponding number (0, 1, 2, 3, etc)
    :type l2num: function

    :param format: File formatting to parse, defaults to "nwchem"
    :type format: str, optional

    :return: Basis sets parsed
    :rtype: dict

    :raises RuntimeError: Unsupported basis file format.
    """

    if (format == "gaussian94" or
        format == "psi4" or
            format == "xtron"):

        return _parse_bases_gbs(basis_set_filepaths, sym2Z, l2num)
    elif (format == "nwchem"):
        return _parse_bases_nw(basis_set_filepaths, sym2Z, l2num)
    else:
        raise RuntimeError("Unsupported basis file format: {}".format(format))


def _find_basis_sets(source_root: str, formats: list = ["nwchem"], recursive: bool = False) -> list:
    """Recursively find all basis set files in the given directory. Basis set
    files are identified using the given extensions list.

    :param source_root: Root directory containing basis set files.
    :type source_root: str

    :param extensions: Possible extensions for basis sets, defaults to ["nwchem"]
    :type extensions: list, optional

    :param recursive: Whether or not to recursively search in subdirectories,
        defaults to False
    :type recursive: bool, optional

    :return: Full paths to basis set files found
    :rtype: list
    """
    basis_sets = []

    for format in formats:
        extension = helpers.lookup_extension(format)

        # Recursively search for basis set files of the given extension
        for dirpath, _, filenames in os.walk(source_root):
            for filename in filenames:
                _, ext = os.path.splitext(filename)

                # Case insensitive comparison of extensions
                if (ext.lower() == extension.lower()):
                    basis_sets.append(os.path.join(dirpath, filename))

            # Stop the recursive search before diving into subdirectories
            # if recursion is not requested
            if (not recursive):
                break

    return basis_sets


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate basis set files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    formats = ["nwchem"]

    # Create some paths
    my_dir = os.path.dirname(os.path.realpath(__file__))
    src_dir = os.path.abspath(args.src_dir)
    name_file = os.path.join(my_dir, "physical_data", "ElementNames.txt")

    # Discover basis set files
    basis_set_dir = os.path.abspath(args.basis_set_source)
    basis_set_filepaths = helpers.find_files(
        basis_set_dir,
        [helpers.lookup_extension(format) for format in formats],
        recursive=args.recursive
    )

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = {ai.sym.lower(): ai.Z for ai in atoms.values()}
    def l2num(l): return "spdfghijklmnoqrtuvwxyzabce".find(l.lower())

    basis_sets = {}
    for format in formats:
        extension = helpers.lookup_extension(format)

        # NOTE: Format order CAN matter!
        #       If the same basis set exists in basis_sets
        #       and the new dict returned from parse_bases(), the
        #       basis_sets version will be replaced by the
        #       parse_bases() version.
        basis_sets.update(_parse_bases(
            basis_set_filepaths[extension], sym2Z, l2num, format=format
        ))

    print("Writing basis sets to file...", end='')
    sys.stdout.flush()

    _write_bases(src_dir, basis_sets)

    print("complete")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser(description=
        "This script will loop over a series of basis sets and write out a "
        "file that will fill them in. The format of the resulting basis sets "
        "is suitable for use with the BasisSetExchange class."
    )

    parser.add_argument('basis_set_source', type=str,
                        help="""Source directory for basis set files. If combined
                             with the \"-r\" flag, this directory will be
                             recursively searched for basis sets.""")

    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")

    parser.add_argument('-r', '--recursive', action="store_true",
                        help="""Toggle on recursive search through the basis
                             set source directory. Default OFF.""")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    main(args)
