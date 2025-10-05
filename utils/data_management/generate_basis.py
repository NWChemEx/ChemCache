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
that will fill them in.

Usage
-----

::

   usage: generate_basis.py [OPTIONS]... basis_set_source src_dir

   positional arguments:
   basis_set_source      Source directory for basis set files. If combined with
                         the "-r" flag, this directory will be recursively
                         searched for basis sets.
   src_dir               Destination directory for generated source files.

   options:
   -h, --help            show this help message and exit
   -r, --recursive       Toggle on recursive search through the basis set
                         source directory. Default OFF.
   -a ATOMS_DIR, --atoms_dir ATOMS_DIR
                         The path to where ElementNames.txt can be found.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

::

   <src_dir>
   ├── bases
   │   └── <all_basis_set_files>
   ├── basis_set_list.hpp
   └── load_basis_sets.cpp
"""

import argparse
import os
import re
import sys
from typing import Callable

import data_management.helper_fxns as helpers
from data_management.generate_atomicinfo import parse_symbols


class Shell:
    """Class representing a shell for an element."""

    def __init__(
        self, ls: list, num_format: str = ".10e", is_pure=True
    ) -> None:
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
        self.shelltype = "pure" if is_pure else "cartesian"

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

    def cxxify(self, center: str, tab: str = "    ") -> str:
        """Create a C++ source representation of the shell.

        :param center: chemist::Center to add the shell to
        :type center: str

        :param tab: String representing a tab, defaults to "    "
        :type tab: str, optional
        """

        add_shell = (
            "{t}{c}.emplace_back(make_shell(\n{t}  pure_t::{shelltype}, {l},"
        )

        lines = []
        for i in range(self.gen):
            l = self.ls[i]
            lines.append(
                add_shell.format(
                    t=tab * 3, c=center, l=l, shelltype=self.shelltype
                )
            )
            cs = "  doubles_t{"
            es = "  doubles_t{"
            for j, ai in enumerate(self.exp):
                ci = format(
                    float(
                        self.coefs[j][i].replace("D", "E").replace("E", "e")
                    ),
                    self.number_format,
                )
                ai_f = format(
                    float(ai.replace("D", "E").replace("E", "e")),
                    self.number_format,
                )
                cs += ci
                es += ai_f
                if j < len(self.exp) - 1:
                    cs += ", "
                    es += ", "
                else:
                    cs += "}"
                    es += "}"
            lines.append("{}{},".format(tab * 3, cs))
            lines.append("{}{}));".format(tab * 3, es))
        return "\n".join(lines)


def _write_basis_files(
    src_dir: str, bs_name: str, basis_set: dict, tab: str = "    "
) -> None:
    source_template = """
#include "../bases.hpp"
#include "{s_name}.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using atomic_basis_pt = simde::AtomicBasisSetFromZ;
using abs_t           = simde::type::atomic_basis_set;
using shell_t         = simde::type::shell;
using center_t        = simde::type::point;
using shells_t        = std::vector<shell_t>;
using doubles_t       = std::vector<double>;
using pure_t          = chemist::ShellType;

static constexpr auto module_desc = R"(
{d_name} atomic basis set
---------------------------------

This module returns the atomic basis set associated with an atomic number.
This module was autogenerated.
)";

MODULE_CTOR({s_name}_atom_basis) {{
    description(module_desc);
    satisfies_property_type<atomic_basis_pt>();
}}

MODULE_RUN({s_name}_atom_basis) {{
    const auto& [Z] = atomic_basis_pt::unwrap_inputs(inputs);
    auto rv         = results();
    auto z_string = std::to_string(Z);

    switch(Z) {{
{cases}
        default: {{
            throw std::out_of_range("Basis Set not available for Z: " + z_string);
        }}
    }}
}}

}} // namespace chemcache
"""

    cases_template = """{t}{t}case({Z}): {{
{t}{t}{t}return atomic_basis_pt::wrap_results(rv, {s_name}_{Z}());
{t}{t}}}"""

    s_name = helpers.sanitize_basis_name(bs_name)
    d_name = helpers.desanitize_basis_name(bs_name)

    basis_dir = os.path.join(src_dir, s_name)
    if not os.path.exists(basis_dir):
        os.mkdir(basis_dir)

    cases = []
    headers = []
    for z in sorted([int(x) for x in basis_set.keys()]):
        shells = []
        for shell in basis_set[str(z)]:
            shells.append(shell.cxxify("shells", tab))
        headers.append(
            "simde::type::atomic_basis_set {s_name}_{Z}();".format(
                s_name=s_name, Z=z
            )
        )
        _write_atomic_basis(
            basis_dir, tab, d_name, s_name, z, shells="\n".join(shells)
        )
        cases.append(cases_template.format(t=tab, Z=z, s_name=s_name))

    out_file = os.path.join(basis_dir, s_name + ".cpp")
    with open(out_file, "w") as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(
            source_template.format(
                d_name=d_name, s_name=s_name, cases="\n".join(cases)
            )
        )

    with open(os.path.join(basis_dir, s_name + ".hpp"), "w") as fout:
        fout.write("#pragma once\n")
        fout.write("#include <simde/types.hpp>\n")
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write("namespace chemcache {\n")
        fout.write("\n".join(headers))
        fout.write("\n}\n")


def _write_atomic_basis(
    src_dir: str, tab: str, d_name: str, s_name: str, z: str, shells: str
) -> None:
    source_template = """
#include "{basis_name}.hpp"
#include <simde/basis_set/atomic_basis_set.hpp>
#include <simde/types.hpp>

namespace chemcache {{

using abs_t     = simde::type::atomic_basis_set;
using shell_t   = simde::type::shell;
using center_t  = simde::type::point;
using shells_t  = std::vector<shell_t>;
using doubles_t = std::vector<double>;
using pure_t    = chemist::ShellType;

abs_t {basis_name}_{Z}(){{

{t}// Basis Set name and origin point
{t}std::string name("{d_name}");
{t}center_t r0(0.0, 0.0, 0.0);

{t}auto make_shell = [&r0](auto pure, auto l, const doubles_t& cs,
                            const doubles_t& es) {{
{t}{t}return shell_t(pure, l, cs.begin(), cs.end(), es.begin(), es.end(), r0);
{t}}};

{t}shells_t shells;
{shells}
{t} return abs_t(name, {Z}, r0, shells.begin(), shells.end());
}} // {basis_name}_{Z}

}} // chemcache
"""
    out_file = os.path.join(src_dir, s_name + "_" + str(z) + ".cpp")
    with open(out_file, "w") as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(
            source_template.format(
                Z=z, d_name=d_name, basis_name=s_name, t=tab, shells=shells
            )
        )


def _write_bases(src_dir: str, bases: dict, tab="    ") -> None:
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
DECLARE_MODULE(molecular_basis);
{d}

namespace bases_mods {{

inline void set_defaults(pluginplay::ModuleManager& mm) {{
    // Default submodules within this subcollection will be set here
    {s}
}}

inline void load_modules(pluginplay::ModuleManager& mm) {{
    // Modules will be added to the ModuleManager here
    {m}
    {b}

    set_defaults(mm);
}}

}} // namespace bases_mods

}} // namespace chemcache
"""

    d_template = "DECLARE_MODULE({}_atom_basis);"
    s_template = 'mm.change_submod("{}", "Atomic Basis", "{} atomic basis");'
    m_template = 'mm.add_module<molecular_basis>("{}");'
    a_template = 'mm.add_module<{}_atom_basis>("{} atomic basis");'

    ntab = "\n" + tab

    ds = []
    ss = []
    ms = []
    bs = []

    for bs_name, basis_set in sorted(bases.items()):
        _write_basis_files(src_dir, bs_name, basis_set)

        s_name = helpers.sanitize_basis_name(bs_name)
        d_name = helpers.desanitize_basis_name(bs_name)
        ds.append(d_template.format(s_name))
        ss.append(s_template.format(d_name, d_name))
        ms.append(m_template.format(d_name))
        bs.append(a_template.format(s_name, d_name))

    bases_file = os.path.join(src_dir, "bases.hpp")
    with open(bases_file, "w") as fout:
        helpers.write_warning(fout, os.path.basename(__file__))
        fout.write(
            bases_template.format(
                d="\n".join(ds),
                s=ntab.join(ss),
                m=ntab.join(ms),
                b=ntab.join(bs),
            )
        )


def _parse_bases_gbs(
    basis_set_filenames: list, sym2Z: dict, l2num: Callable[[str], int]
) -> dict:
    """Parses basis set files from the filepaths given.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list

    :param sym2Z: Dictionary associating atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Function associating orbital letters with a number
    :type l2num: Callable[[str], int]

    :return: Collection of basis sets and the supported elements of each.
    :rtype: dict
    """

    new_atom = re.compile(r"^\s*\D{1,2}\s*0\s*$")
    new_shell = re.compile(r"^\s*[a-zA-Z]+\s*\d+\s*1.00\s*$")
    same_shell = re.compile(r"^\s*(?:-?\d+.\d+(?:(E|e)(\+|-)\d\d)*\s*)+")
    bases = {}
    for filepath in basis_set_filenames:
        # Extract file name without extension
        bs = os.path.splitext(os.path.basename(filepath))[0]

        # Let the user know which basis set is being parsed
        print("Parsing {}...".format(os.path.basename(filepath)), end="")
        # Print immediately
        sys.stdout.flush()

        bases[bs] = {}
        with open(filepath, "r") as f:
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


def _parse_bases_nw(
    basis_set_filepaths: list, sym2Z: dict, l2num: Callable[[str], int]
) -> dict:
    """Parses basis set files from the filepaths given.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list

    :param sym2Z: Dictionary associating atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Function associating orbital letters with a number
    :type l2num: Callable[[str], int]

    :return: Collection of basis sets and the supported elements of each.
    :rtype: dict
    """

    atom_shell_line = re.compile(r"^[a-zA-Z]{1,3}\s+[a-zA-Z]+\s*$")
    shell_data = re.compile(
        r"^\s*(?:-?\d*\.\d+(?:(?:E|e|D|d)(?:\+|-)\d+)*\s*)+$"
    )
    basis_start = re.compile(r'^BASIS "ao basis" (?P<shelltype>\w+) PRINT')
    block_end = re.compile(r"^END$")

    bases = {}
    for filepath in basis_set_filepaths:
        # Extract file name without extension
        basis_set = os.path.splitext(os.path.basename(filepath))[0]

        # Let the user know which basis set is being parsed
        print("Parsing {}...".format(os.path.basename(filepath)), end="")
        # Print immediately
        sys.stdout.flush()

        # Read data from the file into memory
        tmp_basis = {}
        is_pure = True
        with open(filepath, "r") as fin:
            atom_z = 0  # Atomic number of current element
            ls = ""  # Angular momenta of current shells

            basis_block = False

            for line in fin:
                if basis_block:
                    if re.search(atom_shell_line, line):
                        atom_sym, ls = line.split()

                        # Update current atomic number and angular momentum
                        atom_z = sym2Z[atom_sym.lower()]

                        # Add the new element to basis set list
                        if atom_z not in tmp_basis.keys():
                            tmp_basis[atom_z] = {}
                        if ls not in tmp_basis[atom_z].keys():
                            tmp_basis[atom_z][ls] = []

                        tmp_basis[atom_z][ls].append([])

                    elif re.search(shell_data, line):
                        # Build up shell data to generate primitives
                        tmp_basis[atom_z][ls][-1].append(line.split())

                    elif re.search(block_end, line):
                        # Done reading basis; ignore the rest
                        break

                # Check for basis block start
                elif bstart_match := re.search(basis_start, line):
                    # Indicate that we are in a basis block and should parse
                    basis_block = True
                    is_pure = bstart_match["shelltype"] == "SPHERICAL"

        # Did we get anything out of the file?
        if not tmp_basis:
            print("no basis information parsed")
            sys.stdout.flush()
            continue

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
                        shell = Shell(
                            [l2num(l.lower()) for l in ls], is_pure=is_pure
                        )

                        for i in range(len(exp)):
                            coef = [row[i] for row in coefs[j : j + rowspan :]]

                            # No coefficients should be zero
                            if all([float(c) != 0.0 for c in coef]):
                                shell.add_prim(exp[i], coef)

                        bases[basis_set][element].append(shell)

        # Let the user know the basis set is done being parsed
        print("complete")
        # Print immediately
        sys.stdout.flush()

    return bases


def _parse_bases(
    basis_set_filepaths: list,
    sym2Z: dict,
    l2num: Callable[[str], int],
    format: str = "nwchem",
) -> dict:
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
    :type l2num: Callable[[str], int]

    :param format: File formatting to parse, defaults to "nwchem"
    :type format: str, optional

    :return: Basis sets parsed
    :rtype: dict

    :raises RuntimeError: Unsupported basis file format.
    """

    if format == "gaussian94" or format == "psi4" or format == "xtron":
        return _parse_bases_gbs(basis_set_filepaths, sym2Z, l2num)
    elif format == "nwchem":
        return _parse_bases_nw(basis_set_filepaths, sym2Z, l2num)
    else:
        raise RuntimeError("Unsupported basis file format: {}".format(format))


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate basis set files.

    :param args: Command line argument namespace
    :type args: argparse.Namespace
    """

    formats = ["nwchem"]

    # Create some paths
    src_dir = os.path.abspath(args.src_dir)
    name_file = os.path.abspath(
        os.path.join(args.atoms_dir, "ElementNames.txt")
    )

    # Discover basis set files
    basis_set_dir = os.path.abspath(args.basis_set_source)
    basis_set_filepaths = helpers.find_files(
        basis_set_dir,
        [helpers.lookup_extension(format) for format in formats],
        recursive=args.recursive,
    )

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = {ai.sym.lower(): ai.Z for ai in atoms.values()}

    def l2num(l: str) -> int:
        return "spdfghijklmnoqrtuvwxyzabce".find(l.lower())

    basis_sets = {}
    for format in formats:
        extension = helpers.lookup_extension(format)

        # NOTE: Format order CAN matter!
        #       If the same basis set exists in basis_sets
        #       and the new dict returned from parse_bases(), the
        #       basis_sets version will be replaced by the
        #       parse_bases() version.
        basis_sets.update(
            _parse_bases(
                basis_set_filepaths[extension], sym2Z, l2num, format=format
            )
        )

    print("Writing basis sets to file...", end="")
    sys.stdout.flush()

    _write_bases(src_dir, basis_sets)

    print("complete")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: argparse.Namespace
    """

    parser = argparse.ArgumentParser(
        description=(
            "This script will loop over a series of basis sets and write out "
            "a file that will fill them in. The format of the resulting "
            "basis sets is suitable for use with the BasisSetExchange class."
        )
    )

    parser.add_argument(
        "basis_set_source",
        type=str,
        help="""Source directory for basis set files. If combined
                             with the \"-r\" flag, this directory will be
                             recursively searched for basis sets.""",
    )

    parser.add_argument(
        "src_dir",
        type=str,
        help="Destination directory for generated source files.",
    )

    parser.add_argument(
        "-r",
        "--recursive",
        action="store_true",
        help="""Toggle on recursive search through the basis
                             set source directory. Default OFF.""",
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
