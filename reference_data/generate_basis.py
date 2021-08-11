#!/usr/bin/env python3

"""This script will loop over a series of basis sets and write out a file 
that will fill them in.  The format of the resulting basis sets is suitable 
for use with the BasisSetExchange class.

This script creates the following files based on the include and source
directories given. The directories are not created by this script and must
be present before running it.

+---include
|       nwx_basis_list.hpp
|
+---src
|   \---bases
|           <all_basis_set_files>
|       nwx_basis_set_manager_pimpl.cpp
"""

import argparse
import os
import re
import sys

from generate_atomicinfo import parse_symbols
import helper_fxns as helpers

class Shell:
    """Class representing a shell for an element.
    """

    def __init__(self, ls, num_format=".10e"):
        """Initialization function.

        :param ls: List of angular momenta
        :type ls: list of int
        """

        self.ls = ls
        self.exp=[]
        self.coefs=[]
        self.gen = 0
        self.number_format = num_format

    def add_prim(self, exp, coefs):
        """Add a primitive for the shell.

        :param exp: Primitive exponent
        :type exp: str
        :param coefs: Primitive contraction coefficients
        :type coefs: list of str
        """

        self.exp.append(exp)
        self.coefs.append(coefs)
        self.gen = max(len(coefs), self.gen)

    def cxxify(self, tab, fout):
        """Create a C++ source representation of the shell.

        :param tab: String representing a tab.
        :type tab: str

        :param fout: C++ source file opened for writing
        :type fout: class: _io.TextIOWrapper
        """

        for i in range(self.gen):
            l = self.ls[i]
            fout.write("{}rv.add_shell(ShellType::pure, {},\n".format(tab, l))
            cs = "std::vector<double>{"
            es = "std::vector<double>{"
            for j,ai in enumerate(self.exp):
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
            fout.write("{}    {},\n".format(tab, cs))
            fout.write("{}    {});\n".format(tab, es))

def print_pimpl_header(f):
    helpers.write_warning(f, os.path.basename(__file__))

    f.write(
"""#include <libchemist/managers/detail_/basis_set_manager_pimpl.hpp>
#include \"chemcache/nwx_basis_list.hpp\"
        
namespace chemcache::detail_ {
       
class HardCodedBSMan : public BasisSetManagerPIMPL {
public:
    HardCodedBSMan() = default;
    using ao_basis_type = typename BasisSetManagerPIMPL::ao_basis_type;
protected:
    HardCodedBSMan(const HardCodedBSMan& rhs) = default;
private:
    std::unique_ptr<BasisSetManagerPIMPL> clone_() const override {
        return std::unique_ptr<HardCodedBSMan>(new HardCodedBSMan(*this));
    }
            
    ao_basis_type get_basis_(const std::string& name, 
                             size_type Z) const override {         
""")

def print_pimpl_footer(f):
    f.write(
"""throw std::out_of_range(\"Unrecognized basis name\");
    }//end get_basis_
};
        
std::unique_ptr<BasisSetManagerPIMPL> nwx_default_bs() {
    return std::make_unique<HardCodedBSMan>();
}
        
} // namespace chemcache::detail_
""")

def print_basis_header(f, bs_name):
    helpers.write_warning(f, os.path.basename(__file__))

    f.write(
"""#include \"chemcache/nwx_basis_list.hpp\"
 
namespace chemcache::detail_ {{
 
Center<double> {}(std::size_t Z) {{
    switch(Z) {{         
""".format(bs_name))

def print_basis_list(f):
    helpers.write_warning(f, os.path.basename(__file__))

    f.write(
"""#include <libchemist/basis_set/basis_set.hpp>

namespace chemcache::detail_ {

""")

def print_basis_footer(f):
    tab = "    "
    f.write(
"""{}default : {{ 
{}throw std::out_of_range(\"Basis not available for Z\");
{}}}\n{}}} // end switch\n
}} //end function
}} //end namespace chemcache::detail_""".format(tab*2, tab*3, tab*2, tab))

def print_atom_basis(f, z, atom):
    tab = "    "
    f.write("{}case({}) : {{\n".format(tab*2, z))
    f.write("{}Center<double> rv(0.0, 0.0, 0.0);\n".format(tab*3))
    for s in atom:
        s.cxxify(tab*3,f)
    f.write("{}return rv;\n".format(tab*3))
    f.write("{}}} //End case\n".format(tab*2))

def write_bases(inc_dir, src_dir, bases):
    """Writes basis set data to C++ files.

    :param inc_dir: Include directory for header files.
    :type inc_dir: str

    :param src_dir: Source directory for source files.
    :type src_dir: str

    :param bases: Collection of basis sets parsed from files
    :type bases: dict
    """

    tab = "    "
    with open(os.path.join(src_dir,"nwx_basis_set_manager_pimpl.cpp"),'w') as f:
        print_pimpl_header(f)
        with open(os.path.join(inc_dir, "nwx_basis_list.hpp"), 'w') as g:
            print_basis_list(g)
            f.write("{}".format(tab*2))
            for bs_name, bs in sorted(bases.items()):
                s_name = helpers.sanitize_basis_name(bs_name)
                d_name = helpers.desanitize_basis_name(bs_name)
                f.write("if(name == \"{}\") {{ ".format(d_name))
                f.write("return {}(Z); ".format(s_name))
                g.write("Center<double> {}(std::size_t Z);\n".format(s_name))
                bs_file_name = "{}.cpp".format(bs_name)
                bs_path = os.path.join(src_dir,"bases", bs_file_name)
                with open(bs_path, 'w') as h:
                    print_basis_header(h, s_name)
                    for z in sorted([int(x) for x in bs.keys()]):
                        print_atom_basis(h, z, bs[str(z)])
                    print_basis_footer(h)
                f.write("}}\n{}else ".format(tab*2))
            g.write("\n} // namespace chemcache::basis_set\n")
        print_pimpl_footer(f)

def parse_bases_gbs(basis_set_filenames, sym2Z, l2num):
    """Parses basis set files from the filepaths given.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list

    :param sym2Z: Dictionary associating atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Lambda function associating orbital letters with a number
    :type l2num: lambda

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
        with open(filepath,'r') as f:
            atom_z = 0
            for line in f:
                if re.search(new_atom, line):
                    atom_z = sym2Z[line.split()[0].lower()]
                    bases[bs][atom_z] = []
                elif re.search(new_shell, line):
                    ls = [ l2num(l.lower()) for l in line.split()[0]]
                    bases[bs][atom_z].append(Shell(ls))
                elif re.search(same_shell, line):
                    prim = line.split()
                    bases[bs][atom_z][-1].add_prim(prim[0], prim[1:])

        # Let the user know which basis set is being parsed
        print("complete")
        # Print immediately
        sys.stdout.flush()

    return bases

def parse_bases_nw(basis_set_filepaths, sym2Z, l2num):
    """Parses basis set files from the filepaths given.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list

    :param sym2Z: Dictionary associating atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Lambda function associating orbital letters with a number
    :type l2num: lambda

    :return: Collection of basis sets and the supported elements of each.
    :rtype: dict
    """

    atom_shell_line = re.compile("^[a-zA-Z]{1,3}\s+[a-zA-Z]+\s*$")
    shell_data = re.compile(
        "^\s*(?:-?\d+\.\d+(?:(?:E|e|D|d)(?:\+|-)\d\d)*\s*)+$")

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

            for line in fin:
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

                    coefs = [ row for row in exp_coefs_t[1:] ]

                    for j in range(len(coefs) - rowspan + 1):
                        shell = Shell([ l2num(l.lower()) for l in ls ])

                        for i in range(len(exp)):
                            coef = [ row[i] for row in coefs[j:j + rowspan:] ]

                            # No coefficients should be zero
                            if (all([ float(c) != 0.0 for c in coef ])):
                                shell.add_prim(exp[i], coef)
                        
                        bases[basis_set][element].append(shell)

        # Let the user know which basis set is being parsed
        print("complete")
        # Print immediately
        sys.stdout.flush()

    return bases

def parse_bases(basis_set_filepaths, sym2Z, l2num, format="nwchem"):
    """Parse basis set files of the specified format.

    :param basis_set_filepaths: Full paths to basis set files.
    :type basis_set_filepaths: list of str

    :param sym2Z: Mapping from lowercased atomic symbols to atomic numbers
    :type sym2Z: dict

    :param l2num: Conversion function from shell symbol (s, p, d, f, etc)
        to the corresponding number (0, 1, 2, 3, etc)
    :type l2num: lambda

    :param format: File formatting to parse, defaults to "nwchem"
    :type format: str, optional

    :raises RuntimeError: Unsupported basis file format.
    """

    if (format == "gaussian94" or 
        format == "psi4"       or
        format == "xtron"):
        
        return parse_bases_gbs(basis_set_filepaths, sym2Z, l2num)
    elif (format == "nwchem"):
        return parse_bases_nw(basis_set_filepaths, sym2Z, l2num)
    else:
        raise RuntimeError("Unsupported basis file format: {}".format(format))

def find_basis_sets(source_root, formats=["nwchem"], recursive=False):
    """Recursively find all basis set files in the given directory. Basis set
    files are identified using the given extensions list.

    :param source_root: Root directory containing basis set files.
    :type source_root: str
    
    :param extensions: Possible extensions for basis sets, defaults to [".nw"]
    :type extensions: list, optional

    :param recursive: Whether or not to recursively search in subdirectories,
        defaults to False
    :type recursive: bool
    
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

def main(args):
    """Entry point function to generate basis set files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    formats = [ "nwchem" ]

    # Set include directory to default src_dir path if no path is specified
    if (args.inc == ""):
        args.inc = args.src_dir

    # Create some paths
    my_dir    = os.path.dirname(os.path.realpath(__file__))
    inc_dir   = os.path.abspath(args.inc)
    src_dir   = os.path.abspath(args.src_dir)
    test_dir  = os.path.abspath(args.test_path)
    name_file = os.path.join(my_dir, "physical_data", "ElementNames.txt")

    # Discover basis set files
    basis_set_dir       = os.path.abspath(args.basis_set_source)
    basis_set_filepaths = helpers.find_files(
        basis_set_dir, 
        [ helpers.lookup_extension(format) for format in formats ],
        recursive=args.recursive
    )

    # Parse element information
    atoms = {}
    parse_symbols(name_file, atoms)

    sym2Z = { ai.sym.lower() : ai.Z for ai in atoms.values() }
    l2num = lambda l: "spdfghijklmnoqrtuvwxyzabce".find(l.lower())

    basis_sets = {}
    for format in formats:
        extension = helpers.lookup_extension(format)

        # NOTE: Format order CAN matter!
        #       If the same basis set exists in basis_sets
        #       and the new dict returned from parse_bases(), the 
        #       basis_sets version will be replaced by the 
        #       parse_bases() version.
        basis_sets.update(parse_bases(
            basis_set_filepaths[extension], sym2Z, l2num, format=format
        ))

    write_bases(inc_dir, src_dir, basis_sets)

def parse_args():
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """
    
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('basis_set_source', type=str,
                        help="""Source directory for basis set files. If combined
                             with the \"-r\" flag, this directory will be
                             recursively searched for basis sets.""")
    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")
    parser.add_argument('test_path', type=str,
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
