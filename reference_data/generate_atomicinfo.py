#!/usr/bin/env python3

"""This script is used to create the experimental data look up tables for the
atom class.

Original author: Ben Pritchard

In order to run this script simply needs to know where you want the generated
files to live.

For readability and convenience we use a few abbreviations throughout this
script:

- Z the atomic number of an atom
- Sym the atomic symbol of an atom (e.g. H for hydrogen, He for helium)
 
This script creates the following files:

chemcache_root
|
+---src
|       nwx_periodic_table_pimpl.cpp
|
+---tests
|       periodic_table.cpp
"""

import argparse
import os
import re

class AtomicData:
    def __init__(self):
        self.sym = ""
        self.name = ""
        self.Z = 0
        self.mass = 0.0
        self.isotopes = []       # List of isotope mass numbers
        self.isotope_masses = {} # Isotope mass values, indexed by mass number

    def add_isotope(self, num, mass):
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
            rv +="{}: {}".format(x, self.isotope_masses[x])
            if x!= self.isotopes[-1]:
                rv+=", "
        
        rv += "]"
        
        return rv


def parse_symbols(name_file, atoms):
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

def parse_ciaaw_isotopes(iso_file, atoms):
    """Parses an isotope mass file from the Commission on Isotopic Abundances
    and Atomic Weights (CIAAW) and adds it to a given atomic collection.

    :param iso_file: CIAAW isotope mass file
    :type iso_file: str
    :param atoms: Collection of atoms. Loaded isotopes will be added here.
    :type atoms: dict of class:`AtomicData`
    """

    new_atom =  r"(\d+)\s+[a-zA-Z]{1,2}\s+[a-zA-Z]+\s+"
    new_iso = r"(\d+)\**\s+(\d+\.\d+)\s((\d+\s?)+)+"
    new_atom += new_iso

    def par_fxn(match, atom):
        """Parse the isotope match and add it to the given atom.

        :param match: Regex match groups for an isotope
        :type match: tuple
        :param atom: Atom to add the isotope to
        :type atom: class:`AtomicData`
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

def parse_ciaww_mass(mass_file, atoms):
    """Parses a mass file from the Commission on Isotopic Abundances
    and Atomic Weights (CIAAW) and adds it to a given atomic collection.

    :param iso_file: CIAAW mass file
    :type iso_file: str
    :param atoms: Collection of atoms. Loaded masses will be added here.
    :type atoms: dict of class:`AtomicData`
    """

    mass = r"((?:\d+\.\d+(?:\s\d+)*,*\s?)+)"

    with open(mass_file, 'r') as fin:
        for line in fin:
            if re.search(mass, line):
                Z = line.split()[0]

                masses = re.search(mass, line).groups()[0].replace(" ","")

                # Split the mass string and convert to floats
                masses = [ float(mass) for mass in masses.split(',') ]

                # Average masses since it can be single mass or min/max from 
                # error bar
                atoms[Z].mass = sum(masses) / len(masses)

def write_ptable(out_dir, amu2me, atoms):
    """Generate periodic table PIMPL class with atom lookup tables.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of :class:`AtomicData`
    """

    out_file = os.path.join(out_dir, "nwx_periodic_table_pimpl.cpp")
    in_file = "libchemist/detail_/periodic_table_pimpl.hpp"
    sorted_keys = sorted([int(x) for x in atoms.keys()])
    max_Z = sorted_keys[-1]
    tab = "    "
    with open(out_file, 'w') as fout:
        fout.write(
"""/* 
 * This file has been autogenerated by generate_atomicinfo.py.  Any 
 * modifications made to this file will be lost next time generate_atomicinfo.py
 * is run.
 */

#include <algorithm> // For std::transform
#include <{}>
        
namespace chemcache::detail_ {{

class NWXPeriodicTablePIMPL : public libchemist::detail_::PeriodicTablePIMPL {{
public:
    NWXPeriodicTablePIMPL() = default;
protected:
    using my_type = NWXPeriodicTablePIMPL;
    NWXPeriodicTablePIMPL(const my_type& rhs) = default;
private:    
    std::unique_ptr<PeriodicTablePIMPL> clone_() const override {{
        return std::unique_ptr<my_type>(new my_type(*this));
    }}            
            
    size_type max_Z_() const noexcept override {{
        return {};
    }}
            
    Atom get_atom_(size_type Z) const override {{
        switch(Z) {{
""".format(in_file, max_Z))

        for k in sorted_keys:
            ai = atoms[str(k)]
            fout.write("{}case({}): {{ return Atom(".format(tab*3, k))
            fout.write("{}ul, {}, \"{}\");}}\n".format(k, ai.mass*amu2me, ai.sym))

        fout.write(
"""            default : {{ throw std::out_of_range(\" Z > {}\"); }}
        }}
    }}
            
    Atom get_isotope_(size_type Z, size_type mass_num) const override {{
        switch(Z) {{
""".format(max_Z))

        for k in sorted_keys:
            fout.write("{}case({}): {{\n".format(tab*3, k)) #Start case 1
            fout.write("{}switch(mass_num) {{\n".format(tab*4)) #Start switch 2
            ai = atoms[str(k)]

            for mn in sorted([int(x) for x in ai.isotopes]):
                mi = ai.isotope_masses[str(mn)] * amu2me
                fout.write("{}case({}): {{ return Atom(".format(tab*5, mn))
                fout.write("{}ul, {}, \"{}\");}}\n".format(k, mi, ai.sym))
            
            fout.write("{}default : {{ throw std::out_of_range(\"No isotope "
                    "data\"); }}\n".format(tab*5))
            fout.write("{}}}\n".format(tab*4)) #Close switch 2
            fout.write("{}}}\n".format(tab*3)) #Close case 1

        fout.write(
"""            default : {{ throw std::out_of_range(\"Z > {}\"); }}
        }}
    }}
    
    isotope_list isotopes_(size_type Z) const override {{
        switch(Z) {{
""".format(max_Z))

        for k in sorted_keys:
            fout.write("{}case({}) : {{ return {{".format(tab*3, k))
            ai = atoms[str(k)]

            for mn in sorted([int(x) for x in ai.isotopes]):
                fout.write("{}, ".format(mn))

            fout.write("}; }\n")

        fout.write(
"""            default : {{ throw std::out_of_range(\" Z > {}\"); }}
        }}
    }}
    
    size_type sym_2_Z_(const std::string& sym) const override {{
        auto ci_sym = sym;
        std::transform(ci_sym.begin(), ci_sym.end(), ci_sym.begin(), ::tolower);
        """.format(max_Z))

        for k in sorted_keys:
            ai = atoms[str(k)]

            fout.write("if(ci_sym == \"{}\") {{ return {}; }}\n".format(
                ai.sym.lower(), ai.Z))

            fout.write("        else ")

        fout.write(
"""{{ throw std::out_of_range(\"Unrecognized atomic symbol\"); }}
    }}
}};
        
std::unique_ptr<PeriodicTablePIMPL> nwx_default_ptable() {{
    return std::make_unique<NWXPeriodicTablePIMPL>();
}}
}} // namespace chemcache::detail_
""".format())

def write_tests(out_dir, amu2me, atoms):
    """Generate unit tests for periodic table PIMPL class.

    :param out_dir: Output directory for the generated header file.
    :type out_dir: str

    :param amu2me: Ratio of mass of electron to a Dalton.
    :type amu2me: float

    :param atoms: Collection of atoms.
    :type atoms: dict of :class:`AtomicData`
    """

    sorted_keys = sorted([int(x) for x in atoms.keys()])
    max_Z = sorted_keys[-1]

    with open(os.path.join(out_dir, "periodic_table.cpp"), 'w') as fout:
        fout.write(
"""/*
* This file has been autogenerated by generate_atomicinfo.py.  Any changes made
* to it will be lost next time generate_atomicinfo.py is run.
*/

#include <libchemist/PeriodicTable.hpp>
#include <catch/catch.hpp>

using namespace libchemist;
using size_type = typename PeriodicTable::size_type;
using isotope_list = typename PeriodicTable::isotope_list;

void test_ptable(const PeriodicTable& ptable) {{
    REQUIRE(ptable.max_Z() == {});
""".format(max_Z))

        for Z in range(1, max_Z + 1):
            ai = atoms[str(Z)]
            fout.write("    REQUIRE(ptable.get_atom({}) == Atom{{{}, \"{}\", "
                    "{}ul}});\n".format(Z, ai.mass*amu2me, ai.sym, Z))
            fout.write("    REQUIRE(ptable.sym_2_Z(\"{}\") == {});\n".format(
                ai.sym, Z))
            fout.write("    REQUIRE(ptable.isotopes({}) == isotope_list{{"
                    "".format(Z))
            sorted_mn = sorted([int(mn) for mn in ai.isotopes])
            
            for mn in sorted_mn:
                fout.write("{}, ".format(mn))
            
            fout.write("});\n")
            
            for mn in sorted_mn:
                mi = ai.isotope_masses[str(mn)]*amu2me
                fout.write("    REQUIRE(ptable.get_isotope({}, {}) == Atom{{{}ul, "
                        "\"{}\", {}}});\n".format(Z, mn, Z, ai.sym, mi))
        
        fout.write(
"""}

TEST_CASE("PeriodicTable Class") {
    SECTION("Typedefs") {
        REQUIRE(std::is_same_v<size_type, std::size_t>);
        REQUIRE(std::is_same_v<isotope_list, std::vector<size_type>>);
    }
    
    PeriodicTable ptable;
    SECTION("Default CTor") { 
        test_ptable(ptable); 
    }
    SECTION("Copy CTor") {
        PeriodicTable ptable2(ptable); 
        test_ptable(ptable2);  
    }
    SECTION("Copy assign") {
        PeriodicTable ptable2;
        auto& pptable = (ptable2 = ptable);
        test_ptable(ptable2);
        REQUIRE(&pptable == &ptable2);
    }
    SECTION("Move CTor") {
        PeriodicTable ptable2(std::move(ptable)); 
        test_ptable(ptable2); 
    }
    SECTION("Move assign") {
        PeriodicTable ptable2;
        auto& pptable = (ptable2 = std::move(ptable));
        test_ptable(ptable2);
        REQUIRE(&pptable == &ptable2);
    }
}
""")

def main(args):
    """Entry point function to generate atomic info files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    # Get and set some paths
    my_dir    = os.path.dirname(os.path.realpath(__file__))
    data_dir  = os.path.join(my_dir, "physical_data") #Dir w/ files
    out_dir   = os.path.abspath(args.src_dir)
    test_dir  = os.path.abspath(args.test_path)
    name_file = os.path.join(data_dir, "ElementNames.txt")
    iso_file  = os.path.join(data_dir, "CIAAW-ISOTOPEMASSES.txt")
    mass_file = os.path.join(data_dir, "CIAAW-MASSES.txt")
    #cov_file = os.path.join(data_dir, "CovRadii.txt")
    #vdw_file = os.path.join(data_dir, "VanDerWaalRadius.txt")
    #mult_file = os.path.join(data_dir, "NIST-ATOMICION.formatted.txt")

    # Parse atomic data
    atoms = dict()
    parse_symbols(name_file, atoms)
    parse_ciaaw_isotopes(iso_file, atoms)
    parse_ciaww_mass(mass_file, atoms)

    write_ptable(out_dir, args.amu2me, atoms)
    write_tests(test_dir, args.amu2me, atoms)

def parse_args():
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser()
    
    parser.add_argument('src_dir', type=str,
                        help="Destination directory for generated source files.")
    parser.add_argument('test_path', type=str,
                        help="Destination directory for generated unit tests.")
    parser.add_argument('--amu2me', type=float,
                        default=1822.888486192,
                        help="""Ratio of mass of electron to a Dalton. 
                             (Default: 1822.888486192)""")

    return parser.parse_args()

if __name__ == '__main__' :
    args = parse_args()

    main(args)