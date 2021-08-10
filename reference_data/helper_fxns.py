"""
This file contains helper functions used throughout the various parsing scripts
in this directory
"""

import os
from os import TMP_MAX

def find_files(source_root, extensions, recursive=False):
    """Recursively find all files in the given directory with the provided
    extensions.

    :param source_root: Root directory containing relevant files.
    :type source_root: str

    :param extensions: File extensions to search for.
    :type extensions: list of str

    :param recursive: Whether or not to recursively search in subdirectories,
        defaults to False
    :type recursive: bool
    
    :return: Dictionary of file extensions mapping to lists of full paths to 
        files found with the respective extension
    :rtype: dict of lists
    """

    found_files = {}

    for extension in extensions:
        found_files[extension] = []

        # Recursively search for atomic density files of the given extension
        for dirpath, _, filenames in os.walk(source_root):
            for filename in filenames:
                _, ext = os.path.splitext(filename)

                # Case insensitive comparison of extensions
                if (ext.lower() == extension.lower()):
                    found_files[extension].append(
                        os.path.join(dirpath, filename)
                    )

            # Stop the recursive search before diving into subdirectories
            # if recursion is not requested
            if (not recursive):
                break

    return found_files

def lookup_extension(format):
    """Looks up the extension for a basis set file format.

    :param format: BSE basis set file format identifier.
    :type format: str
    :return: Extension for the file format.
    :rtype: str 
    """

    known_formats = {
        "acesii"    : ".acesii", 
        "bdf"       : ".bdf", 
        "bsedebug"  : ".bse", 
        "cfour"     : ".c4bas", 
        "cp2k"      : ".cp2k", 
        "dalton"    : ".dalton", 
        "demon2k"   : ".d2k", 
        "gamess_uk" : ".bas", 
        "gamess_us" : ".bas", 
        "gaussian94": ".gbs", 
        "json"      : ".json", 
        "molcas"    : ".molcas", 
        "molpro"    : ".mpro", 
        "nwchem"    : ".nw", 
        "orca"      : ".orca", 
        "pqs"       : ".pqs", 
        "psi4"      : ".gbs", 
        "qchem"     : ".qchem", 
        "qcschema"  : ".json", 
        "turbomole" : ".tm", 
        "xtron"     : ".gbs"
    }

    return known_formats[format.lower()]

def sanitize_basis_name(bs_name):
    """Sanitizes the basis set name. For example, replace numbers with the
    corresponding word and hyphens, plus-signs, and other symbols with 
    underscores.

    :param bs_name: Basis set name
    :type bs_name: str

    :return: Sanitized basis set name
    :rtype: str
    """

    temp = bs_name.replace("6-", "six_dash_")
    temp = temp.replace("3-", "three_dash_")
    temp = temp.replace("-", "_dash_")
    temp = temp.replace("-", "_dash_")
    temp = temp.replace("+", "_plus_")

    return temp

def desanitize_basis_name(bs_name):
    """\"Desanitizes\" the basis set name. For example, replace "_star" with
    an asterisk character '*'.

    :param bs_name: Sanitized basis set name
    :type bs_name: str

    :return: Unsanitized basis set name
    :rtype: str
    """

    temp = bs_name.replace("_star", "*")
    
    return temp

# Writes warning about auto-generation to file
def write_warning(ostream, script_name):
    warning = (
"""/*
 * This file is autogenerated by {0}. Any modifications made in
 * this file will be lost next time {0} is run.
 */

""").format(script_name)
    
    ostream.write(warning)

def write_header(io, script_name, nmspace="libchemist"):
    io.write("#include \"libchemist/nwx_defaults.hpp\"\n")
    write_warning(script_name, io)
    io.write("namespace {} {{\n".format(nmspace))

def write_footer(io, nmspace="libchemist"):
    io.write("    return rv;\n")
    io.write("}}\n}} // End namespace {}\n".format(nmspace))