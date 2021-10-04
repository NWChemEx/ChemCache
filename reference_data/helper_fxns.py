"""This file contains helper functions used throughout the various parsing 
scripts in this directory.
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
    temp = temp.replace(",", "_comma_")
    temp = temp.replace("(", "_oparen_")
    temp = temp.replace(")", "_cparen_")
    temp = temp.replace("!", "_exclamation_")
    temp = temp.replace("?", "_question_")

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
def write_warning(fout, script_name, prefix=""):
    """Writes a warning saying the file was auto-generated by the given script
    and will be overwritten the next time it is run.

    :param fout: Text file opened for output.
    :type fout: :class:io.TextIOBase

    :param script_name: Name of the script that generates the file.
    :type script_name: str

    :param prefix: An optional, additional prefix to be added to every line of
        the warning comment. This is a quick an dirty fix to output this 
        warning to non-C++ files, where comments may be different.
        Defaults to "".
    :type prefix: str
    """

    warning = (
"""{1}/*
{1} * This file is autogenerated by: {0} 
{1} * 
{1} * NOTE: Any modifications made in this file will be lost next time 
{1} *       {0} is run.
{1} */

""").format(script_name, prefix)
    
    fout.write(warning)
