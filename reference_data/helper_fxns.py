"""
This file contains helper functions used throughout the various parsing scripts
in this directory
"""

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
    
# Writes warning about auto-generation to file
def write_warning(script_name, f):
    f.write("/* Warning!!!! \n")
    f.write(" * This file was auto-generated by {} ".format(script_name))
    f.write("and should not be\n")
    f.write(" * modified.  To make modifications to this file please edit\n")
    f.write(" * {} and regenerate this file.\n".format(script_name))
    f.write(" */\n")

def write_header(script_name, f, nmspace="libchemist"):
    f.write("#include \"libchemist/nwx_defaults.hpp\"\n")
    write_warning(script_name, f)
    f.write("namespace {} {{\n".format(nmspace))

def write_footer(f, nmspace="libchemist"):
    f.write("    return rv;\n")
    f.write("}}\n}} // End namespace {}\n".format(nmspace))