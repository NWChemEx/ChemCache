#!/usr/bin/env python3

"""This script uses a web scraper to download all basis sets from the Basis 
Set Exchange (BSE) in the specified output format. 
"""

import argparse
import os # for os.path manipulation
import requests
import sys # for sys.stdout.flush

import helper_fxns as helpers

class BSEBasisSetScraper:
    """Web scraper to download basis sets from the Basis Set Exchange (BSE).
    """

    def __init__(self, base_url="https://www.basissetexchange.org", 
                 user_agent="NWChemEx BSE Basis Set Scraper",
                 email="", format="nwchem",
                 uncontract_general=False, uncontract_segmented=False,
                 uncontract_spdf=False, optimize_general=False,
                 make_general=False, header_toggle=True):
        """Initializer for BSEBasisSetScraper. Many of the parameter 
        descriptions come from the 
        `BSE REST API documentation <https://molssi-bse.github.io/
        basis_set_exchange-apidoc/reference.html#api-formats>`__. More details 
        about each option may be found there.

        :param base_url: Base URL for REST API, defaults to 
            "https://www.basissetexchange.org"
        :type base_url: str, optional

        :param user_agent: User agent description for BSE, defaults to 
            "NWChemEx BSE Basis Set Scraper"
        :type user_agent: str, optional

        :param email: Email to send to BSE (not shared), defaults to ""
        :type email: str, optional

        :param format: BSE basis set file format identifier, 
            defaults to "nwchem"
        :type format: str, optional

        :param uncontract_general: Remove general contractions, 
            defaults to False
        :type uncontract_general: bool, optional

        :param uncontract_segmented: Remove general and segmented contractions,
            defaults to False
        :type uncontract_segmented: bool, optional

        :param uncontract_spdf: Split combined sp, spd, etc, shells, 
            defaults to False
        :type uncontract_spdf: bool, optional

        :param optimize_general: Optimize general contractions, 
            defaults to False
        :type optimize_general: bool, optional

        :param make_general: Make the basis set as generally-contracted as 
            possible, defaults to False
        :type make_general: bool, optional

        :param header_toggle: Toggle BSE information header on or off, 
            defaults to True
        :type header_toggle: bool, optional
        """

        self.base_url = base_url
        self.set_header(user_agent, email)
        self.valid_basis_sets = self.download_valid_basis_sets()
        self.valid_formats = self.download_valid_formats()

        self.set_default_format(format)
        self.default_header_toggle = header_toggle
        self.default_make_general = make_general
        self.default_optimize_general = optimize_general
        self.default_uncontract_general = uncontract_general
        self.default_uncontract_segmented = uncontract_segmented
        self.default_uncontract_spdf = uncontract_spdf

    def create_params(self, elements=""):
        
        params = {
            "elements"            : elements,
            "header"              : self.default_header_toggle,
            "make_general"        : self.default_make_general,
            "optimize_general"    : self.default_optimize_general,
            "uncontract_general"  : self.default_uncontract_general,
            "uncontract_segmented": self.default_uncontract_segmented,
            "uncontract_spdf"     : self.default_uncontract_spdf
        }

        return params

    def download_basis_set(self, basis_name, elements=""):
        """Download a single basis set. An optional string of elements can be
        provided or left empty to get all elements.

        :param basis_name: BSE basis set name identifier.
        :type basis_name: str
        
        :param elements: Comma-separated string of atomic numbers, 
            defaults to ""
        :type elements: str, optional

        :raises RuntimeError: Basis set could not be obtained from BSE.
        
        :return: Basis set name cleaned to be a file name and the text for the 
            basis set file.
        :rtype: tuple
        """
        
        self.validate_basis_set_name(basis_name)

        # Let the user know that the download has started
        print("Downloading {}...".format(basis_name), end='')
        # Print immediately
        sys.stdout.flush()

        bs = basis_name.replace(" ", "%20").lower()

        url  = self.base_url + "/api"
        url += "/basis/{}".format(bs)
        url += "/format/{}".format(self.default_format.lower())

        # Set additional parameters
        params = self.create_params(elements)

        # Request the basis set from BSE
        response = requests.get(url, params=params, headers=self.headers)

        if (response.status_code != 200):
            print(response.text)
            raise RuntimeError("Could not obtain {}.".format(self.basis_set))

        clean_basis_set_name = bs.lower().replace('%20', '_').replace('*', 
            "_star").replace('/','_')
        
        # Notify the user that the download is complete
        print("complete.")
        
        return clean_basis_set_name, response.text

    def download_valid_basis_sets(self):
        """Download the list of basis sets available from BSE.

        :return: Collection of basis set names
        :rtype: list
        """

        # Let the user know that the download has started
        print("Downloading BSE basis set list...", end='')
        # Print immediately
        sys.stdout.flush()

        url = self.base_url + "/api/metadata"

        response = requests.get(url, headers=self.headers)

        metadata = response.json()

        # Notify the user that the download is complete
        print("complete.")

        return metadata.keys()

    def download_valid_formats(self):
        """Download the list of formats available from BSE.

        :return: Collection of format names
        :rtype: list
        """

        # Let the user know that the download has started
        print("Downloading BSE format list...", end='')
        # Print immediately
        sys.stdout.flush()

        url = self.base_url + "/api/formats"

        response = requests.get(url, headers=self.headers)

        self.format_dict = response.json()

        # Notify the user that the download is complete
        print("complete.")

        return self.format_dict.keys()

    def get_extension(self, format=""):
        """Get the extension for the given BSE format identifier. If no 
        format identifier is given, the class default is used.

        :param format: BSE format identifier, defaults to ""
        :type format: str, optional

        :return: Basis set file extension
        :rtype: str
        """

        if (format != ""):
            self.validate_format_name(format)
        else:
            format = self.default_format

        return helpers.lookup_extension(format)

    def set_header(self, user_agent="", email=""):
        """Generates the header to use in requests.

        :param user_agent: Description of who is pinging their API, 
            defaults to ""
        :type user_agent: str, optional

        :param email: Email to send (not shared), defaults to ""
        :type email: str, optional
        """
        
        self.headers = {
            "User-Agent" : user_agent,
            "From" : email
        }

    def set_default_format(self, format):
        """Set the default format for basis sets.

        :param format: Valid BSE format identifier for basis sets.
        :type format: str
        """
        
        self.validate_format_name(format)

        self.default_format = format

    def validate_basis_set_name(self, basis_name):
        """Validate the basis name against the list of valid basis names 
        retrieved from BSE.

        :param basis_name: Name of the basis set
        :type basis_name: str

        :raises RuntimeError: Invalid basis name was given.
        """

        if (not basis_name.lower() in self.valid_basis_sets):
            raise RuntimeError("Invalid basis name: {}".format(basis_name))

    def validate_format_name(self, format):
        """Validate the format name against the list of valid format names 
        retrieved from BSE.

        :param format: Name of the formatting option
        :type format: str

        :raises RuntimeError: Invalid format option was given
        """

        if (not format.lower() in self.valid_formats):
            raise RuntimeError("Invalid format option: {}".format(format))

def write_basis_set(destination, basis_name, basis_data, extension):
    """Write the basis set out to a file.

    :param basis_name: Name of the basis set.
    :type basis_name: str

    :param basis_data: Text data to write to the basis set file
    :type basis_data: str

    :param extension: Extension for the basis set file. Must include the dot 
        '.' separator if one is needed.
    :type extension: str
    """

    # Let the user know that the writing has started
    print("Writing {} to disk...".format(basis_name), end='')
    # Print immediately
    sys.stdout.flush()

    basis_path = os.path.join(destination, basis_name + extension)
    with open(basis_path,'w') as fout:
        fout.write(basis_data)

    # Notify the user that the writing is complete
    print("complete.")

def main(args):
    """Entry point function to generate basis set files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    scraper = BSEBasisSetScraper(format=args.outformat, 
        optimize_general=args.optimize_general
    )

    print("---")

    for name in scraper.valid_basis_sets:
        clean_name, text = scraper.download_basis_set(name)

        write_basis_set(args.destination, clean_name, text, 
                        scraper.get_extension())

        print("---")

def parse_args():
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('destination', type=str,
                        help="Destination directory for basis set files.")
    parser.add_argument('-o','--outformat', type=str,
                        default="NWChem", 
                        help="Output format. (Default: NWChem)")
    parser.add_argument('-g', '--optimize_general', action="store_true",
                        help="""Toggle on optimizing general contractions. 
                             Default OFF.""")

    return parser.parse_args()

if __name__ == '__main__' :
    args = parse_args()

    main(args)