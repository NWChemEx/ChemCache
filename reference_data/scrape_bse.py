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

"""This script uses a web scraper to download all basis sets from the Basis 
Set Exchange (BSE) in the specified output format.

Usage
-----

::

   usage: scrape_bse.py [-h] [-o OUTFORMAT] [-g] destination

   positional arguments:
     destination           Destination directory for basis set files.

   optional arguments:
     -h, --help            show this help message and exit
     -o OUTFORMAT, --outformat OUTFORMAT
                           Output format. (Default: NWChem)
     -g, --optimize_general
                           Toggle on optimizing general contractions. Default OFF.

This script creates the following files in the given destination::

   +---destination
   |       <all_basis_set_files>
"""

import argparse
import copy  # for deepcopy
import os  # for os.path manipulation
import requests
import sys  # for sys.stdout.flush

import helper_fxns as helpers


class BSEBasisSetScraper:
    """Web scraper to download basis sets from the Basis Set Exchange (BSE).
    """

    def __init__(self, base_url: str = "https://www.basissetexchange.org",
                 user_agent: str = "NWChemEx BSE Basis Set Scraper",
                 email: str = "", format: str = "nwchem",
                 uncontract_general: bool = False,
                 uncontract_segmented: bool = False,
                 uncontract_spdf: bool = False, optimize_general: bool = False,
                 make_general: bool = False, header_toggle: bool = True) -> None:
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

        # Fetch valid basis sets
        self.valid_basis_sets, self.metadata = self.download_valid_basis_sets()
        # No filters at the beginning, so these variables are equal
        self.filtered_basis_sets = self.valid_basis_sets
        self.filtered_metadata = self.metadata
        self.filters = {}

        # Get valid formats and set the default
        self.valid_formats = self.download_valid_formats()
        self.set_default_format(format)

        # Set a bunch of other default parameters
        self.default_header_toggle = header_toggle
        self.default_make_general = make_general
        self.default_optimize_general = optimize_general
        self.default_uncontract_general = uncontract_general
        self.default_uncontract_segmented = uncontract_segmented
        self.default_uncontract_spdf = uncontract_spdf

    def add_filter(self, metadata_key: str, values: list) -> None:
        """Add a metadata filter to the basis set list and update the filtered
        basis set list and metadata.

        This function adds filters to the list of valid basis sets contained
        by this class based on metadata values scraped from BSE. If filters 
        already exist for the metadata key given, the new values will be 
        appended to the existing filter value list. Values must match exactly!

        When multiple filter values exist for a metadata key, basis sets are 
        guaranteed to contain at least one of the filter values, but not 
        necessarily all filter values for the metadata key. However, filter
        values of different metadata keys are applied sequentially, so 
        the filtered basis sets must contain at least one of the
        filter values for each metadata key.

        For example::

           scraper.add_filter("family", ["pople", "dunning"])
           scraper.add_filter("role", ["orbital", "optri"])

        will filter to all basis sets that are of either the "pople" or 
        "dunning" families, but only if they have a role of "orbital" or 
        "optri".

        The filtered basis set names can be retrieved using the data member
        `filtered_basis_sets` or the full filtered metadata can be 
        retrieved with `filtered_metadata`. 

        :param metadata_key: Key for the desired value in basis set metadata.
        :type metadata_key: str

        :param values: Values of the metadata to filter by.
        :type values: list
        """

        if (metadata_key in self.filters):
            self.filters[metadata_key].extend(values)
        else:
            self.filters[metadata_key] = values

        self.filtered_basis_sets, self.filtered_metadata = self.get_filtered_basis_sets()

    def download_basis_set(self, basis_name: str, elements: str = "") -> tuple:
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

        url = self.base_url + "/api"
        url += "/basis/{}".format(bs)
        url += "/format/{}".format(self.default_format.lower())

        # Set additional parameters
        params = self._create_params(elements)

        # Request the basis set from BSE
        response = requests.get(url, params=params, headers=self.headers)

        if (response.status_code != 200):
            print(response.text)
            raise RuntimeError("Could not obtain {}.".format(self.basis_set))

        clean_basis_set_name = bs.lower().replace(
            '%20', '_').replace('*', "_star").replace('/', '_')

        # Notify the user that the download is complete
        print("complete.")

        return clean_basis_set_name, response.text

    def download_valid_basis_sets(self) -> tuple:
        """Download the list of basis sets available from BSE.

        :return: Collections of basis set names and metadata
        :rtype: tuple of list and dict
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

        return list(metadata.keys()), metadata

    def download_valid_formats(self) -> list:
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

    def get_extension(self, format: str = "") -> str:
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

    def set_header(self, user_agent: str = "", email: str = "") -> None:
        """Generates the header to use in requests.

        :param user_agent: Description of who is pinging the BSE API, 
            defaults to ""
        :type user_agent: str, optional

        :param email: Email to send to BSE (not shared), defaults to ""
        :type email: str, optional
        """

        self.headers = {
            "User-Agent": user_agent,
            "From": email
        }

    def set_default_format(self, format: str) -> None:
        """Set the default format for basis sets.

        :param format: Valid BSE format identifier for basis sets.
        :type format: str
        """

        self.validate_format_name(format)

        self.default_format = format

    def get_filtered_basis_sets(self) -> tuple:
        """Filter the existing valid basis sets based on metadata filters 
        currently set in the class. This function does not change the class.

        :return: Returns a filtered list of basis set names and the filtered
            metadata dict
        :rtype: tuple of list and dict
        """

        # Exit early if no filters exist
        if (not len(self.filters)):
            return self.valid_basis_sets, self.metadata

        filtered = copy.deepcopy(self.metadata)

        # Filter by metadata keys set by add_filter
        for metadata_key in self.filters:
            tmp = {}

            # Get each value to filter by in that metadata key
            for filter_value in self.filters[metadata_key]:
                for bs_name in filtered.keys():
                    if (filtered[bs_name][metadata_key] == filter_value):
                        tmp[bs_name] = copy.deepcopy(filtered[bs_name])

            filtered = copy.deepcopy(tmp)

        return list(filtered.keys()), filtered

    def validate_basis_set_name(self, basis_name: str) -> None:
        """Validate the basis name against the list of valid basis names 
        retrieved from BSE.

        :param basis_name: Name of the basis set
        :type basis_name: str

        :raises RuntimeError: Invalid basis name was given.
        """

        if (not basis_name.lower() in self.valid_basis_sets):
            raise RuntimeError("Invalid basis name: {}".format(basis_name))

    def validate_format_name(self, format: str) -> None:
        """Validate the format name against the list of valid format names 
        retrieved from BSE.

        :param format: Name of the formatting option
        :type format: str

        :raises RuntimeError: Invalid format option was given
        """

        if (not format.lower() in self.valid_formats):
            raise RuntimeError("Invalid format option: {}".format(format))

    def _create_params(self, elements: str = "") -> dict:
        """Create the parameter dictionary for a BSE request.

        :param elements: Elements to retrieve bases for, defaults to ""
        :type elements: str, optional

        :return: Dictionary of parameter names (keys) and their values
        :rtype: dict
        """

        params = {
            "elements": elements,
            "header": self.default_header_toggle,
            "make_general": self.default_make_general,
            "optimize_general": self.default_optimize_general,
            "uncontract_general": self.default_uncontract_general,
            "uncontract_segmented": self.default_uncontract_segmented,
            "uncontract_spdf": self.default_uncontract_spdf
        }

        return params


def _write_basis_set(destination: str, basis_name: str, basis_data: str,
                    extension: str) -> None:
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
    with open(basis_path, 'w') as fout:
        fout.write(basis_data)

    # Notify the user that the writing is complete
    print("complete.")


def main(args: argparse.Namespace) -> None:
    """Entry point function to generate basis set files.

    :param args: Command line argument namespace
    :type args: Namespace
    """

    scraper = BSEBasisSetScraper(format=args.outformat,
                                 optimize_general=args.optimize_general
                                 )

    print("---")

    # Add metadata value filters at this point
    scraper.add_filter(
        "family", ["pople", "dunning", "dunning_fit", "dunning_pp_fit", "sto"])

    for name in scraper.filtered_basis_sets:
        clean_name, text = scraper.download_basis_set(name)

        _write_basis_set(args.destination, clean_name, text,
                        scraper.get_extension())

        print("---")


def parse_args() -> argparse.Namespace:
    """Parse command line arguments.

    :return: Values of command line arguments.
    :rtype: Namespace
    """
    parser = argparse.ArgumentParser(description=
        "This script uses a web scraper to download all basis sets from the "
        "Basis Set Exchange (BSE) in the specified output format."
    )

    parser.add_argument('destination', type=str,
                        help="Destination directory for basis set files.")
    parser.add_argument('-o', '--outformat', type=str,
                        default="NWChem",
                        help="Output format. (Default: NWChem)")
    parser.add_argument('-g', '--optimize_general', action="store_true",
                        help="""Toggle on optimizing general contractions. 
                             Default OFF.""")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    main(args)
