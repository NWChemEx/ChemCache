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

import os
import unittest

from data_management import scrape_bse


class ScrapeBSETest(unittest.TestCase):
    # This method prevents multiple downloads, greatly reducing the testing
    # time. However, the same scraper object will be used in each test, so
    # care must be taken to reset the state of the scraper in a test if
    # previous state changes will affect the result
    # scraper = scrape_bse.BSEBasisSetScraper()

    # This method will download the BSE basis set list and format list
    # each time a test is run. This may increase the testing time if many
    # tests are set up, which could get expensive.
    def setUp(self):
        self.scraper = scrape_bse.BSEBasisSetScraper()

    def test_add_filter(self):
        corr = ["sto-2g", "sto-3g", "sto-3g_st_", "sto-4g", "sto-5g", "sto-6g"]

        self.scraper.add_filter("family", ["sto"])
        self.scraper.add_filter("role", ["orbital"])
        self.assertEqual(len(self.scraper.filtered_basis_sets), len(corr))
        self.assertEqual(self.scraper.filtered_basis_sets, corr)

    def test_download_basis_set(self):
        basis_name = "sto-3g"
        clean_name, text = self.scraper.download_basis_set(basis_name)
        corr = ""

        # Read golden file
        current_dir = os.path.dirname(os.path.realpath(__file__))
        ref_path = os.path.join(
            "..", "..", "..", "reference_data", "basis_sets", "default"
        )
        ref_dir = os.path.abspath(os.path.join(current_dir, ref_path))
        basis_file = f"{basis_name}.nw"

        with open(os.path.join(ref_dir, basis_file), "r") as fin:
            corr = fin.read()

        self.assertEqual(clean_name, basis_name)
        self.assertEqual(text, corr)
