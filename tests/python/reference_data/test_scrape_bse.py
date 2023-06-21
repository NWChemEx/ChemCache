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

import unittest
import os
import sys

sys.path.append(
    os.path.join(os.path.dirname(os.path.realpath(__file__)),
                 "../../../reference_data"))

from reference_data import scrape_bse


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
        corr = [
            '3-21g', '4-31g', '5-21g', '6-21g', '6-31++g', '6-31++g_st_',
            '6-31++g_st__st_', '6-31+g', '6-31+g_st_', '6-31+g_st__st_',
            '6-311++g', '6-311++g(2d,2p)', '6-311++g(3df,3pd)', '6-311++g_st_',
            '6-311++g_st__st_', '6-311+g', '6-311+g(2d,p)', '6-311+g_st_',
            '6-311+g_st__st_', '6-311g', '6-311g(2df,2pd)', '6-311g(d,p)',
            '6-311g_st_', '6-311g_st__st_', '6-31g', '6-31g(2df,p)',
            '6-31g(3df,3pd)', '6-31g(d,p)', '6-31g_st_', '6-31g_st__st_',
            'aug-cc-pcv5z', 'aug-cc-pcvdz', 'aug-cc-pcvqz', 'aug-cc-pcvtz',
            'aug-cc-pv(5+d)z', 'aug-cc-pv(d+d)z', 'aug-cc-pv(q+d)z',
            'aug-cc-pv(t+d)z', 'aug-cc-pv5z', 'aug-cc-pv6z', 'aug-cc-pv7z',
            'aug-cc-pvdz', 'aug-cc-pvqz', 'aug-cc-pvtz', 'aug-cc-pwcv5z',
            'aug-cc-pwcvdz', 'aug-cc-pwcvqz', 'aug-cc-pwcvtz', 'aug-pv7z',
            'cc-pcv5z', 'cc-pcvdz', 'cc-pcvqz', 'cc-pcvtz', 'cc-pv(5+d)z',
            'cc-pv(d+d)z', 'cc-pv(q+d)z', 'cc-pv(t+d)z', 'cc-pv5z', 'cc-pv6z',
            'cc-pv8z', 'cc-pv9z', 'cc-pvdz', 'cc-pvdz(seg-opt)', 'cc-pvqz',
            'cc-pvqz(seg-opt)', 'cc-pvtz', 'cc-pvtz(seg-opt)', 'cc-pwcv5z',
            'cc-pwcvdz', 'cc-pwcvqz', 'cc-pwcvtz', 'd-aug-cc-pv5z',
            'd-aug-cc-pv6z', 'd-aug-cc-pvdz', 'd-aug-cc-pvqz', 'd-aug-cc-pvtz',
            'pv6z', 'pv7z', 'aug-cc-pv5z-pp-optri', 'aug-cc-pvdz-pp-optri',
            'aug-cc-pvqz-pp-optri', 'aug-cc-pvtz-pp-optri',
            'aug-cc-pwcv5z-pp-optri', 'aug-cc-pwcvdz-pp-optri',
            'aug-cc-pwcvqz-pp-optri', 'aug-cc-pwcvtz-pp-optri'
        ]

        self.scraper.add_filter("family",
                                ["pople", "dunning", "dunning_pp_fit"])
        self.scraper.add_filter("role", ["orbital", "optri"])

        self.assertEqual(len(self.scraper.filtered_basis_sets), len(corr))
        self.assertEqual(self.scraper.filtered_basis_sets, corr)

    def test_download_basis_set(self):
        basis_name = "3-21g"
        clean_name, text = self.scraper.download_basis_set(basis_name)
        corr = ""

        # Read golden file
        current_dir = os.path.dirname(os.path.realpath(__file__))
        with open(
                os.path.join(current_dir,
                             "../../data/{}.nw".format(basis_name)),
                'r') as fin:
            corr = fin.read()

        self.assertEqual(clean_name, basis_name)
        self.assertEqual(text, corr)


if __name__ == "__main__":
    unittest.main()
