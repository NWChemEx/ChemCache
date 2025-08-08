#
# Copyright 2024 NWChemEx-Project
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
#

import unittest

from chemist import Atom, Molecule
from pluginplay import ModuleManager
from simde import MoleculeFromString

from chemcache import load_modules


class TestMolecules(unittest.TestCase):
    def test_nwx_molecules(self):
        rv_mol = self.mm.run_as(MoleculeFromString(), "NWX Molecules", "water")

        o_mass = 29165.122045980286
        o = Atom("O", 8, o_mass, 0.0, -0.1432223429807816, 0.0)

        h_mass = 1837.4260218693814
        h1 = Atom("H", 1, h_mass, 1.6380335020342418, 1.1365568803584036, 0.0)
        h2 = Atom("H", 1, h_mass, -1.6380335020342418, 1.1365568803584036, 0.0)

        water = Molecule()
        water.push_back(o)
        water.push_back(h1)
        water.push_back(h2)

        self.assertEqual(rv_mol, water)

        with self.assertRaises(IndexError):
            self.mm.run_as(MoleculeFromString(), "NWX Molecules", "Nothing")

    def setUp(self):
        self.mm = ModuleManager()
        load_modules(self.mm)
