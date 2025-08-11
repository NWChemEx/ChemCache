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

from chemist import Atom
from pluginplay import ModuleManager
from simde import AtomFromZ, SymbolFromZ, ZFromSymbol

from chemcache import load_modules


class TestAtoms(unittest.TestCase):
    def test_atoms_average(self):
        h = Atom("H", 1, 1837.4260218693814, 0.0, 0.0, 0.0)
        o = Atom("O", 8, 29165.122045980286, 0.0, 0.0, 0.0)

        rv_atom = self.mm.run_as(AtomFromZ(), "Atom", 1)
        self.assertEqual(rv_atom, h)

        rv_atom = self.mm.run_as(AtomFromZ(), "Atom", 8)
        self.assertEqual(rv_atom, o)

    def test_sym_from_Z(self):
        rv_sym = self.mm.run_as(SymbolFromZ(), "Symbol from Z", 1)
        self.assertEqual(rv_sym, "H")

        rv_sym = self.mm.run_as(SymbolFromZ(), "Symbol from Z", 8)
        self.assertEqual(rv_sym, "O")

        with self.assertRaises(IndexError):
            rv_sym = self.mm.run_as(SymbolFromZ(), "Symbol from Z", 1000)

    def test_Z_from_sym(self):
        rv_Z = self.mm.run_as(ZFromSymbol(), "Z From Symbol", "H")
        self.assertEqual(rv_Z, 1)

        rv_Z = self.mm.run_as(ZFromSymbol(), "Z From Symbol", "O")
        self.assertEqual(rv_Z, 8)

        with self.assertRaises(IndexError):
            rv_Z = self.mm.run_as(ZFromSymbol(), "Z From Symbol", "NWX")

    def setUp(self):
        self.mm = ModuleManager()
        load_modules(self.mm)
