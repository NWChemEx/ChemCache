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

from chemcache import load_modules

from chemist import Atom
from chemist import Molecule
from chemist import ShellType
from chemist import PointD
from chemist.basis_set import ShellD
from chemist.basis_set import AtomicBasisSetD
from chemist.basis_set import AOBasisSetD

from pluginplay import ModuleManager

from simde import AtomicBasisSetFromZ
from simde import MolecularBasisSet


class TestBases(unittest.TestCase):

    def test_atomic_basis_set(self):
        rv_bs = self.mm.run_as(AtomicBasisSetFromZ(), "sto-3g atomic basis", 1)
        self.assertEqual(rv_bs, self._h(0.0))

        rv_bs = self.mm.run_as(AtomicBasisSetFromZ(), "sto-3g atomic basis", 8)

        with self.assertRaises(IndexError):
            self.mm.run_as(AtomicBasisSetFromZ(), "sto-3g atomic basis", 1000)

    def test_molecular_basis_set(self):
        h1 = Atom("H", 1, 0.0, 1.0, 1.0, 1.0)
        h2 = Atom("H", 1, 0.0, 2.0, 2.0, 2.0)
        inp = Molecule()
        inp.push_back(h1)
        inp.push_back(h2)

        rv_bs = self.mm.run_as(MolecularBasisSet(), "sto-3g", inp)

        corr_bs = AOBasisSetD()
        corr_bs.add_center(self._h(1.0))
        corr_bs.add_center(self._h(2.0))

        self.assertEqual(rv_bs, corr_bs)

    def _h(self, r):
        cs = [1.5432896730e-01, 5.3532814230e-01, 4.4463454220e-01]
        es = [3.4252509140e+00, 6.2391372980e-01, 1.6885540400e-01]
        shell = ShellD(ShellType.pure, 0, cs, es, r, r, r)
        h = AtomicBasisSetD("sto-3g", 1, PointD(r, r, r), [shell])
        return h

    def setUp(self):
        self.mm = ModuleManager()
        load_modules(self.mm)
