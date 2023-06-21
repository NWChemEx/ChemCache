#
# Copyright 2023 NWChemEx-Project
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

import pluginplay
import chemist
import simde
import chemcache
import unittest


class TestAtomicDensities(unittest.TestCase):

    def test_sto3g_atomic_dm(self):

        mm = pluginplay.ModuleManager()
        chemcache.load_modules(mm)
        mod = mm.at('sto-3g atomic dm')

        # TODO: Test when property type is exposed
