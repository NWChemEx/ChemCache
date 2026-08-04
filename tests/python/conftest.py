#
# Copyright 2026 NWChemEx-Project
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

import pathlib
import sys

import parallelzone as pz
import pytest

# tests/python/utils_tests/test_scrape_bse.py imports "data_management" as a
# top-level package, but it only exists under utils/ at the repo root (not
# installed anywhere, not a pure-Python package under the wheel). Make it
# importable here.
_utils_dir = pathlib.Path(__file__).resolve().parents[2] / "utils"
if _utils_dir.is_dir():
    sys.path.insert(0, str(_utils_dir))


@pytest.fixture(scope="session", autouse=True)
def _session_runtime_view():
    """
    Holds a single RuntimeView for the whole pytest session.

    MPI may only be initialized/finalized once per process. The first
    RuntimeView constructed owns that responsibility; individual test
    modules construct their own RuntimeView per test (e.g. in setUp),
    which is safe only as long as this session-scoped instance is still
    alive to keep MPI initialized in between. Without this, pytest would
    run each test module independently and MPI would be finalized after
    the first module's tests finished, breaking every module after it.
    """
    rv = pz.runtime.RuntimeView()
    yield rv
