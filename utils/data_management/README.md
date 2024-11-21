<!--
  ~ Copyright 2022 NWChemEx-Project
  ~
  ~ Licensed under the Apache License, Version 2.0 (the "License");
  ~ you may not use this file except in compliance with the License.
  ~ You may obtain a copy of the License at
  ~
  ~ http://www.apache.org/licenses/LICENSE-2.0
  ~
  ~ Unless required by applicable law or agreed to in writing, software
  ~ distributed under the License is distributed on an "AS IS" BASIS,
  ~ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ~ See the License for the specific language governing permissions and
  ~ limitations under the License.
-->

Data Management
===============

The scripts in this directory are used to populate the `reference_data` 
directory and then generate source code from that data.
    
To Test generate_basis.py Locally
=================================

1. Get the data. Try running (in the root of the repo):
   ```
   python3 -m venv venv
   .github/scripts/download_reference_data.sh
   ```
2. Run generate_basis.py
   ```
   PYTHONPATH=/path/to/chemcache/utils \
   python3 generate_basis.py \
   -a /path/to/reference_data/physical_data \
   _bases src/chemcache/bases/
   ```