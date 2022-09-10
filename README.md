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

# ChemCache

This repo contains chemical reference data, *e.g.*, atomic basis sets, physical constants, 
molecular geometries, *etc.*.

## Repo structure

ChemCache has two release branches: `master` and `generated_data`. `master` only version
controls the scripts for obtaining the reference data and the scripts for generating the
source files. `generated_data` additionally contains the reference data and generated
source files. Users of the repo can use whichever branch of the repo is more convenient
for them, as the branches are automatically synchronized (*i.e*, changes to `master` will
automatically cause `generated_data` to be updated).
