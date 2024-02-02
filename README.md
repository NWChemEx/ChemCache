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

ChemCache
=========

[Documentation](https://nwchemex.github.io/ChemCache)

This repo contains chemical reference data, *e.g.*, atomic basis sets, physical 
constants, molecular geometries, *etc.*.

# Features

Write me...

# Repo structure

ChemCache has two release branches: `master` and `generated_data`. `master` only
version controls the scripts for obtaining the reference data and the scripts 
for generating the source files. `generated_data` additionally contains the 
reference data and generated source files. Users of the repo can use whichever 
branch of the repo is more convenient for them, as the branches are 
automatically synchronized (*i.e*, changes to `master` will automatically cause 
`generated_data` to be updated).

# Installation

As with the majority of the NWChemEx stack, ChemCache uses CMake and the 
[CMaize](https://github.com/CMakePP/CMaize) library for configuration and 
building. This means that installation is usually achieved via a variation on:

```.sh
git clone https://github.com/NWChemEx/ChemCache
cd ChemCache
cmake -H. -Bbuild -D...
cmake --build build --target install
```
More detailed install instructions can be found
[here](https://nwchemex.github.io/ChemCache/installation.html).

# Contributing

- [Contributor Guidelines](https://github.com/NWChemEx/.github/blob/1a883d64519f62da7c8ba2b28aabda7c6f196b2c/.github/CONTRIBUTING.md)
- [Contributor License Agreement](https://github.com/NWChemEx/.github/blob/master/.github/CONTRIBUTING.md#contributor-license-agreement-cla)
- [Developer Documentation](https://nwchemex.github.io/ChemCache)
- [Code of Conduct](https://github.com/NWChemEx/.github/blob/master/.github/CODE_OF_CONDUCT.md)

# Acknowledgments

This research was supported by the Exascale Computing Project (17-SC-20-SC), a 
collaborative effort of the U.S. Department of Energy Office of Science and the 
National Nuclear Security Administration.
