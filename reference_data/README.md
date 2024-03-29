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

Reference Data Directory
========================

This directory contains externally known data that is needed for NWChemEx to
run.  The files contained in the subdirectories are ASCII lists/tables of the
data.  Python scripts exist to turn this data into source files and should be
used to update files should the need arise.

The reference data within this directory is further subdivided into the 
following types:

- Physical Data
  - Atomic masses, names, abundances, isotopes, physical constants, *etc.*
- Basis Sets
  - Standard Gaussian basis sets like 6-31g* and aug-cc-pVDZ
- Atomic guesses
  - Precomputed atomic SCF results for initial guesses, from select basis sets
- Molecules
  - Pre-defined molecules for quickly and easily building tests or running 
  pre-defined sets
      
    
