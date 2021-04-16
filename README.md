# ReferenceData

This repo contains basically any data that you can just look up, such as 
basis sets, physical constants, geometries for some well known molecules, etc.

There are multiple parts to this repo. The product of this repo is a library
that provides a variety of fixed data structures that can be delivered to other
components that need that data. This library is compiled from source code but as
such source code is tedious to maintain the source is generated from data files.
A set of Python scripts implement the data to source transformers. As a result 
the fully built library contains the following pieces:

- Data files
- Data to source transformer scripts
- Source code
- Binary library files

## Data files

Data files contain the data that the library ultimately will serve to its 
users. The data files are grouped into different categories dependent on the
kind of information they contain:

### ReferenceData/atomic_densities

Historically it has been shown that a good starting point for electronic
structure calculations can be obtained from a
[superposition of atomic densities](https://doi.org/10.1002/jcc.20393).
This collection stores precomputed atomic densities for a number of
basis sets, that can be used to construct the initial guess.

### ReferenceData/basis_sets

There is the potential to store multiple collections of basis sets. At present
the only is the "default" collection. This collection corresponds to the entire
contents of the [Basis Set Exchange](https://www.basissetexchange.org/). 
Each basis set is stored in the [Gaussian'94 basis set format](https://gaussian.com/gen/).

### ReferenceData/molecules

The structures of well known molecules or molecules that are related to the
science cases of the project. The structures are stored as files in the 
[XYZ format](https://en.wikipedia.org/wiki/XYZ_file_format).

### ReferenceData/physical_data

Different constants and physical constants from multiple sources and in case specific formats:
- Isotope abundances (source [CIAAW](http://www.ciaaw.org/isotopic-abundances.htm) 2013)
- Isotope masses (source [CIAAW](http://www.ciaaw.org/atomic-masses.htm) 2012)
- Atomic masses (source [CIAAW](http://www.ciaaw.org/atomic-weights.htm) 2015)
- Element names (source CIAAW, except for [ElementZero (Ez)](https://masseffect.fandom.com/wiki/Element_Zero)
  which was proposed in the "Mass Effect" video game series). Note, for non-atom centers the label
  "Bq" from "Banquo" has been used in the past. See: [NICS](https://doi.org/10.1021/cr030088+), 
  [QLST](https://doi.org/10.1016/S0065-3276\(08\)00402-4), [Tapia](10.1103/PhysRevA.84.012115),
  [Gastreich](http://shakespeare.nowheres.com/qandr/others/3.15.97/messages/199.html),
  [CCL](http://www.ccl.net/chemistry/resources/messages/2009/03/06.011-dir/index.html),
  [NWChem](https://github.com/nwchemgit/nwchem/tree/master/QA/tests/dft_bsse),
  [Gaussian](https://gaussian.com/molspec/). Gaussian allows for an Oxygen ghost atom to be specified
  by "O-Bq".
  The [Psi4](https://psicode.org/psi4manual/master/psithonmol.html) code handles ghost atoms
  differently, requiring either "Gh(O)" or "@O".
  [MOLPRO](https://www.molpro.net/info/release/doc/manual.pdf) uses "Q" for ghost centers.
  [MPQC](https://valeevgroup.github.io/mpqc-docs/v2/mpqcsimp.html) allows setting the charge of a center to 0.
  [GAMESS](https://www.msg.chem.iastate.edu/gamess/GAMESS_Manual/docs-input.txt) uses "Bq".)
- Covalent radii 
- van der Waals radii
- Atomic electron configurations [NIST](http://www.nist.gov/pml/data/ion_energy.cfm)
- Physical constants [NIST](http://physics.nist.gov/constants)

## Transformer scripts

## Source code

## Binary library files
