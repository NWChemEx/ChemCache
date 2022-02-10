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
