## Dependencies
This is an example to build on local PC with scotch-pastix solver and mmg-gmsh for mesh refinement.
To run the code one needs a fortran compiler, openmpi (or other MPI), BLAS, LAPACK, LAPACKE, Xt, HDF5, SCOTCH+PASTIX (or PETSc), MMG+GMSH libraries
### Debian
```zsh
sudo apt-get install gfortran openmpi-bin openmpi-common openmpi-doc libopenmpi-dev libblas-dev liblapacke-dev  liblapack-dev libxt-dev libhdf5-serial-dev 
```
Then install [SCOTCH](https://gitlab.inria.fr/scotch/scotch) [Pastix](https://solverstack.gitlabpages.inria.fr/pastix/md_docs_doxygen_chapters_Pastix_MPI.html). Install int32 versions of both.
Next, install [mmg](https://github.com/MmgTools/Mmg/wiki/Setup-guide) and load gmsh 4.11.1 SDK [pack](https://gmsh.info/bin/Linux/) (only tested with this gmsh version).

## Building
Make sure that you specify properly the paths to needed libraries.
By defualt `arch.make` file contains instructions for local build with gfortran, OpenMPI, SCOTCH+PASTIX, MMG+GMSH libraries.
```zsh
cd lib
source Make.inc/init_vars_libs.sh
make
```

## Mesh generation

The machine wall outline is needed with knowing the locations of puff and pump.
Then a `.geo` file is needed to make `.msh` and `.mesh` files using GMSH GUI (for example).
`.geo` should be put in `test/res/geometries` folder, `.msh` and `.mesh` files to `test/Meshes`
The refinement level for initial mesh should be of order of 1k elements for a machine like TCV, then it will be automatically refined.

