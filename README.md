# MHDG

Plasma-fluid simulations with the hybridizable discontinuous Galerkin method.
Scientific models, input formats and simulation workflows are indexed in
[docs/](docs/README.md); that coverage is preliminary and incremental.

## Dependencies

The default GCC build needs a Fortran compiler and MPI wrapper, BLAS/LAPACK,
HDF5 with Fortran support, the Gmsh SDK, and SCOTCH/PaStiX or PETSc. A Debian
starting point for the system dependencies is:

```bash
sudo apt-get install make pkg-config gfortran openmpi-bin libopenmpi-dev \
  libblas-dev liblapack-dev liblapacke-dev libxt-dev libhdf5-dev libhwloc-dev
```

Install [SCOTCH](https://gitlab.inria.fr/scotch/scotch) and
[PaStiX](https://solverstack.gitlabpages.inria.fr/pastix/), or provide PETSc.
Use LP64 BLAS/LAPACK and 32-bit indices for SCOTCH/PaStiX. PETSc must use
real double precision, 32-bit indices and Fortran bindings. Install a
[Gmsh SDK](https://gmsh.info/bin/Linux/) and set `MHDG_GMSH_DIR` to its root.

## Building

Review [arch.make](lib/Make.inc/arch.make) and the paths in the chosen
initialization script, then build from `lib/`:

```bash
cd lib
source Make.inc/init_vars_libs.sh
make -j
```

The default configuration builds serial, 2D `NGammaTiTeNeutral` with PaStiX
and copies the executable to `test/`. `arch.make` uses `pkg-config` for HDF5
(`hdf5_fortran`, or Debian's `hdf5` metadata), dense algebra and the selected
sparse solver. OpenBLAS is preferred when available; set
`MHDG_BLAS_PACKAGES='blas lapack'` to select those packages explicitly.
Gmsh's shared library is found under the SDK's `lib` or `lib64` directory.

For Pitagora DCGP, source
[init_vars_libs_pitagora.sh](lib/Make.inc/init_vars_libs_pitagora.sh) instead.
It selects `MODE=parall`, `PASTIX=no`, `PETSC=yes` and defaults to the Gmsh
SDK at `$HOME/libs/gmsh-4.14.1`. The local initializer uses a 4.11.1 SDK path.
Other environment examples are
[arch.make_meso](lib/Make.inc/arch.make_meso) and
[init_vars_libs_meso.sh](lib/Make.inc/init_vars_libs_meso.sh).

Explicit `make` arguments override build selections, for example
`make -j MODE=parall PASTIX=no PETSC=yes`. Use `make clean` before switching
model, serial/parallel mode, compiler, build type or solver backend; otherwise
rebuild with `make -j`. Source the same dependency environment in runtime jobs.
PETSc runtime settings are described under
[linear-solver selection](docs/simulation/configuration.md#linear-solver-selection).

## Running

The solver reads `param.txt` in its working directory. The executable and
`positionFeketeNodesTri2D.h5` are available in `test/`, a convenient location
for a basic run. Copy an appropriate parameter template to `param.txt`, then
set your equilibrium/current inputs, sources and an existing output directory.
Mesh and geometry files can be stored anywhere. Create `res/` inside the
run directory for temporary mesh files.
See [simulation configuration](docs/simulation/configuration.md) for template
selection, geometry preparation and the initialization workflow.

For a serial run from the repository root:

```bash
cd test
mkdir -p res
./MHDG-NGammaTiTeNeutral-serial-2D /path/to/mesh.msh
```

Supply a saved solution as the second argument to initialize or continue:

```bash
./MHDG-NGammaTiTeNeutral-serial-2D /path/to/mesh.msh checkpoint.h5
```

Use the parallel executable with your site's MPI launcher, for example:

```bash
mpirun -np 2 ./MHDG-NGammaTiTeNeutral-parall-2D /path/to/mesh.msh checkpoint.h5
```

Mesh and solution arguments accept full filenames or legacy stems. See
[restart and output](docs/simulation/restart-and-output.md) for checkpoint
interpretation and filenames, and [physical time stepping](docs/simulation/time-stepping.md)
for time units, final-time stopping and adaptive intervals.

## Testing and postprocessing

The [regression guide](regression_tests/README.md) and
[unit-test guide](unit_tests/README.md) describe verification and their
requirements. Regression workflows use external case bundles and include
warm/cold runs, MPI/OpenMP comparisons and golden references.

For visualization and example notebooks, see
[HDG_postprocess](https://github.com/wave46/HDG_postprocess).
