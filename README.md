## Dependencies
This is an example to build on local PC with scotch-pastix solver and gmsh for mesh refinement.
To run the code one needs a fortran compiler, openmpi (or other MPI), BLAS, LAPACK, LAPACKE, Xt, HDF5, SCOTCH+PASTIX (or PETSc), GMSH libraries
### Debian
```zsh
sudo apt-get install gfortran openmpi-bin openmpi-common openmpi-doc libopenmpi-dev libblas-dev liblapacke-dev  liblapack-dev libxt-dev libhdf5-serial-dev 
```
Then install [SCOTCH](https://gitlab.inria.fr/scotch/scotch) [Pastix](https://solverstack.gitlabpages.inria.fr/pastix/md_docs_doxygen_chapters_Pastix_MPI.html). Install int32 versions of both.
Next, install and load gmsh 4.11.1 SDK [pack](https://gmsh.info/bin/Linux/) (only tested with this gmsh version), you can also install gmsh from source code

## Building
Make sure that you specify properly the paths to needed libraries.
By defualt `arch.make` file contains instructions for local build with gfortran, OpenMPI, SCOTCH+PASTIX, GMSH libraries.
`arch.make_meso` contains instructions to build on mesocentre cluster with intel 2020 compiler. `init_vars_libs_meso.sh` is another example of environment setup
```zsh
cd lib
source Make.inc/init_vars_libs.sh
make
```

## Regression testing

The tracked regression harness supports external case bundles, isolated
warm/cold workflows, MPI/OpenMP layouts, and golden-reference comparison. See
the [regression test guide](regression_tests/README.md).

## Mesh generation

The machine wall outline is needed with knowing the locations of puff and pump.
Then a `.geo` file is needed to make `.msh` files using GMSH GUI.
`.geo` should be put in `test/res/geometries` folder, `.msh` files to `test/Meshes`
The refinement level for initial mesh should be of order of 1k P8 elements for a machine like TCV, then it will be automatically refined if mesh adaptivity is used.
Without mesh adaptivity typicla P8 mesh for WEST tokamak has ~15-20k elements refined at the divertor targets and at the wall in the far SOL

## Input configuration

The solver reads `param.txt` from the working directory. The four templates in
`test/` use the same input keys: `param_initial.txt` starts from a mesh at high
diffusion, `param_diffred.txt` reduces diffusion through pseudo-steady states,
`param_steady.txt` reconverges at fixed diffusion, and `param_reference.txt`
provides an adaptive steady-state example with separate particle and heat
diffusivities. Set paths, grid dimensions, and case-specific values before use.

The neutral limiter temperature source is `'ti'` (local ion temperature) or
`'fixed'` (requires `neutral_flux_limiter_tn_eV > 0`). Flux regularization and
the minimum cap flux are specified in `m^-2 s^-1`. Set
`neutral_flux_limiter_save_2d = .true.` in `SWITCH_LST` to save nodal limiter
diagnostics; this can also be used with the limiter off.

Enabled balance diagnostics require a 2D neutral-temperature model and
`bohmtypebc = 0`. `summary` reports compact content, inputs/outputs, and total
balances; `equations` reports individual equation balances; `detailed` adds
physical, discrete, and boundary-condition components. Physical wall reporting
supports Bohm, BohmPump, and BohmPuff boundaries; other physical boundaries
produce a warning.

`INPUT_LST` selects the auxiliary namelists, with paths relative to the working
directory unless absolute paths are supplied:

- `transport_model_path = 'transport_model.nml'` supplies
  `TRANSPORT_MODEL_1D_LST` when `transport_1d` is enabled. See the
  [sample](test/transport_model.nml) and
  [transport guide](src/Models/NGammaTiTe/transport_1d/README.md) for topology
  policies, particle-diffusion taper, and pinch settings.
- `impurity_model_path = 'impurity_model.nml'` supplies
  `IMPURITY_RADIATION_LST` and is required when `impurity_radiation` is enabled.
  The [sample](test/impurity_model.nml) includes every supported key. Species
  names (`'N'`, `'W'`) and non-negative concentrations are arrays, with
  concentrations expressed as fractions of electron density. Mixture settings
  and the moving-equilibrium concentration-file controls belong in this file.

The feedback gains apply to moving-equilibrium control. The spelling
`feedback_propotional_gain` (and its `_xpr` counterpart) matches the reader.

## Running 
The pipeline to get solution for a case without meshadaptivity and in serial is described in [demo](https://github.com/wave46/HDG_postprocess/blob/main/demos/hdg_solution_basics_neutrals.ipynb)
First, go to `test` folder, where executable is. It is convinient to create `Meshes` folder and put your mesh there.
Take the `param_initial.txt` and copy its contents to `param.txt`. Do not forget to specify paths to files with magnetic field and plasma current, as well as the dimensions of the grid and output folder where solutions will be saved
To initialize simulation one should start with high diffusion (usually arond 20m^2/s) and several small timesteps (`dt0=1e3` in adimensional values).
Run the following command
```zsh
./MHDG-NGammaTiTeNeutral-serial-2D ./Meshes/Name_of_the_mesh_without_.msh_but_with_P{polynomial_order}_at_the_end
```
Having first initial guess for high diffusion (in a file with shortest name, without _NR00, _000 and so on), now we can reduce diffusion to the desired values (usually around 1m^2/s). This is done automatically in the code.
One copy `param_diffred.txt` to `param.txt`. For each diffusion value a (pseudo)steady state will be achieved and then the diffusion value will be multiplied by `diffred` parameter and then a steady state for lower diffusion will be found.
Usually simulation runs until it crashes for very low diffusion, or you can specify number of (pseudo)timesteps (i.e. number of (pseudo)steady states will be found for decreasing diffusion values) with parameter `nts` (10 in the example).
Or you can also set the value of minimal desired diffusion `diffmin`.
To restart your simulation from initialization solution from previous step, run 
```zsh
./MHDG-NGammaTiTeNeutral-serial-2D ./Meshes/Name_of_the_mesh_without_.msh_but_with_P{polynomial_order}_at_the_end /path/to/initial/solution
```

Quite often the last diffusion value from previous solution is not exactly what you needed. For example, you started from 20 m^2/s, `diffred`=0.6 and you wished to have `diffmin`=0.5. 
In this case your closest diffusion value will be 0.56 m^2/s. Also it was not a "true" steady state solution, but a solution with very big timestep. 
So it's always convenient to find a base solution for your future scan for precise value of diffusion and in steady state.
Copy `param_steady.txt` to `param.txt` and rerun from soltion with closest diffusion to desired one:

```zsh
./MHDG-NGammaTiTeNeutral-serial-2D ./Meshes/Name_of_the_mesh_without_.msh_but_with_P{polynomial_order}_at_the_end /path/to/solution/with/0.560E+00/in/name
```

From now on you can conduct your scans, for example, puff scan (`puff` in parameters), various diffusions (`diff_n`, `diff_u`, `diff_e`, `diff_ee`), recycling (`R`), etc... 
Sometimes if you change parameters too much, the simulation may "crash", i.e. Newton-Raphson algorithm will not converge. 
You may then change the desired values in scan less or you can also play with minimal and maximal Newton-Raphson damping factors (dumpnr_min,dumpnr_max, dumpnr_width, dumpnr_n0).

## Postprocessing
For visualizing the results, see this [python package](https://github.com/wave46/HDG_postprocess).
