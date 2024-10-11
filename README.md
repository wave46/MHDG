$\kappa\text{-}\epsilon$

## Dependencies
### Debian
```zsh
sudo apt-get install gfortran openmpi-bin openmpi-common openmpi-doc libopenmpi-dev libblas-dev liblapack-dev libxt-dev libhdf5-serial-dev libscotch-dev libptscotch-dev
```
Then install [Pastix](https://gitlab.inria.fr/solverstack/pastix)
### Macos
```zsh
brew install hdf5 scotch
git clone https://gitlab.inria.fr/solverstack/pastix.git 
cd pastix/tools/homebrew
brew install --build-from-source pastix.rb
```
## Building
```zsh
cd lib
source Make.inc/init_vars_libs.sh
make
```
On macos, change `arch.make` to `macos.make` in the `Makefile`. Modify `macos` to your needs.
## Running
Using the circular case, do small time steps at high diffusion, prepare the `param.txt`, accordingly.
Run the executable in the same foder as `param.txt`.
```zsh
test/MHDG-NGammaTiTeNeutral-serial-2D PATH_TO_MESH/CircLimAlign_Quads_Nel588_P6
```
Supposing the results are stored in the `init` folder, modify or created a new `param.txt` file such that the diffusion is decreased during the simulation. The simulation runs until it crashes at the lowest diffusion possible.
```zsh
test/MHDG-NGammaTiTeNeutral-serial-2D PATH_TO_MESH/CircLimAlign_Quads_Nel588_P6 init/the_shortest_filename
```
If running the WEST case, rebuild in parallel mode by modifying the `lib/Make.inc/arch.make` file. Change the testcase in `param.txt`. Run the executable in the same folder as `positionFeketeNodesTri2D.h5`. Do the 2 step process like the circular case. The meshe and the ouput is partitioned into 8 files.
```zsh
mpirun -n 8 test/MHDG-NGammaTiTeNeutral-parall-2D PATH_TO_MESHES/West_NoHole_Nel13118_P
```
## Postprocessing
For visualizing the results, see this [python package](https://github.com/wave46/HDG_postprocess).


