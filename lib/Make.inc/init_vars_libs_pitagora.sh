#!/usr/bin/env bash
# Source this file before building or running MHDG on Pitagora DCGP.

module purge || return 1
module load profile/base || return 1
module load openmpi/4.1.6--gcc--12.3.0-ucx1.20 || return 1
module load petsc/3.22.1--openmpi--4.1.6--gcc--12.3.0-ucx1.20-mumps || return 1
module load hdf5/1.14.3--openmpi--4.1.6--gcc--12.3.0-ucx1.20 || return 1

# Expose all include/library paths in metadata queries, even when modules have
# already added them to CPATH or LIBRARY_PATH.
export PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1
export PKG_CONFIG_ALLOW_SYSTEM_LIBS=1
pkg-config --print-errors --exists PETSc openblas hdf5_fortran || return 1

# Build choices consumed by arch.make; command-line make settings take priority.
export MODE=parall
export PASTIX=no
export PETSC=yes
export MHDG_BLAS_PACKAGES=openblas

# The only private installation; lib/lib64 is detected by arch.make.
export MHDG_GMSH_DIR="${MHDG_GMSH_DIR:-$HOME/libs/gmsh-4.14.1}"
