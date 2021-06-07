#-------------------------------------------------------------------------------
# Source directory
#-------------------------------------------------------------------------------
SDIR=$(PWD)/../src/

#-------------------------------------------------------------------------------
# Compiling type: optimized - debugging - profiling
#-------------------------------------------------------------------------------
COMPTYPE_OPT = opt
COMPTYPE_DEB = deb
COMPTYPE_PRO = pro
#COMPTYPE = $(COMPTYPE_DEB)
COMPTYPE = $(COMPTYPE_OPT)
#COMPTYPE = $(COMPTYPE_PRO)

#-------------------------------------------------------------------------------
# Mode: serial or parallel
#-------------------------------------------------------------------------------
MODE_SERIAL = serial
MODE_PARALL = parall
MODE = $(MODE_SERIAL)
#MODE = $(MODE_PARALL)

#-------------------------------------------------------------------------------
# The compiler
#-------------------------------------------------------------------------------
FC = mpiifort

#-------------------------------------------------------------------------------
# Model
#-------------------------------------------------------------------------------
# Available models
MDL_NGAMMA=NGamma
MDL_NGAMMANEUTRAL=NGammaNeutral
MDL_NGAMMATITE=NGammaTiTe
MDL_NGAMMATITENEUTRAL=NGammaTiTeNeutral

# Model chosen
#MDL=$(MDL_NGAMMA)
#MDL=$(MDL_NGAMMANEUTRAL)
MDL=$(MDL_NGAMMATITE)
#MDL=$(MDL_NGAMMATITENEUTRAL)

# Moving equilibrium
MVEQ_TRUE = true
MVEQ_FALS = false
MVEQ = $(MVEQ_FALS)
#MVEQ = $(MVEQ_TRUE)

#-------------------------------------------------------------------------------
# Dimensions 2D/3D
#-------------------------------------------------------------------------------
DIM_3D=3D
DIM_2D=2D
DIM=$(DIM_2D)

#-------------------------------------------------------------------------------
# Libraries for linear system solver
#-------------------------------------------------------------------------------
# Available libraries: put $(LIB_YES) to use the library, $(LIB_NO) to not use it
LIB_YES=yes
LIB_NO=no
PASTIX=$(LIB_YES)
#PASTIX=$(LIB_NO)
#PSBLAS=$(LIB_YES)
PSBLAS=$(LIB_NO)
#PSBLMG=$(LIB_YES)
PSBLMG=$(LIB_NO)


# MACROS FOR MODEL CHOICE
ifeq ($(MDL),$(MDL_NGAMMA))
 RMDL=NGamma
 MACROS$+= -DNGAMMA 
ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMANEUTRAL))
 RMDL=NGamma
 MACROS$+= -DNGAMMA
 MACROS$+= -DNEUTRAL
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMATITE))
 RMDL=NGammaTiTe
 MACROS+= -DN4EQ
 MACROS+= -DNGAMMA 
 MACROS+= -DENERGY
 MACROS+= -DTEMPERATURE
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMATITENEUTRAL))
 RMDL=NGammaTiTe
 MACROS+= -DN4EQ
 MACROS+= -DNGAMMA 
 MACROS+= -DENERGY
 MACROS+= -DTEMPERATURE
 MACROS+= -DNEUTRAL
 ADDMOD+=hdg_LimitingTechniques.o 
else
 abort Unsupported MDL==$(MDL)
 exit
endif



# MACROS FOR MOVING EQUILIBRIUM
ifeq ($(MVEQ),$(MVEQ_FALS))
#Do Nothing
else ifeq ($(MVEQ),$(MVEQ_TRUE))
 MACROS+= -DMOVINGEQUILIBRIUM
 #ADDMOD+=hdg_MagneticDependingMatrices.o
else
 abort Unsupported MVEQ==$(MVEQ)
 exit
endif


# MACROS FOR SERIAL/PARALLEL
ifeq ($(MODE),$(MODE_SERIAL))
else ifeq ($(MODE),$(MODE_PARALL))
  MACROS+= -DPARALL
  ADDMOD+=Communications.o  
endif

ifeq ($(DIM),$(DIM_3D))
  MACROS+= -DTOR3D
endif

# MACROS FOR LINEAR SYSTEM LIBRARIES
ifeq ($(PASTIX),$(LIB_YES))
 MACROS+= -DWITH_PASTIX
 ADDMOD+=solve_pastix.o
endif
ifeq ($(PSBLMG),$(LIB_YES))
 MACROS+= -DWITH_PSBLAS
 MACROS+= -DWITH_MLD2P4
 ADDMOD+=solve_psblas.o
else ifeq ($(PSBLAS),$(LIB_YES))
 MACROS+= -DWITH_PSBLAS
 ADDMOD+=solve_psblas.o
endif


# MACROS FOR COMPILING OPTIONS
ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
 FCFLAGS = -O0 -g -traceback -check bounds
else ifeq ($(COMPTYPE),$(COMPTYPE_PRO))
 FCFLAGS =  -O2 -xHost
else ifeq ($(COMPTYPE),$(COMPTYPE_OPT))
 FCFLAGS =  -O2 -xHost 
endif


FCFLAGS += -fpp -qopenmp -xCORE-AVX512 -mtune=skylake -fpe0
FCFLAGS += -r8
FCFLAGS += -free

DEF = -DTHREAD_FUNNELED


# MKL root folder
MKL_ROOT = /trinity/shared/apps/tr17.10/x86_64/intel-2018.0.128/compilers_and_libraries_2018.0.128/linux/mkl
# BLAS and LAPACK libraries
LIBBLAS = $(MKL_ROOT)/lib/intel64/libmkl_intel_lp64.a $(MKL_ROOT)/lib/intel64/libmkl_sequential.a $(MKL_ROOT)/lib/intel64/libmkl_core.a

# Includes
FCFLAGS += -I/trinity/shared/apps/tr17.10/x86_64/hdf5-icc18-impi-1.10.1/include/
FCFLAGS += -I/home/ptamain/libs/new_meso/pastix_5.2.2.22_64bits_MKL_scotch_6.0.5a_ifort18/src/../install/

# Libraries needed for linking
LIB += -L/home/ptamain/libs/new_meso/pastix_5.2.2.22_64bits_MKL_scotch_6.0.5a_ifort18/src/../install/ -lpastix -lm -lrt -lifcore
LIB += -L/home/ptamain/libs/new_meso/scotch_6.0.5a_64bits_ifort18/lib -lptscotch -lscotch -lptscotcherrexit
LIB += -L/trinity/shared/apps/tr17.10/x86_64/hwloc-1.11.8/lib -lhwloc -lpthread
LIB += $(LIBBLAS) $(LIBBLAS)
LIB += -L/trinity/shared/apps/tr17.10/x86_64/hdf5-icc18-impi-1.10.1/lib -lhdf5_fortran -lhdf5 -lz -L/usr/lib64 -lX11

# Add macros to the compiling flags
FCFLAGS += $(MACROS)
#-------------------------------------------------------------------------------
# Utilitaires
#-------------------------------------------------------------------------------

RM=/bin/rm -f
CP=/bin/cp -f
CMP=cmp
