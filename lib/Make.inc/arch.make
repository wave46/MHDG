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
COMPTYPE = $(COMPTYPE_OPT)
#COMPTYPE = $(COMPTYPE_DEB)
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
FC = mpifort

#-------------------------------------------------------------------------------
# Model
#-------------------------------------------------------------------------------
MDL_LAPLACE=Laplace
MDL_NGAMMA=NGamma
MDL_NGAMMANEUTRAL=NGammaNeutral
MDL_NGAMMATITE=NGammaTiTe
MDL_NGAMMATITENEUTRAL=NGammaTiTeNeutral
MDL_NGAMMATITENEUTRALK=NGammaTiTeNeutralk
MDL_NGAMMAVORT=NGammaVort
# Model chosen
#MDL=$(MDL_NGAMMA)
#MDL=$(MDL_NGAMMANEUTRAL)
#MDL=$(MDL_NGAMMATITE)
MDL=$(MDL_NGAMMATITENEUTRAL)
#MDL=$(MDL_NGAMMATITENEUTRALK)
#MDL=$(MDL_NGAMMAVORT)
#MDL=$(MDL_LAPLACE)

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
#PETSC=$(LIB_YES)
PETSC=$(LIB_NO)


#-------------------------------------------------------------------------------
# MACROS FOR MODEL CHOICE
#-------------------------------------------------------------------------------
ifeq ($(MDL),$(MDL_NGAMMA))
 RMDL=NGamma
 MACROS+= -DNGAMMA
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMANEUTRAL))
 RMDL=NGamma
 MACROS+= -DNGAMMA
 MACROS+= -DNEUTRAL
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMATITE))
 RMDL=NGammaTiTe
 MACROS+= -DNGAMMA 
 MACROS+= -DTEMPERATURE
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMATITENEUTRAL))
 RMDL=NGammaTiTe
 #turns on continuity and momentum equations
 MACROS+= -DNGAMMA
 #turns on ion and electron energy equations
 MACROS+= -DTEMPERATURE
 #turns on neutral contimuity equations
 MACROS+= -DNEUTRAL
 #Uses AMJUEL splines for recombination and ionization, better turn on. If not, Te splines from NRL formulas employed
 MACROS+= -DAMJUELSPLINES
 #Uses AMJUEL splines for three-body recombination spline. If turned off, takes the one without three-body recombination, may be more stable
 MACROS+= -DTHREEBODYREC
 #Takes legacy approximated expression for cx rate. Left only for back-comparison, should not be usually used.
 #MACROS+= -DLEGACYCX
 #The combination of the two flags applies OpenADAS spline for thermal cx coefficient            
 MACROS+= -DEXPANDEDCX
 MACROS+= -DTHERMALCX
 #Applies soft min and max on neutral diffusion coefficient. Should be more stable
 MACROS+= -DDNNSMOOTH
 #Turns on linearization of neutral diffusion. Better to keep it on
 MACROS+= -DDNNLINEARIZED
 #Model with constant neutral diffusion. If turned on, turn off the previous two flags
 #MACROS+= -DCONSTANTNEUTRALDIFF
 #Actually not saves, but monitors in the output particle balance
 MACROS+= -DSAVEFLUX
 #The following 4 flags are development ones, should not be used
 #MACROS+= -DBOHMLIMIT
 #MACROS+= -DNEUTRALCONVECTION
 #MACROS+= -DPINCH
 #MACROS+= -DRHSBC
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMATITENEUTRALK))
 RMDL=NGammaTiTe
 MACROS+= -DNGAMMA 
 MACROS+= -DTEMPERATURE
 MACROS+= -DNEUTRAL
 MACROS+= -DAMJUELSPLINES
 MACROS+= -DTHREEBODYREC
 #MACROS+= -DRHSBC
 MACROS+= -DDNNSMOOTH
 MACROS+= -DDNNLINEARIZED
 #MACROS+= -DNEUTRALPUMP
 #MACROS+= -DLEGACYCX
 MACROS+= -DSAVEFLUX
 MACROS+= -DEXPANDEDCX
 MACROS+= -DTHERMALCX
 #MACROS+= -DCONSTANTNEUTRALDIFF
 MACROS+= -DKEQUATION
 #MACROS+= -DBOHMLIMIT
 #MACROS+= -DNEUTRALCONVECTION
 #MACROS+= -DDKLINEARIZED
 #MACROS+= -DPINCH
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMAVORT))
 RMDL=NGammaVort
 MACROS+= -DNGAMMA 
 MACROS+= -DVORTICITY
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_LAPLACE))
 RMDL=Laplace
 ADDMOD+=hdg_LimitingTechniques.o
else ifeq ($(MDL),$(MDL_NGAMMALAPLACE))
 RMDL=NGammaLaplace
 MACROS+= -DNGAMMA
 ADDMOD+=hdg_LimitingTechniques.o 
else
 abort Unsupported MDL==$(MDL)
 exit
endif


#-------------------------------------------------------------------------------
# MACROS FOR SERIAL/PARALLEL
#-------------------------------------------------------------------------------
ifeq ($(MODE),$(MODE_SERIAL))
else ifeq ($(MODE),$(MODE_PARALL))
  MACROS+= -DPARALL
  ADDMOD+=Communications.o  
endif

ifeq ($(DIM),$(DIM_3D))
  MACROS+= -DTOR3D
endif

#-------------------------------------------------------------------------------
# MACROS FOR LINEAR SYSTEM LIBRARIES
#-------------------------------------------------------------------------------
ifeq ($(PASTIX),$(LIB_YES))
 MACROS+= -DWITH_PASTIX
 ADDMODSOLV+=solve_pastix.o
 ADDMOD+=solve_pastix.o
endif
ifeq ($(PETSC),$(LIB_YES))
 MACROS+= -DWITH_PETSC
 ADDMODSOLV+=solve_petsc.o
 ADDMOD+=solve_petsc.o
endif
ifeq ($(PSBLMG),$(LIB_YES))
 MACROS+= -DWITH_PSBLAS
 ADDMODSOLV+=solve_psblas.o
 ADDMOD+=solve_psblas.o
else ifeq ($(PSBLAS),$(LIB_YES))
 MACROS+= -DWITH_PSBLAS
 ADDMODSOLV+=solve_psblas.o
 ADDMOD+=solve_psblas.o
endif

#-------------------------------------------------------------------------------
# MACROS FOR COMPILING OPTIONS
#-------------------------------------------------------------------------------
####### Begin gfortran #######
ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
 FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid
 FCFLAGS += -Wall -Wextra -Wconversion -fcheck=all -Wuninitialized -Wtabs
else ifeq ($(COMPTYPE),$(COMPTYPE_PRO))
 FCFLAGS = -Og -pg
else ifeq ($(COMPTYPE),$(COMPTYPE_OPT))
 FCFLAGS = -Ofast
endif

FCFLAGS += -cpp  -fopenmp
FCFLAGS += -fdefault-double-8 -fdefault-real-8 
FCFLAGS += -ffree-line-length-none -fimplicit-none -ffree-form -Wno-tabs
######## End gfortran ########

######## Begin ifort #########
#ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
# FCFLAGS = -fpp -O0 -r8 -qopenmp -xHOST -g -traceback -check all -free -check bounds -debug all
#else ifeq ($(COMPTYPE),$(COMPTYPE_PRO))
# FCFLAGS = -fpp -O3 -parallel -fpe0 -qopenmp -xHOST -r8 -free -g
#else ifeq ($(COMPTYPE),$(COMPTYPE_OPT))
# FCFLAGS = -fpp -O3 -parallel -fpe0 -qopenmp -xHOST -r8 -free
#endif
######### End ifort ##########

DEF = -DTHREAD_FUNNELED

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
# HDF5/HWLOC/X11
#Local
FCFLAGS += -I/usr/include
FCFLAGS += -I/usr/include/x86_64-linux-gnu
FCFLAGS += -I/usr/include/hdf5/serial
FCFLAGS += -I/usr/include/hwloc
FCFLAGS += -I/usr/include/X11

#GMSH/MMG
FCFLAGS += -I$(MHDG_MMG_DIR)/build/include
FCFLAGS += -I$(MHDG_GMSH_DIR)/include

# PASTIX
ifeq ($(PASTIX),$(LIB_YES))
 FCFLAGS += $(shell echo `PKG_CONFIG_PATH=${PKG_CONFIG_PATH} pkg-config --cflags pastix pastixf`)
 FCFLAGS += -I$(MHDG_SCOTCH_DIR)/include
endif

# PSBLAS
ifeq ($(PSBLAS),$(LIB_YES))
 FCFLAGS += -I$(MHDG_PSBLAS_DIR)/include/
 FCFLAGS += -I$(MHDG_PSBLAS_DIR)/modules/
endif

# PETSC
ifeq ($(PETSC),$(LIB_YES))
 FCFLAGS += -I$(MHDG_PETSC_DIR)/include/
 FCFLAGS += -I$(MHDG_PETSC_DIR)/$(PETSC_ARCH)/include
endif

# MLD2P4
ifeq ($(PSBLMG),$(LIB_YES))
 FCFLAGS += -I$(MHDG_MLD2P4_DIR)/modules/
 FCFLAGS += -I$(MHDG_MLD2P4_DIR)/include/
endif

#-------------------------------------------------------------------------------
# Libraries needed for linking
#-------------------------------------------------------------------------------
# HDF5/HWLOC/X11
#Local
LIB += -L/usr/lib/x86_64-linux-gnu -lz -lm -lrt -lpthread
LIB += -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5
LIB += -L/usr/lib/x86_64-linux-gnu/hwloc -lhwloc
LIB += -L/usr/lib/x86_64-linux-gnu -lX11
LIB += -L/usr/lib/x86_64-linux-gnu/xtables -lXt

#GMSH/MMG
LIB += -L$(MHDG_GMSH_DIR)/lib -Llib -lgmsh -L. -Wl,-rpath=$(MHDG_GMSH_DIR)/lib 
LIB += -L$(MHDG_MMG_DIR)/build/lib -lmmg


# PASTIX
ifeq ($(PASTIX),$(LIB_YES))
 #former: lifcore is for intel
 #LIB += -L$(HOME)/libs/scotch_6.0.4/lib/ -lscotch -lscotcherrexit  -lptscotchparmetis -lptscotch -lpthread -lhwloc
 #LIB += -L$(MHDG_PASTIX_DIR)/install -lpastix -lm -lrt -lifcore
 #New GNU
 LIB += -L$(MHDG_SCOTCH_DIR)/lib -lptscotch -lscotch -lptscotcherr -lz -lm -lrt -lpthread
 LIB += $(shell echo `PKG_CONFIG_PATH=${PKG_CONFIG_PATH} pkg-config --libs pastix pastixf`)
 #New INTEL
 #LIB += -L$(MHDG_SCOTCH_DIR)/lib -lptscotch -lscotch -lptscotcherr -lz -lm -lrt -lpthread
 #LIB += -L$(MHDG_PASTIX_DIR)/install -lpastix -lm -lrt -lifcore -lpthread -lhwloc -lptscotch -lscotch -lscotcherr
endif

# BLAS/LAPACK
#Local
LIB += -L/usr/lib/x86_64-linux-gnu -lblas -llapack -llapacke


# PSBLAS/MLD2P4
ifeq  ($(PSBLMG),$(LIB_YES)) 
 LIB += -L$(MHDG_PSBLAS_DIR)/lib/ -L$(MHDG_MLD2P4_DIR)/lib/ 
 LIB += -lpsb_krylov -lmld_prec -lpsb_prec -lpsb_krylov -lpsb_prec -lpsb_util -lpsb_base
else ifeq ($(PSBLAS),$(LIB_YES)) 
 LIB += -L$(MHDG_PSBLAS_DIR)/lib/
 LIB += -lpsb_util -lpsb_krylov -lpsb_prec -lpsb_base
endif

ifeq ($(PSBLMG),$(LIB_YES))
 include $(MHDG_MLD2P4_DIR)/include/Make.inc.mld2p4
else ifeq ($(PSBLAS),$(LIB_YES))
 include $(MHDG_PSBLAS_DIR)/include/Make.inc.psblas
endif


# PETSC
#Local
ifeq ($(PETSC),$(LIB_YES))
 LIB += -L$(MHDG_PETSC_DIR)/lib -lpetsc
 LIB += -L$(MHDG_PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc
endif




#-------------------------------------------------------------------------------
# Add macros to the compiling flags
#-------------------------------------------------------------------------------
FCFLAGS += $(MACROS)

#-------------------------------------------------------------------------------
# Utilitaires
#-------------------------------------------------------------------------------
RM=/bin/rm -f
CP=/bin/cp -f
CMP=cmp
