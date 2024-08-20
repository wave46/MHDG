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
# COMPTYPE = $(COMPTYPE_DEB)
# COMPTYPE = $(COMPTYPE_PRO)

#-------------------------------------------------------------------------------
# Mode: serial or parallel
#-------------------------------------------------------------------------------
MODE_SERIAL = serial
MODE_PARALL = parall
# MODE = $(MODE_SERIAL)
MODE = $(MODE_PARALL)

#-------------------------------------------------------------------------------
# The compiler
#-------------------------------------------------------------------------------
FC = mpif90 
#FC = mpif90 -f90=ifort

#-------------------------------------------------------------------------------
# Model
#-------------------------------------------------------------------------------
MDL_LAPLACE=Laplace
MDL_NGAMMA=NGamma
MDL_NGAMMANEUTRAL=NGammaNeutral
MDL_NGAMMATITE=NGammaTiTe
MDL_NGAMMATITENEUTRAL=NGammaTiTeNeutral
MDL_NGAMMAVORT=NGammaVort
# Model chosen
#MDL=$(MDL_NGAMMA)
#MDL=$(MDL_NGAMMANEUTRAL)
#MDL=$(MDL_NGAMMATITE)
MDL=$(MDL_NGAMMATITENEUTRAL)
#MDL=$(MDL_NGAMMAVORT)
#MDL=$(MDL_LAPLACE)

#-------------------------------------------------------------------------------
# Moving equilibrium (not available anymore for this version)
#-------------------------------------------------------------------------------
#MVEQ_TRUE = true
#MVEQ_FALS = false
#MVEQ = $(MVEQ_FALS)
##MVEQ = $(MVEQ_TRUE)

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

# pkg-config path, change the path or add if needed
PKG_CONFIG_PATH:=${PKG_CONFIG_PATH}:/usr/local/opt/lapack/lib/pkgconfig
MHDG_SCOTCH_DIR=/usr/local/opt/scotch




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
 MACROS+= -DNGAMMA 
 MACROS+= -DTEMPERATURE
 MACROS+= -DNEUTRAL
 MACROS+= -DAMJUELSPLINES
 MACROS+= -DTHREEBODYREC
 #MACROS+= -DRHSBC
 MACROS+= -DDNNSMOOTH
 MACROS+= -DDNNLINEARIZED
 #MACROS+= -DLEGACYCX
 MACROS+= -DSAVEFLUX
 MACROS+= -DEXPANDEDCX
 MACROS+= -DTHERMALCX
 #MACROS+= -DBOHMLIMIT
 #MACROS+= -DCONSTANTNEUTRALDIFF
 #MACROS+= -DNEUTRALCONVECTION
 #MACROS+= -DPINCH
 MACROS+= -DKEQUATION
 #MACROS+= -DKDIFFSMOOTH
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
# MACROS FOR MOVING EQUILIBRIUM
#-------------------------------------------------------------------------------
#ifeq ($(MVEQ),$(MVEQ_FALS))
## Do nothing
#else ifeq ($(MVEQ),$(MVEQ_TRUE))
# MACROS+= -DMOVINGEQUILIBRIUM
# ADDMOD+=hdg_MagneticDependingMatrices.o
#else
# abort Unsupported MVEQ==$(MVEQ)
# exit
#endif


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
 MACROS+= -DWITH_PASTIX -DNEW_PASTIX
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

#-------------------------------------------------------------------------------
# MACROS FOR COMPILING OPTIONS
#-------------------------------------------------------------------------------
####### Begin gfortran #######
ifeq ($(COMPTYPE),$(COMPTYPE_DEB))
FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check #-ffpe-trap=zero,overflow,underflow,invalid
 # FCFLAGS = -Og -g -fbounds-check -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,invalid
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
FCFLAGS += `pkg-config --cflags hdf5_fortran`

# PASTIX
ifeq ($(PASTIX),$(LIB_YES))
 # FCFLAGS += -I$(MHDG_PASTIX_DIR)/install
 # FCFLAGS += $(shell pkg-config --cflags $(MHDG_PASTIX_DIR)/spm.pc $(MHDG_PASTIX_DIR)/pastix.pc)
 FCFLAGS += `pkg-config --cflags pastix pastixf`
 FCFLAGS += -I$(MHDG_SCOTCH_DIR)/include
endif

# PSBLAS
ifeq ($(PSBLAS),$(LIB_YES))
  FCFLAGS += `PKG_CONFIG_PATH=${PKG_CONFIG_PATH} pkg-config --cflags blas lapack lapacke`
 # FCFLAGS += -I$(MHDG_PSBLAS_DIR)/include/
 # FCFLAGS += -I$(MHDG_PSBLAS_DIR)/modules/
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

LIB += `pkg-config --libs hdf5_fortran` -lpthread -lz -lm
# macos, change the path if needed
LIB += -L/Library/Developer/CommandLineTools/SDKs/MacOSX12.1.sdk/usr/lib
# LIB +=  -Wl,-rpath,/usr/local/lib



# PASTIX
ifeq ($(PASTIX),$(LIB_YES))
 LIB += -L$(MHDG_SCOTCH_DIR)/lib -lptscotch -lscotch -lptscotcherr
 LIB += `pkg-config --libs pastix pastixf`
 # LIB += -L$(MHDG_PASTIX_DIR)/install -lpastix -lhwloc
endif

# BLAS/LAPACK
LIB += `PKG_CONFIG_PATH=${PKG_CONFIG_PATH} pkg-config --libs blas lapack lapacke`

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
