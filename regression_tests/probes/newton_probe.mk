# Loaded after the solver Makefile in a private copy of its build directory.
# The normal link recipe already includes LIB.
override LIB += newton_linear_probe.o

Main_utils.o: newton_linear_probe.o

newton_linear_probe.o: newton_linear_probe.f90 hdg_Residual.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC) $< -o $@
