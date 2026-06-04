#---------------------------------------
# project: MHDG
# Makefile
# date: 12/02/2016
#---------------------------------------
include Make.inc/arch.make

.PHONY : clean

all : allseq

allseq : MHDG

#allpar : checkmodepar plato geninitsol
#	@echo 'Makefile done'

default : MHDG

MHDG : MHDG-$(MDL)-$(MODE)-$(DIM)
	@echo ' '
	@echo '===================================================='
	@echo ' Executable done : $<'
	@echo '===================================================='
	@echo ' '

Convergence : Convergence-$(MDL)-$(DIM)
	@echo ' '
	@echo '===================================================='
	@echo ' Executable done : $<'
	@echo '===================================================='
	@echo ' '


clean:
	@echo "Cleaning *.o *.mod *.F"
	$(RM) *.o *.mod *.kmo  *.F *.f90 *.F90
	$(RM) MHDG-* Convergence *~

mod_splines.F : $(SDIR)/Spline/mod_splines.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

mod_splines.o : mod_splines.F LinearAlgebra.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

prec_const.F : $(SDIR)/Definitions/prec_const.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

prec_const.o : prec_const.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

types.F : $(SDIR)/Definitions/types.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

types.o : types.F prec_const.o mod_splines.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

matrices_types.F : $(SDIR)/LinearAlgebra/matrices_types.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

matrices_types.o : matrices_types.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

globals.F : $(SDIR)/Definitions/globals.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

globals.o : globals.F types.o matrices_types.o prec_const.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

MPI_OMP.F : $(SDIR)/MPI_OMP/MPI_OMP.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

MPI_OMP.o : MPI_OMP.F globals.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

magnetic_field.F : $(SDIR)/Models/magnetic_field.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

magnetic_field.o : magnetic_field.F globals.o prec_const.o HDF5_io.o MPI_OMP.o interpolation.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

external_heating.F : $(SDIR)/Models/external_heating.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@
external_heating.o : external_heating.F globals.o inout.o prec_const.o HDF5_io.o MPI_OMP.o interpolation.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

#transport_coefficient_map.F : $(SDIR)/Models/transport_coefficient_map.f90
#	$(CMP) -s $< $@ || $(CP) -f $< $@

#transport_coefficient_map.o : transport_coefficient_map.F globals.o
#	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

HDF5_io.F : $(SDIR)/InOut/HDF5_io.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

HDF5_io.o : HDF5_io.F prec_const.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

gmsh_io.F : $(SDIR)/InOut/gmsh_io.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

gmsh_io.o : gmsh_io.F reference_element.o MPI_OMP.o HDF5_io.o reference_element.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

matrices_tools.F : $(SDIR)/LinearAlgebra/matrices_tools.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

matrices_tools.o : matrices_tools.F matrices_types.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

PrintUtils.F : $(SDIR)/Utils/PrintUtils.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

PrintUtils.o : PrintUtils.F matrices_types.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics.F : $(SDIR)/Utils/Diagnostics/diagnostics.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics.o : diagnostics.F MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics_accumulation.F : $(SDIR)/Utils/Diagnostics/diagnostics_accumulation.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics_accumulation.o : diagnostics_accumulation.F diagnostics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics_reduction.F : $(SDIR)/Utils/Diagnostics/diagnostics_reduction.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics_reduction.o : diagnostics_reduction.F diagnostics.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics_quantities.F : $(SDIR)/Utils/Diagnostics/diagnostics_quantities.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics_quantities.o : diagnostics_quantities.F diagnostics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics_wall_sources.F : $(SDIR)/Utils/Diagnostics/diagnostics_wall_sources.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics_wall_sources.o : diagnostics_wall_sources.F diagnostics.o globals.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics_output_terminal.F : $(SDIR)/Utils/Diagnostics/diagnostics_output_terminal.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics_output_terminal.o : diagnostics_output_terminal.F diagnostics.o diagnostics_quantities.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

diagnostics_output_hdf5.F : $(SDIR)/Utils/Diagnostics/diagnostics_output_hdf5.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

diagnostics_output_hdf5.o : diagnostics_output_hdf5.F diagnostics.o diagnostics_quantities.o diagnostics_wall_sources.o HDF5_io.o globals.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

Communications.F : $(SDIR)/MPI_OMP/Communications.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Communications.o : Communications.F MPI_OMP.o preprocess.o PrintUtils.o globals.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

interpolation.F : $(SDIR)/LinearAlgebra/interpolation.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

interpolation.o : interpolation.F PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

adimensionalization.F : $(SDIR)/Models/adimensionalization.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adimensionalization.o : adimensionalization.F globals.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

LinearAlgebra.F : $(SDIR)/LinearAlgebra/LinearAlgebra.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

LinearAlgebra.o : LinearAlgebra.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

inout.F : $(SDIR)/InOut/inout.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

inout.o : inout.F HDF5_io.o flux_surface_transport_data.o transport_models_1d.o diagnostics.o diagnostics_output_hdf5.o globals.o PrintUtils.o LinearAlgebra.o gmsh_io.o MPI_OMP.o Communications.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

Debug.F : $(SDIR)/Utils/Debug.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Debug.o : Debug.F globals.o PrintUtils.o inout.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

physics.F : $(SDIR)/Models/$(RMDL)/physics.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

physics.o : physics.F globals.o magnetic_field.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

ifeq ($(PASTIX),$(LIB_YES))
solve_pastix.F : $(SDIR)/LinearAlgebra/solve_pastix.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_pastix.o : solve_pastix.F matrices_types.o matrices_tools.o globals.o types.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
endif

ifeq ($(PETSC),$(LIB_YES))
solve_petsc.F : $(SDIR)/LinearAlgebra/solve_petsc.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_petsc.o : solve_petsc.F matrices_types.o matrices_tools.o globals.o types.o initialization.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
endif

ifeq ($(PSBLMG),$(LIB_YES))
solve_psblas.F : $(SDIR)/LinearAlgebra/solve_psblas.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_psblas.o : solve_psblas.F matrices_types.o matrices_tools.o globals.o types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
else ifeq ($(PSBLAS),$(LIB_YES))
solve_psblas.F : $(SDIR)/LinearAlgebra/solve_psblas.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_psblas.o : solve_psblas.F matrices_types.o matrices_tools.o globals.o types.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@
endif

read_input.F : $(SDIR)/InOut/read_input.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

read_input.o : read_input.F prec_const.o globals.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

reference_element.F : $(SDIR)/HDG/reference_element.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

reference_element.o : reference_element.F globals.o LinearAlgebra.o PrintUtils.o MPI_OMP.o HDF5_io.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

preprocess.F : $(SDIR)/HDG/preprocess.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

preprocess.o : preprocess.F globals.o PrintUtils.o types.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

flux_surface_transport_data.F : $(SDIR)/Models/$(RMDL)/transport_1d/flux_surface_transport_data.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

flux_surface_transport_data.o : flux_surface_transport_data.F globals.o MPI_OMP.o HDF5_io.o interpolation.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_config.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_config.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_config.o : transport_models_1d_config.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_derived.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_derived.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_derived.o : transport_models_1d_derived.F
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d.o : transport_models_1d.F transport_models_1d_config.o transport_models_1d_derived.o flux_surface_transport_data.o globals.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_update.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_update.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_update.o : transport_models_1d_update.F transport_models_1d.o flux_surface_transport_data.o globals.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_diffusion_bohm_gyrobohm.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_diffusion_bohm_gyrobohm.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_diffusion_bohm_gyrobohm.o : transport_models_1d_diffusion_bohm_gyrobohm.F transport_models_1d.o flux_surface_transport_data.o globals.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_common.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_common.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_common.o : transport_models_1d_common.F transport_models_1d.o flux_surface_transport_data.o globals.o physics.o interpolation.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_io.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_io.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_io.o : transport_models_1d_io.F transport_models_1d.o flux_surface_transport_data.o globals.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_runtime.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_runtime.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_runtime.o : transport_models_1d_runtime.F transport_models_1d.o flux_surface_transport_data.o globals.o physics.o interpolation.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

transport_models_1d_pinch.F : $(SDIR)/Models/$(RMDL)/transport_1d/transport_models_1d_pinch.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

transport_models_1d_pinch.o : transport_models_1d_pinch.F transport_models_1d.o flux_surface_transport_data.o globals.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

analytical.F : $(SDIR)/Models/$(RMDL)/analytical.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

analytical.o : analytical.F globals.o physics.o prec_const.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

initialization.F : $(SDIR)/HDG/initialization.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

initialization.o : initialization.F analytical.o globals.o physics.o MPI_OMP.o PrintUtils.o physics.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

hdg_PrecalculatedFirstEquation.F : $(SDIR)/HDG/hdg_PrecalculatedFirstEquation.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_PrecalculatedFirstEquation.o : hdg_PrecalculatedFirstEquation.F globals.o LinearAlgebra.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

hdg_ComputeJacobian.F : $(SDIR)/HDG/hdg_ComputeJacobian.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_ComputeJacobian.o : hdg_ComputeJacobian.F globals.o analytical.o LinearAlgebra.o physics.o hdg_LimitingTechniques.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

hdg_Mapping.F : $(SDIR)/HDG/hdg_Mapping.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_Mapping.o : hdg_Mapping.F globals.o LinearAlgebra.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

hdg_BC.F : $(SDIR)/HDG/hdg_BC.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_BC.o : hdg_BC.F globals.o analytical.o LinearAlgebra.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

hdg_Assembly.F : $(SDIR)/HDG/hdg_Assembly.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_Assembly.o : hdg_Assembly.F globals.o matrices_types.o LinearAlgebra.o PrintUtils.o Debug.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

hdg_LimitingTechniques.F : $(SDIR)/HDG/hdg_LimitingTechniques.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

hdg_LimitingTechniques.o : hdg_LimitingTechniques.F globals.o PrintUtils.o LinearAlgebra.o MPI_OMP.o reference_element.o physics.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

solve_global_system.F : $(SDIR)/HDG/solve_global_system.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

solve_global_system.o : solve_global_system.F globals.o Communications.o MPI_OMP.o LinearAlgebra.o $(ADDMODSOLV)
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

domain_decomposition.F : $(SDIR)/MPI_OMP/domain_decomposition.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

domain_decomposition.o : domain_decomposition.F globals.o MPI_OMP.o preprocess.o adaptivity_common.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

compute_element_solution.F : $(SDIR)/HDG/compute_element_solution.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

compute_element_solution.o : compute_element_solution.F globals.o LinearAlgebra.o MPI_OMP.o PrintUtils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

Postprocess.F : $(SDIR)/InOut/Postprocess.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Postprocess.o : Postprocess.F globals.o analytical.o PrintUtils.o LinearAlgebra.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

gmsh.F : $(MHDG_GMSH_DIR)/include/gmsh.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@
#-ffree-form for gfortran, -free for intel
gmsh.o : gmsh.F
	$(FC) -c -O3 -ffree-form  $(INC)  $< -o $@

adaptivity_common.F : $(SDIR)/Adaptivity/adaptivity_common.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adaptivity_common.o : adaptivity_common.F physics.o mod_splines.o globals.o gmsh.o reference_element.o LinearAlgebra.o gmsh_io.o inout.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

element_mapping.F : $(SDIR)/Adaptivity/element_mapping.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

element_mapping.o : element_mapping.F globals.o reference_element.o LinearAlgebra.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

adaptivity_projection.F : $(SDIR)/Adaptivity/adaptivity_projection.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adaptivity_projection.o : adaptivity_projection.F globals.o MPI_OMP.o element_mapping.o reference_element.o LinearAlgebra.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

adaptivity_indicator.F : $(SDIR)/Adaptivity/adaptivity_indicator.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adaptivity_indicator.o : adaptivity_indicator.F preprocess.o globals.o gmsh.o adaptivity_common.o reference_element.o inout.o physics.o Communications.o MPI_OMP.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

adaptivity_estimator.F : $(SDIR)/Adaptivity/adaptivity_estimator.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adaptivity_estimator.o : adaptivity_estimator.F globals.o gmsh.o adaptivity_common.o mod_splines.o reference_element.o MPI_OMP.o Communications.o preprocess.o gmsh_io.o physics.o LinearAlgebra.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

adaptivity_general.F : $(SDIR)/Adaptivity/adaptivity_general.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

adaptivity_general.o : adaptivity_general.F MPI_OMP.o Communications.o globals.o adaptivity_common.o adaptivity_indicator.o adaptivity_estimator.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

Main_utils.F : $(SDIR)/Utils/Main_utils.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Main_utils.o : Main_utils.F external_heating.o inout.o gmsh_io.o adaptivity_general.o adaptivity_projection.o reference_element.o preprocess.o MPI_OMP.o PrintUtils.o Debug.o initialization.o Communications.o domain_decomposition.o hdg_LimitingTechniques.o Postprocess.o $(ADDMODSOLV)
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

MHDG.F : $(SDIR)/MHDG.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

MHDG.o : MHDG.F MPI_OMP.o Main_utils.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

Convergence.F : $(SDIR)/Convergence.f90
	$(CMP) -s $< $@ || $(CP) -f $< $@

Convergence.o : Convergence.F inout.o gmsh_io.o reference_element.o MPI_OMP.o PrintUtils.o Debug.o initialization.o Postprocess.o
	$(FC) -c $(FCFLAGS) $(INC)  $< -o $@

MHDG-$(MDL)-$(MODE)-$(DIM) : Main_utils.o MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o diagnostics.o diagnostics_accumulation.o diagnostics_reduction.o diagnostics_quantities.o diagnostics_wall_sources.o diagnostics_output_terminal.o diagnostics_output_hdf5.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o external_heating.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d_config.o transport_models_1d_derived.o transport_models_1d.o transport_models_1d_update.o transport_models_1d_diffusion_bohm_gyrobohm.o transport_models_1d_common.o transport_models_1d_pinch.o transport_models_1d_io.o transport_models_1d_runtime.o inout.o gmsh_io.o gmsh.o reference_element.o preprocess.o analytical.o initialization.o element_mapping.o adaptivity_projection.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o solve_global_system.o compute_element_solution.o adaptivity_general.o adaptivity_common.o mod_splines.o adaptivity_estimator.o adaptivity_indicator.o  domain_decomposition.o MHDG.o
	$(FC) $(FCFLAGS)  Main_utils.o MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o diagnostics.o diagnostics_accumulation.o diagnostics_reduction.o diagnostics_quantities.o diagnostics_wall_sources.o diagnostics_output_terminal.o diagnostics_output_hdf5.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o external_heating.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d_config.o transport_models_1d_derived.o transport_models_1d.o transport_models_1d_update.o transport_models_1d_diffusion_bohm_gyrobohm.o transport_models_1d_common.o transport_models_1d_pinch.o transport_models_1d_io.o transport_models_1d_runtime.o inout.o gmsh_io.o gmsh.o reference_element.o preprocess.o analytical.o initialization.o element_mapping.o adaptivity_projection.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  Postprocess.o $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o  solve_global_system.o compute_element_solution.o mod_splines.o adaptivity_general.o adaptivity_common.o adaptivity_indicator.o adaptivity_estimator.o domain_decomposition.o MHDG.o -o MHDG-$(MDL)-$(MODE)-$(DIM) $(LIB)
	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test/
#	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test_anal/
#	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test_psblas/

Convergence-$(MDL)-$(DIM) : MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o diagnostics.o diagnostics_accumulation.o diagnostics_reduction.o diagnostics_quantities.o diagnostics_wall_sources.o diagnostics_output_terminal.o diagnostics_output_hdf5.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o external_heating.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d_config.o transport_models_1d_derived.o transport_models_1d.o transport_models_1d_update.o transport_models_1d_diffusion_bohm_gyrobohm.o transport_models_1d_common.o transport_models_1d_pinch.o transport_models_1d_io.o transport_models_1d_runtime.o inout.o gmsh_io.o gmsh.o reference_element.o preprocess.o analytical.o initialization.o element_mapping.o adaptivity_projection.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o solve_global_system.o compute_element_solution.o adaptivity_general.o adaptivity_common.o  mod_splines.o Convergence.o
	$(FC) $(FCFLAGS) $(OBJS) MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o diagnostics.o diagnostics_accumulation.o diagnostics_reduction.o diagnostics_quantities.o diagnostics_wall_sources.o diagnostics_output_terminal.o diagnostics_output_hdf5.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o external_heating.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d_config.o transport_models_1d_derived.o transport_models_1d.o transport_models_1d_update.o transport_models_1d_diffusion_bohm_gyrobohm.o transport_models_1d_common.o transport_models_1d_pinch.o transport_models_1d_io.o transport_models_1d_runtime.o inout.o gmsh_io.o gmsh.o reference_element.o preprocess.o analytical.o initialization.o element_mapping.o adaptivity_projection.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  Postprocess.o $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o  solve_global_system.o compute_element_solution.o adaptivity_general.o adaptivity_common.o mod_splines.o Convergence.o -o Convergence-$(MDL)-$(MODE)-$(DIM) $(LIB)
	/bin/cp -f Convergence-$(MDL)-$(MODE)-$(DIM) ../test/
#	/bin/cp -f Convergence-$(MDL)-$(MODE)-$(DIM) ../test_anal/

#MHDG-$(MDL)-$(MODE)-$(DIM) : MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o transport_coefficient_map.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o solve_global_system.o compute_element_solution.o MHDG.o
#	$(FC) $(FCFLAGS)   MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o transport_coefficient_map.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  Postprocess.o $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o  solve_global_system.o compute_element_solution.o MHDG.o -o MHDG-$(MDL)-$(MODE)-$(DIM) $(LIB)
#	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test/
##	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test_anal/
##	/bin/cp -f MHDG-$(MDL)-$(MODE)-$(DIM) ../test_psblas/
#Convergence-$(MDL)-$(DIM) : MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o transport_coefficient_map.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o solve_global_system.o compute_element_solution.o Convergence.o
#	$(FC) $(FCFLAGS) $(OBJS) MPI_OMP.o matrices_types.o interpolation.o PrintUtils.o prec_const.o types.o globals.o adimensionalization.o Debug.o LinearAlgebra.o magnetic_field.o transport_coefficient_map.o physics.o read_input.o HDF5_io.o flux_surface_transport_data.o transport_models_1d.o inout.o reference_element.o preprocess.o analytical.o initialization.o hdg_PrecalculatedFirstEquation.o hdg_ComputeJacobian.o  Postprocess.o $(ADDMOD) hdg_Mapping.o hdg_BC.o hdg_Assembly.o  solve_global_system.o compute_element_solution.o Convergence.o -o Convergence-$(MDL)-$(MODE)-$(DIM) $(LIB)
#	/bin/cp -f Convergence-$(MDL)-$(MODE)-$(DIM) ../test/
##	/bin/cp -f Convergence-$(MDL)-$(MODE)-$(DIM) ../test_anal/
