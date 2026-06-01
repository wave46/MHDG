SUBMODULE (diagnostics) diagnostics_reduction
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE diag_mpi_reduce(this)
    CLASS(diagnostics_type), INTENT(INOUT) :: this

    CALL this%mpi_reduce_boundary_hdg()
    CALL this%mpi_reduce_particles_content()
    CALL diag_reduce_array(this%energy, diag_energy_term_count)
  END SUBROUTINE diag_mpi_reduce

  MODULE SUBROUTINE diag_mpi_reduce_boundary_hdg(this)
    CLASS(diagnostics_type), INTENT(INOUT) :: this

    CALL diag_reduce_array(this%boundary_hdg, diag_boundary_term_count)
  END SUBROUTINE diag_mpi_reduce_boundary_hdg

  MODULE SUBROUTINE diag_mpi_reduce_particles_content(this)
    CLASS(diagnostics_type), INTENT(INOUT) :: this

    CALL diag_reduce_array(this%particles, diag_particle_term_count)
    CALL diag_reduce_array(this%content, diag_content_term_count)
  END SUBROUTINE diag_mpi_reduce_particles_content

  SUBROUTINE diag_reduce_array(values, nvalues)
    INTEGER, INTENT(IN) :: nvalues
    REAL*8, INTENT(INOUT) :: values(nvalues)
#ifdef PARALL
    INTEGER :: ierr

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, values, nvalues, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
  END SUBROUTINE diag_reduce_array

END SUBMODULE diagnostics_reduction
