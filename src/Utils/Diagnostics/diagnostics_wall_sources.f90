SUBMODULE (diagnostics) diagnostics_wall_sources
  USE GLOBALS, ONLY: Mesh, phys
  IMPLICIT NONE

CONTAINS

  MODULE SUBROUTINE diag_prepare_wall_source_nodal_output()
#ifdef PARALL
    INTEGER :: iel, g, ind_local, ind_global, ierr
#endif

    CALL diag_reset_wall_source_nodal_output()
    IF (.NOT. diag_wall_source_nodal_available()) RETURN

#ifdef PARALL
    ALLOCATE(diag_wall_source_puff_nodal(Mesh%Nel_glob*Mesh%Nnodesperelem))
    ALLOCATE(diag_wall_source_pump_nodal(Mesh%Nel_glob*Mesh%Nnodesperelem))
    ALLOCATE(diag_wall_source_net_nodal(Mesh%Nel_glob*Mesh%Nnodesperelem))
    diag_wall_source_puff_nodal = 0.d0
    diag_wall_source_pump_nodal = 0.d0
    diag_wall_source_net_nodal = 0.d0

    DO iel = 1, Mesh%Nelems
       IF (Mesh%ghostElems(iel) .EQ. 0) THEN
          DO g = 1, Mesh%Nnodesperelem
             ind_local = (iel - 1)*Mesh%Nnodesperelem + g
             ind_global = (Mesh%loc2glob_el(iel) - 1)*Mesh%Nnodesperelem + g
             diag_wall_source_puff_nodal(ind_global) = phys%neutral_wall_source_puff_Nod(ind_local)
             diag_wall_source_pump_nodal(ind_global) = phys%neutral_wall_source_pump_Nod(ind_local)
             diag_wall_source_net_nodal(ind_global) = phys%neutral_wall_source_net_Nod(ind_local)
          ENDDO
       ENDIF
    ENDDO

    CALL MPI_Allreduce(MPI_IN_PLACE, diag_wall_source_puff_nodal, &
       &SIZE(diag_wall_source_puff_nodal), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, diag_wall_source_pump_nodal, &
       &SIZE(diag_wall_source_pump_nodal), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, diag_wall_source_net_nodal, &
       &SIZE(diag_wall_source_net_nodal), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    ALLOCATE(diag_wall_source_puff_nodal(SIZE(phys%neutral_wall_source_puff_Nod)))
    ALLOCATE(diag_wall_source_pump_nodal(SIZE(phys%neutral_wall_source_pump_Nod)))
    ALLOCATE(diag_wall_source_net_nodal(SIZE(phys%neutral_wall_source_net_Nod)))
    diag_wall_source_puff_nodal = phys%neutral_wall_source_puff_Nod
    diag_wall_source_pump_nodal = phys%neutral_wall_source_pump_Nod
    diag_wall_source_net_nodal = phys%neutral_wall_source_net_Nod
#endif
    diag_wall_source_nodal_ready = .TRUE.
  END SUBROUTINE diag_prepare_wall_source_nodal_output

  SUBROUTINE diag_reset_wall_source_nodal_output()
    diag_wall_source_nodal_ready = .FALSE.
    IF (ALLOCATED(diag_wall_source_puff_nodal)) DEALLOCATE(diag_wall_source_puff_nodal)
    IF (ALLOCATED(diag_wall_source_pump_nodal)) DEALLOCATE(diag_wall_source_pump_nodal)
    IF (ALLOCATED(diag_wall_source_net_nodal)) DEALLOCATE(diag_wall_source_net_nodal)
  END SUBROUTINE diag_reset_wall_source_nodal_output

  LOGICAL FUNCTION diag_wall_source_nodal_available() RESULT(available)
    INTEGER :: expected_size
#ifdef PARALL
    INTEGER :: available_local, available_global, ierr
#endif

    expected_size = Mesh%Nelems*Mesh%Nnodesperelem
    available = ALLOCATED(phys%neutral_wall_source_puff_Nod) .AND. &
       &ALLOCATED(phys%neutral_wall_source_pump_Nod) .AND. &
       &ALLOCATED(phys%neutral_wall_source_net_Nod)
    IF (available) THEN
       available = (SIZE(phys%neutral_wall_source_puff_Nod) .EQ. expected_size) .AND. &
          &(SIZE(phys%neutral_wall_source_pump_Nod) .EQ. expected_size) .AND. &
          &(SIZE(phys%neutral_wall_source_net_Nod) .EQ. expected_size)
    ENDIF

#ifdef PARALL
    available_local = 0
    IF (available) available_local = 1
    CALL MPI_Allreduce(available_local, available_global, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    available = available_global .EQ. 1
#endif
  END FUNCTION diag_wall_source_nodal_available

END SUBMODULE diagnostics_wall_sources
