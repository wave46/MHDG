MODULE adaptivity_general_module
    USE globals, ONLY: adapt, Mesh
    USE adaptivity_common_module, ONLY: adaptivity_console_output, calculate_h_map_elements, &
         combine_h_target_ind_est, gmsh_create_from_h_target, load_new_mesh_gmsh, save_copy_new_mesh
    USE adaptivity_indicator_module, ONLY: apply_indicator
    USE adaptivity_estimator_module, ONLY: apply_estimator

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: adaptively_refine_mesh

CONTAINS

    SUBROUTINE adaptively_refine_mesh(mesh_name, count_adapt, order)
#ifdef PARALL
        USE Communications, only: gather_mesh, gather_elemental_values
        USE mpi, ONLY: MPI_BARRIER, MPI_COMM_WORLD
        USE MPI_OMP, ONLY: MPIvar
#endif
        CHARACTER(1024), INTENT(IN) :: mesh_name
        INTEGER, INTENT(IN) :: count_adapt
        INTEGER, INTENT(IN) :: order
        REAL*8 :: h_map_elements(SIZE(Mesh%T, 1))
        REAL*8 :: h_target_elements(SIZE(Mesh%T, 1))
#ifdef PARALL
        INTEGER :: ierr
        INTEGER, POINTER :: T_global(:, :)
        REAL*8, POINTER :: X_global(:,:)
        REAL*8, POINTER :: h_target_elements_global(:)
        NULLIFY(T_global, h_target_elements_global, X_global)
#endif
        CALL adaptivity_console_output()
        CALL calculate_h_map_elements(Mesh%X, Mesh%T(:, 1:3), h_map_elements)
        CALL evaluate_adaptivity(h_map_elements, order, h_target_elements)

#ifdef PARALL
        CALL gather_mesh(Mesh,T_global,X_global)
        CALL gather_elemental_values(Mesh, h_target_elements, h_target_elements_global,allgather=.false.)
        IF (MPIvar%glob_id .EQ. 0) THEN
            CALL gmsh_create_from_h_target(h_target_elements_global,X_global, T_global(:,:3), order)
#else
            CALL gmsh_create_from_h_target(h_target_elements,Mesh%X, Mesh%T(:,:3), order)
#endif            
            CALL save_copy_new_mesh(mesh_name, count_adapt)
#ifdef PARALL
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
        CALL load_new_mesh_gmsh(order)

#ifdef PARALL
        IF (ASSOCIATED(T_global)) THEN
            DEALLOCATE(T_global)
            NULLIFY(T_global)
        ENDIF
        IF (ASSOCIATED(h_target_elements_global)) THEN
            DEALLOCATE(h_target_elements_global)
            NULLIFY(h_target_elements_global)
        ENDIF
        IF (ASSOCIATED(X_global)) THEN
            DEALLOCATE(X_global)
            NULLIFY(X_global)
        ENDIF
#endif
    ENDSUBROUTINE adaptively_refine_mesh

    SUBROUTINE evaluate_adaptivity(h_map_elements, order, h_target_elements)
        REAL*8, INTENT(IN) :: h_map_elements(:)
        INTEGER, INTENT(IN) :: order
        REAL*8              :: h_target_elements_ind(SIZE(h_map_elements))
        REAL*8              :: h_target_elements_est(SIZE(h_map_elements))
        REAL*8, INTENT(INOUT) :: h_target_elements(:)
  
        IF (adapt%evaluator .EQ. 0) THEN
            CALL apply_estimator(h_map_elements, order, h_target_elements_est)
            CALL apply_indicator(h_map_elements, h_target_elements_ind)
            CALL combine_h_target_ind_est(h_map_elements,h_target_elements_est, h_target_elements_ind, h_target_elements)
        ELSEIF (adapt%evaluator .EQ. 1) THEN
            CALL apply_indicator(h_map_elements, h_target_elements_ind)
            h_target_elements = h_target_elements_ind
        ELSEIF (adapt%evaluator .EQ. 2) THEN
            CALL apply_estimator(h_map_elements, order, h_target_elements_est)
            h_target_elements = h_target_elements_est
        ENDIF
     ENDSUBROUTINE evaluate_adaptivity

END MODULE adaptivity_general_module
