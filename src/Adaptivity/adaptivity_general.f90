MODULE adaptivity_general_module
    USE globals
    USE adaptivity_common_module
    USE adaptivity_indicator_module
    USE adaptivity_estimator_module
    USE adaptivity_estimator_indicator_module
    USE MPI_OMP

CONTAINS

    SUBROUTINE adaptivity_new(mesh_name,count_adapt,order)
#ifdef PARALL
        USE Communications, only: gather_connectivity, gather_elemental_values
#endif
        CHARACTER(1024), INTENT(IN)                 :: mesh_name
        INTEGER, INTENT(IN)                         :: count_adapt
        INTEGER, INTENT(IN)                         :: order
        REAL*8                                      :: h_map_elements(SIZE(Mesh%T,1))
        REAL*8                                      :: h_target_elements(SIZE(Mesh%T,1))
        REAL*8, POINTER                             :: h_target_vertices(:)
#ifdef PARALL
        INTEGER                                     :: ierr
        INTEGER,POINTER                              :: T_global(:,:)
        REAL*8,POINTER                              :: h_target_elements_global(:)
        NULLIFY(T_global,h_target_elements_global)
#endif
        NULLIFY(h_target_vertices)
        
        
        
        CALL adaptivity_console_output()
        
        CALL calculate_h_map_elements(Mesh%X,Mesh%T(:,1:3),h_map_elements)


        IF((adapt%evaluator .EQ. 2) .or. (adapt%evaluator .EQ. 0)) THEN
           CALL apply_estimator(h_map_elements,order,h_target_elements)
           h_map_elements = h_target_elements
        ENDIF
        
        IF((adapt%evaluator .EQ. 1) .or. (adapt%evaluator .EQ. 0)) THEN
           CALL apply_indicator(h_map_elements,h_target_elements)
        ENDIF

#ifdef PARALL
        CALL gather_connectivity(Mesh,T_global,allgather=.false.)
        CALL gather_elemental_values(Mesh,h_target_elements,h_target_elements_global,allgather=.false.)
        IF(MPIvar%glob_id .EQ. 0) THEN
            CALL get_h_target_vertices(h_target_elements_global,h_target_vertices,T_global)
#else
            CALL get_h_target_vertices(h_target_elements,h_target_vertices,Mesh%T)
#endif

            CALL generate_new_mesh(mesh_name,h_target_vertices,count_adapt)
#ifdef PARALL
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif
        CALL load_new_mesh(order)



        IF (ASSOCIATED(h_target_vertices)) THEN
           DEALLOCATE(h_target_vertices)
           NULLIFY(h_target_vertices)        
        ENDIF
#ifdef PARALL
        IF (ASSOCIATED(T_global)) THEN
           DEALLOCATE(T_global)
           NULLIFY(T_global)
        ENDIF
        IF (ASSOCIATED(h_target_elements_global)) THEN
           DEALLOCATE(h_target_elements_global)
           NULLIFY(h_target_elements_global)
        ENDIF
#endif

    ENDSUBROUTINE adaptivity_new

    SUBROUTINE adaptivity_console_output()
        CHARACTER(1024) :: buffer
  
        IF ((adapt%evaluator .EQ. 2)) THEN
           buffer = "            ADAPTIVITY ESTIMATOR                 "
        ELSEIF ((adapt%evaluator .EQ. 1)) THEN
           buffer = "            ADAPTIVITY INDICATOR                 "
        ELSEIF((adapt%evaluator .EQ. 0)) THEN
           buffer = "        ADAPTIVITY ESTIMATOR-INDICATOR           "
        ENDIF
  
        IF(MPIvar%glob_id .EQ. 0) THEN
           WRITE(*,*) "*************************************************"
           WRITE(*,*) TRIM(buffer)
           WRITE(*,*) "*************************************************"
        ENDIF
  
     ENDSUBROUTINE adaptivity_console_output

END MODULE adaptivity_general_module