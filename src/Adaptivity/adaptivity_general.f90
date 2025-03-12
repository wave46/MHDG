MODULE adaptivity_general_module
    USE globals
    USE adaptivity_common_module
    USE adaptivity_indicator_module
    USE adaptivity_estimator_module
    USE adaptivity_estimator_indicator_module

CONTAINS

    SUBROUTINE adaptivity_new(mesh_name,count_adapt,order)
        CHARACTER(1024), INTENT(IN)                 :: mesh_name
        INTEGER, INTENT(IN)                         :: count_adapt
        INTEGER, INTENT(IN)                         :: order
        REAL*8                                      :: h_map_elements(SIZE(Mesh%T,1))
        REAL*8                                      :: h_target_elements(SIZE(Mesh%T,1))
        REAL*8, ALLOCATABLE                             :: h_target_vertices(:)

        
        
        
        CALL adaptivity_console_output()
        
        CALL calculate_h_map_elements(Mesh%X,Mesh%T(:,1:3),h_map_elements)


        IF((adapt%evaluator .EQ. 2) .or. (adapt%evaluator .EQ. 0)) THEN
           CALL apply_estimator(h_map_elements,order,h_target_elements)
           h_map_elements = h_target_elements
        ENDIF
        
        IF((adapt%evaluator .EQ. 1) .or. (adapt%evaluator .EQ. 0)) THEN
           CALL apply_indicator(h_map_elements,h_target_elements)
        ENDIF

        CALL get_h_target_vertices(h_target_elements,h_target_vertices)

        CALL generate_new_mesh(mesh_name,h_target_vertices,count_adapt,order)

        IF (ALLOCATED(h_target_vertices)) THEN
           DEALLOCATE(h_target_vertices)        
        ENDIF

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