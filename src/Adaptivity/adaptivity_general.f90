MODULE adaptivity_general_module
    USE globals
    USE adaptivity_common_module
    USE adaptivity_indicator_module
    USE adaptivity_estimator_module
    USE adaptivity_estimator_indicator_module

CONTAINS

    SUBROUTINE adaptivity_new(mesh_name,count_adapt,restart_adapt)

        ! mesh_name is the mesh path + mesh name + .msh extension ("./Meshes/CircLim.msh")
        ! mesh_name_npne (mesh name no path no extension) is just the name of the mesh ("CircLim")
        ! new_mesh_name_npne (new mesh name no path no extension) is just the name of the mesh + param_adapt + count_adapt ("CircLim_param2_n1")
        ! buffer is a dummy array to store intermediate mesh names
        CHARACTER(1024), INTENT(IN)                 :: mesh_name
        INTEGER, INTENT(IN)                         :: count_adapt
        LOGICAL, INTENT(IN)                         :: restart_adapt
        REAL*8                                      :: h_map_elements(SIZE(Mesh%T,1))
        REAL*8                                      :: oscillation_element(SIZE(Mesh%T,1))
        
        
        
        CALL adaptivity_console_output(restart_adapt)
        
        CALL calculate_h_map_elements(Mesh%X,Mesh%T(:,1:3),h_map_elements)
        
        CALL calculate_oscillations(oscillation_element)
 
    ENDSUBROUTINE adaptivity_new

    SUBROUTINE adaptivity_console_output(restart_adapt)
        LOGICAL, INTENT(IN) :: restart_adapt
        CHARACTER(1024) :: buffer
  
        IF ((adapt%evaluator .EQ. 2)) THEN
           buffer = "            ADAPTIVITY ESTIMATOR                 "
        ELSEIF ((adapt%evaluator .EQ. 1)) THEN
           buffer = "            ADAPTIVITY INDICATOR                 "
        ELSEIF((adapt%evaluator .EQ. 0) .OR. (restart_adapt)) THEN
           buffer = "        ADAPTIVITY ESTIMATOR-INDICATOR           "
        ENDIF
  
        IF(MPIvar%glob_id .EQ. 0) THEN
           WRITE(*,*) "*************************************************"
           WRITE(*,*) TRIM(buffer)
           WRITE(*,*) "*************************************************"
        ENDIF
  
     ENDSUBROUTINE adaptivity_console_output

END MODULE adaptivity_general_module