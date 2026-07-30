!*****************************************
! project: MHDG
! file: globals.f90
! date: 06/09/2016
! Definition of the global variable structure
!*****************************************
MODULE globals
  USE prec_const
  USE types
  USE matrices_types
  IMPLICIT NONE

  !*******************************************
  ! Declaration of the GLOBAL variables
  !*******************************************

  TYPE(Reference_element_type), TARGET :: refElPol
  TYPE(Reference_element_type), TARGET :: refElTor
  TYPE(Mesh_type), TARGET :: mesh
  TYPE(Splines_DT), ALLOCATABLE, TARGET   :: splines(:)
  TYPE(Physics_type), TARGET :: phys
  TYPE(Geometry_type), TARGET :: geom
  TYPE(Magnetic_type), TARGET :: magn
  TYPE(Switches_type), TARGET :: switch
  TYPE(Inputs_type), TARGET :: input
  TYPE(Transport_model_input_type), TARGET :: transport_model_input
  TYPE(Time_type), TARGET :: time
  TYPE(Numeric_type), TARGET :: numer
  TYPE(Adaptivity_type), TARGET :: adapt
  TYPE(Utils_type), TARGET :: utils
  TYPE(Lssolver_type), TARGET :: lssolver
  TYPE(Elmat_type), TARGET :: elmat
  TYPE(Sol_type), TARGET :: sol
  TYPE(MAT_CSR_TYP), TARGET :: MatK
  TYPE(RHS_TYP), TARGET :: rhs
  TYPE(Simulationparams_type), TARGET :: simpar
  TYPE(Timing_type), TARGET :: timing

CONTAINS

  !***********************************************
  ! Routine to deallocate all previously allocated
  ! structures
  !***********************************************
  SUBROUTINE free_all()

    ! Reference element
    IF (ALLOCATED(refElPol%Face_nodes)) THEN
       DEALLOCATE (refElPol%Face_nodes)
    END IF
    IF (ALLOCATED(refElPol%Bord_nodes)) THEN
       DEALLOCATE (refElPol%Bord_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes)) THEN
       DEALLOCATE (refElPol%inner_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes_face)) THEN
       DEALLOCATE (refElPol%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElPol%coord3D)) THEN
       DEALLOCATE (refElPol%coord3D)
       NULLIFY(refElPol%coord3D)
    END IF
    IF (ASSOCIATED(refElPol%coord2D)) THEN
       DEALLOCATE (refElPol%coord2D)
       NULLIFY(refElPol%coord2D)
    END IF
    IF (ASSOCIATED(refElPol%coord1D)) THEN
       DEALLOCATE (refElPol%coord1D)
       NULLIFY(refElPol%coord1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points3D)) THEN
       DEALLOCATE (refElPol%gauss_points3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points2D)) THEN
       DEALLOCATE (refElPol%gauss_points2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points1D)) THEN
       DEALLOCATE (refElPol%gauss_points1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights3D)) THEN
       DEALLOCATE (refElPol%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights2D)) THEN
       DEALLOCATE (refElPol%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights1D)) THEN
       DEALLOCATE (refElPol%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElPol%N3D)) THEN
       DEALLOCATE (refElPol%N3D)
    END IF
    IF (ALLOCATED(refElPol%Nxi3D)) THEN
       DEALLOCATE (refElPol%Nxi3D)
    END IF
    IF (ALLOCATED(refElPol%Neta3D)) THEN
       DEALLOCATE (refElPol%Neta3D)
    END IF
    IF (ALLOCATED(refElPol%Nzeta3D)) THEN
       DEALLOCATE (refElPol%Nzeta3D)
    END IF
    IF (ALLOCATED(refElPol%N2D)) THEN
       DEALLOCATE (refElPol%N2D)
    END IF
    IF (ALLOCATED(refElPol%Nxi2D)) THEN
       DEALLOCATE (refElPol%Nxi2D)
    END IF
    IF (ALLOCATED(refElPol%Neta2D)) THEN
       DEALLOCATE (refElPol%Neta2D)
    END IF
    IF (ALLOCATED(refElPol%N1D)) THEN
       DEALLOCATE (refElPol%N1D)
    END IF
    IF (ALLOCATED(refElPol%Nxi1D)) THEN
       DEALLOCATE (refElPol%Nxi1D)
    END IF

    ! Reference element
    IF (ALLOCATED(refElTor%Face_nodes)) THEN
       DEALLOCATE (refElTor%Face_nodes)
    END IF
    IF (ALLOCATED(refElTor%Bord_nodes)) THEN
       DEALLOCATE (refElTor%Bord_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes)) THEN
       DEALLOCATE (refElTor%inner_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes_face)) THEN
       DEALLOCATE (refElTor%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElTor%coord3D)) THEN
       DEALLOCATE (refElTor%coord3D)
       NULLIFY(refElTor%coord3D)
    END IF
    IF (ASSOCIATED(refElTor%coord2D)) THEN
       DEALLOCATE (refElTor%coord2D)
       NULLIFY(refElTor%coord2D)
    END IF
    IF (ASSOCIATED(refElTor%coord1D)) THEN
       DEALLOCATE (refElTor%coord1D)
       NULLIFY(refElTor%coord1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points3D)) THEN
       DEALLOCATE (refElTor%gauss_points3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points2D)) THEN
       DEALLOCATE (refElTor%gauss_points2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points1D)) THEN
       DEALLOCATE (refElTor%gauss_points1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights3D)) THEN
       DEALLOCATE (refElTor%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights2D)) THEN
       DEALLOCATE (refElTor%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights1D)) THEN
       DEALLOCATE (refElTor%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElTor%N3D)) THEN
       DEALLOCATE (refElTor%N3D)
    END IF
    IF (ALLOCATED(refElTor%Nxi3D)) THEN
       DEALLOCATE (refElTor%Nxi3D)
    END IF
    IF (ALLOCATED(refElTor%Neta3D)) THEN
       DEALLOCATE (refElTor%Neta3D)
    END IF
    IF (ALLOCATED(refElTor%Nzeta3D)) THEN
       DEALLOCATE (refElTor%Nzeta3D)
    END IF
    IF (ALLOCATED(refElTor%N2D)) THEN
       DEALLOCATE (refElTor%N2D)
    END IF
    IF (ALLOCATED(refElTor%Nxi2D)) THEN
       DEALLOCATE (refElTor%Nxi2D)
    END IF
    IF (ALLOCATED(refElTor%Neta2D)) THEN
       DEALLOCATE (refElTor%Neta2D)
    END IF
    IF (ALLOCATED(refElTor%N1D)) THEN
       DEALLOCATE (refElTor%N1D)
    END IF
    IF (ALLOCATED(refElTor%Nxi1D)) THEN
       DEALLOCATE (refElTor%Nxi1D)
    ENDIF

    ! Mesh
    IF (ASSOCIATED(Mesh%T)) THEN
       DEALLOCATE(Mesh%T)
       NULLIFY(Mesh%T)
    END IF
    IF (ASSOCIATED(Mesh%Tlin)) THEN
       DEALLOCATE(Mesh%Tlin)
       NULLIFY(Mesh%Tlin)
    END IF
    IF (ASSOCIATED(Mesh%Tb)) THEN
       DEALLOCATE(Mesh%Tb)
       NULLIFY(Mesh%Tb)
    END IF
    IF (ASSOCIATED(Mesh%boundaryFlag)) THEN
       DEALLOCATE(Mesh%boundaryFlag)
       NULLIFY(Mesh%boundaryFlag)
    END IF
    IF (ALLOCATED(Mesh%F)) THEN
       DEALLOCATE (Mesh%F)
    END IF
    IF (ALLOCATED(Mesh%N)) THEN
       DEALLOCATE (Mesh%N)
    END IF
    IF (ALLOCATED(Mesh%face_info)) THEN
       DEALLOCATE (Mesh%face_info)
    END IF
    IF (ALLOCATED(Mesh%faces)) THEN
       DEALLOCATE (Mesh%faces)
    END IF
    IF (ALLOCATED(Mesh%extfaces)) THEN
       DEALLOCATE (Mesh%extfaces)
    END IF
    IF (ALLOCATED(Mesh%intfaces)) THEN
       DEALLOCATE (Mesh%intfaces)
    END IF
    IF (ALLOCATED(Mesh%flipFace)) THEN
       DEALLOCATE (Mesh%flipface)
    END IF
    IF (ALLOCATED(Mesh%Fdir)) THEN
       DEALLOCATE (Mesh%Fdir)
    END IF
    IF (ALLOCATED(Mesh%Diric)) THEN
       DEALLOCATE (Mesh%Diric)
    END IF
    IF (ALLOCATED(Mesh%numberbcs)) THEN
       DEALLOCATE (Mesh%numberbcs)
    END IF
    IF (ASSOCIATED(Mesh%X)) THEN
       DEALLOCATE (Mesh%X)
       NULLIFY(Mesh%X)
    END IF
    IF (ASSOCIATED(Mesh%elemSize)) THEN
       DEALLOCATE (Mesh%elemSize)
       NULLIFY(Mesh%elemSize)
    END IF
    IF (ALLOCATED(Mesh%scdiff_nodes)) THEN
       DEALLOCATE (Mesh%scdiff_nodes)
    END IF
    IF (ASSOCIATED(Mesh%toroidal)) THEN
       DEALLOCATE (Mesh%toroidal)
        NULLIFY(Mesh%toroidal)
    END IF
    IF (ALLOCATED(Mesh%periodic_faces)) THEN
       DEALLOCATE (Mesh%periodic_faces)
    END IF
    IF (ALLOCATED(Mesh%flag_elems_sc)) THEN
       DEALLOCATE (Mesh%flag_elems_sc)
    END IF

    ! sol type
    IF (ASSOCIATED(sol%u)) THEN
       DEALLOCATE (sol%u)
       NULLIFY(sol%u)
    END IF
    IF (ASSOCIATED(sol%u_conv)) THEN
       DEALLOCATE (sol%u_conv)
       NULLIFY(sol%u_conv)
    END IF
    IF (ASSOCIATED(sol%q)) THEN
       DEALLOCATE (sol%q)
       NULLIFY(sol%q)
    END IF
    IF (ASSOCIATED(sol%q_conv)) THEN
       DEALLOCATE (sol%q_conv)
       NULLIFY(sol%q_conv)
    END IF
    IF (ASSOCIATED(sol%u_tilde)) THEN
       DEALLOCATE (sol%u_tilde)
       NULLIFY(sol%u_tilde)
    END IF
    IF (ASSOCIATED(sol%u_tilde0)) THEN
       DEALLOCATE (sol%u_tilde0)
       NULLIFY(sol%u_tilde0)
    END IF
    IF (ALLOCATED(sol%u0)) THEN
       DEALLOCATE (sol%u0)
    END IF
    IF (ALLOCATED(sol%tres)) THEN
       DEALLOCATE (sol%tres)
    END IF
    IF (ALLOCATED(sol%time)) THEN
       DEALLOCATE (sol%time)
    END IF

    IF (ALLOCATED(elMat%iAqq)) THEN
       DEALLOCATE (elMat%iAqq)
    END IF
    IF (ALLOCATED(elMat%Aqu)) THEN
       DEALLOCATE (elMat%Aqu)
    END IF
    IF (ALLOCATED(elMat%Aql)) THEN
       DEALLOCATE (elMat%Aql)
    END IF
    IF (ALLOCATED(elMat%Auq)) THEN
       DEALLOCATE (elMat%Auq)
    END IF
    IF (ALLOCATED(elMat%Auu)) THEN
       DEALLOCATE (elMat%Auu)
    END IF
    IF (ALLOCATED(elMat%Aul)) THEN
       DEALLOCATE (elMat%Aul)
    END IF
    IF (ALLOCATED(elMat%Alq)) THEN
       DEALLOCATE (elMat%Alq)
    END IF
    IF (ALLOCATED(elMat%Alu)) THEN
       DEALLOCATE (elMat%Alu)
    END IF
    IF (ALLOCATED(elMat%All)) THEN
       DEALLOCATE (elMat%All)
    END IF
    IF (ALLOCATED(elMat%Aql_dir)) THEN
       DEALLOCATE (elMat%Aql_dir)
    END IF
    IF (ALLOCATED(elMat%Aul_dir)) THEN
       DEALLOCATE (elMat%Aul_dir)
    END IF
    IF (ALLOCATED(elMat%fH)) THEN
       DEALLOCATE (elMat%fH)
    END IF
    IF (ALLOCATED(elMat%LL)) THEN
       DEALLOCATE (elMat%LL)
    END IF
    IF (ALLOCATED(elMat%L0)) THEN
       DEALLOCATE (elMat%L0)
    ENDIF
    IF (ALLOCATED(elMat%S)) THEN
       DEALLOCATE (elMat%S)
    END IF
    IF (ALLOCATED(elMat%UU)) THEN
       DEALLOCATE (elMat%UU)
    END IF
    IF (ALLOCATED(elMat%U0)) THEN
       DEALLOCATE (elMat%U0)
    END IF

    ! MatK
    IF (ASSOCIATED(MatK%rowptr)) THEN
       DEALLOCATE (MatK%rowptr)
       NULLIFY(MatK%rowptr)
    END IF
    IF (ASSOCIATED(MatK%cols)) THEN
       DEALLOCATE (MatK%cols)
       NULLIFY(MatK%cols)
    END IF
    IF (ASSOCIATED(MatK%vals)) THEN
       DEALLOCATE (MatK%vals)
       NULLIFY(MatK%vals)
    END IF
    IF (ASSOCIATED(MatK%loc2glob)) THEN
       DEALLOCATE (MatK%loc2glob)
       NULLIFY(MatK%loc2glob)
    END IF

    ! RHS
    IF (ASSOCIATED(RHS%vals)) THEN
       DEALLOCATE (RHS%vals)
       NULLIFY(RHS%vals)
    END IF
    IF (ASSOCIATED(RHS%loc2glob)) THEN
       DEALLOCATE (RHS%loc2glob)
       NULLIFY(RHS%loc2glob)
    END IF

    ! Physics
    IF (ASSOCIATED(phys%B)) THEN
       DEALLOCATE (phys%B)
       NULLIFY(phys%B)
    END IF
    IF (ASSOCIATED(phys%magnetic_flux)) THEN
       DEALLOCATE (phys%magnetic_flux)
       NULLIFY(phys%magnetic_flux)
    END IF
    IF (ASSOCIATED(phys%magnetic_psi)) THEN
       DEALLOCATE (phys%magnetic_psi)
       NULLIFY(phys%magnetic_psi)
    END IF
    IF (ASSOCIATED(phys%Bperturb)) THEN
       DEALLOCATE (phys%Bperturb)
       NULLIFY(phys%Bperturb)
    END IF
    IF (ASSOCIATED(phys%phyVarNam)) THEN
       DEALLOCATE (phys%phyVarNam)
       NULLIFY(phys%phyVarNam)
    END IF
    IF (ASSOCIATED(phys%conVarNam)) THEN
       DEALLOCATE (phys%conVarNam)
       NULLIFY(phys%conVarNam)
    END IF
    IF (ALLOCATED(phys%impurity_names)) THEN
       DEALLOCATE (phys%impurity_names)
    END IF
    IF (ALLOCATED(phys%impurity_concentrations)) THEN
       DEALLOCATE (phys%impurity_concentrations)
    END IF
    IF (ALLOCATED(phys%alpha_cooling_factor_impurities)) THEN
       DEALLOCATE (phys%alpha_cooling_factor_impurities)
    END IF
    IF (ASSOCIATED(phys%Jtor)) THEN
       DEALLOCATE (phys%Jtor)
       NULLIFY(phys%Jtor)
    END IF
    IF (ALLOCATED(phys%diff_nn_Vol)) THEN
       DEALLOCATE (phys%diff_nn_Vol)
    END IF
    IF (ALLOCATED(phys%diff_nn_Fac)) THEN
       DEALLOCATE (phys%diff_nn_Fac)
    END IF
    IF (ALLOCATED(phys%diff_nn_Bou)) THEN
       DEALLOCATE (phys%diff_nn_Bou)
    END IF
    IF (ALLOCATED(phys%v_nn_Vol)) THEN
       DEALLOCATE (phys%v_nn_Vol)
    END IF
    IF (ALLOCATED(phys%v_nn_Fac)) THEN
       DEALLOCATE (phys%v_nn_Fac)
    END IF
    IF (ALLOCATED(phys%v_nn_Bou)) THEN
       DEALLOCATE (phys%v_nn_Bou)
    END IF
    IF (ASSOCIATED(phys%puff_exp)) THEN
       DEALLOCATE (phys%puff_exp)
       NULLIFY(phys%puff_exp)
    END IF

    IF (ASSOCIATED(phys%omega)) THEN
       DEALLOCATE (phys%omega)
       NULLIFY(phys%omega)
    END IF
    IF (ASSOCIATED(phys%q_cyl)) THEN
       DEALLOCATE (phys%q_cyl)
       NULLIFY(phys%q_cyl)
    END IF


    IF (ASSOCIATED(phys%rho_1D)) THEN
       DEALLOCATE (phys%rho_1D)
       NULLIFY(phys%rho_1D)
    END IF

    IF (ASSOCIATED(phys%diff_n_1D)) THEN
       DEALLOCATE (phys%diff_n_1D)
       NULLIFY(phys%diff_n_1D)
    END IF

    IF (ASSOCIATED(phys%diff_u_1D)) THEN
       DEALLOCATE (phys%diff_u_1D)
       NULLIFY(phys%diff_u_1D)
    END IF

    IF (ASSOCIATED(phys%diff_e_1D)) THEN
       DEALLOCATE (phys%diff_e_1D)
       NULLIFY(phys%diff_e_1D)
    END IF

    IF (ASSOCIATED(phys%diff_ee_1D)) THEN
       DEALLOCATE (phys%diff_ee_1D)
       NULLIFY(phys%diff_ee_1D)
   END IF   

    ! magnetic
    IF (ASSOCIATED(magn%coils_rmp)) THEN
       DEALLOCATE (magn%coils_rmp)
       NULLIFY(magn%coils_rmp)
    END IF
    IF (ASSOCIATED(magn%coils_ripple)) THEN
       DEALLOCATE (magn%coils_ripple)
       NULLIFY(magn%coils_ripple)
    END IF

    IF (ALLOCATED(simpar%physvar_refval)) THEN
       DEALLOCATE (simpar%physvar_refval)
    END IF
    IF (ALLOCATED(simpar%consvar_refval)) THEN
       DEALLOCATE (simpar%consvar_refval)
    END IF

    IF (ASSOCIATED(phys%external_heating_ions)) THEN
      DEALLOCATE (phys%external_heating_ions,phys%external_heating_electrons)
      NULLIFY (phys%external_heating_ions,phys%external_heating_electrons)
    END IF
         

  END SUBROUTINE free_all

  SUBROUTINE free_mesh
    ! Mesh
    IF (ASSOCIATED(Mesh%T)) THEN
       DEALLOCATE (Mesh%T)
       NULLIFY(Mesh%T)
    END IF
    IF (ASSOCIATED(Mesh%Tlin)) THEN
       DEALLOCATE (Mesh%Tlin)
       NULLIFY(Mesh%Tlin)
    END IF
    IF (ASSOCIATED(Mesh%Tb)) THEN
       DEALLOCATE (Mesh%Tb)
       NULLIFY(Mesh%Tb)
    END IF
    IF (ASSOCIATED(Mesh%boundaryFlag)) THEN
       DEALLOCATE (Mesh%boundaryFlag)
       NULLIFY(Mesh%boundaryFlag)
    END IF
    IF (ALLOCATED(Mesh%F)) THEN
       DEALLOCATE (Mesh%F)
    END IF
    IF (ALLOCATED(Mesh%face_info)) THEN
       DEALLOCATE (Mesh%face_info)
    END IF
    IF (ALLOCATED(Mesh%N)) THEN
       DEALLOCATE (Mesh%N)
    END IF
    IF (ALLOCATED(Mesh%faces)) THEN
       DEALLOCATE (Mesh%faces)
    END IF
    IF (ALLOCATED(Mesh%extfaces)) THEN
       DEALLOCATE (Mesh%extfaces)
    END IF
    IF (ALLOCATED(Mesh%intfaces)) THEN
       DEALLOCATE (Mesh%intfaces)
    END IF
    IF (ALLOCATED(Mesh%flipFace)) THEN
       DEALLOCATE (Mesh%flipface)
    END IF
    IF (ALLOCATED(Mesh%Fdir)) THEN
       DEALLOCATE (Mesh%Fdir)
    END IF
    IF (ALLOCATED(Mesh%Diric)) THEN
       DEALLOCATE (Mesh%Diric)
    END IF
    IF (ALLOCATED(Mesh%numberbcs)) THEN
       DEALLOCATE (Mesh%numberbcs)
    END IF
    IF (ASSOCIATED(Mesh%X)) THEN
       DEALLOCATE (Mesh%X)
       NULLIFY(Mesh%X)
    END IF
    IF (ALLOCATED(Mesh%flag_elems_sc)) THEN
       DEALLOCATE (Mesh%flag_elems_sc)
    END IF
    IF (ASSOCIATED(Mesh%elemSize)) THEN
       DEALLOCATE (Mesh%elemSize)
       NULLIFY(Mesh%elemSize)
    END IF
    IF (ALLOCATED(Mesh%scdiff_nodes)) THEN
       DEALLOCATE (Mesh%scdiff_nodes)
    END IF
    IF (ASSOCIATED(Mesh%toroidal)) THEN
       DEALLOCATE (Mesh%toroidal)
       NULLIFY(Mesh%toroidal)
    END IF
    IF (ALLOCATED(Mesh%periodic_faces)) THEN
       DEALLOCATE (Mesh%periodic_faces)
    END IF

#ifdef PARALL
    IF (ASSOCIATED(Mesh%loc2glob_fa)) THEN
       DEALLOCATE(Mesh%loc2glob_fa)
       NULLIFY(Mesh%loc2glob_fa)
    END IF
    IF (ASSOCIATED(Mesh%loc2glob_el)) THEN
       DEALLOCATE(Mesh%loc2glob_el)
       NULLIFY(Mesh%loc2glob_el)
    END IF
    IF (ASSOCIATED(Mesh%loc2glob_nodes)) THEN
       DEALLOCATE(Mesh%loc2glob_nodes)
       NULLIFY(Mesh%loc2glob_nodes)
    END IF
    IF (ASSOCIATED(Mesh%ghostfaces)) THEN
       DEALLOCATE(Mesh%ghostfaces)
       NULLIFY(Mesh%ghostfaces)
    END IF
    IF (ASSOCIATED(Mesh%ghostelems)) THEN
       DEALLOCATE(Mesh%ghostelems)
       NULLIFY(Mesh%ghostelems)
    END IF

    IF (ASSOCIATED(Mesh%ghostflp)) THEN
       DEALLOCATE(Mesh%ghostflp)
       NULLIFY(Mesh%ghostflp)
    END IF
    IF (ASSOCIATED(Mesh%ghostloc)) THEN
       DEALLOCATE(Mesh%ghostloc)
       NULLIFY(Mesh%ghostloc)
    END IF
    IF (ASSOCIATED(Mesh%ghostpro)) THEN
       DEALLOCATE(Mesh%ghostpro)
       NULLIFY(Mesh%ghostpro)
    END IF
    IF (ASSOCIATED(Mesh%ghelsloc)) THEN
       DEALLOCATE(Mesh%ghelsloc)
       NULLIFY(Mesh%ghelsloc)
    END IF
    IF (ASSOCIATED(Mesh%ghelspro)) THEN
       DEALLOCATE(Mesh%ghelspro)
       NULLIFY(Mesh%ghelspro)
    END IF
    IF (ALLOCATED(Mesh%fc2sd)) THEN
       DEALLOCATE (Mesh%fc2sd)
    END IF
    IF (ALLOCATED(Mesh%pr2sd)) THEN
       DEALLOCATE (Mesh%pr2sd)
    END IF
    IF (ALLOCATED(Mesh%fc2rv)) THEN
       DEALLOCATE (Mesh%fc2rv)
    END IF
    IF (ALLOCATED(Mesh%pr2rv)) THEN
       DEALLOCATE (Mesh%pr2rv)
    END IF
    IF (ALLOCATED(Mesh%el2sd)) THEN
       DEALLOCATE (Mesh%el2sd)
    END IF
    IF (ALLOCATED(Mesh%pe2sd)) THEN
       DEALLOCATE (Mesh%pe2sd)
    END IF
    IF (ALLOCATED(Mesh%el2rv)) THEN
       DEALLOCATE (Mesh%el2rv)
    END IF
    IF (ALLOCATED(Mesh%pe2rv)) THEN
       DEALLOCATE (Mesh%pe2rv)
    END IF
#endif

  END SUBROUTINE free_mesh

  SUBROUTINE free_mesh_loc(Mesh_loc)
    TYPE(Mesh_type), INTENT(INOUT) :: Mesh_loc
    ! Mesh
    IF (ASSOCIATED(Mesh_loc%T)) THEN
       DEALLOCATE(Mesh_loc%T)
       NULLIFY(Mesh_loc%T)
    END IF
    IF (ASSOCIATED(Mesh_loc%Tlin)) THEN
       DEALLOCATE(Mesh_loc%Tlin)
       NULLIFY(Mesh_loc%Tlin)
    END IF
    IF (ASSOCIATED(Mesh_loc%Tb)) THEN
       DEALLOCATE(Mesh_loc%Tb)
       NULLIFY(Mesh_loc%Tb)
    END IF
    IF (ASSOCIATED(Mesh_loc%boundaryFlag)) THEN
       DEALLOCATE(Mesh_loc%boundaryFlag)
       NULLIFY(Mesh_loc%boundaryFlag)
    END IF
    IF (ALLOCATED(Mesh_loc%F)) THEN
       DEALLOCATE (Mesh_loc%F)
    END IF
    IF (ALLOCATED(Mesh_loc%face_info)) THEN
       DEALLOCATE (Mesh_loc%face_info)
    END IF
    IF (ALLOCATED(Mesh_loc%N)) THEN
       DEALLOCATE (Mesh_loc%N)
    END IF
    IF (ALLOCATED(Mesh_loc%faces)) THEN
       DEALLOCATE (Mesh_loc%faces)
    END IF
    IF (ALLOCATED(Mesh_loc%extfaces)) THEN
       DEALLOCATE (Mesh_loc%extfaces)
    END IF
    IF (ALLOCATED(Mesh_loc%intfaces)) THEN
       DEALLOCATE (Mesh_loc%intfaces)
    END IF
    IF (ALLOCATED(Mesh_loc%flipFace)) THEN
       DEALLOCATE (Mesh_loc%flipface)
    END IF
    IF (ALLOCATED(Mesh_loc%Fdir)) THEN
       DEALLOCATE (Mesh_loc%Fdir)
    END IF
    IF (ALLOCATED(Mesh_loc%Diric)) THEN
       DEALLOCATE (Mesh_loc%Diric)
    END IF
    IF (ALLOCATED(Mesh_loc%numberbcs)) THEN
       DEALLOCATE (Mesh_loc%numberbcs)
    END IF
    IF (ASSOCIATED(Mesh_loc%X)) THEN
       DEALLOCATE (Mesh_loc%X)
       NULLIFY(Mesh_loc%X)
    END IF
    IF (ASSOCIATED(Mesh_loc%elemSize)) THEN
       DEALLOCATE (Mesh_loc%elemSize)
       NULLIFY(Mesh_loc%elemSize)
    END IF
    IF (ALLOCATED(Mesh_loc%scdiff_nodes)) THEN
       DEALLOCATE (Mesh_loc%scdiff_nodes)
    END IF
    IF (ASSOCIATED(Mesh_loc%toroidal)) THEN
       DEALLOCATE (Mesh_loc%toroidal)
       NULLIFY(Mesh_loc%toroidal)
    END IF
    IF (ALLOCATED(Mesh_loc%periodic_faces)) THEN
       DEALLOCATE (Mesh_loc%periodic_faces)
    END IF

#ifdef PARALL
    IF (ASSOCIATED(Mesh_loc%loc2glob_fa)) THEN
       DEALLOCATE(Mesh_loc%loc2glob_fa)
       NULLIFY(Mesh_loc%loc2glob_fa)
    END IF
    IF (ASSOCIATED(Mesh_loc%loc2glob_el)) THEN
       DEALLOCATE(Mesh_loc%loc2glob_el)
       NULLIFY(Mesh_loc%loc2glob_el)
    END IF
    IF (ASSOCIATED(Mesh_loc%loc2glob_nodes)) THEN
       DEALLOCATE(Mesh_loc%loc2glob_nodes)
       NULLIFY(Mesh_loc%loc2glob_nodes)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghostfaces)) THEN
       DEALLOCATE(Mesh_loc%ghostfaces)
       NULLIFY(Mesh_loc%ghostfaces)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghostelems)) THEN
       DEALLOCATE(Mesh_loc%ghostelems)
       NULLIFY(Mesh_loc%ghostelems)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghostflp)) THEN
       DEALLOCATE(Mesh_loc%ghostflp)
       NULLIFY(Mesh_loc%ghostflp)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghostloc)) THEN
       DEALLOCATE(Mesh_loc%ghostloc)
       NULLIFY(Mesh_loc%ghostloc)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghostpro)) THEN
       DEALLOCATE(Mesh_loc%ghostpro)
       NULLIFY(Mesh_loc%ghostpro)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghelsloc)) THEN
       DEALLOCATE(Mesh_loc%ghelsloc)
       NULLIFY(Mesh_loc%ghelsloc)
    END IF
    IF (ASSOCIATED(Mesh_loc%ghelspro)) THEN
       DEALLOCATE(Mesh_loc%ghelspro)
       NULLIFY(Mesh_loc%ghelspro)
    END IF
    IF (ALLOCATED(Mesh_loc%fc2sd)) THEN
       DEALLOCATE (Mesh_loc%fc2sd)
    END IF
    IF (ALLOCATED(Mesh_loc%pr2sd)) THEN
       DEALLOCATE (Mesh_loc%pr2sd)
    END IF
    IF (ALLOCATED(Mesh_loc%fc2rv)) THEN
       DEALLOCATE (Mesh_loc%fc2rv)
    END IF
    IF (ALLOCATED(Mesh_loc%pr2rv)) THEN
       DEALLOCATE (Mesh_loc%pr2rv)
    END IF
    IF (ALLOCATED(Mesh_loc%el2sd)) THEN
       DEALLOCATE (Mesh_loc%el2sd)
    END IF
    IF (ALLOCATED(Mesh_loc%pe2sd)) THEN
       DEALLOCATE (Mesh_loc%pe2sd)
    END IF
    IF (ALLOCATED(Mesh_loc%el2rv)) THEN
       DEALLOCATE (Mesh_loc%el2rv)
    END IF
    IF (ALLOCATED(Mesh_loc%pe2rv)) THEN
       DEALLOCATE (Mesh_loc%pe2rv)
    END IF
#endif

END SUBROUTINE free_mesh_loc


  SUBROUTINE free_reference_element
    ! Reference element
    IF (ALLOCATED(refElPol%Face_nodes)) THEN
       DEALLOCATE (refElPol%Face_nodes)
    END IF
    IF (ALLOCATED(refElPol%Bord_nodes)) THEN
       DEALLOCATE (refElPol%Bord_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes)) THEN
       DEALLOCATE (refElPol%inner_nodes)
    END IF
    IF (ALLOCATED(refElPol%inner_nodes_face)) THEN
       DEALLOCATE (refElPol%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElPol%coord3D)) THEN
       DEALLOCATE (refElPol%coord3D)
       NULLIFY(refElPol%coord3D)
    END IF
    IF (ASSOCIATED(refElPol%coord2D)) THEN
       DEALLOCATE (refElPol%coord2D)
       NULLIFY(refElPol%coord2D)
    END IF
    IF (ASSOCIATED(refElPol%coord1D)) THEN
       DEALLOCATE (refElPol%coord1D)
       NULLIFY(refElPol%coord1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points3D)) THEN
       DEALLOCATE (refElPol%gauss_points3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points2D)) THEN
       DEALLOCATE (refElPol%gauss_points2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_points1D)) THEN
       DEALLOCATE (refElPol%gauss_points1D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights3D)) THEN
       DEALLOCATE (refElPol%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights2D)) THEN
       DEALLOCATE (refElPol%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElPol%gauss_weights1D)) THEN
       DEALLOCATE (refElPol%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElPol%N3D)) THEN
       DEALLOCATE (refElPol%N3D)
    END IF
    IF (ALLOCATED(refElPol%Nxi3D)) THEN
       DEALLOCATE (refElPol%Nxi3D)
    END IF
    IF (ALLOCATED(refElPol%Neta3D)) THEN
       DEALLOCATE (refElPol%Neta3D)
    END IF
    IF (ALLOCATED(refElPol%Nzeta3D)) THEN
       DEALLOCATE (refElPol%Nzeta3D)
    END IF
    IF (ALLOCATED(refElPol%N2D)) THEN
       DEALLOCATE (refElPol%N2D)
    END IF
    IF (ALLOCATED(refElPol%Nxi2D)) THEN
       DEALLOCATE (refElPol%Nxi2D)
    END IF
    IF (ALLOCATED(refElPol%Neta2D)) THEN
       DEALLOCATE (refElPol%Neta2D)
    END IF
    IF (ALLOCATED(refElPol%N1D)) THEN
       DEALLOCATE (refElPol%N1D)
    END IF
    IF (ALLOCATED(refElPol%Nxi1D)) THEN
       DEALLOCATE (refElPol%Nxi1D)
    END IF
    IF (ALLOCATED(refElPol%Nlin)) THEN
       DEALLOCATE (refElPol%Nlin)
    END IF

    ! Reference element
    IF (ALLOCATED(refElTor%Face_nodes)) THEN
       DEALLOCATE (refElTor%Face_nodes)
    END IF
    IF (ALLOCATED(refElTor%Bord_nodes)) THEN
       DEALLOCATE (refElTor%Bord_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes)) THEN
       DEALLOCATE (refElTor%inner_nodes)
    END IF
    IF (ALLOCATED(refElTor%inner_nodes_face)) THEN
       DEALLOCATE (refElTor%inner_nodes_face)
    END IF
    IF (ASSOCIATED(refElTor%coord3D)) THEN
       DEALLOCATE (refElTor%coord3D)
       NULLIFY(refElTor%coord3D)
    END IF
    IF (ASSOCIATED(refElTor%coord2D)) THEN
       DEALLOCATE (refElTor%coord2D)
       NULLIFY(refElTor%coord2D)
    END IF
    IF (ASSOCIATED(refElTor%coord1D)) THEN
       DEALLOCATE (refElTor%coord1D)
       NULLIFY(refElTor%coord1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points3D)) THEN
       DEALLOCATE (refElTor%gauss_points3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points2D)) THEN
       DEALLOCATE (refElTor%gauss_points2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_points1D)) THEN
       DEALLOCATE (refElTor%gauss_points1D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights3D)) THEN
       DEALLOCATE (refElTor%gauss_weights3D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights2D)) THEN
       DEALLOCATE (refElTor%gauss_weights2D)
    END IF
    IF (ALLOCATED(refElTor%gauss_weights1D)) THEN
       DEALLOCATE (refElTor%gauss_weights1D)
    END IF
    IF (ALLOCATED(refElTor%N3D)) THEN
       DEALLOCATE (refElTor%N3D)
    END IF
    IF (ALLOCATED(refElTor%Nxi3D)) THEN
       DEALLOCATE (refElTor%Nxi3D)
    END IF
    IF (ALLOCATED(refElTor%Neta3D)) THEN
       DEALLOCATE (refElTor%Neta3D)
    END IF
    IF (ALLOCATED(refElTor%Nzeta3D)) THEN
       DEALLOCATE (refElTor%Nzeta3D)
    END IF
    IF (ALLOCATED(refElTor%N2D)) THEN
       DEALLOCATE (refElTor%N2D)
    END IF
    IF (ALLOCATED(refElTor%Nxi2D)) THEN
       DEALLOCATE (refElTor%Nxi2D)
    END IF
    IF (ALLOCATED(refElTor%Neta2D)) THEN
       DEALLOCATE (refElTor%Neta2D)
    END IF
    IF (ALLOCATED(refElTor%N1D)) THEN
       DEALLOCATE (refElTor%N1D)
    END IF
    IF (ALLOCATED(refElTor%Nxi1D)) THEN
       DEALLOCATE (refElTor%Nxi1D)
    ENDIF

  END SUBROUTINE free_reference_element

  SUBROUTINE free_reference_element_pol(RefEl)
    TYPE(Reference_element_type), INTENT(INOUT) :: refEl

    ! Reference element
    IF (ALLOCATED(RefEl%Face_nodes)) THEN
       DEALLOCATE (RefEl%Face_nodes)
    END IF
    IF (ALLOCATED(RefEl%Bord_nodes)) THEN
       DEALLOCATE (RefEl%Bord_nodes)
    END IF
    IF (ALLOCATED(RefEl%inner_nodes)) THEN
       DEALLOCATE (RefEl%inner_nodes)
    END IF
    IF (ALLOCATED(RefEl%inner_nodes_face)) THEN
       DEALLOCATE (RefEl%inner_nodes_face)
    END IF
    IF (ASSOCIATED(RefEl%coord3D)) THEN
       DEALLOCATE (RefEl%coord3D)
       NULLIFY(RefEl%coord3D)
    END IF
    IF (ASSOCIATED(RefEl%coord2D)) THEN
       DEALLOCATE (RefEl%coord2D)
       NULLIFY(RefEl%coord2D)
    END IF
    IF (ASSOCIATED(RefEl%coord1D)) THEN
       DEALLOCATE (RefEl%coord1D)
       NULLIFY(RefEl%coord1D)
    END IF
    IF (ALLOCATED(RefEl%gauss_points3D)) THEN
       DEALLOCATE (RefEl%gauss_points3D)
    END IF
    IF (ALLOCATED(RefEl%gauss_points2D)) THEN
       DEALLOCATE (RefEl%gauss_points2D)
    END IF
    IF (ALLOCATED(RefEl%gauss_points1D)) THEN
       DEALLOCATE (RefEl%gauss_points1D)
    END IF
    IF (ALLOCATED(RefEl%gauss_weights3D)) THEN
       DEALLOCATE (RefEl%gauss_weights3D)
    END IF
    IF (ALLOCATED(RefEl%gauss_weights2D)) THEN
       DEALLOCATE (RefEl%gauss_weights2D)
    END IF
    IF (ALLOCATED(RefEl%gauss_weights1D)) THEN
       DEALLOCATE (RefEl%gauss_weights1D)
    END IF
    IF (ALLOCATED(RefEl%N3D)) THEN
       DEALLOCATE (RefEl%N3D)
    END IF
    IF (ALLOCATED(RefEl%Nxi3D)) THEN
       DEALLOCATE (RefEl%Nxi3D)
    END IF
    IF (ALLOCATED(RefEl%Neta3D)) THEN
       DEALLOCATE (RefEl%Neta3D)
    END IF
    IF (ALLOCATED(RefEl%Nzeta3D)) THEN
       DEALLOCATE (RefEl%Nzeta3D)
    END IF
    IF (ALLOCATED(RefEl%N2D)) THEN
       DEALLOCATE (RefEl%N2D)
    END IF
    IF (ALLOCATED(RefEl%Nxi2D)) THEN
       DEALLOCATE (RefEl%Nxi2D)
    END IF
    IF (ALLOCATED(RefEl%Neta2D)) THEN
       DEALLOCATE (RefEl%Neta2D)
    END IF
    IF (ALLOCATED(RefEl%N1D)) THEN
       DEALLOCATE (RefEl%N1D)
    END IF
    IF (ALLOCATED(RefEl%Nxi1D)) THEN
       DEALLOCATE (RefEl%Nxi1D)
    END IF
    IF (ALLOCATED(RefEl%Nlin)) THEN
       DEALLOCATE (RefEl%Nlin)
    END IF
    IF (ALLOCATED(RefEl%sFTF)) THEN
       DEALLOCATE (RefEl%sFTF)
    END IF
    IF (ALLOCATED(RefEl%faceNodes3)) THEN
       DEALLOCATE (RefEl%faceNodes3)
    END IF

  END SUBROUTINE free_reference_element_pol

  SUBROUTINE free_reference_element_tor(RefEl)
    TYPE(Reference_element_type), INTENT(INOUT) :: refEl
    ! Reference element
    IF (ALLOCATED(RefEl%Face_nodes)) THEN
       DEALLOCATE (RefEl%Face_nodes)
    END IF
    IF (ALLOCATED(RefEl%Bord_nodes)) THEN
       DEALLOCATE (RefEl%Bord_nodes)
    END IF
    IF (ALLOCATED(RefEl%inner_nodes)) THEN
       DEALLOCATE (RefEl%inner_nodes)
    END IF
    IF (ALLOCATED(RefEl%inner_nodes_face)) THEN
       DEALLOCATE (RefEl%inner_nodes_face)
    END IF
    IF (ASSOCIATED(RefEl%coord3D)) THEN
       DEALLOCATE (RefEl%coord3D)
       NULLIFY(RefEl%coord3D)
    END IF
    IF (ASSOCIATED(RefEl%coord2D)) THEN
       DEALLOCATE (RefEl%coord2D)
       NULLIFY(RefEl%coord2D)
    END IF
    IF (ASSOCIATED(RefEl%coord1D)) THEN
       DEALLOCATE (RefEl%coord1D)
       NULLIFY(RefEl%coord1D)
    END IF
    IF (ALLOCATED(RefEl%gauss_points3D)) THEN
       DEALLOCATE (RefEl%gauss_points3D)
    END IF
    IF (ALLOCATED(RefEl%gauss_points2D)) THEN
       DEALLOCATE (RefEl%gauss_points2D)
    END IF
    IF (ALLOCATED(RefEl%gauss_points1D)) THEN
       DEALLOCATE (RefEl%gauss_points1D)
    END IF
    IF (ALLOCATED(RefEl%gauss_weights3D)) THEN
       DEALLOCATE (RefEl%gauss_weights3D)
    END IF
    IF (ALLOCATED(RefEl%gauss_weights2D)) THEN
       DEALLOCATE (RefEl%gauss_weights2D)
    END IF
    IF (ALLOCATED(RefEl%gauss_weights1D)) THEN
       DEALLOCATE (RefEl%gauss_weights1D)
    END IF
    IF (ALLOCATED(RefEl%N3D)) THEN
       DEALLOCATE (RefEl%N3D)
    END IF
    IF (ALLOCATED(RefEl%Nxi3D)) THEN
       DEALLOCATE (RefEl%Nxi3D)
    END IF
    IF (ALLOCATED(RefEl%Neta3D)) THEN
       DEALLOCATE (RefEl%Neta3D)
    END IF
    IF (ALLOCATED(RefEl%Nzeta3D)) THEN
       DEALLOCATE (RefEl%Nzeta3D)
    END IF
    IF (ALLOCATED(RefEl%N2D)) THEN
       DEALLOCATE (RefEl%N2D)
    END IF
    IF (ALLOCATED(RefEl%Nxi2D)) THEN
       DEALLOCATE (RefEl%Nxi2D)
    END IF
    IF (ALLOCATED(RefEl%Neta2D)) THEN
       DEALLOCATE (RefEl%Neta2D)
    END IF
    IF (ALLOCATED(RefEl%N1D)) THEN
       DEALLOCATE (RefEl%N1D)
    END IF
    IF (ALLOCATED(RefEl%Nxi1D)) THEN
       DEALLOCATE (RefEl%Nxi1D)
    ENDIF

  END SUBROUTINE free_reference_element_tor

  SUBROUTINE free_el_mat
    IF (ALLOCATED(elMat%iAqq)) THEN
       DEALLOCATE (elMat%iAqq)
    END IF
    IF (ALLOCATED(elMat%Aqu)) THEN
       DEALLOCATE (elMat%Aqu)
    END IF
    IF (ALLOCATED(elMat%Aql)) THEN
       DEALLOCATE (elMat%Aql)
    END IF
    IF (ALLOCATED(elMat%Auq)) THEN
       DEALLOCATE (elMat%Auq)
    END IF
    IF (ALLOCATED(elMat%Auu)) THEN
       DEALLOCATE (elMat%Auu)
    END IF
    IF (ALLOCATED(elMat%Aul)) THEN
       DEALLOCATE (elMat%Aul)
    END IF
    IF (ALLOCATED(elMat%Alq)) THEN
       DEALLOCATE (elMat%Alq)
    END IF
    IF (ALLOCATED(elMat%Alu)) THEN
       DEALLOCATE (elMat%Alu)
    END IF
    IF (ALLOCATED(elMat%All)) THEN
       DEALLOCATE (elMat%All)
    END IF
    IF (ALLOCATED(elMat%Aql_dir)) THEN
       DEALLOCATE (elMat%Aql_dir)
    END IF
    IF (ALLOCATED(elMat%Aul_dir)) THEN
       DEALLOCATE (elMat%Aul_dir)
    END IF
    IF (ALLOCATED(elMat%fH)) THEN
       DEALLOCATE (elMat%fH)
    END IF
    IF (ALLOCATED(elMat%LL)) THEN
       DEALLOCATE (elMat%LL)
    END IF
    IF (ALLOCATED(elMat%L0)) THEN
       DEALLOCATE (elMat%L0)
    ENDIF
    IF (ALLOCATED(elMat%S)) THEN
       DEALLOCATE (elMat%S)
    END IF
    IF (ALLOCATED(elMat%UU)) THEN
       DEALLOCATE (elMat%UU)
    END IF
    IF (ALLOCATED(elMat%U0)) THEN
       DEALLOCATE (elMat%U0)
    END IF

  ENDSUBROUTINE free_el_mat

  SUBROUTINE free_mat
    IF(ASSOCIATED(MatK%cols)) THEN
      DEALLOCATE (MatK%cols)
      NULLIFY(MatK%cols)
    ENDIF
    IF(ASSOCIATED(MatK%rowptr)) THEN
      DEALLOCATE (MatK%rowptr)
      NULLIFY(MatK%rowptr)
    ENDIF
    IF(ASSOCIATED(MatK%vals)) THEN
      DEALLOCATE (MatK%vals)
      NULLIFY(MatK%vals)
    ENDIF
    IF(ASSOCIATED(MatK%loc2glob)) THEN
      DEALLOCATE (MatK%loc2glob)
      NULLIFY(MatK%loc2glob)
    ENDIF
    IF(ASSOCIATED(rhs%loc2glob)) THEN
      DEALLOCATE (rhs%loc2glob)
      NULLIFY(rhs%loc2glob)
    ENDIF
    IF(ASSOCIATED(rhs%vals)) THEN
      DEALLOCATE (rhs%vals)
      NULLIFY(rhs%vals)
    ENDIF

  END SUBROUTINE free_mat


  SUBROUTINE free_splines(splines_str_array)
    TYPE(Splines_DT), ALLOCATABLE, INTENT(INOUT)  :: splines_str_array(:)
    INTEGER                                       :: i
    IF(ALLOCATED(splines_str_array)) THEN
       DO i = 1, SIZE(splines_str_array)
          IF(ALLOCATED(splines_str_array(i)%points_number)) DEALLOCATE (splines_str_array(i)%points_number)
          IF(ALLOCATED(splines_str_array(i)%points_coord)) DEALLOCATE (splines_str_array(i)%points_coord)
       ENDDO
       DEALLOCATE(splines_str_array)
    ENDIF

  END SUBROUTINE free_splines

END MODULE globals
