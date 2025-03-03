!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_indicator_module
  USE globals
  USE reference_element
  USE gmsh
  USE adaptivity_common_module
  USE MPI_OMP
  IMPLICIT NONE

CONTAINS

  SUBROUTINE adaptivity_indicator(mesh_name,thresh, param_adapt, count_adapt, order)
    USE in_out, ONLY: copy_file
    USE gmsh_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, convert_gmsh_to_hdf5
    USE preprocess
#ifdef PARALL
    USE Communications, ONLY: gather_1D_vector_int,gather_1D_vector_real,gather_mesh
#endif
    TYPE(gmsh_t)                                :: gmsh
    REAL*8, INTENT(IN)                          :: thresh
    INTEGER, INTENT(IN)                         :: param_adapt, count_adapt, order
    REAL*8,  ALLOCATABLE                        :: h(:), h_target(:), error_oscillation(:)
    INTEGER, ALLOCATABLE                        :: vector_nodes_unique(:)
    INTEGER                                     :: i, N_n_vertex, n_el_unstable
    REAL*8                                      :: eps_plot(Mesh%Nnodes)
    REAL*8, POINTER                             :: h_target_on_nodes(:), nodes_glob(:,:)
    INTEGER, POINTER                            :: connectivity_glob(:,:)
#ifdef PARALL
    REAL*8, ALLOCATABLE                         :: h_root(:), error_oscillation_root(:)
    REAL*8, POINTER                             :: h_glob(:), error_oscillation_glob(:)
    INTEGER, POINTER                            :: vector_nodes_unique_glob(:), count_vec_glob(:)
    INTEGER, ALLOCATABLE                        :: count_vec(:), noghost_index(:)
#endif
    CHARACTER(70)                               :: param_adapt_char, count_adapt_char
    ! mesh_name is the mesh path + mesh name + .msh extension ("./Meshes/CircLim.msh")
    ! mesh_name_npne (mesh name no path no extension) is just the name of the mesh ("CircLim")
    ! new_mesh_name_npne (new mesh name no path no extension) is just the name of the mesh + param_adapt + count_adapt ("CircLim_param2_n1")
    ! buffer is a dummy array to store intermediate mesh names
    CHARACTER(1024), INTENT(IN)                 :: mesh_name
    CHARACTER(1024)                             :: mesh_name_npne, new_mesh_name_npne, buffer
    INTEGER                                     :: ierr

    NULLIFY(h_target_on_nodes, nodes_glob, connectivity_glob)

#ifdef PARALL
    NULLIFY(h_glob, error_oscillation_glob, vector_nodes_unique_glob, count_vec_glob)
#endif

  IF(MPIvar%glob_id .EQ. 0) THEN
     WRITE(*,*) "*************************************************"
     WRITE(*,*) "            ADAPTIVITY INDICATOR                 "
     WRITE(*,*) "*************************************************"
  ENDIF

#ifdef PARALL
    ALLOCATE(noghost_index(Mesh%Nelems-Mesh%nghostelems))
    ! only select the indices of the non-ghost elements
    noghost_index = PACK([(i, i=1, Mesh%Nelems)], Mesh%ghostElems(:) .EQ. 0)
    CALL unique_1D(RESHAPE(Mesh%T(noghost_index,1:refElPol%Nvertices), [SIZE(Mesh%T(noghost_index,1:refElPol%Nvertices),1) * SIZE(Mesh%T(noghost_index,1:refElPol%Nvertices),2)]), vector_nodes_unique)
    CALL gather_1D_vector_int(Mesh%loc2glob_nodes(vector_nodes_unique), vector_nodes_unique_glob, allgather = .FALSE.)
#else
    CALL unique_1D(RESHAPE(Mesh%T(:,1:refElPol%Nvertices), [SIZE(Mesh%T(:,1:refElPol%Nvertices),1) * SIZE(Mesh%T(:,1:refElPol%Nvertices),2)]), vector_nodes_unique)
#endif

    N_n_vertex = SIZE(vector_nodes_unique)

    ALLOCATE(error_oscillation(N_n_vertex))
    ALLOCATE(h(N_n_vertex))
    error_oscillation = 0.
    h = 0.

#ifndef PARALL
    ALLOCATE(h_target(N_n_vertex))
    h_target = 0.
#else
    ALLOCATE(count_vec(N_n_vertex))
    IF(MPIvar%glob_id .EQ. 0) THEN
       N_n_vertex = MAXVAL(vector_nodes_unique_glob)
       ALLOCATE(error_oscillation_root(N_n_vertex))
       error_oscillation_root = 0
       N_n_vertex = SIZE(vector_nodes_unique)
    ENDIF
#endif

    !! use error map to create element size map: h_target
#ifndef PARALL
    CALL h_map(N_n_vertex,Mesh%X(vector_nodes_unique,1:2),Mesh%T(:,1:3),vector_nodes_unique, h)
#else
    CALL h_map(N_n_vertex,Mesh%X(vector_nodes_unique,1:2),Mesh%T(:,1:3),vector_nodes_unique, h, count_vec)

    CALL gather_1D_vector_real(h, h_glob, allgather = .FALSE.)
    CALL gather_1D_vector_int(count_vec, count_vec_glob, allgather = .FALSE.)

    IF(MPIvar%glob_id .EQ. 0) THEN
       ALLOCATE(h_root(MAXVAL(vector_nodes_unique_glob)))
       h_root = 0.
       CALL compute_error_on_vertices_root(h_glob, vector_nodes_unique_glob, count_vec_glob, h_root)
    ENDIF

#endif

    CALL hdg_ShockCapturing_adapt(thresh, eps_plot)
    !! error indicator based on the elemental oscillations
#ifndef PARALL
    CALL read_error(eps_plot, error_oscillation)
    h_target = h  ! Start by assigning h to h_target
    WHERE (ABS(error_oscillation) .GT. 1e-10)
       h_target = h * 0.5
    END WHERE
#else
    CALL read_error(eps_plot, error_oscillation, count_vec)
    CALL gather_1D_vector_real(error_oscillation, error_oscillation_glob, allgather = .FALSE.)

    IF(MPIvar%glob_id .EQ. 0) THEN
       CALL compute_error_on_vertices_root(error_oscillation_glob, vector_nodes_unique_glob, count_vec_glob, error_oscillation_root)

       ! create contribution of estimator and indicator for mmg
       ! if oscillations are detected just use original size/2 otherwise use
       ! estimator calculated with Richardson

       h_target = h_root  ! Start by assigning h to h_target
       WHERE (ABS(error_oscillation_root) .GT. 1e-10)
          h_target = h_root * 0.5
       END WHERE

    ENDIF
#endif

#ifndef PARALL
    ALLOCATE(nodes_glob(SIZE(Mesh%X,1),SIZE(Mesh%X,2)))
    ALLOCATE(connectivity_glob(SIZE(Mesh%T,1),SIZE(Mesh%T,2)))
    nodes_glob = Mesh%X
    connectivity_glob = Mesh%T
#else
   CALL gather_mesh(Mesh,connectivity_glob,nodes_glob)
   DEALLOCATE(vector_nodes_unique)
   CALL unique_1D(vector_nodes_unique_glob,vector_nodes_unique)  
   
    IF(MPIvar%glob_id .EQ. 0) THEN
       N_n_vertex = MAXVAL(vector_nodes_unique_glob)
#endif
      ALLOCATE(h_target_on_nodes(SIZE(Mesh%X,1)))
      h_target_on_nodes = 0.5
      DO i=1,SIZE(vector_nodes_unique,1)
         h_target_on_nodes(vector_nodes_unique(i)) = h_target(i)
      ENDDO

      CALL gmsh_create_from_h_target( h_target_on_nodes, Mesh%X, Mesh%T(:,1:3), order)

#ifdef PARALL
    ENDIF
    ! wait for process 0 to finish writing before loading new mesh
    CALL MPI_BARRIER(mpi_comm_world, ierr)
#endif

    CALL free_mesh

    CALL free_reference_element_pol(refElPol)
    CALL create_reference_element(refElPol,2,order, verbose = 0)

    IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
       CALL load_gmsh_mesh("./res/temp",0)
    ELSE
       CALL load_gmsh_mesh("./res/temp",1)
    ENDIF
    
    CALL free_reference_element_pol(refElPol)
    CALL create_reference_element(refElPol,2,order, verbose = 0)
    CALL mesh_preprocess_serial(ierr)
    CALL read_extended_connectivity('./res/temp.msh')

    Mesh%X = Mesh%X*phys%lscale

    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
       Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
    END IF

    n_el_unstable = 0

    DO i = 1, SIZE(error_oscillation)
       IF(error_oscillation(i) .GT. 1e-10) THEN
          n_el_unstable = n_el_unstable + 1
       ENDIF
    ENDDO

    WRITE(*,'(a, F5.2)') "********** Percentage of oscillating elements on previous mesh: ", REAL(n_el_unstable*100)/REAL(N_n_vertex), "%"

    DEALLOCATE(error_oscillation)
    DEALLOCATE(h)
    DEALLOCATE(vector_nodes_unique)

#ifndef PARALL
    DEALLOCATE(h_target)
#else
    DEALLOCATE(h_glob)
    DEALLOCATE(count_vec_glob)
    DEALLOCATE(error_oscillation_glob)
    DEALLOCATE(vector_nodes_unique_glob)
    DEALLOCATE(noghost_index)
    DEALLOCATE(count_vec)
    IF(ASSOCIATED(h_target_on_nodes)) DEALLOCATE(h_target_on_nodes)    
    NULLIFY(h_target_on_nodes)
    DEALLOCATE(nodes_glob)
    IF(ASSOCIATED(nodes_glob)) DEALLOCATE(nodes_glob) 
    NULLIFY(nodes_glob)
    IF(ASSOCIATED(connectivity_glob)) DEALLOCATE(connectivity_glob)
    NULLIFY(connectivity_glob)
    IF(MPIvar%glob_id .EQ. 0) THEN
       DEALLOCATE(h_target)
       DEALLOCATE(h_root)
       DEALLOCATE(error_oscillation_root)
    ENDIF
    NULLIFY(h_glob, error_oscillation_glob, vector_nodes_unique_glob, count_vec_glob)
#endif


  ENDSUBROUTINE adaptivity_indicator

  SUBROUTINE compute_error_oscillations(oscillations, min_osc, max_osc, n_osc, ir, ir_check, Mesh_prec)
    REAL*8, ALLOCATABLE, INTENT(OUT)  :: oscillations(:)
    REAL*8, INTENT(OUT)               :: min_osc, max_osc
    INTEGER, INTENT(IN)               :: ir
    INTEGER, INTENT(OUT)              :: n_osc,  ir_check
    TYPE(Mesh_type), INTENT(INOUT)    :: Mesh_prec
#ifdef PARALL
    INTEGER                           :: ierr
#endif

    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tps1)
       CALL system_CLOCK(timing%cks1, timing%clock_rate1)
    END IF

    IF(.NOT. ALLOCATED(oscillations)) THEN
       ALLOCATE(oscillations(Mesh%Nelems))
    ELSEIF(SIZE(oscillations) .NE. Mesh%Nelems) THEN
       DEALLOCATE(oscillations)
       ALLOCATE(oscillations(Mesh%Nelems))
    ENDIF
    oscillations = -100.

    CALL check_oscillations(adapt%thr_ind, oscillations)

    max_osc = MAXVAL(oscillations)
    min_osc = MINVAL(oscillations)
    n_osc = COUNT((oscillations .LT. 0.) .AND. (oscillations .GT. -100.))

#ifdef PARALL
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, max_osc, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, min_osc, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_Allreduce(MPI_IN_PLACE, n_osc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) "MAX ERROR OSCILLATION:       ", max_osc
       !WRITE(*,*) "MIN ERROR OSCILLATION:       ", min_osc
       WRITE(*,*) "NUMBER OF OSCILLATIONS:      ", n_osc
    ENDIF

    IF(max_osc .LE. adapt%osc_check) THEN
       IF(MPIvar%glob_id .EQ. 0) THEN
          WRITE(*,*) "Solution saved as checkpoint."
       ENDIF
       IF(SIZE(sol%u_conv) .NE. SIZE(sol%u)) THEN
          DEALLOCATE(sol%u_conv)
          DEALLOCATE(sol%q_conv)
          ALLOCATE(sol%u_conv(SIZE(sol%u)))
          ALLOCATE(sol%q_conv(SIZE(sol%q)))
       ENDIF
       sol%u_conv = sol%u
       sol%q_conv = sol%q
       ir_check = ir
       CALL free_mesh_loc(Mesh_prec)
       CALL deep_copy_mesh_struct(Mesh, Mesh_prec)
    ENDIF



    IF (utils%timing) THEN
       CALL cpu_TIME(timing%tpe1)
       CALL system_CLOCK(timing%cke1, timing%clock_rate1)
       timing%runtadapt = timing%runtadapt + (timing%cke1-timing%cks1)/REAL(timing%clock_rate1)
       timing%cputadapt = timing%cputadapt + timing%tpe1-timing%tps1
    END IF

  ENDSUBROUTINE compute_error_oscillations

  SUBROUTINE check_oscillations(thresh, oscillations)

    REAL*8, INTENT(IN)                                  :: thresh
    REAL*8, OPTIONAL, INTENT(OUT)                       :: oscillations(:)
    REAL*8                                              :: eps_plot(Mesh%Nnodes)

    CALL hdg_ShockCapturing_adapt(thresh, eps_plot, oscillations)

  ENDSUBROUTINE check_oscillations

  SUBROUTINE hdg_ShockCapturing_adapt(thresh, eps_plot, oscillations)

    REAL*8, INTENT(OUT)                                 :: eps_plot(:)
    REAL*8, OPTIONAL, INTENT(OUT)                       :: oscillations(:)
    REAL*8, INTENT(IN)                                  :: thresh
    REAL*8                                              :: eps_elem(Mesh%Nelems), eps_nodal_n(Mesh%Nnodesperelem), eps_nodal_u(RefElPol%Nnodes2D), aux_eps_nodal_n(RefElPol%Nvertices)
    REAL*8                                              :: V1(Mesh%Nelems), V0(Mesh%Nelems), tol = 1e-12
    INTEGER                                             :: Ne, Nv, Nvert, Ncon, i, Nx, inod, iel
    INTEGER                                             :: els(SIZE(Mesh%N, 2))
    REAL*8                                              :: Vand(refElPol%Nnodes2D, refElPol%Nnodes2D), invVand(refElPol%Nnodes2D, refElPol%Nnodes2D)

    Ne = Mesh%Nelems                  ! number of elements
    Nx = SIZE(Mesh%X,1)
    Nv = RefElPol%Nnodes2D          ! number of nodes per element
    Nvert  = RefElPol%Nvertices     ! number of vertices
    Ncon = SIZE(RefElPol%N2D)
    eps_plot = 0.
    eps_elem = 0.
    aux_eps_nodal_n = 0.
    eps_nodal_u     = 0.
    V0 = 0
    V1 = 1


    !******* Find shock capturing coefficient in each element
    ! Vandermonde matrix
    IF (refElPol%elemType == 0) THEN
       ! Triangles
       CALL vandermonde_2d(Vand, refElPol)
    ELSEIF (refElPol%elemType == 1) THEN
       ! Quadrilaterals
       CALL vandermonde_qua(Vand, refElPol)
    ELSE
       WRITE (6, *) "Vandermonde matrix for this element type not coded yet"
       STOP
    END IF
    ! Invert Vandermonde matrix
    CALL invert_matrix(Vand, invVand)

    CALL find_coeff_shock_capturing_adapt(thresh, eps_elem, invVand, oscillations)

    DO iel = 1, Ne
       ! Shock capturing parameter in each node
       SELECT CASE (adapt%shockcp_adapt)
       CASE (1)
          IF (eps_elem(iel) .LT. tol) CYCLE

          ! Constant value in each element
          eps_nodal_n = eps_elem(iel) * 1.0d0
          eps_nodal_u = 0.0d0

          ! Update eps_plot
          eps_plot(Mesh%T(iel, :)) = eps_nodal_n

       CASE (2)
          ! Linear interpolation
          DO inod = 1, refElPol%Nvertices
             els = Mesh%N(Mesh%Tlin(iel, inod), :)
             aux_eps_nodal_n(inod) = 0.
             DO i = 1, SIZE(els)
                IF (els(i) .EQ. 0) CYCLE
                aux_eps_nodal_n(inod) = MAX(aux_eps_nodal_n(inod), eps_elem(els(i)))
             END DO
          END DO

          ! Check IF all aux_eps_nodal_n values are zero
          IF (MAXVAL(aux_eps_nodal_n) .LT. tol) THEN
             CYCLE
          ELSE
             ! Calculate eps_nodal_n using shock_st.N
             eps_nodal_n = MATMUL(refElPol%Nlin, aux_eps_nodal_n)
          ENDIF
          ! Update eps_plot
          eps_plot(Mesh%T(iel, :)) = eps_nodal_n

       CASE DEFAULT
          WRITE(*,*) "Option of shockcp_adapt not allowed. STOP."
          STOP
       END SELECT
    END DO

  ENDSUBROUTINE hdg_ShockCapturing_adapt

  SUBROUTINE find_coeff_shock_capturing_adapt(thresh, eps, invV, oscillations)
    USE physics, ONLY: cons2phys
    REAL*8, INTENT(IN)            :: thresh
    REAL*8, INTENT(OUT)           :: eps(:)
    REAL*8, OPTIONAL, INTENT(OUT) :: oscillations(:)
    REAL*8, INTENT(IN)            :: invV(refElPol%Nnodes2D, refElPol%Nnodes2D)
    INTEGER*4                     :: Ndim, Neq, Nel, Np, Npm1, i, j, counter1, counter2, start, ending
    INTEGER, ALLOCATABLE          :: indices(:)
    REAL*8                        :: se(Mesh%Nelems), s0
    REAL*8, ALLOCATABLE           :: up(:, :),grad_mag(:, :), grad(:, :, :), udet(:)
    REAL*8, ALLOCATABLE           :: um(:, :), umho(:, :)
    REAL*8, PARAMETER             :: tol = 1e-12
    REAL*8                        :: uc1, uc2, uc3, uc4

    Ndim = Mesh%ndim
    Neq = phys%Neq
    Nel = Mesh%Nelems
    Np = refElPol%Nnodes2D

    eps = 0

    ALLOCATE (up(Mesh%Nelems*refElPol%Nnodes2D, phys%npv))
    ALLOCATE (grad(Mesh%Nelems*refElPol%Nnodes2D, phys%npv, Mesh%ndim))
    ALLOCATE (grad_mag(Mesh%Nelems*refElPol%Nnodes2D, phys%npv))
    ALLOCATE (udet(Mesh%Nelems*refElPol%Nnodes2D))
    ALLOCATE (um(refElPol%Nnodes2D, Mesh%Nelems))
    ALLOCATE (umho(refElPol%Nnodes2D, Mesh%Nelems))

    CALL cons2phys(TRANSPOSE(RESHAPE(sol%u, (/neq, Nel*Np/))), up)

    IF(adapt%quant_ind .EQ. 1) THEN
       start = 1
       ending = 1
    ELSEIF(adapt%quant_ind .EQ. 2) THEN
       start = 2
       ending = 2
    ELSEIF(adapt%quant_ind .EQ. 3) THEN
       start = 1
       ending = 2
    ELSE
       WRITE(*,*) "quant_ind not valid, must be between 0 and 3. STOP"
       STOP
    ENDIF


    IF((adapt%quant_ind .EQ. 2) .OR. (adapt%quant_ind .EQ. 3)) THEN
       DO i = 1, SIZE(up,1)
          counter1 = (i-1)*phys%neq+1
          counter2 = (i-1)*phys%neq*mesh%ndim+1

          uc1 = sol%u(counter1)
          uc2 = sol%u(counter1+1)
          uc3 = sol%u(counter1+2)
          uc4 = sol%u(counter1+3)

          grad(i,1,1) = sol%q(counter2)
          grad(i,1,2) = sol%q(counter2+1)
          grad_mag(i,1) = NORM2(grad(i,1,:))

          grad(i,2,1) = -uc2/uc1**2*grad(i,1,1) + 1/uc1*sol%q(counter2+2)
          grad(i,2,2) = -uc2/uc1**2*grad(i,1,2) + 1/uc1*sol%q(counter2+3)
          grad_mag(i,2) = NORM2(grad(i,2,:))

          grad(i,3,1) = -uc3/uc1**2*grad(i,1,1) + 1/uc1*sol%q(counter2+4)
          grad(i,3,2) = -uc3/uc1**2*grad(i,1,2) + 1/uc1*sol%q(counter2+5)
          grad_mag(i,3) = NORM2(grad(i,3,:))

          grad(i,4,1) = -uc4/uc1**2*grad(i,1,1) + 1/uc1*sol%q(counter2+6)
          grad(i,4,2) = -uc4/uc1**2*grad(i,1,2) + 1/uc1*sol%q(counter2+7)
          grad_mag(i,4) = NORM2(grad(i,4,:))

          grad(i,5,1) = 2/(3*phys%Mref)*(sol%q(counter2+4)-0.5/uc1**2*(2*uc1*uc2*sol%q(counter2+2)-uc2**2*grad(i,1,1)))
          grad(i,5,2) = 2/(3*phys%Mref)*(sol%q(counter2+5)-0.5/uc1**2*(2*uc1*uc2*sol%q(counter2+3)-uc2**2*grad(i,1,2)))
          grad_mag(i,5) = NORM2(grad(i,5,:))

          grad(i,6,1) = 2/(3*phys%Mref)*sol%q(counter2+6)
          grad(i,6,2) = 2/(3*phys%Mref)*sol%q(counter2+7)
          grad_mag(i,6) = NORM2(grad(i,6,:))

          grad(i,7,1) = (uc1*sol%q(counter2+4)-uc3*sol%q(counter2))/uc1**2
          grad(i,7,1) = grad(i,7,1) - uc2/uc1**3*(uc1*sol%q(counter2+2)-uc2*sol%q(counter2))
          grad(i,7,1) = grad(i,7,1)*2/(3*phys%Mref)
          grad(i,7,2) = (uc1*sol%q(counter2+5)-uc3*sol%q(counter2+1))/uc1**2
          grad(i,7,2) = grad(i,7,2) - uc2/uc1**3*(uc1*sol%q(counter2+3)-uc2*sol%q(counter2+1))
          grad(i,7,2) = grad(i,7,2)*2/(3*phys%Mref)
          grad_mag(i,7) = NORM2(grad(i,7,:))

          grad(i,8,1) =  2/(3*phys%Mref)*(uc1*sol%q(counter2+6)-uc4*grad(i,1,1))/uc1**2
          grad(i,8,2) =  2/(3*phys%Mref)*(uc1*sol%q(counter2+7)-uc4*grad(i,1,2))/uc1**2
          grad_mag(i,8) = NORM2(grad(i,8,:))

          grad(i,9,1) = phys%Mref**(-0.5)*0.5*(up(i,7)+up(i,8))**(-0.5)*(grad(i,7,1)+grad(i,8,1))
          grad(i,9,2) = phys%Mref**(-0.5)*0.5*(up(i,7)+up(i,8))**(-0.5)*(grad(i,7,2)+grad(i,8,2))
          grad_mag(i,9) = NORM2(grad(i,9,:))

          grad(i,10,1) = 1/up(i,9)*grad(i,2,1)-up(i,2)/up(i,9)**2*grad(i,9,1)
          grad(i,10,2) = 1/up(i,9)*grad(i,2,2)-up(i,2)/up(i,9)**2*grad(i,9,2)
          grad_mag(i,10) = NORM2(grad(i,10,:))

       ENDDO
    ENDIF

    IF(adapt%n_quant_ind .EQ. 0) THEN
       ALLOCATE(indices(phys%npv))
       indices = (/(i,i=1,phys%npv)/)
    ELSEIF((adapt%n_quant_ind .GE. 1) .AND. (adapt%n_quant_ind .LE. 10)) THEN
       ALLOCATE(indices(1))
       indices(1) = adapt%n_quant_ind
    ELSE
       WRITE(*,*) "n_quant_ind not valid, must be between 0 and 10. STOP"
       STOP
    ENDIF

    DO counter1 = start,ending
       DO j = 1, SIZE(indices)

          IF(counter1 .EQ. 1) THEN
             udet = up(:,indices(j))
          ELSE
             udet = grad_mag(:,indices(j))
          ENDIF


          ! Convert solution into modal expansion
          um = 0
          um = MATMUL(invV, RESHAPE(udet, (/Np, Nel/)))

          ! Solution with only the ho mode
          Npm1 = refElPol%Ndeg*(refElPol%Ndeg + 1)/2
          umho = 0.
          umho(Npm1 + 1:Np, :) = um(Npm1 + 1:Np, :)

          ! Shock detector
          se = LOG10(tol + SUM(umho**2, 1)/(SUM(um**2, 1) + tol))

          ! coefficients
          s0 = LOG10(1./refElPol%Ndeg**4)

          DO i = 1, Nel
             IF (SUM(um(:, i)**2) .LT. thresh) THEN
                se(i) = -100.
             END IF
             IF (se(i) .GT. s0) THEN
                eps(i) = 1
                IF(PRESENT(oscillations)) THEN
                   oscillations(i) = MAX(oscillations(i),se(i))
                ENDIF
             END IF
          END DO
       ENDDO
    ENDDO

    DEALLOCATE (up, udet, um, umho)
    DEALLOCATE (grad_mag, grad)
    DEALLOCATE(indices)

  END SUBROUTINE find_coeff_shock_capturing_adapt

  SUBROUTINE read_error(eps_plot, error_oscillation, count_vec)
    REAL*8, INTENT(OUT)                                 :: error_oscillation(:)
    REAL*8, INTENT(IN)                                  :: eps_plot(:)
    INTEGER, INTENT(OUT), OPTIONAL                      :: count_vec(:)
    REAL*8, ALLOCATABLE                                 :: error_vec(:)
    INTEGER, ALLOCATABLE                                :: count_vec_local(:)
    INTEGER                                             :: i,j, ind, N_e_real, N_n_max
    INTEGER, ALLOCATABLE                                :: vertex_nodes(:,:)


    N_e_real = Mesh%Nelems
    N_n_max = Mesh%Nelems*Mesh%Nnodesperelem

    ALLOCATE(vertex_nodes(Mesh%Nelems,refElPol%Nvertices))
    ALLOCATE(error_vec(N_n_max))
    ALLOCATE(count_vec_local(N_n_max))
    error_vec = 0.
    count_vec_local = 0
    vertex_nodes = Mesh%T(:,1:refElPol%Nvertices)

    DO i =1,N_e_real
#ifdef PARALL
       IF(Mesh%ghostElems(i) .EQ. 1) CYCLE
#endif
       DO j=1,RefElPol%Nvertices
          ind = vertex_nodes(i,j)
          error_vec(ind) = error_vec(ind) + eps_plot(ind)
          count_vec_local(ind) = count_vec_local(ind) + 1
       ENDDO
    ENDDO

#ifndef PARALL
    j=1
    DO i=1,N_n_max
       IF (count_vec_local(i) .NE. 0) THEN
          error_oscillation(j) = error_vec(i)/count_vec_local(i)
          j=j+1
       ENDIF
    ENDDO
#else
    j=1
    DO i=1,N_n_max
       IF (count_vec_local(i) .NE. 0) THEN
          error_oscillation(j) = error_vec(i)
          count_vec(j) = count_vec_local(i)
          j=j+1
       ENDIF
    ENDDO
#endif

    DEALLOCATE(error_vec)
    DEALLOCATE(count_vec_local)
    DEALLOCATE(vertex_nodes)

  ENDSUBROUTINE read_error

  PURE SUBROUTINE unique_2D(input_matrix, output_matrix)
    INTEGER, INTENT(IN)                :: input_matrix(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT)  :: output_matrix(:,:)
    INTEGER                            :: num_rows, num_cols, i, j, k, count
    LOGICAL, ALLOCATABLE               :: is_unique(:)

    ! Determine the number of rows and columns in the input matrix
    num_rows = SIZE(input_matrix, 1)
    num_cols = SIZE(input_matrix, 2)

    ! Initialize an array to track unique elements along the second dimension
    ALLOCATE(is_unique(num_cols))
    is_unique = .TRUE.

    ! Initialize the output matrix
    ALLOCATE(output_matrix(num_rows, num_cols))

    ! Loop through each row and remove duplicate elements along the second dimension
    DO i = 1, num_rows
       ! Reset is_unique array for each row
       is_unique = .TRUE.
       count = 0

       DO j = 1, num_cols
          IF (is_unique(j)) THEN
             count = count + 1
             output_matrix(i, count) = input_matrix(i, j)

             ! Check for duplicates in the rest of the row
             DO k = j + 1, num_cols
                IF (input_matrix(i, j) == input_matrix(i, k)) THEN
                   is_unique(k) = .FALSE.
                END IF
             END DO
          END IF
       END DO
    END DO

    ! Deallocate the temporary array
    DEALLOCATE(is_unique)

  END SUBROUTINE unique_2D

ENDMODULE adaptivity_indicator_module
