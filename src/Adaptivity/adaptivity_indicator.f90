!************************************************************
! project: MHDG
! file: inout.f90
! date: 06/09/2016
! Module for schock capturing adaptivity
!************************************************************

MODULE adaptivity_indicator_module
  USE globals
  USE reference_element
  USE LinearAlgebra
  USE gmsh
  USE adaptivity_common_module
  IMPLICIT NONE

CONTAINS

  SUBROUTINE adaptivity_indicator(mesh_name,thresh, param_adapt, count_adapt, order)
    USE in_out, ONLY: copy_file
    USE gmsh_io_module, ONLY: load_gmsh_mesh, HDF5_save_mesh, convert_gmsh_to_hdf5
    USE preprocess

    TYPE(gmsh_t)                                :: gmsh
    REAL*8, INTENT(IN)                          :: thresh
    INTEGER, INTENT(IN)                         :: param_adapt, count_adapt, order

    REAL*8,  ALLOCATABLE                        :: two_d_nodes(:,:)
    INTEGER, ALLOCATABLE                        :: two_d_elements(:,:)
    REAL*8,  ALLOCATABLE                        :: h(:), h_target(:), error_oscillation(:)
    INTEGER, ALLOCATABLE                        :: vector_nodes_unique(:)
    INTEGER                                     :: i, N_n_vertex, n_el_unstable
    REAL*8                                      :: eps_plot(Mesh%Nnodes)

    CHARACTER(70)                               :: param_adapt_char, count_adapt_char
    ! mesh_name is the mesh path + mesh name + .msh extension ("./Meshes/CircLim.msh")
    ! mesh_name_npne (mesh name no path no extension) is just the name of the mesh ("CircLim")
    ! new_mesh_name_npne (new mesh name no path no extension) is just the name of the mesh + param_adapt + count_adapt ("CircLim_param2_n1")
    ! buffer is a dummy array to store intermediate mesh names
    CHARACTER(1024), INTENT(IN)                 :: mesh_name
    CHARACTER(1024)                             :: mesh_name_npne, new_mesh_name_npne, buffer
    INTEGER                                     :: ierr


    WRITE(*,*) '*************** Starting refinement procedure with oscillations: ****************'

    CALL hdg_ShockCapturing_adapt(thresh, eps_plot)

    IF(MPIvar%glob_id .EQ. 0) THEN
       WRITE(*,*) '******************************** Unique: *****************************************'
    ENDIF

    CALL unique_1D(RESHAPE(Mesh%T(:,1:RefElPol%Nvertices), [SIZE(Mesh%T(:,1:RefElPol%Nvertices),1) * SIZE(Mesh%T(:,1:RefElPol%Nvertices),2)]), vector_nodes_unique)

    N_n_vertex = SIZE(vector_nodes_unique)

    ALLOCATE(two_d_nodes(N_n_vertex,2))
    ALLOCATE(two_d_elements(SIZE(Mesh%T,1),3))
    ALLOCATE(error_oscillation(N_n_vertex))
    ALLOCATE(h(N_n_vertex))
    ALLOCATE(h_target((N_n_vertex)))

    two_d_nodes = 0.
    two_d_elements = 0.
    error_oscillation = 0.
    h = 0.

    two_d_nodes = Mesh%X(vector_nodes_unique,:)
    two_d_elements = Mesh%T(:,1:3)

    !! error indicator based on the elemental oscillations
    CALL read_error(eps_plot, error_oscillation)

    !! use error map to create element size map: h_target
    CALL h_map(N_n_vertex,two_d_nodes,two_d_elements,vector_nodes_unique, h);!%*1.901*1e-3;

    ! create contribution of estimator and indicator for mmg
    ! if oscillations are detected just use original size/2 otherwise use
    ! estimator calculated with Richardson

    DO i=1, SIZE(h) !here we add the oscillation contribution
       IF (ABS(error_oscillation(i)) .GT. 1e-10) THEN
          ! refine dividing by 2 the element size when find oscillations
          h_target(i) = h(i)*0.5;
       ELSE
          h_target(i) = h(i)
       ENDIF
    ENDDO

    CALL generate_htarget_sol_file(N_n_vertex, h_target)

    CALL extract_mesh_name_from_fullpath_woext(mesh_name, mesh_name_npne)

    WRITE(param_adapt_char, *) param_adapt
    WRITE(count_adapt_char, *) count_adapt
    new_mesh_name_npne = TRIM(ADJUSTL(mesh_name_npne)) // '_param'// TRIM(ADJUSTL(param_adapt_char)) // '_n' // TRIM(ADJUSTL(count_adapt_char))

    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".mesh"
    CALL mmg_create_mesh_from_h_target(buffer)

    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne))
    CALL convert_mesh2msh(buffer)
    CALL convert_msh2mesh(buffer)

    CALL delete_file("./res/temp.mesh")
    CALL delete_file("./res/temp.msh")
    CALL delete_file("./res/ElSizeMap.sol")


    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".mesh"
    CALL copy_file(buffer, "./res/temp.mesh")
    !CALL delete_file(buffer)

    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".msh"


    CALL open_merge_with_geometry(gmsh, buffer)
    CALL copy_file(buffer, "./res/temp.msh")
    !CALL delete_file(buffer)


    buffer = "./res/" // TRIM(ADJUSTL(new_mesh_name_npne)) // ".sol"
    CALL delete_file(buffer)

    !CALL merge_with_geometry(gmsh)

    WRITE(*,*) "********** Loading mesh P1  **********"
    CALL free_mesh

    IF((switch%testcase .GE. 60) .AND. (switch%testcase .LE. 80)) THEN
       CALL load_gmsh_mesh("./res/temp",0)
    ELSE
       CALL load_gmsh_mesh("./res/temp",1)
    ENDIF

    CALL create_reference_element(refElPol,2,1, verbose = 0)
    CALL mesh_preprocess(ierr)

    Mesh%X = Mesh%X*phys%lscale

    IF(ierr .EQ. 0) THEN
       WRITE(*,*) "Error! Corresponding face in Tb not found. STOP"
       STOP
    ENDIF

    CALL HDF5_save_mesh("./newmesh_pre.h5", Mesh%Ndim, mesh%Nelems, mesh%Nextfaces, mesh%Nnodes, mesh%Nnodesperelem, mesh%Nnodesperface, mesh%elemType, mesh%T, mesh%X, mesh%Tb, mesh%boundaryFlag)

    CALL read_extended_connectivity('./res/temp.msh')

    WRITE(*,*) "********** Increasing order mesh **********"
    CALL set_order_mesh(order)
    CALL free_reference_element
    CALL create_reference_element(refElPol,2,order, verbose = 0)
    CALL mesh_preprocess(ierr)

    Mesh%X = Mesh%X*phys%lscale

    IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
       Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
    END IF

    CALL HDF5_save_mesh("./newmesh_notround.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)

    WRITE(*,*) "********** Rounding edges **********"
    !CALL round_edges(Mesh)

    ! overwrite the temp.msh file with the new one with rounded edges (still order 1)
    CALL write_msh_file(Mesh%X,Mesh%T)
    ! convert the mesh to .mesh
    CALL convert_msh2mesh('./res/temp')

    CALL HDF5_save_mesh("./newmesh_round.h5", Mesh%Ndim, Mesh%Nelems, Mesh%Nextfaces, Mesh%Nnodes, Mesh%Nnodesperelem, Mesh%Nnodesperface, Mesh%elemType, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)

    n_el_unstable = 0

    DO i = 1, SIZE(error_oscillation)
       IF(error_oscillation(i) .GT. 1e-10) THEN
          n_el_unstable = n_el_unstable + 1
       ENDIF
    ENDDO

    WRITE(*,'(a, F5.2)') "********** Percentage of oscillating elements on previous mesh: ", REAL(n_el_unstable*100)/REAL(N_n_vertex), "%"

    DEALLOCATE(two_d_nodes)
    DEALLOCATE(two_d_elements)
    DEALLOCATE(error_oscillation)
    DEALLOCATE(h)
    DEALLOCATE(vector_nodes_unique)

  ENDSUBROUTINE adaptivity_indicator

  SUBROUTINE check_oscillations(thresh, error_oscillation, oscillations)

    REAL*8, INTENT(IN)                                  :: thresh
    REAL*8, INTENT(OUT)                                 :: error_oscillation(:)
    REAL*8, OPTIONAL, INTENT(OUT)                       :: oscillations(:)
    REAL*8                                              :: eps_plot(Mesh%Nnodes)


    CALL hdg_ShockCapturing_adapt(thresh, eps_plot, oscillations)
    CALL read_error(eps_plot, error_oscillation)

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

       END SELECT
    END DO

  ENDSUBROUTINE hdg_ShockCapturing_adapt

  SUBROUTINE find_coeff_shock_capturing_adapt(thresh, eps, invV, oscillations)
    USE physics, ONLY: cons2phys
    REAL*8, INTENT(IN)      :: thresh
    REAL*8, INTENT(OUT)     :: eps(:)
    REAL*8, OPTIONAL, INTENT(OUT)                       :: oscillations(:)
    REAL*8, INTENT(IN)      :: invV(refElPol%Nnodes2D, refElPol%Nnodes2D)
    INTEGER*4               :: Ndim, Neq, Nel, Np, Npm1, i, j, counter1, counter2, start, END
    INTEGER, ALLOCATABLE    :: indices(:)
    REAL*8                  :: se(Mesh%Nelems), s0
    REAL*8, ALLOCATABLE     :: up(:, :),grad_mag(:, :), grad(:, :, :), udet(:)
    REAL*8, ALLOCATABLE     :: um(:, :), umho(:, :)
    REAL*8, PARAMETER       :: tol = 1e-12
    REAL*8                  :: uc1, uc2, uc3, uc4

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
      end = 1
    ELSEIF(adapt%quant_ind .EQ. 2) THEN
      start = 2
      end = 2
    ELSEIF(adapt%quant_ind .EQ. 3) THEN
      start = 1
      end = 2
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

    DO counter1 = start,END
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
                se(i) = -100
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

    ! IF(PRESENT(oscillations)) THEN
    !   WRITE(time_buffer, *) time%it
    !   buffer = './fortran_save_'// trim(adjustl(time_buffer)) // '.h5'
    !   call HDF5_create(trim(adjustl(buffer)), file_id, ierr)
    !   call HDF5_array2D_saving(file_id, um, size(um,1), size(um,2), 'sol_f')
    !   call HDF5_array2D_saving(file_id, umho, size(umho,1), size(umho,2), 'sol_ho_f')
    !   call HDF5_array1D_saving(file_id, oscillations, size(oscillations), 'oscillations_f')
    !   call HDF5_array2D_saving_int(file_id, Mesh%T, size(Mesh%T,1), size(Mesh%T,2), 'T_f')
    !   call HDF5_array2D_saving(file_id, Mesh%X, size(Mesh%X,1), size(Mesh%X,2), 'X_f')
    !   call HDF5_array2D_saving(file_id, up, size(up,1), size(up,2), 'up_f')
    !   call HDF5_array2D_saving(file_id, grad_mag, size(grad_mag,1), size(grad_mag,2), 'grad_mag_f')
    !   call HDF5_array1D_saving(file_id, sol%u, size(sol%u,1), 'u_f')
    !   call HDF5_close(file_id)
    ! ENDIF

    DEALLOCATE (up, udet, um, umho)
    DEALLOCATE (grad_mag, grad)
    DEALLOCATE(indices)
  END SUBROUTINE find_coeff_shock_capturing_adapt

  SUBROUTINE read_error(eps_plot, error_oscillation)
    REAL*8, INTENT(OUT)                                 :: error_oscillation(:)
    REAL*8, INTENT(IN)                                  :: eps_plot(:)
    REAL*8, ALLOCATABLE                                 :: error_vec(:)
    INTEGER, ALLOCATABLE                                :: count_vec(:)
    INTEGER                                             :: i,j, ind, N_e_real, N_n_max
    INTEGER                                             :: vertex_nodes(Mesh%Nelems,refElPol%Nvertices)


    N_e_real = Mesh%Nelems
    N_n_max = Mesh%Nelems*Mesh%Nnodesperelem

    ALLOCATE(error_vec(N_n_max))
    ALLOCATE(count_vec(N_n_max))

    error_vec = 0.
    count_vec = 0
    vertex_nodes = Mesh%T(:,1:refElPol%Nvertices)

    DO i =1,N_e_real
       DO j=1,RefElPol%Nvertices
          ind = vertex_nodes(i,j)
          error_vec(ind) = error_vec(ind) + eps_plot(ind)
          count_vec(ind) = count_vec(ind) + 1
       ENDDO
    ENDDO

    j=1
    DO i=1,N_n_max
       IF (count_vec(i) .NE. 0) THEN
          error_oscillation(j) = error_vec(i)/count_vec(i)
          j=j+1;
       ENDIF
    ENDDO

    DEALLOCATE(error_vec)
    DEALLOCATE(count_vec)
  ENDSUBROUTINE read_error


  SUBROUTINE h_map(N_n_vertex, two_d_nodes, two_d_elements, vector_nodes_unique, h)
    INTEGER, INTENT(IN)              :: N_n_vertex
    INTEGER                          :: N_e_real
    REAL*8, INTENT(IN)               :: two_d_nodes(:,:)
    INTEGER, INTENT(IN)              :: two_d_elements(:,:), vector_nodes_unique(:)
    REAL*8, INTENT(OUT)              :: h(N_n_vertex)
    INTEGER                          :: count_vec(N_n_vertex)
    REAL*8                           :: g(N_n_vertex), l(N_n_vertex)
    INTEGER                          :: i, j, A, B, C
    REAL*8, DIMENSION(2, 2)          :: J1, J2, J3
    REAL*8, DIMENSION(2, 2)          :: M1, M2, M3
    REAL*8                           :: detJ1, detJ2, detJ3
    REAL*8                           :: sqrt3

    sqrt3 = SQRT(3.0)
#ifndef PARALL
    N_e_real = SIZE(Mesh%T,1)
#else
    N_e_real = SIZE(Mesh%T,1) - Mesh%nghostelems
#endif

    ! Initialize arrays
    g = 0.
    l = 0.
    count_vec = 0
    J1 = 0.
    J2 = 0.
    J3 = 0.

    DO i = 1, N_e_real

       ! Find indexes of the vertex nodes
       A = 0
       B = 0
       C = 0

       DO j = 1, N_n_vertex
          IF ((two_d_elements(i, 1) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             A = j
          ENDIF
          IF ((two_d_elements(i, 2) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             B = j
          ENDIF
          IF ((two_d_elements(i, 3) + 1) .EQ. (vector_nodes_unique(j) + 1)) THEN
             C = j
          ENDIF
       ENDDO

       IF(A*B*C .EQ. 0) THEN
          WRITE(*,*) "Index not found in h_map. STOP"
          STOP
       ENDIF

       ! Calculate Jacobian matrices
       CALL jacobian(two_d_nodes, A, B, C, J1)
       CALL jacobian(two_d_nodes, A, B, C, J2)
       CALL jacobian(two_d_nodes, A, B, C, J3)

       ! Calculate determinants of Jacobians
       detJ1 = J1(1, 1) * J1(2, 2) - J1(1, 2) * J1(2, 1)
       detJ2 = J2(1, 1) * J2(2, 2) - J2(1, 2) * J2(2, 1)
       detJ3 = J3(1, 1) * J3(2, 2) - J3(1, 2) * J3(2, 1)

       ! Calculate metrics
       M1 = MATMUL(TRANSPOSE(J1), J1)
       M2 = MATMUL(TRANSPOSE(J2), J2)
       M3 = MATMUL(TRANSPOSE(J3), J3)

       ! Using Jacobian
       g(A) = g(A) + SQRT(2. * detJ1 / sqrt3)
       count_vec(A) = count_vec(A) + 1

       g(B) = g(B) + SQRT(2. * detJ2 / sqrt3)
       count_vec(B) = count_vec(B) + 1

       g(C) = g(C) + SQRT(2. * detJ3 / sqrt3)
       count_vec(C) = count_vec(C) + 1
    ENDDO

    ! Calculate h values
    DO i = 1, N_n_vertex
       IF (count_vec(i) .GT. 0) THEN
          h(i) = g(i) / count_vec(i)
       ELSE
          WRITE(*,*) "GOT A DIVISION BY 0 IN h_map ADAPTIVITY, SOMETHING IS WRONG!"
          STOP
       ENDIF
    ENDDO
  END SUBROUTINE h_map

  SUBROUTINE jacobian(two_d_nodes, A, B, C, J)
    REAL*8, INTENT(IN)              :: two_d_nodes(:,:)
    REAL*8, INTENT(OUT)             :: J(2,2)
    INTEGER, INTENT(IN)             :: A, B, C

    ! Calculate Jacobian matrix
    J(1,1) = two_d_nodes(B,1) - two_d_nodes(A,1)
    J(2,1) = two_d_nodes(B,2) - two_d_nodes(A,2)
    J(1,2) = two_d_nodes(C,1) - two_d_nodes(A,1)
    J(2,2) = two_d_nodes(C,2) - two_d_nodes(A,2)
  END SUBROUTINE jacobian


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
