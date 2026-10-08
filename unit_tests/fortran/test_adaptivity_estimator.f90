PROGRAM test_adaptivity_estimator
  USE globals
  USE MPI_OMP, ONLY: MPIvar
  USE reference_element, ONLY: create_reference_element, free_reference_element_pol
  USE adaptivity_estimator_module, ONLY: hdg_post_process_matrix, hdg_postprocess_solution, apply_estimator
  USE mpi, ONLY: MPI_INIT, MPI_FINALIZE, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_COMM_WORLD
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  IMPLICIT NONE

  TYPE(Reference_element_type) :: higher
  REAL*8, ALLOCATABLE :: K(:,:,:), B(:,:,:), mean_matrix(:,:,:)
  REAL*8, ALLOCATABLE :: u(:), q(:), reconstructed(:), interpolated(:), expected(:)
  REAL*8, ALLOCATABLE :: saved_gradient(:)
  REAL*8 :: h(1), target(1), scale, density, gradient(5), coefficients(5)
  INTEGER :: neq, sizes(4), n, ns, i, j, case_id, axisym_case, scale_case, norm_case, ierr

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPIvar%glob_id, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPIvar%glob_size, ierr)
  IF (MPIvar%glob_size /= 1) ERROR STOP 'estimator unit driver requires one MPI rank'
  utils%printint = 0
  Mesh%ndim = 2
  Mesh%elemType = 0
  Mesh%Nelems = 1
  Mesh%Nnodes = 15
  Mesh%Nnodesperelem = 15
  ALLOCATE(Mesh%X(15,2), Mesh%T(1,15))
#ifdef PARALL
  ALLOCATE(Mesh%ghostElems(1))
  Mesh%ghostElems = 0
#endif
  Mesh%T(1,:) = [(i,i=1,15)]
  CALL create_reference_element(refElPol, 2, 4, verbose=0)
  CALL create_reference_element(higher, 2, 5, verbose=0)
  n = refElPol%Nnodes2D
  ns = higher%Nnodes2D
  Mesh%X(:,1) = 2.d0 + 0.0282842712474619d0*(refElPol%coord2D(:,1)+1.d0)
  Mesh%X(:,2) = 0.0282842712474619d0*(refElPol%coord2D(:,2)+1.d0)

  ! Every equation has a distinct affine field and gradient. The two-equation
  ! layout must remain correct, and other equation counts must not mix strides.
  sizes = [5,2,6,1]
  DO axisym_case = 0,1
    switch%axisym = axisym_case == 1
    DO case_id = 1,SIZE(sizes)
      neq = sizes(case_id)
      phys%neq = neq
      ALLOCATE(K(ns*neq,ns*neq,1), B(ns*neq,ns*neq*2,1), mean_matrix(ns*neq,neq,1))
      ALLOCATE(u(n*neq),q(n*neq*2),reconstructed(ns*neq),interpolated(ns*neq),expected(ns*neq))
      DO i=1,n
        DO j=1,neq
          u((i-1)*neq+j) = 2.d0*j + 0.3d0*j*Mesh%X(i,1) - 0.1d0*(j+1)*Mesh%X(i,2)
          q((i-1)*neq*2+(j-1)*2+1) = 0.3d0*j
          q((i-1)*neq*2+(j-1)*2+2) = -0.1d0*(j+1)
        END DO
      END DO
      DO i=1,ns
        DO j=1,neq
          expected((i-1)*neq+j) = 2.d0*j + 0.3d0*j* &
            (2.d0+0.0282842712474619d0*(higher%coord2D(i,1)+1.d0)) - &
            0.1d0*(j+1)*0.0282842712474619d0*(higher%coord2D(i,2)+1.d0)
        END DO
      END DO
      CALL hdg_post_process_matrix(Mesh%X,Mesh%T,higher,4,5,K,B,mean_matrix)
      CALL hdg_postprocess_solution(q,u,K,B,mean_matrix,higher,refElPol,1,reconstructed,interpolated)
      IF (.NOT. ALL(ieee_is_finite(reconstructed))) ERROR STOP 'nonfinite estimator reconstruction'
      IF (MAXVAL(ABS(reconstructed-expected)) > 2.d-10) THEN
        WRITE(*,*) 'Affine reconstruction failed: equations=',neq,' axisymmetric=',switch%axisym, &
          ' max error=',MAXVAL(ABS(reconstructed-expected))
        ERROR STOP 'estimator gradient equation layout'
      ENDIF
      DEALLOCATE(K,B,mean_matrix,u,q,reconstructed,interpolated,expected)
    END DO
  END DO

  ! The production estimator receives X in metres and q with respect to x/L.
  ! This varying density with constant primitive velocity/temperatures is
  ! exactly reconstructible; all requested fields must permit the 10 cm cap.
  phys%neq = 5
  phys%npv = 11
  phys%Mref = 12.d0
  phys%idx_rhon_eq = 5
  phys%idx_rhon_pv = 11
  switch%axisym = .TRUE.
  ALLOCATE(sol%u(n*5),sol%q(n*10),saved_gradient(n*10),adapt%param_est(4))
  adapt%param_est = [1,2,7,8]
  adapt%tol_est = 0.05d0
  coefficients = [1.d0,0.2d0,25.02d0,20.d0,0.01d0]
  gradient = 8.d0*coefficients
  DO i=1,n
    density = 1.d0 + 8.d0*(Mesh%X(i,1)-2.d0) - 4.d0*Mesh%X(i,2)
    sol%u((i-1)*5+1:i*5) = density*coefficients
  END DO
  h = 0.08d0
  DO scale_case=1,3
    SELECT CASE(scale_case)
    CASE(1)
      scale = 1.d0
    CASE(2)
      scale = 0.001901d0
    CASE(3)
      scale = 0.03d0
    END SELECT
    phys%lscale = scale
    DO i=1,n
      DO j=1,5
        sol%q((i-1)*10+(j-1)*2+1) = scale*gradient(j)
        sol%q((i-1)*10+(j-1)*2+2) = -0.5d0*scale*gradient(j)
      END DO
    END DO
    saved_gradient = sol%q
    DO norm_case=0,1
      adapt%difference = norm_case
      CALL apply_estimator(h,4,target)
      IF (.NOT. ALL(ieee_is_finite(target))) ERROR STOP 'nonfinite estimator target'
      IF (ABS(target(1)-0.1d0) > 1.d-12) THEN
        WRITE(*,*) 'Affine target failed: length scale=',scale,' difference=',norm_case,' target=',target
        ERROR STOP 'estimator length-unit consistency and 10 cm cap'
      ENDIF
      IF (ANY(sol%q /= saved_gradient)) ERROR STOP 'estimator changed stored solver gradients'
    END DO
  END DO
  CALL free_reference_element_pol(higher)
  CALL free_reference_element_pol(refElPol)
  CALL MPI_FINALIZE(ierr)
  WRITE(*,*) 'adaptivity estimator: equation layouts, length scales, physical fields and cap PASS'
END PROGRAM test_adaptivity_estimator
