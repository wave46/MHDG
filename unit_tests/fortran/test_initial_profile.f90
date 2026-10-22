#if !defined(TEMPERATURE) || defined(TOR3D)
#error This check requires a 2D temperature model.
#endif
PROGRAM test_initial_profile
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, Mesh, refElPol, switch
  USE analytical, ONLY: analytical_solution
  IMPLICIT NONE
  REAL*8 :: normalized(3), expected(3), x(3), y(3), axis, orientation
  REAL*8, ALLOCATABLE :: u(:,:), boundary(:,:)
  INTEGER :: sign_case, offset_case

  phys%neq = 4
  phys%npv = 10
  phys%idx_rhon_eq = 0
  phys%idx_rhon_pv = 0
  phys%idx_gamman_eq = 0
  phys%idx_un_pv = 0
  phys%idx_k_eq = 0
  phys%idx_k_pv = 0
#ifdef NEUTRAL
  phys%neq = 5
  phys%npv = 11
  phys%idx_rhon_eq = 5
  phys%idx_rhon_pv = 11
#ifdef NEUTRALGAMMA
  phys%neq = 6
  phys%npv = 12
  phys%idx_gamman_eq = 6
  phys%idx_un_pv = 12
#endif
#ifdef KEQUATION
  phys%neq = phys%neq+1
  phys%npv = phys%npv+1
  phys%idx_k_eq = phys%neq
  phys%idx_k_pv = phys%npv
#endif
#endif
  refElPol%Nnodes2d = 3
  ALLOCATE(Mesh%T(1,3), phys%magnetic_flux(3), phys%magnetic_psi(3))
  ALLOCATE(u(3,phys%neq), boundary(1,phys%neq))
  Mesh%T(1,:) = [3,1,2]
  Mesh%xmin = 1.d0
  Mesh%xmax = 2.d0
  Mesh%ymin = -1.d0
  Mesh%ymax = 1.d0
  x = [1.2d0,1.5d0,1.8d0]
  y = 0.d0
  normalized = [0.001d0,0.8d0,1.6d0]
  phys%magnetic_psi = normalized
  expected = EXP(-normalized(Mesh%T(1,:))**2/(2.d0*0.7d0**2))
  switch%testcase = 54

  DO sign_case = 1,2
    orientation = (-1.d0)**sign_case
    DO offset_case = 1,2
      axis = 3.d0+20.d0*(offset_case-1)
      phys%magnetic_flux = axis+orientation*1.2d0*normalized
      phys%Flux2Dmin = axis
      phys%Flux2Dmax = MAXVAL(phys%magnetic_flux)
      CALL analytical_solution(1,x,y,u)
      IF (.NOT. ALL(ieee_is_finite(u))) ERROR STOP 'nonfinite cold start'
      IF (ANY(u(:,1) <= 0.d0)) ERROR STOP 'zero-density cold start'
      IF (MAXVAL(ABS(u(:,1)-expected)) > 1.d-12) ERROR STOP 'normalized Gaussian density'
      IF (MAXVAL(ABS(u(:,3)-18.d0*expected**2)) > 1.d-12 .OR. &
          MAXVAL(ABS(u(:,4)-18.d0*expected**2)) > 1.d-12) ERROR STOP 'Gaussian energies'
    ENDDO
  ENDDO

  ! Boundary evaluation keeps its constant state rather than using nodal flux.
  CALL analytical_solution(1,x(:1),y(:1),boundary)
  IF (ABS(boundary(1,1)-1.d0) > 1.d-12 .OR. &
      ANY(ABS(boundary(1,3:4)-18.d0) > 1.d-12)) ERROR STOP 'boundary profile changed'
  switch%testcase = 60
  CALL analytical_solution(1,x,y,u)
  IF (ANY(ABS(u(:,1)-1.d0) > 1.d-12) .OR. &
      ANY(ABS(u(:,3:4)-18.d0) > 1.d-12)) ERROR STOP 'other initializer changed'
  WRITE(*,*) 'Initial profile: PASS; both flux signs, flux offsets, core/SOL nodes and unchanged boundary/constant cases'
END PROGRAM test_initial_profile
