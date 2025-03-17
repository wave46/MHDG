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

  SUBROUTINE apply_indicator(h_map_elements,h_target_elements)
      REAL*8, INTENT(IN)    :: h_map_elements(:)
      REAL*8, INTENT(INOUT)   :: h_target_elements(:)
      REAL*8                :: eps_element(SIZE(Mesh%T,1))
      REAL*8                :: oscillations(SIZE(Mesh%T,1))

      oscillations = 0.
      CALL find_oscillations_elements(eps_element,oscillations)

      CALL refine_h_map(h_map_elements,eps_element,h_target_elements)
      
   ENDSUBROUTINE apply_indicator
     

  SUBROUTINE find_oscillations_elements(eps_element,oscillations)
   REAL*8, INTENT(OUT)           :: eps_element(:)
   REAL*8, INTENT(OUT), OPTIONAL :: oscillations(:)
   REAL*8                        :: Vand(refElPol%Nnodes2D, refElPol%Nnodes2D), invVand(refElPol%Nnodes2D, refElPol%Nnodes2D)

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

   CALL find_coeff_shock_capturing_adapt(adapt%thr_ind, eps_element, invVand, oscillations)

   END SUBROUTINE find_oscillations_elements

   SUBROUTINE refine_h_map(h_map_elements, eps_element,h_target_elements)
      REAL*8, INTENT(IN) :: h_map_elements(:)
      REAL*8, INTENT(IN)    :: eps_element(:)
      REAL*8, INTENT(INOUT)   :: h_target_elements(:)
      INTEGER               :: unstable_elements
      INTEGER               :: i
  
      unstable_elements = 0
      h_target_elements = h_map_elements
  
      DO i = 1, SIZE(h_target_elements)
          SELECT CASE (adapt%shockcp_adapt)
          CASE (1)
              CALL refine_if_oscillating(h_target_elements(i), eps_element(i), unstable_elements)
          CASE (2)
              CALL refine_if_neighbors_oscillating(h_target_elements(i), eps_element, i, unstable_elements)
          CASE DEFAULT
              WRITE(*,*) "Option of shockcp_adapt not allowed. STOP."
              STOP
          END SELECT
      ENDDO
  
      WRITE(*,'(A, F5.2, A)') "********** Percentage of refined elements on previous mesh by indicator: ", REAL(unstable_elements*100)/REAL(SIZE(h_map_elements)), "%"
  
  END SUBROUTINE refine_h_map
  
  SUBROUTINE refine_if_oscillating(h_map_element, eps_element, unstable_elements)
      REAL*8, INTENT(INOUT) :: h_map_element
      REAL*8, INTENT(IN)    :: eps_element
      INTEGER, INTENT(INOUT) :: unstable_elements
  
      IF (eps_element .GT. 1e-10) THEN
          h_map_element = h_map_element * 0.5
          unstable_elements = unstable_elements + 1
      END IF
  END SUBROUTINE refine_if_oscillating
  
  SUBROUTINE refine_if_neighbors_oscillating(h_map_element, eps_element, i, unstable_elements)
      REAL*8, INTENT(INOUT) :: h_map_element
      REAL*8, INTENT(IN)    :: eps_element(:)
      INTEGER, INTENT(IN)   :: i
      INTEGER, INTENT(INOUT) :: unstable_elements
      INTEGER               :: inod, els(SIZE(Mesh%N, 2))
  
      DO inod = 1, refElPol%Nvertices
          els = Mesh%N(Mesh%Tlin(i, inod), :)
          IF (ANY(eps_element(PACK(els, els /= 0)) .GT. 1e-10)) THEN
              h_map_element = h_map_element * 0.5
              unstable_elements = unstable_elements + 1
              EXIT
          END IF
      END DO
  END SUBROUTINE refine_if_neighbors_oscillating

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

    CALL check_oscillations(oscillations)

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

  SUBROUTINE check_oscillations(oscillations)

    REAL*8, OPTIONAL, INTENT(OUT)                       :: oscillations(:)
    REAL*8                                              :: eps_elem(Mesh%Nelems)

    CALL find_oscillations_elements(eps_elem, oscillations)

  ENDSUBROUTINE check_oscillations

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

ENDMODULE adaptivity_indicator_module
