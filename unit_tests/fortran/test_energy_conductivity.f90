PROGRAM test_energy_conductivity
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys, switch
  USE physics, ONLY: computeAlphai, computeAlphae, compute_dAlpha_dUi, compute_dAlpha_dUe
  IMPLICIT NONE
  REAL*8, PARAMETER :: energy_floor=1.d-20, perturbation=1.d-6
  REAL*8 :: u(5), trial(5), derivative(5), a, expected, plus, minus, step, numerical
  REAL*8 :: energies(8), cap, worst, error, scale
  INTEGER :: species, row, regime, j, checks

  phys%Mref=1.d0
  phys%epn=2.5d0
  switch%testcase=50
  ! Distinct species caps catch use of the wrong configured temperature.
  phys%T_fluxlim_maxi=2.d0
  phys%T_fluxlim_maxe=4.d0
  checks=0
  worst=0.d0

  DO species=1,2
    row=species+2
    j=0
    cap=1.5d0*phys%Mref*phys%T_fluxlim_maxi
    IF (species==2) cap=1.5d0*phys%Mref*phys%T_fluxlim_maxe
    energies=[-energy_floor,0.d0,.5d0*energy_floor,energy_floor, &
         2.d0*energy_floor,1.d0,cap,cap+1.d0]
    DO regime=1,SIZE(energies)
      u=[1.d0,0.d0,1.d0,1.d0,.2d0]
      u(row)=energies(regime)
      CALL coefficient(u,species,a,derivative)
      expected=MIN(MAX(energies(regime),energy_floor),cap)**phys%epn
      CALL close(a,expected,expected,2.d-14,'conductivity value')
      IF (regime<=4 .OR. regime==8) THEN
        IF (ANY(derivative/=0.d0)) ERROR STOP 'constant conductivity branch has nonzero derivative'
      ELSE
        expected=phys%epn*energies(regime)**(phys%epn-1.d0)
        CALL close(derivative(row),expected,expected,2.d-14,'active conductivity derivative')
      ENDIF
      CALL close(DOT_PRODUCT(derivative,u),0.d0,a,2.d-14,'degree-zero contraction')
      CALL coefficient(1.75d0*u,species,plus)
      CALL close(plus,a,a,2.d-14,'degree-zero value')

      ! One-sided checks select the documented constant/uncapped branches at
      ! the floor/cap joins. Other samples admit centered differences.
      step=perturbation*MAX(ABS(energies(regime)),energy_floor)
      trial=u
      IF (regime==4 .OR. regime==7) THEN
        trial(row)=u(row)-step
        CALL coefficient(trial,species,minus)
        numerical=(a-minus)/step
      ELSE
        trial(row)=u(row)+step
        CALL coefficient(trial,species,plus)
        trial(row)=u(row)-step
        CALL coefficient(trial,species,minus)
        numerical=(plus-minus)/(2.d0*step)
      ENDIF
      scale=MAX(a/energy_floor,1.d-100)
      IF (regime>=5) scale=MAX(ABS(derivative(row)),a/MAX(ABS(energies(regime)),energy_floor))
      CALL close(numerical,derivative(row),scale,3.d-6,'conductivity branch finite difference')
    ENDDO

    ! At nonzero flow, perturb every independent conservative component.
    ! The ion-energy case exercises the kinetic-energy subtraction.
    u=[1.7d0,.68d0,2.1d0,2.8d0,.2d0]
    CALL coefficient(u,species,a,derivative)
    DO j=1,SIZE(u)
      step=perturbation*MAX(ABS(u(j)),.2d0)
      trial=u
      trial(j)=u(j)+step
      CALL coefficient(trial,species,plus)
      trial(j)=u(j)-step
      CALL coefficient(trial,species,minus)
      CALL close((plus-minus)/(2.d0*step),derivative(j),a,3.d-7,'coupled conductivity finite difference')
    ENDDO
    CALL close(DOT_PRODUCT(derivative,u),0.d0,a,2.d-14,'flow-state degree-zero contraction')

    ! Testcase 2 intentionally retains the uncapped power law.
    switch%testcase=2
    u=[1.d0,0.d0,cap+1.d0,cap+1.d0,.2d0]
    CALL coefficient(u,species,a,derivative)
    CALL close(a,(cap+1.d0)**phys%epn,a,2.d-14,'convergence cap bypass value')
    CALL close(derivative(row),phys%epn*(cap+1.d0)**(phys%epn-1.d0),a,2.d-14,'convergence cap bypass derivative')
    switch%testcase=50
  ENDDO
  WRITE(*,*) 'Energy conductivity: PASS; checks, worst scaled error:',checks,worst

CONTAINS

  SUBROUTINE coefficient(state,which,value,tangent)
    REAL*8, INTENT(IN) :: state(:)
    INTEGER, INTENT(IN) :: which
    REAL*8, INTENT(OUT) :: value
    REAL*8, OPTIONAL, INTENT(OUT) :: tangent(:)
    IF (which==1) THEN
      value=computeAlphai(state)
      IF (PRESENT(tangent)) CALL compute_dAlpha_dUi(state,tangent)
    ELSE
      value=computeAlphae(state)
      IF (PRESENT(tangent)) CALL compute_dAlpha_dUe(state,tangent)
    ENDIF
  END SUBROUTINE coefficient

  SUBROUTINE close(actual,reference,normalization,tolerance,label)
    REAL*8, INTENT(IN) :: actual,reference,normalization,tolerance
    CHARACTER(*), INTENT(IN) :: label
    error=ABS(actual-reference)/MAX(ABS(actual),ABS(reference),normalization)
    IF (.NOT.ieee_is_finite(error) .OR. error>tolerance) THEN
      WRITE(*,*) label,': species, regime, component:',species,regime,j
      WRITE(*,*) 'actual, expected, scaled error:',actual,reference,error
      ERROR STOP 'conductivity consistency'
    ENDIF
    checks=checks+1
    worst=MAX(worst,error)
  END SUBROUTINE close

END PROGRAM test_energy_conductivity
