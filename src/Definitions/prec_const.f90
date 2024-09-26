!********************************************
! project: MHDG
! file: prec_const.f90
! date: 06/09/2016
! Definition of real precision and real constants
!********************************************

MODULE prec_const
  IMPLICIT NONE

  ! Precision for real
  !*****************
  INTEGER, PARAMETER :: float = SELECTED_REAL_KIND(13, 99)

  ! Some useful constants
  !*********************
  REAL(float), PARAMETER :: PI = 3.141592653589793238462643383279502884197_float
  REAL(float), PARAMETER :: sqrt2 = 1.414213562373095_float

END MODULE prec_const
