#if !defined(NEUTRALGAMMA) || !defined(TEMPERATURE)
#error This check requires a temperature model with NEUTRALGAMMA.
#endif
PROGRAM test_neutral_energy_source
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite
  USE globals, ONLY: phys
  USE physics, ONLY: neutral_rt, compute_fEiN, compute_dfEiN_dU
  IMPLICIT NONE
  REAL*8 :: U(6), zero(6), df(6), Sn(6,6), Sn0(6), densities(3)
  REAL*8 :: f, kiz, kcx, dkiz(6), dkcx(6), realized, expected, tolerance
  INTEGER :: density_case, rate_case, velocity_case, checks

  phys%idx_rhon_eq = 5
  phys%idx_gamman_eq = 6
  phys%ionization_ion_energy_fraction = 0.d0
  neutral_rt%recombination_energy = 0.d0
  zero = 0.d0
  densities = [0.5d-7, 1.d-7, 2.d-7]
  checks = 0
  DO density_case = 1,SIZE(densities)
    DO rate_case = 1,2
      DO velocity_case = 0,1
        U = [2.d0, 0.d0, 3.d0, 4.d0, densities(density_case), &
             DBLE(velocity_case)*3.d-8]
        CALL compute_fEiN(U,f)
        CALL compute_dfEiN_dU(U,df)
        kiz = 2.d0
        kcx = 3.d0
        dkiz = 0.d0
        dkcx = 0.d0
        IF (rate_case == 2) THEN
          ! Synthetic density-dependent ionization and temperature-ratio CX.
          kiz = 2.d0*U(1)**0.2d0
          dkiz(1) = 0.2d0*kiz/U(1)
          kcx = 3.d0*(U(3)/U(1))**0.1d0
          dkcx(1) = -0.1d0*kcx/U(1)
          dkcx(3) = 0.1d0*kcx/U(3)
        ENDIF
        ! Isolate neutral bulk kinetic-energy transfer in the real assembly.
        CALL assemblyNeutral(U=U,niz=0.d0,dniz_dU=zero,nrec=0.d0,dnrec_dU=zero, &
             sigmaviz=kiz,dsigmaviz_dU=dkiz,sigmavrec=0.d0,dsigmavrec_dU=zero, &
             fGammacx=0.d0,dfGammacx_dU=zero,fGammarec=0.d0,dfGammarec_dU=zero, &
             fGammaN=0.d0,dfGammaN_dU=zero,sigmavcx=kcx,dsigmavcx_dU=dkcx, &
             fEiiz=0.d0,dfEiiz_dU=zero,fEirec=0.d0,dfEirec_dU=zero, &
             fEicx=0.d0,dfEicx_dU=zero,fEiN=f,dfEiN_dU=df,Sn=Sn,Sn0=Sn0, &
             sigmavEiz=0.d0,dsigmavEiz_dU=zero,sigmavErec=0.d0,dsigmavErec_dU=zero)
        realized = -DOT_PRODUCT(Sn(3,:),U)-Sn0(3)
        expected = 0.5d0*U(1)*U(6)**2/U(5)*(kiz+kcx)
        tolerance = 128.d0*EPSILON(1.d0)*ABS(expected)
        IF (.NOT. ieee_is_finite(realized) .OR. ABS(realized-expected) > tolerance) THEN
          WRITE(*,*) 'neutral density, rate case, velocity case:',U(5),rate_case,velocity_case
          WRITE(*,*) 'assembled source, physical source:',realized,expected
          ERROR STOP 'neutral kinetic-energy source reconstruction'
        ENDIF
        checks = checks+1
      ENDDO
    ENDDO
  ENDDO
  WRITE(*,*) 'Neutral kinetic-energy source: PASS; checks:',checks

CONTAINS

#include "neutral_source_assembly_routine.F"

END PROGRAM test_neutral_energy_source
