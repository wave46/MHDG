#ifdef TEMPERATURE
! Model coefficients only; HDG owns jumps, matrix placement and boundary laws.
SUBMODULE (physics) physics_stabilization
  IMPLICIT NONE
CONTAINS

  MODULE PROCEDURE evaluate_stab5
    REAL*8 :: advective,bn,bnorm,Etan
    INTEGER :: i

    bn = DOT_PRODUCT(b(1:SIZE(normal)),normal)
    bnorm = NORM2(b(1:SIZE(normal)))
    advective = MAX(ABS(5.d0/3.d0*up(2)*bn), &
         ABS(0.3d0*bn*(3.d0*uc(1)+ &
         SQRT(ABS(10.d0*uc(3)*uc(1)+10.d0*uc(4)*uc(1)-5.d0*uc(2)**2)))/uc(1)))
    tau = 0.d0
    DO i = 1,SIZE(uc)
      tau(i,i) = advective+diff_iso(i,i)*length_factor
    ENDDO
    tau(3,3) = tau(3,3)+thermal_stabilization(up(7),phys%T_fluxlim_maxi,phys%diff_pari)
    tau(4,4) = tau(4,4)+thermal_stabilization(up(8),phys%T_fluxlim_maxe,phys%diff_pare)
#ifdef NEUTRALGAMMA
    IF (phys%idx_gamman_eq > 0 .AND. phys%idx_rhon_eq > 0) THEN
      CALL computeEtan(uc,Etan)
      tau(phys%idx_gamman_eq,phys%idx_gamman_eq) = &
           advective+Etan/uc(phys%idx_rhon_eq)*length_factor
    ENDIF
#endif

  CONTAINS

    FUNCTION thermal_stabilization(temperature,temperature_cap,parallel_diffusion) RESULT(value)
      REAL*8, INTENT(IN) :: temperature,temperature_cap,parallel_diffusion
      REAL*8 :: value
      value = ABS(bn)*parallel_diffusion*MIN(temperature_cap,temperature)**2.5d0* &
           bnorm/uc(1)*length_factor
    END FUNCTION thermal_stabilization

  END PROCEDURE evaluate_stab5

END SUBMODULE physics_stabilization
#endif
