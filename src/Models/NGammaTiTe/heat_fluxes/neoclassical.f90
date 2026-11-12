#ifdef TEMPERATURE
! Hydrogenic-majority SOLPS-ITER neoclassical heat-conductivity fits.
! Source conventions and applicability: docs/models/heat-flux-limiter.md.
SUBMODULE (physics:heat_fluxes) neoclassical
  USE prec_const, ONLY: PI
  IMPLICIT NONE
  REAL*8, PARAMETER :: charge_si = 1.602176634d-19, epsilon0_si = 8.8541878128d-12
CONTAINS

  MODULE SUBROUTINE apply_neoclassical_correction(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this

    CALL correct_species_conductivity(this%state(1),this%energy(1),this%energy_derivative(:,1), &
         this%length_m,this%epsilon,3,this%ion_heat_flux%conductivity)
    CALL correct_species_conductivity(this%state(1),this%energy(2),this%energy_derivative(:,2), &
         this%length_m,this%epsilon,4,this%electron_heat_flux%conductivity)
  END SUBROUTINE apply_neoclassical_correction

  SUBROUTINE correct_species_conductivity(density,energy,energy_derivative,length_m,epsilon,row,conductivity)
    REAL*8, INTENT(IN) :: density,energy,energy_derivative(4),length_m,epsilon
    INTEGER, INTENT(IN) :: row
    TYPE(heat_conductivity), INTENT(INOUT) :: conductivity
    REAL*8 :: temperature_eV,factor,slope

    IF (energy > energy_floor) THEN
      temperature_eV = 2.d0/(3.d0*phys%Mref)*energy*simpar%refval_temperature
      CALL neoclassical_heat_factor(temperature_eV,density*simpar%refval_density, &
           length_m,epsilon,row,factor,slope)
      conductivity%value = conductivity%value*factor
      conductivity%state_derivative = conductivity%value/energy*(phys%epn-2.d0*slope)*energy_derivative
      conductivity%state_derivative(1) = conductivity%state_derivative(1)+conductivity%value/density*slope
      conductivity%density_degree = slope
    ENDIF
  END SUBROUTINE correct_species_conductivity

  PURE SUBROUTINE neoclassical_heat_factor(temperature_eV,density_m3,length_m,epsilon, &
       energy_row,factor,collisionality_slope)
    REAL*8, INTENT(IN) :: temperature_eV,density_m3,length_m,epsilon
    INTEGER, INTENT(IN) :: energy_row
    REAL*8, INTENT(OUT) :: factor,collisionality_slope
    REAL*8 :: collision_scale,x,amplitude,a,b,c,d,denominator,fit,dfit_dlogx,prefactor

    ! x = epsilon**(3/2)*nu_star = L/(2*pi*tau_s*sqrt(T_J/m_s)).
    ! The species mass cancels for a singly charged hydrogenic majority.
    collision_scale = (4.d0*PI*epsilon0_si/charge_si**2)**2*charge_si**2
    IF (energy_row == 4) THEN
      collision_scale = 3.d0/(4.d0*SQRT(2.d0*PI))*collision_scale
      amplitude=1.d0; a=2.01d0; b=1.53d0; c=.52d0; d=.89d0; prefactor=.58d0
    ELSE
      collision_scale = 3.d0/(4.d0*SQRT(PI))*collision_scale
      amplitude=.66d0+1.88d0*SQRT(epsilon)-1.54d0*epsilon
      a=1.03d0; b=.31d0; c=1.17d0; d=.74d0; prefactor=1.60d0
    ENDIF
    x = length_m*density_m3*coulomb_logarithm/(2.d0*PI*collision_scale*temperature_eV**2)
    ! Divide K by epsilon**(3/2) analytically; this also gives the axis limit.
    denominator = epsilon**1.5d0+a*epsilon**.75d0*SQRT(x)+b*x
    fit = amplitude/denominator+c*x/(1.d0+d*x)
    factor = prefactor/fit
    dfit_dlogx = -amplitude*(.5d0*a*epsilon**.75d0*SQRT(x)+b*x)/denominator**2 &
         +c*x/(1.d0+d*x)**2
    collisionality_slope = -dfit_dlogx/fit
  END SUBROUTINE neoclassical_heat_factor
END SUBMODULE neoclassical
#endif
