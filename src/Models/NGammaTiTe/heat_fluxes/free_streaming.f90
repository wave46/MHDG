#ifdef TEMPERATURE
! Open-field free-streaming saturation with a lagged multiplier.
SUBMODULE (physics:heat_fluxes) free_streaming
  IMPLICIT NONE
  REAL*8, PARAMETER :: heat_flux_floor = 1.d-10, temperature_floor_eV = 1.d-5
CONTAINS

  MODULE SUBROUTINE apply_free_streaming_saturation(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this
    REAL*8 :: free_streaming_energy_floor,ion_free_streaming_flux,electron_free_streaming_flux

    free_streaming_energy_floor = temperature_floor_eV/simpar%refval_temperature*3*phys%Mref/2.
    ion_free_streaming_flux = (2./3.)**1.5*MAX(this%energy(1),free_streaming_energy_floor)**1.5*this%state(1)
    electron_free_streaming_flux = (2./3.)**1.5*SQRT(ion_mass_si/electron_mass_si) &
         *MAX(this%energy(2),free_streaming_energy_floor)**1.5*this%state(1)
    CALL saturate_species_flux(ion_free_streaming_flux,phys%c_fli,this%ion_heat_flux)
    CALL saturate_species_flux(electron_free_streaming_flux,phys%c_fle,this%electron_heat_flux)
  END SUBROUTINE apply_free_streaming_saturation

  SUBROUTINE saturate_species_flux(free_streaming_flux,fraction,heat_flux)
    REAL*8, INTENT(IN) :: free_streaming_flux,fraction
    TYPE(conductive_heat_flux_result), INTENT(INOUT) :: heat_flux
    REAL*8 :: capped_free_streaming_flux

    capped_free_streaming_flux = MAX(ABS(free_streaming_flux),heat_flux_floor)
    heat_flux%limiter = 1.d0/(1.d0+ABS(heat_flux%diffusive_flux)/(fraction*capped_free_streaming_flux))
    heat_flux%diffusive_flux = heat_flux%limiter*heat_flux%diffusive_flux
    heat_flux%state_derivative = heat_flux%limiter*heat_flux%state_derivative
    heat_flux%gradient_derivative = heat_flux%limiter*heat_flux%gradient_derivative
    heat_flux%affine = heat_flux%limiter*heat_flux%affine
  END SUBROUTINE saturate_species_flux
END SUBMODULE free_streaming
#endif
