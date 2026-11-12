#ifdef TEMPERATURE
! Classical parallel conductivity and conductive flux derivatives.
SUBMODULE (physics:heat_fluxes) classical
  IMPLICIT NONE
CONTAINS

  MODULE SUBROUTINE compute_classical_conductivities(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this

    this%energy(1) = this%state(3)/this%state(1)-0.5d0*(this%state(2)/this%state(1))**2
    this%energy(2) = this%state(4)/this%state(1)
    CALL computeVi(this%state,this%energy_derivative(:,1))
    CALL computeVe(this%state,this%energy_derivative(:,2))
    CALL compute_species_conductivity(this%energy(1),this%energy_derivative(:,1),this%ion_heat_flux%conductivity)
    CALL compute_species_conductivity(this%energy(2),this%energy_derivative(:,2),this%electron_heat_flux%conductivity)
  END SUBROUTINE compute_classical_conductivities

  MODULE SUBROUTINE compute_conductive_fluxes(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this
    REAL*8 :: parallel_state_gradient(4),ion_hessian(4,4),electron_hessian(4,4)

    parallel_state_gradient = MATMUL(this%field_direction,this%gradient)
    CALL compute_dV_dUi(this%state,ion_hessian)
    CALL compute_dV_dUe(this%state,electron_hessian)
    CALL compute_species_flux(this%energy_derivative(:,1),ion_hessian,parallel_state_gradient, &
         this%field_direction,phys%diff_pari,this%ion_heat_flux)
    CALL compute_species_flux(this%energy_derivative(:,2),electron_hessian,parallel_state_gradient, &
         this%field_direction,phys%diff_pare,this%electron_heat_flux)
  END SUBROUTINE compute_conductive_fluxes

  SUBROUTINE compute_species_conductivity(energy,energy_derivative,conductivity)
    REAL*8, INTENT(IN) :: energy,energy_derivative(4)
    TYPE(heat_conductivity), INTENT(OUT) :: conductivity

    conductivity%value = MAX(energy,energy_floor)**phys%epn
    conductivity%state_derivative = 0.d0
    conductivity%density_degree = 0.d0
    IF (energy > energy_floor) &
         conductivity%state_derivative = phys%epn*energy**(phys%epn-1.d0)*energy_derivative
  END SUBROUTINE compute_species_conductivity

  SUBROUTINE compute_species_flux(energy_derivative,energy_hessian,parallel_state_gradient,b,parallel_diffusion,heat_flux)
    REAL*8, INTENT(IN) :: energy_derivative(4),energy_hessian(4,4),parallel_state_gradient(4),b(3),parallel_diffusion
    TYPE(conductive_heat_flux_result), INTENT(INOUT) :: heat_flux
    REAL*8 :: parallel_energy_gradient,coefficient,conductivity_factor
    INTEGER :: k

    coefficient = parallel_diffusion*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
    conductivity_factor = heat_flux%conductivity%value
    parallel_energy_gradient = DOT_PRODUCT(parallel_state_gradient,energy_derivative)
    heat_flux%diffusive_flux = coefficient*conductivity_factor*parallel_energy_gradient
    heat_flux%state_derivative = coefficient*(parallel_energy_gradient*heat_flux%conductivity%state_derivative &
         +conductivity_factor*MATMUL(parallel_state_gradient,energy_hessian))
    DO k=1,3
      heat_flux%gradient_derivative(k,:) = coefficient*conductivity_factor*b(k)*energy_derivative
    ENDDO
    ! Energy has degree zero in U; its state derivative has degree minus one.
    heat_flux%affine = heat_flux%diffusive_flux*(heat_flux%conductivity%density_degree-1.d0)
  END SUBROUTINE compute_species_flux
END SUBMODULE classical
#endif
