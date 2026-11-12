#ifdef TEMPERATURE
! Parallel conductive heat fluxes: conductivity, regional corrections and saturation.
SUBMODULE (physics) heat_fluxes
  USE flux_surface_geometry, ONLY: core_surface_geometry
  USE magnetic_topology, ONLY: magnetic_region_core
  IMPLICIT NONE
  REAL*8, PARAMETER :: energy_floor = 1.d-20
CONTAINS

  ! init captures one local assembly state and its geometry.
  ! Flow: classical conductivity -> closed-core correction -> conductivity ceiling
  !       -> conductive flux -> open-field free-streaming saturation.
  MODULE SUBROUTINE init(this,U,topology_region,rho,Q,b)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this
    REAL*8, INTENT(IN) :: U(:)
    INTEGER, OPTIONAL, INTENT(IN) :: topology_region
    REAL*8, OPTIONAL, INTENT(IN) :: rho,Q(:,:),b(:)
    LOGICAL :: closed_surface

    this%state = U(:4)
    this%ion_heat_flux = conductive_heat_flux_result()
    this%electron_heat_flux = conductive_heat_flux_result()
    this%gradient = 0.d0
    this%field_direction = 0.d0
    this%has_flux_inputs = PRESENT(Q) .AND. PRESENT(b)
    IF (this%has_flux_inputs) THEN
      this%gradient(:SIZE(Q,1),:) = Q(:,:4)
      this%field_direction(:SIZE(Q,1)) = b(:SIZE(Q,1))
    ENDIF
    this%length_m = 0.d0
    this%epsilon = 0.d0
    closed_surface = .FALSE.
    IF (PRESENT(topology_region)) closed_surface = topology_region == magnetic_region_core
    this%use_neoclassical_correction = phys%neoclassical_heat_conduction .AND. closed_surface
    this%use_free_streaming_saturation = switch%flux_limiter .AND. .NOT.closed_surface
    IF (phys%neoclassical_heat_conduction) THEN
      IF (.NOT.PRESENT(topology_region)) ERROR STOP 'Neoclassical conductivity requires magnetic topology'
      IF (closed_surface) THEN
        IF (.NOT.PRESENT(rho)) ERROR STOP 'Neoclassical conductivity requires rho'
        CALL core_surface_geometry%interp(rho,this%length_m,this%epsilon)
      ENDIF
    ENDIF
  END SUBROUTINE init

  MODULE SUBROUTINE evaluate_fluxes(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this

    IF (.NOT.this%has_flux_inputs) ERROR STOP 'Flux evaluation requires Q and b in init'
    CALL this%evaluate_conductivities()
    CALL this%compute_conductive_fluxes()
    IF (this%use_free_streaming_saturation) CALL this%apply_free_streaming_saturation()
  END SUBROUTINE evaluate_fluxes

  MODULE SUBROUTINE evaluate_conductivities(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this

    this%ion_heat_flux = conductive_heat_flux_result()
    this%electron_heat_flux = conductive_heat_flux_result()
    CALL this%compute_classical_conductivities()
    IF (this%use_neoclassical_correction) CALL this%apply_neoclassical_correction()
    CALL this%apply_conductivity_ceiling()
  END SUBROUTINE evaluate_conductivities

  MODULE SUBROUTINE apply_conductivity_ceiling(this)
    CLASS(parallel_heat_conduction_flux), INTENT(INOUT) :: this

    IF (switch%testcase /= 2) THEN
      CALL cap_species_conductivity(this%energy(1),phys%T_fluxlim_maxi, &
           this%use_neoclassical_correction,this%ion_heat_flux%conductivity)
      CALL cap_species_conductivity(this%energy(2),phys%T_fluxlim_maxe, &
           this%use_neoclassical_correction,this%electron_heat_flux%conductivity)
    ENDIF
  END SUBROUTINE apply_conductivity_ceiling

  SUBROUTINE cap_species_conductivity(energy,temperature_cap,corrected,conductivity)
    REAL*8, INTENT(IN) :: energy,temperature_cap
    LOGICAL, INTENT(IN) :: corrected
    TYPE(heat_conductivity), INTENT(INOUT) :: conductivity
    REAL*8 :: ceiling
    LOGICAL :: capped

    ceiling = (1.5d0*phys%Mref*temperature_cap)**phys%epn
    IF (corrected) THEN
      capped = conductivity%value >= ceiling
    ELSE
      capped = 2.d0/(3.d0*phys%Mref)*energy > temperature_cap
    ENDIF
    IF (capped) THEN
      conductivity%value = ceiling
      conductivity%state_derivative = 0.d0
      conductivity%density_degree = 0.d0
    ENDIF
  END SUBROUTINE cap_species_conductivity
END SUBMODULE heat_fluxes
#endif
