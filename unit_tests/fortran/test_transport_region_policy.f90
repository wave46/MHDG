PROGRAM test_transport_region_policy
  USE magnetic_topology, ONLY: magnetic_region_undefined, magnetic_region_core, &
       magnetic_region_main_sol, magnetic_region_private_flux
  USE transport_models_1d_config
  IMPLICIT NONE

  INTEGER :: policy
  LOGICAL :: valid
  REAL*8 :: velocity(2)

  CALL tm1d_region_policy_from_name('legacy_all_regions', policy, valid)
  CALL assert_true(valid, 'legacy policy parsing')
  CALL assert_true(policy == transport_region_legacy_all_regions, &
       'legacy policy identifier')
  CALL assert_true(TRIM(tm1d_region_policy_name(policy)) == &
       'legacy_all_regions', 'legacy policy name round trip')
  CALL assert_true(tm1d_region_is_included(policy, magnetic_region_undefined), &
       'legacy undefined inclusion')
  CALL assert_true(tm1d_region_is_included(policy, magnetic_region_private_flux), &
       'legacy private-flux inclusion')

  CALL tm1d_region_policy_from_name('CORE_AND_MAIN_SOL', policy, valid)
  CALL assert_true(valid, 'case-insensitive core-and-SOL parsing')
  CALL assert_true(TRIM(tm1d_region_policy_name(policy)) == &
       'core_and_main_sol', 'core-and-SOL name round trip')
  CALL assert_true(tm1d_region_is_included(policy, magnetic_region_core), &
       'core-and-SOL core inclusion')
  CALL assert_true(tm1d_region_is_included(policy, magnetic_region_main_sol), &
       'core-and-SOL main-SOL inclusion')
  CALL assert_true(.NOT. tm1d_region_is_included(policy, magnetic_region_private_flux), &
       'core-and-SOL private-flux exclusion')
  CALL assert_true(.NOT. tm1d_region_is_included(policy, magnetic_region_undefined), &
       'core-and-SOL undefined-region exclusion')

  CALL tm1d_region_policy_from_name('core_only', policy, valid)
  CALL assert_true(valid, 'core-only parsing')
  CALL assert_true(TRIM(tm1d_region_policy_name(policy)) == &
       'core_only', 'core-only name round trip')
  CALL assert_true(tm1d_region_is_included(policy, magnetic_region_core), &
       'core-only core inclusion')
  CALL assert_true(.NOT. tm1d_region_is_included(policy, magnetic_region_main_sol), &
       'core-only main-SOL exclusion')
  CALL assert_true(.NOT. tm1d_region_is_included(policy, magnetic_region_private_flux), &
       'core-only private-flux exclusion')

  CALL tm1d_region_policy_from_name('not_a_policy', policy, valid)
  CALL assert_true(.NOT. valid, 'invalid policy rejection')
  CALL assert_true(TRIM(tm1d_region_policy_name(-1)) == 'unknown', &
       'invalid policy formatting')

  CALL tm1d_region_policy_from_name('core_and_main_sol', policy, valid)
  CALL tm1d_build_pinch_velocity(2.d0, policy, [0.d0, 1.d0], &
       [1.d0, 0.d0], velocity)
  CALL assert_close(velocity(1), 2.d0, 'positive pinch is outward')
  CALL assert_close(velocity(2), 0.d0, 'positive pinch has no tangential component')

  CALL tm1d_build_pinch_velocity(-2.d0, policy, [0.d0, 1.d0], &
       [1.d0, 0.d0], velocity)
  CALL assert_close(velocity(1), -2.d0, 'negative pinch is inward')

  CALL tm1d_build_pinch_velocity(2.d0, policy, [0.d0, 1.d0], &
       [0.d0, 0.d0], velocity)
  CALL assert_close(NORM2(velocity), 0.d0, 'null topology normal suppresses pinch')

  CALL tm1d_build_pinch_velocity(2.d0, transport_region_legacy_all_regions, &
       [0.d0, 1.d0], [0.d0, 0.d0], velocity)
  CALL assert_close(velocity(1), 2.d0, 'legacy magnetic-field fallback')
  WRITE (*, '(A)') 'transport region policy checks: PASS'

CONTAINS

  SUBROUTINE assert_true(condition, label)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (.NOT. condition) THEN
       WRITE (*, '(2A)') 'FAIL: ', TRIM(label)
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_true

  SUBROUTINE assert_close(value, expected, label)
    REAL*8, INTENT(IN) :: value, expected
    CHARACTER(LEN=*), INTENT(IN) :: label

    IF (ABS(value - expected) > 1.d-12) THEN
       WRITE (*, '(A,1X,A,2(1X,ES12.4))') 'FAIL:', TRIM(label), value, expected
       ERROR STOP 1
    ENDIF
  END SUBROUTINE assert_close

END PROGRAM test_transport_region_policy
