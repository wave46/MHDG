PROGRAM test_equilibrium_series
  USE prec_const, ONLY: float
  USE heating_series, ONLY: heating_series_t, heating_sample_t
  USE equilibrium_series, ONLY: equilibrium_series_t, equilibrium_sample_t
  IMPLICIT NONE
  TYPE(equilibrium_series_t) :: series, other
  TYPE(equilibrium_sample_t) :: sample
  REAL(float), ALLOCATABLE :: r(:), z(:), times(:), jr(:), jz(:)
  REAL(float) :: requests(10), t, x, y
  INTEGER :: ierr, i, j, k
  CHARACTER(LEN=1000) :: filename, message, current_filename
  CHARACTER(LEN=32) :: mode

  CALL get_command_argument(1, filename)
  CALL get_command_argument(2, mode)
  IF (INDEX(mode, 'heating') == 1) THEN
     CALL check_heating_input()
     STOP
  ENDIF
  CALL get_command_argument(3, current_filename)
  CALL open_input(series)
  CALL require(ierr == 0, message)
  CALL series%grid(r, z)
  CALL series%time_points(times)
  CALL series%current_grid(jr, jz)
  IF (mode == 'inspect') THEN
     CALL series%sample_at(times(1), sample, ierr, message)
     CALL require(ierr == 0, message)
     WRITE(*,'(A,3(1X,I0))') 'EQUILIBRIUM_INSPECT_PASS', SIZE(r), SIZE(z), SIZE(times)
     CALL series%clear()
     STOP
  ENDIF
  CALL require(SIZE(r) == 4 .AND. SIZE(z) == 3, 'non-square grid sizes')
  CALL require(MAXVAL(ABS(r-[2._float,2.1_float,2.3_float,2.6_float])) < 1.e-12_float, 'R coordinates')
  CALL require(MAXVAL(ABS(z-[-1._float,-0.1_float,0.5_float])) < 1.e-12_float, 'Z coordinates')
  requests = [0._float,0.03_float,0.08_float,0.12_float,0.56_float,0.69_float,0.8_float, &
       0.2_float,0.01_float,0.56_float]
  DO k = 1, SIZE(requests)
     CALL series%sample_at(requests(k), sample, ierr, message)
     CALL require(ierr == 0, message)
     t = requests(k)
     IF (INDEX(mode, 'scale') > 0) THEN
        CALL require(ALLOCATED(sample%poloidal_scale), 'declared scale absent from sample')
        CALL require(ABS(sample%poloidal_scale+1._float/(2*ACOS(-1._float))) < 1.e-14_float, &
             'constant declared scale through cache/interpolation')
        sample%poloidal_scale = 999._float
     ENDIF
     IF (mode == 'stationary') t = 0
     DO i = 1, SIZE(r)
        DO j = 1, SIZE(z)
           x = r(i)
           y = z(j)
           CALL require(ABS(sample%psi(j,i)-(2*x*x+3*y+7*t)) < 1.e-12_float, 'psi interpolation/orientation')
           CALL require(ABS(sample%bphi(j,i)-(5+x-y+11*t)) < 1.e-12_float, 'Bphi interpolation')
           IF (mode /= 'optional') THEN
              CALL require(ABS(sample%br(j,i)-(0.2_float+x+2*y+3*t)) < 1.e-12_float, 'Br interpolation')
              CALL require(ABS(sample%bz(j,i)-(-0.4_float+2*x-y-4*t)) < 1.e-12_float, 'Bz interpolation')
           ENDIF
        ENDDO
     ENDDO
     IF (ALLOCATED(jr)) THEN
        DO i = 1, SIZE(jr)
           DO j = 1, SIZE(jz)
              CALL require(ABS(sample%jtor(j,i)-(1+3*jr(i)+4*jz(j)+13*t)) < 1.e-12_float, &
                   'Jtor interpolation on its own grid')
           ENDDO
        ENDDO
     ENDIF
     CALL require(ABS(sample%psi_sep-(40+17*t)) < 1.e-12_float, 'psiSep interpolation')
     CALL require(ABS(sample%time-requests(k)) < 1.e-12_float, 'sample physical time')
     IF (mode == 'optional') CALL require(.NOT. ALLOCATED(sample%br) .AND. &
          .NOT. ALLOCATED(sample%bz) .AND. .NOT. ALLOCATED(sample%jtor), 'optional references')
  ENDDO
  IF (mode /= 'stationary') THEN
     CALL require(SIZE(times) == 5, 'trajectory timestamps')
     CALL series%sample_at(0.03_float, sample, ierr, message)
     CALL require(ierr == 0, message)
     CALL require(sample%frame_left == 0 .AND. sample%frame_right == 1 .AND. &
          ABS(sample%fraction-0.375_float) < 1.e-12_float, 'interpolation position')
     CALL series%sample_at(0.8_float+EPSILON(t), sample, ierr, message)
     CALL require(ierr == 0, message)
     CALL require(sample%frame_left == 4 .AND. sample%frame_right == 4, 'final endpoint position')
     CALL series%sample_at(0.9_float, sample, ierr, message)
     CALL require(ierr /= 0 .AND. INDEX(message,'outside trajectory coverage') > 0, 'coverage rejection')
  ELSE
     CALL series%sample_at(10._float, sample, ierr, message)
     CALL require(ierr == 0, message)
     CALL require(ABS(sample%psi_sep-40._float) < 1.e-12_float .AND. &
          sample%time == 10._float, 'fixed equilibrium independent of simulation time')
  ENDIF
  ! Closing/reopening one provider must leave another provider usable.
  CALL open_input(other)
  CALL require(ierr == 0, message)
  CALL series%clear()
  CALL series%sample_at(0._float, sample, ierr, message)
  CALL require(ierr /= 0, 'closed provider rejection')
  CALL other%sample_at(0.08_float, sample, ierr, message)
  CALL require(ierr == 0, message)
  CALL open_input(other)
  CALL require(ierr == 0, message)
  CALL other%sample_at(0._float, sample, ierr, message)
  CALL require(ierr == 0, message)
  CALL other%clear()
  WRITE(*,'(A)') 'EQUILIBRIUM_SERIES_PASS'
CONTAINS
  SUBROUTINE check_heating_input()
    TYPE(heating_series_t) :: heating
    TYPE(heating_sample_t) :: power
    REAL(float), ALLOCATABLE :: hr(:), hz(:), ht(:)
    REAL(float) :: targets(8), actual, expected
    INTEGER :: m, ir, iz

    CALL heating%open(filename, ierr, message)
    CALL require(ierr == 0, message)
    CALL heating%grid(hr, hz)
    CALL heating%time_points(ht)
    CALL require(SIZE(hr) == 3 .AND. SIZE(hz) == 2, 'independent heating grid dimensions')
    targets = [0._float,0.03_float,0.08_float,0.17_float,0.8_float,0.45_float,0.69_float,0.01_float]
    DO m = 1, SIZE(targets)
       CALL heating%sample_at(targets(m), power, ierr, message)
       CALL require(ierr == 0, message)
       actual = targets(m)
       IF (mode == 'heating_fixed') actual = 0
       DO ir = 1, SIZE(hr)
          DO iz = 1, SIZE(hz)
             CALL require(ABS(power%ions(iz,ir)-(20+2*hr(ir)+3*hz(iz)+7*actual)) < 1.e-12_float, &
                  'ion heating orientation/interpolation')
             CALL require(ABS(power%electrons(iz,ir)-(30+hr(ir)**2+hz(iz)-4*actual)) < 1.e-12_float, &
                  'electron heating orientation/interpolation')
          ENDDO
       ENDDO
       CALL require(ABS(power%time-targets(m)) < 1.e-12_float, 'heating physical time')
    ENDDO
    IF (mode == 'heating_fixed') THEN
       CALL heating%sample_at(10._float, power, ierr, message)
       CALL require(ierr == 0, message)
    ELSE
       CALL require(SIZE(ht) == 5 .AND. ABS(ht(2)-0.05_float) < 1.e-12_float, 'independent heating times')
       CALL heating%sample_at(0.8_float+EPSILON(actual), power, ierr, message)
       CALL require(ierr == 0, message)
       CALL require(power%frame_left == 4 .AND. power%frame_right == 4, 'heating final endpoint')
       CALL heating%sample_at(0.9_float, power, ierr, message)
       CALL require(ierr /= 0, 'heating coverage rejection')
    ENDIF
    CALL heating%sample_at(0._float, power, ierr, message)
    CALL require(ierr == 0, message)
    expected = power%ions(1,1)
    power%ions = -999
    CALL heating%sample_at(0._float, power, ierr, message)
    CALL require(ierr == 0, message)
    CALL require(power%ions(1,1) == expected, 'heating samples must not mutate cached primitives')
    CALL heating%clear()
    CALL heating%sample_at(0._float, power, ierr, message)
    CALL require(ierr /= 0, 'closed heating input')
    CALL heating%open(filename, ierr, message)
    CALL require(ierr == 0, message)
    CALL heating%sample_at(0._float, power, ierr, message)
    CALL require(ierr == 0, message)
    CALL heating%clear()
    WRITE(*,'(A)') 'EQUILIBRIUM_SERIES_PASS'
  END SUBROUTINE check_heating_input

  SUBROUTINE open_input(input)
    TYPE(equilibrium_series_t), INTENT(INOUT) :: input
    IF (INDEX(mode, 'numbered') == 1) THEN
       CALL input%open_numbered(filename, 4, ierr, message, current_filename)
    ELSE
       CALL input%open(filename, ierr, message, current_filename)
    ENDIF
  END SUBROUTINE open_input

  SUBROUTINE require(condition, detail)
    LOGICAL, INTENT(IN) :: condition
    CHARACTER(LEN=*), INTENT(IN) :: detail
    IF (.NOT. condition) THEN
       WRITE(*,'(A)') 'PROVIDER_ERROR: '//TRIM(detail)
       STOP 1
    ENDIF
  END SUBROUTINE require
END PROGRAM test_equilibrium_series
