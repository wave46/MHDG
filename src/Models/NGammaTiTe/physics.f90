!*****************************************
! project: MHDG
! file: physics.f90
! date: 20/09/2017
! Define the physics of the model
!  ******** N-Gamma-Ti-Te system     ****
!*****************************************
MODULE physics
  USE globals
  USE magnetic_field
  IMPLICIT NONE

  REAL*8, PARAMETER :: eirene_rate_te_min_phys = 1.d-1
  REAL*8, PARAMETER :: eirene_rate_te_max_phys = 2.d4
  REAL*8, PARAMETER :: eirene_rate_ti_min_phys = 1.d-1
  REAL*8, PARAMETER :: eirene_rate_ti_max_phys = 3.d3
  REAL*8, PARAMETER :: eirene_rate_ne_min_phys = 1.d14
  REAL*8, PARAMETER :: eirene_rate_ne_max_phys = 1.d22
  REAL*8, PARAMETER :: neutral_state_tol_default = 1.d-20
  REAL*8, PARAMETER :: neutral_te_floor_phys = 1.d-10
  REAL*8, PARAMETER :: neutral_ti_floor_phys = 1.d-10
  REAL*8, PARAMETER :: neutral_ne_floor = 1.d-20
  REAL*8, PARAMETER :: neutral_sigmavnn_prefactor_phys = 5.2958d-11*1.d-6
  REAL*8, PARAMETER :: neutral_rydberg_energy_phys = 13.6d0
  REAL*8, PARAMETER :: neutral_iz_te_floor_phys = 0.05d0
  REAL*8, PARAMETER :: neutral_rec_te_floor_phys = 0.1d0
  REAL*8, PARAMETER :: neutral_cx_te_floor_phys = 0.05d0
  REAL*8, PARAMETER :: neutral_tloss_offset_phys = 25.d0
  REAL*8, PARAMETER :: neutral_tloss_amplitude_phys = 170.d0
  REAL*8, PARAMETER :: neutral_tlossrec_prefactor_phys = 8.d0
  REAL*8, PARAMETER :: neutral_tlossrec_cap_phys = 250.d0
  REAL*8, PARAMETER :: neutral_legacy_cx_prefactor_phys = 2.5d-15/EXP(-0.5d0)
  REAL*8, PARAMETER :: neutral_manuelcx_coeffs(5) = (/-6.837d-4, 8.004d-3, -3.581d-2, 4.518d-1, -3.259d1/)
  REAL*8, PARAMETER :: boltzmann_constant_si = 1.38064852d-23
  REAL*8, PARAMETER :: elementary_charge_si = 1.60217662d-19

  TYPE neutral_runtime_constants_t
    REAL*8 :: rate_scale
    REAL*8 :: energy_weight
    REAL*8 :: eirene_te_min
    REAL*8 :: eirene_te_max
    REAL*8 :: eirene_ti_min
    REAL*8 :: eirene_ti_max
    REAL*8 :: eirene_ne_min
    REAL*8 :: eirene_ne_max
    REAL*8 :: state_tol
    REAL*8 :: te_floor
    REAL*8 :: ti_floor
    REAL*8 :: sigmavnn_prefactor
    REAL*8 :: transport_ti_floor
    REAL*8 :: transport_ti_supp
    REAL*8 :: rydberg_energy
    REAL*8 :: iz_te_floor
    REAL*8 :: rec_te_floor
    REAL*8 :: cx_te_floor
    REAL*8 :: log_temperature_ref
    REAL*8 :: tloss_offset
    REAL*8 :: tloss_amplitude
    REAL*8 :: tloss_decay
    REAL*8 :: tlossrec_prefactor
    REAL*8 :: tlossrec_cap
    REAL*8 :: tlossrec_growth
    REAL*8 :: tlossrec_cap_log
    REAL*8 :: legacy_cx_prefactor
    REAL*8 :: recombination_energy
  END TYPE neutral_runtime_constants_t

  TYPE(neutral_runtime_constants_t), SAVE :: neutral_rt

CONTAINS

  !*******************************************
  ! Set model layout bookkeeping
  !*******************************************
  SUBROUTINE set_model_layout()

    phys%idx_rhon_eq = 0
    phys%idx_gamman_eq = 0
    phys%idx_k_eq = 0
    phys%idx_rhon_pv = 0
    phys%idx_un_pv = 0
    phys%idx_k_pv = 0

    phys%Neq = 4
    phys%npv = 10
    simpar%model = 'N-Gamma-Ti-Te'

#ifdef NEUTRAL
    phys%idx_rhon_eq = 5
    phys%idx_rhon_pv = 11
    phys%Neq = phys%idx_rhon_eq
    phys%npv = phys%idx_rhon_pv
    simpar%model = 'N-Gamma-Ti-Te-Neutral'
#endif

#ifdef NEUTRALGAMMA
    if (phys%idx_rhon_eq == 0) then
      phys%idx_rhon_eq = 5
      phys%idx_rhon_pv = 11
      phys%Neq = phys%idx_rhon_eq
      phys%npv = phys%idx_rhon_pv
    end if
    phys%idx_gamman_eq = phys%idx_rhon_eq + 1
    phys%idx_un_pv = phys%idx_rhon_pv + 1
    phys%Neq = phys%idx_gamman_eq
    phys%npv = phys%idx_un_pv
    simpar%model = 'N-Gamma-Ti-Te-NeutralGamma'
#endif

#ifdef KEQUATION
    phys%idx_k_eq = phys%Neq + 1
    phys%idx_k_pv = phys%npv + 1
    phys%Neq = phys%idx_k_eq
    phys%npv = phys%idx_k_pv
#ifdef NEUTRALGAMMA
    simpar%model = 'N-Gamma-Ti-Te-NeutralGamma-k'
#else
#ifdef NEUTRAL
    simpar%model = 'N-Gamma-Ti-Te-Neutral-k'
#endif
#endif
#endif

  END SUBROUTINE set_model_layout

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE initPhys()

    call set_model_layout()

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

    phys%phyVarNam = ""
    phys%conVarNam = ""

    ! Set the name of the physical variables
    phys%phyVarNam(1) = "rho" ! density
    phys%phyVarNam(2) = "u"   ! parallel velocity
    phys%phyVarNam(3) = "Ei"  ! total energy of ions
    phys%phyVarNam(4) = "Ee"  ! total energy of electrons
    phys%phyVarNam(5) = "pi"  ! pressure of ions
    phys%phyVarNam(6) = "pe"  ! pressure of electrons
    phys%phyVarNam(7) = "Ti"  ! temperature of ions
    phys%phyVarNam(8) = "Te"  ! temperature of electrons
    phys%phyVarNam(9) = "Csi" ! sound speed
    phys%phyVarNam(10)= "M"   ! Mach
    if (phys%idx_rhon_pv > 0) phys%phyVarNam(phys%idx_rhon_pv)= "rhon" ! density neutral
    if (phys%idx_un_pv > 0) phys%phyVarNam(phys%idx_un_pv)= "un" ! neutral parallel velocity
    if (phys%idx_k_pv > 0) phys%phyVarNam(phys%idx_k_pv)= "k" ! turbulent energy

    ! Set the name of the conservative variables
    phys%conVarNam(1) = "rho"   ! U1 = rho
    phys%conVarNam(2) = "Gamma" ! U2 = rho*u
    phys%conVarNam(3) = "nEi"   ! U3 = rho*Ei
    phys%conVarNam(4) = "nEe"   ! U4 = rho*Ee
    if (phys%idx_rhon_eq > 0) phys%conVarNam(phys%idx_rhon_eq) = "rhon" ! neutral density
    if (phys%idx_gamman_eq > 0) phys%conVarNam(phys%idx_gamman_eq) = "Gamman" ! neutral momentum
    if (phys%idx_k_eq > 0) phys%conVarNam(phys%idx_k_eq) = "k" ! turbulent energy

    simpar%Ndim = 2
#ifdef TOR3D
    simpar%Ndim = 3
#endif
    simpar%Neq = phys%Neq

    ALLOCATE (simpar%physvar_refval(phys%npv))
    ALLOCATE (simpar%consvar_refval(phys%Neq))
    simpar%physvar_refval = 0.
    simpar%consvar_refval = 0.

    simpar%physvar_refval(1) = simpar%refval_density
    simpar%physvar_refval(2) = simpar%refval_speed
    simpar%physvar_refval(3) = simpar%refval_specenergy
    simpar%physvar_refval(4) = simpar%refval_specenergy
    simpar%physvar_refval(5) = simpar%refval_specpress
    simpar%physvar_refval(6) = simpar%refval_specpress
    simpar%physvar_refval(7) = simpar%refval_temperature
    simpar%physvar_refval(8) = simpar%refval_temperature
    simpar%physvar_refval(9) = simpar%refval_speed
    simpar%physvar_refval(10) = 1.
    if (phys%idx_rhon_pv > 0) simpar%physvar_refval(phys%idx_rhon_pv) = simpar%refval_neutral
    if (phys%idx_un_pv > 0) simpar%physvar_refval(phys%idx_un_pv) = simpar%refval_speed
#ifdef KEQUATION
    if (phys%idx_k_pv > 0) simpar%physvar_refval(phys%idx_k_pv) = simpar%refval_k
#endif
    simpar%consvar_refval(1) = simpar%refval_density
    simpar%consvar_refval(2) = simpar%refval_momentum
    simpar%consvar_refval(3) = simpar%refval_specenergydens
    simpar%consvar_refval(4) = simpar%refval_specenergydens
    if (phys%idx_rhon_eq > 0) simpar%consvar_refval(phys%idx_rhon_eq) = simpar%refval_neutral
    if (phys%idx_gamman_eq > 0) simpar%consvar_refval(phys%idx_gamman_eq) = simpar%refval_momentum
#ifdef KEQUATION
    if (phys%idx_k_eq > 0) simpar%consvar_refval(phys%idx_k_eq) = simpar%refval_k
#endif
#ifdef EXPANDEDCX
#ifdef AMJUELCX
    ! coefficients for AMJUEL spline 3.1.8 FJ
    phys%alpha_cx = (/-1.841756e+01,  5.282950e-01, -2.200477e-01,  9.750192e-02,&
    -1.749183e-02,  4.954298e-04,  2.174910e-04, -2.530206e-05,&
    8.230751e-07/)
#endif
#ifdef THERMALCX
    phys%alpha_cx = (/-1.87744894e+01,  4.51800000e-01, -3.58100000e-02,  8.00400000e-03, -6.83700000e-04/)
#endif
#endif
#ifdef AMJUELSPLINES

    ! coefficients for AMJUEL 2.1.5JH
    phys%alpha_iz(:,1) = (/-3.29264710e+01,  1.42397767e+01, -6.51943873e+00,&
                            2.00999615e+00, -4.28959442e-01,  6.04783461e-02,&
                           -5.30473797e-03,  2.60694695e-04, -5.46790307e-06/)
    phys%alpha_iz(:,2) = (/1.29348138e-02, -1.17314396e-02, -7.18982575e-03,&
                           1.27597974e-02, -5.34086632e-03,  9.62490059e-04,&
                          -7.85487245e-05,  2.31744225e-06,  6.07738004e-09/)
    phys%alpha_iz(:,3) = (/5.51756251e-03,  1.06344011e-03,  9.24737741e-04,&
                          -4.69347962e-03,  2.32458236e-03, -4.18298118e-04,&
                           2.73582380e-05,  5.14889078e-08, -4.71289307e-08/)
    phys%alpha_iz(:,4) = (/-7.85381632e-04, -1.60005353e-03,  2.03702675e-03,&
                           -2.38922414e-05, -3.21722808e-04,  7.95723018e-05,&
                           -5.91534856e-06, -7.14418252e-09,  1.08685876e-08/)
    phys%alpha_iz(:,5) = (/1.43612850e-04,  1.13655464e-05, -3.66871720e-04,&
                           1.35806992e-04,  6.66058141e-06, -7.44704256e-06,&
                           8.66630287e-07, -2.54019475e-08, -3.44841725e-10/)
    phys%alpha_iz(:,6) = (/-3.88375028e-07,  5.17766228e-05,  5.36863032e-06,&
                           -1.45489756e-05,  2.39653187e-06,  1.84915526e-07,&
                           -6.11551482e-08,  4.09785784e-09, -8.71418322e-11/)
    phys%alpha_iz(:,7) = (/-1.48977436e-06, -7.94799990e-06,  3.71395891e-06,&
                            4.21203150e-08, -1.78520832e-07,  1.61823364e-08,&
                            1.07547317e-09, -2.04865734e-10,  8.02366070e-12/)
    phys%alpha_iz(:,8) = (/1.41636143e-07,  4.50850568e-07, -3.12576437e-07,&
                           5.50604467e-08,  6.09564957e-10, -8.18292830e-10,&
                           4.11800067e-11,  3.02791637e-12, -2.39651850e-13/)
    phys%alpha_iz(:,9) = (/-3.89093208e-09, -8.95261409e-09,  7.45121322e-09,&
                           -1.85267764e-09,  1.47020423e-10,  4.83578962e-12,&
                           -1.08932309e-12,  1.15585402e-14,  2.17364528e-15/)
    ! coefficients for AMJUEL 2.1.5JH
    phys%alpha_energy_iz(:,1) = (/-2.50812402e+01,  9.96163441e+00, -4.77618017e+00,&
    1.63071304e+00, -3.86224646e-01,  5.90834812e-02,&
   -5.50214904e-03,  2.82569314e-04, -6.12637364e-06/)
    phys%alpha_energy_iz(:,2) = (/1.73410814e-02, -1.57330788e-02,  2.97091760e-04,&
    3.45781992e-03, -1.35470702e-03,  2.46767178e-04,&
   -2.55003960e-05,  1.47937419e-06, -3.76891493e-08/)
    phys%alpha_energy_iz(:,3) = (/-1.89177715e-02,  1.84373426e-02, -3.80775886e-03,&
    -1.18284624e-03,  5.75833588e-04, -9.70777644e-05,&
     9.18691271e-06, -5.35799299e-07,  1.51653965e-08/)
    phys%alpha_energy_iz(:,4) = (/7.82341508e-03, -7.52650697e-03,  2.10882029e-03,&
    -1.06633298e-05, -6.05390376e-05,  7.78447320e-06,&
    -5.44204893e-07,  4.79775889e-08, -2.40410533e-09/)
    phys%alpha_energy_iz(:,5) = (/-1.63154981e-03,  1.44548221e-03, -4.15664835e-04,&
    2.94249910e-05, -4.22728571e-07,  7.20587415e-07,&
   -8.99016149e-08, -1.63703056e-09,  3.75320669e-10/)
    phys%alpha_energy_iz(:,6) = (/ 1.88643572e-04, -1.43008955e-04,  3.40709882e-05,&
    -8.97623559e-07, -2.35247387e-09, -9.88325233e-08,&
     8.93347501e-09,  5.73579811e-10, -6.00200317e-11/)
    phys%alpha_energy_iz(:,7) = (/-1.22470024e-05,  7.22137655e-06, -8.30018517e-07,&
    -3.00628338e-07,  4.54414015e-08,  2.57712484e-09,&
    -1.93796015e-10, -6.97804260e-11,  5.01115698e-12/)
    phys%alpha_energy_iz(:,8) = (/ 4.17042724e-07, -1.63374822e-07, -2.67173046e-08,&
    2.40019846e-08, -2.85457989e-09, -1.00767541e-10,&
    1.66441312e-11,  1.83602098e-12, -1.64199582e-13/)
    phys%alpha_energy_iz(:,9) = (/-5.77555689e-09,  1.02821579e-09,  1.14006243e-09,&
    -4.90129725e-10,  3.86861210e-11,  6.21520030e-12,&
    -9.52124215e-13,  1.88492217e-14,  1.25729669e-15/)
    !AMJUEL energy losses
    phys%alpha_energy_rec(:,1) = (/-2.59245035e+01, -7.29067024e-01,  2.36392587e-02,&
    3.64533393e-03,  1.59418465e-03, -1.21666803e-03,&
    2.37611590e-04, -1.93097764e-05,  5.59925778e-07/)
    phys%alpha_energy_rec(:,2) = (/1.22209727e-02, -1.54032393e-02,  1.16445335e-02,&
    -1.00582079e-03, -1.58223801e-05, -3.50307014e-04,&
     1.17270978e-04, -1.31840149e-05,  4.97782332e-07/)
    phys%alpha_energy_rec(:,3) = (/4.27849940e-05, -3.40609378e-03, -5.84520933e-03,&
    6.95635227e-04,  4.07369562e-04,  1.04350030e-04,&
   -6.69518205e-05,  8.84802545e-06, -3.61501382e-07/)
    phys%alpha_energy_rec(:,4) = (/ 1.94396774e-03,  1.53224343e-03,  2.85414587e-03,&
    -9.30505637e-04, -9.37916924e-05,  9.53616277e-06,&
     1.18818401e-05, -2.07237071e-06,  9.46698931e-08/)
    phys%alpha_energy_rec(:,5) = (/-7.12347460e-04, -4.65842377e-04, -5.07748529e-04,&
    2.58489629e-04,  1.49089050e-06, -6.90868188e-06,&
   -4.38151436e-07,  2.05591999e-07, -1.14648523e-08/)
    phys%alpha_energy_rec(:,6) = (/1.30352340e-04,  5.97244875e-05,  4.21110664e-05,&
    -3.29464390e-05,  2.24529287e-06,  8.23201901e-07,&
    -6.93626717e-08, -7.48963265e-09,  6.77233892e-10/)
    phys%alpha_energy_rec(:,7) = (/-1.18656075e-05, -4.07084329e-06, -1.25143662e-06,&
    2.11292402e-06, -3.15090101e-07, -2.90533105e-08,&
    6.59224926e-09, -7.07379703e-11, -1.77649634e-11/)
    phys%alpha_energy_rec(:,8) = (/5.33445563e-07,  1.37870988e-07, -1.62655575e-08,&
    -6.54468284e-08,  1.63196564e-08, -3.16903852e-10,&
    -1.77888796e-10,  1.04708751e-11,  7.19919506e-14/)
    phys%alpha_energy_rec(:,9) = (/-9.34985789e-09, -1.81807973e-09,  1.07345881e-09,&
    7.81029308e-10, -2.98409303e-10,  2.44276577e-11,&
    1.16076211e-12, -1.87744627e-13,  3.92930028e-15/)
    !Cooling factor for Nitrogen. 1D fit in loglog space for ADAS data in coronal limit, fitted in the range of 0.2 eV to 4e3 eV
    IF (phys%impurity_name == 'N') THEN
      phys%alpha_cooling_factor = (/-2.49348163e+01,  9.52628451e+00, -4.39511346e+00,  2.00446916e+00,&
      -2.27166819e+00,  9.95587194e-01,  1.11209167e+00, -1.13089140e+00,&
       2.22902131e-01,  1.37491696e-01, -9.74320916e-02,  2.88337222e-02,&
      -5.03419501e-03,  5.54252119e-04, -3.80153983e-05,  1.49061146e-06,&
      -2.56122449e-08/)
    ELSEIF (phys%impurity_name == 'W') THEN
      phys%alpha_cooling_factor = (/-3.29798537e+01,  3.38045169e+01, -3.81202398e+01,  2.47450333e+01,&
      -9.04921504e+00,  1.91417658e+00, -2.32405139e-01,  1.46022765e-02,&
      -2.46177645e-04, -1.89159754e-05,  7.73790432e-07,  0.,&
      0.,  0., 0.,  0.,&
      0./)
    ELSE 
      WRITE(6,*) 'Warning: cooling factor not defined for impurity ', TRIM(phys%impurity_name)
      STOP
    ENDIF
    ! coefficients for AMJUEL 2.1.8JH
#ifdef THREEBODYREC
    phys%alpha_rec(:,1) = (/-2.85572848e+01, -7.66404261e-01, -4.93042400e-03,&
                            -5.38683098e-03, -1.62603924e-04,  6.08090765e-06,&
                             2.10110205e-05, -2.77071760e-06,  1.03823594e-07/)
    phys%alpha_rec(:,2) = (/3.48856323e-02, -3.58323337e-03, -3.62024535e-03,&
                           -9.53284048e-04,  1.88804863e-04, -1.01489068e-05,&
                            2.24567656e-05, -4.69598237e-06,  2.52316661e-07/)
    phys%alpha_rec(:,3) = (/-2.79964439e-02, -7.45251429e-03,  6.95871196e-03,&
                             4.63175381e-04,  1.28857769e-04, -1.14502889e-04,&
                            -2.24562427e-06,  3.25087887e-06, -2.14539040e-07/)
    phys%alpha_rec(:,4) = (/1.20954532e-02,  2.70929976e-03, -2.13925730e-03,&
                           -5.37117970e-04, -1.63458052e-05,  5.94219398e-05,&
                           -2.94487376e-06, -9.38729079e-07,  7.38143524e-08/)
    phys%alpha_rec(:,5) = (/-2.43663080e-03, -7.74512977e-04,  4.60388371e-04,&
                            1.54335050e-04, -9.60103695e-06, -1.21185172e-05,&
                            1.00210510e-06,  1.39239163e-07, -1.29971368e-08/)
    phys%alpha_rec(:,6) = (/2.83789372e-04,  1.14244470e-04, -5.99163684e-05,&
                           -2.25756584e-05,  3.42526239e-06,  1.11896550e-06,&
                           -1.29132080e-07, -1.13909329e-08,  1.26518958e-09/)
    phys%alpha_rec(:,7) = (/-1.88651117e-05, -9.38278352e-06,  4.72926255e-06,&
                             1.73078295e-06, -4.07701994e-07, -4.27532157e-08,&
                             7.78615546e-09,  5.17850560e-10, -6.85420397e-11/)
    phys%alpha_rec(:,8) = (/6.75215560e-07,  3.90280010e-07, -1.99348540e-07,&
                           -6.61824078e-08,  2.04204110e-08,  3.70861611e-10,&
                           -2.44112778e-10, -9.45240216e-12,  1.83661503e-12/)
    phys%alpha_rec(:,9) = (/-1.00589386e-08, -6.38741159e-09,  3.35258987e-09,&
                             1.01336428e-09, -3.70797772e-10,  7.06845011e-12,&
                             3.77320848e-12, -4.67272402e-14, -1.64049236e-14/)
    ! coefficients for AMJUEL 2.1.8a
#else
     phys%alpha_rec(:,1) = (/-2.86177956e+01, -7.25199707e-01, -1.73502332e-02,&
     -3.55775280e-03, -2.77788226e-04,  2.06029540e-05,&
      1.59323839e-05, -2.11658076e-06,  7.66599010e-08/)
     phys%alpha_rec(:,2) = (/-1.78616692e-02,  3.21096605e-03, -3.11251743e-03,&
     1.55896611e-03, -9.32993286e-05, -1.28371165e-04,&
     3.70550340e-05, -3.85417246e-06,  1.40078912e-07/)
     phys%alpha_rec(:,3) = (/6.39155334e-04,  4.55025150e-03,  1.07786335e-03,&
     -1.03733153e-03,  1.09633177e-04,  7.31231189e-05,&
     -2.40723586e-05,  2.66239203e-06, -1.00895147e-07/)
     phys%alpha_rec(:,4) = (/-4.50941526e-04, -1.88230646e-03, -2.61695897e-04,&
     2.81723717e-04, -4.56748839e-05, -1.06480515e-05,&
     4.91521392e-06, -6.12084620e-07,  2.49521491e-08/)
     phys%alpha_rec(:,5) = (/7.09545902e-05,  3.98313304e-04,  5.45933281e-05,&
     -4.40781517e-05,  8.49578724e-06, -1.49877643e-07,&
     -3.34660940e-07,  5.66372822e-08, -2.67848413e-09/)
     phys%alpha_rec(:,6) = (/-5.66030993e-06, -4.85183529e-05, -8.63530868e-06,&
     4.64601735e-06, -7.26107627e-07,  1.19908760e-07,&
    -4.91275369e-09, -1.47422116e-09,  1.17013833e-10/)
     phys%alpha_rec(:,7) = (/1.16018663e-07,  3.40483450e-06,  8.38310637e-07,&
     -3.36565455e-07,  2.32699294e-08, -5.66807913e-09,&
      1.30239368e-09, -7.37309518e-11, -1.58825470e-13/)
     phys%alpha_rec(:,8) = (/7.56498607e-09, -1.28083999e-07, -4.13335200e-08,&
     1.42835079e-08,  2.20808955e-10, -1.01855404e-10,&
    -3.16901361e-11,  4.31445723e-12, -1.22634522e-13/)
     phys%alpha_rec(:,9) = (/-2.96981503e-10,  1.98283997e-09,  7.87249173e-10,&
     -2.52215335e-10, -1.98997939e-11,  7.76657896e-12,&
     -1.78376276e-13, -4.79167750e-14,  2.32940245e-15/)
#endif
#endif
#ifdef NEUTRAL
    ! Allocate atomic data
    phys%E = (/0, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 250, 270, 300, 350, 400, 500, 600, 700, 1000/)
    phys%theta = (/0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90/)
    
    phys%RN_DW(1,:) = (/0.7987, 0.8017, 0.8047, 0.8078, 0.8169, 0.8260, 0.8351, 0.8482, 0.8612, 0.8743, 0.8928,&
                        0.9114, 0.9261, 0.9409, 0.9514, 0.9620, 0.9682, 0.9716, 0.9751/)
    phys%RN_DW(2,:) = (/0.7949, 0.7980, 0.8010, 0.8040, 0.8131, 0.8223, 0.8314, 0.8446, 0.8578, 0.8710, 0.8898,&
                        0.9085, 0.9238, 0.9391, 0.9504, 0.9618, 0.9687, 0.9725, 0.9763/)
    phys%RN_DW(3,:) = (/0.7614, 0.7643, 0.7672, 0.7702, 0.7793, 0.7885, 0.7977, 0.8125, 0.8272, 0.8420, 0.8622,&
                        0.8824, 0.9029, 0.9234, 0.9417, 0.9601, 0.9730, 0.9800, 0.9871 /)
    phys%RN_DW(4,:) = (/0.7242, 0.7270, 0.7297, 0.7325, 0.7418, 0.7511, 0.7604, 0.7768, 0.7933, 0.8097, 0.8316,&
                        0.8534, 0.8797, 0.9059, 0.9321, 0.9582, 0.9778, 0.9885, 0.9991 /)
    phys%RN_DW(5,:) = (/0.7091, 0.7119, 0.7146, 0.7174, 0.7265, 0.7357, 0.7448, 0.7607, 0.7766, 0.7925, 0.8148,&
                        0.8371, 0.8645, 0.8919, 0.9220, 0.9520, 0.9765, 0.9900, 1.0000 /)
    phys%RN_DW(6,:) = (/0.6941, 0.6968, 0.6995, 0.7023, 0.7113, 0.7203, 0.7293, 0.7446, 0.7599, 0.7753, 0.7981,&
                        0.8208, 0.8494, 0.8780, 0.9118, 0.9457, 0.9752, 0.9916, 1.0000 /)
    phys%RN_DW(7,:) = (/0.6790, 0.6817, 0.6844, 0.6872, 0.6960, 0.7049, 0.7138, 0.7285, 0.7433, 0.7580, 0.7813,&
                        0.8045, 0.8343, 0.8640, 0.9017, 0.9395, 0.9738, 0.9932, 1.0000 /)
    phys%RN_DW(8,:) = (/0.6730, 0.6755, 0.6781, 0.6806, 0.6894, 0.6983, 0.7072, 0.7220, 0.7369, 0.7518, 0.7740,&
                        0.7962, 0.8260, 0.8558, 0.8945, 0.9333, 0.9711, 0.9932, 1.0000 /)
    phys%RN_DW(9,:) = (/0.6671, 0.6694, 0.6717, 0.6740, 0.6828, 0.6917, 0.7005, 0.7156, 0.7306, 0.7456, 0.7667,& 
                        0.7879, 0.8178, 0.8476, 0.8874, 0.9271, 0.9683, 0.9933, 1.0000 /)
    phys%RN_DW(10,:) = (/0.6611, 0.6632, 0.6653, 0.6674, 0.6762, 0.6851, 0.6939, 0.7091, 0.7242, 0.7394, 0.7595,&
                         0.7795, 0.8095, 0.8395, 0.8802, 0.9209, 0.9655, 0.9933, 1.0000 /)
    phys%RN_DW(11,:) = (/0.6552, 0.6570, 0.6589, 0.6608, 0.6696, 0.6785, 0.6873, 0.7026, 0.7179, 0.7332, 0.7522,&
                         0.7712, 0.8013, 0.8313, 0.8730, 0.9147, 0.9627, 0.9933, 1.0000 /)
    phys%RN_DW(12,:) = (/0.6492, 0.6509, 0.6525, 0.6542, 0.6630, 0.6719, 0.6807, 0.6961, 0.7115, 0.7269, 0.7449,&
                         0.7629, 0.7930, 0.8232, 0.8659, 0.9085, 0.9599, 0.9934, 1.0000 /)
    phys%RN_DW(13,:) = (/0.6117, 0.6157, 0.6197, 0.6237, 0.6300, 0.6364, 0.6427, 0.6584, 0.6741, 0.6898, 0.7107,&
                         0.7315, 0.7588, 0.7861, 0.8273, 0.8685, 0.9315, 0.9904, 1.0000 /)
    phys%RN_DW(14,:) = (/0.6024, 0.6053, 0.6081, 0.6110, 0.6196, 0.6282, 0.6368, 0.6513, 0.6658, 0.6804, 0.7006,&
                         0.7209, 0.7475, 0.7740, 0.8144, 0.8547, 0.9199, 0.9882, 1.0000 /)
    phys%RN_DW(15,:) = (/0.5981, 0.6010, 0.6039, 0.6069, 0.6156, 0.6243, 0.6330, 0.6476, 0.6622, 0.6768, 0.6969,&
                         0.7170, 0.7436, 0.7702, 0.8102, 0.8502, 0.9153, 0.9872, 1.0000 /)
    phys%RN_DW(16,:) = (/0.5924, 0.5954, 0.5985, 0.6015, 0.6102, 0.6189, 0.6276, 0.6423, 0.6569, 0.6715, 0.6917,&
                         0.7120, 0.7382, 0.7644, 0.8041, 0.8438, 0.9091, 0.9857, 1.0000 /)
    phys%RN_DW(17,:) = (/0.5839, 0.5869, 0.5900, 0.5930, 0.6019, 0.6108, 0.6196, 0.6344, 0.6492, 0.6640, 0.6840,&
                         0.7041, 0.7303, 0.7565, 0.7954, 0.8344, 0.8997, 0.9829, 1.0000 /)
    phys%RN_DW(18,:) = (/0.5767, 0.5794, 0.5821, 0.5848, 0.5939, 0.6031, 0.6122, 0.6273, 0.6423, 0.6574, 0.6774,&
                         0.6974, 0.7235, 0.7497, 0.7882, 0.8267, 0.8908, 0.9799, 1.0000 /)
    phys%RN_DW(19,:) = (/0.5639, 0.5669, 0.5698, 0.5728, 0.5819, 0.5911, 0.6002, 0.6150, 0.6298, 0.6445, 0.6653,&
                         0.6860, 0.7119, 0.7378, 0.7758, 0.8137, 0.8766, 0.9739, 1.0000 /)
    phys%RN_DW(20,:) = (/0.5528, 0.5555, 0.5582, 0.5609, 0.5700, 0.5791, 0.5881, 0.6039, 0.6197, 0.6355, 0.6566,&
                         0.6777, 0.7032, 0.7288, 0.7658, 0.8028, 0.8649, 0.9671, 1.0000 /)
    phys%RN_DW(21,:) = (/0.5425, 0.5451, 0.5477, 0.5503, 0.5599, 0.5694, 0.5790, 0.5943, 0.6097, 0.6250, 0.6466,&
                         0.6683, 0.6952, 0.7222, 0.7581, 0.7941, 0.8554, 0.9609, 1.0000 /)
    phys%RN_DW(22,:) = (/0.5174, 0.5210, 0.5247, 0.5283, 0.5379, 0.5474, 0.5570, 0.5730, 0.5890, 0.6050, 0.6253,&
                         0.6456, 0.6741, 0.7026, 0.7381, 0.7737, 0.8327, 0.9419, 1.0000/)
#endif

  ! Bohm-GyroBohm definitions
  phys%c_bohm_i = 1.6e-4
  phys%c_bohm_e = 8.e-5
  phys%c_gyroBohm_i = 1.75e-2
  phys%c_gyroBohm_e = 3.5e-2
  phys%prandtl = 1.
  phys%c_bohm_n = 1.

  call initialize_neutral_rate_runtime_constants()
  call adimensionalize_neutral_rate_coefficients()

  ENDSUBROUTINE initPhys

  SUBROUTINE initialize_neutral_rate_runtime_constants()
    REAL*8 :: density_scale

    density_scale = simpar%refval_density/1.d14
    neutral_rt%rate_scale = simpar%refval_density*simpar%refval_time
    neutral_rt%energy_weight = phys%Mref/simpar%refval_temperature

    neutral_rt%eirene_te_min = eirene_rate_te_min_phys/simpar%refval_temperature
    neutral_rt%eirene_te_max = eirene_rate_te_max_phys/simpar%refval_temperature
    neutral_rt%eirene_ti_min = eirene_rate_ti_min_phys/simpar%refval_temperature
    neutral_rt%eirene_ti_max = eirene_rate_ti_max_phys/simpar%refval_temperature
    neutral_rt%eirene_ne_min = eirene_rate_ne_min_phys/1.d14/density_scale
    neutral_rt%eirene_ne_max = eirene_rate_ne_max_phys/1.d14/density_scale

    neutral_rt%state_tol = neutral_state_tol_default
    neutral_rt%te_floor = neutral_te_floor_phys/simpar%refval_temperature
    neutral_rt%ti_floor = neutral_ti_floor_phys/simpar%refval_temperature
    neutral_rt%sigmavnn_prefactor = neutral_rt%rate_scale*neutral_sigmavnn_prefactor_phys* &
      &(simpar%refval_temperature*elementary_charge_si/boltzmann_constant_si)**0.25d0
    neutral_rt%transport_ti_floor = 1.d-6/simpar%refval_temperature
    neutral_rt%transport_ti_supp = neutral_rt%transport_ti_floor
    neutral_rt%rydberg_energy = neutral_rydberg_energy_phys/simpar%refval_temperature
    neutral_rt%iz_te_floor = neutral_iz_te_floor_phys/simpar%refval_temperature
    neutral_rt%rec_te_floor = neutral_rec_te_floor_phys/simpar%refval_temperature
    neutral_rt%cx_te_floor = neutral_cx_te_floor_phys/simpar%refval_temperature
    neutral_rt%log_temperature_ref = LOG(simpar%refval_temperature)
    neutral_rt%tloss_offset = neutral_tloss_offset_phys*neutral_rt%energy_weight
    neutral_rt%tloss_amplitude = neutral_tloss_amplitude_phys*neutral_rt%energy_weight
    neutral_rt%tloss_decay = 0.5d0*simpar%refval_temperature
    neutral_rt%tlossrec_prefactor = neutral_tlossrec_prefactor_phys*neutral_rt%energy_weight
    neutral_rt%tlossrec_cap = neutral_tlossrec_cap_phys*neutral_rt%energy_weight
    neutral_rt%tlossrec_growth = simpar%refval_temperature/9.d0
    neutral_rt%tlossrec_cap_log = LOG(neutral_tlossrec_cap_phys/neutral_tlossrec_prefactor_phys)
    neutral_rt%legacy_cx_prefactor = neutral_rt%rate_scale*neutral_legacy_cx_prefactor_phys
    neutral_rt%recombination_energy = neutral_rydberg_energy_phys*neutral_rt%energy_weight
  ENDSUBROUTINE initialize_neutral_rate_runtime_constants

  SUBROUTINE adimensionalize_neutral_rate_coefficients()
    REAL*8 :: log_temp_shift, log_density_shift
    REAL*8 :: log_rate_scale, log_energy_rate_scale

    log_temp_shift = LOG(simpar%refval_temperature)
    log_density_shift = LOG(simpar%refval_density/1.d14)
    log_rate_scale = LOG(simpar%refval_density*simpar%refval_time)
    log_energy_rate_scale = LOG(simpar%refval_density*simpar%refval_time* &
      &simpar%refval_charge/simpar%refval_mass*simpar%refval_time**2/simpar%refval_length**2)

#ifdef AMJUELSPLINES
    call shift_logpoly_2d_9x9(phys%alpha_iz, log_temp_shift, log_density_shift, log_rate_scale)
    call shift_logpoly_2d_9x9(phys%alpha_rec, log_temp_shift, log_density_shift, log_rate_scale)
    call shift_logpoly_2d_9x9(phys%alpha_energy_iz, log_temp_shift, log_density_shift, log_energy_rate_scale)
    call shift_logpoly_2d_9x9(phys%alpha_energy_rec, log_temp_shift, log_density_shift, log_energy_rate_scale)
#endif
#ifdef EXPANDEDCX
#ifdef AMJUELCX
    call shift_logpoly_1d_9(phys%alpha_cx, log_temp_shift, log_rate_scale)
#endif
#ifdef THERMALCX
    call shift_logpoly_1d_5(phys%alpha_cx, log_temp_shift, log_rate_scale)
#endif
    call shift_logpoly_1d_17(phys%alpha_cooling_factor, log_temp_shift, log_energy_rate_scale)
#endif
  ENDSUBROUTINE adimensionalize_neutral_rate_coefficients

  SUBROUTINE shift_logpoly_2d_9x9(alpha, shift_x, shift_y, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(9,9)
    REAL*8, INTENT(IN)    :: shift_x, shift_y, log_scale

    call shift_logpoly_2d(alpha, shift_x, shift_y, log_scale)
  ENDSUBROUTINE shift_logpoly_2d_9x9

  SUBROUTINE shift_logpoly_1d_9(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(9)
    REAL*8, INTENT(IN)    :: shift_x, log_scale

    call shift_logpoly_1d(alpha, shift_x, log_scale)
  ENDSUBROUTINE shift_logpoly_1d_9

  SUBROUTINE shift_logpoly_1d_5(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(5)
    REAL*8, INTENT(IN)    :: shift_x, log_scale

    call shift_logpoly_1d(alpha, shift_x, log_scale)
  ENDSUBROUTINE shift_logpoly_1d_5

  SUBROUTINE shift_logpoly_1d_17(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(17)
    REAL*8, INTENT(IN)    :: shift_x, log_scale

    call shift_logpoly_1d(alpha, shift_x, log_scale)
  ENDSUBROUTINE shift_logpoly_1d_17

  SUBROUTINE shift_logpoly_2d(alpha, shift_x, shift_y, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(:,:)
    REAL*8, INTENT(IN)    :: shift_x, shift_y, log_scale
    REAL*8                :: shifted(size(alpha,1), size(alpha,2))
    REAL*8                :: x_factor, y_factor
    INTEGER               :: i, j, ip, jp, pow_x, pow_y

    shifted = 0.d0
    do j = 1, size(alpha,2)
      do i = 1, size(alpha,1)
        do jp = 1, j
          pow_y = (j - jp)
          y_factor = 1.d0
          if (pow_y > 0) y_factor = shift_y**pow_y
          do ip = 1, i
            pow_x = (i - ip)
            x_factor = 1.d0
            if (pow_x > 0) x_factor = shift_x**pow_x
            shifted(ip,jp) = shifted(ip,jp) + alpha(i,j)* &
              &binomial_coefficient(i - 1, ip - 1)*binomial_coefficient(j - 1, jp - 1)* &
              &x_factor*y_factor
          end do
        end do
      end do
    end do

    shifted(1,1) = shifted(1,1) + log_scale
    alpha = shifted
  ENDSUBROUTINE shift_logpoly_2d

  SUBROUTINE shift_logpoly_1d(alpha, shift_x, log_scale)
    REAL*8, INTENT(INOUT) :: alpha(:)
    REAL*8, INTENT(IN)    :: shift_x, log_scale
    REAL*8                :: shifted(size(alpha))
    REAL*8                :: x_factor
    INTEGER               :: i, ip, pow_x

    shifted = 0.d0
    do i = 1, size(alpha)
      do ip = 1, i
        pow_x = (i - ip)
        x_factor = 1.d0
        if (pow_x > 0) x_factor = shift_x**pow_x
        shifted(ip) = shifted(ip) + alpha(i)*binomial_coefficient(i - 1, ip - 1)*x_factor
      end do
    end do

    shifted(1) = shifted(1) + log_scale
    alpha = shifted
  ENDSUBROUTINE shift_logpoly_1d

  REAL*8 FUNCTION binomial_coefficient(n, k)
    INTEGER, INTENT(IN) :: n, k
    INTEGER             :: i, kk

    if (k < 0 .or. k > n) then
      binomial_coefficient = 0.d0
      return
    end if

    if (k == 0 .or. k == n) then
      binomial_coefficient = 1.d0
      return
    end if

    kk = MIN(k, n - k)
    binomial_coefficient = 1.d0
    do i = 1, kk
      binomial_coefficient = binomial_coefficient*DBLE(n - kk + i)/DBLE(i)
    end do
  ENDFUNCTION binomial_coefficient

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  PURE SUBROUTINE phys2cons(up, ua)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: up
    REAL*8, DIMENSION(:, :), INTENT(out) :: ua

    ua(:, 1) = ABS(up(:, 1))
    ua(:, 2) = up(:, 1)*up(:, 2)
    ua(:, 3) = up(:, 1)*up(:, 3)
    ua(:, 4) = up(:, 1)*up(:, 4)
#ifdef NEUTRAL
    if (phys%idx_rhon_eq > 0 .and. phys%idx_rhon_pv > 0) ua(:,phys%idx_rhon_eq) = ABS(up(:,phys%idx_rhon_pv))
#ifdef NEUTRALGAMMA
    if (phys%idx_gamman_eq > 0 .and. phys%idx_un_pv > 0 .and. phys%idx_rhon_pv > 0) then
      ua(:,phys%idx_gamman_eq) = ABS(MAX(up(:,phys%idx_rhon_pv), 1.d-7))*up(:,phys%idx_un_pv)
    end if
#endif
#ifdef KEQUATION
    if (phys%idx_k_eq > 0 .and. phys%idx_k_pv > 0) ua(:,phys%idx_k_eq) = ABS(up(:,phys%idx_k_pv))
#endif
#endif

  ENDSUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  PURE SUBROUTINE cons2phys(ua, up)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: ua
    REAL*8, DIMENSION(:, :), INTENT(out) :: up
    REAL*8,  DIMENSION(SIZE(ua,1))       :: U1
    REAL*8,  DIMENSION(SIZE(ua,1))       :: U5

    U1 = max(ua(:,1),1e-20)
    U5 = 0.d0
#ifdef NEUTRAL
    if (phys%idx_rhon_eq > 0) U5 = ABS(ua(:,phys%idx_rhon_eq))
#endif

    up(:, 1) = ABS(U1)                                                           ! density
    up(:, 2) = ua(:, 2)/U1                                            ! u parallel
    up(:, 3) = ua(:, 3)/U1                                            ! total energy of ions
    up(:, 4) = ua(:, 4)/U1                                            ! total energy of electrons
    up(:, 5) = (2./(3.*phys%Mref)*(ua(:, 3) - 0.5*ua(:, 2)**2/U1)) ! pressure of ions
    up(:, 6) = (2./(3.*phys%Mref)*ua(:, 4))                              ! pressure of electrons
    up(:, 7) = max(up(:, 5)/U1,1e-20)                                           ! temperature of ions
    up(:, 8) = max(up(:, 6)/U1,1e-20)                                         ! temperature of electrons
    up(:, 9) = SQRT(max((up(:, 5) + up(:, 6))/U1*phys%Mref,1e-20))                        ! sound speed
    up(:, 10) = up(:, 2)/up(:, 9)                                           ! Mach
#ifdef NEUTRAL
    if (phys%idx_rhon_eq > 0 .and. phys%idx_rhon_pv > 0) up(:,phys%idx_rhon_pv) = ABS(ua(:,phys%idx_rhon_eq)) ! density neutral
#ifdef NEUTRALGAMMA
    if (phys%idx_gamman_eq > 0 .and. phys%idx_un_pv > 0) then
      up(:,phys%idx_un_pv) = ua(:,phys%idx_gamman_eq)/ABS(MAX(U5, 1.d-7))
    end if
#endif
#ifdef KEQUATION
    if (phys%idx_k_eq > 0 .and. phys%idx_k_pv > 0) up(:,phys%idx_k_pv) = ABS(ua(:,phys%idx_k_eq)) ! turbulent energy
#endif
#endif


  ENDSUBROUTINE cons2phys

    ! ******************************
    ! Split diffusion terms
    ! ******************************
    SUBROUTINE compute_W2(U,W2,diff_n,diff_u)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_u
    REAL*8             :: W2(:)
      W2 = 0.
      W2(1) = (diff_n-diff_u)*U(2)/U(1)
    ENDSUBROUTINE compute_W2

    SUBROUTINE compute_W3(U,W3,diff_n,diff_u,diff_e)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_u, diff_e
    REAL*8             :: W3(:)
    REAL*8 		:: rhovar
    REAL*8 		:: sigmavar

      rhovar = (diff_e-diff_u)
      sigmavar = (diff_n-diff_e)

      W3 = 0.
      W3(1) = sigmavar*U(3)/U(1) + rhovar*(U(2)/U(1))**2
      W3(2) = -rhovar*U(2)/U(1)
    ENDSUBROUTINE compute_W3

    SUBROUTINE compute_W4(U,W4,diff_n,diff_ee)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_ee
    REAL*8             :: W4(:)

      W4 = 0.
      W4(1) = (diff_n-diff_ee)*U(4)/U(1)
    ENDSUBROUTINE compute_W4



    SUBROUTINE compute_dW2_dU(U,dW2_dU,diff_n,diff_u)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_u
    REAL*8             :: dW2_dU(:,:)
      dW2_dU = 0.
      dW2_dU(1,1) = -U(2)/(U(1)**2)
      dW2_dU(1,2) = 1./U(1)

      dW2_dU = (diff_n-diff_u)*dW2_dU
    ENDSUBROUTINE compute_dW2_dU

    SUBROUTINE compute_dW3_dU(U,res,diff_n,diff_u,diff_e)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_u,diff_e
    REAL*8             :: res(:,:)
    REAL*8 		:: rhovar
    REAL*8 		:: sigmavar

      rhovar = (diff_e-diff_u)
      sigmavar = (diff_n-diff_e)

      res = 0.
      res(1,1) = -sigmavar*U(3)/(U(1)**2)-2*rhovar*(U(2)**2)/(U(1)**3)
      res(1,2) = 2*rhovar*U(2)/(U(1)**2)
      res(1,3) = sigmavar*1./U(1)

      res(2,1) = rhovar*U(2)/(U(1)**2)
      res(2,2) = -rhovar*1./U(1)
    ENDSUBROUTINE compute_dW3_dU

    SUBROUTINE compute_dW4_dU(U,res,diff_n,diff_ee)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_ee
    REAL*8             :: res(:,:)
    res = 0.
    res(1,1) = -U(4)*(diff_n-diff_ee)/(U(1)**2)
    res(1,4) = 1.*(diff_n-diff_ee)/U(1)
  ENDSUBROUTINE compute_dW4_dU

  SUBROUTINE compute_Ti(U, Ti)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Ti

    Ti = 2.d0/(3.d0*phys%Mref)*(U(3)/U(1) - 0.5d0*(U(2)/U(1))**2)
  ENDSUBROUTINE compute_Ti

  SUBROUTINE compute_dTi_dU(U, dTi_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dTi_dU(:)

    CALL computeVi(U, dTi_dU)
    dTi_dU = dTi_dU*2.d0/(3.d0*phys%Mref)
  ENDSUBROUTINE compute_dTi_dU

  SUBROUTINE compute_Te(U, Te)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Te

    Te = 2.d0/(3.d0*phys%Mref)*U(4)/U(1)
  ENDSUBROUTINE compute_Te

  SUBROUTINE compute_dTe_dU(U, dTe_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dTe_dU(:)

    dTe_dU = 0.d0
    dTe_dU(1) = -2.d0/(3.d0*phys%Mref)*U(4)/U(1)**2
    dTe_dU(4) = 2.d0/(3.d0*phys%Mref)/U(1)
  ENDSUBROUTINE compute_dTe_dU

  SUBROUTINE compute_limited_Ti(U, Ti_limited)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Ti_limited

    CALL compute_Ti(U, Ti_limited)
    CALL softplus(Ti_limited, neutral_rt%transport_ti_floor)
  ENDSUBROUTINE compute_limited_Ti

  SUBROUTINE compute_dlimited_Ti_dU(U, dTi_limited_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dTi_limited_dU(:)
    REAL*8              :: Ti, soft_deriv
    REAL*8              :: dTi_dU(size(U))

    CALL compute_Ti(U, Ti)
    CALL compute_dTi_dU(U, dTi_dU)
    CALL softplus_deriv(Ti, neutral_rt%transport_ti_floor, soft_deriv)
    dTi_limited_dU = dTi_dU*soft_deriv
  ENDSUBROUTINE compute_dlimited_Ti_dU

  SUBROUTINE compute_neutral_transport_prefactor(U, coeff)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: coeff
    REAL*8              :: Ti_limited

    CALL compute_limited_Ti(U, Ti_limited)
    coeff = 0.5d0*phys%a*Ti_limited
  ENDSUBROUTINE compute_neutral_transport_prefactor

  SUBROUTINE compute_dneutral_transport_prefactor_dU(U, dcoeff_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dcoeff_dU(:)
    REAL*8              :: dTi_limited_dU(size(U))

    CALL compute_dlimited_Ti_dU(U, dTi_limited_dU)
    dcoeff_dU = 0.5d0*phys%a*dTi_limited_dU
  ENDSUBROUTINE compute_dneutral_transport_prefactor_dU

  SUBROUTINE compute_neutral_diffusion_denominator(U, denom)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: denom
    REAL*8              :: sigmaviz, sigmavnn, sigmavcx
    INTEGER             :: inn

    inn = phys%idx_rhon_eq
    CALL compute_sigmaviz(U, sigmaviz)
    CALL compute_sigmavnn(U, sigmavnn)
    CALL compute_sigmavcx(U, sigmavcx)
    denom = U(1)*(sigmaviz + sigmavcx) + U(inn)*sigmavnn
  ENDSUBROUTINE compute_neutral_diffusion_denominator

  SUBROUTINE compute_dneutral_diffusion_denominator_dU(U, ddenom_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: ddenom_dU(:)
    REAL*8              :: sigmaviz, sigmavnn, sigmavcx
    REAL*8              :: dsigmaviz_dU(size(U)), dsigmavnn_dU(size(U)), dsigmavcx_dU(size(U))
    INTEGER             :: inn

    inn = phys%idx_rhon_eq
    CALL compute_sigmaviz(U,sigmaviz)
    CALL compute_sigmavnn(U,sigmavnn)
    CALL compute_sigmavcx(U,sigmavcx)
    CALL compute_dsigmaviz_dU(U,dsigmaviz_dU)
    CALL compute_dsigmavnn_dU(U,dsigmavnn_dU)
    CALL compute_dsigmavcx_dU(U,dsigmavcx_dU)

    ddenom_dU = U(1)*(dsigmaviz_dU + dsigmavcx_dU) + U(inn)*dsigmavnn_dU
    ddenom_dU(1) = ddenom_dU(1) + sigmaviz + sigmavcx
    ddenom_dU(inn) = ddenom_dU(inn) + sigmavnn
  ENDSUBROUTINE compute_dneutral_diffusion_denominator_dU

#ifdef NEUTRALGAMMA
  SUBROUTINE compute_neutral_gamma_denominator(U, denom)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: denom
    REAL*8              :: sigmavcx, sigmavnn
    INTEGER             :: inn

    inn = phys%idx_rhon_eq
    CALL compute_sigmavcx(U, sigmavcx)
    CALL compute_sigmavnn(U, sigmavnn)
    denom = U(1)*sigmavcx + U(inn)*sigmavnn
  ENDSUBROUTINE compute_neutral_gamma_denominator

  SUBROUTINE compute_dneutral_gamma_denominator_dU(U, ddenom_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: ddenom_dU(:)
    REAL*8              :: sigmavcx, sigmavnn
    REAL*8              :: dsigmavcx_dU(size(U)), dsigmavnn_dU(size(U))
    INTEGER             :: inn

    inn = phys%idx_rhon_eq
    CALL compute_sigmavcx(U, sigmavcx)
    CALL compute_sigmavnn(U, sigmavnn)
    CALL compute_dsigmavcx_dU(U, dsigmavcx_dU)
    CALL compute_dsigmavnn_dU(U, dsigmavnn_dU)
    ddenom_dU = U(1)*dsigmavcx_dU + U(inn)*dsigmavnn_dU
    ddenom_dU(1) = ddenom_dU(1) + sigmavcx
    ddenom_dU(inn) = ddenom_dU(inn) + sigmavnn
  ENDSUBROUTINE compute_dneutral_gamma_denominator_dU
#endif

  SUBROUTINE compute_Dnn(U, Dnn)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Dnn
    REAL*8              :: coeff, denom

#ifdef CONSTANTNEUTRALDIFF
    Dnn = MAX(phys%diff_nn, phys%diff_nn_min)
#else
    CALL compute_neutral_transport_prefactor(U, coeff)
    CALL compute_neutral_diffusion_denominator(U, denom)
    Dnn = coeff/denom
    CALL double_softplus(Dnn, phys%diff_nn_min, phys%diff_nn)
#endif
  ENDSUBROUTINE compute_Dnn

  SUBROUTINE compute_neutral_free_streaming_speed(U, cn)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: cn
    REAL*8              :: Ti_limited

    CALL compute_limited_Ti(U, Ti_limited)
    cn = SQRT(MAX(phys%Mref*Ti_limited, 0.d0))
  ENDSUBROUTINE compute_neutral_free_streaming_speed

  SUBROUTINE compute_neutral_flux_limiter(U, Q, phi, gamma_unlim, gamma_max, ratio)
    REAL*8, INTENT(IN)  :: U(:), Q(:)
    REAL*8, INTENT(OUT) :: phi, gamma_unlim(:), gamma_max, ratio
    REAL*8              :: Qpr(simpar%Ndim, simpar%Neq)
    REAL*8              :: Dnn, cn, gamma_abs_eps
#ifdef NEUTRALPNEW
    REAL*8              :: W5p(simpar%Neq)
#endif
    INTEGER             :: inn

    inn = phys%idx_rhon_eq
    gamma_unlim = 0.d0
    gamma_max = 0.d0
    ratio = 0.d0
    phi = 1.d0
    IF (inn <= 0) RETURN

    Qpr = RESHAPE(Q, (/simpar%Ndim, simpar%Neq/))
    CALL compute_Dnn(U, Dnn)
    gamma_unlim = -Dnn*Qpr(:,inn)
#ifdef NEUTRALPNEW
    CALL compute_W5p(U, W5p)
    gamma_unlim = gamma_unlim - MATMUL(Qpr, W5p)
#endif

    CALL compute_neutral_free_streaming_speed(U, cn)
    gamma_max = MAX(phys%neutral_flux_limiter_fs_fraction*MAX(U(inn), 0.d0)*cn, &
      &phys%neutral_flux_limiter_fs_flux_min)
    gamma_abs_eps = SQRT(DOT_PRODUCT(gamma_unlim, gamma_unlim) + phys%neutral_flux_limiter_eps**2)

    IF (gamma_max > 0.d0) THEN
      ratio = gamma_abs_eps/gamma_max
    ELSEIF (gamma_abs_eps > 0.d0) THEN
      ratio = HUGE(1.d0)
    ELSE
      ratio = 0.d0
    ENDIF

    IF (ratio < HUGE(1.d0)) THEN
      phi = 1.d0/(1.d0 + ratio)
    ELSE
      phi = 0.d0
    ENDIF
    IF (Dnn > 0.d0) phi = MAX(phi, MIN(1.d0, phys%diff_nn_min/Dnn))
  ENDSUBROUTINE compute_neutral_flux_limiter

  !*****************************************
  ! Jacobian matrices
  !****************************************
  SUBROUTINE jacobianMatrices(U, A)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: A(:, :)
    REAL*8              :: Unn
    INTEGER             :: ik, ign, inn
    ![ 0,                                           1,                              0,                0; ...
    ! -2/3*U(2)**2/U(1)**2                          4/3*U(2)/U(1)                   2/3,              2/3; ...
    ! -5/3*U(2)*U(3)/U(1)**2+2/3*U(2)**3/U(1)**3    5/3*U(3)/U(1)-U(2)**2/U(1)**2   5/3*U(2)/U(1),    0 ;   ...
    ! -5/3*U(4)*U(2)/U(1)**2,                       5/3*U(4)/U(1),                  0,                5/3*U(2)/U(1)]
    ! k equation line
    ! -U(6)*U(2)/U(1)**2,                           U(6)/U(1),                      0,                0,            0,        U(2)/U(1)]
    A = 0.d0
    ik = phys%idx_k_eq
    ign = phys%idx_gamman_eq
    inn = phys%idx_rhon_eq
    IF (switch%decoup) THEN
      A(1, 2) = 1.

      A(2, 1) = (-U(2)**2/U(1)**2 + phys%Mref)
      A(2, 2) = 2.*U(2)/U(1)

      A(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      A(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      A(3, 3) = 5./3.*U(2)/U(1)

      A(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      A(4, 2) = 5./3.*U(4)/U(1)
      A(4, 4) = 5./3.*U(2)/U(1)
#ifdef KEQUATION
      if (ik > 0) then
        A(ik, 1) = -U(ik)*U(2)/U(1)**2
        A(ik, 2) = U(ik)/U(1)
        A(ik, ik) = U(2)/U(1)
      end if
#endif
#ifdef NEUTRAL
#ifdef NEUTRALGAMMA
      if (inn > 0 .and. ign > 0) then
        Unn = MAX(1.d-7,U(inn))
        A(inn, ign) = 1.d0
        A(ign, 1) = 2.d0/3.d0*Unn*(-U(3)/U(1)**2 + U(2)**2/U(1)**3)
        A(ign, 2) = -2.d0/3.d0*Unn*U(2)/U(1)**2
        A(ign, 3) = 2.d0/3.d0*Unn/U(1)
        A(ign, inn) = -U(ign)**2/Unn**2 + 2.d0/3.d0*(U(3)/U(1) - 0.5d0*U(2)**2/U(1)**2)
        A(ign, ign) = 2.d0*U(ign)/Unn
      end if
#endif
#endif
    ELSE

      A(1, 2) = 1.

      A(2, 1) = -2./3.*U(2)**2/U(1)**2
      A(2, 2) = 4./3.*U(2)/U(1)
      A(2, 3) = 2./3.
      A(2, 4) = 2./3.

      A(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      A(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      A(3, 3) = 5./3.*U(2)/U(1)

      A(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      A(4, 2) = 5./3.*U(4)/U(1)
      A(4, 4) = 5./3.*U(2)/U(1)
#ifdef KEQUATION
      if (ik > 0) then
        A(ik, 1) = -U(ik)*U(2)/U(1)**2
        A(ik, 2) = U(ik)/U(1)
      !IF (U(6) >= 0.) THEN
          A(ik, ik) = U(2)/U(1)
      !ENDIF
      end if
#endif
#ifdef NEUTRAL
#ifdef NEUTRALGAMMA
      if (inn > 0 .and. ign > 0) then
        Unn = MAX(1.d-7,U(inn))
        A(inn, ign) = 1.d0
        A(ign, 1) = 2.d0/3.d0*Unn*(-U(3)/U(1)**2 + U(2)**2/U(1)**3)
        A(ign, 2) = -2.d0/3.d0*Unn*U(2)/U(1)**2
        A(ign, 3) = 2.d0/3.d0*Unn/U(1)
        A(ign, inn) = -U(ign)**2/Unn**2 + 2.d0/3.d0*(U(3)/U(1) - 0.5d0*U(2)**2/U(1)**2)
        A(ign, ign) = 2.d0*U(ign)/Unn
      end if
#endif
#endif
    END IF
  ENDSUBROUTINE jacobianMatrices

  !*****************************************
  ! Jacobian matrix for face computations
  !****************************************
  SUBROUTINE jacobianMatricesFace(U, bn, An)
    REAL*8, INTENT(in)  :: U(:), bn
    REAL*8, INTENT(out) :: An(:, :)
    REAL*8              :: Unn
    INTEGER             :: ik, ign, inn
    An = 0.d0
    ik = phys%idx_k_eq
    ign = phys%idx_gamman_eq
    inn = phys%idx_rhon_eq
    IF (switch%decoup) THEN
      An(1, 2) = 1.

      An(2, 1) = (-U(2)**2/U(1)**2 + phys%Mref)
      An(2, 2) = 2.*U(2)/U(1)

      An(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      An(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      An(3, 3) = 5./3.*U(2)/U(1)

      An(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      An(4, 2) = 5./3.*U(4)/U(1)
      An(4, 4) = 5./3.*U(2)/U(1)
#ifdef KEQUATION
      if (ik > 0) then
        An(ik, 1) = -U(ik)*U(2)/U(1)**2
        An(ik, 2) = U(ik)/U(1)
        An(ik, ik) = U(2)/U(1)
      end if
#endif
#ifdef NEUTRAL
#ifdef NEUTRALGAMMA
      if (inn > 0 .and. ign > 0) then
        Unn = MAX(1.d-7,U(inn))
        An(inn, ign) = 1.d0
        An(ign, 1) = 2.d0/3.d0*Unn*(-U(3)/U(1)**2 + U(2)**2/U(1)**3)
        An(ign, 2) = -2.d0/3.d0*Unn*U(2)/U(1)**2
        An(ign, 3) = 2.d0/3.d0*Unn/U(1)
        An(ign, inn) = -U(ign)**2/Unn**2 + 2.d0/3.d0*(U(3)/U(1) - 0.5d0*U(2)**2/U(1)**2)
        An(ign, ign) = 2.d0*U(ign)/Unn
      end if
#endif
#endif
    ELSE
      An(1, 2) = 1.

      An(2, 1) = -2./3.*U(2)**2/U(1)**2
      An(2, 2) = 4./3.*U(2)/U(1)
      An(2, 3) = 2./3.
      An(2, 4) = 2./3.

      An(3, 1) = -5./3.*U(2)*U(3)/U(1)**2 + 2./3.*U(2)**3/U(1)**3
      An(3, 2) = 5./3.*U(3)/U(1) - U(2)**2/U(1)**2
      An(3, 3) = 5./3.*U(2)/U(1)

      An(4, 1) = -5./3.*U(4)*U(2)/U(1)**2
      An(4, 2) = 5./3.*U(4)/U(1)
      An(4, 4) = 5./3.*U(2)/U(1)
#ifdef KEQUATION
      if (ik > 0) then
        An(ik, 1) = -U(ik)*U(2)/U(1)**2
        An(ik, 2) = U(ik)/U(1)
        An(ik, ik) = U(2)/U(1)
      end if
#endif
#ifdef NEUTRAL
#ifdef NEUTRALGAMMA
      if (inn > 0 .and. ign > 0) then
        Unn = MAX(1.d-7,U(inn))
        An(inn, ign) = 1.d0
        An(ign, 1) = 2.d0/3.d0*Unn*(-U(3)/U(1)**2 + U(2)**2/U(1)**3)
        An(ign, 2) = -2.d0/3.d0*Unn*U(2)/U(1)**2
        An(ign, 3) = 2.d0/3.d0*Unn/U(1)
        An(ign, inn) = -U(ign)**2/Unn**2 + 2.d0/3.d0*(U(3)/U(1) - 0.5d0*U(2)**2/U(1)**2)
        An(ign, ign) = 2.d0*U(ign)/Unn
      end if
#endif
#endif
    ENDIF
    An = bn*An
  ENDSUBROUTINE jacobianMatricesFace

  !*****************************************
  ! Jacobian matrix for the Bohm BC
  !****************************************
  SUBROUTINE jacobianMatricesBohm(U, A)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: A(:, :)
    REAL*8              :: auxi, auxe
    INTEGER             :: inn

    A = 0.
    inn = phys%idx_rhon_eq
    auxi = (5.-2.*phys%Gmbohm)/3.
    auxe = (5.-2.*phys%Gmbohme)/3.
    A(1, 1) = -auxi*(U(2)**3/U(1)**3 - U(2)*U(3)/U(1)**2)
    A(1, 2) = auxi*(U(3)/U(1) - 3./2.*U(2)**2/U(1)**2)
    A(1, 3) = auxi*U(2)/U(1)

    A(2, 1) = -auxe*U(2)*U(4)/U(1)**2
    A(2, 2) = auxe*U(4)/U(1)
    A(2, 4) = auxe*U(2)/U(1)

  ENDSUBROUTINE jacobianMatricesBohm

#ifdef NEUTRAL
  !*****************************************
  ! Jacobian matrices for Neutrals
  !****************************************
  SUBROUTINE jacobianMatricesN(U, Up, Q, b, sigmaviz, sigmavcx, Ax, Ay)
    REAL*8, INTENT(in)  :: U(:), Up(:), Q(:,:), b(:), sigmaviz, sigmavcx
    REAL*8, INTENT(out) :: Ax(:,:),Ay(:,:)
    REAL*8              :: GradTi(simpar%Ndim), GradTimod, Vnn, Csnn
    INTEGER             :: inn

    inn = phys%idx_rhon_eq

  GradTi(1) = 2./(3*phys%Mref)*( (U(2)**2/U(1)**3 - U(3)/U(1)**2)*Q(1,1) - (U(2)/U(1)**2)*Q(1,2) + 1./U(1)*Q(1,3) )
  GradTi(2) = 2./(3*phys%Mref)*( (U(2)**2/U(1)**3 - U(3)/U(1)**2)*Q(2,1) - (U(2)/U(1)**2)*Q(2,2) + 1./U(1)*Q(2,3) )
  GradTi = simpar%refval_temperature/simpar%refval_length*GradTi

    GradTimod = SQRT(GradTi(1)**2 + GradTi(2)**2)

    Vnn = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTimod/(U(1)*(ABS(sigmaviz) + ABS(sigmavcx)))
    Csnn = SQRT(simpar%refval_charge*simpar%refval_temperature/simpar%refval_mass*Up(7))

  Ax = 0.d0
  !Neutral convective velocity
!  if (Csnn .ge. Vnn) then
!     Ax(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(1)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))
!  else
!     Ax(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(1)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))*Csnn/Vnn
!  endif
!  Ax = Ax/simpar%refval_speed
  !Neutral velocity parallel to magnetic field
  Ax(inn,inn) = Ax(inn,inn) - Up(2)*b(1)

  Ay = 0.d0
!  if (Csnn .ge. Vnn) then
!     Ay(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(2)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))
!  else
!     Ay(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(2)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))*Csnn/Vnn
!  endif
!  Ay = Ay/simpar%refval_speed
  !Neutral velocity parallel to magnetic field
  Ay(inn,inn) = Ay(inn,inn) - Up(2)*b(2)

  Ax = 0.
  Ay = 0.

  ENDSUBROUTINE jacobianMatricesN
#endif

  

  SUBROUTINE add_1D_diff(rho,d_iso,d_ani)
    REAL*8, INTENT(IN)  :: rho(:)
    REAL*8, INTENT(INOUT) :: d_iso(:, :, :), d_ani(:, :, :)
    REAL*8              :: diff_1D(size(rho),4)

    CALL interpolate_1D_diff(rho, diff_1D)

    d_iso(1, 1, :) = d_iso(1, 1, :) + diff_1D(:,1)
    d_iso(2, 2, :) = d_iso(2, 2, :) + diff_1D(:,2)
    d_iso(3, 3, :) = d_iso(3, 3, :) + diff_1D(:,3)
    d_iso(4, 4, :) = d_iso(4, 4, :) + diff_1D(:,4)

    d_ani(1, 1, :) = d_ani(1, 1, :) + diff_1D(:,1)
    d_ani(2, 2, :) = d_ani(2, 2, :) + diff_1D(:,2)
    d_ani(3, 3, :) = d_ani(3, 3, :) + diff_1D(:,3)
    d_ani(4, 4, :) = d_ani(4, 4, :) + diff_1D(:,4)

  ENDSUBROUTINE add_1D_diff

 SUBROUTINE interpolate_1D_diff(rho,diff_1D)
   USE interpolation
   REAL*8, INTENT(IN)  :: rho(:)
   REAL*8, INTENT(INOUT) :: diff_1D(size(rho),4)
   REAL*8            :: rho_loc 
   INTEGER :: i,idx
   
   !linear interpolation of the diffusion coefficients
   DO i=1,size(rho)
       rho_loc = rho(i)
       IF (rho_loc .LT. phys%rho_1D_min) rho_loc = phys%rho_1D_min+1e-10
       IF (rho_loc.GT. phys%rho_1D_max) rho_loc = phys%rho_1D_max-1e-10

       idx = binarySearch(phys%rho_1D_size, phys%rho_1D, rho_loc,1e-12)
       diff_1D(i,1) = phys%diff_n_1D(idx) + (phys%diff_n_1D(idx+1)-phys%diff_n_1D(idx))*(rho_loc-phys%rho_1D(idx))/(phys%rho_1D(idx+1)-phys%rho_1D(idx))
       diff_1D(i,2) = phys%diff_u_1D(idx) + (phys%diff_u_1D(idx+1)-phys%diff_u_1D(idx))*(rho_loc-phys%rho_1D(idx))/(phys%rho_1D(idx+1)-phys%rho_1D(idx))
       diff_1D(i,3) = phys%diff_e_1D(idx) + (phys%diff_e_1D(idx+1)-phys%diff_e_1D(idx))*(rho_loc-phys%rho_1D(idx))/(phys%rho_1D(idx+1)-phys%rho_1D(idx))
       diff_1D(i,4) = phys%diff_ee_1D(idx) + (phys%diff_ee_1D(idx+1)-phys%diff_ee_1D(idx))*(rho_loc-phys%rho_1D(idx))/(phys%rho_1D(idx+1)-phys%rho_1D(idx))

   END DO

  END SUBROUTINE interpolate_1D_diff

  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
#ifndef KEQUATION
  SUBROUTINE setLocalDiff(xy, u, d_iso, d_ani)
#else
  SUBROUTINE setLocalDiff(xy, u, d_iso, d_ani, q_cyl)
#endif
    real*8, intent(in)  		:: xy(:, :)
    real*8, intent(in)  		:: u(:,:)
#ifdef KEQUATION
    real*8, intent(in)  		:: q_cyl(:)
#endif
    real*8, intent(out)		 :: d_iso(:, :, :), d_ani(:, :, :)
    real*8		              :: iperdiff(size(xy, 1))
    integer                     :: inn
#ifdef NEUTRAL
    integer             		:: i
    real*8, dimension(size(u,1))	:: Dnn
#ifdef KEQUATION
    real*8, dimension(size(u,1))          :: D_k,U6,c_s
    real*8                         :: r
    integer                        :: ik
#endif
#endif


    inn = phys%idx_rhon_eq
    ! d_iso(Neq,Neq,Ngauss),d_ani(Neq,Neq,Ngauss)
    ! first index corresponds to the equation
    ! second index corresponds to the unknown
    ! third index correspond to the gauss point
    ! example d_iso(2,3,ig) is the diffusion in the non-diagonal
    ! diffusion term in the second equation/third variable
    ! d_iso>0,d_ani=0 is an isotropic diffusion
    ! d_iso>0,d_ani=d_iso is a perpendicular diffusion
    ! d_iso=0,d_ani<0 is a (positive) parallel diffusion

    d_iso = 0.
    d_ani = 0.
    !*****************************
    ! Diagonal terms
    !*****************************
    d_iso(1, 1, :) = phys%diff_n
    d_iso(2, 2, :) = phys%diff_u
    d_iso(3, 3, :) = phys%diff_e
    d_iso(4, 4, :) = phys%diff_ee
    IF ((switch%ME .EQV. .TRUE.) ) THEN
       IF (switch%testcase .EQ. 85) THEN !Iter core-edge with evolving equilibria plus diffusion decrease
          d_iso(1, 1, :) = phys%diff_n - (phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          d_iso(2, 2, :) = phys%diff_u - (phys%diff_u - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          d_iso(3, 3, :) = phys%diff_e - (phys%diff_e - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
          d_iso(4, 4, :) = phys%diff_ee - (phys%diff_ee - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
       ELSE IF (switch%testcase .EQ. 86) THEN
          d_iso(1,1,:) = MAX(switch%diffmin,phys%diff_n - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
          d_iso(2,2,:) = MAX(switch%diffmin,phys%diff_u - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
          d_iso(3,3,:) = MAX(switch%diffmin,phys%diff_e - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
          d_iso(4,4,:) = MAX(switch%diffmin,phys%diff_ee - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
       ELSE IF (switch%testcase .EQ. 87) THEN
          d_iso(1,1,:) = phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
          d_iso(2,2,:) = phys%diff_u - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
          d_iso(3,3,:) = phys%diff_e - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
          d_iso(4,4,:) = phys%diff_ee - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
       ENDIF

       IF (switch%diff_reverse_Ip) THEN
         d_iso(1,1,:) = phys%diff_n*phys%I_0/phys%I_p
         d_iso(2,2,:) = phys%diff_u*phys%I_0/phys%I_p
         d_iso(3,3,:) = phys%diff_e*phys%I_0/phys%I_p
         d_iso(4,4,:) = phys%diff_ee*phys%I_0/phys%I_p
       ENDIF
      phys%ME_diff_n = d_iso(1,1,1)
      phys%ME_diff_u = d_iso(2,2,1)
      phys%ME_diff_e = d_iso(3,3,1)
      phys%ME_diff_ee = d_iso(4,4,1)
    ENDIF
#ifdef NEUTRAL
#ifdef KEQUATION
    ik = phys%idx_k_eq
    U6 = u(:,ik)
#endif
#ifndef CONSTANTNEUTRALDIFF
    DO i=1,SIZE(u,1)
       CALL compute_Dnn(u(i,:), Dnn(i))
       d_iso(inn,inn,i) = Dnn(i)
    END DO
#else
    d_iso(inn,inn,:)=phys%diff_nn
#endif

#ifdef KEQUATION
    DO i= 1,SIZE(c_s, 1)
       CALL compute_cs(u(i,:), c_s(i))
      ! for all equations
#ifndef DKLINEARIZED
      if (c_s(i)<=1.e-20) then
        D_k(i) = phys%diff_k_min
        !WRITE(6,*) 'NEGATIVE C_S ', c_s(i)
        !stop
      else
#endif
        if ((switch%testcase .ge. 50) .and.(switch%testcase .le. 59)) then
          r = xy(i,1)
        elseif ((switch%testcase .ge. 60) .and.(switch%testcase .le. 69)) then
          r = xy(i,1) + geom%R0/simpar%refval_length
        endif
        D_k(i) = r*U6(i)/c_s(i)
        if (switch%testcase == 60) then
          D_k(i) = D_k(i)*geom%q*2.*PI
        else
          D_k(i) =  D_k(i)*q_cyl(i)*2.*PI
        endif

#ifndef DKLINEARIZED

        D_k(i) = max(phys%diff_k_min,min(phys%diff_k_max,D_k(i) ))
      endif
#else
        !for circular case q_cyl assume constant

          CALL double_softplus(D_k(i),phys%diff_k_min,phys%diff_k_max)
#endif


    enddo
    d_iso(ik,ik,:) = D_k+phys%diff_n
    d_ani(ik,ik,:) = d_iso(ik,ik,:)
    d_iso(1,1,:) = d_iso(1,1,:) + D_k
    d_iso(2,2,:) = d_iso(2,2,:) + D_k
    d_iso(3,3,:) = d_iso(3,3,:) + D_k
    d_iso(4,4,:) = d_iso(4,4,:) + D_k
    !WRITE(6,*) d_iso(6,6,:)*simpar%refval_length**2/simpar%refval_time
    !WRITE(6,*) d_iso(1,1,:)*simpar%refval_length**2/simpar%refval_time
#endif
    !Dnn = sum(Dnn)/size(u,1)
#else
    if (inn > 0) d_iso(inn,inn,:) = 0.
#endif
    d_ani(1, 1, :) = d_iso(1,1,:)
    d_ani(2, 2, :) = d_iso(2,2,:)
    d_ani(3, 3, :) = d_iso(3,3,:)
    d_ani(4, 4, :) = d_iso(4,4,:)
    IF ((switch%ME .EQV. .TRUE.) .AND. (switch%testcase .GT. 84)) THEN !Iter core-edge with evolving equilibria plus diffusion decrease
       IF (switch%testcase .EQ. 85) THEN
         d_ani(1, 1, :) = phys%diff_n - (phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
         d_ani(2, 2, :) = phys%diff_u - (phys%diff_u - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
         d_ani(3, 3, :) = phys%diff_e - (phys%diff_e - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
         d_ani(4, 4, :) = phys%diff_ee - (phys%diff_ee - 0.5*simpar%refval_time/simpar%refval_length**2)/14.65*(phys%I_p - 0.35)
       ELSE IF (switch%testcase .EQ. 86) THEN
          d_ani(1,1,:) = MAX(switch%diffmin,phys%diff_n - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
          d_ani(2,2,:) = MAX(switch%diffmin,phys%diff_u - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
          d_ani(3,3,:) = MAX(switch%diffmin,phys%diff_e - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
          d_ani(4,4,:) = MAX(switch%diffmin,phys%diff_ee - 2.18*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.)))
       ELSE IF (switch%testcase .EQ. 87) THEN
          d_ani(1,1,:) = phys%diff_n - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
          d_ani(2,2,:) = phys%diff_u - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
          d_ani(3,3,:) = phys%diff_e - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
          d_ani(4,4,:) = phys%diff_ee - 0.5*simpar%refval_time/simpar%refval_length**2*(TANH((phys%I_p - 0.35)/5.))
       ENDIF
    ENDIF
#ifdef NEUTRAL
    d_ani(inn,inn,:) = 0.
#endif

    !*****************************
    ! Non diagonal terms
    !*****************************
    ! No non-diagonal terms defined for this model
    CALL computeIperDiffusion(xy,u, iperdiff)
    d_iso(1, 1, :) = d_iso(1, 1, :)+iperdiff
    d_iso(2, 2, :) = d_iso(2, 2, :)+iperdiff
    d_iso(3, 3, :) = d_iso(3, 3, :)+iperdiff
    d_iso(4, 4, :) = d_iso(4, 4, :)+iperdiff
  ENDSUBROUTINE setLocalDiff

  !*******************************************
  ! Compute local diffusion in points
  !*******************************************
  SUBROUTINE computeIperDiffusion(X, u, ipdiff)
    REAL*8, INTENT(IN)  :: X(:, :), u(:,:)
    REAL*8, INTENT(OUT) :: ipdiff(:)
    REAL*8             :: d,dref,maxamp
    REAL*8,ALLOCATABLE :: xcorn(:), ycorn(:)
    REAL*8             :: rad(SIZE(X, 1))
    REAL*8             :: h,rhog
    INTEGER            :: i,g, opt

    ipdiff = 0.

    IF (switch%difcor .GT. 0) THEN

						 SELECT CASE (switch%difcor)
										CASE (1)
												! Circular case with infinitely small limiter
          ALLOCATE(xcorn(1),ycorn(1))
												xcorn = geom%R0
												ycorn = -0.75
										CASE (2)
												! Circular case with infinitely small limiter
          ALLOCATE(xcorn(1),ycorn(1))
												xcorn = geom%R0
												ycorn = -0.287
										CASE (3)
												! West
          ALLOCATE(xcorn(1),ycorn(1))
												xcorn = 2.7977
												ycorn = -0.5128
       CASE(4)
										 ! ITER
          ALLOCATE(xcorn(4),ycorn(4))
												xcorn(1) =  4.023150000000000
												ycorn(1) = -2.544840000000000
												xcorn(2) =  6.23648
												ycorn(2) = -3.23689
												xcorn(3) =  4.02837
												ycorn(3) =  3.588
												xcorn(4) =  5.74627
												ycorn(4) =  4.51401
										CASE DEFAULT
												WRITE (6, *) "Case not valid"
												STOP
						 END SELECT

						 !!**********************************************************
						 !! Gaussian around the corner
						 !!**********************************************************
						 h = 1e-1
       DO i=1,SIZE(xcorn)
          rad = SQRT((X(:, 1)*phys%lscale - xcorn(i))**2 + (X(:, 2)*phys%lscale - ycorn(i))**2)
          ipdiff = ipdiff + numer%dc_coe*EXP(-(2*rad/h)**2)

       END DO
       DO i=1,SIZE(ipdiff)
          IF (ipdiff(i)<1e-5) ipdiff(i)=0.
       END DO
       DEALLOCATE(xcorn,ycorn)
    END IF

    maxamp = 4.
    opt = 2
    IF (switch%limrho .EQ. 2 .AND. MINVAL(u(:,1)) .LT. numer%minrho    ) THEN
       DO g = 1,SIZE(u,1)
          rhog = u(g,1)
          IF (rhog<0.) rhog=0.
          IF (rhog<numer%minrho) THEN
						       d = numer%minrho-rhog ! 0 < d < minrho
             IF (opt.EQ.1) THEN
						          dref = maxamp * d/numer%minrho  ! 0 < dref < maxamp
!			            ipdiff(g) = exp(dref) ! 1 < ipdiff(g) < exp(maxamp)
                ipdiff(g) = ipdiff(g) + EXP(dref) - 1 ! 0 < ipdiff(g) < exp(maxamp)-1
             ELSE IF (opt==2) THEN
						          ipdiff(g) = ipdiff(g) +  1./((1.-d/numer%minrho)*2+1./50.) - 1.
             ENDIF
          ENDIF
       END DO
    ENDIF


  ENDSUBROUTINE computeIperDiffusion

  !*****************************************
  ! Curvature term matrix
  !****************************************
  SUBROUTINE GimpMatrix(U, divb, G)
    REAL*8, INTENT(in)  :: U(:), divb
    REAL*8, INTENT(out) :: G(:, :)

    !G = divb*[              0,                            0,                            0                              0; ...
    !                  1/3*U(2)**2/U(1)**2            -2/3*U(2)/U(1)                      2/3,                           2/3;...
    !                        0                             0                             0                              0; ...
    !                        0                             0                             0                              0];

    G = 0.d0
    IF (switch%decoup) THEN
      G(2, 1) = phys%Mref
    ELSE
      G(2, 1) = 1./3.*U(2)**2/U(1)**2
      G(2, 2) = -2./3.*U(2)/U(1)
      G(2, 3) = 2./3.
      G(2, 4) = 2./3.
    END IF
    G = divb*G
  ENDSUBROUTINE GimpMatrix

#ifdef NEUTRALGAMMA
  SUBROUTINE GimpMatrixN(U, divb, Gn)
    REAL*8, INTENT(in)  :: U(:), divb
    REAL*8, INTENT(out) :: Gn(:, :)
    REAL*8              :: U1, U2, U3, Unn
    INTEGER             :: inn, ign

    Gn = 0.d0
    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq

    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    Unn = MAX(U(inn),1.d-7)

    Gn(ign, 1) = 2.d0/3.d0*Unn*(-U3/U1**2 + U2**2/U1**3)
    Gn(ign, 2) = -2.d0/3.d0*Unn*U2/U1**2
    Gn(ign, 3) = 2.d0/3.d0*Unn/U1
    Gn(ign, inn) = 2.d0/3.d0*(U3/U1 - 0.5d0*U2**2/U1**2)

    Gn = divb*Gn
  ENDSUBROUTINE GimpMatrixN
#endif

  !*****************************************
  ! Parallel diffusion terms
  !****************************************
  SUBROUTINE computeVi(U, V)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: V(:)
    V = 0.d0
    V(1) = U(2)**2/U(1)**3 - U(3)/U(1)**2
    V(2) = -U(2)/U(1)**2
    V(3) = 1./U(1)
  ENDSUBROUTINE computeVi

  SUBROUTINE computeVe(U, V)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: V(:)
    V = 0.d0
    V(1) = -U(4)/U(1)**2
    V(4) = 1./U(1)
  ENDSUBROUTINE computeVe

  !                                SUBROUTINE computeVe(U,V)
  !                                real*8, intent(IN)  :: U(:)
  !                                real*8, intent(OUT) :: V(:)
  !                                V = 0.d0
  !                                V(1) = U(2)**2/U(1)**3 - U(4)/U(1)**2
  !                                V(2) = -U(2)/U(1)**2
  !                                V(4) = 1./U(1)
  !                                ENDSUBROUTINE computeVe

  SUBROUTINE compute_dV_dUi(U, dV_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dV_dU(:, :)
    dV_dU = 0.
    dV_dU(1, 1) = 2*U(3)/U(1)**3 - 3*U(2)**2/U(1)**4
    dV_dU(1, 2) = 2*U(2)/U(1)**3
    dV_dU(1, 3) = -1/U(1)**2
    dV_dU(2, 1) = 2*U(2)/U(1)**3
    dV_dU(2, 2) = -1/U(1)**2
    dV_dU(3, 1) = -1/U(1)**2
  ENDSUBROUTINE compute_dV_dUi

  SUBROUTINE compute_dV_dUe(U, dV_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dV_dU(:, :)
    dV_dU = 0.
    dV_dU(1, 1) = 2*U(4)/U(1)**3
    dV_dU(1, 4) = -1/U(1)**2
    dV_dU(4, 1) = -1/U(1)**2
  ENDSUBROUTINE compute_dV_dUe

#ifdef NEUTRALGAMMA
  SUBROUTINE computeEtan(U,Etan)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Etan
    REAL*8, PARAMETER   :: tol = 1.d-7
    REAL*8              :: coeff, denom, eta_coeff, double_soft_deriv, Unn
    INTEGER             :: inn

    inn = phys%idx_rhon_eq
    Unn = U(inn)
    IF (Unn < tol) Unn = tol

#ifdef CONSTANTNEUTRALDIFF
    Etan = Unn*phys%diff_nn
#else
    CALL compute_neutral_transport_prefactor(U, coeff)
    CALL compute_neutral_gamma_denominator(U, denom)
    eta_coeff = coeff/denom
    CALL double_softplus_deriv(eta_coeff, 10.d0*phys%diff_n, phys%diff_nn, double_soft_deriv)
    CALL double_softplus(eta_coeff, 10.d0*phys%diff_n, phys%diff_nn)
    Etan = Unn*eta_coeff
#endif
  ENDSUBROUTINE computeEtan


  SUBROUTINE compute_dEtan_dU(U,dEtan_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dEtan_dU(:)
    REAL*8, PARAMETER   :: tol = 1.d-7
    REAL*8              :: coeff, denom, eta_coeff, double_soft_deriv, Unn
    REAL*8              :: dcoeff_dU(SIZE(U)), ddenom_dU(SIZE(U)), deta_dU(SIZE(U))
    INTEGER             :: inn

    dEtan_dU = 0.d0
    inn = phys%idx_rhon_eq
    Unn = U(inn)
    IF (Unn < tol) Unn = tol

#ifdef CONSTANTNEUTRALDIFF
    IF (U(inn) >= tol) dEtan_dU(inn) = phys%diff_nn
#else
    CALL compute_neutral_transport_prefactor(U, coeff)
    CALL compute_dneutral_transport_prefactor_dU(U, dcoeff_dU)
    CALL compute_neutral_gamma_denominator(U, denom)
    CALL compute_dneutral_gamma_denominator_dU(U, ddenom_dU)

    eta_coeff = coeff/denom
    CALL double_softplus_deriv(eta_coeff, 10.d0*phys%diff_n, phys%diff_nn, double_soft_deriv)
    deta_dU = (dcoeff_dU/denom - coeff*ddenom_dU/denom**2)*double_soft_deriv
    CALL double_softplus(eta_coeff, 10.d0*phys%diff_n, phys%diff_nn)

    dEtan_dU = Unn*deta_dU
    IF (U(inn) >= tol) dEtan_dU(inn) = dEtan_dU(inn) + eta_coeff
#endif
  ENDSUBROUTINE compute_dEtan_dU


  SUBROUTINE computeVun(U,Vun)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Vun(:)
    REAL*8, PARAMETER   :: tol = 1.d-7
    REAL*8              :: Unn
    INTEGER             :: inn, ign

    Vun = 0.d0
    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq
    Unn = U(inn)
    IF (Unn < tol) Unn = tol

    Vun(inn) = -U(ign)/Unn**2
    Vun(ign) = 1.d0/Unn
  ENDSUBROUTINE computeVun


  SUBROUTINE compute_dVun_dU(U,dVun_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dVun_dU(:, :)
    REAL*8, PARAMETER   :: tol = 1.d-7
    REAL*8              :: Unn
    INTEGER             :: inn, ign

    dVun_dU = 0.d0
    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq
    Unn = U(inn)
    IF (Unn < tol) Unn = tol

    IF (U(inn) >= tol) THEN
      dVun_dU(inn, inn) = 2.d0*U(ign)/Unn**3
      dVun_dU(ign, inn) = -1.d0/Unn**2
    END IF
    dVun_dU(inn, ign) = -1.d0/Unn**2
  ENDSUBROUTINE compute_dVun_dU
#endif

  FUNCTION computeAlphai(U) RESULT(res)
    REAL*8 :: U(:)
    REAL*8 :: res, aux
    REAL*8, PARAMETER :: tol = 1.e-20
    aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2
    IF ((2./(3.*phys%Mref)*aux > phys%T_fluxlim_maxi) .AND. (switch%testcase .NE. 2)) THEN
      res = (3.*phys%Mref/2*phys%T_fluxlim_maxi)**(phys%epn)
    ELSE
       IF (aux<tol) aux = tol
      res = aux**phys%epn
    ENDIF
  END FUNCTION computeAlphai

  FUNCTION computeAlphae(U) RESULT(res)
    REAL*8 :: U(:)
    REAL*8 :: res, aux
    REAL*8, PARAMETER :: tol = 1.e-20
    aux = U(4)/U(1)
    IF ((2./(3.*phys%Mref)*aux > phys%T_fluxlim_maxe) .AND. (switch%testcase .NE. 2)) THEN
      res = (3.*phys%Mref/2*phys%T_fluxlim_maxe)**(phys%epn)
    ELSE
       IF (aux<tol) aux = tol
      res = aux**phys%epn
    ENDIF
  END FUNCTION computeAlphae

  SUBROUTINE compute_dAlpha_dUi(U, res)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: res(:)
    real*8             :: aux
    REAL*8, PARAMETER :: tol = 1.e-20

    aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2

    IF ((2./(3.*phys%Mref)*aux > phys%T_fluxlim_maxi) .AND. (switch%testcase .NE. 2)) THEN !! don't apply flux limiter if it is a convergence test
      res = 0.
    ELSE
       IF (aux<0) aux = tol
      res = 0.d0
      res(1) = -U(3)/U(1)**2+U(2)**2/U(1)**3
      res(2) = -U(2)/U(1)**2
      res(3) = 1./U(1)
      res=phys%epn*aux**(phys%epn-1)*res
    ENDIF
  ENDSUBROUTINE compute_dAlpha_dUi

  SUBROUTINE compute_dAlpha_dUe(U, res)
    real*8, intent(IN) :: U(:)
    real*8, intent(OUT):: res(:)
    real*8             :: aux
    REAL*8, PARAMETER :: tol = 1.e-20

    aux = U(4)/U(1)

    IF ((2./(3.*phys%Mref)*aux > phys%T_fluxlim_maxe) .AND. (switch%testcase .NE. 2)) THEN !! don't apply flux limiter if it is a convergence test
      res = 0.
    ELSE
       IF (aux<0) aux = tol
      res = 0.d0
      res(1) = -U(4)/U(1)**2
      res(4) = 1./U(1)
      res=phys%epn*aux**(phys%epn-1)*res
    ENDIF
  ENDSUBROUTINE compute_dAlpha_dUe

  SUBROUTINE compute_free_streaming_heat_flux_electrons(U, qfs)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: qfs
    REAL*8, PARAMETER :: tmin = 1e-5 ! eV    
    REAL*8            :: tmin_cons
    REAL*8, PARAMETER :: mi = 3.35e-27         ! Ionic mass [kg]
    REAL*8, PARAMETER :: me = 9.109e-31        ! Electronic mass [kg]
    REAL*8             :: t, n

    tmin_cons = tmin/simpar%refval_temperature*3*phys%Mref/2.
    

    t = MAX(U(4)/U(1), tmin_cons)
    n = U(1)
    

    qfs = (2./3.)**1.5*SQRT(mi/me)*t**1.5*n


  END SUBROUTINE

  SUBROUTINE compute_dfree_streaming_heat_flux_electrons_dU(U, res)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: res(:)
    REAL*8, PARAMETER :: mi = 3.35e-27         ! Ionic mass [kg]
    REAL*8, PARAMETER :: me = 9.109e-31        ! Electronic mass [kg]
    REAL*8, PARAMETER :: tmin = 1e-5 ! eV    
    REAL*8            :: tmin_cons

    res = 0.d0

    tmin_cons = tmin/simpar%refval_temperature*3*phys%Mref/2.
    
    IF (U(4)/U(1) > tmin_cons) THEN
      res(1) = -0.5/U(1)
      res(4) = 3./2./U(4)
      res = (2./3.*U(4)/U(1))**1.5*U(1)*SQRT(mi/me)*res
    ENDIF
  END SUBROUTINE

  SUBROUTINE compute_free_streaming_heat_flux_ions(U,qfs)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: qfs
    REAL*8, PARAMETER :: tmin = 1e-5 ! eV    
    REAL*8            :: tmin_cons
    REAL*8             :: t, n

    tmin_cons = tmin/simpar%refval_temperature*3*phys%Mref/2.

    t = MAX(U(3)/U(1) - 0.5*U(2)**2/U(1)**2, tmin_cons)
    n = U(1)



    qfs = (2./3.)**1.5*t**1.5*n


  END SUBROUTINE

  SUBROUTINE compute_dfree_streaming_heat_flux_ions_dU(U, res)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: res(:)
    REAL*8, PARAMETER :: tmin = 1e-5 ! eV    
    REAL*8            :: tmin_cons

    res = 0.d0
    
    tmin_cons = tmin/simpar%refval_temperature*3*phys%Mref/2.
    
    IF (U(3)/U(1) - 0.5*U(2)**2/U(1)**2  > tmin_cons) THEN

      res(1) = -0.5*(U(3)-2.*U(2)**2/U(1))/(U(1)*(U(3)-1./2.*U(2)**2/U(1)))
      res(2) = -3./2.*U(2)/(U(1)*(U(3)-1./2.*U(2)**2/U(1)))
      res(3) = 3./2./(U(3)-1./2.*U(2)**2/U(1))
      res = (2./3.*(U(3)-1./2.*U(2)**2/U(1))/U(1))**1.5*U(1)*res
    ENDIF
  END SUBROUTINE

  SUBROUTINE compute_spitzer_harm_flux_ions(U, Q, b, q_sh_i)
    REAL*8, INTENT(IN) :: U(:), Q(:,:), b(:)
    REAL*8, INTENT(OUT):: q_sh_i
    REAL*8             :: Alphai, gmi, coefi
    REAL*8             :: Vveci(SIZE(U))

    coefi = phys%diff_pari*(2./(3.*phys%Mref))**(1 + phys%epn)
    Alphai = computeAlphai(U)
    CALL computeVi(U, Vveci)
    gmi = dot_PRODUCT(MATMUL(Q,Vveci),b)
    q_sh_i = -coefi*Alphai*gmi
  
  END SUBROUTINE compute_spitzer_harm_flux_ions

  SUBROUTINE compute_spitzer_harm_flux_electrons(U, Q, b, q_sh_e)
    REAL*8, INTENT(IN) :: U(:), Q(:,:), b(:)
    REAL*8, INTENT(OUT):: q_sh_e
    REAL*8             :: Alphae, gme, coefe
    REAL*8             :: Vvece(SIZE(U))

    coefe = phys%diff_pare*(2./(3.*phys%Mref))**(1 + phys%epn)
    Alphae = computeAlphae(U)
    CALL computeVe(U, Vvece)
    gme = dot_PRODUCT(MATMUL(Q,Vvece),b)
    q_sh_e = -coefe*Alphae*gme

  END SUBROUTINE compute_spitzer_harm_flux_electrons

  SUBROUTINE compute_flux_limiter(qfs,qsh,c_fl, flux_limiter)
    REAL*8, INTENT(IN)  :: qfs, qsh, c_fl
    REAL*8, INTENT(OUT) :: flux_limiter
    REAL*8, PARAMETER :: tol = 1.e-10
    REAL*8 :: qfs_clamped
    

    qfs_clamped = MAX(ABS(qfs), tol)

    flux_limiter = 1.0 / (1.0 + ABS(qsh) / (c_fl * qfs_clamped))


  END SUBROUTINE compute_flux_limiter
   

  ! ******************************
  ! Parallel electric field terms
  ! ******************************
  SUBROUTINE compute_W(U, W)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: W

    W = 2./3.*U(2)/U(1)
  ENDSUBROUTINE compute_W

  SUBROUTINE compute_dW_dU(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:, :)
    res = 0.
    res(4, 1) = -2./3.*U(2)/U(1)**2
    res(4, 2) = 2./3./U(1)
  ENDSUBROUTINE compute_dW_dU

  ! ******************************
  ! Temperature exchange terms
  ! ******************************
  SUBROUTINE compute_s(U, s)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: s, U1, U4, U3
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U4 = U(4)
    U3 = U(3)
    IF (U4 < tol) U4 = tol
    IF (U1 < tol) U1 = tol
    IF (U3 < tol) U3 = tol
    ! keeping this thing for high diffusion, but take care for low values
#ifndef KEQUATION
    if ((phys%diff_ee .gt. 0.0380) .and. (switch%testcase .ne. 2) .and. (switch%psdtime)) then
      s = 1./(phys%tie*0.0380/phys%diff_ee)*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
#else
    if (((phys%diff_ee+phys%diff_k_min) .gt. 0.0380) .and. (switch%testcase .ne. 2) .and. (switch%psdtime)) then
      s = 1./(phys%tie*0.0380/(phys%diff_ee+phys%diff_k_min))*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
#endif
    else
      s = 1./(phys%tie)*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
    ENDIF
    !s = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4 - U3 + 0.5*(U(2)**2/U1))
  ENDSUBROUTINE compute_s

  SUBROUTINE compute_ds_dU(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), U1, U4, U3
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U4 = U(4)
    U3 = U(3)
    IF (U4 < tol) U4 = tol
    IF (U1 < tol) U1 = tol
    IF (U3 < tol) U3 = tol
    res = 0.
    res(1) = 2.5*(U1/U4)**1.5*(U4 - U3 + 0.5*(U(2)**2/U1)) - 0.5*U1**0.5*U(2)**2/U4**1.5
    res(2) = U(2)*(U1/U4)**1.5
    res(3) = -U1**2.5/U4**1.5
    res(4) = -1.5*(U1/U4)**2.5*(U4 - U3 + 0.5*U(2)**2/U1) + U1**2.5/U4**1.5
    ! keeping this thing for high diffusion, but take care for low values
#ifndef KEQUATION
    if ((phys%diff_ee .gt. 0.0380) .and. (switch%testcase .ne. 2) .and. (switch%psdtime)) then
      res = 1./(phys%tie*0.0380/phys%diff_ee)*(2./3./phys%Mref)**(-0.5)*res
#else
    if (((phys%diff_ee+phys%diff_k_min) .gt. 0.0380) .and. (switch%testcase .ne. 2) .and. (switch%psdtime)) then
      res = 1./(phys%tie*0.0380/(phys%diff_ee+phys%diff_k_min))*(2./3./phys%Mref)**(-0.5)*res

#endif
    else
     res = 1./(phys%tie)*(2./3./phys%Mref)**(-0.5)*res
    ENDIF
    !res = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*res
  ENDSUBROUTINE compute_ds_dU


  !*****************************
  !Ohmic Heating Source
  !*****************************
  SUBROUTINE compute_Sohmic(U,Sohmic)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: Sohmic,U1,U4
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U4 = U(4)
    IF (U4<tol) U4=tol
    IF (U1<tol) U1=tol
    Sohmic = phys%Zeff*phys%Pohmic*(((3*phys%Mref)/2)**1.5)*((U1/U4)**1.5)
  ENDSUBROUTINE compute_Sohmic


  SUBROUTINE compute_dSohmic_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U4 = U(4)
    IF (U4<tol) U4=tol
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = (1.5*U1**0.5)/(U4**1.5)
    res(4) = -(1.5*U1**1.5)/(U4**2.5)
    res = phys%Zeff*phys%Pohmic*(((3*phys%Mref)/2)**1.5)*res
  ENDSUBROUTINE compute_dSohmic_dU


  ! ******************************
  ! Neutral Source terms
  ! ******************************
#ifdef NEUTRAL
  SUBROUTINE compute_RN(E,theta,RN)
    ! Compute the recycling coefficeint RN(E,theta) interpolating the TRIM data
    USE interpolation
    real*8, intent(IN)   :: E,theta
    real*8, intent(OUT)  :: RN
    REAL*8               :: E_clipped, theta_clipped
    integer              :: ip, jp
  
    RN = 1.
  
    ip = size(phys%E)
    jp = size(phys%theta)

    E_clipped = max(1e-10,min(1e3-1e-10,E))
    theta_clipped = max(1e-10,min(90-1e-10,theta))
  
    RN = interpolate(ip, phys%E, jp, phys%theta, phys%RN_DW, E_clipped, theta_clipped, 1e-12)
  
  END SUBROUTINE compute_RN

  SUBROUTINE compute_niz(U,niz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: niz,U1,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U1 = U(1)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    niz = U1*U5
  ENDSUBROUTINE compute_niz


  SUBROUTINE compute_dniz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U1 = U(1)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.
    res(1) = U5
    res(inn) = U1
  ENDSUBROUTINE compute_dniz_dU


  SUBROUTINE compute_nrec(U,nrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: nrec,U1
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    IF (U1<tol) U1=tol
    nrec = U1**2
  ENDSUBROUTINE compute_nrec


  SUBROUTINE compute_dnrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = 2.*U1
  ENDSUBROUTINE compute_dnrec_dU


  SUBROUTINE compute_fGammacx(U,fGammacx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fGammacx,U2,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U2 = U(2)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    IF (U5<tol) U5=tol
    fGammacx = U2*U5
  ENDSUBROUTINE compute_fGammacx


  SUBROUTINE compute_dfGammacx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U2,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U2 = U(2)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    IF (U5<tol) U5=tol
    res = 0.
    res(2) = U5
    res(inn) = U2
  ENDSUBROUTINE compute_dfGammacx_dU


  SUBROUTINE compute_fGammarec(U,fGammarec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fGammarec,U1,U2
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U2 = U(2)
    IF (U1<tol) U1=tol
    fGammarec = U1*U2
  ENDSUBROUTINE compute_fGammarec


  SUBROUTINE compute_dfGammarec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U2
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U2 = U(2)
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = U2
    res(2) = U1
  ENDSUBROUTINE compute_dfGammarec_dU

#ifdef NEUTRALGAMMA
  SUBROUTINE compute_fGammaN(U,fGammaN)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: fGammaN
    INTEGER             :: ign

    fGammaN = 0.d0
    ign = phys%idx_gamman_eq
    fGammaN = U(1)*U(ign)
  ENDSUBROUTINE compute_fGammaN


  SUBROUTINE compute_dfGammaN_dU(U,res)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: res(:)
    INTEGER             :: ign

    res = 0.d0
    ign = phys%idx_gamman_eq

    res(1) = U(ign)
    res(ign) = U(1)
  ENDSUBROUTINE compute_dfGammaN_dU
#endif


#ifdef TEMPERATURE
#ifndef AMJUELSPLINES
  SUBROUTINE compute_sigmaviz(U,sigmaviz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmaviz, te, e0

    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
    ELSE
      te = neutral_rt%iz_te_floor
    ENDIF

    te = MAX(te, neutral_rt%iz_te_floor)
    e0 = te/neutral_rt%rydberg_energy
    sigmaviz = neutral_rt%rate_scale*1.d-11*SQRT(e0)/(neutral_rydberg_energy_phys**1.5d0*(6.d0 + e0))*EXP(-1.d0/e0)
  ENDSUBROUTINE compute_sigmaviz


  SUBROUTINE compute_dsigmaviz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), te, e0, sigmaviz, dlograte_dte
    REAL*8             :: dte_dU(size(U))

    res = 0.d0
    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
      IF (te>neutral_rt%iz_te_floor) THEN
        CALL compute_dTe_dU(U, dte_dU)
        CALL compute_sigmaviz(U, sigmaviz)
        e0 = te/neutral_rt%rydberg_energy
        dlograte_dte = (0.5d0/e0 - 1.d0/(6.d0 + e0) + 1.d0/e0**2)/neutral_rt%rydberg_energy
        res = sigmaviz*dlograte_dte*dte_dU
      END IF
    END IF
  ENDSUBROUTINE compute_dsigmaviz_dU


  SUBROUTINE compute_sigmavrec(U,sigmavrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavrec, te, e0, g

    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
    ELSE
      te = neutral_rt%rec_te_floor
    ENDIF

    te = MAX(te, neutral_rt%rec_te_floor)
    e0 = neutral_rt%rydberg_energy/te
    g = 0.43d0 + 0.5d0*LOG(e0) + 0.469d0*e0**(-1.d0/3.d0)
    sigmavrec = neutral_rt%rate_scale*5.2d-20*SQRT(e0)*g
  ENDSUBROUTINE compute_sigmavrec


  SUBROUTINE compute_dsigmavrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), te, e0, g, dg_de0, sigmavrec, dlograte_dte
    REAL*8             :: dte_dU(size(U))

    res = 0.d0
    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
      IF (te>neutral_rt%rec_te_floor) THEN
        CALL compute_dTe_dU(U, dte_dU)
        CALL compute_sigmavrec(U, sigmavrec)
        e0 = neutral_rt%rydberg_energy/te
        g = 0.43d0 + 0.5d0*LOG(e0) + 0.469d0*e0**(-1.d0/3.d0)
        dg_de0 = 0.5d0/e0 - 0.469d0/3.d0*e0**(-4.d0/3.d0)
        dlograte_dte = -(0.5d0/e0 + dg_de0/g)*e0/te
        res = sigmavrec*dlograte_dte*dte_dU
      END IF
    END IF
  ENDSUBROUTINE compute_dsigmavrec_dU

  SUBROUTINE compute_sigmavEiz(U,sigmavEiz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavEiz, sigmaviz, Tloss

    CALL compute_sigmaviz(U, sigmaviz)
    CALL compute_Tloss(U, Tloss)
    sigmavEiz = sigmaviz*Tloss
  ENDSUBROUTINE compute_sigmavEiz

  SUBROUTINE compute_dsigmavEiz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), sigmaviz, Tloss
    REAL*8             :: dsigmaviz_dU(size(U)), dTloss_dU(size(U))

    CALL compute_sigmaviz(U, sigmaviz)
    CALL compute_dsigmaviz_dU(U, dsigmaviz_dU)
    CALL compute_Tloss(U, Tloss)
    CALL compute_dTloss_dU(U, dTloss_dU)
    res = dsigmaviz_dU*Tloss + sigmaviz*dTloss_dU
  ENDSUBROUTINE compute_dsigmavEiz_dU

  SUBROUTINE compute_sigmavErec(U,sigmavErec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavErec, sigmavrec, Tlossrec

    CALL compute_sigmavrec(U, sigmavrec)
    CALL compute_Tlossrec(U, Tlossrec)
    sigmavErec = sigmavrec*Tlossrec
  ENDSUBROUTINE compute_sigmavErec

  SUBROUTINE compute_dsigmavErec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), sigmavrec, Tlossrec
    REAL*8             :: dsigmavrec_dU(size(U)), dTlossrec_dU(size(U))

    CALL compute_sigmavrec(U, sigmavrec)
    CALL compute_dsigmavrec_dU(U, dsigmavrec_dU)
    CALL compute_Tlossrec(U, Tlossrec)
    CALL compute_dTlossrec_dU(U, dTlossrec_dU)
    res = dsigmavrec_dU*Tlossrec + sigmavrec*dTlossrec_dU
  ENDSUBROUTINE compute_dsigmavErec_dU
#else


  ! Routines for 2D splines in te, ne loglogspace
  SUBROUTINE compute_2D_eirene_rate(te,ne,alpha,rate)
    ! this routine calculates eirene rate in te, ne space
    ! if in the region of applicability, then just takes the values according to the splines
    ! if ne>ne_max (LTE limit) or ne<ne_min (Corona limit) then no more dependancy on ne
    ! if te<te_min or te>te_max then extrapolates in log space taking the derivative of the edge and linearly expanding
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    REAL*8             :: dlograte_dlogte

    rate = 0.
    ! region 1, where ne is applicable
    IF ((ne>=neutral_rt%eirene_ne_min) .AND. (ne<=neutral_rt%eirene_ne_max)) THEN
       IF ((te>=neutral_rt%eirene_te_min) .AND.(te<=neutral_rt%eirene_te_max)) THEN
          CALL compute_2D_logeirene_rate(te,ne,alpha,rate)
       ELSEIF (te<neutral_rt%eirene_te_min) THEN
          CALL compute_2D_logeirene_rate(neutral_rt%eirene_te_min,ne,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_min,ne,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(neutral_rt%eirene_te_min))
       ELSEIF (te>neutral_rt%eirene_te_max) THEN
          CALL compute_2D_logeirene_rate(neutral_rt%eirene_te_max,ne,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_max,ne,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(neutral_rt%eirene_te_max))
       ENDIF
    ! beyond range of ne applicability
    ELSEIF(ne<neutral_rt%eirene_ne_min) THEN
      ! if te is still applicable
       IF ((te>=neutral_rt%eirene_te_min) .AND.(te<=neutral_rt%eirene_te_max)) THEN
          CALL compute_2D_logeirene_rate(te,neutral_rt%eirene_ne_min,alpha,rate)
      ! if te < te_min
       ELSEIF (te<neutral_rt%eirene_te_min) THEN
          CALL compute_2D_logeirene_rate(neutral_rt%eirene_te_min,neutral_rt%eirene_ne_min,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_min,neutral_rt%eirene_ne_min,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(neutral_rt%eirene_te_min))
      ! if te>te_max
       ELSEIF (te>neutral_rt%eirene_te_max) THEN
          CALL compute_2D_logeirene_rate(neutral_rt%eirene_te_max,neutral_rt%eirene_ne_min,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_max,neutral_rt%eirene_ne_min,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(neutral_rt%eirene_te_max))
       ENDIF
    ! if ne>ne_max
    ELSEIF (ne>neutral_rt%eirene_ne_max) THEN
       IF ((te>=neutral_rt%eirene_te_min) .AND.(te<=neutral_rt%eirene_te_max)) THEN
          CALL compute_2D_logeirene_rate(te,neutral_rt%eirene_ne_max,alpha,rate)
       ELSEIF (te<neutral_rt%eirene_te_min) THEN
          CALL compute_2D_logeirene_rate(neutral_rt%eirene_te_min,neutral_rt%eirene_ne_max,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_min,neutral_rt%eirene_ne_max,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(neutral_rt%eirene_te_min))
       ELSEIF (te>neutral_rt%eirene_te_max) THEN
          CALL compute_2D_logeirene_rate(neutral_rt%eirene_te_max,neutral_rt%eirene_ne_max,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_max,neutral_rt%eirene_ne_max,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(neutral_rt%eirene_te_max))
       ENDIF
    ENDIF
    ! rate is in cm^3/s in EIRENE

    IF (rate < -100) THEN
        rate = 0.0
    ELSE
        rate = EXP(rate)/1.e6
    END IF
  ENDSUBROUTINE compute_2D_eirene_rate

  SUBROUTINE compute_2D_eirene_rate_du(U1,U4,te,ne,alpha,rate_du)
    REAL*8, INTENT(IN) :: U1,U4,te,ne,alpha(:,:)
    REAL*8             :: rate
    REAL*8             :: dlograte_dlogne, dlograte_dlogte
    REAL*8, INTENT(OUT):: rate_du(:)

    CALL compute_2D_eirene_rate(te,ne,alpha,rate)
    rate_du = 0.
    dlograte_dlogne = 0.
    dlograte_dlogte = 0.
    IF ((ne>=neutral_rt%eirene_ne_min) .AND. (ne<=neutral_rt%eirene_ne_max)) THEN
       IF ((te>=neutral_rt%eirene_te_min) .AND.(te<=neutral_rt%eirene_te_max)) THEN
          CALL compute_dlogeirene_2D_dlogne_rate(te,ne,alpha,dlograte_dlogne)
          CALL compute_dlogeirene_2D_dlogte_rate(te,ne,alpha,dlograte_dlogte)
       ELSEIF (te<neutral_rt%eirene_te_min) THEN
          CALL compute_dlogeirene_2D_dlogne_rate(neutral_rt%eirene_te_min,ne,alpha,dlograte_dlogne)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_min,ne,alpha,dlograte_dlogte)
       ELSEIF (te>neutral_rt%eirene_te_max) THEN
          CALL compute_dlogeirene_2D_dlogne_rate(neutral_rt%eirene_te_max,ne,alpha,dlograte_dlogne)
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_max,ne,alpha,dlograte_dlogte)
       ENDIF
    ! beyond range of ne applicability (ne derivative is now zero)
    ELSEIF(ne<neutral_rt%eirene_ne_min) THEN
      ! if te is still applicable
       IF ((te>=neutral_rt%eirene_te_min) .AND.(te<=neutral_rt%eirene_te_max)) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te,neutral_rt%eirene_ne_min,alpha,dlograte_dlogte)
      ! if te < te_min
       ELSEIF (te<neutral_rt%eirene_te_min) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_min,neutral_rt%eirene_ne_min,alpha,dlograte_dlogte)
      ! if te>te_max
       ELSEIF (te>neutral_rt%eirene_te_max) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_max,neutral_rt%eirene_ne_min,alpha,dlograte_dlogte)
       ENDIF
    ! if ne>ne_max (ne derivative is now zero)
    ELSEIF (ne>neutral_rt%eirene_ne_max) THEN
       IF ((te>=neutral_rt%eirene_te_min) .AND.(te<=neutral_rt%eirene_te_max)) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te,neutral_rt%eirene_ne_max,alpha,dlograte_dlogte)
       ELSEIF (te<neutral_rt%eirene_te_min) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_min,neutral_rt%eirene_ne_max,alpha,dlograte_dlogte)
       ELSEIF (te>neutral_rt%eirene_te_max) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(neutral_rt%eirene_te_max,neutral_rt%eirene_ne_max,alpha,dlograte_dlogte)
       ENDIF
    ENDIF
    rate_du(1) = rate_du(1) + dlograte_dlogte*(-1./U1)
    rate_du(1) = rate_du(1) + dlograte_dlogne*(1./U1)
    rate_du(4) = rate_du(4) + dlograte_dlogte*(1./U4)
    rate_du = rate_du*rate
  ENDSUBROUTINE compute_2D_eirene_rate_du

  SUBROUTINE compute_2D_logeirene_rate(te,ne,alpha,rate)
    ! this routine calculate log (eirene_rate) for given te, ne in log log space
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    REAL*8             :: logte, logne  
    INTEGER            :: i,j
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    logte = LOG(te)
    logne = LOG(ne)
    DO j=1,SIZE(alpha,2)
       DO i = 1,SIZE(alpha,1)
          rate = rate + alpha(i,j)*logne**(j-1)*logte**(i-1)
       END DO
    END DO
  ENDSUBROUTINE compute_2D_logeirene_rate

  SUBROUTINE compute_dlogeirene_2D_dlogte_rate(te,ne,alpha,rate)
    ! this routines calculate derivative dlog (eirene_rate)/dlog(te) for given te, ne in log log space
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    REAL*8             :: logte, logne 
    INTEGER            :: i,j
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    logte = LOG(te)
    logne = LOG(ne)
    DO j=1,SIZE(alpha,2)
       DO i = 2,SIZE(alpha,1)
          rate = rate + alpha(i,j)*(i-1)*logne**(j-1)*logte**(i-2)
       END DO
    END DO
  ENDSUBROUTINE compute_dlogeirene_2D_dlogte_rate

  SUBROUTINE compute_dlogeirene_2D_dlogne_rate(te,ne,alpha,rate)
    ! this routine calculate derivative dlog (eirene_rate)/dlog(ne) for given te, ne in log log space
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    REAL*8             :: logte, logne
    INTEGER            :: i,j
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    logte = LOG(te)
    logne = LOG(ne)
    DO j=2,SIZE(alpha,2)
       DO i = 1,SIZE(alpha,1)
          rate = rate + alpha(i,j)*(j-1)*logne**(j-2)*logte**(i-1)
       END DO
    END DO
  ENDSUBROUTINE compute_dlogeirene_2D_dlogne_rate

  SUBROUTINE compute_sigmaviz(U,sigmaviz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmaviz,U1,U4,te,ne

    U1 = U(1)
    U4 = U(4)

    IF ((U1>neutral_rt%state_tol) .AND. (U4>neutral_rt%state_tol)) THEN ! basically it's a below zero check
      CALL compute_Te(U, te)
      ne = U1
    ELSE!some low values
      ne = neutral_ne_floor
      te = neutral_rt%te_floor
    ENDIF

    sigmaviz = 0.

    CALL compute_2D_eirene_rate(te,ne,phys%alpha_iz,sigmaviz)
  ENDSUBROUTINE compute_sigmaviz

  SUBROUTINE compute_dsigmaviz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,te,ne

    U1 = U(1)
    U4 = U(4)

    res = 0.
    IF ((U1>neutral_rt%state_tol) .AND. (U4>neutral_rt%state_tol)) THEN ! basically it's a below zero check
      ne = U1
      CALL compute_Te(U, te)
       CALL compute_2D_eirene_rate_du(U1,U4,te,ne,phys%alpha_iz,res)
    ENDIF !let non-linear part as zero if negative solutions

  ENDSUBROUTINE compute_dsigmaviz_dU

  SUBROUTINE compute_sigmavEiz(U,sigmavEiz)
    real*8, intent(IN) :: U(:)
    real*8             :: sigmavEiz,U1,U4,te,ne
    U1 = U(1)
    U4 = U(4)

    if ((U1>neutral_rt%state_tol) .and. (U4>neutral_rt%state_tol)) then ! basically it's a below zero check
      CALL compute_Te(U, te)
      ne = U1
    else!some low values
      ne = neutral_ne_floor
      te = neutral_rt%te_floor
    endif

    sigmavEiz = 0.

    call compute_2D_eirene_rate(te,ne,phys%alpha_energy_iz,sigmavEiz)
  ENDSUBROUTINE compute_sigmavEiz

  SUBROUTINE compute_dsigmavEiz_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4,te,ne
    U1 = U(1)
    U4 = U(4)

    res = 0.
    if ((U1>neutral_rt%state_tol) .and. (U4>neutral_rt%state_tol)) then ! basically it's a below zero check
      ne = U1
      CALL compute_Te(U, te)
      call compute_2D_eirene_rate_du(U1,U4,te,ne,phys%alpha_energy_iz,res)
    endif !let non-linear part as zero if negative solutions

  ENDSUBROUTINE compute_dsigmavEiz_dU

  ! Neutral-neutral collision reaction rate
  SUBROUTINE compute_sigmavnn(U,sigmavnn)
    real*8, intent(IN) :: U(:)
    real*8             :: sigmavnn,ti

    CALL compute_Ti(U, ti)
    if (ti .LT. neutral_rt%ti_floor) then ! basically it's a below zero check
      !some low values
      ti = neutral_rt%ti_floor
    endif

    sigmavnn = neutral_rt%sigmavnn_prefactor*ti**0.25d0
  ENDSUBROUTINE compute_sigmavnn

  SUBROUTINE compute_dsigmavnn_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),ti
    real*8             :: dti_dU(size(U))

    res = 0.
    CALL compute_Ti(U, ti)
    if (ti>neutral_rt%ti_floor) then ! basically it's a below zero check
      CALL compute_dTi_dU(U, dti_dU)
      res = (0.25d0*neutral_rt%sigmavnn_prefactor/ti**0.75d0) * dti_dU
    endif !let non-linear part as zero if negative solutions
  ENDSUBROUTINE compute_dsigmavnn_dU
  SUBROUTINE compute_sigmavrec(U,sigmavrec)
    real*8, intent(IN) :: U(:)
    real*8             :: sigmavrec,U1,U4,te,ne
    U1 = U(1)
    U4 = U(4)

    if ((U1>neutral_rt%state_tol) .and. (U4>neutral_rt%state_tol)) then ! basically it's a below zero check
      CALL compute_Te(U, te)
      ne = U1
    else!some low values
      ne = neutral_ne_floor
      te = neutral_rt%te_floor
    endif
    call compute_2D_eirene_rate(te,ne,phys%alpha_rec,sigmavrec)
  ENDSUBROUTINE compute_sigmavrec

  SUBROUTINE compute_dsigmavrec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4,te,ne
    U1 = U(1)
    U4 = U(4)

    res = 0.
    if ((U1>neutral_rt%state_tol) .and. (U4>neutral_rt%state_tol)) then ! basically it's a below zero check
      ne = U1
      CALL compute_Te(U, te)
      call compute_2D_eirene_rate_du(U1,U4,te,ne,phys%alpha_rec,res)
    endif !let non-linear part as zero if negative solutions
  ENDSUBROUTINE compute_dsigmavrec_dU

  SUBROUTINE compute_sigmavErec(U,sigmavErec)
    real*8, intent(IN) :: U(:)
    real*8             :: sigmavErec,U1,U4,te,ne
    U1 = U(1)
    U4 = U(4)

    if ((U1>neutral_rt%state_tol) .and. (U4>neutral_rt%state_tol)) then ! basically it's a below zero check
      CALL compute_Te(U, te)
      ne = U1
    else!some low values
      ne = neutral_ne_floor
      te = neutral_rt%te_floor
    endif
    call compute_2D_eirene_rate(te,ne,phys%alpha_energy_rec,sigmavErec)
  ENDSUBROUTINE compute_sigmavErec

  SUBROUTINE compute_dsigmavErec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U4,te,ne
    U1 = U(1)
    U4 = U(4)

    res = 0.
    if ((U1>neutral_rt%state_tol) .and. (U4>neutral_rt%state_tol)) then ! basically it's a below zero check
      ne = U1
      CALL compute_Te(U, te)
      call compute_2D_eirene_rate_du(U1,U4,te,ne,phys%alpha_energy_rec,res)
    endif !let non-linear part as zero if negative solutions
  ENDSUBROUTINE compute_dsigmavErec_dU
#endif
#ifdef MANUELCX
!ADAS truncated CX
  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavcx, te, logte_dim

    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
    ELSE
      te = neutral_rt%cx_te_floor
    ENDIF

    te = MAX(te, neutral_rt%cx_te_floor)
    logte_dim = LOG(te) + neutral_rt%log_temperature_ref
    sigmavcx = neutral_rt%rate_scale*EXP(neutral_manuelcx_coeffs(1)*logte_dim**4 + neutral_manuelcx_coeffs(2)*logte_dim**3 + &
      &neutral_manuelcx_coeffs(3)*logte_dim**2 + neutral_manuelcx_coeffs(4)*logte_dim + neutral_manuelcx_coeffs(5))
  ENDSUBROUTINE compute_sigmavcx


  SUBROUTINE compute_dsigmavcx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), te, logte_dim, sigmavcx, dlograte_dte
    REAL*8             :: dte_dU(size(U))

    res = 0.d0
    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
      IF (te>neutral_rt%cx_te_floor) THEN
        CALL compute_dTe_dU(U, dte_dU)
        CALL compute_sigmavcx(U, sigmavcx)
        logte_dim = LOG(te) + neutral_rt%log_temperature_ref
        dlograte_dte = (4.d0*neutral_manuelcx_coeffs(1)*logte_dim**3 + 3.d0*neutral_manuelcx_coeffs(2)*logte_dim**2 + &
          &2.d0*neutral_manuelcx_coeffs(3)*logte_dim + neutral_manuelcx_coeffs(4))/te
        res = sigmavcx*dlograte_dte*dte_dU
      END IF
    END IF
  ENDSUBROUTINE compute_dsigmavcx_dU
#endif
#ifdef LEGACYCX
  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavcx, te, e0

    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
    ELSE
      te = neutral_rt%te_floor
    ENDIF

    te = MAX(te, neutral_rt%te_floor)
    e0 = 0.5d0/(simpar%refval_temperature*te)
    sigmavcx = neutral_rt%legacy_cx_prefactor*EXP(-e0)
  ENDSUBROUTINE compute_sigmavcx

  SUBROUTINE compute_dsigmavcx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), te, e0, sigmavcx, dlograte_dte
    REAL*8             :: dte_dU(size(U))

    res = 0.d0
    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
      IF (te>neutral_rt%te_floor) THEN
        CALL compute_dTe_dU(U, dte_dU)
        CALL compute_sigmavcx(U, sigmavcx)
        e0 = 0.5d0/(simpar%refval_temperature*te)
        dlograte_dte = e0/te
        res = sigmavcx*dlograte_dte*dte_dU
      END IF
    END IF
  ENDSUBROUTINE compute_dsigmavcx_dU
#endif
#ifdef EXPANDEDCX
! These routines use AMUJUEL splines
  SUBROUTINE compute_eirene_1D_rate(t,alpha,rate)
    ! This routine calculates extrapolated AMJUEL 1D rate (here on ion temperature) for given temperature and coefficients
    real*8, intent(IN) :: t,alpha(:)
    real*8, intent(OUT):: rate
    real*8             :: dlograte_dlogt
    rate = 0.
    if ((t>=neutral_rt%eirene_ti_min) .AND. (t<=neutral_rt%eirene_ti_max)) then
      call compute_logeirene_1D_rate(t,alpha,rate)
    elseif(t<neutral_rt%eirene_ti_min) then
      call compute_logeirene_1D_rate(neutral_rt%eirene_ti_min,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(neutral_rt%eirene_ti_min,alpha,dlograte_dlogt)
      rate = rate + dlograte_dlogt*(log(t)- log(neutral_rt%eirene_ti_min))
    else
      call compute_logeirene_1D_rate(neutral_rt%eirene_ti_max,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(neutral_rt%eirene_ti_max,alpha,dlograte_dlogt)
      rate = rate + dlograte_dlogt*(log(t)- log(neutral_rt%eirene_ti_max))
    endif
    ! rates are not higher than 1 m^3/s, if rate is higher than that value, then there is something weird
    if (rate>6.*log(10.)) then
      WRITE(6,*) "Something weird in compute_eirene_rate, probably, solution is not good already"
      WRITE(6,*) " t equal to", t
      WRITE(6,*) " rate equal to", rate
      stop
    endif
    if (rate < -100) then
      rate = 0.0
    else
      rate = exp(rate)/1.e6
    endif
  ENDSUBROUTINE compute_eirene_1D_rate

  SUBROUTINE compute_eirene_1D_rate_vs_ti_du(ti,dti_dU,alpha,res)
    ! This routine calculates extrapolated AMJUEL 1D rate (typically on temperature) for given temperature and coefficients
    real*8, intent(IN) :: ti,dti_dU(:),alpha(:)
    real*8, intent(OUT):: res(:)
    real*8             :: dlograte_dlogte,rate
    res = 0.

    if ((ti>=neutral_rt%eirene_ti_min) .AND. (ti<=neutral_rt%eirene_ti_max)) then
      call compute_eirene_1D_rate(ti,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(ti,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dti_dU(1)/ti)
      res(2) = res(2) + dlograte_dlogte*(dti_dU(2)/ti)
      res(3) = res(3) + dlograte_dlogte*(dti_dU(3)/ti)
      res = rate*res
    elseif(ti<neutral_rt%eirene_ti_min) then
      call compute_eirene_1D_rate(ti,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(neutral_rt%eirene_ti_min,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dti_dU(1)/ti)
      res(2) = res(2) + dlograte_dlogte*(dti_dU(2)/ti)
      res(3) = res(3) + dlograte_dlogte*(dti_dU(3)/ti)
      res = rate*res
    else
      call compute_eirene_1D_rate(ti,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(neutral_rt%eirene_ti_max,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dti_dU(1)/ti)
      res(2) = res(2) + dlograte_dlogte*(dti_dU(2)/ti)
      res(3) = res(3) + dlograte_dlogte*(dti_dU(3)/ti)
      res = rate*res
    endif

  ENDSUBROUTINE compute_eirene_1D_rate_vs_ti_du

  SUBROUTINE compute_eirene_1D_rate_vs_te_du(te,dte_dU,alpha,res)
    ! This routine calculates extrapolated AMJUEL 1D rate (typically on temperature) for given temperature and coefficients
    real*8, intent(IN) :: te,dte_dU(:),alpha(:)
    real*8, intent(OUT):: res(:)
    real*8             :: dlograte_dlogte,rate
    res = 0.

    if ((te>=neutral_rt%eirene_ti_min) .AND. (te<=neutral_rt%eirene_ti_max)) then
      call compute_eirene_1D_rate(te,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(te,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dte_dU(1)/te)
      res(4) = res(4) + dlograte_dlogte*(dte_dU(4)/te)
      res = rate*res
    elseif(te<neutral_rt%eirene_ti_min) then
      call compute_eirene_1D_rate(te,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(neutral_rt%eirene_ti_min,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dte_dU(1)/te)
      res(4) = res(4) + dlograte_dlogte*(dte_dU(4)/te)
      res = rate*res
    else
      call compute_eirene_1D_rate(te,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogt(neutral_rt%eirene_ti_max,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dte_dU(1)/te)
      res(4) = res(4) + dlograte_dlogte*(dte_dU(4)/te)
      res = rate*res
    endif

  ENDSUBROUTINE compute_eirene_1D_rate_vs_te_du


  SUBROUTINE compute_logeirene_1D_rate(ti,alpha,rate)
    ! Calculates 1D AMJUEL spline in loglog space
    real*8, intent(IN) :: ti,alpha(:)
    real*8, intent(OUT):: rate
    integer            :: i
    rate = 0.

    do i = 1,size(alpha,1)
      rate = rate + alpha(i)*log(ti)**(i-1)
    end do

  ENDSUBROUTINE compute_logeirene_1D_rate


  SUBROUTINE compute_d_logeirene_1D_rate_dlogt(t,alpha,d_log_rate_dt)
    ! calculates derivative of AMJUEL 1D spline in loglog space
    real*8, intent(IN) :: t, alpha(:)
    real*8, intent(OUT):: d_log_rate_dt
    integer            :: i
    d_log_rate_dt = 0.

    do i = 2,size(alpha,1)
      d_log_rate_dt = d_log_rate_dt + (i-1)*alpha(i)*log(t)**(i-2)
    end do
  ENDSUBROUTINE compute_d_logeirene_1D_rate_dlogt
  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    ! calculates AMJUEL CX rate
    real*8, intent(IN)  :: U(:)
    real*8              :: sigmavcx,ti

    CALL compute_Ti(U, ti)
    if (ti<neutral_rt%ti_floor) then ! basically it's a below zero check
      !some low values
      ti = neutral_rt%ti_floor
    endif
    sigmavcx = 0.

    call compute_eirene_1D_rate(ti, phys%alpha_cx, sigmavcx)
  ENDSUBROUTINE compute_sigmavcx


  SUBROUTINE compute_dsigmavcx_dU(U,res)
    ! calculates derivative of AMJUEL CX rate for linearization
    real*8, intent(IN) :: U(:)
    real*8             :: res(:), ti
    real*8             :: dti_dU(size(U))
    res = 0.
    CALL compute_Ti(U, ti)

    if (ti>neutral_rt%ti_floor) then
      CALL compute_dTi_dU(U, dti_dU)
      call compute_eirene_1D_rate_vs_ti_dU(ti,dti_dU,phys%alpha_cx,res)
    endif

  ENDSUBROUTINE compute_dsigmavcx_dU

  SUBROUTINE compute_cooling_factor(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res,U1,U4,te

    U1 = U(1)
    U4 = U(4)
    res = 0.
    IF ((U1>neutral_rt%state_tol) .AND. (U4>neutral_rt%state_tol)) THEN ! basically it's a below zero check
      CALL compute_Te(U, te)
      CALL compute_eirene_1D_rate(te,phys%alpha_cooling_factor,res)

    ENDIF
  ENDSUBROUTINE compute_cooling_factor

  SUBROUTINE compute_dcooling_factor_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,te
    REAL*8             :: dte_dU(size(U))

    U1 = U(1)
    U4 = U(4)
    res = 0.
    IF ((U1>neutral_rt%state_tol) .AND. (U4>neutral_rt%state_tol)) THEN ! basically it's a below zero check
      CALL compute_Te(U, te)
      CALL compute_dTe_dU(U, dte_dU)
      CALL compute_eirene_1D_rate_vs_te_du(te,dte_dU,phys%alpha_cooling_factor,res)
      res = res!/simpar%refval_charge
    ENDIF
  ENDSUBROUTINE compute_dcooling_factor_dU


#endif

  SUBROUTINE compute_Dnn_dU(U, Dnn_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: Dnn_dU(:)
    REAL*8              :: Dnn, double_soft_deriv, denom
    REAL*8              :: dcoeff_dU(size(U))
    REAL*8              :: ddenom_dU(size(U))

    Dnn_dU = 0.d0
#ifndef CONSTANTNEUTRALDIFF
    CALL compute_neutral_transport_prefactor(U, Dnn)
    CALL compute_dneutral_transport_prefactor_dU(U, dcoeff_dU)
    CALL compute_neutral_diffusion_denominator(U, denom)
    CALL compute_dneutral_diffusion_denominator_dU(U, ddenom_dU)

    Dnn = Dnn/denom
    CALL double_softplus_deriv(Dnn, phys%diff_nn_min, phys%diff_nn, double_soft_deriv)

    Dnn_dU = dcoeff_dU/denom - (Dnn*denom)*ddenom_dU/denom**2
    Dnn_dU = Dnn_dU*double_soft_deriv
#endif
  ENDSUBROUTINE compute_Dnn_dU

  SUBROUTINE compute_W5p(U, W5p)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: W5p(:)
    REAL*8              :: Dnn, Ti_limited, alpha, supp, ti_factor
    INTEGER             :: inn

    ti_factor = 2.d0/(3.d0*phys%Mref)
    inn = phys%idx_rhon_eq
    CALL compute_Dnn(U, Dnn)
    CALL compute_limited_Ti(U, Ti_limited)
    supp = Ti_limited/(Ti_limited + neutral_rt%transport_ti_supp)

    alpha = numer%neutralp_lambda*ti_factor*U(inn)*Dnn/Ti_limited
    CALL computeVi(U, W5p)
    W5p = alpha*supp*W5p
  ENDSUBROUTINE compute_W5p

  SUBROUTINE compute_dW5p_dU(U, dW5p_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dW5p_dU(:, :)
    REAL*8              :: Vi(size(U)), dVi_dU(size(U), size(U))
    REAL*8              :: Dnn, Dnn_dU(size(U))
    REAL*8              :: Ti_limited
    REAL*8              :: dTi_limited_dU(size(U))
    REAL*8              :: alpha, dalpha_dU(size(U)), ti_factor
    REAL*8              :: supp, dsupp_dU(size(U))
    INTEGER             :: inn, j

    dW5p_dU = 0.d0
    ti_factor = 2.d0/(3.d0*phys%Mref)
    inn = phys%idx_rhon_eq

    CALL computeVi(U, Vi)
    CALL compute_dV_dUi(U, dVi_dU)
    CALL compute_Dnn(U, Dnn)
    CALL compute_Dnn_dU(U, Dnn_dU)
    CALL compute_limited_Ti(U, Ti_limited)
    CALL compute_dlimited_Ti_dU(U, dTi_limited_dU)
    supp = Ti_limited/(Ti_limited + neutral_rt%transport_ti_supp)
    dsupp_dU = neutral_rt%transport_ti_supp*dTi_limited_dU/(Ti_limited + neutral_rt%transport_ti_supp)**2

    alpha = numer%neutralp_lambda*ti_factor*U(inn)*Dnn/Ti_limited
    dalpha_dU = numer%neutralp_lambda*ti_factor*(U(inn)*Dnn_dU/Ti_limited - U(inn)*Dnn*dTi_limited_dU/Ti_limited**2)
    dalpha_dU(inn) = dalpha_dU(inn) + numer%neutralp_lambda*ti_factor*Dnn/Ti_limited

    DO j = 1, SIZE(U)
      dW5p_dU(:,j) = Vi*(supp*dalpha_dU(j) + alpha*dsupp_dU(j)) + alpha*supp*dVi_dU(:,j)
    END DO
  ENDSUBROUTINE compute_dW5p_dU

  SUBROUTINE compute_Tloss(U,Tloss)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: Tloss, te

    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
    ELSE
      te = neutral_rt%te_floor
    ENDIF

    te = MAX(te, neutral_rt%te_floor)
    Tloss = neutral_rt%tloss_offset + neutral_rt%tloss_amplitude*EXP(-neutral_rt%tloss_decay*te)
  ENDSUBROUTINE compute_Tloss


  SUBROUTINE compute_dTloss_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), te
    REAL*8             :: dte_dU(size(U))

    res = 0.d0
    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
      IF (te>neutral_rt%te_floor) THEN
        CALL compute_dTe_dU(U, dte_dU)
        res = -neutral_rt%tloss_amplitude*neutral_rt%tloss_decay*EXP(-neutral_rt%tloss_decay*te)*dte_dU
      END IF
    END IF
  ENDSUBROUTINE compute_dTloss_dU


  SUBROUTINE compute_Tlossrec(U,Tlossrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: Tlossrec, te, logTlossrec

    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
    ELSE
      te = neutral_rt%te_floor
    ENDIF

    te = MAX(te, neutral_rt%te_floor)
    logTlossrec = neutral_rt%tlossrec_growth*te
    IF (logTlossrec .LE. neutral_rt%tlossrec_cap_log) THEN
       Tlossrec = neutral_rt%tlossrec_prefactor*EXP(logTlossrec)
    ELSE
       Tlossrec = neutral_rt%tlossrec_cap
    ENDIF
  ENDSUBROUTINE compute_Tlossrec


  SUBROUTINE compute_dTlossrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), te, Tlossrec
    REAL*8             :: dte_dU(size(U))

    res = 0.d0
    IF ((U(1)>neutral_rt%state_tol) .AND. (U(4)>neutral_rt%state_tol)) THEN
      CALL compute_Te(U, te)
      IF (te>neutral_rt%te_floor) THEN
        IF (neutral_rt%tlossrec_growth*te .LE. neutral_rt%tlossrec_cap_log) THEN
          CALL compute_dTe_dU(U, dte_dU)
          CALL compute_Tlossrec(U, Tlossrec)
          res = Tlossrec*neutral_rt%tlossrec_growth*dte_dU
        END IF
      END IF
    END IF
  ENDSUBROUTINE compute_dTlossrec_dU


  SUBROUTINE compute_fEiiz(U,fEiiz)
    real*8, intent(IN) :: U(:)
    real*8             :: fEiiz,U1,U3,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U3 = U(3)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    U1 = U(1)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    if (U5<tol) U5=tol
    !PSI review
    fEiiz = U5*(U3-0.5*U(2)**2/U1)
  ENDSUBROUTINE compute_fEiiz


  SUBROUTINE compute_dfEiiz_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U3,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U1 = U(1)
    U3 = U(3)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    if (U5<tol) U5=tol
    res = 0.
    res(3) = U5
    res(inn) = U3
    !PSI review
    res(1) = res(1)+0.5*U(2)**2/U1**2*U5
    res(2) = -1.*U(2)/U1*U5
    res(inn) = res(inn) - 0.5*U(2)**2/U1
  ENDSUBROUTINE compute_dfEiiz_dU


  SUBROUTINE compute_fEirec(U,fEirec)
    real*8, intent(IN) :: U(:)
    real*8             :: fEirec,U1,U3
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U3 = U(3)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    fEirec = U1*U3

    fEirec = fEirec
  ENDSUBROUTINE compute_fEirec


  SUBROUTINE compute_dfEirec_dU(U,res)
    real*8, intent(IN) :: U(:)
    real*8             :: res(:),U1,U3
    REAL*8, PARAMETER :: tol = 1.e-20
    U1 = U(1)
    U3 = U(3)
    if (U1<tol) U1=tol
    if (U3<tol) U3=tol
    res = 0.
    res(1) = U3
    res(3) = U1
  ENDSUBROUTINE compute_dfEirec_dU


  SUBROUTINE compute_fEicx(U,fEicx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fEicx,U1,U2,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U1 = U(1)
    U2 = U(2)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    fEicx = (U5*U2**2)/U1*0.5
  ENDSUBROUTINE compute_fEicx


  SUBROUTINE compute_dfEicx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U2,U5
    REAL*8, PARAMETER :: tol = 1.e-20
    INTEGER            :: inn
    U1 = U(1)
    U2 = U(2)
    inn = phys%idx_rhon_eq
    U5 = U(inn)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.
    res(1) = -U5*(U2/U1)**2
    res(2) = 2.*U5*U2/U1
    res(inn) = (U2**2)/U1
    res(:) = res(:)*0.5
  ENDSUBROUTINE compute_dfEicx_dU

#ifdef NEUTRALGAMMA
  SUBROUTINE compute_fEiN(U,fEiN)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: fEiN
    REAL*8, PARAMETER   :: tol = 1.d-7
    REAL*8              :: Unn
    INTEGER             :: inn, ign

    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq
    Unn = U(inn)
    IF (Unn < tol) Unn = tol

    fEiN = 0.5d0*U(1)*U(ign)**2/Unn
  ENDSUBROUTINE compute_fEiN


  SUBROUTINE compute_dfEiN_dU(U,res)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: res(:)
    REAL*8, PARAMETER   :: tol = 1.d-7
    REAL*8              :: Unn
    INTEGER             :: inn, ign

    res = 0.d0
    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq
    Unn = U(inn)
    IF (Unn < tol) Unn = tol

    res(1) = U(ign)**2/Unn
    IF (U(inn) >= tol) res(inn) = -U(1)*(U(ign)/Unn)**2
    res(ign) = 2.d0*U(1)*U(ign)/Unn

    res = 0.5d0*res
  ENDSUBROUTINE compute_dfEiN_dU
#endif


#ifdef NEUTRAL
#ifdef KEQUATION
#ifdef DKLINEARIZED
SUBROUTINE compute_ddk_du(U,xy,q_cyl,ddk_du)
    ! Routine that computes linearization of turbulent diffusion
    real*8, intent(IN) :: U(:), xy(:),q_cyl
    real*8, intent(OUT) :: ddk_du(:)
    real*8              :: cs,dk
    real*8              :: dcs_du(size(U,1))
    !softplus stuff
    real*8              :: double_soft_deriv
    REAL*8, PARAMETER :: tol = 1.e-20

    ddk_du(:) = 0.
    !modification softplus dk
    !call compute_cs(U,cs)
    !!if (cs>tol) then
    !  dk = 2.*PI*q_cyl*xy(1)*U(6)/cs
    !  call double_softplus_deriv(dk,phys%diff_k_min,phys%diff_k_max,double_soft_deriv)
    !  !if ((dk>phys%diff_k_min) .and.(dk<phys%diff_k_max)) then
    !    call compute_dcs_du(U,dcs_du)
    !    ddk_du(:) = -1*dk/cs*dcs_du(:)
    !    ddk_du(6) = ddk_du(6) + dk/U(6)
    !    ddk_du(:) = ddk_du(:)*double_soft_deriv
    !  !endif
    !
    !!endif



ENDSUBROUTINE compute_ddk_du
#endif
SUBROUTINE compute_gamma_I(U,Q, Btor, gradBtor, R, gamma_I)
  ! growth rate for turbulent energy
    REAL*8, INTENT(IN) :: U(:), Q(:,:), gradBtor(:), Btor, R
    REAL*8             :: U1,U2,U3,U4, ti, te, gr_p_gr_b, cs, p, theta, ti_te
    REAL*8, INTENT(OUT) :: gamma_I
    REAL*8, PARAMETER :: tol = 1.e-20
  U1 = U(1)
  U2 = U(2)
  U3 = U(3)
  U4 = U(4)
  gamma_I=0.
    CALL compute_cs(U, cs)
    ! grad(pi) x gradB
    p = U3-1./2.*U2**2/U1

    ti = MAX(p/U1,tol)
    te = MAX(U4/U1,tol)
    IF (p<tol) p = tol
    theta = 5.*(1.+ti/te)
    gr_p_gr_b = gradBtor(1)*(Q(1,3)-Q(1,2)*U2/U1+1./2.*U2**2/U1**2*Q(1,1))+gradBtor(2)*(Q(2,3)-Q(2,2)*U2/U1+1./2.*U2**2/U1**2*Q(2,1))
    gr_p_gr_b = gr_p_gr_b/Btor/p-theta/R**2
    IF (gr_p_gr_b >= 0) THEN
       gamma_I = cs*SQRT(gr_p_gr_b)
    ELSE
      !gamma_I = -1.*cs*sqrt(-1.*gr_p_gr_b)
      gamma_I=0.
    ENDIF
ENDSUBROUTINE compute_gamma_I

SUBROUTINE compute_gamma_ke(U, Q, B, gradB, q_cyl, omega, gamma_ke)
  ! growth rate for turbulent energy
  real*8, intent(IN) :: U(:), Q(:, :), gradB(:), B, q_cyl, omega
  logical :: is_core
  real*8             :: n, v, ti, te, V0, nB, nu_e, DB, D_perp, nu_perp, d_star, rho_L, nu_star, L_para, dn_dr, dn_dz, C_Omega, tau_para, tau, C_star, aa, an, a_phi, b_nr, b_phir, b_ni, b_phii, gr, gi
  real*8, intent(OUT) :: gamma_ke
  REAL*8, PARAMETER :: tol = 1.e-20, m_ratio = sqrt(3670.4829678537167), coulomb_log = 15.

  n = max(tol, U(1))
  v = U(2)/n
  Ti = max(tol, 2./3./phys%Mref*(U(3)/n - v**2))
  Te = max(tol, 2./3./phys%Mref*U(4)/n)
  call compute_cs(U, V0)
  V0 = V0*m_ratio
  nB = n

  nu_e = 2.91e-12*n*simpar%refval_density*coulomb_log*(Te*simpar%refval_temperature)**(-1.5)*simpar%refval_time
  DB = Te/abs(B)
  D_perp = 1e-2*DB
  nu_perp = 1e-2*DB
  L_para = PI*q_cyl*geom%R0/simpar%refval_length
  rho_L = V0/omega
  nu_star = L_para/V0*nu_e
  d_star = sqrt((D_perp + nu_perp)/DB)
  dn_dr = Q(1, 1)
  dn_dz = Q(2, 1)

  C_star = V0/Omega/L_para
  C_Omega = m_ratio/nu_star
  is_core = (nu_star<m_ratio)
  C_Omega = merge(C_Omega, min(1., C_Omega), is_core)*C_star

  an = D_perp/DB/d_star + sqrt(C_Omega)
  a_phi = nu_perp/DB/d_star + d_star
  b_nr = rho_L/sqrt(2*d_star)*(dn_dr - dn_dz)/nB
  b_ni = C_Omega**0.75
  b_phir = -sqrt(2*d_star)*rho_L/abs(B)*gradB(1)
  b_phii = merge(d_star*C_Omega**0.75, 0., is_core)

  aa = (an + a_phi)/2
  tau_para = L_para/V0
  tau = tau_para/sqrt(C_Omega)*C_star

  gi = -(b_nr * b_phii + b_ni * b_phir) / C_Omega
  gr = aa**2 - an * a_phi - (b_nr * b_phir - b_ni * b_phii) / C_Omega

  gamma_ke = (sqrt((gr + norm2([gr, gi], dim=1))/2) - aa)/tau

ENDSUBROUTINE compute_gamma_ke

SUBROUTINE compute_ce(U,Q, Btor, gradBtor, r,omega_c,q_cyl, ce)
  ! dissipation rate for turbulent energy
    REAL*8, INTENT(IN) :: U(:), Q(:,:), gradBtor(:), r, Btor,omega_c,q_cyl
    REAL*8             :: U1,U2,U3, gamma_I, rhoL, cs, gamma_e , k_loc
    REAL*8, INTENT(OUT) :: ce
    REAL*8, PARAMETER :: tol = 1.e-20
  U1 = U(1)
  U2 = U(2)
  U3 = U(3)

  gamma_e = 4.5
    CALL compute_cs(U, cs)
    IF (cs < tol) cs = tol
    CALL compute_rhoL(U, r, omega_c, rhoL)
    IF (rhoL < tol) rhoL = tol
    CALL compute_gamma_I(U,Q,Btor,gradBtor,r,gamma_I)

  ce = gamma_I*(PI**2/8./gamma_e/rhoL**2/cs**2+1./phys%k_max)

ENDSUBROUTINE compute_ce
SUBROUTINE compute_rhoL(U, R,omega_c, rhoL)
  ! Larmor radii
    REAL*8, INTENT(IN) :: U(:), R,omega_c
    REAL*8             :: cs
    REAL*8, INTENT(OUT) :: rhoL
    CALL compute_cs(U, cs)

  rhoL = cs/omega_c/R

ENDSUBROUTINE compute_rhoL

SUBROUTINE compute_dissip(U, dissip)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: U6
    REAL*8             :: dissip
    REAL*8, PARAMETER :: tol = 1.e-10
  U6 = U(phys%idx_k_eq)
  dissip = U6**2
ENDSUBROUTINE  compute_dissip

SUBROUTINE compute_ddissip_du(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: U6
    REAL*8             :: res(:)
    REAL*8, PARAMETER :: tol = 1.e-10
  U6 = U(phys%idx_k_eq)
  !if (U6 < tol) U6 = tol
  res = 0.
  res(phys%idx_k_eq) = 2.*U6
ENDSUBROUTINE  compute_ddissip_du

#endif
#endif

#endif
!TEMPERATURE

#endif
!NEUTRAL

  !***********************************************************************
  !
  !    COMPUTATION OF THE STABILIZATION PARAMETER
  !
  !***********************************************************************

  !*******************************************
  ! Compute the stabilization tensor tau
  !*******************************************
  SUBROUTINE computeTauGaussPoints(up, uc, q, b, n, iel, isext, xy, tau,diff_iso,diff_ani)
    real*8, intent(in)  :: up(:), uc(:), q(:), b(:), n(:), xy(:)
    REAL*8, intent(in)    :: isext
    integer, intent(in) ::  iel
    real*8, intent(out) :: tau(:, :)
#ifdef NEUTRAL
    REAL*8              :: tau_aux(size(uc))
    REAL*8,INTENT(IN)   :: diff_iso(:, :),diff_ani(:, :)
#else
    REAL*8              :: tau_aux(size(uc))
    REAL*8,INTENT(IN)   :: diff_iso(:, :),diff_ani(:, :)
#endif
    integer             :: ndim
    integer             :: ik, inn, ign
    real*8              :: bn, bnorm,xyd(1,size(xy)),uu(1,size(uc)),qq(1,size(q))
    REAL*8              :: q_fs_i, q_fs_e, flux_limiter_i, flux_limiter_e, q_sh_i, q_sh_e
    REAL*8              :: Qpr(simpar%Ndim,simpar%Neq)
#ifdef NEUTRALGAMMA
    REAL*8              :: Etan
#endif

    real*8              :: U1, U2, U3, U4
    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    tau = 0.
    ik = phys%idx_k_eq
    inn = phys%idx_rhon_eq
    ign = phys%idx_gamman_eq
    ndim = SIZE(n)
    bn = dot_PRODUCT(b(1:ndim), n)
    bnorm = NORM2(b(1:ndim))
    xyd(1,:) = xy(:)
    uu(1,:) = uc(:)
    qq(1,:) = q(:)
    Qpr = RESHAPE(q,(/simpar%Ndim,simpar%Neq/))
    ! Compute flux limiters at the face
    IF (switch%flux_limiter) THEN
      call compute_free_streaming_heat_flux_electrons(uc,q_fs_e)
      call compute_free_streaming_heat_flux_ions(uc,q_fs_i)
      call compute_spitzer_harm_flux_electrons(uc,Qpr,b,q_sh_e)
      call compute_spitzer_harm_flux_ions(uc,Qpr,b,q_sh_i)
      call compute_flux_limiter(q_fs_e,q_sh_e,phys%c_fle,flux_limiter_e)
      call compute_flux_limiter(q_fs_i,q_sh_i,phys%c_fli,flux_limiter_i)

    ELSE
      flux_limiter_e = 1.
      flux_limiter_i = 1.
    END IF

#ifdef NEUTRALGAMMA
    IF (ign > 0) CALL computeEtan(uc,Etan)
#endif

    IF (numer%stab == 2) THEN
       IF (ABS(isext - 1.) .LT. 1e-12) THEN
        ! exterior faces
          tau_aux = ABS((4*uc(2)*bn)/uc(1))
       ELSE
          tau_aux = MAX(ABS(5./3.*up(2)*bn), ABS(0.3*bn*(3*uc(1) + SQRT(ABS(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)))
       ENDIF
#ifdef TOR3D
       IF (ABS(n(3)) > 0.1) THEN
        ! Poloidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(3) = tau_aux(3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(4) = tau_aux(4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifdef NEUTRAL
        tau_aux(inn) = tau_aux(inn) + phys%diff_nn*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#endif
       ELSE
#endif
        ! Toroidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(3) = tau_aux(3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(4) = tau_aux(4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
        tau_aux(inn) = tau_aux(inn) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif
#ifdef TOR3D
       ENDIF
#endif

    ELSEIF (numer%stab == 3) THEN
       IF (ABS(isext - 1.) .LT. 1e-12) THEN
        ! exterior faces
          tau_aux = MAX(ABS((5*U2 - 2*U2*phys%Gmbohme)/(3*U1)), ABS((5*U2 - 2*U2*phys%Gmbohm)/(3*U1)))
       ELSE
          tau_aux = MAX(ABS((3*U2*bn + 5**(0.5)*bn*(-U2**2 + 2*U1*U3 + 2*U1*U4)**(0.5))/(3*U1)),&
          &abs((3*U2*bn - 5**(0.5)*bn*(-U2**2 + 2*U1*U3 + 2*U1*U4)**(0.5))/(3*U1)),&
          &abs((U2*bn)/U1),&
          &abs((5*U2*bn)/(3*U1)))
       ENDIF
#ifdef TOR3D
       IF (ABS(n(3)) > 0.1) THEN
        ! Poloidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(3) = tau_aux(3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(4) = tau_aux(4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifdef NEUTRAL
        tau_aux(inn) = tau_aux(inn) + phys%diff_nn*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#endif
       ELSE
#endif
        ! Toroidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
        tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(3) = tau_aux(3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(4) = tau_aux(4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
        tau_aux(inn) = tau_aux(inn) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif
#ifdef TOR3D
       ENDIF
#endif

    ELSEIF (numer%stab == 4) THEN
       tau_aux = MAX(ABS(5./3.*up(2)*bn), ABS(0.3*bn*(3*uc(1) + SQRT(ABS(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)), &
            phys%lscale/geom%R0*ABS(bn)*phys%diff_pari*up(7)**2.5, phys%lscale/geom%R0*ABS(bn)*phys%diff_pare*up(8)**2.5)

    ELSEIF (numer%stab == 5) THEN
       !IF (ABS(isext - 1.) .LT. 1e-12) THEN
       ! ! exterior faces
       !   tau_aux = ABS((4*uc(2)*bn)/uc(1))
       !ELSE
          tau_aux = MAX(ABS(5./3.*up(2)*bn), ABS(0.3*bn*(3*uc(1) + SQRT(ABS(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)))
#ifdef TOR3D
       IF (ABS(n(3)) > 0.1) THEN
        ! Poloidal face
        tau_aux(1) = tau_aux(1) + phys%diff_n
        tau_aux(2) = tau_aux(2) + phys%diff_u
          tau_aux(3) = tau_aux(3) + phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1)*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(4) = tau_aux(4) + phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1)*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifdef NEUTRAL
        tau_aux(inn) = phys%diff_nn!numer%tau(5) !tau_aux(inn) + diff_iso(inn,inn,1)
#endif
       ELSE
#endif
        ! Toroidal face
        tau_aux(1) = tau_aux(1) + diff_iso(1,1)*refElPol%ndeg/Mesh%elemSize(iel)
        tau_aux(2) = tau_aux(2) + diff_iso(2,2)*refElPol%ndeg/Mesh%elemSize(iel)
          tau_aux(3) = tau_aux(3) + diff_iso(3,3)*refElPol%ndeg/Mesh%elemSize(iel) + flux_limiter_i*ABS(bn)*phys%diff_pari*(MIN(phys%T_fluxlim_maxi,up(7)))**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)!/phys%lscale
          tau_aux(4) = tau_aux(4) + diff_iso(4,4)*refElPol%ndeg/Mesh%elemSize(iel) + flux_limiter_e*ABS(bn)*phys%diff_pare*(MIN(phys%T_fluxlim_maxe,up(8)))**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)!/phys%lscale
#ifdef NEUTRAL
        tau_aux(inn) = tau_aux(inn) + diff_iso(inn,inn)*refElPol%ndeg/Mesh%elemSize(iel) !! !numer%tau(5) diff_iso(inn,inn,1)
#ifdef KEQUATION
        if (ik > 0) tau_aux(ik) = tau_aux(ik) + diff_iso(ik,ik)*refElPol%ndeg/Mesh%elemSize(iel)
#endif
#ifdef NEUTRALGAMMA
        if (ign > 0 .and. inn > 0) tau_aux(ign) = MAX(numer%tau(ign), tau_aux(ign) + Etan/MAX(uc(inn),1.d-7))*refElPol%ndeg/Mesh%elemSize(iel)
#endif
#endif
!        ! Toroidal face
!        tau_aux(1) = tau_aux(1) +  diff_iso(1,1,1)
!        tau_aux(2) = tau_aux(2) +  diff_iso(2,2,1)
!        tau_aux(3) = tau_aux(3) +  diff_iso(3,3,1) + abs(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!        tau_aux(4) = tau_aux(4) +  diff_iso(4,4,1) + abs(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
!#ifdef NEUTRAL
!        tau_aux(5) = tau_aux(5) +  diff_iso(5,5,1)
!#endif
#ifdef TOR3D
       ENDIF
#endif

    ELSE
       WRITE (6, *) "Wrong stabilization type: ", numer%stab
       STOP
    ENDIF
    tau(1, 1) = tau_aux(1)
    tau(2, 2) = tau_aux(2)
    tau(3, 3) = tau_aux(3)
    tau(4, 4) = tau_aux(4)
#ifdef NEUTRAL
    tau(inn, inn) = tau_aux(inn)
#ifdef NEUTRALGAMMA
    if (ign > 0) tau(ign, ign) = tau_aux(ign)
#endif
#ifdef KEQUATION
    if (ik > 0) tau(ik,ik) = tau_aux(ik)
#endif
#endif
  ENDSUBROUTINE computeTauGaussPoints

  !!

  SUBROUTINE computeTauGaussPoints_matrix(up, uc, b, n, xy, isext, iel, tau)

    REAL*8, INTENT(in)  :: up(:), uc(:), b(:), n(:), xy(:), isext
    REAL*8, INTENT(out) :: tau(:, :)
    INTEGER, INTENT(in) :: iel
    REAL*8, PARAMETER :: eps = 1e-12
    REAL*8              :: bn, bnorm

    REAL*8 :: U1, U2, U3, U4
    REAL*8 :: t2, t3, t4, t5, t6, t7, t8, t9
    REAL*8 :: t10, t11, t12, t13, t14, t15, t16, t17, t18, t19
    REAL*8 :: t20, t21, t22, t23, t24, t25, t26, t27, t28, t29
    REAL*8 :: t30, t31, t32, t33, t34, t35, t36, t37, t38, t39
    REAL*8 :: t40, t41, t42, t43, t44, t45, t46, t47, t48, t49
    REAL*8 :: t50, t51, t52, t53, t54, t55, t56, t57, t58, t59
    REAL*8 :: t60, t61, t62, t63, t64, t65, t66, t67, t68, t69
    REAL*8 :: t70, t71, t72, t73, t74, t75, t76, t77, t78, t79, t80
    REAL*8 :: x, y

    tau = 0.
    bn = dot_PRODUCT(b, n)
    bnorm = NORM2(b)

    x = xy(1)
    y = xy(2)

    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    !************************************
    !
    ! *****     CONVECTIVE PART  ********
    !
    !************************************
    IF (ABS(isext - 1.) .LT. 1e-12) THEN

      !************************************
      !   EXTERIOR FACES
      !************************************
       tau(3, 1) = (1.0D0/U1**2*ABS(U2*bn)*(U1*U3 - U2**2)*(-4.0D0))/ABS(U1)
       tau(3, 2) = (ABS(U2*bn)*(U1*U3*2.0D0 - U2**2*3.0D0)*2.0D0)/(U1*(U2 + eps)*ABS(U1))
       tau(3, 3) = (ABS(U2*bn)*4.0D0)/ABS(U1)

    ELSE
      !************************************
      !   INTERIOR FACES
      !************************************
       t2 = ABS(U1)
      t3 = U1*U3*2.0D0
      t4 = U1*U4*2.0D0
      t5 = U2**2
      !          t6 = t3+t4-t5
       t6 = ABS(t3 + t4 - t5)
      t7 = U2*bn
       t8 = ABS(t7)
      t9 = t5**2
      t10 = U2*bn*3.0D0
       t11 = SQRT(5.0D0)
       t12 = SQRT(t6)
      t13 = bn*t11*t12
      t14 = t10 + t13
       t15 = ABS(t14)
      t16 = t10 - t13
       t17 = ABS(t16)
      t18 = t6**(3.0D0/2.0D0)
      t19 = 1.0D0/t6
      t20 = 1.0D0/t2
      t21 = U2*3.0D0
      t22 = t11*t12
      t23 = U1**2
      t24 = t21 + t22
      t25 = bn*t24
       t26 = ABS(t25)
      t27 = t21 - t22
      t28 = bn*t27
       t29 = ABS(t28)
      t30 = t8*(-6.0D0) + t26 + t29
      t31 = t19*t20*t23*t30*(1.0D0/5.0D0)
      t32 = 1.0D0/U1
      t33 = U1*U3*1.0D1
      t34 = U1*U4*1.0D1
      t51 = t5*9.0D0
      t35 = t33 + t34 - t51
      t36 = 1.0D0/t35
      t37 = U2*t12*5.0D0
      t38 = U1*U3*t11*1.0D1
      t39 = U1*U4*t11*1.0D1
      t40 = U2*t11*t12
      t41 = t3 + t4
      t42 = t2*t6*1.35D2
      t43 = t42 - t2*t41*6.0D1
      t44 = 1.0D0/t43
      t45 = U1*U2*t5*t8*7.2D1
      t46 = U1*U2*t6*t15*1.5D1
      t47 = U1*U2*t6*t17*1.5D1
      t48 = U1*t11*t15*t18*5.0D0
      t49 = U1*t5*t11*t12*t17*4.0D0
      t50 = t19*t44*(t45+t46+t47+t48+t49-U1*t11*t17*t18*5.0D0-U1*U2*t6*t8*9.0D1-U1*U2*t5*t15*1.2D1-U1*U2*t5*t17*1.2D1-U1*t5*t11*t12*t15*4.0D0)
      t52 = U1*U3*5.0D0
      t53 = U1*U4*5.0D0
      t54 = t5*(-4.0D0) + t52 + t53
      t55 = 1.0D0/U1**2
      t56 = t5*7.0D0
      t57 = 1.0D0/t6**(3.0D0/2.0D0)
      t58 = t33 + t34 - t40 - t56
      t59 = t9*1.5D1
      t60 = U3**2
      t61 = t23*t60*5.0D1
      t62 = U3*U4*t23*5.0D1
      t63 = U2*t5*t11*t12*3.0D0
      t65 = U1*U3*t5*5.5D1
      t66 = U1*U4*t5*3.0D1
      t64 = t59 + t61 + t62 + t63 - t65 - t66
      t67 = U2*2.0D0
      t68 = U1*U3*5.0D1
      t69 = U1*U4*5.0D1
      t73 = t5*4.5D1
      t70 = t68 + t69 - t73
      t71 = 1.0D0/t70
      t72 = t22 + t67
      t74 = t59 + t61 + t62 - t63 - t65 - t66
      t75 = t11*t20*t26*t57*t71*t72*t74*(1.0D0/1.5D1)
      t76 = t11*t20*t29*t57*t64*t71*(t22 - t67)*(1.0D0/1.5D1)
       t77 = 1.0D0/SQRT(t6)
      t78 = U1*U4*t12*t26*5.0D0
      t79 = U1*U4*t12*t29*5.0D0
      t80 = U1*U2*U4*t11*t26*2.0D0
      tau(1,1) = (t19*(t8*t9*2.4D1-t9*t15*4.0D0-t9*t17*4.0D0+t6**2*t8*5.0D1-t5*t6*t8*7.0D1+t5*t6*t15*5.0D0+t5*t6*t17*5.0D0-U2*t11*t15*t18*5.0D0+U2*t11*t17*t18*5.0D0+U2*t5*t11*t12*t15*4.0D0-U2*t5*t11*t12*t17*4.0D0))/(t2*t6*9.0D1-t2*t41*4.0D1)
      tau(1, 2) = t19*t20*(U1*U2*t8*(-1.2D1) + U1*U2*t15*2.0D0 + U1*U2*t17*2.0D0 - U1*t11*t12*t15 + U1*t11*t12*t17)*(-1.0D0/1.0D1)
      tau(1, 3) = t31
      tau(1, 4) = t31
      tau(2,1) = U2*t8*t19*t20*t32*t54*(2.0D0/5.0D0)+U2*t11*t19*t20*t29*t32*t36*t58*(t37-t38-t39+t5*t11*1.1D1)*(1.0D0/1.5D2)-U2*t11*t19*t20*t26*t32*t36*(t37+t38+t39-t5*t11*1.1D1)*(t5*(-7.0D0)+t33+t34+t40)*(1.0D0/1.5D2)
      tau(2,2) = t19*t20*(t5*t8*3.6D1-t5*t15*6.0D0+t6*t15*5.0D0-t5*t17*6.0D0+t6*t17*5.0D0+U2*t11*t12*t15-U2*t11*t12*t17)*(1.0D0/3.0D1)
      tau(2, 3) = t50
      tau(2, 4) = t50
      tau(3,1) = t5*t8*t19*t20*t54*t55*(1.0D0/5.0D0)-U4*t5*t8*t20*t32*t36*(2.5D1/3.0D0)+U2*t11*t20*t29*t36*t55*t57*t58*t64*(1.0D0/1.5D2)-U2*t11*t20*t26*t36*t55*t57*(t33+t34+t40-t56)*(t59+t61+t62-U1*U3*t5*5.5D1-U1*U4*t5*3.0D1-U2*t5*t11*t12*3.0D0)*(1.0D0/1.5D2)
      tau(3,2) = t11*t20*t29*t32*t57*t64*(-1.0D0/1.5D2)+t11*t20*t26*t32*t57*(t59+t61+t62-t63-U1*U3*t5*5.5D1-U1*U4*t5*3.0D1)*(1.0D0/1.5D2)+U2*t5*t8*t19*t20*t32*(3.0D0/5.0D0)
      tau(3, 3) = t75 + t76 - t5*t8*t19*t20*(3.0D0/5.0D0) + U1*U4*t8*t20*t36*(5.0D1/3.0D0)
      tau(3, 4) = t75 + t76 - t5*t8*t19*t20*(3.0D0/5.0D0) - t8*t20*t36*(t33 - t51)*(5.0D0/3.0D0)
      tau(4,1) = (t11*t77*(U2*U4*t5*t15*(-1.0D1)+U2*U4*t6*t15*2.5D1+U2*U4*t5*t17*1.0D1-U2*U4*t6*t17*2.5D1-U4*t5*t8*t11*t12*5.0D1+U4*t5*t11*t12*t15*5.0D0+U4*t5*t11*t12*t17*5.0D0)*(-1.0D0/5.0D0))/(U1*t2*t6*5.4D1-U1*t2*t41*2.4D1)
      tau(4, 2) = t20*t77*(U4*t11*t15 - U4*t11*t17)*(1.0D0/6.0D0)
      tau(4, 3) = t20*t36*t77*(t78 + t79 + t80 - U1*U4*t8*t12*5.0D1 - U1*U2*U4*t11*t29*2.0D0)*(1.0D0/3.0D0)
      tau(4, 4) = t20*t36*t77*(t78 + t79 + t80 - t5*t8*t12*4.5D1 + U1*U3*t8*t12*5.0D1 - U1*U2*U4*t11*t29*2.0D0)*(1.0D0/3.0D0)
    END IF

    !************************************
    !
    ! *****     DIFFUSIVE PART  ********
    !
    !************************************
    tau(1, 1) = tau(1, 1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
    tau(2, 2) = tau(2, 2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
    tau(3, 3) = tau(3, 3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
    tau(4, 4) = tau(4, 4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale

  ENDSUBROUTINE computeTauGaussPoints_matrix

#ifdef TEMPERATURE
  !!!! Routines to apply smoothening on limiting values of neutral diffusion
    SUBROUTINE double_softplus(x, xmin, xmax)
      ! this routine constrains value x between xmin and xmax
      ! using paradigm of softplus function
      ! for xmin it is a typical softplus
      ! f(x) = xmin+width*ln(1+exp((x-xmin)/w)
      ! w here and after = w*xmin(or max), where w is defined inside the function
      ! parameter width states for the region where smoothening is applied xmax+-width*w
      ! for xmax it is somewhat inversed softplus:
      ! f(x) = width*ln(1+exp(xmax/width))-width*ln(1+exp(-(x-xmax)/width))
      ! for x>= xmax+width*w*xmax : f(x)=xmax
      ! for xmax-width*w*xmax<=x<xmax+width*w*xmax : f(x) = w*xmax*ln(1+exp(1/w))-width*w*xmax*ln(1+exp(-(x-xmax)/(w*xmax))
      ! for xmin+width*w*xmin<=x<xmax-width*w*xmax : f(x) = x
      ! for xmin-width*w*xmin<=x<xmin+width*w*xmin : f(x) = xmin + w*xmin*ln(1+exp((x-xmin)/(w*xmin))
      ! x<xmin-width*w*xmin : f(x) = xmin
      REAL*8, INTENT(IN) :: xmin, xmax
      REAL*8, INTENT(INOUT):: x
      REAL*8             :: w,width
      w = 0.01
      width = 10
      IF (x>=xmax+w*width*xmax) THEN
        x = xmax
      ELSEIF ((x>=xmax-w*width*xmax) .AND. (x<xmax+w*width*xmax)) THEN
         x = xmax-w*xmax*LOG(1+EXP(-(x-xmax)/(w*xmax)))
      !elseif ((x>=xmin+w*width*xmin) .and. (x<xmax-w*width*xmax)) then
        ! do nothing
      ELSEIF ((x>=xmin-w*width*xmin) .AND. (x<xmin+w*width*xmin)) THEN
         x = xmin + w*xmin*LOG(1+EXP((x-xmin)/(w*xmin)))
      ELSEIF (x<xmin-w*width*xmin) THEN
        x = xmin
      ENDIF
    ENDSUBROUTINE double_softplus

    SUBROUTINE double_softplus_deriv(x, xmin, xmax,deriv)
      ! this calculates dervitive of double_softplus
      ! for x>= xmax+width*w*xmax : f'(x)=0
      ! for xmax-width*w*xmax<=x<xmax+width*w*xmax : f'(x) = 1/(1+exp((x-xmax)/(w*xmax)))
      ! for xmin+width*w*xmin<=x<xmax-width*w*xmax : f'(x) = 1.
      ! for xmin-width*w*xmin<=x<xmin+width*w*xmin : f'(x) = 1/(1+exp(-(x-xmin)/(w*xmin)))
      ! x<xmin-width*w*xmin : f'(x) = 0
      REAL*8, INTENT(IN) :: x,xmin, xmax
      REAL*8, INTENT(OUT):: deriv
      REAL*8             :: w, width
      w = 0.01
      width = 10
      IF (x>=xmax+w*width*xmax) THEN
        deriv = 0.
      ELSEIF ((x>=xmax-w*width*xmax) .AND. (x<xmax+w*width*xmax)) THEN
         deriv = 1./(1.+EXP((x-xmax)/(w*xmax)))
      ELSEIF ((x>=xmin+w*width*xmin) .AND. (x<xmax-w*width*xmax)) THEN
        deriv = 1.
      ELSEIF ((x>=xmin-w*width*xmin) .AND. (x<xmin+w*width*xmin)) THEN
         deriv = 1./(1.+EXP(-1.*(x-xmin)/(w*xmin)))
        !WRITE(6,*) 'Low diffusion ', x*simpar%refval_diffusion
        !stop

      ELSEIF (x<xmin-w*width*xmin) THEN
        deriv = 0.
      ENDIF
    ENDSUBROUTINE double_softplus_deriv

    SUBROUTINE softplus(x, xmin)
      ! this routine limits value x with xmin
      ! using paradigm of softplus function
      ! f(x) = xmin+width*ln(1+exp((x-xmin)/width)
      ! w here and after = w*xmin(or max), where w is defined inside the function
      ! parameter width states for the region where smoothening is applied xmax+-width*w
      ! for x>=xmin-width*w*xmin : f(x) = xmin + w*xmin*ln(1+exp((x-xmin)/(w*xmin))
      ! x<xmin-width*w*xmin : f(x) = xmin
      REAL*8, INTENT(IN) :: xmin
      REAL*8, INTENT(INOUT):: x
      REAL*8             :: w, width
      w = 0.01
      width = 10
      !if (x>=xmin+w*width.xmin) then
      !  x = x !do nothing
      IF ((x>=xmin-w*width*xmin) .AND. (x<xmin+w*width*xmin)) THEN
         x = xmin + w*xmin*LOG(1+EXP((x-xmin)/(w*xmin)))
      ELSEIF (x<xmin-w*width*xmin) THEN
        x = xmin
      ENDIF
    ENDSUBROUTINE softplus

    SUBROUTINE softplus_deriv(x, xmin,deriv)
      ! this routine calculates derivtiv of softplus
      ! x>=xmin+width*w*xmin: f'(x) = 1.
      ! for xmin-width*w*xmin<=x<xmin+width*w*xmin : f'(x) = 1/(1+exp(-(x-xmin)/(w*xmin)))
      ! x<xmin-width*w*xmin : f'(x) = 0
      REAL*8, INTENT(IN) :: x, xmin
      REAL*8, INTENT(OUT):: deriv
      REAL*8             :: w, width
      w = 0.01
      width = 10
      IF (x>=xmin+w*width*xmin) THEN
        deriv = 1.
      ELSEIF ((x>=xmin-w*width*xmin) .AND. (x<xmin+w*width*xmin)) THEN
         deriv =  1./(1.+EXP(-1.*(x-xmin)/(w*xmin)))
      ELSEIF (x<xmin-w*width*xmin) THEN
        deriv = 0.
      ENDIF
    ENDSUBROUTINE softplus_deriv
    !*******************************************
    ! Compute the terms relative to k equations
    !*******************************************
    SUBROUTINE compute_cs(U, cs)
      ! Sound speed of plasma
      REAL*8, INTENT(IN) :: U(:)
      REAL*8             :: U1,U2,U3,U4
      REAL*8, INTENT(OUT) :: cs
      REAL :: tol
    tol = 1.e-20
      U1 = U(1)
      U2 = U(2)
      U3 = U(3)
      U4 = U(4)
      cs = 2./3./U1*(U3+U4-1./2.*U2**2/U1)
      !modification softplus dk
      !IF (cs<0.) cs = tol**2
      call softplus(cs,tol)
      cs = sqrt(cs)
    ENDSUBROUTINE compute_cs
    SUBROUTINE compute_dcs_du(U, dcs_du)
      ! Sound speed derivative
      real*8, intent(IN) :: U(:)
      real*8             :: U1,U2,U3,U4,cs,cs_real,soft_deriv
      real*8, intent(OUT) :: dcs_du(:)
      REAL :: tol
      tol = 1.e-20
      U1 = U(1)
      U2 = U(2)
      U3 = U(3)
      U4 = U(4)
      !if (U4 < tol) U4 = tol
      !if (U1 < tol) U1 = tol
      !if (U3 < tol) U3 = tol
      dcs_du = 0.
      !modification softplus dk
      cs_real = 2./3./U1*(U3+U4-1./2.*U2**2/U1)
      call softplus_deriv(cs_real,tol,soft_deriv)
      call compute_cs(U,cs)
      !if (cs>tol) then
        dcs_du(1) = -1.*(U3+U4-U2**2/U1)/U1**2
        dcs_du(2) = -1.*U2/U1**2
        dcs_du(3) = 1./U1
        dcs_du(4) = 1./U1


        dcs_du = dcs_du/3./cs*soft_deriv
      !endif
    ENDSUBROUTINE compute_dcs_du
#endif

END MODULE physics
