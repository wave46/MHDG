!*****************************************
! project: MHDG
! file: physics.f90
! date: 20/09/2017
! Define the physics of the model
!  ******** N-Gamma-Ti-Te system     ****
!*****************************************
MODULE physics

  USE globals
  USE printUtils
  USE magnetic_field
  IMPLICIT NONE

CONTAINS

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE initPhys()

    ! number of equation of the problem
    phys%Neq = 4
#ifdef NEUTRAL
    phys%Neq = 5
#ifdef KEQUATION
    !so far we only use k equation with neutrals and the convention is that the k-equation is always the last
    phys%Neq = 6
#endif
#endif

    ! number of physical variables
    phys%npv = 10
#ifdef NEUTRAL
    phys%npv = 11
#ifdef KEQUATION
    phys%npv = 12
#endif
#endif

    ALLOCATE (phys%phyVarNam(phys%npv))
    ALLOCATE (phys%conVarNam(phys%Neq))

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
#ifdef NEUTRAL
    phys%phyVarNam(11)= "rhon"   ! density neutral
#ifdef KEQUATION
    phys%phyVarNam(12)= "k"   ! turbulent energy
#endif
#endif

    ! Set the name of the conservative variables
    phys%conVarNam(1) = "rho"   ! U1 = rho
    phys%conVarNam(2) = "Gamma" ! U2 = rho*u
    phys%conVarNam(3) = "nEi"   ! U3 = rho*Ei
    phys%conVarNam(4) = "nEe"   ! U4 = rho*Ee
#ifdef NEUTRAL
    phys%conVarNam(5) = "rhon"  ! U5 = rhon
#ifdef KEQUATION
    phys%conVarNam(6) = "k"  ! U6 = k
#endif
#endif

    simpar%model = 'N-Gamma-Ti-Te'
#ifdef NEUTRAL
    simpar%model = 'N-Gamma-Ti-Te-Neutral'
#ifdef KEQUATION
    simpar%model = 'N-Gamma-Ti-Te-Neutral-k'
#endif
#endif
    simpar%Ndim = 2
#ifdef TOR3D
    simpar%Ndim = 3
#endif
    simpar%Neq = phys%Neq

    ALLOCATE (simpar%physvar_refval(phys%npv))
    ALLOCATE (simpar%consvar_refval(phys%Neq))
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
#ifdef NEUTRAL
    simpar%physvar_refval(11) = simpar%refval_neutral
#ifdef KEQUATION
    simpar%physvar_refval(12) = simpar%refval_k
#endif
#endif
    simpar%consvar_refval(1) = simpar%refval_density
    simpar%consvar_refval(2) = simpar%refval_momentum
    simpar%consvar_refval(3) = simpar%refval_specenergydens
    simpar%consvar_refval(4) = simpar%refval_specenergydens
#ifdef NEUTRAL
    simpar%consvar_refval(5) = simpar%refval_neutral
#ifdef KEQUATION
    simpar%consvar_refval(6) = simpar%refval_k
#endif
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
  END SUBROUTINE initPhys

  !*******************************************
  ! Convert physical variable to conservative
  ! variables
  !*******************************************
  SUBROUTINE phys2cons(up, ua)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: up
    REAL*8, DIMENSION(:, :), INTENT(out) :: ua

    ua(:, 1) = ABS(up(:, 1))
    ua(:, 2) = up(:, 1)*up(:, 2)
    ua(:, 3) = up(:, 1)*up(:, 3)
    ua(:, 4) = up(:, 1)*up(:, 4)
#ifdef NEUTRAL
    ua(:,5) = ABS(up(:,11))
#ifdef KEQUATION
    ua(:,6) = ABS(up(:,12))
#endif
#endif

  END SUBROUTINE phys2cons

  !*******************************************
  ! Convert conservative variable to physical
  ! variables
  !*******************************************
  SUBROUTINE cons2phys(ua, up)
    REAL*8, DIMENSION(:, :), INTENT(in)  :: ua
    REAL*8, DIMENSION(:, :), INTENT(out) :: up
    REAL*8,  DIMENSION(SIZE(ua,1))       :: U1

    U1 = ABS(ua(:,1))

    up(:, 1) = ABS(U1)                                                           ! density
    up(:, 2) = ua(:, 2)/U1                                            ! u parallel
    up(:, 3) = ua(:, 3)/U1                                            ! total energy of ions
    up(:, 4) = ua(:, 4)/U1                                            ! total energy of electrons
    up(:, 5) = ABS(2./(3.*phys%Mref)*(ua(:, 3) - 0.5*ua(:, 2)**2/U1)) ! pressure of ions
    up(:, 6) = ABS(2./(3.*phys%Mref)*ua(:, 4))                              ! pressure of electrons
    up(:, 7) = up(:, 5)/U1                                            ! temperature of ions
    up(:, 8) = up(:, 6)/U1                                            ! temperature of electrons
    up(:, 9) = SQRT((up(:, 7) + up(:, 8))*phys%Mref)                        ! sound speed
    up(:, 10) = up(:, 2)/up(:, 9)                                           ! Mach
#ifdef NEUTRAL
    up(:,11) = ABS(ua(:,5))                                                 ! density neutral
#ifdef KEQUATION
    up(:,12) = ABS(ua(:,6))                                                ! turbulent energy
#endif
#endif


  END SUBROUTINE cons2phys

  ! ******************************
  ! Split diffusion terms
  ! ******************************
  SUBROUTINE compute_W2(U,W2,diff_n,diff_u)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_u
    REAL*8             :: W2(:)
    W2 = 0.
    W2(1) = (diff_n-diff_u)*U(2)/U(1)
  END SUBROUTINE compute_W2

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
  END SUBROUTINE compute_W3

  SUBROUTINE compute_W4(U,W4,diff_n,diff_ee)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_ee
    REAL*8             :: W4(:)

    W4 = 0.
    W4(1) = (diff_n-diff_ee)*U(4)/U(1)
  END SUBROUTINE compute_W4



  SUBROUTINE compute_dW2_dU(U,dW2_dU,diff_n,diff_u)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_u
    REAL*8             :: dW2_dU(:,:)
    dW2_dU = 0.
    dW2_dU(1,1) = -U(2)/(U(1)**2)
    dW2_dU(1,2) = 1./U(1)

    dW2_dU = (diff_n-diff_u)*dW2_dU
  END SUBROUTINE compute_dW2_dU

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
  END SUBROUTINE compute_dW3_dU

  SUBROUTINE compute_dW4_dU(U,res,diff_n,diff_ee)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(IN) :: diff_n, diff_ee
    REAL*8             :: res(:,:)
    res = 0.
    res(1,1) = -U(4)*(diff_n-diff_ee)/(U(1)**2)
    res(1,4) = 1.*(diff_n-diff_ee)/U(1)
  END SUBROUTINE compute_dW4_dU


  !*****************************************
  ! Jacobian matrices
  !****************************************
  SUBROUTINE jacobianMatrices(U, A)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: A(:, :)
    ![ 0,                                           1,                              0,                0; ...
    ! -2/3*U(2)**2/U(1)**2                          4/3*U(2)/U(1)                   2/3,              2/3; ...
    ! -5/3*U(2)*U(3)/U(1)**2+2/3*U(2)**3/U(1)**3    5/3*U(3)/U(1)-U(2)**2/U(1)**2   5/3*U(2)/U(1),    0 ;   ...
    ! -5/3*U(4)*U(2)/U(1)**2,                       5/3*U(4)/U(1),                  0,                5/3*U(2)/U(1)]
    ! k equation line
    ! -U(6)*U(2)/U(1)**2,                           U(6)/U(1),                      0,                0,            0,        U(2)/U(1)]
    A = 0.d0
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
       A(6, 1) = -U(6)*U(2)/U(1)**2
       A(6, 2) = U(6)/U(1)
       A(6, 6) = U(2)/U(1)
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
       A(6, 1) = -U(6)*U(2)/U(1)**2
       A(6, 2) = U(6)/U(1)
       !IF (U(6) >= 0.) THEN
       A(6, 6) = U(2)/U(1)
       !ENDIF
#endif
#ifdef NEUTRAL
#ifdef NEUTRALCONVECTION
       A(5, 1) = -U(5)*U(2)/U(1)**2
       A(5, 2) = U(5)/U(1)
       A(5, 5) = U(2)/U(1)
       !A(5, :) = simpar%refval_time/(simpar%refval_length**2*phys%diff_n)*A(5,:)
#endif
#endif
    END IF
  END SUBROUTINE jacobianMatrices

#ifdef NEUTRALP
  SUBROUTINE jacobianMatricesNP(U, Anp)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: Anp(:)
    REAL*8              :: cs_n

    Anp = 0.d0
    cs_n = SQRT(ABS(2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2)))

    Anp(1) = 1./(3.*cs_n)*(U(2)**2/U(1)**3 - U(3)/U(1)**2)
    Anp(2) = - 1./(3.*cs_n)*U(2)*U(5)/U(1)**2
    Anp(3) = 1./(3.*cs_n)*U(5)/U(1)
    Anp(5) = cs_n
  END SUBROUTINE jacobianMatricesNP
#endif

  !*****************************************
  ! Jacobian matrix for face computations
  !****************************************
  SUBROUTINE jacobianMatricesFace(U, bn, An)
    REAL*8, INTENT(in)  :: U(:), bn
    REAL*8, INTENT(out) :: An(:, :)
    An = 0.d0
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
       An(6, 1) = -U(6)*U(2)/U(1)**2
       An(6, 2) = U(6)/U(1)
       An(6, 6) = U(2)/U(1)
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
       An(6, 1) = -U(6)*U(2)/U(1)**2
       An(6, 2) = U(6)/U(1)
       An(6, 6) = U(2)/U(1)
#endif
#ifdef NEUTRAL
#ifdef NEUTRALCONVECTION
       An(5, 1) = -U(5)*U(2)/U(1)**2
       An(5, 2) = U(5)/U(1)
       An(5, 5) = U(2)/U(1)
       !An(5, :) = simpar%refval_time/(simpar%refval_length**2*phys%diff_n)*An(5,:)
#endif
#endif
    ENDIF
    An = bn*An
  END SUBROUTINE jacobianMatricesFace

#ifdef NEUTRALP
  SUBROUTINE jacobianMatricesFaceNP(U, bn, Anpn)
    REAL*8, INTENT(in)  :: U(:), bn
    REAL*8, INTENT(out) :: Anpn(:)
    REAL*8              :: cs_n

    Anpn = 0.d0
    cs_n = SQRT(2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2))

    Anpn(1) = 1./(3.*cs_n)*(U(2)**2/U(1)**3 - U(3)/U(1)**2)
    Anpn(2) = - 1./(3.*cs_n)*U(2)*U(5)/U(1)**2
    Anpn(3) = 1./(3.*cs_n)*U(5)/U(1)
    Anpn(5) = cs_n

    Anpn = bn*Anpn
  END SUBROUTINE jacobianMatricesFaceNP
#endif

  !*****************************************
  ! Jacobian matrix for the Bohm BC
  !****************************************
  SUBROUTINE jacobianMatricesBohm(U, A)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: A(:, :)
    REAL*8              :: auxi, auxe

    A = 0.
    auxi = (5.-2.*phys%Gmbohm)/3.
    auxe = (5.-2.*phys%Gmbohme)/3.
    A(1, 1) = -auxi*(U(2)**3/U(1)**3 - U(2)*U(3)/U(1)**2)
    A(1, 2) = auxi*(U(3)/U(1) - 3./2.*U(2)**2/U(1)**2)
    A(1, 3) = auxi*U(2)/U(1)

    A(2, 1) = -auxe*U(2)*U(4)/U(1)**2
    A(2, 2) = auxe*U(4)/U(1)
    A(2, 4) = auxe*U(2)/U(1)

#ifdef NEUTRAL
#ifdef NEUTRALCONVECTION
    A(5, 1) = -U(5)*U(2)/U(1)**2
    A(5, 2) = U(5)/U(1)
    A(5, 5) = U(2)/U(1)
    !A(5, :) = simpar%refval_time/(simpar%refval_length**2*phys%diff_n)*A(5,:)
#endif
#endif
  END SUBROUTINE jacobianMatricesBohm

#ifdef NEUTRALP
  SUBROUTINE jacobianMatricesBohmNP(U, Anp)
    REAL*8, INTENT(in)  :: U(:)
    REAL*8, INTENT(out) :: Anp(:)
    REAL*8              :: cs_n

    Anp = 0.d0
    cs_n = SQRT(ABS(2./3.*(U(3)/U(1) - 1./2.*U(2)**2/U(1)**2)))

    Anp(1) = 1./(3.*cs_n)*(U(2)**2/U(1)**3 - U(3)/U(1)**2)
    Anp(2) = - 1./(3.*cs_n)*U(2)*U(5)/U(1)**2
    Anp(3) = 1./(3.*cs_n)*U(5)/U(1)
    Anp(5) = cs_n
  END SUBROUTINE jacobianMatricesBohmNP
#endif

#ifdef NEUTRAL
  !*****************************************
  ! Jacobian matrices for Neutrals
  !****************************************
  SUBROUTINE jacobianMatricesN(U, Up, Q, b, sigmaviz, sigmavcx, Ax, Ay)
    REAL*8, INTENT(in)  :: U(:), Up(:), Q(:,:), b(:), sigmaviz, sigmavcx
    REAL*8, INTENT(out) :: Ax(:,:),Ay(:,:)
    REAL*8              :: GradTi(simpar%Ndim), GradTimod, Vnn, Csnn

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
    Ax(5,5) = Ax(5,5) - Up(2)*b(1)

    Ay = 0.d0
    !  if (Csnn .ge. Vnn) then
    !     Ay(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(2)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))
    !  else
    !     Ay(5,5) = simpar%refval_charge/(simpar%refval_mass*simpar%refval_density)*GradTi(2)/(U(1)*(abs(sigmaviz) + abs(sigmavcx)))*Csnn/Vnn
    !  endif
    !  Ay = Ay/simpar%refval_speed
    !Neutral velocity parallel to magnetic field
    Ay(5,5) = Ay(5,5) - Up(2)*b(2)

    Ax = 0.
    Ay = 0.

  END SUBROUTINE jacobianMatricesN
#endif

  !*****************************************
  ! Set the perpendicular diffusion
  !****************************************
  SUBROUTINE setLocalDiff(xy, u, d_iso, d_ani)
    REAL*8, INTENT(in)  		:: xy(:, :)
    REAL*8, INTENT(in)  		:: u(:,:)
    REAL*8, INTENT(out)		  :: d_iso(:, :, :), d_ani(:, :, :)
    REAL*8		              :: iperdiff(SIZE(xy, 1))
#ifdef NEUTRAL
    INTEGER             		:: i
    REAL*8				          :: ti_min=1e-6,ti
    REAL*8, DIMENSION(SIZE(u,1))	:: U1, U2, U3, U4, U5, sigmaviz, sigmavcx, Dnn
#ifdef KEQUATION
    REAL*8, DIMENSION(SIZE(u,1))          :: D_k,U6,c_s
    REAL*8                         :: r
#endif
#endif


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
    IF ((switch%ME .EQV. .TRUE.)  .AND. (switch%testcase .GT. 84)) THEN
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
       phys%ME_diff_n = d_iso(1,1,1)
       phys%ME_diff_u = d_iso(2,2,1)
       phys%ME_diff_e = d_iso(3,3,1)
       phys%ME_diff_ee = d_iso(4,4,1)
    ENDIF
#ifndef NEUTRALP
#ifdef NEUTRAL
    !d_iso(5, 5, :) = phys%diff_nn
    U1 = u(:,1)
    U2 = u(:,2)
    U3 = u(:,3)
    U4 = u(:,4)
    U5 = u(:,5)
#ifdef KEQUATION
    U6 = u(:,6)
#endif
#ifndef CONSTANTNEUTRALDIFF
    DO i=1,SIZE(u,1)
       CALL compute_sigmaviz(u(i,:),sigmaviz(i))
       CALL compute_sigmavcx(u(i,:),sigmavcx(i))
    END DO
    !ti = max(simpar%refval_temperature*2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2),0.1)
    !Dnn = simpar%refval_charge*ti/(simpar%refval_mass*simpar%refval_density*U1*(sigmaviz + sigmavcx))
    !Dnn = Dnn*simpar%refval_time/simpar%refval_length**2

    !Set a threshold on Dnn
    DO i=1,SIZE(Dnn,1)
#ifndef DNNSMOOTH
       ti = MAX(simpar%refval_temperature*2./(3.*phys%Mref)*(U3(i)/U1(i) - 1./2.*(U2(i)/U1(i))**2),ti_min)
       Dnn(i) = simpar%refval_charge*ti/(simpar%refval_mass*simpar%refval_density*U1(i)*(sigmaviz(i) + sigmavcx(i)))
       Dnn(i) = Dnn(i)*simpar%refval_time/simpar%refval_length**2
       IF (Dnn(i) .GT.  phys%diff_nn) THEN
          d_iso(5,5,i) = phys%diff_nn
       ELSEIF (Dnn(i)<10.*d_iso(1,1,i))THEN
          d_iso(5,5,i) = 10.*d_iso(1,1,i)
       ELSE
          d_iso(5,5,i) = Dnn(i)
       ENDIF
#else
       ti = simpar%refval_temperature*2./(3.*phys%Mref)*(U3(i)/U1(i) - 1./2.*(U2(i)/U1(i))**2)
       CALL softplus(ti, ti_min)
       Dnn(i) = simpar%refval_charge*ti/(simpar%refval_mass*simpar%refval_density*U1(i)*(sigmaviz(i) + sigmavcx(i)))
       Dnn(i) = Dnn(i)*simpar%refval_time/simpar%refval_length**2

       CALL double_softplus(Dnn(i),10.*phys%diff_n,phys%diff_nn)
       d_iso(5,5,i) = Dnn(i)


#endif

    END DO
#else
    d_iso(5,5,:)=phys%diff_nn
#endif

#ifdef KEQUATION
    DO i= 1,SIZE(c_s, 1)
       CALL compute_cs(u(i,:), c_s(i))
       ! for all equations
       IF (c_s(i)<=1.e-20) THEN
          D_k(i) = phys%diff_k_min
          !WRITE(6,*) 'NEGATIVE C_S ', c_s(i)
          !stop
       ELSE
          IF ((switch%testcase .GE. 50) .AND.(switch%testcase .LE. 59)) THEN
             r = xy(i,1)
          ELSEIF ((switch%testcase .GE. 60) .AND.(switch%testcase .LE. 69)) THEN
             r = xy(i,1) + geom%R0/simpar%refval_length
          ENDIF
          D_k(i) = r*U6(i)/c_s(i)
          IF (switch%testcase == 60) THEN
             D_k(i) = D_k(i)*geom%q
          ENDIF

#ifndef KDIFFSMOOTH

          D_k(i) = MAX(phys%diff_k_min,MIN(phys%diff_k_max,D_k(i) ))
#else
          !for circular case q_cyl assume constant

          CALL double_softplus(D_k(i),phys%diff_k_min,phys%diff_k_max)
#endif

       ENDIF
    ENDDO
    d_iso(6,6,:) = D_k+phys%diff_n
    d_ani(6,6,:) = d_iso(6,6,:)
    d_iso(1,1,:) = d_iso(1,1,:) + D_k
    d_iso(2,2,:) = d_iso(2,2,:) + D_k
    d_iso(3,3,:) = d_iso(3,3,:) + D_k
    d_iso(4,4,:) = d_iso(4,4,:) + D_k
    !WRITE(6,*) d_iso(6,6,:)*simpar%refval_length**2/simpar%refval_time
    !WRITE(6,*) d_iso(1,1,:)*simpar%refval_length**2/simpar%refval_time
#endif
    !Dnn = sum(Dnn)/size(u,1)
#endif
#else
    d_iso(5,5,:) = 0.
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
    d_ani(5,5,:) = 0.
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
  END SUBROUTINE setLocalDiff

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


  END SUBROUTINE computeIperDiffusion

  !*****************************************
  ! Pinch term
  !***************************************
  SUBROUTINE computePinch(b,psi,APinch)
    REAL*8, INTENT(IN)     :: b(:),psi
    REAL*8, INTENT(OUT)    :: APinch(:,:)
    REAL*8                 :: v_p,bnorm(2)

    APinch = 0.

    bnorm = b(:)/NORM2(b)
    v_p = phys%v_p*(psi**2 + psi**2*TANH((0.95 - psi)/0.02))
    !IF (v_p .lt. 1.e-4/simpar%refval_speed) v_p = 0.

    APinch(1,1) = v_p*bnorm(2)
    APinch(1,2) = v_p*(-bnorm(1))

    !WRITE(6,*) 'psi = ', psi
    !WRITE(6,*) 'v_p = ', v_p
  END SUBROUTINE computePinch

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
  END SUBROUTINE GimpMatrix

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
  END SUBROUTINE computeVi

  SUBROUTINE computeVe(U, V)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: V(:)
    V = 0.d0
    V(1) = -U(4)/U(1)**2
    V(4) = 1./U(1)
  END SUBROUTINE computeVe

  !                                SUBROUTINE computeVe(U,V)
  !                                real*8, intent(IN)  :: U(:)
  !                                real*8, intent(OUT) :: V(:)
  !                                V = 0.d0
  !                                V(1) = U(2)**2/U(1)**3 - U(4)/U(1)**2
  !                                V(2) = -U(2)/U(1)**2
  !                                V(4) = 1./U(1)
  !                                END SUBROUTINE computeVe

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
  END SUBROUTINE compute_dV_dUi

  SUBROUTINE compute_dV_dUe(U, dV_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dV_dU(:, :)
    dV_dU = 0.
    dV_dU(1, 1) = 2*U(4)/U(1)**3
    dV_dU(1, 4) = -1/U(1)**2
    dV_dU(4, 1) = -1/U(1)**2
  END SUBROUTINE compute_dV_dUe

  FUNCTION computeAlphai(U) RESULT(res)
    REAL*8 :: U(:)
    REAL*8 :: res, aux
    REAL, PARAMETER :: tol = 1e-5
    aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2
    IF ((2./(3.*phys%Mref)*aux > 1.) .AND. (switch%testcase .NE. 2)) THEN
       res = (3.*phys%Mref/2)**(phys%epn)
    ELSE
       IF (aux<tol) aux = tol
       res = aux**phys%epn
    ENDIF
    !! applying softplus instead strong limit
    !if (switch%testcase .ne. 2) then
    !  call double_softplus(aux, tol, 3.*phys%Mref/2)
    !endif
    !res = aux**phys%epn
  END FUNCTION computeAlphai

  FUNCTION computeAlphae(U) RESULT(res)
    REAL*8 :: U(:)
    REAL*8 :: res, aux
    REAL, PARAMETER :: tol = 1e-5
    aux = U(4)/U(1)
    IF ((2./(3.*phys%Mref)*aux > 1.) .AND. (switch%testcase .NE. 2)) THEN
       res = (3.*phys%Mref/2)**(phys%epn)
    ELSE
       IF (aux<tol) aux = tol
       res = aux**phys%epn
    ENDIF
    !! applying softplus instead strong limit
    !if (switch%testcase .ne. 2) then
    !  call double_softplus(aux, tol, 3.*phys%Mref/2)
    !endif
    !res = aux**phys%epn
  END FUNCTION computeAlphae

  SUBROUTINE compute_dAlpha_dUi(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(OUT):: res(:)
    REAL*8             :: aux
    REAL, PARAMETER :: tol = 1e-5
    ! applying softplus instead strong limit
    aux = U(3)/U(1) - 0.5*U(2)**2/U(1)**2
    !double_soft_deriv = 1.
    !if (switch%testcase .ne. 2) then  !! don't apply flux limiter if it is a convergence test
    !  call double_softplus_deriv(aux, tol,3.*phys%Mref/2,double_soft_deriv)
    !  call double_softplus(aux, tol, 3.*phys%Mref/2)
    !endif
    !res = 0.
    !res(1) = -U(3)/U(1)**2+U(2)**2/U(1)**3
    !res(2) = -U(2)/U(1)**2
    !res(3) = 1./U(1)
    !res=phys%epn*aux**(phys%epn-1)*res
    !if (switch%testcase .ne. 2) then
    !  res = res*double_soft_deriv
    !endif
    IF ((2./(3.*phys%Mref)*aux > 1.) .AND. (switch%testcase .NE. 2)) THEN !! don't apply flux limiter if it is a convergence test
       res = 0.
    ELSE
       IF (aux<0) aux = tol
       res = 0.d0
       res(1) = -U(3)/U(1)**2+U(2)**2/U(1)**3
       res(2) = -U(2)/U(1)**2
       res(3) = 1./U(1)
       res=phys%epn*aux**(phys%epn-1)*res
    ENDIF
  END SUBROUTINE compute_dAlpha_dUi

  SUBROUTINE compute_dAlpha_dUe(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(OUT):: res(:)
    REAL*8             :: aux
    REAL, PARAMETER :: tol = 1e-5
    ! applying softplus instead strong limit
    aux = U(4)/U(1)
    !double_soft_deriv = 1.
    !if (switch%testcase .ne. 2) then  !! don't apply flux limiter if it is a convergence test
    !  call double_softplus_deriv(aux, tol,3.*phys%Mref/2,double_soft_deriv)
    !  call double_softplus(aux, tol, 3.*phys%Mref/2)
    !endif
    !res = 0.
    !res(1) = -U(4)/U(1)**2
    !res(4) = 1./U(1)
    !res=phys%epn*aux**(phys%epn-1)*res
    !if (switch%testcase .ne. 2) then
    !  res = res*double_soft_deriv
    !endif
    IF ((2./(3.*phys%Mref)*aux > 1.) .AND. (switch%testcase .NE. 2)) THEN !! don't apply flux limiter if it is a convergence test
       res = 0.
    ELSE
       IF (aux<0) aux = tol
       res = 0.d0
       res(1) = -U(4)/U(1)**2
       res(4) = 1./U(1)
       res=phys%epn*aux**(phys%epn-1)*res
    ENDIF
  END SUBROUTINE compute_dAlpha_dUe

  ! ******************************
  ! Parallel electric field terms
  ! ******************************
  SUBROUTINE compute_W(U, W)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: W

    W = 2./3.*U(2)/U(1)
  END SUBROUTINE compute_W

  SUBROUTINE compute_dW_dU(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:, :)
    res = 0.
    res(4, 1) = -2./3.*U(2)/U(1)**2
    res(4, 2) = 2./3./U(1)
  END SUBROUTINE compute_dW_dU

  ! ******************************
  ! Temperature exchange terms
  ! ******************************
  SUBROUTINE compute_s(U, s)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: s, U1, U4, U3
    REAL, PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    U3 = U(3)
    IF (U4 < tol) U4 = tol
    IF (U1 < tol) U1 = tol
    IF (U3 < tol) U3 = tol
    ! keeping this thing for high diffusion, but take care for low values
    IF ((phys%diff_ee .GT. 0.0380) .AND. (switch%testcase .NE. 2)) THEN
       s = 1./(phys%tie*0.0380/phys%diff_ee)*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
    ELSE
       s = 1./(phys%tie)*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4-U3+0.5*(U(2)**2/U1))
    ENDIF
    !s = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*(U1**(2.5)/U4**1.5)*(U4 - U3 + 0.5*(U(2)**2/U1))
  END SUBROUTINE compute_s

  SUBROUTINE compute_ds_dU(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), U1, U4, U3
    REAL, PARAMETER :: tol = 1e-20
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
    IF ((phys%diff_ee .GT. 0.0380) .AND. (switch%testcase .NE. 2)) THEN
       res = 1./(phys%tie*0.0380/phys%diff_ee)*(2./3./phys%Mref)**(-0.5)*res
    ELSE
       res = 1./(phys%tie)*(2./3./phys%Mref)**(-0.5)*res
    ENDIF
    !res = 1./phys%tie*(2./3./phys%Mref)**(-0.5)*res
  END SUBROUTINE compute_ds_dU


  !*****************************
  !Ohmic Heating Source
  !*****************************
  SUBROUTINE compute_Sohmic(U,Sohmic)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: Sohmic,U1,U4
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    IF (U4<tol) U4=tol
    IF (U1<tol) U1=tol
    Sohmic = phys%Pohmic*(((3*phys%Mref)/2)**1.5)*((U1/U4)**1.5)
  END SUBROUTINE compute_Sohmic


  SUBROUTINE compute_dSohmic_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    IF (U4<tol) U4=tol
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = (1.5*U1**0.5)/(U4**1.5)
    res(4) = -(1.5*U1**1.5)/(U4**2.5)
    res = phys%Pohmic*(((3*phys%Mref)/2)**1.5)*res
  END SUBROUTINE compute_dSohmic_dU


  ! ******************************
  ! Neutral Source terms
  ! ******************************
#ifdef NEUTRAL

  SUBROUTINE compute_niz(U,niz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: niz,U1,U5
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    niz = U1*U5
  END SUBROUTINE compute_niz


  SUBROUTINE compute_dniz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U5
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.
    res(1) = U5
    res(5) = U1
  END SUBROUTINE compute_dniz_dU


  SUBROUTINE compute_nrec(U,nrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: nrec,U1
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    IF (U1<tol) U1=tol
    nrec = U1**2
  END SUBROUTINE compute_nrec


  SUBROUTINE compute_dnrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = 2.*U1
  END SUBROUTINE compute_dnrec_dU


  SUBROUTINE compute_fGammacx(U,fGammacx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fGammacx,U2,U5
    REAL,PARAMETER :: tol = 1e-20
    U2 = U(2)
    U5 = U(5)
    IF (U5<tol) U5=tol
    fGammacx = U2*U5
  END SUBROUTINE compute_fGammacx


  SUBROUTINE compute_dfGammacx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U2,U5
    REAL,PARAMETER :: tol = 1e-20
    U2 = U(2)
    U5 = U(5)
    IF (U5<tol) U5=tol
    res = 0.
    res(2) = U5
    res(5) = U2
  END SUBROUTINE compute_dfGammacx_dU


  SUBROUTINE compute_fGammarec(U,fGammarec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fGammarec,U1,U2
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    IF (U1<tol) U1=tol
    fGammarec = U1*U2
  END SUBROUTINE compute_fGammarec


  SUBROUTINE compute_dfGammarec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U2
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    IF (U1<tol) U1=tol
    res = 0.
    res(1) = U2
    res(2) = U1
  END SUBROUTINE compute_dfGammarec_dU


#ifdef TEMPERATURE
#ifndef AMJUELSPLINES
  SUBROUTINE compute_sigmaviz(U,sigmaviz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmaviz,U1,U4,T0,Ery,E0
    REAL, PARAMETER    :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    Ery = 13.6
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    E0 = (2.*T0*U4)/(3.*phys%Mref*Ery*U1)
    !Threshold on Te >= 0.2 eV
    IF (E0*Ery .LE. 0.05) E0 = 0.05/Ery
    sigmaviz = 1.e-11*SQRT(E0)*(1./((Ery**1.5)*(6. + E0)))*EXP(-1./E0)
  END SUBROUTINE compute_sigmaviz


  SUBROUTINE compute_dsigmaviz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0,Ery,E0
    REAL, PARAMETER    :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    Ery = 13.6
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    res = 0.
    E0 = (2.*T0*U4)/(3.*phys%Mref*Ery*U1)
    !Threshold on Te >= 0.2 eV
    !if (E0*Ery .le. 0.05)  E0 = 0.05/Ery
    IF (E0*Ery .GE. 0.05) THEN
       res(1) = (1./U1)*( - 1./2 + 1./((6./E0) + 1)   - 1./E0 )
       res(4) = 1./(2.*U4) - 1./(6.*(U4/E0) + U4) + 1./(E0*U4)
       res = 1.e-11*SQRT(E0)*(1./((Ery**1.5)*(6. + E0)))*EXP(-1./E0)*res
    END IF
  END SUBROUTINE compute_dsigmaviz_dU


  SUBROUTINE compute_sigmavrec(U,sigmavrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavrec,U1,U4,T0,Ery,E0
    REAL, PARAMETER    :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    Ery = 13.6
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    E0 = (Ery*3.*phys%Mref*U1)/(T0*2.*U4)
    !Threshold on Te >= 0.1 eV
    IF (E0/Ery .GE. 10.) E0 = 10.*Ery
    sigmavrec = 5.2e-20*SQRT(E0)*(0.43 + 0.5*LOG(E0) + 0.469*(E0**(-1./3)))
  END SUBROUTINE compute_sigmavrec


  SUBROUTINE compute_dsigmavrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0,Ery,E0
    REAL, PARAMETER    :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    Ery = 13.6
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    E0 = (Ery*3.*phys%Mref*U1)/(T0*2.*U4)
    res = 0.
    !Threshold on Te >= 0.1 eV
    IF (E0/Ery .LE. 10.) THEN
       res(1) = 5.2e-20*SQRT(E0)*(0.43 + 0.5*LOG(E0) + 0.469*(E0**(-1./3)))/(2.*U1) + 5.2e-20*SQRT(E0)*(0.5/U1 +&
            &0.469*(E0**(-1./3))*(-1./(3.*U1)))
       res(4) = - 5.2e-20*SQRT(E0)*(0.43 + 0.5*LOG(E0) + 0.469*(E0**(-1./3)))/(2.*U4) + 5.2e-20*SQRT(E0)*( - 0.5/U4 +&
            &0.469*(E0**(-1./3))*(1./(3.*U4)))
    ENDIF
  END SUBROUTINE compute_dsigmavrec_dU
#else


  ! Routines for 2D splines in te, ne loglogspace
  SUBROUTINE compute_2D_eirene_rate(te,ne,alpha,rate)
    ! this routine calculates eirene rate in te, ne space
    ! if in the region of applicability, then just takes the values according to the splines
    ! if ne>ne_max (LTE limit) or ne<ne_min (Corona limit) then no more dependancy on ne
    ! if te<te_min or te>te_max then extrapolates in log space taking the derivative of the edge and linearly expanding
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    REAL*8             :: te_min = 0.1
    REAL*8             :: te_max = 2.e4
    REAL*8             :: ne_min = 1.
    REAL*8             :: ne_max = 1.e8
    REAL*8             :: dlograte_dlogte

    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    ! region 1, where ne is applicable
    IF ((ne>=ne_min) .AND. (ne<=ne_max)) THEN
       IF ((te>=te_min) .AND.(te<=te_max)) THEN
          CALL compute_2D_logeirene_rate(te,ne,alpha,rate)
       ELSEIF (te<te_min) THEN
          CALL compute_2D_logeirene_rate(te_min,ne,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(te_min,ne,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(te_min))
       ELSEIF (te>te_max) THEN
          CALL compute_2D_logeirene_rate(te_max,ne,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(te_max,ne,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(te_max))
       ENDIF
       ! beyond range of ne applicability
    ELSEIF(ne<ne_min) THEN
       ! if te is still applicable
       IF ((te>=te_min) .AND.(te<=te_max)) THEN
          CALL compute_2D_logeirene_rate(te,ne_min,alpha,rate)
          ! if te < te_min
       ELSEIF (te<te_min) THEN
          CALL compute_2D_logeirene_rate(te_min,ne_min,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(te_min,ne_min,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(te_min))
          ! if te>te_max
       ELSEIF (te>te_max) THEN
          CALL compute_2D_logeirene_rate(te_max,ne_min,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(te_max,ne_min,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(te_max))
       ENDIF
       ! if ne>ne_max
    ELSEIF (ne>ne_max) THEN
       IF ((te>=te_min) .AND.(te<=te_max)) THEN
          CALL compute_2D_logeirene_rate(te,ne_max,alpha,rate)
       ELSEIF (te<te_min) THEN
          CALL compute_2D_logeirene_rate(te_min,ne_max,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(te_min,ne_max,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(te_min))
       ELSEIF (te>te_max) THEN
          CALL compute_2D_logeirene_rate(te_max,ne_max,alpha,rate)
          CALL compute_dlogeirene_2D_dlogte_rate(te_max,ne_max,alpha,dlograte_dlogte)
          rate = rate+dlograte_dlogte*(LOG(te)-LOG(te_max))
       ENDIF
    ENDIF
    ! rate is in cm^3/s in EIRENE
    rate = EXP(rate)/1.e6
  END SUBROUTINE compute_2D_eirene_rate

  SUBROUTINE compute_2D_eirene_rate_du(U1,U4,te,ne,alpha,rate_du)
    REAL*8, INTENT(IN) :: U1,U4,te,ne,alpha(:,:)
    REAL*8             :: rate
    REAL*8             :: te_min = 0.1
    REAL*8             :: te_max = 2.e4
    REAL*8             :: ne_min = 1.
    REAL*8             :: ne_max = 1.e8
    REAL*8             :: dlograte_dlogne, dlograte_dlogte
    REAL*8, INTENT(OUT):: rate_du(:)

    CALL compute_2D_eirene_rate(te,ne,alpha,rate)
    ! In EIRENE the density is scaled to 1.e14
    rate_du = 0.
    dlograte_dlogne = 0.
    dlograte_dlogte = 0.
    IF ((ne>=ne_min) .AND. (ne<=ne_max)) THEN
       IF ((te>=te_min) .AND.(te<=te_max)) THEN
          CALL compute_dlogeirene_2D_dlogne_rate(te,ne,alpha,dlograte_dlogne)
          CALL compute_dlogeirene_2D_dlogte_rate(te,ne,alpha,dlograte_dlogte)
       ELSEIF (te<te_min) THEN
          CALL compute_dlogeirene_2D_dlogne_rate(te_min,ne,alpha,dlograte_dlogne)
          CALL compute_dlogeirene_2D_dlogte_rate(te_min,ne,alpha,dlograte_dlogte)
       ELSEIF (te>te_max) THEN
          CALL compute_dlogeirene_2D_dlogne_rate(te_max,ne,alpha,dlograte_dlogne)
          CALL compute_dlogeirene_2D_dlogte_rate(te_max,ne,alpha,dlograte_dlogte)
       ENDIF
       ! beyond range of ne applicability (ne derivative is now zero)
    ELSEIF(ne<ne_min) THEN
       ! if te is still applicable
       IF ((te>=te_min) .AND.(te<=te_max)) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te,ne_min,alpha,dlograte_dlogte)
          ! if te < te_min
       ELSEIF (te<te_min) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te_min,ne_min,alpha,dlograte_dlogte)
          ! if te>te_max
       ELSEIF (te>te_max) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te_max,ne_min,alpha,dlograte_dlogte)
       ENDIF
       ! if ne>ne_max (ne derivative is now zero)
    ELSEIF (ne>ne_max) THEN
       IF ((te>=te_min) .AND.(te<=te_max)) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te,ne_max,alpha,dlograte_dlogte)
       ELSEIF (te<te_min) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te_min,ne_max,alpha,dlograte_dlogte)
       ELSEIF (te>te_max) THEN
          CALL compute_dlogeirene_2D_dlogte_rate(te_max,ne_max,alpha,dlograte_dlogte)
       ENDIF
    ENDIF
    rate_du(1) = rate_du(1) + dlograte_dlogte*(-1./U1)
    rate_du(1) = rate_du(1) + dlograte_dlogne*(1./U1)
    rate_du(4) = rate_du(4) + dlograte_dlogte*(1./U4)
    rate_du = rate_du*rate
  END SUBROUTINE compute_2D_eirene_rate_du

  SUBROUTINE compute_2D_logeirene_rate(te,ne,alpha,rate)
    ! this routine calculate log (eirene_rate) for given te, ne in log log space
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    INTEGER            :: i,j
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    DO j=1,SIZE(alpha,2)
       DO i = 1,SIZE(alpha,1)
          rate = rate + alpha(i,j)*LOG(ne)**(j-1)*LOG(te)**(i-1)
       END DO
    END DO
  END SUBROUTINE compute_2D_logeirene_rate

  SUBROUTINE compute_dlogeirene_2D_dlogte_rate(te,ne,alpha,rate)
    ! this routines calculate derivative dlog (eirene_rate)/dlog(te) for given te, ne in log log space
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    INTEGER            :: i,j
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    DO j=1,SIZE(alpha,2)
       DO i = 2,SIZE(alpha,1)
          rate = rate + alpha(i,j)*(i-1)*LOG(ne)**(j-1)*LOG(te)**(i-2)
       END DO
    END DO
  END SUBROUTINE compute_dlogeirene_2D_dlogte_rate

  SUBROUTINE compute_dlogeirene_2D_dlogne_rate(te,ne,alpha,rate)
    ! this routine calculate derivative dlog (eirene_rate)/dlog(ne) for given te, ne in log log space
    REAL*8, INTENT(IN) :: te,ne,alpha(:,:)
    REAL*8, INTENT(OUT):: rate
    INTEGER            :: i,j
    ! In EIRENE the density is scaled to 1.e14
    rate = 0.
    DO j=2,SIZE(alpha,2)
       DO i = 1,SIZE(alpha,1)
          rate = rate + alpha(i,j)*(j-1)*LOG(ne)**(j-2)*LOG(te)**(i-1)
       END DO
    END DO
  END SUBROUTINE compute_dlogeirene_2D_dlogne_rate




  SUBROUTINE compute_sigmaviz(U,sigmaviz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmaviz,U1,U4,T0,te,ne,n0
    REAL, PARAMETER    :: tol = 1.e-20  !tolerance for U4 = 3/2*Mref*U1min*te_min/T0

    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19

    IF ((U1>tol) .AND. (U4>tol)) THEN ! basically it's a below zero check
       te = T0*2/3./phys%Mref*U4/U1
       ne = n0*U1/1.e14
    ELSE!some low values
       ne = n0*1.e-20/1.e14
       te = 1.e-10
    ENDIF

    sigmaviz = 0.

    CALL compute_2D_eirene_rate(te,ne,phys%alpha_iz,sigmaviz)
  END SUBROUTINE compute_sigmaviz

  SUBROUTINE compute_dsigmaviz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0, te,ne,n0
    REAL, PARAMETER    :: tol = 1.e-20 !tolerance for U4 = 3/2*Mref*U1min*te_min/T0

    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19

    res = 0.
    IF ((U1>tol) .AND. (U4>tol)) THEN ! basically it's a below zero check
       ne = n0*U1/1.e14
       te = T0*2/3./phys%Mref*U4/U1
       CALL compute_2D_eirene_rate_du(U1,U4,te,ne,phys%alpha_iz,res)
    ENDIF !let non-linear part as zero if negative solutions

  END SUBROUTINE compute_dsigmaviz_dU

  SUBROUTINE compute_sigmavrec(U,sigmavrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavrec,U1,U4,T0,te,ne,n0
    REAL, PARAMETER    :: tol = 1.e-20  !tolerance for U4 = 3/2*Mref*U1min*te_min/T0

    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19

    IF ((U1>tol) .AND. (U4>tol)) THEN ! basically it's a below zero check
       te = T0*2/3./phys%Mref*U4/U1
       ne = n0*U1/1.e14
    ELSE!some low values
       ne = n0*1.e-20/1.e14
       te = 1.e-10
    ENDIF
    CALL compute_2D_eirene_rate(te,ne,phys%alpha_rec,sigmavrec)
  END SUBROUTINE compute_sigmavrec

  SUBROUTINE compute_dsigmavrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:), U1, U4, T0, te, ne, n0
    REAL, PARAMETER    :: tol = 1.e-20 !tolerance for U4 = 3/2*Mref*U1min*te_min/T0

    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    n0 = 1.e19

    res = 0.
    IF ((U1>tol) .AND. (U4>tol)) THEN ! basically it's a below zero check
       ne = n0*U1/1.e14
       te = T0*2/3./phys%Mref*U4/U1
       CALL compute_2D_eirene_rate_du(U1,U4,te,ne,phys%alpha_rec,res)
    ENDIF !let non-linear part as zero if negative solutions
  END SUBROUTINE compute_dsigmavrec_dU
#endif
#ifdef MANUELCX
  !ADAS truncated CX
  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavcx,U1,U4,T0,E0
    REAL*8              :: p1,p2,p3,p4,p5
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    !Analytical expression from Hugo Bufferand
    !E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    !sigmavcx = (2.5e-15/exp(-0.5))*exp(-E0)
    !Fit log log from ADAS database
    p1 = -0.0006837
    p2 = 0.008004
    p3 = -0.03581
    p4 = 0.4518
    p5 = -32.59
    E0 = T0*2./(3.*phys%Mref)*U4/U1

    !Threshold on Te >= 0.2 eV
    IF (E0 .LE. 0.05) E0 = 0.05
    sigmavcx = EXP(p1*LOG(E0)**4 + p2*LOG(E0)**3 + p3*LOG(E0)**2 + p4*LOG(E0) + p5)
  END SUBROUTINE compute_sigmavcx


  SUBROUTINE compute_dsigmavcx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0,E0
    REAL*8             :: p1,p2,p3,p4,p5
    REAL, PARAMETER    :: tol = 1e-20
    T0 = 50.
    U1 = U(1)
    U4 = U(4)
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    !Analytical expression from Hugo Bufferand
    !E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    !res = 0.
    !res(1) = -1./U4
    !res(4) = U1/(U4**2)
    !res = (1.5*phys%Mref/(2.*T0))*(2.5e-15/exp(-0.5))*exp(-E0)*res
    !Fit log log from ADAS database
    p1 = -0.0006837
    p2 = 0.008004
    p3 = -0.03581
    p4 = 0.4518
    p5 = -32.59
    E0 = T0*2./(3.*phys%Mref)*U4/U1
    res = 0.
    !Threshold on Te >= 0.2 eV
    !if (E0 .le. 0.05) E0 = 0.05
    IF (E0 .GE. 0.05) THEN
       res(1) = (4.*p1*LOG(E0)**3 + 3.*p2*LOG(E0)**2 + 2.*p3*LOG(E0) + p4)*(-U4/U1**2)
       res(4) = (4.*p1*LOG(E0)**3 + 3.*p2*LOG(E0)**2 + 2.*p3*LOG(E0) + p4)*1./U1
       res = EXP(p1*LOG(E0)**4 + p2*LOG(E0)**3 + p3*LOG(E0)**2 + p4*LOG(E0) + p5)*U1/U4*res
    END IF
  END SUBROUTINE compute_dsigmavcx_dU
#endif
#ifdef LEGACYCX
  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: sigmavcx,U1,U4,T0,E0
    REAL*8              :: p1,p2,p3,p4,p5
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    !Analytical expression from Hugo Bufferand
    !E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    !sigmavcx = (2.5e-15/exp(-0.5))*exp(-E0)

    E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    sigmavcx = (2.5e-15/EXP(-0.5))*EXP(-E0)
  END SUBROUTINE compute_sigmavcx

  SUBROUTINE compute_dsigmavcx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0,E0
    REAL*8             :: p1,p2,p3,p4,p5
    REAL, PARAMETER    :: tol = 1e-20
    T0 = 50.
    U1 = U(1)
    U4 = U(4)
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    !Analytical expression from Hugo Bufferand

    E0 = (0.5*3.*phys%Mref*U1)/(2.*T0*U4)
    res = 0.
    res(1) = -1./U4
    res(4) = U1/(U4**2)
    res = (1.5*phys%Mref/(2.*T0))*(2.5e-15/EXP(-0.5))*EXP(-E0)*res
  END SUBROUTINE compute_dsigmavcx_dU
#endif
#ifdef EXPANDEDCX
! These routines use AMUJUEL splines
  SUBROUTINE compute_eirene_1D_rate(ti,alpha,rate)
    ! This routine calculates extrapolated AMJUEL 1D rate (here on ion temperature) for given temperature and coefficients
    real*8, intent(IN) :: ti,alpha(:)          
    real*8, intent(OUT):: rate
    real*8             :: ti_min=0.1
    real*8             :: dlograte_dlogti
    integer            :: i
    rate = 0.
    if (ti>=ti_min) then 
      call compute_logeirene_1D_rate(ti,alpha,rate)
    else
      call compute_logeirene_1D_rate(ti_min,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogti(ti_min,alpha,dlograte_dlogti)
      rate = rate + dlograte_dlogti*(log(ti)- log(ti_min))
    endif
    ! rates are not higher than 1 m^3/s, if rate is higher than that value, then there is something weird
    if (rate>6.*log(10.)) then
      WRITE(6,*) "Something weird in compute_eirene_te_rate, probably, solution is not good already"
      WRITE(6,*) " ti equal to", ti
      WRITE(6,*) " rate equal to", rate
      stop
    endif
    rate = exp(rate)/1.e6
  END SUBROUTINE compute_eirene_1D_rate

  SUBROUTINE compute_eirene_1D_rate_du(U1,U2,U3,ti,dti_dU,alpha,res)
    ! This routine calculates extrapolated AMJUEL 1D rate (typically on temperature) for given temperature and coefficients
    real*8, intent(IN) :: U1,U2,U3,ti,dti_dU(:),alpha(:)          
    real*8, intent(OUT):: res(:)
    real*8             :: ti_min=0.1
    real*8             :: dlograte_dlogte,rate
    integer            :: i
    res = 0.

    if (ti>ti_min) then
      call compute_eirene_1D_rate(ti,alpha,rate)
      call compute_d_logeirene_1D_rate_dlogti(ti,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dti_dU(1)/ti)
      res(2) = res(2) + dlograte_dlogte*(dti_dU(2)/ti)
      res(3) = res(3) + dlograte_dlogte*(dti_dU(3)/ti)
      res = rate*res
    else
      call compute_eirene_1D_rate(ti,phys%alpha_cx,rate)
      call compute_d_logeirene_1D_rate_dlogti(ti_min,alpha,dlograte_dlogte)
      res(1) = res(1) + dlograte_dlogte*(dti_dU(1)/ti)
      res(2) = res(2) + dlograte_dlogte*(dti_dU(2)/ti)
      res(3) = res(3) + dlograte_dlogte*(dti_dU(3)/ti)
      res = rate*res
    endif

  END SUBROUTINE compute_eirene_1D_rate_du


  SUBROUTINE compute_logeirene_1D_rate(ti,alpha,rate)
    ! Calculates 1D AMJUEL spline in loglog space
    real*8, intent(IN) :: ti,alpha(:)          
    real*8, intent(OUT):: rate
    integer            :: i
    rate = 0.

    do i = 1,size(alpha,1)    
      rate = rate + alpha(i)*log(ti)**(i-1)
    end do

  END SUBROUTINE compute_logeirene_1D_rate


  SUBROUTINE compute_d_logeirene_1D_rate_dlogti(ti,alpha,d_log_rate_dti)
    ! calculates derivative of AMJUEL 1D spline in loglog space
    real*8, intent(IN) :: ti, alpha(:)          
    real*8, intent(OUT):: d_log_rate_dti
    integer            :: i
    d_log_rate_dti = 0.

    do i = 2,size(alpha,1)    
      d_log_rate_dti = d_log_rate_dti + (i-1)*alpha(i)*log(ti)**(i-2)
    end do
  END SUBROUTINE compute_d_logeirene_1D_rate_dlogti
  SUBROUTINE compute_sigmavcx(U,sigmavcx)
    ! calculates AMJUEL CX rate
    real*8, intent(IN) :: U(:)
    real*8             :: sigmavcx,U1,U2,U3,T0,E0,ti
    integer             :: i
    real,parameter :: tol = 1.e-20 !tolerance for U4 = 3/2*Mref*U1min*te_min/T0
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    T0 = 50.
    
    !if (U1<tol) U1=tol
    !if (U4<tol) then
    if ((U1>tol) .and. (U3>tol)) then ! basically it's a below zero check
      ti = T0*2/3. /phys%Mref * (U3/U1 - 1/2 *U2**2/U1**2)
    else!some low values  
      ti = 1.e-10
    endif
    sigmavcx = 0.
  
    call compute_eirene_1D_rate(ti, phys%alpha_cx, sigmavcx)
  END SUBROUTINE compute_sigmavcx
  
  
  SUBROUTINE compute_dsigmavcx_dU(U,res)
    ! calculates derivative of AMJUEL CX rate for linearization
    real*8, intent(IN) :: U(:)
    real*8             :: res(:), U1,U2,U3,T0,ti, ti_min = 0.1
    real*8, allocatable :: dti_dU(:)
    real*8             :: sigmavcx, sigmavcx_dte
    real, parameter    :: tol = 1.e-20  !tolerance for U4 = 3/2*Mref*U1min*te_min/T0
    integer            :: i
    
    allocate(dti_dU(size(U)))
    
    T0 = 50.
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    res = 0.
    dti_dU = 0.
    
    if ((U1>tol) .and. (U3>tol)) then
      ti = T0*2/3. /phys%Mref * (U3/U1 - 1/2 *U2**2/U1**2)
      dti_dU(1) = dti_dU(1) + 1.*(-U3 + U2**2/U1) / U1**2
      dti_dU(2) = dti_dU(2) - 1.*U2/U1**2
      dti_dU(3) = dti_dU(3) + 1./U1
      dti_dU(:) = dti_dU(:) * T0*2/3. /phys%Mref
      call compute_eirene_1D_rate_dU(U1,U2,U3,ti,dti_dU,phys%alpha_cx,res)      
    endif
    
  
  END SUBROUTINE compute_dsigmavcx_dU
#endif
#ifdef DNNSMOOTH
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
  END SUBROUTINE double_softplus

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
  END SUBROUTINE double_softplus_deriv

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
  END SUBROUTINE softplus

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
  END SUBROUTINE softplus_deriv


#endif
#ifdef DNNLINEARIZED
  SUBROUTINE compute_Dnn_dU(U, Dnn_dU)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(OUT) :: Dnn_dU(:)
    REAL*8              :: Dnn, ti,  ti_min=1e-6
#ifdef DNNSMOOTH
    REAL*8              :: double_soft_deriv, soft_deriv
#else
    REAL*8              :: ti_real
#endif
    REAL*8              :: sigmaviz, sigmavcx
    REAL*8              :: dti_du(SIZE(U,1)), dsigmaviz_dU(SIZE(U,1)), dsigmavcx_dU(SIZE(U,1))
    Dnn_dU(:) = 0.
    !if ((U(3)>=tol) .and. (U(1)>=tol) .and. (U(5)>=tol)) then
    ! calculation of atomic rates
    CALL compute_sigmaviz(U,sigmaviz)
    CALL compute_sigmavcx(U,sigmavcx)
    ! calculation of temperature before limitation
    ti = simpar%refval_temperature*2./(3.*phys%Mref)*(U(3)/U(1) - 1./2.*(U(2)/U(1))**2)
#ifdef DNNSMOOTH
    CALL softplus_deriv(ti, ti_min,soft_deriv)
    CALL softplus(ti,ti_min)
#else
    ti_real = ti
    ti = MAX(ti_min,ti)
#endif
    ! calculation of Dnn before limitation
    Dnn = simpar%refval_charge*ti/(simpar%refval_mass*simpar%refval_density*U(1)*(sigmaviz + sigmavcx))*simpar%refval_time/simpar%refval_length**2
#ifdef DNNSMOOTH
    CALL double_softplus_deriv(Dnn,10.*phys%diff_n,phys%diff_nn,double_soft_deriv)   !to check the mulptiplier for Dnn_min
#else
    IF ((Dnn>10.*phys%diff_n) .AND. (Dnn<phys%diff_nn)) THEN
#endif
       ! ti derivative
       dti_du(:) = 0.
#ifndef DNNSMOOTH
       IF (ti_real>ti) THEN
#endif
          dti_du(1) = -U(3)/U(1)**2+U(2)**2/U(1)**3
          dti_du(2) = -U(2)/U(1)**2
          dti_du(3) = 1./U(1)
          dti_du(:) = dti_du(:)*simpar%refval_temperature*2./(3.*phys%Mref)
#ifndef DNNSMOOTH
       ENDIF
#endif
       ! atomic rates derivatives
       CALL compute_dsigmaviz_dU(U,dsigmaviz_dU)
       CALL compute_dsigmavcx_dU(U,dsigmavcx_dU)

       ! arrange all ingredients
       Dnn_dU(:) = 0.
       ! ti part
       Dnn_dU(:) = Dnn_dU(:)+dti_du(:)*simpar%refval_charge/(simpar%refval_mass*simpar%refval_density*U(1)*(sigmaviz + sigmavcx))
#ifdef DNNSMOOTH
       Dnn_dU(:) = Dnn_dU(:)*soft_deriv
#endif
       ! n part
       Dnn_dU(1) = Dnn_dU(1)-ti*simpar%refval_charge/(simpar%refval_mass*simpar%refval_density*U(1)**2*(sigmaviz + sigmavcx))
       ! atomic rates part
       Dnn_dU(:) = Dnn_dU(:)-ti*simpar%refval_charge/(simpar%refval_mass*simpar%refval_density*U(1)*(sigmaviz + sigmavcx)**2)*(dsigmaviz_dU(:)+dsigmavcx_dU(:))

       Dnn_dU(:) = Dnn_dU(:)*simpar%refval_time/simpar%refval_length**2
#ifdef DNNSMOOTH
       Dnn_dU(:) = Dnn_dU(:)*double_soft_deriv
#endif
       !endif
#ifndef DNNSMOOTH
    ENDIF
#endif


  END SUBROUTINE  compute_Dnn_dU
#endif
  SUBROUTINE compute_Tloss(U,Tloss)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: Tloss,U1,U4,T0
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    Tloss = 25. + 170.*EXP(-(T0*U4)/(3.*phys%Mref*U1))
  END SUBROUTINE compute_Tloss


  SUBROUTINE compute_dTloss_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    res = 0.
    res(1) = U4/(U1**2)
    res(4) = -1./U1
    res = 170.*EXP(-(T0*U4)/(3.*phys%Mref*U1))*(T0/(3.*phys%Mref))*res
  END SUBROUTINE compute_dTloss_dU


  SUBROUTINE compute_Tlossrec(U,Tlossrec)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: Tlossrec,U1,U4,T0
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    !Tlossrec = min(250.,8.*exp((2.*U4*T0)/(3.*phys%Mref*U1*9.)))
    Tlossrec = (2.*U4*T0)/(3.*phys%Mref*U1*9.)
    IF (Tlossrec .LE. LOG(250./8.)) THEN
       Tlossrec = 8.*EXP(Tlossrec)
    ELSE
       Tlossrec = 250.
    ENDIF
  END SUBROUTINE compute_Tlossrec


  SUBROUTINE compute_dTlossrec_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U4,T0,Tlossrec
    REAL, PARAMETER    :: tol = 1e-20
    U1 = U(1)
    U4 = U(4)
    T0 = 50.
    IF (U1<tol) U1=tol
    IF (U4<tol) U4=tol
    res = 0.
    !Tlossrec = 8.*exp((2.*U4*T0)/(3.*phys%Mref*U1*9.))
    !if (Tlossrec .le. 250.) then
    !  res(1) = -U4/(U1**2)
    !  res(4) = 1./U1
    !  res = Tlossrec*((2.*T0)/(27.*phys%Mref))*res
    !endif
    Tlossrec = (2.*U4*T0)/(3.*phys%Mref*U1*9.)
    IF (Tlossrec .LE. LOG(250./8.)) THEN
       res(1) = -U4/(U1**2)
       res(4) = 1./U1
       res = 8.*EXP(Tlossrec)*((2.*T0)/(27.*phys%Mref))*res
    ENDIF
  END SUBROUTINE compute_dTlossrec_dU


  SUBROUTINE compute_fEiiz(U,fEiiz)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fEiiz,U3,U5
    REAL,PARAMETER :: tol = 1e-20
    U3 = U(3)
    U5 = U(5)
    IF (U3<tol) U3=tol
    IF (U5<tol) U5=tol
    fEiiz = U3*U5
  END SUBROUTINE compute_fEiiz


  SUBROUTINE compute_dfEiiz_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U3,U5
    REAL,PARAMETER :: tol = 1e-20
    U3 = U(3)
    U5 = U(5)
    IF (U3<tol) U3=tol
    IF (U5<tol) U5=tol
    res = 0.
    res(3) = U5
    res(5) = U3
  END SUBROUTINE compute_dfEiiz_dU


  SUBROUTINE compute_fEirec(U,fEirec)
   real*8, intent(IN) :: U(:)
   real*8             :: fEirec,U1,U3
   real,parameter :: tol = 1e-20
   U1 = U(1)
   U3 = U(3)
   if (U1<tol) U1=tol
   if (U3<tol) U3=tol
   fEirec = U1*U3
   !PSI review
   fEirec = fEirec-0.5*U(2)**2/U1*U3
 END SUBROUTINE compute_fEirec


 SUBROUTINE compute_dfEirec_dU(U,res)
   real*8, intent(IN) :: U(:)
   real*8             :: res(:),U1,U3
   real,parameter :: tol = 1e-20
   U1 = U(1)
   U3 = U(3)
   if (U1<tol) U1=tol
   if (U3<tol) U3=tol
   res = 0.
   res(1) = U3
   res(3) = U1
   !PSI review
   res(1) = res(1)+0.5*U(2)**2/U1**2*U3
   res(2) = -1.*U(2)/U1*U3
   res(3) = res(3) - 0.5*U(2)**2/U1
 END SUBROUTINE compute_dfEirec_dU


  SUBROUTINE compute_fEicx(U,fEicx)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: fEicx,U1,U2,U5
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    fEicx = (U5*U2**2)/U1
  END SUBROUTINE compute_fEicx


  SUBROUTINE compute_dfEicx_dU(U,res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: res(:),U1,U2,U5
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.
    res(1) = -U5*(U2/U1)**2
    res(2) = 2.*U5*U2/U1
    res(5) = (U2**2)/U1
  END SUBROUTINE compute_dfEicx_dU
#ifdef NEUTRAL
  !*******************************************
  ! Compute the terms relative to k equations
  !*******************************************
#ifdef KEQUATION

  SUBROUTINE compute_cs(U, cs)
    ! Sound speed of plasma
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: U1,U2,U3,U4
    REAL*8, INTENT(OUT) :: cs
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    cs = 2./3./U1*(U3+U4-1./2.*U2**2/U1)
    IF (cs<0.) cs = tol**2
    cs = SQRT(cs)
  END SUBROUTINE compute_cs
  SUBROUTINE compute_dcs_du(U, dcs_du)
    ! Sound speed derivative
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: U1,U2,U3,U4,cs
    REAL*8, INTENT(OUT) :: dcs_du(:)
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    IF (U4 < tol) U4 = tol
    IF (U1 < tol) U1 = tol
    IF (U3 < tol) U3 = tol
    dcs_du = 0.
    CALL compute_cs(U,cs)
    IF (cs>tol) THEN
       dcs_du(1) = -1.*(U3+U4-U2**2/U1)/U1**2
       dcs_du(2) = -1.*U2/U1**2
       dcs_du(3) = 1./U1
       dcs_du(4) = 1./U1


       dcs_du = dcs_du/3./cs
    ENDIF
  END SUBROUTINE compute_dcs_du
  SUBROUTINE compute_gamma_I(U,Q, Btor, gradBtor, R, gamma_I)
    ! growth rate for turbulent energy
    REAL*8, INTENT(IN) :: U(:), Q(:,:), gradBtor(:), Btor, R
    REAL*8             :: U1,U2,U3,U4, ti, te, gr_p_gr_b, cs, p, theta, ti_te
    REAL*8, INTENT(OUT) :: gamma_I
    REAL,PARAMETER :: tol = 1e-20
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
  END SUBROUTINE compute_gamma_I

  SUBROUTINE compute_ce(U,Q, Btor, gradBtor, r,omega_c,q_cyl, ce)
    ! dissipation rate for turbulent energy
    REAL*8, INTENT(IN) :: U(:), Q(:,:), gradBtor(:), r, Btor,omega_c,q_cyl
    REAL*8             :: U1,U2,U3, gamma_I, rhoL, cs, gamma_e , k_loc
    REAL*8, INTENT(OUT) :: ce
    REAL,PARAMETER :: tol = 1e-20
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)

    gamma_e = 4.5
    CALL compute_cs(U, cs)
    IF (cs < tol) cs = tol
    CALL compute_rhoL(U, r, omega_c, rhoL)
    IF (rhoL < tol) rhoL = tol
    CALL compute_gamma_I(U,Q,Btor,gradBtor,r,gamma_I)

    ce = gamma_I*(2.*PI/16./1./gamma_e/rhoL**2/q_cyl/cs**2+1./phys%k_max)

  END SUBROUTINE compute_ce
  SUBROUTINE compute_rhoL(U, R,omega_c, rhoL)
    ! Larmor radii
    REAL*8, INTENT(IN) :: U(:), R,omega_c
    REAL*8             :: cs
    REAL*8, INTENT(OUT) :: rhoL
    CALL compute_cs(U, cs)

    rhoL = cs/omega_c/R

  END SUBROUTINE compute_rhoL

  SUBROUTINE compute_dissip(U, dissip)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: U6
    REAL*8             :: dissip
    REAL,PARAMETER :: tol = 1e-10
    U6 = U(6)
    dissip = U6**2
  END SUBROUTINE  compute_dissip

  SUBROUTINE compute_ddissip_du(U, res)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8             :: U6
    REAL*8             :: res(:)
    REAL,PARAMETER :: tol = 1e-10
    U6 = U(6)
    !if (U6 < tol) U6 = tol
    res = 0.
    res(6) = 2.*U6
  END SUBROUTINE  compute_ddissip_du

#endif
#endif

#ifdef NEUTRALP
  SUBROUTINE computeDpn(U,Q,Vpn,Dpn)
    REAL*8, INTENT(IN) :: U(:),Q(:,:),Vpn(:)
    REAL*8, INTENT(OUT):: Dpn
    REAL*8             :: U1,U2,U3,U4,U5
    REAL*8             ::sigmaviz,sigmavcx,cs_n,Dpn_th
    REAL*8             :: Grad_Pn(simpar%Ndim)
    REAL,PARAMETER :: tol = 1e-10
    Dpn = 0.
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U2<tol) U2=tol
    IF (U3<tol) U3=tol
    IF (U4<tol) U4=tol
    IF (U5<tol) U5=tol
    CALL compute_sigmaviz(U,sigmaviz)
    CALL compute_sigmavcx(U,sigmavcx)
    Dpn = 1./(simpar%refval_time*simpar%refval_density*U1*(sigmaviz + sigmavcx))
    Dpn = 2./3.*Dpn
    !Set a threshold Dpn*|grad(Pn)| <= cs_n*n_n
    !cs_n = 0.
    !Grad_Pn = 0.
    !cs_n = sqrt(phys%Mref*abs(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    !Grad_Pn = 2./(3.*phys%Mref)*matmul(Q,Vpn)
    !Dpn_th = abs(cs_n*U5/norm2(Grad_Pn))
    !if (Dpn .gt. Dpn_th) then
    IF (Dpn .GT. phys%diff_nn) THEN
       Dpn = phys%diff_nn
       !fth(3) = 5.e-7*abs(U5)!*abs(U3/U1)
       !fth(4) = 5.e-7*abs(U5)!*abs(U4/U1)
    END IF
    !Dpn*exp(abs(Dpn - Dpn_th)/Dpn_th) + Dpn_th
  END SUBROUTINE computeDpn

  SUBROUTINE compute_dDpn_dU(U,Q,Vpn,res)
    REAL*8, INTENT(IN) :: U(:),Q(:,:),Vpn(:)
    REAL*8, INTENT(OUT):: res(:)
    REAL*8             :: sigmaviz,sigmavcx,Dpn
    REAL*8             :: dsigmaviz_dU(phys%Neq),dsigmavcx_dU(phys%Neq)
    REAL*8             :: U1,U2,U3,U4,U5
    REAL,PARAMETER     :: tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.

    CALL compute_sigmaviz(U,sigmaviz)
    CALL compute_sigmavcx(U,sigmavcx)
    CALL compute_dsigmaviz_dU(U,dsigmaviz_dU)
    CALL compute_dsigmavcx_dU(U,dsigmavcx_dU)
    CALL computeDpn(U,Q,Vpn,Dpn)

    IF (Dpn .LT. phys%diff_nn) THEN
       res(1) = sigmaviz + sigmavcx + U1*(dsigmaviz_dU(1) + dsigmavcx_dU(1))
       res(4) = U1*(dsigmaviz_dU(4) + dsigmavcx_dU(4))
       res = -3./2*Dpn**2*res
    END IF
  END SUBROUTINE compute_dDpn_dU

  SUBROUTINE computeVpn(U,Vpn)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(OUT):: Vpn(:)
    REAL*8             :: U1,U2,U3,U5
    REAL,PARAMETER :: tol = 1e-12
    Vpn = 0.
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    Vpn(1) = U5*(U2**2/U1**3 - U3/U1**2)
    Vpn(2) = - U5*U2/U1**2
    Vpn(3) = U5/U1
    Vpn(5) = U3/U1 - 1./2.*U2**2/U1**2
  END SUBROUTINE computeVpn

  SUBROUTINE compute_dVpn_dU(U,dVpn_dU)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: dVpn_dU(:, :)
    REAL*8              :: U1,U2,U3,U5
    REAL,PARAMETER :: tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    dVpn_dU = 0.

    dVpn_dU(1, 1) = (2.*U3/U1**3 - 3.*U2**2/U1**4)*U5
    dVpn_dU(1, 2) = (2.*U2/U1**3)*U5
    dVpn_dU(1, 3) = -U5/U1**2
    dVpn_dU(1, 5) = U2**2/U1**3 - U3/U1**2

    dVpn_dU(2, 1) = (2.*U2/U1**3)*U5
    dVpn_dU(2, 2) = - U5/U1**2
    dVpn_dU(2, 5) = - U2/U1**2

    dVpn_dU(3, 1) = -U5/U1**2
    dVpn_dU(3, 5) = 1./U1

    dVpn_dU(5, 1) = U2**2/U1**3 - U3/U1**2
    dVpn_dU(5, 2) = -U2/U1**2
    dVpn_dU(5, 3) = 1./U1
  END SUBROUTINE compute_dVpn_dU

  SUBROUTINE computeGammared(U,res)
    REAL*8, INTENT(IN)  :: U(:)
    REAL*8, INTENT(OUT) :: res
    REAL*8              :: U1,U2,U3,U5,Tmin,T
    REAL,PARAMETER      :: tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol

    res = 0.
    Tmin = 0.2/simpar%refval_temperature  ! Threshold at 0.2 eV
    T = 2./(3.*phys%Mref)*(U3/U1 - 1./2.*U2**2/U1**2)

    !WRITE(6,*) 'Tmin/T = ', Tmin/T

    IF (Tmin/T .LE. 1) THEN
       res = 0.*3./2.*(phys%Mref**2)*(Tmin/T)*U5
    ELSE
       res = 0.*3./2.*(phys%Mref**2)*U5
    END IF

  END SUBROUTINE computeGammared

  SUBROUTINE computeAlphaCoeff(U,Q,Vpn,res)
    REAL*8, INTENT(IN)    :: U(:),Q(:,:),Vpn(:)
    REAL*8, INTENT(OUT)   :: res
    REAL*8                :: U1,U2,U3,U4,U5
    REAL*8                :: cs_n,Dpn,t
    REAL*8                :: Grad_Pn(simpar%Ndim)
    REAL,PARAMETER        :: gamma = 4.,tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 1.

    cs_n = SQRT(phys%Mref*ABS(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    Grad_Pn = 2./(3.*phys%Mref)*MATMUL(Q,Vpn)
    CALL computeDpn(U,Q,Vpn,Dpn)
    !t = abs(Dpn*norm2(Grad_Pn))/abs(cs_n*U5)
    t = ABS(Dpn)/phys%diff_nn

    IF (t .GE. 1) THEN
       res = EXP((1-t)/gamma)
    END IF
  END SUBROUTINE computeAlphaCoeff

  SUBROUTINE computeBetaCoeff(U,Q,Vpn,res)
    REAL*8, INTENT(IN)    :: U(:),Q(:,:),Vpn(:)
    REAL*8, INTENT(OUT)   :: res
    REAL*8                :: U1,U2,U3,U4,U5
    REAL*8                :: cs_n,Dpn,t
    REAL*8                :: Grad_Pn(simpar%Ndim)
    REAL,PARAMETER        :: gamma = 4.,tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.

    cs_n = SQRT(phys%Mref*ABS(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    Grad_Pn = 2./(3.*phys%Mref)*MATMUL(Q,Vpn)
    CALL computeDpn(U,Q,Vpn,Dpn)
    !t = abs(Dpn*norm2(Grad_Pn))/abs(cs_n*U5)
    t = ABS(Dpn)/phys%diff_nn

    IF (t .GE. 1) THEN
       res = 1 - EXP((1-t)/gamma)
    END IF
  END SUBROUTINE computeBetaCoeff

  SUBROUTINE computeGammaLim(U,Q,Vpn,res)
    REAL*8, INTENT(IN)    :: U(:),Q(:,:),Vpn(:)
    REAL*8, INTENT(OUT)   :: res
    REAL*8                :: U1,U2,U3,U4,U5
    REAL*8                :: cs_n,Dpn,GammaDpn,GammaLim
    REAL*8                :: Grad_Pn(simpar%Ndim)
    REAL,PARAMETER        :: tol = 1e-12
    U1 = U(1)
    U2 = U(2)
    U3 = U(3)
    U4 = U(4)
    U5 = U(5)
    IF (U1<tol) U1=tol
    IF (U5<tol) U5=tol
    res = 0.

    cs_n = SQRT(phys%Mref*ABS(2./(3.*phys%Mref)*(U3/U1 - 1./2.*(U2/U1)**2)))
    Grad_Pn = 2./(3.*phys%Mref)*MATMUL(Q,Vpn)
    CALL computeDpn(U,Q,Vpn,Dpn)
    GammaDpn = ABS(Dpn*NORM2(Grad_Pn))
    GammaLim = ABS(cs_n*U5)

    res = 1./(1. + GammaDpn/GammaLim)
  END SUBROUTINE computeGammaLim
#endif
  !NEUTRAL PRESSURE

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
  SUBROUTINE computeTauGaussPoints(up, uc, b, n, iel, isext, xy, tau)
    REAL*8, INTENT(in)  :: up(:), uc(:), b(:), n(:), xy(:)
    REAL, INTENT(in)    :: isext
    INTEGER, INTENT(in) :: iel
    REAL*8, INTENT(out) :: tau(:, :)
#ifdef NEUTRAL
#ifndef KEQUATION
    REAL*8              :: tau_aux(5),diff_iso(5,5,1),diff_ani(5,5,1)
#else
    REAL*8              :: tau_aux(6),diff_iso(6,6,1),diff_ani(6,6,1)
#endif
#ifdef NEUTRALP
    REAL*8              :: Dpn
    REAL*8              :: Vpn(simpar%Neq),Qpr(simpar%Ndim,simpar%Neq)
#endif
#else
    REAL*8              :: tau_aux(4),diff_iso(4,4,1),diff_ani(4,4,1)
#endif
    INTEGER             :: ndim
    REAL*8              :: bn, bnorm,xyd(1,SIZE(xy)),uu(1,SIZE(uc))
    REAL*8              :: U1, U2, U3, U4
    U1 = uc(1)
    U2 = uc(2)
    U3 = uc(3)
    U4 = uc(4)

    tau = 0.
    ndim = SIZE(n)
    bn = dot_PRODUCT(b(1:ndim), n)
    bnorm = NORM2(b(1:ndim))
    xyd(1,:) = xy(:)
    uu(1,:) = uc(:)

    CALL setLocalDiff(xyd, uu, diff_iso, diff_ani)

#ifdef NEUTRALP
    ! Compute Vpn(U^(k-1))
    CALL computeVpn(uc,Vpn)
    ! Compute Dpn(U^(k-1))
    Qpr = RESHAPE(q,(/simpar%Ndim,simpar%Neq/))
    !CALL computeDpn(uc,Qpr,Vpn,Dpn)
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
          tau_aux(5) = tau_aux(5) + phys%diff_nn*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#endif
       ELSE
#endif
          ! Toroidal face
          tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(3) = tau_aux(3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(4) = tau_aux(4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
          tau_aux(5) = tau_aux(5) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
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
          tau_aux(5) = tau_aux(5) + phys%diff_nn*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#endif
       ELSE
#endif
          ! Toroidal face
          tau_aux(1) = tau_aux(1) + phys%diff_n*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(2) = tau_aux(2) + phys%diff_u*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(3) = tau_aux(3) + (phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(4) = tau_aux(4) + (phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1))*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifdef NEUTRAL
          tau_aux(5) = tau_aux(5) + phys%diff_nn*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#endif
#ifdef TOR3D
       ENDIF
#endif

    ELSEIF (numer%stab == 4) THEN
       tau_aux = MAX(ABS(5./3.*up(2)*bn), ABS(0.3*bn*(3*uc(1) + SQRT(ABS(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)), &
            phys%lscale/geom%R0*ABS(bn)*phys%diff_pari*up(7)**2.5, phys%lscale/geom%R0*ABS(bn)*phys%diff_pare*up(8)**2.5)

    ELSEIF (numer%stab == 5) THEN
       IF (ABS(isext - 1.) .LT. 1e-12) THEN
          ! exterior faces
          tau_aux = ABS((4*uc(2)*bn)/uc(1))
       ELSE
          tau_aux = MAX(ABS(5./3.*up(2)*bn), ABS(0.3*bn*(3*uc(1) + SQRT(ABS(10*uc(3)*uc(1) + 10*uc(4)*uc(1) - 5*uc(2)**2)))/uc(1)))
          !#ifdef NEUTRALP
          !        tau_aux(5) = max(abs(5./3.*up(2)*bn), abs(0.3*bn*(3*uc(5) + sqrt(abs(10*uc(3)/uc(1)*uc(5)**2 - 5*(uc(5)*uc(2)/uc(1))**2)))/uc(5)))
          !#endif
       ENDIF
#ifdef TOR3D
       IF (ABS(n(3)) > 0.1) THEN
          ! Poloidal face
          tau_aux(1) = tau_aux(1) + phys%diff_n
          tau_aux(2) = tau_aux(2) + phys%diff_u
          tau_aux(3) = tau_aux(3) + phys%diff_e + ABS(bn)*phys%diff_pari*up(7)**2.5*bnorm/uc(1)*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
          tau_aux(4) = tau_aux(4) + phys%diff_ee + ABS(bn)*phys%diff_pare*up(8)**2.5*bnorm/uc(1)*refElTor%Ndeg/(numer%tmax*xy(1)/numer%ntor)/phys%lscale
#ifndef NEUTRALP
#ifdef NEUTRAL
          tau_aux(5) = phys%diff_nn!numer%tau(5) !tau_aux(5) + diff_iso(5,5,1)
#endif
#else
          tau_aux(5) = tau_aux(5) + phys%diff_nn!phys%diff_nn
#endif
       ELSE
#endif
          ! Toroidal face
          tau_aux(1) = tau_aux(1) + 6*diff_iso(1,1,1)
          tau_aux(2) = tau_aux(2) + 6*diff_iso(2,2,1)
          tau_aux(3) = tau_aux(3) + 6*diff_iso(3,3,1) + ABS(bn)*phys%diff_pari*(MIN(1.,up(7)))**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
          tau_aux(4) = tau_aux(4) + 6*diff_iso(4,4,1) + ABS(bn)*phys%diff_pare*(MIN(1.,up(8)))**2.5*bnorm/uc(1)*refElPol%ndeg/Mesh%elemSize(iel)/phys%lscale
#ifndef NEUTRALP
#ifdef NEUTRAL
          tau_aux(5) = tau_aux(5) + phys%diff_nn !! !numer%tau(5) diff_iso(5,5,1)
#ifdef KEQUATION
          tau_aux(6) = tau_aux(6) + 6.*diff_iso(6,6,1)
#endif
#endif
#else
          tau_aux(5) = tau_aux(5) + numer%tau(5)
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
    tau(5, 5) = tau_aux(5)
#ifdef KEQUATION
    tau(6,6) = tau_aux(6)
#endif
#endif
  END SUBROUTINE computeTauGaussPoints

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

  END SUBROUTINE computeTauGaussPoints_matrix

  !  !***********************************************************************
  !  !
  !  !                           MAGNETIC FIELD
  !  !
  !  !***********************************************************************
  !  #ifdef TOR3D
  !  SUBROUTINE defineMagneticField(x,y,t,b,divb,drift)
  !    real*8, intent(in)      :: x(:),y(:),t(:)
  !    real*8, intent(out),optional     :::: b(:,:),divb(:),drift(:,:)
  !    real*8                  :: R0,q,r
  !    real*8                  :: xmax,xmin,ymax,ymin,xm,ym,p,divbp
  !    real*8                  :: xx,yy,tt,Bp,Bt,Br,Bz,BB,dmax,B0,xr,yr
  !    integer*4               :: i,j,ind,N2D,N1D
  !
  !
  !    N2d = size(X,1)
  !    N1d = size(t,1)
  !    xmax = Mesh%xmax
  !    xmin = Mesh%xmin
  !    ymax = Mesh%ymax
  !    ymin = Mesh%ymin
  !    xm = 0.5*(xmax+xmin)
  !    ym = 0.5*(ymax+ymin)
  !    !  xm = -0.5
  !    !  ym = -0.5
  !    ! Initialization
  !    ! Initialization
  !    if (present(b)) then
  !      b     = 0.
  !    endif
  !    if (present(divb)) then
  !      divb  = 0.
  !    endif
  !    if (present(drift)) then
  !      drift = 0.
  !    endif
  !    if (present(Bmod)) then
  !      Bmod  = 0.
  !    endif
  !
  !    DO i=1,N2d
  !      DO j=1,N1d
  !        xx = x(i)
  !        yy = y(i)
  !        tt = t(j)
  !        ind = (j-1)*N2d+i
  !
  !        SELECT CASE(switch%testcase)
  !        CASE(1)
  !          IF (switch%axisym) THEN
  !            WRITE(6,*) "This is NOT an axisymmetric test case!"
  !            stop
  !          END IF
  !          ! Cartesian case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
  !          Br = (yy-ym)
  !          Bz = (-xx+xm)
  !          Bt = 1.
  !          divbp = 0.
  !        CASE(2)
  !          IF (.not.switch%axisym) THEN
  !            WRITE(6,*) "This is an axisymmetric test case!"
  !            stop
  !          END IF
  !          ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
  !          Br = (yy-ym)/xx
  !          Bz = (-xx+xm)/xx
  !          Bt = 1.
  !          divbp = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
  !        CASE(3)
  !          IF (.not.switch%axisym) THEN
  !            WRITE(6,*) "This is an axisymmetric test case!"
  !            stop
  !          END IF
  !          ! Axysimmetric case, circular field centered in [xm, ym] in the poloidal plane, Bt = 1
  !          Br = (yy-ym)/xx
  !          Bz = (-xx+xm)/xx
  !          Bt = 1.
  !          divbp = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
  !
  !        CASE(50:59)
  !          write(6,*) "Error in defineMagneticField: you should not be here!"
  !          STOP
  !        CASE(60:69)
  !
  !          ! Circular case with limiter
  !          R0 = geom%R0
  !          q  = geom%q
  !          B0 = 2 ! not influential
  !          xr = xx*phys%lscale
  !          yr = yy*phys%lscale
  !
  !          r  = sqrt((xr-R0)**2+yr**2)
  !          Br = -B0*yr/(xr*q*sqrt(1- (r/R0)**2 ) )
  !          Bz = B0*(xr-R0)/(xr*q*sqrt(1- (r/R0)**2 ) )
  !          Bt = B0*R0/xr
  !
  !          if (present(divb)) then
  !            IF (switch%axisym) THEN
  !              divbp = -yy/xx/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
  !            ELSE
  !              WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
  !              STOP
  !            END IF
  !          endif
  !          if (present(divb)) then
  !            IF (switch%driftdia) THEN
  !              drift(:,2) =  -1./R0*phys%lscale
  !            END IF
  !          endif
  !
  !
  !
  !
  !
  !
  !        CASE DEFAULT
  !          WRITE(6,*) "Error! Test case not valid"
  !          STOP
  !        END SELECT
  !
  !        Bp = sqrt(Br**2+Bz**2)
  !        BB = sqrt(Bp**2+Bt**2)
  !        b(ind,1) = Br/BB
  !        b(ind,2) = Bz/BB
  !        b(ind,3) = Bt/BB
  !        if (present(divb)) then
  !          divb(ind) = divbp
  !        endif
  !      END DO
  !    END DO
  !  END SUBROUTINE defineMagneticField
  !  #else
  !  SUBROUTINE defineMagneticField(x,y,b,divb,drift,Bmod)
  !    real*8, intent(in)      :: x(:),y(:)
  !    real*8, intent(out)     :: b(:,:)
  !    real*8, intent(out),optional::divb(:),drift(:,:),Bmod(:)
  !    real*8                  :: R0,q,r(size(x)),auxdivb(size(x)),auxdrift(size(b,1),size(b,2))
  !    real*8                  :: xmax,xmin,ymax,ymin,xm,ym,xx,yy
  !    real*8                  :: Br,Bz,Bt,BB,Bp
  !    integer :: i
  !    ! Initialization
  !    b     = 0.
  !    auxdivb(:) = 0.
  !    xmax = Mesh%xmax
  !    xmin = Mesh%xmin
  !    ymax = Mesh%ymax
  !    ymin = Mesh%ymin
  !    xm = 0.5*(xmax+xmin)
  !    ym = 0.5*(ymax+ymin)
  !    if (present(Bmod)) then
  !      Bmod  = 0.
  !    endif
  !    !  xm = -0.5
  !    !  ym = -0.5
  !    SELECT CASE(switch%testcase)
  !    CASE(1)
  !      ! Cartesian case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 20+sin(wx*x)*cos(wy*y)
  !      !            b(:,1) = 1./30.*(x-y**2+2)
  !      !            b(:,2) = 1./30.*(x*y+y)
  !      !            auxdivb(:) = 1./15.+(1./30.)*x
  !      DO i=1,size(x)
  !        xx = x(i)
  !        yy = y(i)
  !        Br = (yy-ym)
  !        Bz = (-xx+xm)
  !        Bt = 1.
  !        Bp = sqrt(Br**2+Bz**2)
  !        BB = sqrt(Bp**2+Bt**2)
  !        b(i,1) = Br/BB
  !        b(i,2) = Bz/BB
  !        auxdivb(i) = 0.
  !      END DO
  !
  !    CASE(2)
  !      ! Axisymmetric case with div(b)~=0, n = 2+sin(wx*x )*sin(wy*y),  u = cos(wx*x)*cos(wy*y), Ei = 20+cos(wx*x)*sin(wy*y), Ee = 20+sin(wx*x)*cos(wy*y)
  !      !            b(:,1) = 1./30.*(x-y**2+2)
  !      !            b(:,2) = 1./30.*(x*y+y)
  !      !            auxdivb(:) = 1./30.+((1./30.)*x-(1./30.)*y**2+1./15.)/x+1./30.*(x+1)
  !      DO i=1,size(x)
  !        xx = x(i)
  !        yy = y(i)
  !        Br = (yy-ym)/xx
  !        Bz = (-xx+xm)/xx
  !        Bt = 1.
  !        Bp = sqrt(Br**2+Bz**2)
  !        BB = sqrt(Bp**2+Bt**2)
  !        b(i,1) = Br/BB
  !        b(i,2) = Bz/BB
  !        auxdivb(i) = (xx**2*yy-xx**2*ym+xm**2*yy-xm**2*ym+3*yy*ym**2-3*yy**2*ym+yy**3-ym**3-2*xx*xm*yy+2*xx*xm*ym)/(xx**4*((xx-xm)**2/xx**2+(yy-ym)**2/xx**2+1)**(1.5))
  !      END DO
  !    CASE(5)
  !      ! Cartesian case, square mesh, horizontal field
  !      b(:,1) = 0.1
  !      b(:,2) = 0.
  !    CASE(6)
  !      ! Cartesian case, square mesh, horizontal field
  !      DO i=1,size(x)
  !        IF (y(i).ge.0.5) THEN
  !          b(i,1) = 0.1
  !        ELSE
  !          b(i,1) = -0.1
  !        END IF
  !      END DO
  !    CASE(50:59)
  !      write(6,*) "Error in defineMagneticField: you should not be here!"
  !      STOP
  !    CASE(60:69)
  !
  !      ! Circular case with limiter
  !      R0 = geom%R0
  !      q  = geom%q
  !      r  = phys%lscale*sqrt((x-R0/phys%lscale)**2+y**2)
  !      b(:,1) = -phys%lscale*y/sqrt(R0**2*q**2+(1-q**2)*r**2)
  !      b(:,2) = phys%lscale*(x-R0/phys%lscale)/sqrt(R0**2*q**2+(1-q**2)*r**2)
  !
  !      IF (switch%axisym) THEN
  !        auxdivb(:) = -y/x/sqrt(R0**2*q**2+(1-q**2)*r**2)*phys%lscale
  !      ELSE
  !        WRITE(6,*) "Not coded: usually here you should have an axisym simulation"
  !        STOP
  !      END IF
  !
  !      IF (switch%driftdia) THEN
  !        auxdrift(:,2) =  -1./R0*phys%lscale
  !      END IF
  !    CASE DEFAULT
  !      WRITE(6,*) "Error! Test case not valid"
  !      STOP
  !    END SELECT
  !
  !    IF (present(divb)) THEN
  !      divb = auxdivb
  !    ENDIF
  !    IF (present(drift)) THEN
  !      drift = auxdrift
  !    ENDIF
  !  END SUBROUTINE defineMagneticField
  !  #endif
  !
  !  SUBROUTINE loadMagneticField()
  !    USE interpolation
  !    USE HDF5
  !    USE HDF5_io_module
  !    integer        :: i,ierr,ip,jp
  !    integer(HID_T) :: file_id
  !    real*8,pointer,dimension(:,:) :: r2D,z2D,flux2D,Br2D,Bz2D,Bphi2D
  !    real*8,allocatable,dimension(:,:) :: bx,by,bmod,divb,bmodx,bmody,driftx,drifty
  !    real*8,allocatable,dimension(:)   :: xvec,yvec
  !    real*8                            :: x,y
  !
  !    WRITE(6,*) "******* Loading magnetic field *******"
  !    ! Allocate storing space in phys
  !    ALLOCATE(phys%b(Mesh%Nnodes,Mesh%Ndim))
  !    ALLOCATE(phys%divb(Mesh%Nnodes))
  !    ALLOCATE(phys%drift(Mesh%Nnodes,Mesh%Ndim))
  !    ALLOCATE(phys%flux2d(Mesh%Nnodes))
  !
  !    ! Dimensions of the file storing the magnetic field for West
  !    ip = 541
  !    jp = 391
  !    ALLOCATE(r2D(ip,jp))
  !    ALLOCATE(z2D(ip,jp))
  !    ALLOCATE(flux2D(ip,jp))
  !    ALLOCATE(Br2D(ip,jp))
  !    ALLOCATE(Bz2D(ip,jp))
  !    ALLOCATE(Bphi2D(ip,jp))
  !    ALLOCATE(bx(ip,jp))
  !    ALLOCATE(by(ip,jp))
  !    ALLOCATE(bmod(ip,jp))
  !    ALLOCATE(bmodx(ip,jp))
  !    ALLOCATE(bmody(ip,jp))
  !    ALLOCATE(divb(ip,jp))
  !    ALLOCATE(driftx(ip,jp))
  !    ALLOCATE(drifty(ip,jp))
  !
  !    ! Read file
  !    CALL HDF5_open('WEST_far_465.h5',file_id,IERR)
  !    CALL HDF5_array2D_reading(file_id,r2D,'r2D')
  !    CALL HDF5_array2D_reading(file_id,z2D,'z2D')
  !    CALL HDF5_array2D_reading(file_id,flux2D,'flux2D')
  !    CALL HDF5_array2D_reading(file_id,Br2D,'Br2D')
  !    CALL HDF5_array2D_reading(file_id,Bz2D,'Bz2D')
  !    CALL HDF5_array2D_reading(file_id,Bphi2D,'Bphi2D')
  !    CALL HDF5_close(file_id)
  !
  !    ! Apply length scale
  !    r2D = r2D/phys%lscale
  !    z2D = z2D/phys%lscale
  !
  !    ! Compute b
  !    bmod = sqrt(Br2D**2+Bz2D**2+Bphi2D**2)
  !    bx = -Br2D/bmod
  !    by = -Bz2D/bmod
  !
  !    ! Compute divergence of b
  !    divb = 0.
  !    IF (switch%axisym) THEN
  !      ! 1/r*(d r*br/dr)+dbz/dz
  !      divb(2:ip-1,2:jp-1) = 1./r2D(2:ip-1,2:jp-1)*(r2D(2:ip-1,3:jp)*bx(2:ip-1,3:jp)- &
  !        r2D(2:ip-1,1:jp-2)*bx(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-  &
  !        r2D(2:ip-1,1:jp-2))+(by(3:ip,2:jp-1)-by(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
  !
  !    ELSE
  !      ! dbr/dr+dbz/dz
  !      divb(2:ip-1,2:jp-1) = (bx(2:ip-1,3:jp-1)-bx(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-r2D(2:ip-1,1:jp-2))+ &
  !        (by(3:ip,2:jp-1)-by(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
  !    END IF
  !
  !    ! Compute drift velocity
  !    driftx = 0.
  !    drifty = 0.
  !    IF (switch%driftdia) THEN
  !      bmodx = (bmod(2:ip-1,3:jp)-bmod(2:ip-1,1:jp-2))/(r2D(2:ip-1,3:jp)-r2D(2:ip-1,1:jp-2))
  !      bmody = (bmod(3:ip,2:jp-1)-bmod(1:ip-2,2:jp-1))/(z2D(3:ip,2:jp-1)-z2D(1:ip-2,2:jp-1))
  !      driftx(2:ip-1,2:jp-1) =  -Bphi2D(2:ip-1,2:jp-1)*bmody/bmod(2:ip-1,2:jp-1)**3
  !      drifty(2:ip-1,2:jp-1) =   Bphi2D(2:ip-1,2:jp-1)*bmodx/bmod(2:ip-1,2:jp-1)**3
  !    END IF
  !
  !    ! Interpolate
  !    ALLOCATE(xvec(jp))
  !    ALLOCATE(yvec(ip))
  !    xvec = r2D(1,:)
  !    yvec = z2D(:,1)
  !    DO i = 1,Mesh%Nnodes
  !      x = Mesh%X(i,1)
  !      y = Mesh%X(i,2)
  !      phys%b(i,1) = interpolate( ip, yvec,jp, xvec, bx, y,x, 1e-12)
  !      phys%b(i,2) = interpolate( ip, yvec,jp, xvec, by, y,x, 1e-12)
  !      phys%divb(i) = interpolate( ip, yvec,jp, xvec, divb, y,x, 1e-12)
  !      phys%drift(i,1) = interpolate( ip, yvec,jp, xvec,driftx, y,x, 1e-12)
  !      phys%drift(i,2) = interpolate( ip, yvec,jp, xvec,drifty, y,x, 1e-12)
  !      phys%flux2D(i) = interpolate( ip, yvec,jp, xvec,flux2D, y,x, 1e-12)
  !    END DO
  !
  !    ! Free memory
  !    DEALLOCATE(Br2D,Bz2D,Bphi2D,xvec,yvec)
  !    DEALLOCATE(r2D,z2D,flux2D,bx,by,bmod,bmodx,bmody,divb,driftx,drifty)
  !
  !  END SUBROUTINE loadMagneticField
  !
  !
  !
  !  SUBROUTINE loadMagneticFieldTemporalEvolution()
  !    USE HDF5
  !    USE HDF5_io_module
  !    USE MPI_OMP
  !    integer        :: ierr,k
  !    character(LEN=20) :: fname = 'Evolving_equilibrium'
  !    character(10)  :: npr,nid,nit
  !    character(len=1000) :: fname_complete
  !    integer(HID_T) :: file_id
  !
  !    WRITE(6,*) "******* Loading magnetic field *******"
  !
  !    ! Allocate storing space in phys
  !    ALLOCATE(phys%Br(Mesh%Nnodes))
  !    ALLOCATE(phys%Bz(Mesh%Nnodes))
  !    ALLOCATE(phys%Bt(Mesh%Nnodes))
  !    ALLOCATE(phys%flux2d(Mesh%Nnodes))
  !
  !    ! File name
  !    write(nit, "(i10)") time%it
  !    nit = trim(adjustl(nit))
  !    k = INDEX(nit, " ") -1
  !
  !    IF (MPIvar%glob_size.GT.1) THEN
  !      write(nid,*) MPIvar%glob_id+1
  !      write(npr,*) MPIvar%glob_size
  !      fname_complete = trim(adjustl(fname))//'_'//trim(adjustl(nid))//'_'//trim(adjustl(npr))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
  !    ELSE
  !      fname_complete = trim(adjustl(fname))//'_'//REPEAT("0", 4 - k)//trim(ADJUSTL(nit))//'.h5'
  !    END IF
  !
  !    write(6,*) 'Magnetic field loaded from file: ', trim(adjustl(fname_complete))
  !
  !    ! Read file
  !    CALL HDF5_open(fname_complete,file_id,IERR)
  !    CALL HDF5_array1D_reading(file_id,phys%Br,'Br')
  !    CALL HDF5_array1D_reading(file_id,phys%Bz,'Bz')
  !    CALL HDF5_array1D_reading(file_id,phys%Bt,'Bt')
  !    CALL HDF5_array1D_reading(file_id,phys%flux2D,'flux')
  !    CALL HDF5_close(file_id)
  !
  !  END SUBROUTINE loadMagneticFieldTemporalEvolution


END MODULE physics
