!**********************************
! project: MHDG
! file: READ_input.f90
! date: 04/09/2016
! Subroutine for input file loading
!**********************************

!********************************
! Loads input file
!********************************
SUBROUTINE READ_input()
  USE prec_const
  USE globals
  USE MPI_OMP
  IMPLICIT NONE

  LOGICAL               :: driftdia,driftexb, axisym, steady,dotiming,psdtime,decoup,bxgradb, read_gmsh,readMeshFromSol, set_2d_order, gmsh2h5,igz, adaptivity, time_adapt, NR_adapt, div_adapt, rest_adapt,osc_adapt
  LOGICAL               :: ckeramp,saveNR,filter,saveTau,transport_1d,lstiming,fixdPotLim,dirivortcore,dirivortlim,convvort,logrho
  INTEGER               :: thresh, difcor, tis, stab,pertini,init,order_2d
  INTEGER               :: itmax, itrace, rest, istop, sollib, kspitrace,rprecond, Nrprecond, kspitmax, kspnorm, gmresres,mglevels,mgtypeform
  INTEGER               :: uinput, printint, testcase, nrp
  INTEGER               :: nts, tsw, freqdisp, freqsave, shockcp, limrho
  INTEGER               :: shockcp_adapt, evaluator, difference, freq_t_adapt,freq_NR_adapt, quant_ind
  INTEGER,ALLOCATABLE,DIMENSION(:) :: n_quant_ind,param_est
  INTEGER               :: num_param_est, num_n_quant_ind
  REAL*8                :: thr_ind, tol_est, osc_tol, osc_check
  INTEGER               :: bcflags(1:10), ntor, ptor, npartor,bohmtypebc
  REAL*8                :: dt0, R0, diff_n, diff_u, v_p, tau(1:5), tNr, tTM, div, Tbg
  REAL*8                :: tfi, a, bohmth,bohm_energy_thresh, q, diffred, diffmin
  REAL*8                :: sc_coe, so_coe, df_coe, thr, thrpre, minrho, dc_coe, sc_sen
  REAL*8                :: epn, Mref, diff_pari, diff_e, Gmbohm, Gmbohme
  REAL*8                :: diff_pare, diff_ee, tie, dumpnr_min,dumpnr_max,dumpnr_width,dumpnr_n0, tmax, tol, rtol, atol
  REAL*8                :: diff_vort, diff_pot, etapar, Potfloat,diagsource(10)
  CHARACTER(100)        :: msg
  CHARACTER(20)         :: kmethd, ptype, kspmethd, pctype

  CHARACTER(len=20)     :: smther, smther2, prol, restr, solve, restr2, prol2, solve2, mlcycle
  CHARACTER(len=20)     :: aggr_prol, par_aggr_alg, aggr_ord, aggr_filter, csolve, csbsolve, cmat
  INTEGER               :: jsweeps, novr, fill, jsweeps2, novr2, fill2, outer_sweeps, maxlevs, csize, cfill, cjswp
  REAL*8                :: thrsol, thrsol2, mncrratio, athres, cthres
  REAL*8                :: heating_power, heating_dr,heating_dz,heating_sigmar,heating_sigmaz
  INTEGER               :: heating_equation
  REAL*8                :: exbdump, part_source,ener_source, density_source, ener_source_e, ener_source_ee, sigma_source, fluxg_trunc

  ! Info for input and output
  CHARACTER(len = 1000) :: field_path, jtor_path,save_folder, geometry_path,puff_path, target_density_path,target_density_xpr_path, impurity_concentration_path, zeff_path
  CHARACTER(len = 1000) :: transport_model_path = 'transport_model.nml'
  INTEGER               :: field_dimensions(1:2), jtor_dimensions(1:2),puff_dimension, target_density_dimension,target_density_xpr_dimension,impurity_concentration_dimension, zeff_dimension
  LOGICAL               :: field_from_grid, compute_from_flux, divide_by_2pi

  ! RMP and Ripple
  LOGICAL               :: RMP, Ripple
  REAL*8                :: amp_rmp, torElongCoils_rmp, amp_ripple, triang, ellip
  INTEGER               :: nbCoils_rmp, parite, nbRow, nbCoils_ripple

  ! Neutral and Ohmic heating
  LOGICAL               :: OhmicSrc, apply_trim
  REAL*8                :: Zeff,Pohmic,diff_nn,Re,Re_pump,puff,feedback_propotional_gain,feedback_integral_gain,feedback_derivative_gain,cryopump_power,puff_slope
  REAL*8                :: feedback_propotional_gain_xpr, feedback_integral_gain_xpr, feedback_derivative_gain_xpr
#ifdef KEQUATION
  ! k equation
  REAL*8                :: diff_k_min, diff_k_max, k_max
#endif
  ! Movin Equilibrium
  LOGICAL               :: ME
  LOGICAL               :: diff_reverse_Ip
  REAL*8                :: I_0


  ! impurity radiation
  LOGICAL               :: impurity_radiation
  CHARACTER(1000)       :: impurity_name
  REAL*8                :: impurity_concentration

  ! target density
  INTEGER               :: target_variable
  !external heating
  LOGICAL               :: external_heating, external_heating_from_grid
  CHARACTER(1000)       :: external_heating_path

  ! flux limiter
  LOGICAL               :: flux_limiter
  REAL*8                :: c_fli, c_fle
  REAL*8                :: T_fluxlim_maxi, T_fluxlim_maxe

  ! 1D diffusion
  LOGICAL               :: import_diffusion_1D
  CHARACTER(1000)       :: diffusion_1D_path

  !bohm-gyrobohm
  LOGICAL               :: bohm_gyrobohm

  ! Defining the variables to READ from the file
  NAMELIST /SWITCH_LST/ steady,read_gmsh, readMeshFromSol, set_2d_order, order_2d, gmsh2h5, axisym,external_heating, impurity_radiation, init, driftdia, driftexb, testcase, OhmicSrc, ME,diff_reverse_Ip, target_variable, RMP, Ripple, psdtime, diffred, diffmin, &
       & shockcp, limrho, difcor, thresh, filter, decoup, ckeramp, saveNR, saveTau, transport_1d, fixdPotLim, dirivortcore,dirivortlim, convvort,pertini,&
       & logrho,bxgradb,flux_limiter,import_diffusion_1D,bohm_gyrobohm
  NAMELIST /INPUT_LST/ field_path, field_dimensions,field_from_grid,compute_from_flux,divide_by_2pi, jtor_path, jtor_dimensions,external_heating_path,external_heating_from_grid, save_folder,puff_path,puff_dimension,target_density_path,target_density_dimension,target_density_xpr_path,target_density_xpr_dimension,impurity_concentration_path,impurity_concentration_dimension,zeff_path,zeff_dimension, diffusion_1D_path, transport_model_path
  NAMELIST /NUMER_LST/ tau,nrp,tNR,tTM,div,sc_coe,sc_sen,minrho,so_coe,df_coe,dc_coe,thr,thrpre,stab,dumpnr_min,dumpnr_max,dumpnr_width,dumpnr_n0,ntor,ptor,tmax,npartor,bohmtypebc,exbdump
  NAMELIST /ADAPT_LST/ adaptivity,shockcp_adapt, evaluator, param_est, thr_ind, quant_ind, n_quant_ind,tol_est, difference, time_adapt, NR_adapt, freq_t_adapt, freq_NR_adapt, div_adapt, rest_adapt, osc_adapt, osc_tol, osc_check, geometry_path
  NAMELIST /GEOM_LST/ R0, q
  NAMELIST /MAGN_LST/ amp_rmp,nbCoils_rmp,torElongCoils_rmp,parite,nbRow,amp_ripple,nbCoils_ripple,triang,ellip ! RMP and Ripple
  NAMELIST /TIME_LST/ dt0, nts, tfi, tsw, tis
#ifndef KEQUATION
  NAMELIST /PHYS_LST/ diff_n, diff_u, diff_e, diff_ee, diff_vort, v_p, diff_nn,I_0, heating_power, heating_dr,heating_dz,heating_sigmar,heating_sigmaz,heating_equation,&
  & Re, Re_pump, apply_trim, puff,impurity_name,impurity_concentration,feedback_propotional_gain,feedback_integral_gain,feedback_derivative_gain,& 
  & feedback_propotional_gain_xpr, feedback_integral_gain_xpr, feedback_derivative_gain_xpr, cryopump_power,puff_slope, density_source, ener_source_e, ener_source_ee, sigma_source, fluxg_trunc, part_source,ener_source,Zeff, Pohmic, Tbg, bcflags, bohmth,&
    &bohm_energy_thresh,Gmbohm, Gmbohme, a, Mref, tie, diff_pari, diff_pare, diff_pot, epn, etapar, Potfloat,diagsource, c_fli, c_fle, T_fluxlim_maxi, T_fluxlim_maxe
#else
  NAMELIST /PHYS_LST/ diff_n, diff_u, diff_e, diff_ee, diff_vort, v_p, diff_nn,I_0,heating_power, heating_dr,heating_dz,heating_sigmar,heating_sigmaz,heating_equation, Re, Re_pump, apply_trim, puff,impurity_name,impurity_concentration,feedback_propotional_gain,feedback_integral_gain,feedback_derivative_gain,feedback_propotional_gain_xpr, feedback_integral_gain_xpr, feedback_derivative_gain_xpr,cryopump_power,puff_slope, density_source, ener_source_e, ener_source_ee, sigma_source, fluxg_trunc, part_source,ener_source,&
  & diff_k_min, diff_k_max, k_max, Zeff,Pohmic, Tbg, bcflags, bohmth,&
    &bohm_energy_thresh,Gmbohm, Gmbohme, a, Mref, tie, diff_pari, diff_pare, diff_pot, epn, etapar, Potfloat,diagsource, c_fli, c_fle, T_fluxlim_maxi, T_fluxlim_maxe
#endif
  NAMELIST /UTILS_LST/ PRINTint, dotiming, freqdisp, freqsave
  NAMELIST /LSSOLV_LST/ sollib, lstiming, kspitrace, rtol, atol, kspitmax, igz, rprecond,Nrprecond, kspnorm, kspmethd, pctype, gmresres,mglevels, mgtypeform,itmax, itrace, rest, istop, tol, kmethd, ptype,&
       &smther, jsweeps,&
       &novr, restr, prol, solve, fill, thrsol, smther2, jsweeps2, novr2, restr2, prol2, solve2, fill2, thrsol2, mlcycle,&
       &outer_sweeps, maxlevs, csize, aggr_prol, par_aggr_alg, aggr_ord, aggr_filter, mncrratio, athres,&
       &csolve, csbsolve, cmat, cfill, cthres, cjswp

  ! preallocate adaptivity arrays after specification statements
  ALLOCATE(param_est(1000))
  ALLOCATE(n_quant_ind(1000))
  param_est = -1
  n_quant_ind = -1

  ! Reading the file
  uinput = 100
  diagsource = 0.
  OPEN (uinput, file='param.txt', status='unknown')
  READ (uinput, SWITCH_LST)
  READ (uinput, INPUT_LST)
  READ (uinput, NUMER_LST)
  READ (uinput, ADAPT_LST)
  READ (uinput, GEOM_LST)
  READ (uinput, MAGN_LST)
  READ (uinput, TIME_LST)
  READ (uinput, PHYS_LST)
  READ (uinput, UTILS_LST)
  READ (uinput, LSSOLV_LST)
  CLOSE (uinput)

  IF (impurity_radiation) THEN
     IF ((TRIM(ADJUSTL(impurity_name)) .NE. 'N') ) THEN
         IF ((TRIM(ADJUSTL(impurity_name)) .NE. 'W')) THEN
            PRINT *, 'Only nitrogen or tungsten is allowed for impurity radiation so far. Stopping'
            STOP
         ENDIF
     ENDIF
  ENDIF

  ! Storing at the right place
  switch%steady           = steady
  switch%read_gmsh        = read_gmsh
  switch%readMeshFromSol  = readMeshFromSol
  switch%set_2d_order     = set_2d_order
  switch%order_2d         = order_2d
  switch%axisym           = axisym
  switch%init             = init
  switch%driftdia         = driftdia
  switch%driftexb         = driftexb
  switch%testcase         = testcase
  switch%ohmicsrc         = OhmicSrc
  switch%ME               = ME
  switch%diff_reverse_Ip = diff_reverse_Ip
  switch%target_variable   = target_variable
  switch%RMP              = RMP
  switch%Ripple           = Ripple
  switch%psdtime          = psdtime
  switch%diffred          = diffred
  switch%diffmin          = diffmin
  switch%shockcp          = shockcp
  switch%limrho           = limrho
  switch%difcor           = difcor
  switch%thresh           = thresh
  switch%filter           = filter
  switch%decoup           = decoup
  switch%ckeramp          = ckeramp
  switch%saveNR           = saveNR
  switch%saveTau          = saveTau
  switch%transport_1d = transport_1d
  switch%gmsh2h5          = gmsh2h5
  switch%fixdPotLim       = fixdPotLim
  switch%dirivortcore     = dirivortcore
  switch%dirivortlim      = dirivortlim
  switch%convvort         = convvort
  switch%pertini          = pertini
  switch%logrho           = logrho
  switch%flux_limiter     = flux_limiter
  switch%bxgradb          = bxgradb
  switch%external_heating = external_heating
  switch%impurity_radiation = impurity_radiation
  switch%import_diffusion_1D = import_diffusion_1D
  switch%bohm_gyrobohm    = bohm_gyrobohm
  input%field_path        = TRIM(ADJUSTL(field_path))
  input%field_dimensions  = field_dimensions
  input%field_from_grid   = field_from_grid
  input%compute_from_flux = compute_from_flux
  input%divide_by_2pi     = divide_by_2pi
  input%jtor_path         = TRIM(ADJUSTL(jtor_path))
  input%jtor_dimensions   = jtor_dimensions
  input%external_heating_path = TRIM(ADJUSTL(external_heating_path))
  input%external_heating_from_grid = external_heating_from_grid
  input%save_folder       = TRIM(ADJUSTL(save_folder))
  input%puff_path        = TRIM(ADJUSTL(puff_path))
  input%target_density_path = TRIM(ADJUSTL(target_density_path))
  input%target_density_xpr_path = TRIM(ADJUSTL(target_density_xpr_path))
  input%target_density_dimension = target_density_dimension
  input%target_density_xpr_dimension = target_density_xpr_dimension
  input%impurity_concentration_path = TRIM(ADJUSTL(impurity_concentration_path))
  input%impurity_concentration_dimension = impurity_concentration_dimension
  input%zeff_path         = TRIM(ADJUSTL(zeff_path))
  input%zeff_dimension    = zeff_dimension
  input%puff_dimension   = puff_dimension
  input%diffusion_1D_path = TRIM(ADJUSTL(diffusion_1D_path))
  input%transport_model_path = TRIM(ADJUSTL(transport_model_path))
  numer%tau               = tau
  numer%nrp               = nrp
  numer%tNR               = tNR
  numer%tTM               = tTM
  numer%div               = div
  numer%sc_coe            = sc_coe
  numer%sc_sen            = sc_sen
  numer%minrho            = minrho
  numer%so_coe            = so_coe
  numer%df_coe            = df_coe
  numer%dc_coe            = dc_coe
  numer%thr               = thr
  numer%thrpre            = thrpre
  numer%stab              = stab
  numer%dumpnr            = dumpnr_min
  numer%dumpnr_min        = dumpnr_min
  numer%dumpnr_max        = dumpnr_max
  numer%dumpnr_width      = dumpnr_width
  numer%dumpnr_n0         = dumpnr_n0
  numer%ntor              = ntor
  numer%ptor              = ptor
  numer%tmax              = tmax
  numer%npartor           = npartor
  numer%bohmtypebc        = bohmtypebc
  numer%exbdump           = exbdump
  adapt%adaptivity        = adaptivity
  adapt%shockcp_adapt     = shockcp_adapt
  num_param_est        = COUNT(param_est /= -1.0)
  ALLOCATE(adapt%param_est(num_param_est))
  adapt%param_est         = param_est(1:num_param_est)
  adapt%evaluator         = evaluator
  adapt%difference        = difference
  adapt%thr_ind           = thr_ind
  adapt%quant_ind         = quant_ind
  num_n_quant_ind    = COUNT(n_quant_ind /= -1.0)
  ALLOCATE(adapt%n_quant_ind(num_n_quant_ind))
  adapt%n_quant_ind       = n_quant_ind(1:num_n_quant_ind)
  adapt%tol_est           = tol_est
  adapt%time_adapt        = time_adapt
  adapt%NR_adapt          = NR_adapt
  adapt%freq_t_adapt      = freq_t_adapt
  adapt%freq_NR_adapt     = freq_NR_adapt
  adapt%div_adapt         = div_adapt
  adapt%rest_adapt        = rest_adapt
  adapt%osc_adapt         = osc_adapt
  adapt%osc_tol           = osc_tol
  adapt%osc_check         = osc_check

  IF (switch%transport_1d) CALL read_transport_model_input()
  adapt%geometry_path     = TRIM(ADJUSTL(geometry_path))
  geom%R0                 = R0
  geom%q                  = q
  magn%amp_rmp            = amp_rmp
  magn%nbCoils_rmp        = nbCoils_rmp
  magn%torElongCoils_rmp  = torElongCoils_rmp
  magn%parite             = parite
  magn%nbRow              = nbRow
  magn%amp_ripple         = amp_ripple
  magn%nbCoils_ripple     = nbCoils_ripple
  magn%triang             = triang
  magn%ellip              = ellip
  time%dt0                = dt0
  time%nts                = nts
  time%tfi                = tfi
  time%tsw                = tsw
  time%tis                = tis
  phys%diff_n             = diff_n
  phys%diff_u             = diff_u
  phys%diff_e             = diff_e
  phys%diff_ee            = diff_ee
  phys%diff_vort          = diff_vort
  phys%v_p                = v_p
  phys%diff_nn            = diff_nn
  phys%I_0                = I_0
  phys%heating_power      = heating_power
  phys%heating_dr         = heating_dr
  phys%heating_dz         = heating_dz
  phys%heating_sigmar     = heating_sigmar
  phys%heating_sigmaz     = heating_sigmaz
  phys%heating_equation   = heating_equation
  phys%Re                 = Re
  phys%Re_pump            = Re_pump
  phys%apply_trim         = apply_trim
  phys%cryopump_power     = cryopump_power
  phys%puff               = puff
  phys%impurity_name      = TRIM(ADJUSTL(impurity_name))
  phys%impurity_concentration = impurity_concentration
  phys%feedback_propotional_gain = feedback_propotional_gain
  phys%feedback_integral_gain = feedback_integral_gain
  phys%feedback_derivative_gain = feedback_derivative_gain
  phys%feedback_propotional_gain_xpr = feedback_propotional_gain_xpr
  phys%feedback_integral_gain_xpr = feedback_integral_gain_xpr
  phys%feedback_derivative_gain_xpr = feedback_derivative_gain_xpr
  phys%puff_slope         = puff_slope
  phys%density_source     = density_source
  phys%ener_source_e      = ener_source_e
  phys%ener_source_ee     = ener_source_ee
#ifdef KEQUATION
  phys%diff_k_min         = diff_k_min
  phys%diff_k_max         = diff_k_max
  phys%k_max              = k_max
#endif
  phys%sigma_source       = sigma_source
  phys%fluxg_trunc        = fluxg_trunc
  phys%part_source        = part_source
  phys%ener_source        = ener_source
  phys%Zeff               = Zeff
  phys%Pohmic             = Pohmic
  phys%Tbg                = Tbg
  phys%bcflags            = bcflags
  phys%bohmth             = bohmth
  phys%bohm_energy_thresh = bohm_energy_thresh
  phys%Gmbohm             = Gmbohm
  phys%Gmbohme            = Gmbohme
  phys%a                  = a
  phys%Mref               = Mref
  phys%tie                = tie
  phys%diff_pari          = diff_pari
  phys%diff_pare          = diff_pare
  phys%diff_pot           = diff_pot
  phys%c_fli              = c_fli
  phys%c_fle              = c_fle
  phys%T_fluxlim_maxi     = T_fluxlim_maxi
  phys%T_fluxlim_maxe     = T_fluxlim_maxe
  phys%epn                = epn
  phys%etapar             = etapar
  phys%Potfloat           = Potfloat
  phys%diagsource         = diagsource
  utils%PRINTint          = printint
  utils%timing            = dotiming
  utils%freqdisp          = freqdisp
  utils%freqsave          = freqsave
  lssolver%sollib         = sollib
  lssolver%timing         = lstiming
  lssolver%kspitrace      = kspitrace
  lssolver%rtol           = rtol
  lssolver%atol           = atol
  lssolver%kspitmax       = kspitmax
  lssolver%igz            = igz
  lssolver%rprecond       = rprecond
  lssolver%Nrprecond      = Nrprecond
  lssolver%kspnorm        = kspnorm
  lssolver%kspmethd       = kspmethd
  lssolver%pctype         = pctype
  lssolver%gmresres       = gmresres
  lssolver%mglevels       = mglevels
  lssolver%mgtypeform     = mgtypeform
  lssolver%itmax          = itmax
  lssolver%itrace         = itrace
  lssolver%rest           = rest
  lssolver%istop          = istop
  lssolver%tol            = tol
  lssolver%kmethd         = kmethd
  lssolver%ptype          = ptype
  lssolver%smther         = smther
  lssolver%jsweeps        = jsweeps
  lssolver%novr           = novr
  lssolver%restr          = restr
  lssolver%prol           = prol
  lssolver%solve          = solve
  lssolver%fill           = fill
  lssolver%thr            = thrsol
  lssolver%smther2        = smther2
  lssolver%jsweeps2       = jsweeps2
  lssolver%novr2          = novr2
  lssolver%restr2         = restr2
  lssolver%prol2          = prol2
  lssolver%solve2         = solve2
  lssolver%fill2          = fill2
  lssolver%thr2           = thrsol2
  lssolver%mlcycle        = mlcycle
  lssolver%outer_sweeps   = outer_sweeps
  lssolver%maxlevs        = maxlevs
  lssolver%csize          = csize
  lssolver%aggr_prol      = aggr_prol
  lssolver%par_aggr_alg   = par_aggr_alg
  lssolver%aggr_ord       = aggr_ord
  lssolver%aggr_filter    = aggr_filter
  lssolver%mncrratio      = mncrratio
  lssolver%athres         = athres
  lssolver%csolve         = csolve
  lssolver%csbsolve       = csbsolve
  lssolver%cmat           = cmat
  lssolver%cfill          = cfill
  lssolver%cthres         = cthres
  lssolver%cjswp          = cjswp

  IF (switch%steady) THEN
     msg = 'Steady state simulation'
  ELSEIF (switch%psdtime) THEN
     msg = 'Pseudotime simulation for reducing diffusion'
  ELSE
     msg = 'Time advancing simulation'
  END IF
  
  ! Some checking of the inputs
  IF (.NOT. (switch%read_gmsh .OR. switch%readMeshFromSol)) THEN
       WRITE (6, *) "Error: you must choose between reading the mesh from a Gmsh file or from a solution file"
       STOP
  END IF
  IF (time%tis .GT. 6) THEN
     WRITE (6, *) "Error: wrong time integration scheme in parameters: tis=", time%tis
     STOP
  END IF
  IF (numer%dumpnr < 0. .OR. numer%dumpnr > 1.) THEN
     WRITE (6, *) "Error: wrong dumping factor for Newton-Raphson: dumpnr=", numer%dumpnr
     STOP
  END IF
  IF (switch%testcase.GE.50 .AND. .NOT.switch%axisym) THEN
     WRITE (6, *) "Error: this should be an axisymmetric simulation"
     STOP
  ENDIF
  ! Adimensionalize units and store reference values
  CALL adimensionalization()

  ! A little message for the user...
  IF (MPIvar%glob_id .EQ. 0) THEN
     PRINT *, '                                                                      '
     PRINT *, '                                                                      '
#ifdef TOR3D
     PRINT *, '------------------------  MHDG SIMULATION IN 3D ----------------------'
#else
     PRINT *, '------------------------  MHDG SIMULATION IN 2D ----------------------'
#endif
     PRINT *, '                                                                      '
#ifndef TEMPERATURE
#ifdef NEUTRAL
#ifdef KEQUATION
     PRINT *, ' MODEL: N-Gamma isothermal with neutral with k equation               '
#else
     PRINT *, ' MODEL: N-Gamma isothermal with neutral                               '
#endif
#else
     PRINT *, ' MODEL: N-Gamma isothermal                                            '
#endif
#else
#ifdef NEUTRAL
     PRINT *, ' MODEL: N-Gamma-Ti-Te with neutral                                    '
#else
     PRINT *, ' MODEL: N-Gamma-Ti-Te                                                 '
#endif
#endif
     PRINT *, '                                                                      '
     PRINT *, '----------------------------------------------------------------------'
     PRINT *, '  Simulation type: ', ADJUSTL(TRIM(msg))
     PRINT *, '----------------------------------------------------------------------'
     PRINT *, 'Parameter file loaded:'
     PRINT *, '        ***************** Geometry ****************************'
     PRINT *, '                - R0:                                                 ', R0
     PRINT *, '                - Security factor:                                    ', q
     PRINT *, '      	***************** Magnetic ****************************'
     PRINT *, '		            - RMP amplitude:                                      ', amp_rmp
     PRINT *, '		            - RMP coil number:                                    ', nbCoils_rmp
     PRINT *, '		            - RMP coil toroidal elongation:                       ', torElongCoils_rmp
     PRINT *, '		            - RMP parity:                                         ', parite
     PRINT *, '		            - RMP row number:                                     ', nbRow
     PRINT *, '		            - Ripple amplitude:                                   ', amp_ripple
     PRINT *, '		            - Ripple coil number:                                 ', nbCoils_ripple
     PRINT *, '		            - triangularity:                                      ', triang
     PRINT *, '		            - ellipticity:                                        ', ellip
     PRINT *, '        ***************** Time stepping ************************'
     PRINT *, '                - dt0:                                                ', time%dt0
     PRINT *, '                - dt modification:                                    ', time%tsw
     PRINT *, '                - final time:                                         ', time%tfi
     PRINT *, '                - max number of time steps:                           ', time%nts
     PRINT *, '                - time integration scheme:                            ', time%tis
     PRINT *, '        ***************** Physics *****************************'
     PRINT *, '                - perp. diffusion in the continuity equation:         ', phys%diff_n
     PRINT *, '                - perp. diffusion in the momentum equation:           ', phys%diff_u
     PRINT *, '                - pinch velocity in continuity equation:              ', phys%v_p
#ifdef TEMPERATURE
     PRINT *, '                - perp. diffusion in the ions energy equation:        ', phys%diff_e
     PRINT *, '                - perp. diffusion in the electrons energy equation:   ', phys%diff_ee
     PRINT *, '                - paral. diffusion for ions temperature:              ', phys%diff_pari
     PRINT *, '                - paral. diffusion for electrons temperature:         ', phys%diff_pare
     PRINT *, '                - temperature exchange coefficient ions/electrons:    ', phys%tie
     PRINT *, '                - temperature diffusion exponential:                  ', phys%epn
     PRINT *, '                - reference Mach number:                              ', phys%Mref
     PRINT *, '                - gamma for Bohm boundary condition on ions:          ', phys%Gmbohm
     PRINT *, '                - gamma for Bohm boundary condition for electrons:    ', phys%Gmbohme
     IF(switch%ohmicsrc) THEN
       PRINT *, '             - Zeff                                                ', phys%Zeff
       PRINT *, '             - Ohmic heating                                       ', phys%Pohmic
     ENDIF
#endif
#ifdef VORTICITY
     PRINT *, '                - perp. diffusion in the vorticity equation:          ', phys%diff_vort
     PRINT *, '                - perp. diffusion in the potential equation:          ', phys%diff_pot
#endif
#ifdef NEUTRAL
     PRINT *, '                - diffusion in the neutral equation:                  ', phys%diff_nn
     PRINT *, '                - recycling coefficient in the neutral equation:      ', phys%Re
     PRINT *, '                - recycling coefficient pump in the neutral equation: ', phys%Re_pump
     PRINT *, '                - applying trim:                                      ', phys%apply_trim
     PRINT *, '                - puff coefficient in the neutral equation:           ', phys%puff
     PRINT *, '                - cryopump power coefficient in the neutral equation: ', phys%cryopump_power
     PRINT *, '                - flux limiter applied:                               ', switch%flux_limiter
     IF (switch%flux_limiter) THEN
       PRINT *, '                - c_fli (ions):                                      ', phys%c_fli
       PRINT *, '                - c_fle (electrons):                                 ', phys%c_fle
     ENDIF
     PRINT *, '                - T_maxi for the flux limiter:                       ', phys%T_fluxlim_maxi
     PRINT *, '                - T_maxe for the flux limiter:                       ', phys%T_fluxlim_maxe
     PRINT *, '                - impurity name:                                     ', TRIM(ADJUSTL(phys%impurity_name))
     PRINT *, '                - impurity concentration:                            ', phys%impurity_concentration
     PRINT *, '                - import 1D diffusion profile:                       ', switch%import_diffusion_1D
     PRINT *, '                - Bohm-gyro-Bohm model applied:                      ', switch%bohm_gyrobohm
     IF (switch%ME) THEN
        PRINT *, '                - I_0 for moving equilibrium:                         ', phys%I_0
        PRINT *, '             - puff increment slope:                               ', phys%puff_slope
          IF (switch%target_variable /= 0) THEN
             PRINT *, '             - feedback_propotional_gain:                               ', phys%feedback_propotional_gain
             PRINT *, '             - feedback_integral_gain:                                  ', phys%feedback_integral_gain
             PRINT *, '             - feedback_derivative_gain:                                 ', phys%feedback_derivative_gain
             IF (switch%target_variable == 3) THEN
                PRINT *, '             - feedback_propotional_gain_xpr:                            ', phys%feedback_propotional_gain_xpr
                PRINT *, '             - feedback_integral_gain_xpr:                               ', phys%feedback_integral_gain_xpr
                PRINT *, '             - feedback_derivative_gain_xpr:                              ', phys%feedback_derivative_gain_xpr
             ENDIF
          ENDIF
     ENDIF
     PRINT *, '                - particle source at core:                            ', part_source
     PRINT *, '                - energy source at core:                              ', ener_source
#endif
#ifdef KEQUATION
     PRINT *, '                - minimum perp diffusion in the k equation:           ', phys%diff_k_min
     PRINT *, '                - maximum perp diffusion in the k equation:           ', phys%diff_k_max
     PRINT *, '                - maximum k:                                          ', phys%k_max
#endif
     PRINT *, '                - constant for the momentum equation (isoth)          ', phys%a
     PRINT *, '                - diagonal implicit sources                           ', phys%diagsource
     PRINT *, '        ***************** Switches ****************************'
     PRINT *, '                - stady state simulation:                             ', switch%steady
     PRINT *, '                - axisym:                                             ', switch%axisym
     IF (switch%init.EQ.1) THEN
        PRINT *, '  Initializing with analytical solution at nodes                    '
     ELSE
        PRINT *, '  Initializing with L2 projection                                   '
     ENDIF
     PRINT *, '                - driftdia:                                           ', driftdia
     PRINT *, '                - driftexb:                                           ', driftexb
     PRINT *, '                - test case:                                          ', testcase
     PRINT *, '                - Ohmic heating:                                      ', OhmicSrc
     PRINT *, '                - Moving equilibrium:                                 ', ME
     PRINT *, '                - Diffusion adjusted to plasma current:               ', diff_reverse_Ip
     PRINT *, '                - RMP:                                                ', RMP
     PRINT *, '                - Ripple:                                             ', Ripple
     PRINT *, '                - shockcp:                                            ', shockcp
     PRINT *, '                - minrho:                                             ', minrho
     PRINT *, '                - logrho:                                             ', logrho
     PRINT *, '                - thresh:                                             ', thresh
     PRINT *, '                - filter:                                             ', filter
#ifdef TEMPERATURE
     PRINT *, '                - decoup:                                             ', decoup
#endif
     PRINT *, '                - ckeramp:                                            ', ckeramp
     PRINT *, '                - saveNR:                                             ', saveNR
     PRINT *, '                - saveTau:                                            ', saveTau
     PRINT *, '                - transport_1d:                          ', transport_1d
     PRINT *, '                - transport_model_path:                  ', TRIM(ADJUSTL(input%transport_model_path))
     PRINT *, '                - rho_core (transport model):            ', transport_model_input%rho_core
     PRINT *, '                - rho_edge (transport model):            ', transport_model_input%rho_edge
     PRINT *, '                - rho_diffusion_model_max (transport model):       ', transport_model_input%rho_diffusion_model_max
     PRINT *, '                - c_bohm_i (transport model):           ', transport_model_input%c_bohm_i
     PRINT *, '                - c_gyrobohm_i (transport model):       ', transport_model_input%c_gyrobohm_i
     PRINT *, '                - c_bohm_e (transport model):           ', transport_model_input%c_bohm_e
     PRINT *, '                - c_gyrobohm_e (transport model):       ', transport_model_input%c_gyrobohm_e
     PRINT *, '                - c_bohm_n (transport model):           ', transport_model_input%c_bohm_n
     PRINT *, '                - prandtl (transport model):            ', transport_model_input%prandtl
     PRINT *, '                - pinch_model (transport model):        ', transport_model_input%pinch_model
     PRINT *, '                - c_pinch (transport model):            ', transport_model_input%c_pinch
     PRINT *, '                - nu_th (transport model):              ', transport_model_input%nu_th
     PRINT *, '        ***************** Numerics ****************************'
     PRINT *, '                - stabilization type:                                 ', numer%stab
     PRINT *, '                - tau(1):                                             ', numer%tau(1)
     PRINT *, '                - tau(2):                                             ', numer%tau(2)
     PRINT *, '                - tau(3):                                             ', numer%tau(3)
     PRINT *, '                - tau(4):                                             ', numer%tau(4)
     PRINT *, '                - tau(5):                                             ', numer%tau(5)
     PRINT *, '                - max number of N-R iterations:                       ', numer%nrp
     PRINT *, '                - tolerance for the N-R scheme:                       ', numer%tNR
     PRINT *, '                - tolerance for the steady state achievement:         ', numer%tTM
     IF (switch%shockcp .GT. 0) THEN
        PRINT *, '                - shock capturing coeff:                            ', numer%sc_coe
        PRINT *, '                - shock capturing sensibility:                      ', numer%sc_sen
     END IF
     IF (switch%limrho .GT. 0) THEN
        PRINT *, '                - applying limiting of rho at value:                ', numer%minrho
     END IF
     IF (switch%limrho .EQ. 1 .OR. switch%limrho .EQ. 3) THEN
        PRINT *, '                - coefficient of the source for limiting rho:       ', numer%so_coe
     END IF
     IF (switch%limrho .EQ. 2 .OR. switch%limrho .EQ. 3) THEN
        PRINT *, '                - coefficient of the diffusion for limiting rho:    ', numer%df_coe
     END IF
     IF (switch%difcor .GT. 0) THEN
        PRINT *, '                - adding diffusion in corners, in position:         ', switch%difcor
        PRINT *, '                - diffusion coefficient in corners:                 ', numer%dc_coe
     END IF
     IF (switch%thresh .GT. 0) THEN
        PRINT *, '                - using a threshold at rho:                         ', numer%thr
        PRINT *, '                - using a threshold at pressure:                    ', numer%thrpre
     END IF
     PRINT *, '                - using a dumping factor for Newton-Raphson:          ', numer%dumpnr
#ifdef TOR3D
     PRINT *, '                - number of elements in the toroidal direction:       ', numer%ntor
     PRINT *, '                - polynomial degree in the toroidal direction:        ', numer%ptor
#ifdef PARALL
     PRINT *, '                - number of MPI partitions in the toroidal direction :', numer%npartor
#endif
     PRINT *, '                - max extention in the toroidal direction:            ', numer%tmax
#endif
     PRINT *, '                - dumping for ExB term:       ', numer%exbdump
     PRINT *, '                - type of bohm boundary condition:                     ',  numer%bohmtypebc
     IF(adapt%adaptivity) THEN
        PRINT *, '        ***************** Adaptivity ****************************'
        PRINT *, '                - adaptivity :                                         ', adapt%adaptivity
        PRINT *, '                - adaptivity evaluator:                                ', adapt%evaluator
        PRINT *, '                - option for shock capturing:                          ', adapt%shockcp_adapt
        PRINT *, '                - threshold for the indicator:                         ', adapt%thr_ind
        PRINT *, '                - tollerance estimator:                                ', adapt%tol_est
        PRINT *, '                - parameter for estimator:                             ', adapt%param_est
        PRINT *, '                - adaptivity in time:                                  ', adapt%time_adapt
        PRINT *, '                - adaptivity in NR:                                    ', adapt%NR_adapt
        PRINT *, '                - frequence of time adaptivity:                        ', adapt%freq_t_adapt
        PRINT *, '                - frequence of NR adaptivity:                          ', adapt%freq_NR_adapt
        PRINT *, '                - NR divergence adaptivity:                            ', adapt%div_adapt
        PRINT *, '                - restart adaptivity:                                  ', adapt%rest_adapt
        PRINT *, '                - oscillations on the solution adaptivity:             ', adapt%osc_adapt
        PRINT *, '                - oscillations tollerance adaptivity:                  ', adapt%osc_tol
        PRINT *, '                - threshold oscillations for saving checkpoint:        ', adapt%osc_check

        IF(adapt%difference .EQ. 0) THEN
           PRINT *, '                - difference type for the estimator:              ', 'relative'
        ELSE
           PRINT *, '                - difference type for the estimator:              ', 'absolute'
        ENDIF
     ENDIF
     PRINT *, '        ***************** Linear solver params******************'
     IF (lssolver%sollib == 1) THEN
        PRINT *, '                - Library used for the linear system:   PASTIX      '
     ELSEIF (lssolver%sollib == 2) THEN
        PRINT *, '                - Library used for the linear system:   PSBLAS      '
        PRINT *, '                - Iterative method:                                 ', lssolver%kmethd
        PRINT *, '                - Preconditioner:                                   ', lssolver%ptype
        PRINT *, '                - Stopping criterion type                           ', lssolver%istop
        PRINT *, '                - Stopping criterion tolerance                      ', lssolver%tol
        PRINT *, '                - Restart:                                          ', lssolver%rest
     ELSEIF (lssolver%sollib == 3) THEN
        PRINT *, '                - Library used for the linear system:   PETSc      '
        PRINT *, '                - Iterative method:                                 ', lssolver%kspmethd
        PRINT *, '                - Preconditioner:                                   ', lssolver%pctype
        PRINT *, '                - Relative tollerance                               ', lssolver%rtol
        PRINT *, '                - Absolute tollerance                               ', lssolver%atol
        PRINT *, '                - Max number of iteration:                          ', lssolver%kspitmax
        PRINT *, '                - Set to zero the initial guess                     ', lssolver%igz
        PRINT *, '                - Recompute preconditioner at each NR               ', lssolver%rprecond
        PRINT *, '                - Norm type                                         ', lssolver%kspnorm
        PRINT *, '                - Gmres restart value                               ', lssolver%gmresres
        PRINT *, '                - Display convergence at each iteration             ', lssolver%kspitrace
        PRINT *, '                - MultiGrid (MG) levels                             ', lssolver%mglevels
        PRINT *, '                - MultiGrid (MG) type form                          ', lssolver%mgtypeform
     ENDIF

     PRINT *, '        '
  END IF

END SUBROUTINE READ_input

SUBROUTINE read_transport_model_input()
  USE globals
  USE MPI_OMP
  IMPLICIT NONE

  INTEGER :: utransport, ios
  REAL*8 :: rho_core, rho_edge, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys
  INTEGER :: pinch_model
  NAMELIST /TRANSPORT_MODEL_1D_LST/ rho_core, rho_edge, rho_diffusion_model_max, c_bohm_i, c_gyrobohm_i, c_bohm_e, c_gyrobohm_e, c_bohm_n, prandtl, pinch_model, c_pinch, nu_th, vpinch_const_phys, rho_pinch_axis_width, rho_pinch_model_max, rho_pinch_edge_width, rho_blend_width, diff_n_min_phys, diff_u_min_phys, diff_e_min_phys, diff_ee_min_phys

  rho_core = transport_model_input%rho_core
  rho_edge = transport_model_input%rho_edge
  rho_diffusion_model_max = transport_model_input%rho_diffusion_model_max
  c_bohm_i = transport_model_input%c_bohm_i
  c_gyrobohm_i = transport_model_input%c_gyrobohm_i
  c_bohm_e = transport_model_input%c_bohm_e
  c_gyrobohm_e = transport_model_input%c_gyrobohm_e
  c_bohm_n = transport_model_input%c_bohm_n
  prandtl = transport_model_input%prandtl
  pinch_model = transport_model_input%pinch_model
  c_pinch = transport_model_input%c_pinch
  nu_th = transport_model_input%nu_th
  vpinch_const_phys = transport_model_input%vpinch_const_phys
  rho_pinch_axis_width = transport_model_input%rho_pinch_axis_width
  rho_pinch_model_max = transport_model_input%rho_pinch_model_max
  rho_pinch_edge_width = transport_model_input%rho_pinch_edge_width
  rho_blend_width = transport_model_input%rho_blend_width
  diff_n_min_phys = transport_model_input%diff_n_min_phys
  diff_u_min_phys = transport_model_input%diff_u_min_phys
  diff_e_min_phys = transport_model_input%diff_e_min_phys
  diff_ee_min_phys = transport_model_input%diff_ee_min_phys

  IF (LEN_TRIM(input%transport_model_path) == 0) RETURN

  utransport = 101
  OPEN(utransport, file=TRIM(ADJUSTL(input%transport_model_path)), status='old', iostat=ios)
  IF (ios /= 0) THEN
     IF (MPIvar%glob_id == 0) WRITE(6,*) 'Warning: could not open transport model settings file: ', TRIM(ADJUSTL(input%transport_model_path))
     RETURN
  END IF

  READ(utransport, nml=TRANSPORT_MODEL_1D_LST, iostat=ios)
  CLOSE(utransport)
  IF (ios /= 0) THEN
     IF (MPIvar%glob_id == 0) WRITE(6,*) 'Warning: could not read TRANSPORT_MODEL_1D_LST from file: ', TRIM(ADJUSTL(input%transport_model_path))
     RETURN
  END IF

  transport_model_input%rho_core = rho_core
  transport_model_input%rho_edge = rho_edge
  transport_model_input%rho_diffusion_model_max = rho_diffusion_model_max
  transport_model_input%c_bohm_i = c_bohm_i
  transport_model_input%c_gyrobohm_i = c_gyrobohm_i
  transport_model_input%c_bohm_e = c_bohm_e
  transport_model_input%c_gyrobohm_e = c_gyrobohm_e
  transport_model_input%c_bohm_n = c_bohm_n
  transport_model_input%prandtl = prandtl
  transport_model_input%pinch_model = pinch_model
  transport_model_input%c_pinch = c_pinch
  transport_model_input%nu_th = nu_th
  transport_model_input%vpinch_const_phys = vpinch_const_phys
  transport_model_input%rho_pinch_axis_width = rho_pinch_axis_width
  transport_model_input%rho_pinch_model_max = rho_pinch_model_max
  transport_model_input%rho_pinch_edge_width = rho_pinch_edge_width
  transport_model_input%rho_blend_width = rho_blend_width
  transport_model_input%diff_n_min_phys = diff_n_min_phys
  transport_model_input%diff_u_min_phys = diff_u_min_phys
  transport_model_input%diff_e_min_phys = diff_e_min_phys
  transport_model_input%diff_ee_min_phys = diff_ee_min_phys
END SUBROUTINE read_transport_model_input
