!************************************************************
! project: MHDG
! file: types.f90
! date: 06/09/2016
! Declaration of all the types used
! in the code
!************************************************************

MODULE types
  USE prec_const
  USE mod_splines
  IMPLICIT NONE

  !*******************************************************
  ! Boundary condition flags
  !*******************************************************
  INTEGER, PARAMETER, PUBLIC :: max_num_diff_bc = 100
  ! 1-10 Dirichlet type
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet = 1
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet_weak_form = 2
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet_weak_form_oldvalues = 3
  INTEGER, PARAMETER, PUBLIC :: bc_Transmission = 4
  INTEGER, PARAMETER, PUBLIC :: bc_dirichlet_and_Neumann = 5
   INTEGER, PARAMETER, PUBLIC :: bc_iter_core = 6
  ! 20-.. In out type
  INTEGER, PARAMETER, PUBLIC :: bc_inout = 20
  ! 30-.. Neumann type
  INTEGER, PARAMETER, PUBLIC :: bc_NeumannH = 30
  ! 50-.. Bohm type
  INTEGER, PARAMETER, PUBLIC :: bc_Bohm = 50
  INTEGER, PARAMETER, PUBLIC :: bc_BohmNM = 51
  INTEGER, PARAMETER, PUBLIC :: bc_BohmSimp = 52
  INTEGER, PARAMETER, PUBLIC :: bc_BohmPump = 55 ! bohm for plasma, pump for neutrals
  INTEGER, PARAMETER, PUBLIC :: bc_BohmPuff = 56 ! bohm for plasma, puff for neutrals

  ! 60-.. Slip wall type
  INTEGER, PARAMETER, PUBLIC :: bc_slip_wall = 60
  ! 70-.. Periodic
  INTEGER, PARAMETER, PUBLIC :: bc_periodic = 70


  ! Boundary conditions flags and names
  CHARACTER(100) ::  bc_flag_type(10), bc_flag_name(max_num_diff_bc)

  !*******************************************************
  ! Reference element
  !*******************************************************
  TYPE :: Reference_element_type
     INTEGER*4 :: ndim          ! Number of dimensions of the reference element
     INTEGER*4 :: elemType      ! type of mesh elements (0 = triangles, 1 = quadrilaterals, 2 = thetrahedra, 3 = hexahedra)
     INTEGER*4 :: Nvertices     ! Number of vertices of the element: 3 for triangles, 4 for quads and thetraedra, 8 for hexahedra
     INTEGER*4 :: Ndeg          ! degree of interpolation
     INTEGER*4 :: Nnodes3D      ! number of 3D nodes in the element
     INTEGER*4 :: Nnodes2D      ! number of 2D nodes in the element
     INTEGER*4 :: Nnodes1D      ! number of 1D nodes in the element
     INTEGER*4 :: Nfaces        ! number of faces in the element: a face in a 2D mesh is a segment, in a 3D mesh is a 2D entity
     INTEGER*4 :: Nbords        ! number of bords (only in 3D), is the number of 1D entities of the element
     INTEGER*4 :: Nextnodes     ! number of exterior nodes of the element
     INTEGER*4 :: Nfacenodes    ! number of nodes in each face of the element
     INTEGER*4 :: Nfacenodeslin ! number of nodes in each linear face of the element
     INTEGER*4 :: Nbordnodes    ! number of nodes in each bord of the element (only in 3D)
     INTEGER*4 :: Ninnernodes   ! number of inner nodes of the element
     INTEGER*4 :: Ninnernodesface  ! number of inner nodes of the 2D faces (only for 3D meshes)
     INTEGER*4 :: NGauss1D      ! number of Gauss points of the 1D quadrature
     INTEGER*4 :: NGauss2D      ! number of Gauss points of the 2D quadrature
     INTEGER*4 :: NGauss3D      ! number of Gauss points of the 3D quadrature
     INTEGER*4, ALLOCATABLE :: Face_nodes(:, :)  ! numbering of the face nodes
     INTEGER*4, ALLOCATABLE :: Bord_nodes(:)    ! numbering of the 1D edge nodes (for 3D meshes)
     INTEGER*4, ALLOCATABLE :: inner_nodes(:)   ! numbering of the inner nodes
     INTEGER*4, ALLOCATABLE :: inner_nodes_face(:)   ! numbering of the inner nodes for 2D faces (for 3D meshes)
     REAL*8, POINTER :: coord3d(:, :) => NULL()       ! spatial coordinates of the 3D nodes
     REAL*8, POINTER :: coord2d(:, :) => NULL()       ! spatial coordinates of the 2D nodes
     REAL*8, POINTER :: coord1d(:) => NULL()       ! spatial coordinates of the 1D nodes
     REAL*8, ALLOCATABLE :: gauss_points3D(:, :) ! Gauss integration points for the 3D quadrature
     REAL*8, ALLOCATABLE :: gauss_points2D(:, :) ! Gauss integration points for the 2D quadrature
     REAL*8, ALLOCATABLE :: gauss_points1D(:)   ! Gauss integration points for the 1D quadrature
     REAL*8, ALLOCATABLE :: gauss_weights3D(:)  ! Weights of the integration points for the 3D quadrature
     REAL*8, ALLOCATABLE :: gauss_weights2D(:)  ! Weights of the integration points for the 2D quadrature
     REAL*8, ALLOCATABLE :: gauss_weights1D(:)  ! Weights of the integration points for the 1D quadrature
     REAL*8, ALLOCATABLE :: N3D(:, :)      ! Shape functions 3D at the Gauss points
     REAL*8, ALLOCATABLE :: Nxi3D(:, :)    ! Derivative with respect to xi of the shape functions 3D at the Gauss points
     REAL*8, ALLOCATABLE :: Neta3D(:, :)   ! Derivative with respect to eta of the shape functions 3D at the Gauss points
     REAL*8, ALLOCATABLE :: Nzeta3D(:, :)  ! Derivative with respect to zeta of the shape functions 3D at the Gauss points
     REAL*8, ALLOCATABLE :: N2D(:, :)      ! Shape functions 2D at the Gauss points
     REAL*8, ALLOCATABLE :: Nxi2D(:, :)    ! Derivative with respect to xi of the shape functions 2D at the Gauss points
     REAL*8, ALLOCATABLE :: Neta2D(:, :)   ! Derivative with respect to eta of the shape functions 2D at the Gauss points
     REAL*8, ALLOCATABLE :: N1D(:, :)      ! Shape functions 1D at the Gauss points
     REAL*8, ALLOCATABLE :: Nxi1D(:, :)    ! Derivative of the shape functions 1D at the Gauss points
     REAL*8, ALLOCATABLE :: Nlin(:, :)     ! Linear shape functions 2D at the nodes of the element (only used in shock capturing so far...)
     REAL*8, ALLOCATABLE :: sFTF(:, :)     ! Shape functions for the toroidal faces
     INTEGER*4, ALLOCATABLE :: faceNodes3(:, :) ! Face numbering for toroidal faces
     INTEGER*4          :: Nfl           ! Number of nodes in lateral faces for toroidal 3d computations
     INTEGER*4          :: Ngl           ! Number of Gauss points in lateral faces for toroidal 3d computations
     INTEGER*4          :: Nft           ! Number of nodes in all faces for toroidal 3d computations
  END TYPE Reference_element_type

  !*******************************************************
  ! Mesh
  !*******************************************************
  TYPE :: Mesh_type
     INTEGER*4 :: Ndim                ! Number of dimensions of the mesh
     INTEGER*4 :: Nnodes              ! Number of nodes in the mesh
     INTEGER*4 :: Nnodesperelem       ! Number of nodes per element in the mesh (dimension 2 of the connectvity matrix)
     INTEGER*4 :: Nnodesperface       ! Number of nodes per face in the mesh (dimension 2 of boundary connectvity matrix )
     INTEGER*4 :: Nelems              ! Number of elements in the mesh
     INTEGER*4 :: Nfaces              ! Number of faces in the mesh
     INTEGER*4 :: Nextfaces           ! Number of exterior faces
     INTEGER*4 :: Nintfaces           ! Number of interior faces
     INTEGER*4 :: elemType            ! 0 for quads - 1 for triangles - 2 for thetra - 3 for hexa
     INTEGER*4 :: ndir                ! number of Dirichlet faces
     INTEGER*4 :: ukf                 ! number of faces in the mesh minus the number of Dirichlet faces
     INTEGER*4, POINTER :: T(:, :) => NULL()    ! Elements connectivity matrix
     INTEGER*4, POINTER :: Tlin(:, :) => NULL()    ! Elements linear connectivity matrix
     INTEGER*4, POINTER :: Tb(:, :) => NULL()    ! Outer faces connectivity matrix
     INTEGER*4, POINTER :: boundaryFlag(:) => NULL()    ! Flag for the boundary condition for each external face (set in the mesh generator)
     INTEGER*4, ALLOCATABLE :: F(:, :)          ! Faces connectivity matrix
     INTEGER*4, ALLOCATABLE :: N(:, :)          ! Nodes connectivity matrix
     INTEGER*4, ALLOCATABLE :: face_info(:, :)          ! Elemental face info
     INTEGER*4, ALLOCATABLE :: faces(:, :, :)    ! for each triangle i, stores info k on each face j: faces(i,j,1) = # of neighbouring triangle (0 if external
    ! boundary), faces(i,j,2) = type of boundary (), faces(i,j,3) = type of boundary condition
     INTEGER*4, ALLOCATABLE :: extfaces(:, :)   ! for each exterior face, stores the number of the triangle, the number of the face, and the type of BC
     INTEGER*4, ALLOCATABLE :: intfaces(:, :)   ! for each interior face, stores the number of the triangle, the number of the face, the number of the
    ! neighboring triangle, the number of its face and the number of the node of the neighboring triangle that
    ! matches the first knot of the triangle
     LOGICAL, ALLOCATABLE :: flipface(:, :)    ! for each triangle, and for each face, 0 if the order of the numbering in the face is to be kept, 1 if the
    ! order is to be reversed
     LOGICAL, ALLOCATABLE    :: Fdir(:, :)         ! for each element states if each local face is of Dirichlet type
     INTEGER*4, ALLOCATABLE  :: periodic_faces(:)  ! Mapping for periodic faces
     INTEGER*4, ALLOCATABLE  :: Diric(:)
     INTEGER*4, ALLOCATABLE  :: numberbcs(:)
     REAL*8, ALLOCATABLE     :: elemSize(:)   ! element size (area in 2D, volume in 3D) [Number of elements]
     REAL*8, POINTER         :: X(:, :) => NULL()    ! nodes coordinates
     INTEGER*4              :: Nnodes_toroidal     ! Number of nodes in the toroidal direction
     REAL*8, POINTER         :: toroidal(:) => NULL()    ! nodes coordinates in the toroidal direction
    ! Limiting & shock capturing stuff
     INTEGER, ALLOCATABLE    :: flag_elems_rho(:)      ! Flagged elements for limiting rho [Number of elements]
     INTEGER, ALLOCATABLE    :: flag_elems_sc(:)     !  Flagged elements for shock-capturing [Number of elements]
     REAL*8, ALLOCATABLE     :: minrho_elems(:)        ! Minimum value of density in the flagged elements[Number of elements]
     REAL*8, ALLOCATABLE     :: sour_elems(:)          ! Source to limit rho in the flagged elements[Number of elements]
     REAL*8, ALLOCATABLE     :: diff_elems(:)          ! Diffusion to limit rho in the flagged elements[Number of elements]
     REAL*8, ALLOCATABLE     :: scdiff_nodes(:, :)      ! Shock capturing diffusion in each node [Number of elements,Number of nodes per element]
     REAL*8                 :: xmax, xmin, ymax, ymin    ! Limit of the GLOBAL matrix, across mpi partitions
     REAL*8                  :: puff_area         ! area of the puff bounday condition
     REAL*8                  :: pump_area         ! area of the pump bounday condition
     REAL*8                  :: core_area         ! area of the core bounday condition
     REAL*8,ALLOCATABLE      :: Xg(:,:)               ! 2D Gauss point coordinates
     REAL*8,ALLOCATABLE      :: Xgf(:,:)              ! 1D Gauss point coordinates at interior faces
     REAL*8,ALLOCATABLE      :: Xgb(:,:)          ! 1D Gauss point  coordinates at boundary faces
     INTEGER*4, POINTER	    :: T_gmsh(:, :) => NULL()    ! Elements connectivity matrix to write the msh file
     INTEGER*4, POINTER      :: Tb_gmsh(:, :) => NULL()    ! Outer faces connectivity matrix to write the msh file
     REAL*8, POINTER         :: X_P1(:,:) => NULL()         ! coordinates of the nodes on the P1 mesh to write to the msh file
#ifdef PARALL
     INTEGER, POINTER       :: loc2glob_fa(:)     ! mapping number of the faces for creating the global matrix [number of faces in the mesh]
     INTEGER, POINTER       :: loc2glob_el(:)     ! mapping number of the elements from local to global [number of elements in the mesh]
     INTEGER, POINTER        :: loc2glob_nodes(:)    ! mapping number of the nodes from local to global [number of nodes in the mesh]
     INTEGER, POINTER       :: ghostfaces(:)      ! integer that states if a face is to be assembled locally or not [number of faces in the mesh]
     INTEGER, POINTER       :: ghostelems(:)      ! integer that states if an element is to be assembled locally or not (for 3D) [number of elements in the mesh]
     INTEGER               :: nghostfaces        ! number of ghost faces
     INTEGER               :: nghostelems        ! number of ghost elements (used only for 3D)
     INTEGER               :: Nel_glob           ! number of elements of global mesh
     INTEGER               :: Nfa_glob           ! number of faces of global mesh
     INTEGER               :: Nno_glob           ! number of nodes of global mesh
     INTEGER               :: Ndir_glob          ! number of Dirichlet faces in the global mesh
     INTEGER               :: Ngho_glob          ! number of ghost faces in the global mesh
    ! readed from input
     INTEGER, POINTER       :: ghostflp(:)        ! flipFaces for the ghost faces [number of ghost faces]
     INTEGER, POINTER       :: ghostloc(:)        ! local numbering of the ghost face in the process that assemble it [number of ghost faces]
     INTEGER, POINTER       :: ghostpro(:)        ! the process that assemble the ghost face [number of ghost faces]
     INTEGER, POINTER       :: ghelsloc(:)        ! local numbering of the ghost element in the process that assemble it (only 3D) [number of ghost elements]
     INTEGER, POINTER       :: ghelspro(:)        ! the process that assemble the ghost element (only 3D) [number of ghost elements]
    ! built after reading from input
     INTEGER, ALLOCATABLE   :: fc2sd(:)          ! face 2 send: faces computed locally that the local process has to send (vector)
     INTEGER, ALLOCATABLE   :: pr2sd(:)          ! process 2 send: to which process the faces computed locally need to be sent (vector)
     INTEGER, ALLOCATABLE   :: fc2rv(:)          ! face 2 receive: ghost faces computed by other processes that the local process need to receive[number of ghost faces]
     INTEGER, ALLOCATABLE   :: pr2rv(:)          ! process 2 receive: from which process the faces computed externally need to be received [number of ghost faces] (it is the same as ghostpro)

     INTEGER, ALLOCATABLE   :: el2sd(:)          ! element 2 send: elements computed locally that the local process has to send (only 3D) (vector)
     INTEGER, ALLOCATABLE   :: pe2sd(:)          ! process 2 send: to which process the elements computed locally need to be sent (only 3D) (vector)
     INTEGER, ALLOCATABLE   :: el2rv(:)          ! element 2 receive: ghost elements computed by other processes that the local process need to receive (only 3D) [number of ghost elements]
     INTEGER, ALLOCATABLE   :: pe2rv(:)          ! process 2 receive: from which process the elements computed externally need to be received (only 3D) [number of ghost elements] (it is the same as ghelspro)

    !     integer,allocatable   :: connpro(:)        ! processes connected to the local process
#endif
  END TYPE Mesh_type

  TYPE Splines_DT
     INTEGER                       :: boundaryFlag
     INTEGER                       :: spline_number
     INTEGER                       :: tot_n_splines
     INTEGER                       :: n_points
     REAL*8                        :: angle, xmax, xmin, ymax, ymin, xcenter, ycenter
     REAL*8                        :: RotMat(2,2), AntiRotMat(2,2)
     INTEGER, ALLOCATABLE          :: points_number(:)
     REAL*8, ALLOCATABLE           :: points_coord(:,:)
     TYPE(spline)                  :: spline

  ENDTYPE Splines_DT

  !*******************************************************
  ! Physics: type for physical model info
  !*******************************************************
  TYPE Physics_type
     INTEGER         :: neq            ! Number of equations
     INTEGER         :: npv            ! Number of physical variables
     REAL*8          :: diff_n, diff_u ! Perpendicular diffusion in the continuity and momentum equation
     REAL*8          :: v_p            ! Pinch velocity in the continuity equation
     REAL*8          :: a              ! Proportionality constant between pressure and density for isothermal model (p = a*rho)
     REAL*8          :: dfcoef         ! Constant related to the diamagnetic drift velocity
     REAL*8          :: dexbcoef         ! Constant related to the ExB drift velocity
     INTEGER*4       :: bcflags(1:10)          ! Set the correspondence between mesh boundary flag (Mesh%boundaryFlag) and the boundary condition
     REAL*8          :: diagsource(1:10)       ! Diagonal implicit sources
     CHARACTER(LEN=20), POINTER:: phyVarNam(:) => NULL() ! Names of the physical variables (set in initPhys)
     CHARACTER(LEN=20), POINTER:: conVarNam(:) => NULL() ! Names of the conservative variables (set in initPhys)
     REAL*8          :: lscale          ! Length scale for the non-dimensionalization of the equations
    ! Magnetic field defined for each node of the mesh.
     REAL*8, POINTER :: B(:, :)            ! Magnetic field, Br,Bz,Bphi [n of nodes  x 3]
     REAL*8          :: B0                 ! Reference value for the magnetic field [Tesla]
     REAL*8, POINTER :: magnetic_flux(:)   ! Magnetic flux   [n of nodes]
     REAL*8, POINTER :: magnetic_psi(:)    ! Magnetic flux normalized to separatrix magnetic flux   [n of nodes]
     REAL*8          :: Flux2Dmin          ! Minimum of the magnetic flux, across the MPI partitions
     REAL*8          :: Flux2Dmax          ! Maximum of the magnetic flux, across the MPI partitions
#ifdef KEQUATION
     REAL*8, POINTER :: omega(:)           ! larmor frequency   [n of nodes]
     REAL*8, POINTER :: q_cyl(:)           ! q cylindrical   [n of nodes]
#endif
     REAL*8          :: r_axis             ! R-coordinate of magnetic axis
     REAL*8          :: z_axis             ! Z-coordinate of magnetic axis

     REAL*8, POINTER :: Bperturb(:, :)     ! Magnetic perturbation, Br,Bz,Bphi [n of nodes  x 3]
     REAL*8          :: Tbg                ! Background temperature in the isothermal model
     REAL*8, POINTER :: Jtor(:)            ! Toroidal Current
     REAL*8          :: I_p                ! Total plasma current
     REAL*8          :: bohmth             ! Threshold for imposing the Bohm boundary condition
    ! Energy equation coefficients
     REAL*8          :: diff_e             ! Perpendicular diffusion in the energy equation
     REAL*8          :: epn                ! Exponential of the parallel diffusion (usually 5/2)
     REAL*8          :: Mref               ! Reference Mach number
     REAL*8          :: diff_pari          ! Parallel diffusion for the temperature (usually 1e7)
     REAL*8          :: Gmbohm             ! gamma for Bohm boundary condition on energy:
    ! Temperature equations coefficients (the ions coefficients are the ones defined previously)
     REAL*8          :: diff_ee            ! Perpendicular diffusion in the elcetron energy equation
     REAL*8          :: diff_pare          ! Parallel diffusion for the electron temperature
     REAL*8          :: tie                ! Temperature exchange coefficient between ions and electrons
     REAL*8          :: Gmbohme            ! gamma for Bohm boundary condition on electron energy:
     REAL*8          :: Pohmic             ! Ohmic heating power
    ! Coefficients for the vorticity equations
     REAL*8          :: diff_vort
     REAL*8          :: diff_pot
     REAL*8          :: etapar
     REAL*8          :: c1, c2             ! coefficients coming from the adimensionalization
     REAL*8          :: Potfloat
    ! Coefficients for ion heating
    real*8          :: heating_power      ! total power of ion heating source
    real*8          :: heating_amplitude  ! adimentionalized gaussian amplitude including integration coefficients
    real*8          :: heating_dr         ! displacement of the soruce from magnetic axis in r direction
    real*8          :: heating_dz         ! displacement of the soruce from magnetic axis in z direction
    real*8          :: heating_sigmar     ! width of the soruce from magnetic axis in r direction
    real*8          :: heating_sigmaz     ! width of the soruce from magnetic axis in z direction
    integer         :: heating_equation   ! Equation to which additional heating is applied (3 for ions, 4 for  electrons)
    ! Coefficients for the neutral equations
     REAL*8          :: diff_nn            ! Diffusion in the neutral equation
     REAL*8,ALLOCATABLE:: diff_nn_Vol(:)   ! Diffusion in the neutral equation at 2D Gauss points
     REAL*8,ALLOCATABLE:: diff_nn_Fac(:)   ! Diffusion in the neutral equation at 1D Gauss points on interior faces
     REAL*8,ALLOCATABLE:: diff_nn_Bou(:)   ! Diffusion in the neutral equation at 1D Gauss points on boundary faces
     REAL*8,ALLOCATABLE:: v_nn_Vol(:,:)    ! Convective velocity in the neutral equation at 2D Gauss points
     REAL*8,ALLOCATABLE:: v_nn_Fac(:,:)    ! Convective velocity in the neutral equation at 1D Gauss points on interior faces
     REAL*8,ALLOCATABLE:: v_nn_Bou(:,:)    ! Convective velocity in the neutral equation at 1D Gauss points on boundary faces
     REAL*8          :: Re                 ! Recycling for the neutral equation
     REAL*8          :: Re_pump            ! Recycling for the neutral equation in the pump region
     REAL*8          :: puff               ! Puff coefficient
    real*8          :: cryopump_power     ! Cryopump power in [m^3/s] coefficient
     REAL*8          :: puff_slope         ! Puff increment coefficient (only for moving equilibrium for ITER)
     REAL*8,POINTER :: puff_exp(:)    ! Puff experimental coefficient (only for moving equilibriums)
     REAL*8          :: part_source        ! Particle source for ITER
     REAL*8          :: ener_source        ! Particle source for ITER
     REAL*8          :: density_source     ! Density source for WEST (2D, case 52)
     REAL*8          :: ener_source_e      ! Ion energy source for WEST and ITER (2D, case 52 and 81)
     REAL*8          :: ener_source_ee     ! Electron source for WEST and ITER (2D, case 52 and 81)
     REAL*8          :: sigma_source       ! Sigma for the gaussian sources for WEST and ITER (2D, case 52 and 81)
     REAL*8          :: fluxg_trunc        ! Value of the NORMALISED magnetic flux at which to truncate the gaussian sources for WEST (2D, case 52), refer to Source_shape.m file
    ! Diffusion coefficients ITER evolving equilibrium
     REAL*8          :: ME_diff_n
     REAL*8          :: ME_diff_u
     REAL*8          :: ME_diff_e
     REAL*8          :: ME_diff_ee
#ifdef EXPANDEDCX
#ifdef AMJUELCX
     REAL*8, DIMENSION(9):: alpha_cx       ! Coefficients for charge exchange coefficients spline
#endif
#ifdef THERMALCX
     REAL*8, DIMENSION(5):: alpha_cx       ! Coefficients for charge exchange coefficients spline
#endif
#endif
#ifdef AMJUELSPLINES
    ! Atomic rates coefficients
     REAL*8, DIMENSION(9,9):: alpha_iz     ! Coefficients for ionization coefficients spline from EIRENE, (te,ne) grid
     REAL*8, DIMENSION(9,9):: alpha_rec     ! Coefficients for recompination coefficients spline from EIRENE, (te,ne) grid
    real*8, dimension(9,9):: alpha_energy_iz     ! Coefficients for radiation losses due to ionization  spline from EIRENE, (te,ne) grid
    real*8, dimension(9,9):: alpha_energy_rec     ! Coefficients for radiation losses due to recombination  spline from EIRENE, (te,ne) grid
#endif
#ifdef KEQUATION
    ! Coefficients for the k equation
     REAL*8          :: diff_k_min         ! Mininmum diffusion in the k equation
     REAL*8          :: diff_k_max         ! Maximum diffusion in the k equation
     REAL*8          :: k_max              ! Maximum k
#endif
  END TYPE Physics_type

  !*******************************************************
  ! Geometry: type for geometry info
  !*******************************************************
  TYPE Geometry_type
     INTEGER     :: ndim     ! Number of dimensions of the problem (can be different from the ndim of the mesh)
     REAL*8      :: R0       ! Major radius at the magnetic axis
     REAL*8      :: q        ! Safety factor
  END TYPE Geometry_type

  !*******************************************************
  ! Magnetic: type for 3D magnetic perturbation
  !*******************************************************
  TYPE Magnetic_type
     REAL*8          :: amp_rmp            ! amplitude RMP
     INTEGER         :: nbCoils_rmp        ! number coils RMP (full torus)
     REAL*8          :: torElongCoils_rmp  ! Toroidal elongation of RMP coils (rectangular coils)
     INTEGER         :: parite             ! parite RMP (for 2 row, -1 even parite, 1 odd)
     INTEGER         :: nbRow              ! number of rows (1,2 or 3) of RMP coils, default is 2
     REAL*8          :: amp_ripple         ! amplitude Ripple
     INTEGER         :: nbCoils_ripple     ! number coils Ripple (full torus)
     REAL*8          :: triang             ! triangularity (0: None)
     REAL*8          :: ellip              ! ellipticity (1: None)
     REAL*8, POINTER :: coils_rmp(:, :, :) ! Coil coordinates for RMP (nbCoils*4*Discr,start-stop*(xyz)=6,rowNb) (4 for square coils)
     REAL*8, POINTER :: coils_ripple(:, :) ! Coil coordinates for Ripple (nbCoils*Discr,start-stop*(xyz)=6)
  END TYPE Magnetic_type

  !*******************************************************
  ! Switches: type for main code switches
  !*******************************************************
  TYPE Switches_type
     LOGICAL :: axisym ! Is it an axisymmetric simulation?
    ! true =
    ! false =
     LOGICAL :: rmp      ! To activate resonant magnetic perturbation
     LOGICAL :: ripple   ! To activate ripple
     LOGICAL :: ohmicsrc ! Set to TRUE to consider ohmic source of energy
     LOGICAL :: ME       ! Set to TRUE to allow magnetic equilibrium evolution in time
     LOGICAL :: driftdia ! Set to TRUE to consider diamagnetic drift
     LOGICAL :: driftexb ! Set to TRUE to consider ExB drift
     LOGICAL :: steady
     LOGICAL :: read_gmsh ! read a mesh file in gmsh format
     LOGICAL :: readMeshFromSol ! read a mesh file from the solution
     LOGICAL :: set_2d_order ! if read_gmsh = .true., set order
     INTEGER :: order_2d    ! if read_gmsh = .true., set order = .true., what 2d polynomial order?
     LOGICAL :: gmsh2h5
     LOGICAL :: saveMeshSol
     LOGICAL :: time_init ! true if it is a time initialization simulation. The time counter "it" does not increment  (i.e. when the analitical initialisation is not good enough). Used for moving equilibrium (case 59)
     INTEGER :: init     ! 1-init. analy. solution at nodes; 2-L2 projection
    ! Set to TRUE for a steady state computation
    ! Set to FALSE for a transient computation
     INTEGER :: testcase  ! Define the testcase ( which analytical solution, body force, magnetic field...)
     LOGICAL :: psdtime ! Reduce the diffusion every time we reach the steady state
    ! condition (only works if steady=.false.)
     REAL*8  :: diffred ! Reduction factor of the diffusion for psdtime simulation
     REAL*8  :: diffmin ! Minimum value of the diffusion for a psdtime simulation
     INTEGER :: shockcp ! Shock capturing option
     INTEGER :: limrho  ! Add a source for limiting the min value of rho
     INTEGER :: difcor  ! Add diffusion in corners
     INTEGER :: thresh  ! Use a threshold for limiting the min value of rho
    ! (rho-rho*E in case of N-Gamma-Energy model, rho-rho*Ei-rho*Ee in case of N-Gamma-Ti-Te model)
     LOGICAL :: filter  ! Filter solution to avoid oscillation in empty zones
     LOGICAL :: decoup  ! Decouple N-Gamma from Te-Ti (only used for N-Gamma-Ti-Te model)
     LOGICAL :: ckeramp ! Chech the error amplification in the linear system solution (for very ill-conditioned matrices)
     LOGICAL :: saveNR  ! Save solution at each NR iteration
     LOGICAL :: saveTau ! Save tau on faces
     LOGICAL :: fixdPotLim
     LOGICAL :: dirivortcore
     LOGICAL :: dirivortlim
     LOGICAL :: convvort ! consider the convective term in the vorticity equation
     LOGICAL :: bxgradb  ! consider the term in BxGradB in the vorticity equation
     INTEGER :: pertini  ! add perturbation in the initial solution
    ! 1 -add sinusoidal perturbation
    ! 2 -add density blob
     LOGICAL :: logrho   ! solve for the density logarithm instead of density
  END TYPE Switches_type

  !***************************************************************
  ! Paths: type for storing paths to load inputs and sotre outputs
  !***************************************************************
  TYPE Inputs_type
     CHARACTER(len=1000) :: field_path ! where do we read magnetic field from (WEST cases so far)
     CHARACTER(len=1000) :: jtor_path ! where do we read plasma current from (WEST cases so far)
     CHARACTER(len=1000) :: save_folder ! where to save last solution
     LOGICAL             :: field_from_grid !if true, then reads equilibrium file n rectangular grid; if false - on nodes of the mesh
     LOGICAL             :: compute_from_flux ! if components B_R, B_Z are computed from flux or not
     LOGICAL             :: divide_by_2pi     ! correspondng to flux definition if it is needed to divide by 2pi or not
     INTEGER             :: field_dimensions(1:2) ! dimensions of magnetic field files (2D WEST cases so far)
     INTEGER             :: jtor_dimensions(1:2) ! dimensions of magnetic field files (2D WEST cases so far)
  END TYPE Inputs_type

  !*******************************************************
  ! Time: type for the time stepping information
  !*******************************************************
  TYPE Time_type
     REAL*8      :: dt0  ! initial time step
     REAL*8      :: dt   ! current time step
    real*8      :: dt_ME   ! time step from equilibrium !stored as dimensional value
    real*8      :: t_ME   ! time step from equilibrium (to get puff rate) !stored as dimensional value
     REAL*8      :: tfi  ! final time of the simulation
     INTEGER     :: it   ! the number of the current time step
     INTEGER     :: ik   ! same as it but always incrementing (also in case of pseudotime..)
     INTEGER     :: ndt  ! max number of time steps to do in the current session
     INTEGER     :: tsw  ! switch to modify the time step
     INTEGER     :: nts  ! max number of time iterations to do in the current session (only for transient simulations)
     INTEGER     :: tis  ! time integration scheme
    ! 1 - first order
    ! 2 - second order
     REAL*8      :: t    ! time of the simulation (initialized to finish time of previous simulation if restart, to 0 if new simulation)
  END TYPE Time_type

  !*******************************************************
  ! Numerics: type for numeric scheme parameters
  !*******************************************************
  TYPE Numeric_type
     INTEGER        :: nrp      ! Max number of Newton-Raphson iterations
     REAL*8         :: tNR      ! Tolerance of the Newton-Raphson scheme
     REAL*8         :: tTM      ! Tolerance for the steady state achievement
     REAL*8         :: div      ! Divergence detector
     REAL*8         :: tau(1:5) ! Stabilization parameter for each equation (4 values max for now...)
     REAL*8         :: sc_coe   ! Shock capturing coefficient
     REAL*8         :: sc_sen   ! Shock capturing sensibility
     REAL*8         :: minrho   ! Value of rho to start applying limiting
     REAL*8         :: so_coe   ! Source coefficient for limiting rho
     REAL*8         :: df_coe   ! Diffusion coefficient for limiting rho
     REAL*8         :: dc_coe   ! Diffusion coefficient in corners
     REAL*8         :: thr      ! Threshold to limit rho
     REAL*8         :: thrpre   ! Threshold to limit pressure
     INTEGER        :: stab     ! Stabilization type
    ! 1 - constant tau (one for each equation) in the whole domain
    ! 2 -
    ! 3 -
     REAL*8         :: dumpnr   ! dumping factor for Newton-Raphson. 0<dumpnr<1
     REAL*8         :: dumpnr_min   ! dumping factor minimum for Newton-Raphson. 0<dumpnr<1
     REAL*8         :: dumpnr_max   ! dumping factor maximum for Newton-Raphson. 0<dumpnr<1
     REAL*8         :: dumpnr_width   ! dumping factor width of hyperbolic tangential for Newton-Raphson. 0<dumpnr<1
     REAL*8         :: dumpnr_n0   ! dumping factor x0 for hyperbolic tangential for Newton-Raphson. 0<dumpnr<1
     INTEGER        :: ntor     ! Number of elements in the toroidal direction
     INTEGER        :: ptor     ! Polynomial degree in the toroidal direction
     REAL*8         :: tmax     ! Max extention in the toroidal direction
     INTEGER        :: npartor  ! Number of MPI divisions in the toroidal direction
     INTEGER        :: bohmtypebc ! Implementation of the Bohm bc for Gamma
     REAL*8         :: exbdump ! Dumping for ExB drifts
  END TYPE Numeric_type

  !*******************************************************
  ! Adaptivity: parameters for the adaptivity procedures
  !*******************************************************
  TYPE adaptivity_type
     LOGICAL :: adaptivity       ! Adaptivity enabled or not (1 or 0)
     INTEGER :: shockcp_adapt
     INTEGER :: evaluator        ! 1,2,3 (indicator, estimator, both)
     REAL*8  :: thr_ind          ! Threshold  for the detection of oscillations (indicator)
     INTEGER :: quant_ind        ! 1,2,3 (physical variable(s), gradient of variable(s), both)
     INTEGER :: n_quant_ind      ! 1,2,3 for iso (n,u, Mach). 1,2,3,4,5 for aniso (n,u,Ti,Te,Mach), 10 for all of them
     REAL*8  :: tol_est          ! Relative difference between the p+1 and p solution (estimator)
     INTEGER :: param_est        ! 1,2,3 ... physical variable used in the estimator (-1 for all of them)
     INTEGER :: difference       ! 0,1 (relative, absolute)
     LOGICAL :: time_adapt       ! refine at times steps (pseudo or not)
     LOGICAL :: NR_adapt         ! refine at NR steps
     INTEGER :: freq_t_adapt     ! Frequency of time refinement, only if time_adapt = .true.
     INTEGER :: freq_NR_adapt    ! Frequency of NR refinement, only if NR_adapt = .true.
     LOGICAL :: div_adapt        ! refine if NR divergence
     LOGICAL :: rest_adapt       ! call adaptivity at the very beginning, only for restart simulations
     LOGICAL :: osc_adapt        ! refine if oscillations are lower than threshold
     REAL*8  :: osc_tol          ! refine if oscillations are lower than this threshold, only if osc_adapt = .true.
     REAL*8  :: osc_check        ! save the solution as checkpoint if the maximum value of oscillations are lower than this threshold
  END TYPE adaptivity_type

  !*******************************************************
  ! Utilities: type for printing/debugging/saving...
  !*******************************************************
  TYPE Utils_type
     INTEGER :: printint       ! Integer for printing
     LOGICAL :: timing         ! Timing of the code
     INTEGER :: freqdisp       ! Frequency of results display
     INTEGER :: freqsave       ! Frequency of solution save
  END TYPE Utils_type

  !*******************************************************
  ! Linear system solver parameters
  !*******************************************************
  TYPE Lssolver_type
     INTEGER           :: sollib    ! Solver library to be used
    ! 1-Pastix
    ! 2-PSBLAS
     ! 3-PETSc

     LOGICAL           :: timing    ! timing of the linear solver
     ! Parameters relative to the library PETSc
     INTEGER           :: kspitrace     ! Display convergence at each iteration
     REAL*8            :: rtol       ! Relative tolerance
     REAL*8            :: atol       ! Absolute tolerance
     INTEGER           :: kspitmax      ! Max number of iterations
     CHARACTER(len=20) :: kspmethd     ! Krylov method (see list on the library manual)
     LOGICAL           :: igz        ! Set initial guess of the iterative method to zeros
     INTEGER           :: rprecond   ! Recompute preconditioner
     INTEGER           :: Nrprecond  ! Recompute preconditioner each Rrprecond iterations
     INTEGER           :: kspnorm    ! norm type to be used
     CHARACTER(len=20) :: ksptype    ! Krylov solver type
     CHARACTER(len=20) :: pctype     ! Preconditioner type
     INTEGER           :: gmresres   ! Restart value for GMRES
     INTEGER           :: mglevels   ! Number of levels for the MultiGrid preconditioner
     INTEGER           :: mgtypeform   ! Form type of the MultiGrid preconditioner


     ! Parameters relative to the library PSBLAS
     CHARACTER(len=20) :: kmethd    ! Krylov method (see list on the library manual)
     INTEGER           :: istop     ! Stopping criterion (see spec on the library manual)
     INTEGER           :: itmax     ! Max number of iterations
     INTEGER           :: itrace    ! Display convergence at each iteration
     INTEGER           :: rest      ! Restart
     REAL*8            :: tol       ! Stopping tolerance
     CHARACTER(len=20) :: ptype     ! Preconditioner type
    ! Parameters relative to the library MLD2P4
    ! First smoother / 1-lev preconditioner
     CHARACTER(len=20) :: smther       ! smoother type
     INTEGER           :: jsweeps      ! (pre-)smoother / 1-lev prec sweeps
     INTEGER           :: novr         ! Number of olverlap layers, for Additive Schwartz only
     CHARACTER(len=20) :: restr        ! restriction  over application for Additive Schwartz only
     CHARACTER(len=20) :: prol         ! Type of prolongation operator for Additive Schwartz only
     CHARACTER(len=20) :: solve        ! local subsolver
     INTEGER           :: fill         ! fill-in level p of the ILU factorizations
     REAL              :: thr          ! threshold for ILUT
    ! Second smoother/ AMG post-smoother (if NONE ignored in main)
     CHARACTER(len=20) :: smther2      ! smoother type
     INTEGER           :: jsweeps2     ! (post-)smoother sweeps
     INTEGER           :: novr2        ! number of overlap layers
     CHARACTER(len=20) :: restr2       ! restriction  over application of AS
     CHARACTER(len=20) :: prol2        ! prolongation over application of AS
     CHARACTER(len=20) :: solve2       ! local subsolver
     INTEGER           :: fill2        ! fill-in for incomplete LU
     REAL              :: thr2         ! threshold for ILUT
    ! general AMG data
     CHARACTER(len=20) :: mlcycle      ! multi-level cycle
     INTEGER           :: outer_sweeps ! number of multilevel cycles
     INTEGER           :: maxlevs      ! Maximum number of levels
     INTEGER           :: csize        ! Coarse size threshold
    ! aggregation
     REAL              :: mncrratio    ! Minimum coarsening ratio
     REAL              :: athres       ! Aggregation threshold
     CHARACTER(len=20) :: aggr_prol    ! Prolongator used by the aggregation algorithm
     CHARACTER(len=20) :: par_aggr_alg ! Parallel aggregation algorithm
     CHARACTER(len=20) :: aggr_ord     ! Initial ordering of indices for the aggregation algorithm
     CHARACTER(len=20) :: aggr_filter  ! Matrix used in computing the smoothed prolongator
    ! coasest-level solver
     CHARACTER(len=20) :: csolve       ! coarsest-lev solver
     CHARACTER(len=20) :: csbsolve     ! coarsest-lev solver
     CHARACTER(len=20) :: cmat         ! coarsest mat layout
     INTEGER           :: cfill        ! fill-in for incompl LU
     REAL              :: cthres       ! Threshold for ILUT
     INTEGER           :: cjswp        ! sweeps for GS/JAC subsolver
  END TYPE Lssolver_type

  !**********************************************************
  ! Solution: contains the solution at the current time step
  !**********************************************************
  TYPE Sol_type
     REAL*8, POINTER :: u(:) => NULL()  ! Elemental solution
     REAL*8, POINTER :: u_conv(:) => NULL()
     REAL*8, POINTER :: u_init(:) => NULL()
     REAL*8, POINTER :: u0_init(:) => NULL()
     REAL*8, POINTER :: u_tilde(:) => NULL()  ! Face solution
     REAL*8, POINTER :: u_tilde0(:) => NULL()  ! Face solution
     REAL*8, POINTER :: q(:) => NULL()  ! Elemental solution for the gradient
     REAL*8, POINTER :: q_conv(:) => NULL()
     REAL*8, POINTER :: q_init(:) => NULL()
     REAL*8, POINTER :: q0_init(:) => NULL()
     REAL*8, ALLOCATABLE :: u0(:, :)           ! Elemental solution at previous time steps
     REAL*8, ALLOCATABLE :: tres(:)           ! Time residual
     REAL*8, ALLOCATABLE :: time(:)           ! Time evolution
     INTEGER             :: Nt                ! Number of time-steps
  END TYPE Sol_type

  !**********************************************************
  ! Simulation parameters: for saving purpose
  !**********************************************************
  TYPE Simulationparams_type
     CHARACTER(len=50) :: model
     INTEGER   :: Ndim
     INTEGER   :: Neq
     REAL, ALLOCATABLE :: consvar_refval(:)
     REAL, ALLOCATABLE :: physvar_refval(:)
    ! Reference values
     REAL*8    :: refval_length
     REAL*8    :: refval_mass
     REAL*8    :: refval_charge
     REAL*8    :: refval_time
     REAL*8    :: refval_temperature
     REAL*8    :: refval_density
     REAL*8    :: refval_neutral
#ifdef KEQUATION
     REAL*8    :: refval_k
#endif
     REAL*8    :: refval_speed
     REAL*8    :: refval_potential
     REAL*8    :: refval_vorticity
     REAL*8    :: refval_magfield
     REAL*8    :: refval_current
     REAL*8    :: refval_diffusion
     REAL*8    :: refval_momentum
     REAL*8    :: refval_specpress
     REAL*8    :: refval_specenergy
     REAL*8    :: refval_specenergydens
    ! Adimesional isothermal compressibility coefficient
     REAL*8    :: compress_coeff

    ! Dimensions used in the simulation
     CHARACTER(len=20)    :: refval_length_dimensions
     CHARACTER(len=20)    :: refval_mass_dimensions
     CHARACTER(len=20)    :: refval_charge_dimensions
     CHARACTER(len=20)    :: refval_time_dimensions
     CHARACTER(len=20)    :: refval_temperature_dimensions
     CHARACTER(len=20)    :: refval_density_dimensions
     CHARACTER(len=20)    :: refval_neutral_dimensions
#ifdef KEQUATION
     CHARACTER(len=20)    :: refval_k_dimensions
#endif
     CHARACTER(len=20)    :: refval_speed_dimensions
     CHARACTER(len=20)    :: refval_potential_dimensions
     CHARACTER(len=20)    :: refval_vorticity_dimensions
     CHARACTER(len=20)    :: refval_magfield_dimensions
     CHARACTER(len=20)    :: refval_current_dimensions
     CHARACTER(len=20)    :: refval_diffusion_dimensions
     CHARACTER(len=20)    :: refval_momentum_dimensions
     CHARACTER(len=20)    :: refval_specpress_dimensions
     CHARACTER(len=20)    :: refval_specenergy_dimensions
     CHARACTER(len=20)    :: refval_specenergydens_dimensions

    ! Physical parameters used in the computation
     REAL*8    :: a
     REAL*8    :: Mref
     REAL*8    :: c1
     REAL*8    :: c2
     REAL*8    :: diff_pari
     REAL*8    :: diff_pare
     REAL*8    :: diff_n
     REAL*8    :: diff_u
     REAL*8    :: diff_e
     REAL*8    :: diff_ee
  END TYPE Simulationparams_type

  !**********************************************************
  ! Elemental matrices: type to store the elemental matrices
  ! used during the computation
  !**********************************************************
  TYPE :: elmat_type
     REAL*8, ALLOCATABLE :: iAqq(:, :, :)
     REAL*8, ALLOCATABLE :: Aqu(:, :, :)
     REAL*8, ALLOCATABLE :: Aql(:, :, :)
     REAL*8, ALLOCATABLE :: Auq(:, :, :)
     REAL*8, ALLOCATABLE :: Auu(:, :, :)
     REAL*8, ALLOCATABLE :: Aul(:, :, :)
     REAL*8, ALLOCATABLE :: Alq(:, :, :)
     REAL*8, ALLOCATABLE :: Alu(:, :, :)
     REAL*8, ALLOCATABLE :: ALL(:, :, :)
     REAL*8, ALLOCATABLE :: Aql_dir(:, :)
     REAL*8, ALLOCATABLE :: Aul_dir(:, :)

     REAL*8, ALLOCATABLE :: M(:, :, :)
     REAL*8, ALLOCATABLE :: Cv(:, :, :)
     REAL*8, ALLOCATABLE :: H(:, :, :)
     REAL*8, ALLOCATABLE :: Hdir(:, :)
     REAL*8, ALLOCATABLE :: D(:, :, :)
     REAL*8, ALLOCATABLE :: E(:, :, :)
     REAL*8, ALLOCATABLE :: Edir(:, :)
     REAL*8, ALLOCATABLE :: S(:, :)
     REAL*8, ALLOCATABLE :: UU(:, :, :)
     REAL*8, ALLOCATABLE :: U0(:, :)
     REAL*8, ALLOCATABLE :: Hf(:, :, :)
     REAL*8, ALLOCATABLE :: Df(:, :, :)
     REAL*8, ALLOCATABLE :: Ef(:, :, :)
     REAL*8, ALLOCATABLE :: fH(:, :)
     REAL*8, ALLOCATABLE :: B(:, :, :)
     REAL*8, ALLOCATABLE :: C(:, :, :)
     REAL*8, ALLOCATABLE :: Cdir(:, :)
     REAL*8, ALLOCATABLE :: P(:, :, :)
     REAL*8, ALLOCATABLE :: G(:, :, :)
     REAL*8, ALLOCATABLE :: IL(:, :, :)
     REAL*8, ALLOCATABLE :: Lf(:, :, :)
     REAL*8, ALLOCATABLE :: Qf(:, :, :)
     REAL*8, ALLOCATABLE :: LL(:, :, :)
     REAL*8, ALLOCATABLE :: L0(:, :)
    ! Limiting rho
     REAL*8, ALLOCATABLE :: S_lrho(:, :)
     REAL*8, ALLOCATABLE :: P_lrho(:, :, :)
    ! Shock capturing
     REAL*8, ALLOCATABLE :: P_sc(:, :, :)
     REAL*8, ALLOCATABLE :: Lf_sc(:, :, :)
     REAL*8, ALLOCATABLE :: TQ(:, :, :)
     REAL*8, ALLOCATABLE :: TQhf(:, :, :)
     REAL*8, ALLOCATABLE :: Tfhf(:, :)
     REAL*8, ALLOCATABLE :: Tf(:, :)
     REAL*8, ALLOCATABLE :: Thdir(:, :)

  END TYPE elmat_type

  !**********************************************************
  ! Timing: structure used to store execution times
  !**********************************************************
  TYPE Timing_type
     INTEGER :: cks1, clock_rate1, cke1, clock_start1, clock_end1
     INTEGER :: cks2, clock_rate2, cke2, clock_start2, clock_end2
     REAL*8  :: time_start1, time_finish1, tps1, tpe1
     REAL*8  :: time_start2, time_finish2, tps2, tpe2
     REAL*8  :: cputpre, cputmap, cputass, cputbcd, cputsol, cputjac, cputglb, cputadapt
     REAL*8  :: runtpre, runtmap, runtass, runtbcd, runtsol, runtjac, runtglb, runtadapt
     REAL(8) :: clstime1, clstime2, clstime3, clstime4, clstime5, clstime6
     REAL(8) :: rlstime1, rlstime2, rlstime3, rlstime4, rlstime5, rlstime6
     REAL(8) :: cputcom
     REAL(8) :: runtcom
  END TYPE Timing_type


CONTAINS

  SUBROUTINE set_boundary_flag_names()
    bc_flag_type(1) = 'Tb_Dirichlet'
    bc_flag_type(2) = 'Tb_LEFT'
    bc_flag_type(3) = 'Tb_RIGHT'
    bc_flag_type(4) = 'Tb_UP'
    bc_flag_type(5) = 'Tb_PUMP'
    bc_flag_type(6) = 'Tb_PUFF'
    bc_flag_type(7) = 'Tb_LIM'
    bc_flag_type(8) = 'Tb_IN'
    bc_flag_type(9) = 'Tb_OUT'
    bc_flag_type(10) = 'Tb_ULIM'

    bc_flag_name(1) = 'Dirichlet strong form'
    bc_flag_name(2) = 'Dirichlet weak form'
    bc_flag_name(3) = 'Dirichlet weak form using old values of the variables'
    bc_flag_name(4) = 'Transmission boundary conditions'
    bc_flag_name(5) = 'Dirichlet in some equations and Neumann in others '
    bc_flag_name(6) = 'ITER core BC: plasma flux = plasma flux, energy flux imposed '
    bc_flag_name(20) = 'Inlet-Outlet'
    bc_flag_name(30) = 'Neumann homogeneus'
    bc_flag_name(50) = 'Bohm'
    bc_flag_name(51) = 'Bohm for the velocity, Neumann homogeneous for the other variables'
    bc_flag_name(52) = 'Simplified Bohm bc: normal Bohm for the velocity, Grad//T=0 for Ti, Te'
    bc_flag_name(55) = 'Bohm for plasma, pump for neutrals'
    bc_flag_name(56) = 'Bohm for plasma, puff for neutrals'
    bc_flag_name(60) = 'Slip wall'
    bc_flag_name(70) = 'Periodic'

  END SUBROUTINE set_boundary_flag_names

  SUBROUTINE deep_copy_mesh_struct(Mesh1, Mesh2)
    TYPE(Mesh_type), INTENT(IN)               :: Mesh1
    TYPE(Mesh_type), INTENT(OUT)              :: Mesh2

    Mesh2%Ndim = Mesh1%Ndim
    Mesh2%Nnodes = Mesh1%Nnodes
    Mesh2%Nnodesperelem = Mesh1%Nnodesperelem
    Mesh2%Nnodesperface = Mesh1%Nnodesperface
    Mesh2%Nelems = Mesh1%Nelems
    Mesh2%Nfaces = Mesh1%Nfaces
    Mesh2%Nextfaces = Mesh1%Nextfaces
    Mesh2%Nintfaces = Mesh1%Nintfaces
    Mesh2%elemType = Mesh1%elemType
    Mesh2%ndir   = Mesh1%ndir
    Mesh2%ukf = Mesh1%ukf

    IF(ASSOCIATED(Mesh1%T)) THEN
       ALLOCATE(Mesh2%T(SIZE(Mesh1%T, 1), SIZE(Mesh1%T, 2)))
       Mesh2%T = Mesh1%T
    ENDIF

    IF(ASSOCIATED(Mesh1%Tlin)) THEN
       ALLOCATE(Mesh2%Tlin(SIZE(Mesh1%Tlin, 1), SIZE(Mesh1%Tlin, 2)))
       Mesh2%Tlin = Mesh1%Tlin
    ENDIF

    IF(ASSOCIATED(Mesh1%Tb)) THEN
       ALLOCATE(Mesh2%Tb(SIZE(Mesh1%Tb, 1), SIZE(Mesh1%Tb, 2)))
       Mesh2%Tb = Mesh1%Tb
    ENDIF

    IF(ASSOCIATED(Mesh1%boundaryFlag)) THEN
       ALLOCATE(Mesh2%boundaryFlag(SIZE(Mesh1%boundaryFlag)))
       Mesh2%boundaryFlag = Mesh1%boundaryFlag
    ENDIF

    IF(ALLOCATED(Mesh1%F)) THEN
       ALLOCATE(Mesh2%F(SIZE(Mesh1%F, 1), SIZE(Mesh1%F, 2)))
       Mesh2%F = Mesh1%F
    ENDIF

    IF(ALLOCATED(Mesh1%N)) THEN
       ALLOCATE(Mesh2%N(SIZE(Mesh1%N, 1), SIZE(Mesh1%N, 2)))
       Mesh2%N = Mesh1%N
    ENDIF

    IF(ALLOCATED(Mesh1%faces)) THEN
       ALLOCATE(Mesh2%faces(SIZE(Mesh1%faces, 1), SIZE(Mesh1%faces, 2), SIZE(Mesh1%faces, 3)))
       Mesh2%faces = Mesh1%faces
    ENDIF

    IF(ALLOCATED(Mesh1%extfaces)) THEN
       ALLOCATE(Mesh2%extfaces(SIZE(Mesh1%extfaces, 1), SIZE(Mesh1%extfaces, 2)))
       Mesh2%extfaces = Mesh1%extfaces
    ENDIF

    IF(ALLOCATED(Mesh1%intfaces)) THEN
       ALLOCATE(Mesh2%intfaces(SIZE(Mesh1%intfaces, 1), SIZE(Mesh1%intfaces, 2)))
       Mesh2%intfaces = Mesh1%intfaces
    ENDIF

    IF(ALLOCATED(Mesh1%flipface)) THEN
       ALLOCATE(Mesh2%flipface(SIZE(Mesh1%flipface, 1), SIZE(Mesh1%flipface, 2)))
       Mesh2%flipface = Mesh1%flipface
    ENDIF

    IF(ALLOCATED(Mesh1%Fdir)) THEN
       ALLOCATE(Mesh2%Fdir(SIZE(Mesh1%Fdir, 1), SIZE(Mesh1%Fdir, 2)))
       Mesh2%Fdir = Mesh1%Fdir
    ENDIF

    IF(ALLOCATED(Mesh1%periodic_faces)) THEN
       ALLOCATE(Mesh2%periodic_faces(SIZE(Mesh1%periodic_faces)))
       Mesh2%periodic_faces = Mesh1%periodic_faces
    ENDIF

    IF(ALLOCATED(Mesh1%Diric)) THEN
       ALLOCATE(Mesh2%Diric(SIZE(Mesh1%Diric)))
       Mesh2%Diric = Mesh1%Diric
    ENDIF

    IF(ALLOCATED(Mesh1%numberbcs)) THEN
       ALLOCATE(Mesh2%numberbcs(SIZE(Mesh1%numberbcs)))
       Mesh2%numberbcs = Mesh1%numberbcs
    ENDIF

    IF(ALLOCATED(Mesh1%elemSize)) THEN
       ALLOCATE(Mesh2%elemSize(SIZE(Mesh1%elemSize)))
       Mesh2%elemSize = Mesh1%elemSize
    ENDIF

    IF(ASSOCIATED(Mesh1%X)) THEN
       ALLOCATE(Mesh2%X(SIZE(Mesh1%X, 1), SIZE(Mesh1%X, 2)))
       Mesh2%X = Mesh1%X
    ENDIF

    IF(ASSOCIATED(Mesh1%toroidal)) THEN
       ALLOCATE(Mesh2%toroidal(SIZE(Mesh1%toroidal)))
       Mesh2%toroidal = Mesh1%toroidal
    ENDIF

    IF(ALLOCATED(Mesh1%flag_elems_rho)) THEN
       ALLOCATE(Mesh2%flag_elems_rho(SIZE(Mesh1%flag_elems_rho)))
       Mesh2%flag_elems_rho = Mesh1%flag_elems_rho
    ENDIF

    IF(ALLOCATED(Mesh1%flag_elems_sc)) THEN
       ALLOCATE(Mesh2%flag_elems_sc(SIZE(Mesh1%flag_elems_sc)))
       Mesh2%flag_elems_sc = Mesh1%flag_elems_sc
    ENDIF

    IF(ALLOCATED(Mesh1%minrho_elems)) THEN
       ALLOCATE(Mesh2%minrho_elems(SIZE(Mesh1%minrho_elems)))
       Mesh2%minrho_elems = Mesh1%minrho_elems
    ENDIF

    IF(ALLOCATED(Mesh1%sour_elems)) THEN
       ALLOCATE(Mesh2%sour_elems(SIZE(Mesh1%sour_elems)))
       Mesh2%sour_elems = Mesh1%sour_elems
    ENDIF

    IF(ALLOCATED(Mesh1%diff_elems)) THEN
       ALLOCATE(Mesh2%diff_elems(SIZE(Mesh1%diff_elems)))
       Mesh2%diff_elems = Mesh1%diff_elems
    ENDIF

    IF(ALLOCATED(Mesh1%scdiff_nodes)) THEN
       ALLOCATE(Mesh2%scdiff_nodes(SIZE(Mesh1%scdiff_nodes, 1), SIZE(Mesh1%scdiff_nodes, 2)))
       Mesh2%scdiff_nodes = Mesh1%scdiff_nodes
    ENDIF

    Mesh2%Nnodes_toroidal  = Mesh1%Nnodes_toroidal
    Mesh2%xmax = Mesh1%xmax
    Mesh2%xmin = Mesh1%xmin
    Mesh2%ymax = Mesh1%ymax
    Mesh2%ymin = Mesh1%ymin
    Mesh2%puff_area = Mesh1%puff_area
    Mesh2%core_area = Mesh1%core_area

#ifdef PARALL
    Mesh2%nghostfaces = Mesh1%nghostfaces
    Mesh2%nghostelems  = Mesh1%nghostelems
    Mesh2%Nel_glob = Mesh1%Nel_glob
    Mesh2%Nfa_glob = Mesh1%Nfa_glob
    Mesh2%Nno_glob = Mesh1%Nno_glob
    Mesh2%Ndir_glob = Mesh1%Ndir_glob
    Mesh2%Ngho_glob = Mesh1%Ngho_glob

    IF (ASSOCIATED(Mesh1%loc2glob_fa)) THEN
       ALLOCATE(Mesh2%loc2glob_fa(SIZE(Mesh1%loc2glob_fa)))
       Mesh2%loc2glob_fa = Mesh1%loc2glob_fa
    ENDIF
    IF (ASSOCIATED(Mesh1%loc2glob_el)) THEN
       ALLOCATE(Mesh2%loc2glob_el(SIZE(Mesh1%loc2glob_el)))
       Mesh2%loc2glob_el = Mesh1%loc2glob_el
    ENDIF
    IF (ASSOCIATED(Mesh1%loc2glob_nodes)) THEN
       ALLOCATE(Mesh2%loc2glob_nodes(SIZE(Mesh1%loc2glob_nodes)))
       Mesh2%loc2glob_nodes = Mesh1%loc2glob_nodes
    ENDIF
    IF (ASSOCIATED(Mesh1%ghostfaces)) THEN
       ALLOCATE(Mesh2%ghostfaces(SIZE(Mesh1%ghostfaces)))
       Mesh2%ghostfaces = Mesh1%ghostfaces
    ENDIF
    IF (ASSOCIATED(Mesh1%ghostelems)) THEN
       ALLOCATE(Mesh2%ghostelems(SIZE(Mesh1%ghostelems)))
       Mesh2%ghostelems = Mesh1%ghostelems
    ENDIF
    IF (ASSOCIATED(Mesh1%ghostflp)) THEN
       ALLOCATE(Mesh2%ghostflp(SIZE(Mesh1%ghostflp)))
       Mesh2%ghostflp = Mesh1%ghostflp
    ENDIF
    IF (ASSOCIATED(Mesh1%ghostloc)) THEN
       ALLOCATE(Mesh2%ghostloc(SIZE(Mesh1%ghostloc)))
       Mesh2%ghostloc = Mesh1%ghostloc
    ENDIF
    IF (ASSOCIATED(Mesh1%ghostpro)) THEN
       ALLOCATE(Mesh2%ghostpro(SIZE(Mesh1%ghostpro)))
       Mesh2%ghostpro = Mesh1%ghostpro
    ENDIF
    IF (ASSOCIATED(Mesh1%ghelsloc)) THEN
       ALLOCATE(Mesh2%ghelsloc(SIZE(Mesh1%ghelsloc)))
       Mesh2%ghelsloc = Mesh1%ghelsloc
    ENDIF
    IF (ASSOCIATED(Mesh1%ghelspro)) THEN
       ALLOCATE(Mesh2%ghelspro(SIZE(Mesh1%ghelspro)))
       Mesh2%ghelspro = Mesh1%ghelspro
    ENDIF
    IF (ALLOCATED(Mesh1%fc2sd)) THEN
       ALLOCATE(Mesh2%fc2sd(SIZE(Mesh1%fc2sd)))
       Mesh2%fc2sd = Mesh1%fc2sd
    ENDIF
    IF (ALLOCATED(Mesh1%pr2sd)) THEN
       ALLOCATE(Mesh2%pr2sd(SIZE(Mesh1%pr2sd)))
       Mesh2%pr2sd = Mesh1%pr2sd
    ENDIF
    IF (ALLOCATED(Mesh1%fc2rv)) THEN
       ALLOCATE(Mesh2%fc2rv(SIZE(Mesh1%fc2rv)))
       Mesh2%fc2rv = Mesh1%fc2rv
    ENDIF
    IF (ALLOCATED(Mesh1%pr2rv)) THEN
       ALLOCATE(Mesh2%pr2rv(SIZE(Mesh1%pr2rv)))
       Mesh2%pr2rv = Mesh1%pr2rv
    ENDIF
    IF (ALLOCATED(Mesh1%el2sd)) THEN
       ALLOCATE(Mesh2%el2sd(SIZE(Mesh1%el2sd)))
       Mesh2%el2sd = Mesh1%el2sd
    ENDIF
    IF (ALLOCATED(Mesh1%el2rv)) THEN
       ALLOCATE(Mesh2%el2rv(SIZE(Mesh1%el2rv)))
       Mesh2%el2rv = Mesh1%el2rv
    ENDIF
    IF (ALLOCATED(Mesh1%pe2rv)) THEN
       ALLOCATE(Mesh2%pe2rv(SIZE(Mesh1%pe2rv)))
       Mesh2%pe2rv = Mesh1%pe2rv
    ENDIF


#endif

  ENDSUBROUTINE deep_copy_mesh_struct

  SUBROUTINE deep_copy_refel_struct(refEl1, refEl2)
    TYPE(Reference_element_type), INTENT(IN)               :: refEl1
    TYPE(Reference_element_type), INTENT(OUT)              :: refEl2

    refEl2%ndim = refEl1%ndim
    refEl2%elemType = refEl1%elemType
    refEl2%Nvertices = refEl1%Nvertices
    refEl2%Ndeg = refEl1%Ndeg
    refEl2%Nnodes3D = refEl1%Nnodes3D
    refEl2%Nnodes2D = refEl1%Nnodes2D
    refEl2%Nnodes1D = refEl1%Nnodes1D
    refEl2%Nfaces = refEl1%Nfaces
    refEl2%Nbords = refEl1%Nbords
    refEl2%Nextnodes   = refEl1%Nextnodes
    refEl2%Nfacenodes = refEl1%Nfacenodes
    refEl2%Nfacenodeslin = refEl1%Nfacenodeslin
    refEl2%Nbordnodes = refEl1%Nbordnodes
    refEl2%Ninnernodes = refEl1%Ninnernodes
    refEl2%Ninnernodesface = refEl1%Ninnernodesface
    refEl2%NGauss1D = refEl1%NGauss1D
    refEl2%NGauss2D = refEl1%NGauss2D
    refEl2%NGauss3D = refEl1%NGauss3D
    refEl2%Nfl = refEl1%Nfl
    refEl2%Ngl = refEl1%Ngl
    refEl2%Nft = refEl1%Nft

    IF(ALLOCATED(refEl1%Face_nodes)) THEN
       ALLOCATE(refEl2%Face_nodes, mold = refEl1%Face_nodes)
       refEl2%Face_nodes = refEl1%Face_nodes
    ENDIF
    IF(ALLOCATED(refEl1%Bord_nodes)) THEN
       ALLOCATE(refEl2%Bord_nodes, mold = refEl1%Bord_nodes)
       refEl2%Bord_nodes = refEl1%Bord_nodes
    ENDIF
    IF(ALLOCATED(refEl1%inner_nodes)) THEN
       ALLOCATE(refEl2%inner_nodes, mold = refEl1%inner_nodes)
       refEl2%inner_nodes = refEl1%inner_nodes
    ENDIF
    IF(ALLOCATED(refEl1%inner_nodes_face)) THEN
       ALLOCATE(refEl2%inner_nodes_face, mold = refEl1%inner_nodes_face)
       refEl2%inner_nodes_face = refEl1%inner_nodes_face
    ENDIF
    IF(ASSOCIATED(refEl1%coord3d)) THEN
       ALLOCATE(refEl2%coord3d, mold = refEl1%coord3d)
       refEl2%coord3d = refEl1%coord3d
    ENDIF
    IF(ASSOCIATED(refEl1%coord2d)) THEN
       ALLOCATE(refEl2%coord2d, mold = refEl1%coord2d)
       refEl2%coord2d = refEl1%coord2d
    ENDIF
    IF(ASSOCIATED(refEl1%coord1d)) THEN
       ALLOCATE(refEl2%coord1d, mold = refEl1%coord1d)
       refEl2%coord1d = refEl1%coord1d
    ENDIF
    IF(ALLOCATED(refEl1%gauss_points3D)) THEN
       ALLOCATE(refEl2%gauss_points3D, mold = refEl1%gauss_points3D)
       refEl2%gauss_points3D = refEl1%gauss_points3D
    ENDIF
    IF(ALLOCATED(refEl1%gauss_points2D)) THEN
       ALLOCATE(refEl2%gauss_points2D, mold = refEl1%gauss_points2D)
       refEl2%gauss_points2D = refEl1%gauss_points2D
    ENDIF
    IF(ALLOCATED(refEl1%gauss_points1D)) THEN
       ALLOCATE(refEl2%gauss_points1D, mold = refEl1%gauss_points1D)
       refEl2%gauss_points1D = refEl1%gauss_points1D
    ENDIF
    IF(ALLOCATED(refEl1%gauss_weights3D)) THEN
       ALLOCATE(refEl2%gauss_weights3D, mold = refEl1%gauss_weights3D)
       refEl2%gauss_weights3D = refEl1%gauss_weights3D
    ENDIF
    IF(ALLOCATED(refEl1%gauss_weights2D)) THEN
       ALLOCATE(refEl2%gauss_weights2D, mold = refEl1%gauss_weights2D)
       refEl2%gauss_weights2D = refEl1%gauss_weights2D
    ENDIF
    IF(ALLOCATED(refEl1%gauss_weights1D)) THEN
       ALLOCATE(refEl2%gauss_weights1D, mold = refEl1%gauss_weights1D)
       refEl2%gauss_weights1D = refEl1%gauss_weights1D
    ENDIF
    IF(ALLOCATED(refEl1%N3D)) THEN
       ALLOCATE(refEl2%N3D, mold = refEl1%N3D)
       refEl2%N3D = refEl1%N3D
    ENDIF
    IF(ALLOCATED(refEl1%Nxi3D)) THEN
       ALLOCATE(refEl2%Nxi3D, mold = refEl1%Nxi3D)
       refEl2%Nxi3D = refEl1%Nxi3D
    ENDIF
    IF(ALLOCATED(refEl1%Neta3D)) THEN
       ALLOCATE(refEl2%Neta3D, mold = refEl1%Neta3D)
       refEl2%Neta3D = refEl1%Neta3D
    ENDIF
    IF(ALLOCATED(refEl1%Nzeta3D)) THEN
       ALLOCATE(refEl2%Nzeta3D, mold = refEl1%Nzeta3D)
       refEl2%Nzeta3D = refEl1%Nzeta3D
    ENDIF
    IF(ALLOCATED(refEl1%N2D)) THEN
       ALLOCATE(refEl2%N2D, mold = refEl1%N2D)
       refEl2%N2D = refEl1%N2D
    ENDIF
    IF(ALLOCATED(refEl1%Nxi2D)) THEN
       ALLOCATE(refEl2%Nxi2D, mold = refEl1%Nxi2D)
       refEl2%Nxi2D = refEl1%Nxi2D
    ENDIF
    IF(ALLOCATED(refEl1%Neta2D)) THEN
       ALLOCATE(refEl2%Neta2D, mold = refEl1%Neta2D)
       refEl2%Neta2D = refEl1%Neta2D
    ENDIF
    IF(ALLOCATED(refEl1%N1D)) THEN
       ALLOCATE(refEl2%N1D, mold = refEl1%N1D)
       refEl2%N1D = refEl1%N1D
    ENDIF
    IF(ALLOCATED(refEl1%Nxi1D)) THEN
       ALLOCATE(refEl2%Nxi1D, mold = refEl1%Nxi1D)
       refEl2%Nxi1D = refEl1%Nxi1D
    ENDIF
    IF(ALLOCATED(refEl1%Nlin)) THEN
       ALLOCATE(refEl2%Nlin, mold = refEl1%Nlin)
       refEl2%Nlin = refEl1%Nlin
    ENDIF
    IF(ALLOCATED(refEl1%sFTF)) THEN
       ALLOCATE(refEl2%sFTF, mold = refEl1%sFTF)
       refEl2%sFTF = refEl1%sFTF
    ENDIF
    IF(ALLOCATED(refEl1%faceNodes3)) THEN
       ALLOCATE(refEl2%faceNodes3, mold = refEl1%faceNodes3)
       refEl2%faceNodes3 = refEl1%faceNodes3
    ENDIF


  ENDSUBROUTINE deep_copy_refel_struct


END MODULE types
