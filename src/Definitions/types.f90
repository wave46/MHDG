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
    integer*4 :: ndim          ! Number of dimensions of the reference element
    integer*4 :: elemType      ! type of mesh elements (0 = triangles, 1 = quadrilaterals, 2 = thetrahedra, 3 = hexahedra)
    integer*4 :: Nvertices     ! Number of vertices of the element: 3 for triangles, 4 for quads and thetraedra, 8 for hexahedra
    integer*4 :: Ndeg          ! degree of interpolation
    integer*4 :: Nnodes3D      ! number of 3D nodes in the element
    integer*4 :: Nnodes2D      ! number of 2D nodes in the element
    integer*4 :: Nnodes1D      ! number of 1D nodes in the element
    integer*4 :: Nfaces        ! number of faces in the element: a face in a 2D mesh is a segment, in a 3D mesh is a 2D entity
    integer*4 :: Nbords        ! number of bords (only in 3D), is the number of 1D entities of the element
    integer*4 :: Nextnodes     ! number of exterior nodes of the element
    integer*4 :: Nfacenodes    ! number of nodes in each face of the element
    integer*4 :: Nfacenodeslin ! number of nodes in each linear face of the element
    integer*4 :: Nbordnodes    ! number of nodes in each bord of the element (only in 3D)
    integer*4 :: Ninnernodes   ! number of inner nodes of the element
    integer*4 :: Ninnernodesface  ! number of inner nodes of the 2D faces (only for 3D meshes)
    integer*4 :: NGauss1D      ! number of Gauss points of the 1D quadrature
    integer*4 :: NGauss2D      ! number of Gauss points of the 2D quadrature
    integer*4 :: NGauss3D      ! number of Gauss points of the 3D quadrature
    integer*4, allocatable :: Face_nodes(:, :)  ! numbering of the face nodes
    integer*4, allocatable :: Bord_nodes(:)    ! numbering of the 1D edge nodes (for 3D meshes)
    integer*4, allocatable :: inner_nodes(:)   ! numbering of the inner nodes
    integer*4, allocatable :: inner_nodes_face(:)   ! numbering of the inner nodes for 2D faces (for 3D meshes)
    real*8, pointer :: coord3d(:, :) => null()       ! spatial coordinates of the 3D nodes
    real*8, pointer :: coord2d(:, :) => null()       ! spatial coordinates of the 2D nodes
    real*8, pointer :: coord1d(:) => null()       ! spatial coordinates of the 1D nodes
    real*8, allocatable :: gauss_points3D(:, :) ! Gauss integration points for the 3D quadrature
    real*8, allocatable :: gauss_points2D(:, :) ! Gauss integration points for the 2D quadrature
    real*8, allocatable :: gauss_points1D(:)   ! Gauss integration points for the 1D quadrature
    real*8, allocatable :: gauss_weights3D(:)  ! Weights of the integration points for the 3D quadrature
    real*8, allocatable :: gauss_weights2D(:)  ! Weights of the integration points for the 2D quadrature
    real*8, allocatable :: gauss_weights1D(:)  ! Weights of the integration points for the 1D quadrature
    real*8, allocatable :: N3D(:, :)      ! Shape functions 3D at the Gauss points
    real*8, allocatable :: Nxi3D(:, :)    ! Derivative with respect to xi of the shape functions 3D at the Gauss points
    real*8, allocatable :: Neta3D(:, :)   ! Derivative with respect to eta of the shape functions 3D at the Gauss points
    real*8, allocatable :: Nzeta3D(:, :)  ! Derivative with respect to zeta of the shape functions 3D at the Gauss points
    real*8, allocatable :: N2D(:, :)      ! Shape functions 2D at the Gauss points
    real*8, allocatable :: Nxi2D(:, :)    ! Derivative with respect to xi of the shape functions 2D at the Gauss points
    real*8, allocatable :: Neta2D(:, :)   ! Derivative with respect to eta of the shape functions 2D at the Gauss points
    real*8, allocatable :: N1D(:, :)      ! Shape functions 1D at the Gauss points
    real*8, allocatable :: Nxi1D(:, :)    ! Derivative of the shape functions 1D at the Gauss points
    real*8, allocatable :: Nlin(:, :)     ! Linear shape functions 2D at the nodes of the element (only used in shock capturing so far...)
    real*8, allocatable :: sFTF(:, :)     ! Shape functions for the toroidal faces
    integer*4, allocatable :: faceNodes3(:, :) ! Face numbering for toroidal faces
    integer*4          :: Nfl           ! Number of nodes in lateral faces for toroidal 3d computations
    integer*4          :: Ngl           ! Number of Gauss points in lateral faces for toroidal 3d computations
    integer*4          :: Nft           ! Number of nodes in all faces for toroidal 3d computations
  END TYPE Reference_element_type

  !*******************************************************
  ! Mesh
  !*******************************************************
  TYPE :: Mesh_type
    integer*4 :: Ndim                ! Number of dimensions of the mesh
    integer*4 :: Nnodes              ! Number of nodes in the mesh
    integer*4 :: Nnodesperelem       ! Number of nodes per element in the mesh (dimension 2 of the connectvity matrix)
    integer*4 :: Nnodesperface       ! Number of nodes per face in the mesh (dimension 2 of boundary connectvity matrix )
    integer*4 :: Nelems              ! Number of elements in the mesh
    integer*4 :: Nfaces              ! Number of faces in the mesh
    integer*4 :: Nextfaces           ! Number of exterior faces
    integer*4 :: Nintfaces           ! Number of interior faces
    integer*4 :: elemType            ! 0 for quads - 1 for triangles - 2 for thetra - 3 for hexa
    integer*4 :: ndir                ! number of Dirichlet faces
    integer*4 :: ukf                 ! number of faces in the mesh minus the number of Dirichlet faces
    integer*4, pointer :: T(:, :) => null()    ! Elements connectivity matrix
    integer*4, pointer :: Tlin(:, :) => null()    ! Elements linear connectivity matrix
    integer*4, pointer :: Tb(:, :) => null()    ! Outer faces connectivity matrix
     INTEGER*4, POINTER :: boundaryFlag(:) => NULL()    ! Flag for the boundary condition for each external face (set in the mesh generator)
    integer*4, allocatable :: F(:, :)          ! Faces connectivity matrix
    integer*4, allocatable :: N(:, :)          ! Nodes connectivity matrix
     INTEGER*4, ALLOCATABLE :: face_info(:, :)          ! Elemental face info
    integer*4, allocatable :: faces(:, :, :)    ! for each triangle i, stores info k on each face j: faces(i,j,1) = # of neighbouring triangle (0 if external
    ! boundary), faces(i,j,2) = type of boundary (), faces(i,j,3) = type of boundary condition
    integer*4, allocatable :: extfaces(:, :)   ! for each exterior face, stores the number of the triangle, the number of the face, and the type of BC
    integer*4, allocatable :: intfaces(:, :)   ! for each interior face, stores the number of the triangle, the number of the face, the number of the
    ! neighboring triangle, the number of its face and the number of the node of the neighboring triangle that
    ! matches the first knot of the triangle
    logical, allocatable :: flipface(:, :)    ! for each triangle, and for each face, 0 if the order of the numbering in the face is to be kept, 1 if the
    ! order is to be reversed
    logical, allocatable    :: Fdir(:, :)         ! for each element states if each local face is of Dirichlet type
    integer*4, allocatable  :: periodic_faces(:)  ! Mapping for periodic faces
    integer*4, allocatable  :: Diric(:)
    integer*4, allocatable  :: numberbcs(:)
    real*8, allocatable     :: elemSize(:)   ! element size (area in 2D, volume in 3D) [Number of elements]
    real*8, pointer         :: X(:, :) => null()    ! nodes coordinates
    integer*4              :: Nnodes_toroidal     ! Number of nodes in the toroidal direction
    real*8, pointer         :: toroidal(:) => null()    ! nodes coordinates in the toroidal direction
    ! Limiting & shock capturing stuff
    integer, allocatable    :: flag_elems_rho(:)      ! Flagged elements for limiting rho [Number of elements]
    integer, allocatable    :: flag_elems_sc(:)     !  Flagged elements for shock-capturing [Number of elements]
    real*8, allocatable     :: minrho_elems(:)        ! Minimum value of density in the flagged elements[Number of elements]
    real*8, allocatable     :: sour_elems(:)          ! Source to limit rho in the flagged elements[Number of elements]
    real*8, allocatable     :: diff_elems(:)          ! Diffusion to limit rho in the flagged elements[Number of elements]
    real*8, allocatable     :: scdiff_nodes(:, :)      ! Shock capturing diffusion in each node [Number of elements,Number of nodes per element]
    real*8                 :: xmax, xmin, ymax, ymin    ! Limit of the GLOBAL matrix, across mpi partitions
    real*8                  :: puff_area         ! area of the puff bounday condition
    real*8                  :: pump_area         ! area of the pump bounday condition
    real*8                  :: core_area         ! area of the core bounday condition
    real*8,allocatable      :: Xg(:,:)               ! 2D Gauss point coordinates
    real*8,allocatable      :: Xgf(:,:)              ! 1D Gauss point coordinates at interior faces
    real*8,allocatable      :: Xgb(:,:)          ! 1D Gauss point  coordinates at boundary faces
     INTEGER*4, POINTER	    :: T_gmsh(:, :) => NULL()    ! Elements connectivity matrix to write the msh file
     INTEGER*4, POINTER      :: Tb_gmsh(:, :) => NULL()    ! Outer faces connectivity matrix to write the msh file
     REAL*8, POINTER         :: X_P1(:,:) => NULL()         ! coordinates of the nodes on the P1 mesh to write to the msh file
#ifdef PARALL
    integer, pointer       :: loc2glob_fa(:)     ! mapping number of the faces for creating the global matrix [number of faces in the mesh]
    integer, pointer       :: loc2glob_el(:)     ! mapping number of the elements from local to global [number of elements in the mesh]
     INTEGER, POINTER        :: loc2glob_nodes(:)    ! mapping number of the nodes from local to global [number of nodes in the mesh]
    integer, pointer       :: ghostfaces(:)      ! integer that states if a face is to be assembled locally or not [number of faces in the mesh]
    integer, pointer       :: ghostelems(:)      ! integer that states if an element is to be assembled locally or not (for 3D) [number of elements in the mesh]
    integer               :: nghostfaces        ! number of ghost faces
    integer               :: nghostelems        ! number of ghost elements (used only for 3D)
    integer               :: Nel_glob           ! number of elements of global mesh
    integer               :: Nfa_glob           ! number of faces of global mesh
    integer               :: Nno_glob           ! number of nodes of global mesh
    integer               :: Ndir_glob          ! number of Dirichlet faces in the global mesh
    integer               :: Ngho_glob          ! number of ghost faces in the global mesh
    ! readed from input
    integer, pointer       :: ghostflp(:)        ! flipFaces for the ghost faces [number of ghost faces]
    integer, pointer       :: ghostloc(:)        ! local numbering of the ghost face in the process that assemble it [number of ghost faces]
    integer, pointer       :: ghostpro(:)        ! the process that assemble the ghost face [number of ghost faces]
    integer, pointer       :: ghelsloc(:)        ! local numbering of the ghost element in the process that assemble it (only 3D) [number of ghost elements]
    integer, pointer       :: ghelspro(:)        ! the process that assemble the ghost element (only 3D) [number of ghost elements]
    ! built after reading from input
    integer, allocatable   :: fc2sd(:)          ! face 2 send: faces computed locally that the local process has to send (vector)
    integer, allocatable   :: pr2sd(:)          ! process 2 send: to which process the faces computed locally need to be sent (vector)
    integer, allocatable   :: fc2rv(:)          ! face 2 receive: ghost faces computed by other processes that the local process need to receive[number of ghost faces]
    integer, allocatable   :: pr2rv(:)          ! process 2 receive: from which process the faces computed externally need to be received [number of ghost faces] (it is the same as ghostpro)

    integer, allocatable   :: el2sd(:)          ! element 2 send: elements computed locally that the local process has to send (only 3D) (vector)
    integer, allocatable   :: pe2sd(:)          ! process 2 send: to which process the elements computed locally need to be sent (only 3D) (vector)
    integer, allocatable   :: el2rv(:)          ! element 2 receive: ghost elements computed by other processes that the local process need to receive (only 3D) [number of ghost elements]
    integer, allocatable   :: pe2rv(:)          ! process 2 receive: from which process the elements computed externally need to be received (only 3D) [number of ghost elements] (it is the same as ghelspro)

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
    integer         :: neq            ! Number of equations
    integer         :: npv            ! Number of physical variables
    real*8          :: diff_n, diff_u ! Perpendicular diffusion in the continuity and momentum equation
    real*8          :: v_p            ! Pinch velocity in the continuity equation
    real*8          :: a              ! Proportionality constant between pressure and density for isothermal model (p = a*rho)
    real*8          :: dfcoef         ! Constant related to the diamagnetic drift velocity
    real*8          :: dexbcoef         ! Constant related to the ExB drift velocity
    integer*4       :: bcflags(1:10)          ! Set the correspondence between mesh boundary flag (Mesh%boundaryFlag) and the boundary condition
    real*8          :: diagsource(1:10)       ! Diagonal implicit sources
    character(LEN=20), pointer:: phyVarNam(:) => Null() ! Names of the physical variables (set in initPhys)
    character(LEN=20), pointer:: conVarNam(:) => Null() ! Names of the conservative variables (set in initPhys)
    real*8          :: lscale          ! Length scale for the non-dimensionalization of the equations
    ! Magnetic field defined for each node of the mesh.
    real*8, pointer :: B(:, :)            ! Magnetic field, Br,Bz,Bphi [n of nodes  x 3]
    real*8          :: B0                 ! Reference value for the magnetic field [Tesla]
    real*8, pointer :: magnetic_flux(:)   ! Magnetic flux   [n of nodes]
    real*8, pointer :: magnetic_psi(:)    ! Magnetic flux normalized to separatrix magnetic flux   [n of nodes]
    real*8          :: Flux2Dmin          ! Minimum of the magnetic flux, across the MPI partitions
    real*8          :: Flux2Dmax          ! Maximum of the magnetic flux, across the MPI partitions
#ifdef KEQUATION
    real*8, pointer :: omega(:)           ! larmor frequency   [n of nodes]
    real*8, pointer :: q_cyl(:)           ! q cylindrical   [n of nodes]
#endif
    real*8          :: r_axis             ! R-coordinate of magnetic axis
    real*8          :: z_axis             ! Z-coordinate of magnetic axis

    real*8, pointer :: Bperturb(:, :)     ! Magnetic perturbation, Br,Bz,Bphi [n of nodes  x 3]
    real*8          :: Tbg                ! Background temperature in the isothermal model
    real*8, pointer :: Jtor(:)            ! Toroidal Current
    real*8          :: I_p                ! Total plasma current
    real*8          :: bohmth             ! Threshold for imposing the Bohm boundary condition
    ! Energy equation coefficients
    real*8          :: diff_e             ! Perpendicular diffusion in the energy equation
    real*8          :: epn                ! Exponential of the parallel diffusion (usually 5/2)
    real*8          :: Mref               ! Reference Mach number
    real*8          :: diff_pari          ! Parallel diffusion for the temperature (usually 1e7)
    real*8          :: Gmbohm             ! gamma for Bohm boundary condition on energy:
    ! Temperature equations coefficients (the ions coefficients are the ones defined previously)
    real*8          :: diff_ee            ! Perpendicular diffusion in the elcetron energy equation
    real*8          :: diff_pare          ! Parallel diffusion for the electron temperature
    real*8          :: tie                ! Temperature exchange coefficient between ions and electrons
    real*8          :: Gmbohme            ! gamma for Bohm boundary condition on electron energy:
    real*8          :: Pohmic             ! Ohmic heating power
    ! Coefficients for the vorticity equations
    real*8          :: diff_vort
    real*8          :: diff_pot
    real*8          :: etapar
    real*8          :: c1, c2             ! coefficients coming from the adimensionalization
    real*8          :: Potfloat
    ! Coefficients for ion heating
    real*8          :: heating_power      ! total power of ion heating source
    real*8          :: heating_amplitude  ! adimentionalized gaussian amplitude including integration coefficients
    real*8          :: heating_dr         ! displacement of the soruce from magnetic axis in r direction
    real*8          :: heating_dz         ! displacement of the soruce from magnetic axis in z direction
    real*8          :: heating_sigmar     ! width of the soruce from magnetic axis in r direction
    real*8          :: heating_sigmaz     ! width of the soruce from magnetic axis in z direction
    integer         :: heating_equation   ! Equation to which additional heating is applied (3 for ions, 4 for  electrons)
    ! Coefficients for the neutral equations
    real*8          :: diff_nn            ! Diffusion in the neutral equation
    real*8,allocatable:: diff_nn_Vol(:)   ! Diffusion in the neutral equation at 2D Gauss points
    real*8,allocatable:: diff_nn_Fac(:)   ! Diffusion in the neutral equation at 1D Gauss points on interior faces
    real*8,allocatable:: diff_nn_Bou(:)   ! Diffusion in the neutral equation at 1D Gauss points on boundary faces 
    real*8,allocatable:: v_nn_Vol(:,:)    ! Convective velocity in the neutral equation at 2D Gauss points
    real*8,allocatable:: v_nn_Fac(:,:)    ! Convective velocity in the neutral equation at 1D Gauss points on interior faces
    real*8,allocatable:: v_nn_Bou(:,:)    ! Convective velocity in the neutral equation at 1D Gauss points on boundary faces 
    real*8          :: Re                 ! Recycling for the neutral equation
    real*8          :: Re_pump            ! Recycling for the neutral equation in the pump region
    real*8          :: puff               ! Puff coefficient
    real*8          :: cryopump_power     ! Cryopump power in [m^3/s] coefficient
    real*8          :: puff_slope         ! Puff increment coefficient (only for moving equilibrium for ITER)
     real*8,pointer :: puff_exp(:)    ! Puff experimental coefficient (only for moving equilibriums)
    real*8          :: part_source        ! Particle source for ITER
    real*8          :: ener_source        ! Particle source for ITER
    real*8          :: density_source     ! Density source for WEST (2D, case 52)
    real*8          :: ener_source_e      ! Ion energy source for WEST and ITER (2D, case 52 and 81)
    real*8          :: ener_source_ee     ! Electron source for WEST and ITER (2D, case 52 and 81)
    real*8          :: sigma_source       ! Sigma for the gaussian sources for WEST and ITER (2D, case 52 and 81)
    real*8          :: fluxg_trunc        ! Value of the NORMALISED magnetic flux at which to truncate the gaussian sources for WEST (2D, case 52), refer to Source_shape.m file
    ! Diffusion coefficients ITER evolving equilibrium
    real*8          :: ME_diff_n
    real*8          :: ME_diff_u
    real*8          :: ME_diff_e
    real*8          :: ME_diff_ee
#ifdef EXPANDEDCX
#ifdef AMJUELCX
    real*8, dimension(9):: alpha_cx       ! Coefficients for charge exchange coefficients spline
#endif
#ifdef THERMALCX
    real*8, dimension(5):: alpha_cx       ! Coefficients for charge exchange coefficients spline
#endif
#endif
#ifdef AMJUELSPLINES
    ! Atomic rates coefficients
    real*8, dimension(9,9):: alpha_iz     ! Coefficients for ionization coefficients spline from EIRENE, (te,ne) grid
    real*8, dimension(9,9):: alpha_rec     ! Coefficients for recombination coefficients spline from EIRENE, (te,ne) grid
    real*8, dimension(9,9):: alpha_energy_iz     ! Coefficients for radiation losses due to ionization  spline from EIRENE, (te,ne) grid
    real*8, dimension(9,9):: alpha_energy_rec     ! Coefficients for radiation losses due to recombination  spline from EIRENE, (te,ne) grid
#endif
#ifdef KEQUATION
    ! Coefficients for the k equation
    real*8          :: diff_k_min         ! Mininmum diffusion in the k equation
    real*8          :: diff_k_max         ! Maximum diffusion in the k equation
    real*8          :: k_max              ! Maximum k
#endif
  END TYPE Physics_type

  !*******************************************************
  ! Geometry: type for geometry info
  !*******************************************************
  TYPE Geometry_type
    integer     :: ndim     ! Number of dimensions of the problem (can be different from the ndim of the mesh)
    real*8      :: R0       ! Major radius at the magnetic axis
    real*8      :: q        ! Safety factor
  END TYPE Geometry_type

  !*******************************************************
  ! Magnetic: type for 3D magnetic perturbation
  !*******************************************************
  TYPE Magnetic_type
    real*8          :: amp_rmp            ! amplitude RMP
    integer         :: nbCoils_rmp        ! number coils RMP (full torus)
    real*8          :: torElongCoils_rmp  ! Toroidal elongation of RMP coils (rectangular coils)
    integer         :: parite             ! parite RMP (for 2 row, -1 even parite, 1 odd)
    integer         :: nbRow              ! number of rows (1,2 or 3) of RMP coils, default is 2
    real*8          :: amp_ripple         ! amplitude Ripple
    integer         :: nbCoils_ripple     ! number coils Ripple (full torus)
    real*8          :: triang             ! triangularity (0: None)
    real*8          :: ellip              ! ellipticity (1: None)
    real*8, pointer :: coils_rmp(:, :, :) ! Coil coordinates for RMP (nbCoils*4*Discr,start-stop*(xyz)=6,rowNb) (4 for square coils)
    real*8, pointer :: coils_ripple(:, :) ! Coil coordinates for Ripple (nbCoils*Discr,start-stop*(xyz)=6)
  END TYPE Magnetic_type

  !*******************************************************
  ! Switches: type for main code switches
  !*******************************************************
  TYPE Switches_type
    logical :: axisym ! Is it an axisymmetric simulation?
    ! true =
    ! false =
    logical :: rmp      ! To activate resonant magnetic perturbation
    logical :: ripple   ! To activate ripple
    logical :: ohmicsrc ! Set to TRUE to consider ohmic source of energy
    logical :: ME       ! Set to TRUE to allow magnetic equilibrium evolution in time
    logical :: driftdia ! Set to TRUE to consider diamagnetic drift
    logical :: driftexb ! Set to TRUE to consider ExB drift
    logical :: steady
     LOGICAL :: read_gmsh ! read a mesh file in gmsh format
     LOGICAL :: readMeshFromSol ! read a mesh file from the solution
     LOGICAL :: set_2d_order ! if read_gmsh = .true., set order
     INTEGER :: order_2d    ! if read_gmsh = .true., set order = .true., what 2d polynomial order?
     LOGICAL :: gmsh2h5
     LOGICAL :: saveMeshSol
    logical :: time_init ! true if it is a time initialization simulation. The time counter "it" does not increment  (i.e. when the analitical initialisation is not good enough). Used for moving equilibrium (case 59)
    integer :: init     ! 1-init. analy. solution at nodes; 2-L2 projection
    ! Set to TRUE for a steady state computation
    ! Set to FALSE for a transient computation
    integer :: testcase  ! Define the testcase ( which analytical solution, body force, magnetic field...)
    logical :: psdtime ! Reduce the diffusion every time we reach the steady state
    ! condition (only works if steady=.false.)
    real*8  :: diffred ! Reduction factor of the diffusion for psdtime simulation
    real*8  :: diffmin ! Minimum value of the diffusion for a psdtime simulation
    integer :: shockcp ! Shock capturing option
    integer :: limrho  ! Add a source for limiting the min value of rho
    integer :: difcor  ! Add diffusion in corners
    integer :: thresh  ! Use a threshold for limiting the min value of rho
    ! (rho-rho*E in case of N-Gamma-Energy model, rho-rho*Ei-rho*Ee in case of N-Gamma-Ti-Te model)
    logical :: filter  ! Filter solution to avoid oscillation in empty zones
    logical :: decoup  ! Decouple N-Gamma from Te-Ti (only used for N-Gamma-Ti-Te model)
    logical :: ckeramp ! Chech the error amplification in the linear system solution (for very ill-conditioned matrices)
    logical :: saveNR  ! Save solution at each NR iteration
    logical :: saveTau ! Save tau on faces
    logical :: fixdPotLim
    logical :: dirivortcore
    logical :: dirivortlim
    logical :: convvort ! consider the convective term in the vorticity equation
    logical :: bxgradb  ! consider the term in BxGradB in the vorticity equation
    integer :: pertini  ! add perturbation in the initial solution
    ! 1 -add sinusoidal perturbation
    ! 2 -add density blob
    logical :: logrho   ! solve for the density logarithm instead of density
  END TYPE Switches_type

  !***************************************************************
  ! Paths: type for storing paths to load inputs and sotre outputs
  !***************************************************************
  TYPE Inputs_type
    character(len=1000) :: field_path ! where do we read magnetic field from (WEST cases so far)
    character(len=1000) :: jtor_path ! where do we read plasma current from (WEST cases so far)
    character(len=1000) :: save_folder ! where to save last solution
    logical             :: field_from_grid !if true, then reads equilibrium file n rectangular grid; if false - on nodes of the mesh
    logical             :: compute_from_flux ! if components B_R, B_Z are computed from flux or not
    logical             :: divide_by_2pi     ! correspondng to flux definition if it is needed to divide by 2pi or not
    integer             :: field_dimensions(1:2) ! dimensions of magnetic field files (2D WEST cases so far)
    integer             :: jtor_dimensions(1:2) ! dimensions of magnetic field files (2D WEST cases so far)
  END TYPE Inputs_type

  !*******************************************************
  ! Time: type for the time stepping information
  !*******************************************************
  TYPE Time_type
    real*8      :: dt0  ! initial time step
    real*8      :: dt   ! current time step
    real*8      :: dt_ME   ! time step from equilibrium !stored as dimensional value
    real*8      :: t_ME   ! time step from equilibrium (to get puff rate) !stored as dimensional value
    real*8      :: tfi  ! final time of the simulation
    integer     :: it   ! the number of the current time step
    integer     :: ik   ! same as it but always incrementing (also in case of pseudotime..)
    integer     :: ndt  ! max number of time steps to do in the current session
    integer     :: tsw  ! switch to modify the time step
    integer     :: nts  ! max number of time iterations to do in the current session (only for transient simulations)
    integer     :: tis  ! time integration scheme
    ! 1 - first order
    ! 2 - second order
    real*8      :: t    ! time of the simulation (initialized to finish time of previous simulation if restart, to 0 if new simulation)
  END TYPE Time_type

  !*******************************************************
  ! Numerics: type for numeric scheme parameters
  !*******************************************************
  TYPE Numeric_type
    integer        :: nrp      ! Max number of Newton-Raphson iterations
    real*8         :: tNR      ! Tolerance of the Newton-Raphson scheme
    real*8         :: tTM      ! Tolerance for the steady state achievement
    real*8         :: div      ! Divergence detector
    real*8         :: tau(1:5) ! Stabilization parameter for each equation (4 values max for now...)
    real*8         :: sc_coe   ! Shock capturing coefficient
    real*8         :: sc_sen   ! Shock capturing sensibility
    real*8         :: minrho   ! Value of rho to start applying limiting
    real*8         :: so_coe   ! Source coefficient for limiting rho
    real*8         :: df_coe   ! Diffusion coefficient for limiting rho
    real*8         :: dc_coe   ! Diffusion coefficient in corners
    real*8         :: thr      ! Threshold to limit rho
    real*8         :: thrpre   ! Threshold to limit pressure
    integer        :: stab     ! Stabilization type
    ! 1 - constant tau (one for each equation) in the whole domain
    ! 2 -
    ! 3 -
    real*8         :: dumpnr   ! dumping factor for Newton-Raphson. 0<dumpnr<1
    real*8         :: dumpnr_min   ! dumping factor minimum for Newton-Raphson. 0<dumpnr<1
    real*8         :: dumpnr_max   ! dumping factor maximum for Newton-Raphson. 0<dumpnr<1
    real*8         :: dumpnr_width   ! dumping factor width of hyperbolic tangential for Newton-Raphson. 0<dumpnr<1
    real*8         :: dumpnr_n0   ! dumping factor x0 for hyperbolic tangential for Newton-Raphson. 0<dumpnr<1
    integer        :: ntor     ! Number of elements in the toroidal direction
    integer        :: ptor     ! Polynomial degree in the toroidal direction
    real*8         :: tmax     ! Max extention in the toroidal direction
    integer        :: npartor  ! Number of MPI divisions in the toroidal direction
    integer        :: bohmtypebc ! Implementation of the Bohm bc for Gamma
    real*8         :: exbdump ! Dumping for ExB drifts
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
    integer :: printint       ! Integer for printing
    logical :: timing         ! Timing of the code
    integer :: freqdisp       ! Frequency of results display
    integer :: freqsave       ! Frequency of solution save
  END TYPE Utils_type

  !*******************************************************
  ! Linear system solver parameters
  !*******************************************************
  TYPE Lssolver_type
    integer           :: sollib    ! Solver library to be used
    ! 1-Pastix
    ! 2-PSBLAS
     ! 3-PETSc

    logical           :: timing    ! timing of the linear solver
    ! Parameters relative to the library PSBLAS
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
    character(len=20) :: kmethd    ! Krylov method (see list on the library manual)
    integer           :: istop     ! Stopping criterion (see spec on the library manual)
    integer           :: itmax     ! Max number of iterations
    integer           :: itrace    ! Display convergence at each iteration
    integer           :: rest      ! Restart
    real*8            :: tol       ! Stopping tolerance
    character(len=20) :: ptype     ! Preconditioner type
    ! Parameters relative to the library MLD2P4
    ! First smoother / 1-lev preconditioner
    character(len=20) :: smther       ! smoother type
    integer           :: jsweeps      ! (pre-)smoother / 1-lev prec sweeps
    integer           :: novr         ! Number of olverlap layers, for Additive Schwartz only
    character(len=20) :: restr        ! restriction  over application for Additive Schwartz only
    character(len=20) :: prol         ! Type of prolongation operator for Additive Schwartz only
    character(len=20) :: solve        ! local subsolver
    integer           :: fill         ! fill-in level p of the ILU factorizations
    real              :: thr          ! threshold for ILUT
    ! Second smoother/ AMG post-smoother (if NONE ignored in main)
    character(len=20) :: smther2      ! smoother type
    integer           :: jsweeps2     ! (post-)smoother sweeps
    integer           :: novr2        ! number of overlap layers
    character(len=20) :: restr2       ! restriction  over application of AS
    character(len=20) :: prol2        ! prolongation over application of AS
    character(len=20) :: solve2       ! local subsolver
    integer           :: fill2        ! fill-in for incomplete LU
    real              :: thr2         ! threshold for ILUT
    ! general AMG data
    character(len=20) :: mlcycle      ! multi-level cycle
    integer           :: outer_sweeps ! number of multilevel cycles
    integer           :: maxlevs      ! Maximum number of levels
    integer           :: csize        ! Coarse size threshold
    ! aggregation
    real              :: mncrratio    ! Minimum coarsening ratio
    real              :: athres       ! Aggregation threshold
    character(len=20) :: aggr_prol    ! Prolongator used by the aggregation algorithm
    character(len=20) :: par_aggr_alg ! Parallel aggregation algorithm
    character(len=20) :: aggr_ord     ! Initial ordering of indices for the aggregation algorithm
    character(len=20) :: aggr_filter  ! Matrix used in computing the smoothed prolongator
    ! coasest-level solver
    character(len=20) :: csolve       ! coarsest-lev solver
    character(len=20) :: csbsolve     ! coarsest-lev solver
    character(len=20) :: cmat         ! coarsest mat layout
    integer           :: cfill        ! fill-in for incompl LU
    real              :: cthres       ! Threshold for ILUT
    integer           :: cjswp        ! sweeps for GS/JAC subsolver
  END TYPE Lssolver_type

  !**********************************************************
  ! Solution: contains the solution at the current time step
  !**********************************************************
  TYPE Sol_type
    real*8, pointer :: u(:) => null()  ! Elemental solution
     REAL*8, POINTER :: u_conv(:) => NULL()
     REAL*8, POINTER :: u_init(:) => NULL()
     REAL*8, POINTER :: u0_init(:) => NULL()
    real*8, pointer :: u_tilde(:) => null()  ! Face solution
     REAL*8, POINTER :: u_tilde0(:) => NULL()  ! Face solution
    real*8, pointer :: q(:) => null()  ! Elemental solution for the gradient
     REAL*8, POINTER :: q_conv(:) => NULL()
     REAL*8, POINTER :: q_init(:) => NULL()
     REAL*8, POINTER :: q0_init(:) => NULL()
    real*8, allocatable :: u0(:, :)           ! Elemental solution at previous time steps
    real*8, allocatable :: tres(:)           ! Time residual
    real*8, allocatable :: time(:)           ! Time evolution
    integer             :: Nt                ! Number of time-steps
  END TYPE Sol_type

  !**********************************************************
  ! Simulation parameters: for saving purpose
  !**********************************************************
  TYPE Simulationparams_type
    character(len=50) :: model
    integer   :: Ndim
    integer   :: Neq
    real, allocatable :: consvar_refval(:)
    real, allocatable :: physvar_refval(:)
    ! Reference values
    real*8    :: refval_length
    real*8    :: refval_mass
    real*8    :: refval_charge
    real*8    :: refval_time
    real*8    :: refval_temperature
    real*8    :: refval_density
    real*8    :: refval_neutral
#ifdef KEQUATION
    real*8    :: refval_k
#endif
    real*8    :: refval_speed
    real*8    :: refval_potential
    real*8    :: refval_vorticity
    real*8    :: refval_magfield
    real*8    :: refval_current
    real*8    :: refval_diffusion
    real*8    :: refval_momentum
    real*8    :: refval_specpress
    real*8    :: refval_specenergy
    real*8    :: refval_specenergydens
    ! Adimesional isothermal compressibility coefficient
    real*8    :: compress_coeff

    ! Dimensions used in the simulation
    character(len=20)    :: refval_length_dimensions
    character(len=20)    :: refval_mass_dimensions
    character(len=20)    :: refval_charge_dimensions
    character(len=20)    :: refval_time_dimensions
    character(len=20)    :: refval_temperature_dimensions
    character(len=20)    :: refval_density_dimensions
    character(len=20)    :: refval_neutral_dimensions
#ifdef KEQUATION
    character(len=20)    :: refval_k_dimensions
#endif
    character(len=20)    :: refval_speed_dimensions
    character(len=20)    :: refval_potential_dimensions
    character(len=20)    :: refval_vorticity_dimensions
    character(len=20)    :: refval_magfield_dimensions
    character(len=20)    :: refval_current_dimensions
    character(len=20)    :: refval_diffusion_dimensions
    character(len=20)    :: refval_momentum_dimensions
    character(len=20)    :: refval_specpress_dimensions
    character(len=20)    :: refval_specenergy_dimensions
    character(len=20)    :: refval_specenergydens_dimensions

    ! Physical parameters used in the computation
    real*8    :: a
    real*8    :: Mref
    real*8    :: c1
    real*8    :: c2
    real*8    :: diff_pari
    real*8    :: diff_pare
    real*8    :: diff_n
    real*8    :: diff_u
    real*8    :: diff_e
    real*8    :: diff_ee
  END TYPE Simulationparams_type

  !**********************************************************
  ! Elemental matrices: type to store the elemental matrices
  ! used during the computation
  !**********************************************************
  TYPE :: elmat_type
    real*8, allocatable :: iAqq(:, :, :)
    real*8, allocatable :: Aqu(:, :, :)
    real*8, allocatable :: Aql(:, :, :)
    real*8, allocatable :: Auq(:, :, :)
    real*8, allocatable :: Auu(:, :, :)
    real*8, allocatable :: Aul(:, :, :)
    real*8, allocatable :: Alq(:, :, :)
    real*8, allocatable :: Alu(:, :, :)
    real*8, allocatable :: All(:, :, :)
    real*8, allocatable :: Aql_dir(:, :)
    real*8, allocatable :: Aul_dir(:, :)

    real*8, allocatable :: M(:, :, :)
    real*8, allocatable :: Cv(:, :, :)
    real*8, allocatable :: H(:, :, :)
    real*8, allocatable :: Hdir(:, :)
    real*8, allocatable :: D(:, :, :)
    real*8, allocatable :: E(:, :, :)
    real*8, allocatable :: Edir(:, :)
    real*8, allocatable :: S(:, :)
    real*8, allocatable :: UU(:, :, :)
    real*8, allocatable :: U0(:, :)
    real*8, allocatable :: Hf(:, :, :)
    real*8, allocatable :: Df(:, :, :)
    real*8, allocatable :: Ef(:, :, :)
    real*8, allocatable :: fH(:, :)
    real*8, allocatable :: B(:, :, :)
    real*8, allocatable :: C(:, :, :)
    real*8, allocatable :: Cdir(:, :)
    real*8, allocatable :: P(:, :, :)
    real*8, allocatable :: G(:, :, :)
    real*8, allocatable :: IL(:, :, :)
    real*8, allocatable :: Lf(:, :, :)
    real*8, allocatable :: Qf(:, :, :)
    real*8, allocatable :: LL(:, :, :)
    real*8, allocatable :: L0(:, :)
    ! Limiting rho
    real*8, allocatable :: S_lrho(:, :)
    real*8, allocatable :: P_lrho(:, :, :)
    ! Shock capturing
    real*8, allocatable :: P_sc(:, :, :)
    real*8, allocatable :: Lf_sc(:, :, :)
    real*8, allocatable :: TQ(:, :, :)
    real*8, allocatable :: TQhf(:, :, :)
    real*8, allocatable :: Tfhf(:, :)
    real*8, allocatable :: Tf(:, :)
    real*8, allocatable :: Thdir(:, :)

  END TYPE elmat_type

  !**********************************************************
  ! Timing: structure used to store execution times
  !**********************************************************
  TYPE Timing_type
    integer :: cks1, clock_rate1, cke1, clock_start1, clock_end1
    integer :: cks2, clock_rate2, cke2, clock_start2, clock_end2
    real*8  :: time_start1, time_finish1, tps1, tpe1
    real*8  :: time_start2, time_finish2, tps2, tpe2
     REAL*8  :: cputpre, cputmap, cputass, cputbcd, cputsol, cputjac, cputglb, cputadapt
     REAL*8  :: runtpre, runtmap, runtass, runtbcd, runtsol, runtjac, runtglb, runtadapt
    real(8) :: clstime1, clstime2, clstime3, clstime4, clstime5, clstime6
    real(8) :: rlstime1, rlstime2, rlstime3, rlstime4, rlstime5, rlstime6
    real(8) :: cputcom
    real(8) :: runtcom
  END TYPE Timing_type


CONTAINS

  SUBROUTINE set_boundary_flag_names()
    bc_flag_type(1) = 'Tb_Dirichlet'
    bc_flag_type(2) = 'Tb_LEFT'
    bc_flag_type(3) = 'Tb_RIGHT'
    bc_flag_type(4) = 'Tb_UP'
    bc_flag_type(5) = 'Tb_DOWN'
    bc_flag_type(6) = 'Tb_WALL'
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
