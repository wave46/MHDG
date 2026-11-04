!*****************************************
! project: MHDG
! file: Main_utils.f90
! date: 25/11/2024
! Set of functions and subroutine used in MHDG.f90

MODULE Main_utils

   USE in_out
   USE GMSH_io_module, ONLY: load_gmsh_mesh, gmsh_mesh2d_write
   USE reference_element
   USE preprocess
   USE MPI_OMP
   USE printutils
   USE debug
   USE initialization
#ifdef WITH_PASTIX
   USE solve_pastix, ONLY: terminate_mat_PASTIX
#endif
#ifdef WITH_PSBLAS
   USE solve_psblas, ONLY: terminate_PSBLAS
#endif
#ifdef WITH_PETSC
   USE solve_petsc, ONLY: matPETSC, InitPETSC, terminate_PETSC, FinalizePETSC
#endif
   USE adaptivity_general_module
   USE adaptivity_projection_module, ONLY: projectSolutionDifferentMeshes_general, projectSolutionDifferentMeshes_general_arrays
   USE Postprocess, ONLY: computeL2ErrorAnalyticSol
#ifdef PARALL
   USE Communications
   USE domain_decomposition_module
#else
   USE domain_decomposition_module, ONLY: reorder_mesh, reorder_solution
#endif
   USE HDG_LimitingTechniques

   IMPLICIT NONE

   INTEGER                      :: Np, Nel, Nfp, Nf, Ndim, Nthreads
   INTEGER                       :: it, ir, ir_check, it0, nts, nu, nut, nb_args, IERR, k, i,is, ir_count_adapt = 0, count_adapt = 0, ir_adapt = 0, checkpoint = 1, divergence_counter_adapt = 0, convergence_counter = 1, order

   LOGICAL, ALLOCATABLE         :: mkelms(:)
   REAL*8                       :: dt, dt0, errNR, errNR_adapt, errlstime
   REAL*8, ALLOCATABLE          :: uiter(:), L2err(:), u0_temp(:, :)
   REAL*8, POINTER              :: uiter_best(:) => NULL(), qiter_best(:) => NULL()
   CHARACTER(LEN=1024)          :: mesh_name, mesh_name_proj, save_name
   CHARACTER(len=255)           :: gmsh_filename
   REAL*8                       :: cputtot, runttot, time_start, time_end
   INTEGER                      :: OMP_GET_MAX_THREADS
   INTEGER                      :: clock_rate, clock_start, clock_end
   REAL*8, ALLOCATABLE           :: xs(:, :)
   REAL*8, ALLOCATABLE          :: oscillations(:)
   REAL*8                       :: max_osc = -100., min_osc = -100.
   INTEGER                      :: n_osc
   TYPE(Mesh_type)              :: Mesh_prec
   LOGICAL                      :: restart_adapt, time_adapt, reorder = .TRUE.
   INTEGER, ALLOCATABLE         :: perm_elems(:), perm_faces(:)
#ifdef PARALL
   REAL*8, POINTER              :: X_glob(:, :) => NULL(), B_glob(:, :) => NULL()
   REAL*8, POINTER              :: u_glob(:) => NULL(), q_glob(:) => NULL(), magnetic_flux_glob(:) => NULL(), u_tilde_glob(:) => NULL()
#endif
#ifdef TOR3D
   INTEGER                      :: Neq, ntorloc, Nfl, Np1dPol, Np1dTor, Ng1dPol, Ng1dTor
#endif

CONTAINS

#ifndef TOR3D
#ifndef PARALL
   SUBROUTINE adaptivity_serial_project_solution_2D()

      Mesh%X = Mesh%X*phys%lscale
      Mesh%X = Mesh%X/phys%lscale

      ! Project the check-point solution to the new mesh
      CALL projectSolutionDifferentMeshes_general(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X, sol%u_conv, sol%q_conv, sol%u, sol%q)

   END SUBROUTINE adaptivity_serial_project_solution_2D
#else
   SUBROUTINE adaptivity_parallel_project_solution_2D()
      INTEGER, POINTER          :: T_glob(:,:)
      REAL*8, POINTER           :: X_glob(:,:)

      NULLIFY(T_glob, X_glob)

      Mesh%X = Mesh%X/phys%lscale

      ! Domain decomposition on the new mesh
      CALL domain_decomposition(mesh_name)

      Mesh%X = Mesh%X*phys%lscale

      ! Restart communication
      CALL init_com()

      ! Gather checkpoint solution and project it to the new mesh
      CALL gather_mesh(Mesh_prec, T_glob=T_glob, X_glob=X_glob)
      CALL gather_solution(Mesh_in = Mesh_prec, Nnodesperelem = Mesh_prec%Nnodesperelem, u_in = sol%u_conv, q_in = sol%q_conv, u_glob = u_glob, q_glob = q_glob)
      CALL projectSolutionDifferentMeshes_general(T_glob, X_glob, Mesh%T, Mesh%X, u_glob, q_glob, sol%u, sol%q)

      DEALLOCATE (X_glob, T_glob, u_glob, q_glob)
      NULLIFY (X_glob, T_glob, u_glob, q_glob)

      ! Gather projected solution on decomposed mesh
      CALL gather_mesh(Mesh_in=Mesh, T_glob=T_glob, X_glob=X_glob)
      CALL gather_solution(Mesh_in=Mesh, Nnodesperelem=Mesh%Nnodesperelem, u_in=sol%u, q_in=sol%q, u_glob=u_glob, q_glob=q_glob)

      DEALLOCATE (X_glob, T_glob, u_glob, q_glob)
      NULLIFY (X_glob, T_glob, u_glob, q_glob)

   END SUBROUTINE adaptivity_parallel_project_solution_2D
#endif
#else
#ifndef PARALL
   SUBROUTINE adaptivity_serial_project_solution_3D()
      INTEGER                   :: j, itor, ptor, iel, Ntor, Nelems, Nelems_prec, Nnodesperelem, Nnodesperelem_tor1d, Neq, start_range_u, end_range_u, start_range_q, end_range_q
      INTEGER, ALLOCATABLE      :: u_tor_slice_indices_elem(:), q_tor_slice_indices_elem(:), u_tor_slice_indices_prec(:), q_tor_slice_indices_prec(:), u_tor_slice_indices(:), q_tor_slice_indices(:)
      REAL*8, POINTER           :: u_tor_slice(:), q_tor_slice(:), u_conv_tor_slice(:), q_conv_tor_slice(:)

      NULLIFY (u_tor_slice, q_tor_slice, u_conv_tor_slice, q_conv_tor_slice)
      DEALLOCATE (sol%u, sol%q)
      NULLIFY (sol%u, sol%q)

      Ntor = numer%ntor
      Nelems_prec = Mesh_prec%Nelems
      Nelems = Mesh%Nelems
      Nnodesperelem = Mesh%Nnodesperelem
      Nnodesperelem_tor1d = numer%ptor + 1
      Neq = phys%Neq
      Ndim = mesh%Ndim + 1

      ALLOCATE (sol%u(Nelems*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (sol%q(Nelems*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*Ndim))

      sol%u = 0.
      sol%q = 0.

      Mesh%X = Mesh%X*phys%lscale
      Mesh%X = Mesh%X/phys%lscale

      ! indices corresponding to one element DoFs
      ALLOCATE (u_tor_slice_indices_elem(Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice_indices_elem(Nnodesperelem*Neq*Ndim))

      ! indices corresponding to the whole toroidal slice for the converged solution on the previous mesh
      ALLOCATE (u_tor_slice_indices_prec(Nelems_prec*Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice_indices_prec(Nelems_prec*Nnodesperelem*Neq*Ndim))

      ! solution for one toroidal slice on the previous mesh
      ALLOCATE (u_conv_tor_slice(Nelems_prec*Nnodesperelem*Neq))
      ALLOCATE (q_conv_tor_slice(Nelems_prec*Nnodesperelem*Neq*Ndim))

      ! indices corresponding to the whole toroidal slice for the new mesh
      ALLOCATE (u_tor_slice_indices(Nelems*Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice_indices(Nelems*Nnodesperelem*Neq*Ndim))

      ! solution for one toroidal slice on the new mesh
      ALLOCATE (u_tor_slice(Nelems*Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice(Nelems*Nnodesperelem*Neq*Ndim))

      DO itor = 1, Ntor
         DO ptor = 1, Nnodesperelem_tor1d

            start_range_u = 1
            end_range_u = 1
            start_range_q = 1
            end_range_q = 1

            ! Extract the solution corresponding to the current toroidal slice
            DO iel = 1, Nelems_prec
               u_tor_slice_indices_elem = (itor-1)*Nelems_prec*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]
               q_tor_slice_indices_elem = (itor-1)*Nelems_prec*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim + (ptor-1)*Neq*Nnodesperelem*Ndim + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d*Ndim + [(i, i = 1, Neq*Nnodesperelem*Ndim)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1
               end_range_q = start_range_q + Neq*Nnodesperelem*Ndim - 1

               u_tor_slice_indices_prec(start_range_u:end_range_u) = u_tor_slice_indices_elem
               q_tor_slice_indices_prec(start_range_q:end_range_q) = q_tor_slice_indices_elem

               start_range_u = end_range_u + 1
               start_range_q = end_range_q + 1
            END DO

            u_conv_tor_slice = sol%u_conv(u_tor_slice_indices_prec)
            q_conv_tor_slice = sol%q_conv(q_tor_slice_indices_prec)

            ! Project the check-point solution to the new mesh
            CALL projectSolutionDifferentMeshes_general(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X, u_conv_tor_slice, q_conv_tor_slice, u_tor_slice, q_tor_slice)

            start_range_u = 1
            end_range_u = 1
            start_range_q = 1
            end_range_q = 1

            ! Fill in the toroidal slice projected solution in the correct location of the solution vector on the new mesh
            DO iel = 1, Nelems
               u_tor_slice_indices_elem = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]
               q_tor_slice_indices_elem = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim + (ptor-1)*Neq*Nnodesperelem*Ndim + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d*Ndim + [(i, i = 1, Neq*Nnodesperelem*Ndim)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1
               end_range_q = start_range_q + Neq*Nnodesperelem*Ndim - 1

               u_tor_slice_indices(start_range_u:end_range_u) = u_tor_slice_indices_elem
               q_tor_slice_indices(start_range_q:end_range_q) = q_tor_slice_indices_elem

               start_range_u = end_range_u + 1
               start_range_q = end_range_q + 1
            END DO

            sol%u(u_tor_slice_indices) = u_tor_slice
            sol%q(q_tor_slice_indices) = q_tor_slice

         END DO
      END DO

      DEALLOCATE (u_tor_slice, q_tor_slice, u_conv_tor_slice, q_conv_tor_slice, u_tor_slice_indices_elem, q_tor_slice_indices_elem, u_tor_slice_indices_prec, q_tor_slice_indices_prec, u_tor_slice_indices, q_tor_slice_indices)
      NULLIFY (u_tor_slice, q_tor_slice, u_conv_tor_slice, q_conv_tor_slice)

   END SUBROUTINE adaptivity_serial_project_solution_3D
#else
   SUBROUTINE adaptivity_parallel_project_solution_3D()
      INTEGER                   :: j, itor, ptor, iel, Ntor, Nelems, Nelems_prec,Nelems_glob, Nelems_prec_glob, Nnodesperelem, Nnodesperelem_tor1d, Neq, start_range_u, end_range_u, start_range_q, end_range_q
      INTEGER, ALLOCATABLE      :: u_tor_slice_indices_elem(:), q_tor_slice_indices_elem(:), u_tor_slice_indices_prec(:), q_tor_slice_indices_prec(:), u_tor_slice_indices(:), q_tor_slice_indices(:)
      REAL*8, POINTER           :: u_tor_slice(:), q_tor_slice(:), u_conv_tor_slice(:), q_conv_tor_slice(:)
      REAL*8, POINTER           :: X_glob(:, :)
      INTEGER, POINTER           :: T_glob(:, :)

      Mesh%X = Mesh%X/phys%lscale

      NULLIFY (u_tor_slice, q_tor_slice, u_conv_tor_slice, q_conv_tor_slice, X_glob, T_glob)
      DEALLOCATE (sol%u, sol%q)
      NULLIFY (sol%u, sol%q)

      CALL domain_decomposition(mesh_name)

      Mesh%X = Mesh%X*phys%lscale

      ! Restart communication
      CALL init_com()

      Ntor = numer%ntor
      Nelems_prec_glob = Mesh_prec%Nel_glob
      Nelems_prec = Mesh_prec%Nelems
      Nelems_glob = Mesh%Nel_glob
      Nelems = Mesh%Nelems
      Nnodesperelem = Mesh%Nnodesperelem
      Nnodesperelem_tor1d = numer%ptor + 1
      Neq = phys%Neq
      Ndim = mesh%Ndim + 1

      ALLOCATE (sol%u(Nelems*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq))
      ALLOCATE (sol%q(Nelems*Ntor*Nnodesperelem_tor1d*Nnodesperelem*Neq*Ndim))

      sol%u = 0.
      sol%q = 0.

      Mesh%X = Mesh%X*phys%lscale
      Mesh%X = Mesh%X/phys%lscale

      ! indices corresponding to one element DoFs
      ALLOCATE (u_tor_slice_indices_elem(Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice_indices_elem(Nnodesperelem*Neq*Ndim))

      ! indices corresponding to the whole toroidal slice for the converged solution on the previous mesh
      ALLOCATE (u_tor_slice_indices_prec(Nelems_prec_glob*Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice_indices_prec(Nelems_prec_glob*Nnodesperelem*Neq*Ndim))

      ! solution for one toroidal slice on the previous mesh
      ALLOCATE (u_conv_tor_slice(Nelems_prec_glob*Nnodesperelem*Neq))
      ALLOCATE (q_conv_tor_slice(Nelems_prec_glob*Nnodesperelem*Neq*Ndim))
      ! indices corresponding to the whole toroidal slice for the new mesh
      ALLOCATE (u_tor_slice_indices(Nelems*Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice_indices(Nelems*Nnodesperelem*Neq*Ndim))

      ! solution for one toroidal slice on the new mesh
      ALLOCATE (u_tor_slice(Nelems*Nnodesperelem*Neq))
      ALLOCATE (q_tor_slice(Nelems*Nnodesperelem*Neq*Ndim))

      ! Gather checkpoint solution and project it to the new mesh
      CALL gather_mesh(Mesh_prec, T_glob=T_glob, X_glob=X_glob)
      CALL gather_solution_3D(Mesh_in=Mesh_prec, u_in=sol%u_conv, q_in=sol%q_conv, u_glob=u_glob, q_glob=q_glob)

      DO itor = 1, Ntor
         DO ptor = 1, Nnodesperelem_tor1d

            start_range_u = 1
            end_range_u = 1
            start_range_q = 1
            end_range_q = 1

            ! Extract the solution corresponding to the current toroidal slice
            DO iel = 1, Nelems_prec_glob
               u_tor_slice_indices_elem = (itor-1)*Nelems_prec_glob*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]
               q_tor_slice_indices_elem = (itor-1)*Nelems_prec_glob*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim + (ptor-1)*Neq*Nnodesperelem*Ndim + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d*Ndim + [(i, i = 1, Neq*Nnodesperelem*Ndim)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1
               end_range_q = start_range_q + Neq*Nnodesperelem*Ndim - 1

               u_tor_slice_indices_prec(start_range_u:end_range_u) = u_tor_slice_indices_elem
               q_tor_slice_indices_prec(start_range_q:end_range_q) = q_tor_slice_indices_elem

               start_range_u = end_range_u + 1
               start_range_q = end_range_q + 1
            END DO

            u_conv_tor_slice = u_glob(u_tor_slice_indices_prec)
            q_conv_tor_slice = q_glob(q_tor_slice_indices_prec)

            CALL projectSolutionDifferentMeshes_general(T_glob, X_glob, Mesh%T, Mesh%X, u_conv_tor_slice, q_conv_tor_slice, u_tor_slice, q_tor_slice)

            start_range_u = 1
            end_range_u = 1
            start_range_q = 1
            end_range_q = 1

            ! Fill in the toroidal slice projected solution in the correct location of the solution vector on the new mesh
            DO iel = 1, Nelems
               u_tor_slice_indices_elem = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]
               q_tor_slice_indices_elem = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq*Ndim + (ptor-1)*Neq*Nnodesperelem*Ndim + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d*Ndim + [(i, i = 1, Neq*Nnodesperelem*Ndim)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1
               end_range_q = start_range_q + Neq*Nnodesperelem*Ndim - 1

               u_tor_slice_indices(start_range_u:end_range_u) = u_tor_slice_indices_elem
               q_tor_slice_indices(start_range_q:end_range_q) = q_tor_slice_indices_elem

               start_range_u = end_range_u + 1
               start_range_q = end_range_q + 1
            END DO

            sol%u(u_tor_slice_indices) = u_tor_slice
            sol%q(q_tor_slice_indices) = q_tor_slice

         END DO
      END DO

      DEALLOCATE (u_tor_slice, q_tor_slice, u_conv_tor_slice, q_conv_tor_slice, u_tor_slice_indices_elem, q_tor_slice_indices_elem, u_tor_slice_indices_prec, q_tor_slice_indices_prec, u_tor_slice_indices, q_tor_slice_indices)
      NULLIFY (u_tor_slice, q_tor_slice, u_conv_tor_slice, q_conv_tor_slice)

   END SUBROUTINE adaptivity_parallel_project_solution_3D
#endif
#endif

   SUBROUTINE adaptivity()

      CHARACTER(70)               :: param_adapt_char, count_adapt_char
      CHARACTER(1024)             :: file_name

      ! Start timing
      IF (utils%timing) THEN
         CALL cpu_TIME(timing%tps1)
         CALL system_CLOCK(timing%cks1, timing%clock_rate1)
      END IF

      IF (restart_adapt) THEN
         time%t = time%t + time%dt
         time%it = time%it + 1
         time%ik = time%ik + 1
         sol%Nt = sol%Nt + 1
      END IF

      !! Start refining procedure
      count_adapt = count_adapt + 1

      ! Deep copy previous Mesh and reference element
      CALL free_mesh_loc(Mesh_prec)
      CALL deep_copy_mesh_struct(Mesh,Mesh_prec)
      CALL free_reference_element_pol(refElPol_prec)
      CALL deep_copy_refel_struct(refElPol,refElPol_prec)

      ! Rescale to dimensional values
      Mesh%X = Mesh%X*phys%lscale
      Mesh_prec%X = Mesh_prec%X*phys%lscale

      ! Free elmat, mat, magnetic field, Jtor, puff (to save some memory for the adaptivity)
      CALL free_before_adaptivity()

      CALL adaptively_refine_mesh(mesh_name, count_adapt, order)

#ifdef TOR3D
      CALL define_toroidal_discretization()
#endif

      ! save new mesh
      WRITE (param_adapt_char, *) adapt%param_est
      WRITE (count_adapt_char, *) count_adapt
      file_name = TRIM(ADJUSTL(mesh_name)) // '_param' //TRIM(ADJUSTL(param_adapt_char)) // '_n' //TRIM(ADJUSTL(count_adapt_char))
      Mesh%X = Mesh%X/phys%lscale
      CALL HDF5_save_mesh_serial(Mesh, file_name)
      Mesh%X = Mesh%X*phys%lscale

#ifndef TOR3D
#ifndef PARALL
      CALL adaptivity_serial_project_solution_2D()
#else
      CALL adaptivity_parallel_project_solution_2D()
#endif
#else
#ifndef PARALL
      CALL adaptivity_serial_project_solution_3D()
#else
      CALL adaptivity_parallel_project_solution_3D()
#endif
#endif

#ifndef PARALL
      WRITE (*, *) "Number of elements previous mesh: ", Mesh_prec%Nelems
      WRITE (*, *) "Number of elements current mesh:  ", Mesh%Nelems
#else
      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "Number of elements previous mesh: ", Mesh_prec%Nel_glob
         WRITE (*, *) "Number of elements current mesh:  ", Mesh%Nel_glob
      END IF
#endif

#ifdef TOR3D
      CALL define_toroidal_discretization()
#endif

      ! Extract new trace solution
      CALL update_solution_arrays()

      ! Rescale back to adimensional values
      Mesh%X = Mesh%X/phys%lscale
      Mesh_prec%X = Mesh_prec%X/phys%lscale

      ! reset parameters
      CALL set_parameters()

      ! Re-initialize magnetic field (the Mesh is needed)
      CALL initialize_magnetic_field()

      ! Re-load magnetic field and Jtor
      CALL load_magnetic_field_Jtor()


      ! Re-Allocation and initialization of the elemental matrices
      CALL init_elmat()

      ! restart to first ever iteration
      matK%start = .TRUE.

      ! Re-Initialize shock capturing
      IF ((switch%shockcp .GT. 0) .OR. ((adapt%adaptivity) .AND. (adapt%shockcp_adapt .GT. 0)))  THEN
         CALL initializeShockCapturing()
      END IF

      ! Re-Compute first equation (definition of the gradient)
      CALL HDG_precalculatedfirstequation()

      !! UPDATE VARIABLES
      errNR_adapt = 1e10
      IF(restart_adapt) THEN
         time%t = time%t - time%dt
         time%it = time%it - 1
         time%ik = time%ik - 1
         sol%Nt = sol%Nt - 1
         restart_adapt = .FALSE.
      END IF

      IF (utils%timing) THEN
         CALL cpu_TIME(timing%tpe1)
         CALL system_CLOCK(timing%cke1, timing%clock_rate1)
         timing%runtadapt = timing%runtadapt + (timing%cke1-timing%cks1)/REAL(timing%clock_rate1)
         timing%cputadapt = timing%cputadapt + timing%tpe1-timing%tps1
      END IF

   END SUBROUTINE adaptivity


   SUBROUTINE free_before_adaptivity()
      IF(ASSOCIATED(phys%B)) DEALLOCATE(phys%B)
      IF(ASSOCIATED(phys%magnetic_flux)) DEALLOCATE(phys%magnetic_flux)
      IF(ASSOCIATED(phys%Bperturb)) DEALLOCATE(phys%Bperturb)
      IF(ASSOCIATED(phys%magnetic_psi)) DEALLOCATE(phys%magnetic_psi)
      IF(ASSOCIATED(phys%Jtor))     DEALLOCATE(phys%Jtor)
      IF(ASSOCIATED(phys%puff_exp)) DEALLOCATE(phys%puff_exp)

      NULLIFY(phys%magnetic_flux, phys%Bperturb, phys%magnetic_psi, phys%Jtor, phys%puff_exp)
      CALL free_el_mat()
      CALL free_mat()
   END SUBROUTINE free_before_adaptivity


   !************************************************
   ! Display results
   !************************************************
   SUBROUTINE displayResults()
      INTEGER                   :: ieq
      REAL*8, ALLOCATABLE       :: uphy(:, :)
      REAL*8                    :: Vmax(phys%npv), Vmin(phys%npv)
      nu = SIZE(sol%u)

      ALLOCATE (uphy(nu/phys%Neq, phys%npv))

      ! Compute physical variables
      CALL cons2phys(TRANSPOSE(RESHAPE(sol%u, (/phys%Neq, nu/phys%Neq/))), uphy)
      DO ieq = 1, phys%npv
         Vmax(ieq) = MAXVAL(uphy(:, ieq))
         Vmin(ieq) = MINVAL(uphy(:, ieq))
      END DO

#ifdef PARALL
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmax, phys%npv, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, Vmin, phys%npv, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (6, '(" * ", 60("-"), "*")')
         WRITE (6, '(" * Time (adimensional) = ", E12.5, 27X, " *")') time%t
         WRITE (6, '(" * Dt (adimensional)        = ", E12.5, 27X, " *")') time%dt
         WRITE (6, '(" * ", 45("^"), 14X, " *")')
         WRITE (6, '(" * ", 10("_"), "      Minimum   ", 4X, "     Maximum   ", 14X, " *")')
         DO ieq = 1, phys%npv
            WRITE (6, '(" * ", A7, " -->", ES16.8, 3X, ES16.8, 13X, " *")') &
            & TRIM(phys%phyVarNam(ieq)), Vmin(ieq), Vmax(ieq)
         END DO
         WRITE (6, '(" * ", 60("-"), "*")')
         WRITE (6, '(" * Time residual  = ", 1X, 2(E16.8, 2X), 13X, " *")') sol%tres(it)
         WRITE (6, *) '  '
         WRITE (6, *) '  '
         WRITE (6, *) '  '
      END IF
#ifdef PARALL
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      DEALLOCATE (uphy)

   END SUBROUTINE displayResults

   !************************************************
   ! Compute the residual
   !************************************************
   FUNCTION computeResidual(u, uref, coeff) RESULT(res)
      REAL*8                    :: u(:), uref(:)
      INTEGER                   :: nglo
      REAL*8                    :: res, sum2, coeff
#ifdef PARALL
      INTEGER                   :: ierr
#endif

      sum2 = SUM((u - uref)**2)
      nglo = SIZE(u)

#ifdef PARALL
      CALL mpi_allreduce(MPI_IN_PLACE, sum2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL mpi_allreduce(nu, nglo, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      res = SQRT(sum2)/SQRT(REAL(nglo))/coeff
   END FUNCTION computeResidual

   !************************************************
   ! Set the name of the solution
   !************************************************
   SUBROUTINE setSolName(save_name, mesh_name, it, convNR, convT)
      CHARACTER(LEN=1024), INTENT(INOUT):: save_name
      CHARACTER(LEN=1024), INTENT(IN)   :: mesh_name
      INTEGER, INTENT(IN)                  :: it
      LOGICAL, INTENT(IN)                  :: convNR, convT
      CHARACTER(LEN=20)                 :: Num
      INTEGER                             :: l, i

      ! At the beginning, the save name is the mesh name..
      save_name = TRIM(ADJUSTL(mesh_name))

      ! Eliminate path info
      l = LEN(save_name)
      i = INDEX(save_name, '/', .TRUE.)
      save_name = save_name(i + 1:l)
      ! look for the P of the order of the polynomial (P4 for example)
      i = INDEX(save_name, 'P', .TRUE.)
      WRITE (Num, "(i10)") refElPol%nDeg
      save_name(i+1:i+2) = TRIM(ADJUSTL(Num))


#ifdef TOR3D
      ! Add the number of toroidal elements
      WRITE (Num, "(i10)") numer%ntor
      Num = TRIM(ADJUSTL(Num))
      save_name = TRIM(ADJUSTL(save_name))//"_Ntor"//Num

      ! Add the poloidal interpolation in the toroidal direction
      WRITE (Num, "(i10)") numer%ptor
      Num = TRIM(ADJUSTL(Num))
      save_name = TRIM(ADJUSTL(save_name))//"Ptor"//Num

#endif
      ! Diffusion
      WRITE (Num, "(E10.3)") phys%diff_n*simpar%refval_diffusion
      save_name = TRIM(ADJUSTL(save_name))//"_DPe"//TRIM(ADJUSTL(Num))
#ifdef TEMPERATURE
      WRITE (Num, "(E10.3)") phys%diff_pari
      save_name = TRIM(ADJUSTL(save_name))//"_DPai"//TRIM(ADJUSTL(Num))

      WRITE (Num, "(E10.3)") phys%diff_pare
      save_name = TRIM(ADJUSTL(save_name))//"_DPae"//TRIM(ADJUSTL(Num))
#endif
      ! Complete the save name: if not converged the NR, I put NR + the iteration number
      IF (.NOT. convNR) THEN
         WRITE (Num, "(i10)") it
         Num = TRIM(ADJUSTL(Num))
         k = INDEX(Num, " ") - 1
         save_name = TRIM(ADJUSTL(save_name))//"_NR"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
      END IF

      ! Complete the save name: if not converged the time scheme, I put the iteration number
      IF (.NOT. convT) THEN
         WRITE (Num, "(i10)") it/utils%freqsave
         Num = TRIM(ADJUSTL(Num))
         k = INDEX(Num, " ") - 1
         save_name = TRIM(ADJUSTL(save_name))//"_"//REPEAT("0", 4 - k)//TRIM(ADJUSTL(Num))
      END IF

      IF (switch%decoup) THEN
         save_name = TRIM(ADJUSTL(save_name))//'_UNCP'
      END IF

      ! Add "Sol_"
#ifdef TOR3D
      save_name = 'Sol3D_'//TRIM(ADJUSTL(save_name))
#else
      save_name = 'Sol2D_'//TRIM(ADJUSTL(save_name))
#endif
      ! Add save_folder
      save_name = TRIM(ADJUSTL(input%save_folder))//save_name
   END SUBROUTINE setSolName

   SUBROUTINE compute_dt(errlstime)
      REAL*8, INTENT(IN)        :: errlstime

      IF ((errlstime*time%dt) < 1e-3) THEN
         WRITE (6, *) "******** Changing time step ***********"
         time%dt = time%dt*2.
      END IF
   END SUBROUTINE compute_dt


#ifdef TOR3D
   !**********************************************
   ! Definition of the toroidal discretization
   !**********************************************
   SUBROUTINE define_toroidal_discretization
      INTEGER                   :: i, ntorloc, itor, itorg, nnodes_toroidal
      INTEGER                   :: ind(numer%ptor + 1)
      REAL*8                    :: tdiv(numer%ntor + 1), tel(numer%ptor + 1), htor

      ! Toroidal discretization
      htor = numer%tmax/numer%ntor
      tdiv = 0.
      DO i = 1, numer%ntor
         tdiv(i + 1) = i*htor
      END DO
#ifdef PARALL
      IF (MPIvar%ntor .GT. 1) THEN
         ntorloc = numer%ntor/MPIvar%ntor + 1
      ELSE
         ntorloc = numer%ntor
      END IF
#else
      ntorloc = numer%ntor
#endif

      nnodes_toroidal = ntorloc + 1 + (numer%ptor - 1)*ntorloc
      Mesh%Nnodes_toroidal = nnodes_toroidal
      ALLOCATE (Mesh%toroidal(nnodes_toroidal))

      DO itor = 1, ntorloc
#ifdef PARALL
         itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
         IF (itorg == numer%ntor + 1) itorg = 1
#else
         itorg = itor
#endif
         tel = tdiv(itorg) + 0.5*(refElTor%coord1d+1)*(tdiv(itorg + 1) - tdiv(itorg))
         ind = (itor - 1)*numer%ptor + (/(i, i=1, numer%ptor + 1)/)
         Mesh%toroidal(ind) = tel
      END DO
   END SUBROUTINE define_toroidal_discretization
#endif

   SUBROUTINE update_solution_arrays()
      DEALLOCATE(sol%u_tilde)
      DEALLOCATE(sol%u_tilde0)
      ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
      ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
      CALL extractFaceSolution()
      sol%u_tilde0 = sol%u_tilde

      IF(restart_adapt) THEN
         DEALLOCATE(sol%u0)
         ALLOCATE (sol%u0(SIZE(sol%u), time%tis))
         sol%u0 = 0.
         sol%u0(:, 1) = sol%u
      END IF

      IF(.NOT. ((adapt%NR_adapt) .AND. (MOD(ir,adapt%freq_NR_adapt) .EQ. 0))) THEN
         ! Update check-point solution to the one projected on the new mesh
         DEALLOCATE(sol%u_conv)
         ALLOCATE(sol%u_conv(SIZE(sol%u)))
         sol%u_conv = sol%u
         DEALLOCATE(sol%q_conv)
         ALLOCATE(sol%q_conv(SIZE(sol%q)))
         sol%q_conv = sol%q
      END IF

   END SUBROUTINE update_solution_arrays

   SUBROUTINE update_solution()
      IF (time%tis .GT. 1 .AND. it .LT. time%tis) THEN
         DO is = it, 1, -1
            sol%u0(:, is + 1) = sol%u0(:, is)
         END DO
      ELSEIF (time%tis .GT. 1) THEN
         DO is = time%tis, 2, -1
            sol%u0(:, is) = sol%u0(:, is - 1)
         END DO
      END IF
      sol%u0(:, 1) = sol%u
   END SUBROUTINE update_solution

   SUBROUTINE update_uiter()
      IF(SIZE(uiter) .NE. SIZE(sol%u0(:,1))) THEN
         DEALLOCATE(uiter)
         ALLOCATE(uiter(SIZE(sol%u0(:,1))))
      END IF
      uiter = sol%u0(:, 1)
   END SUBROUTINE update_uiter

   SUBROUTINE update_uiter_qiter_best(uiter_best, qiter_best, u, q)
      REAL*8, POINTER, INTENT(OUT)  :: uiter_best(:), qiter_best(:)
      REAL*8, INTENT(IN)           :: u(:), q(:)

      IF(ASSOCIATED(uiter_best)) THEN
         IF(SIZE(uiter_best) .NE. SIZE(u)) THEN
            DEALLOCATE(uiter_best)
            DEALLOCATE(qiter_best)
            ALLOCATE(uiter_best(SIZE(u)))
            ALLOCATE(qiter_best(SIZE(q)))
            uiter_best = u
            qiter_best = q
         END IF
      ELSE
         ALLOCATE(uiter_best(SIZE(u)))
         ALLOCATE(qiter_best(SIZE(q)))
         uiter_best = u
         qiter_best = q
      END IF


   END SUBROUTINE update_uiter_qiter_best


   SUBROUTINE set_parameters()
      order = refElPol%nDeg
#ifdef TOR3D
      Ndim = 3                                               ! Number of dimensions
#ifdef PARALL
      ntorloc = numer%ntor/MPIvar%ntor
#else
      ntorloc = numer%ntor
#endif
      Nel = Mesh%Nelems*ntorloc                             ! Number of 3D elements
      Np = refElPol%Nnodes2D*refElTor%Nnodes1D             ! Number of nodes for each 3D element
      Nfl = refElPol%Nnodes1D*refElTor%Nnodes1D             ! Number of nodes in the lateral faces
      Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
      Nf = Mesh%Nfaces                                     ! Number of faces in the 2D mesh
#ifdef PARALL
      IF (MPIvar%ntor .GT. 1) THEN
         nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) + phys%Neq*refElPol%Nnodes2D*Mesh%Nelems ! Size of utilde
      ELSE
         nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
      END IF
#else
      nut = phys%Neq*ntorloc*(Nfl*Nf + refElPol%Nnodes2D*Mesh%Nelems) ! Size of utilde
#endif
      Nfp = refElPol%Nnodes2D*2 + refElPol%Nfaces*Nfl         ! Number of nodes in all the faces of a 3D element
      nu = phys%Neq*Nel*Np
#else
      Ndim = Mesh%ndim
      Nel = Mesh%Nelems
      Np = refElPol%Nnodes2D
      Nf = refElPol%Nfaces
      Nfp = refElPol%Nfacenodes*Nf
      nut = Mesh%Nfaces*Mesh%Nnodesperface*phys%Neq
      nu = Mesh%Nelems*Mesh%Nnodesperelem*phys%Neq
#endif
   END SUBROUTINE set_parameters

   SUBROUTINE initialize_solu0_uiter_uconv()
      ALLOCATE(uiter(SIZE(sol%u)))
      ALLOCATE(uiter_best(SIZE(sol%u)))
      ALLOCATE(qiter_best(SIZE(sol%q)))
      ALLOCATE(sol%u0(SIZE(sol%u), time%tis))
      ALLOCATE (sol%u_conv(SIZE(sol%u)))
      ALLOCATE (sol%q_conv(SIZE(sol%q)))
      sol%u0 = 0.
      sol%u0(:, 1) = sol%u
      uiter = sol%u
      uiter_best = sol%u
      qiter_best = sol%q
      sol%u_conv = sol%u
      sol%q_conv = sol%q
   END SUBROUTINE initialize_solu0_uiter_uconv

   SUBROUTINE update_u0(u)
      REAL*8, INTENT(IN)        :: u(:)

      IF(SIZE(sol%u0,1) .NE. SIZE(u)) THEN
         DEALLOCATE(sol%u0)
         ALLOCATE (sol%u0(SIZE(u), time%tis))
      END IF
      sol%u0 = 0.
      sol%u0(:, 1) = u
   END SUBROUTINE update_u0

   SUBROUTINE update_uconv_qconv(u_conv, q_conv)
      REAL*8, INTENT(IN)        :: u_conv(:), q_conv(:)
      IF(SIZE(sol%u_conv) .NE. SIZE(u_conv)) THEN
         DEALLOCATE(sol%u_conv)
         ALLOCATE(sol%u_conv(SIZE(u_conv)))
      END IF
      IF(SIZE(sol%q_conv) .NE. SIZE(q_conv)) THEN
         DEALLOCATE(sol%q_conv)
         ALLOCATE(sol%q_conv(SIZE(q_conv)))
      END IF
      sol%u_conv = u_conv
      sol%q_conv = q_conv
   END SUBROUTINE update_uconv_qconv

   SUBROUTINE update_dumpnr()
      numer%dumpnr = numer%dumpnr_min+(numer%dumpnr_max-numer%dumpnr_min)/2.*(1+TANH((ir-numer%dumpnr_n0)/numer%dumpnr_width))
   END SUBROUTINE update_dumpnr

   SUBROUTINE initialize_solution()
      IF (nb_args .EQ. 3) THEN
         ALLOCATE(sol%u_tilde(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
         ALLOCATE(sol%u_tilde0(phys%neq*Mesh%Nfaces*Mesh%Nnodesperface))
         CALL extractFaceSolution()
      ELSE IF (nb_args .EQ. 2) THEN
         ! restart simulation: load solution from file (the name is given in argument)
         CALL HDF5_load_solution(save_name)
         ALLOCATE(sol%u_tilde0(SIZE(sol%u_tilde)))
      ELSE
         CALL init_sol()
      END IF
   END SUBROUTINE initialize_solution

   SUBROUTINE load_mesh()

      IF (nb_args .EQ. 3) THEN
         ! find projection points
         CALL load_mesh_h5(mesh_name)
         ALLOCATE(xs(SIZE(Mesh%T,1)*SIZE(Mesh%T,2),2))
         xs = Mesh%X(colint(TRANSPOSE(Mesh%T)),:)
         CALL free_mesh()

         ! Load solution to project from
         CALL load_mesh_h5(mesh_name_proj)
         CALL create_reference_element(refElPol, 2, verbose = 1)
         CALL mesh_preprocess(ierr)
         CALL HDF5_load_solution(save_name)

         ! Project the solution
         !call projectSolutionDifferentMeshes(xs)
         DEALLOCATE(xs)
         CALL free_mesh()
         CALL free_reference_element()
      END IF

      IF(switch%readMeshFromSol) THEN
         CALL HDF5_load_mesh_from_solution(save_name)
      ELSE

         ! Load the mesh file from gmsh input or h5
         CALL load_gmsh_mesh(mesh_name, 1)

         CALL create_reference_element(refElPol, 2, verbose = 1)

         IF((switch%set_2d_order) .AND. (refElPol%nDeg .NE. switch%order_2d) ) THEN
            CALL mesh_preprocess_serial(ierr)
            CALL set_order_mesh(switch%order_2d)
            CALL free_reference_element_pol(refElPol)
            Mesh%X = Mesh%X*phys%lscale
            Mesh%xmax = Mesh%xmax*phys%lscale
            Mesh%xmin = Mesh%xmin*phys%lscale
            Mesh%ymax = Mesh%ymax*phys%lscale
            Mesh%ymin = Mesh%ymin*phys%lscale
         END IF

         IF(switch%gmsh2h5) THEN
            IF(MPIvar%glob_id .EQ. 0) THEN
               h5_filename = TRIM(ADJUSTL(mesh_name)) // '.h5'
               CALL convert_gmsh_to_hdf5(h5_filename, SIZE(Mesh%X,2), SIZE(Mesh%T,1), SIZE(Mesh%Tb,1), SIZE(Mesh%X,1), SIZE(Mesh%T,2), refElPol%nDeg+1, 0, Mesh%T, Mesh%X, Mesh%Tb, Mesh%boundaryFlag)
            END IF
         END IF


         IF ((switch%axisym .AND. switch%testcase .GE. 60 .AND. switch%testcase .LT. 80)) THEN
            Mesh%X(:,1) = Mesh%X(:,1) - geom%R0
         END IF
         IF (adapt%shockcp_adapt .GT. 0) THEN
            gmsh_filename      = TRIM(ADJUSTL(mesh_name))//'.msh'
            CALL copy_file(gmsh_filename,"./res/temp.msh")
         END IF

      END IF
   END SUBROUTINE load_mesh

   SUBROUTINE check_input_arguments()

      ! Check the number of input arguments
      nb_args = iargc()

      IF (nb_args .LT. 1) THEN
         PRINT *, " Error: the mesh name is needed"
         STOP
      END IF

      CALL getarg(1, mesh_name)
      mesh_name = ADJUSTL(TRIM(mesh_name))

      IF (nb_args .GT. 1) THEN
         CALL getarg(2, save_name)
         save_name = ADJUSTL(TRIM(save_name))
         IF(MPIvar%glob_id .EQ. 0) THEN
            PRINT *, " Restart simulation with solution: ", save_name
         END IF
      END IF

      IF (nb_args .GT. 2) THEN
         CALL getarg(3, mesh_name_proj)
         mesh_name_proj = ADJUSTL(TRIM(mesh_name_proj))
         PRINT *, " Projecting solution from: ", mesh_name_proj
      END IF

      IF (nb_args .GT. 3) THEN
         PRINT *, " Too many arguments "
         STOP
      END IF
   END SUBROUTINE check_input_arguments

   SUBROUTINE check_for_NaNs()
      DO i = 1, SIZE(sol%u)
         IF (ISNAN(sol%u(i))) THEN
            WRITE (6, *) "NaN detected. STOPPING."
            STOP
         END IF
      END DO
   END SUBROUTINE check_for_NaNs

#ifndef TOR3D
   SUBROUTINE project_u0_newmesh()
      INTEGER, POINTER          :: T_glob(:,:)
      REAL*8, POINTER           :: X_glob(:,:)

      NULLIFY(T_glob, X_glob)

      IF (.NOT. ((adapt%NR_adapt) .AND. (MOD(ir, adapt%freq_NR_adapt) .EQ. 0))) THEN
         ALLOCATE (u0_temp(SIZE(sol%u0, 1), SIZE(sol%u0, 2)))
         u0_temp = sol%u0
      ELSE
         ALLOCATE (u0_temp(SIZE(sol%u_conv), SIZE(sol%u0, 2)))
         u0_temp = 0
         u0_temp(:, 1) = sol%u_conv
      END IF

      DEALLOCATE (sol%u0)
      ALLOCATE (sol%u0(Mesh%Nelems*Mesh%Nnodesperelem*phys%neq, time%tis))
      sol%u0 = 0.

      Mesh%X = Mesh%X*phys%lscale
      Mesh_prec%X = Mesh_prec%X*phys%lscale

#ifndef PARALL
      CALL projectSolutionDifferentMeshes_general_arrays(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X, u1=u0_temp(:, 1), u2=sol%u0(:, 1))
#else
      ! Ghost cells are removed from the mesh and the solution and then the result is gathered over the processes

      CALL gather_mesh(Mesh_in=Mesh_prec, T_glob=T_glob, X_glob=X_glob)
      CALL gather_solution(Mesh_in=Mesh_prec, Nnodesperelem=Mesh_prec%Nnodesperelem, u_in=u0_temp(:, 1), u_glob=u_glob)

      CALL projectSolutionDifferentMeshes_general_arrays(T_glob, X_glob, Mesh%T, Mesh%X, u1=u_glob, u2=sol%u0(:, 1))
      DEALLOCATE (T_glob, X_glob, u_glob)
      NULLIFY (T_glob, X_glob, u_glob)
#endif

      Mesh%X = Mesh%X/phys%lscale
      Mesh_prec%X = Mesh_prec%X/phys%lscale

      DEALLOCATE (u0_temp)
   END SUBROUTINE project_u0_newmesh
#else
#ifndef PARALL
   SUBROUTINE project_u0_newmesh()
      INTEGER                   :: itor, ptor, iel, Ntor, Nelems, Nelems_prec, Nnodesperelem, Nnodesperelem_tor1d, Neq, start_range_u, end_range_u, size_utilde, Nfaces, Nnodesperface
      INTEGER, ALLOCATABLE      :: u0_tor_slice_indices_elem(:), u0_tor_slice_indices_prec(:), u0_tor_slice_indices(:)
      REAL*8, ALLOCATABLE           :: u0_tor_slice(:), u0_prec_tor_slice(:)

      Mesh%X = Mesh%X*phys%lscale
      Mesh_prec%X = Mesh_prec%X*phys%lscale

      Ntor = numer%ntor
      Nfaces = Mesh%Nfaces
      Nnodesperface = Mesh%Nnodesperface
      Nelems_prec = Mesh_prec%Nelems
      Nelems = Mesh%Nelems
      Nnodesperelem = Mesh%Nnodesperelem
      Nnodesperelem_tor1d = numer%ptor + 1
      Neq = phys%Neq
      Ndim = mesh%Ndim + 1

      sol%u = 0.
      sol%q = 0.

      ! indices corresponding to one element DoFs
      ALLOCATE (u0_tor_slice_indices_elem(Nnodesperelem*Neq))

      ! indices corresponding to the whole toroidal slice for the converged solution on the previous mesh
      ALLOCATE (u0_tor_slice_indices_prec(Nelems_prec*Nnodesperelem*Neq))

      ! solution for one toroidal slice on the previous mesh
      ALLOCATE (u0_prec_tor_slice(Nelems_prec*Nnodesperelem*Neq))

      ! indices corresponding to the whole toroidal slice for the new mesh
      ALLOCATE (u0_tor_slice_indices(Nelems*Nnodesperelem*Neq))

      ! solution for one toroidal slice on the new mesh
      ALLOCATE (u0_tor_slice(Nelems*Nnodesperelem*Neq))

      IF (.NOT. ((adapt%NR_adapt) .AND. (MOD(ir, adapt%freq_NR_adapt) .EQ. 0))) THEN
         ALLOCATE (u0_temp(SIZE(sol%u0, 1), SIZE(sol%u0, 2)))
         u0_temp = sol%u0
      ELSE
         ALLOCATE (u0_temp(SIZE(sol%u_conv), SIZE(sol%u0, 2)))
         u0_temp = 0
         u0_temp(:, 1) = sol%u_conv
      END IF

      DEALLOCATE (sol%u0)

      ALLOCATE (sol%u0(neq*Nelems*ntorloc*Nnodesperelem_tor1d*Nnodesperelem, 1)) ! Size of u, time%tis
      sol%u0 = 0.

      DO itor = 1, Ntor
         DO ptor = 1, Nnodesperelem_tor1d

            start_range_u = 1
            end_range_u = 1

            ! Extract the solution corresponding to the current toroidal slice
            DO iel = 1, Nelems_prec
               u0_tor_slice_indices_elem = (itor-1)*Nelems_prec*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1

               u0_tor_slice_indices_prec(start_range_u:end_range_u) = u0_tor_slice_indices_elem

               start_range_u = end_range_u + 1
            END DO

            u0_prec_tor_slice = u0_temp(u0_tor_slice_indices_prec, 1)

            CALL projectSolutionDifferentMeshes_general_arrays(Mesh_prec%T, Mesh_prec%X, Mesh%T, Mesh%X, u1=u0_prec_tor_slice, u2=u0_tor_slice)

            ! Put the projected solution in the right place in sol%u0
            start_range_u = 1
            end_range_u = 1

            DO iel = 1, Nelems
               u0_tor_slice_indices_elem = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1

               u0_tor_slice_indices(start_range_u:end_range_u) = u0_tor_slice_indices_elem

               start_range_u = end_range_u + 1
            END DO

            sol%u0(u0_tor_slice_indices, 1) = u0_tor_slice

         END DO
      END DO

      Mesh%X = Mesh%X/phys%lscale
      Mesh_prec%X = Mesh_prec%X/phys%lscale

      DEALLOCATE (u0_temp)
      DEALLOCATE (u0_tor_slice_indices_elem, u0_tor_slice_indices_prec, u0_tor_slice_indices, u0_prec_tor_slice, u0_tor_slice)
   END SUBROUTINE project_u0_newmesh
#else
   SUBROUTINE project_u0_newmesh()
      INTEGER                   :: itor, ptor, iel, Ntor, Nelems, Nelems_prec, Nelems_glob, Nelems_prec_glob, Nnodesperelem, Nnodesperelem_tor1d, Neq, start_range_u, end_range_u
      INTEGER, ALLOCATABLE      :: u0_tor_slice_indices_elem(:), u0_tor_slice_indices_prec(:), u0_tor_slice_indices(:)
      REAL*8, ALLOCATABLE       :: u0_tor_slice(:), u0_prec_tor_slice(:), u0_temp(:, :)
      REAL*8, POINTER           :: u0_glob(:)
      REAL*8, POINTER           :: X_glob(:, :)
      INTEGER, POINTER           :: T_glob(:, :)

      Mesh%X = Mesh%X*phys%lscale
      Mesh_prec%X = Mesh_prec%X*phys%lscale

      NULLIFY (u0_glob, X_glob, T_glob)

      ! Store u0 data before decomposition
      IF (.NOT. ((adapt%NR_adapt) .AND. (MOD(ir, adapt%freq_NR_adapt) .EQ. 0))) THEN
         ALLOCATE (u0_temp(SIZE(sol%u0, 1), SIZE(sol%u0, 2)))
         u0_temp = sol%u0
      ELSE
         ALLOCATE (u0_temp(SIZE(sol%u_conv), SIZE(sol%u0, 2)))
         u0_temp = 0
         u0_temp(:, 1) = sol%u_conv
      END IF

      Ntor = numer%ntor
      Nelems_prec_glob = Mesh_prec%Nel_glob
      Nelems_prec = Mesh_prec%Nelems
      Nelems_glob = Mesh%Nel_glob
      Nelems = Mesh%Nelems
      Nnodesperelem = Mesh%Nnodesperelem
      Nnodesperelem_tor1d = numer%ptor + 1
      Neq = phys%Neq

      DEALLOCATE (sol%u0)
      ALLOCATE (sol%u0(neq*Nelems*ntorloc*Nnodesperelem_tor1d*Nnodesperelem, 1))
      sol%u0 = 0.

      ! indices corresponding to one element DoFs
      ALLOCATE (u0_tor_slice_indices_elem(Nnodesperelem*Neq))

      ! indices corresponding to the whole toroidal slice for the converged solution on the previous mesh
      ALLOCATE (u0_tor_slice_indices_prec(Nelems_prec_glob*Nnodesperelem*Neq))

      ! solution for one toroidal slice on the previous mesh
      ALLOCATE (u0_prec_tor_slice(Nelems_prec_glob*Nnodesperelem*Neq))

      ! indices corresponding to the whole toroidal slice for the new mesh
      ALLOCATE (u0_tor_slice_indices(Nelems*Nnodesperelem*Neq))

      ! solution for one toroidal slice on the new mesh
      ALLOCATE (u0_tor_slice(Nelems*Nnodesperelem*Neq))

      ! Gather checkpoint solution and project it to the new mesh
      CALL gather_mesh(Mesh_prec, T_glob=T_glob, X_glob=X_glob)
      CALL gather_solution_3D(Mesh_in=Mesh_prec, u_in=u0_temp(:, 1), u_glob=u0_glob)

      DO itor = 1, Ntor
         DO ptor = 1, Nnodesperelem_tor1d

            start_range_u = 1
            end_range_u = 1

            ! Extract the solution corresponding to the current toroidal slice
            DO iel = 1, Nelems_prec_glob
               u0_tor_slice_indices_elem = (itor-1)*Nelems_prec_glob*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1

               u0_tor_slice_indices_prec(start_range_u:end_range_u) = u0_tor_slice_indices_elem

               start_range_u = end_range_u + 1
            END DO

            u0_prec_tor_slice = u0_glob(u0_tor_slice_indices_prec)

            CALL projectSolutionDifferentMeshes_general_arrays(T_glob, X_glob, Mesh%T, Mesh%X, u1=u0_prec_tor_slice, u2=u0_tor_slice)

            ! Put the projected solution in the right place in sol%u0
            start_range_u = 1
            end_range_u = 1

            DO iel = 1, Nelems
               u0_tor_slice_indices_elem = (itor-1)*Nelems*Nnodesperelem*Nnodesperelem_tor1d*Neq + (ptor-1)*Neq*Nnodesperelem + (iel-1)*Neq*Nnodesperelem*Nnodesperelem_tor1d + [(i, i = 1, Neq*Nnodesperelem)]

               end_range_u = start_range_u + Neq*Nnodesperelem - 1

               u0_tor_slice_indices(start_range_u:end_range_u) = u0_tor_slice_indices_elem

               start_range_u = end_range_u + 1
            END DO

            sol%u0(u0_tor_slice_indices, 1) = u0_tor_slice

         END DO
      END DO

      DEALLOCATE (u0_glob, X_glob, T_glob)
      NULLIFY (u0_glob, X_glob, T_glob)

      DEALLOCATE (u0_tor_slice_indices_elem, u0_tor_slice_indices_prec, u0_tor_slice_indices, u0_prec_tor_slice, u0_tor_slice)

      Mesh%X = Mesh%X/phys%lscale
      Mesh_prec%X = Mesh_prec%X/phys%lscale

      DEALLOCATE (u0_temp)
   END SUBROUTINE project_u0_newmesh
#endif
#endif

#ifdef TOR3D
   !**********************************************
   ! Definition of the toroidal discretization
   !**********************************************
   SUBROUTINE define_toroidal_discretization
      INTEGER                   :: i, ntorloc, itor, itorg, nnodes_toroidal
      INTEGER                   :: ind(numer%ptor + 1)
      REAL*8                    :: tdiv(numer%ntor + 1), tel(numer%ptor + 1), htor

      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (*, *) "*************************************************"
         WRITE (*, *) "*         DEFINE TOROIDAL DISCRETIZATION        *"
         WRITE (*, *) "*************************************************"
      END IF

      ! Toroidal discretization
      htor = numer%tmax/numer%ntor
      tdiv = 0.
      DO i = 1, numer%ntor
         tdiv(i + 1) = i*htor
      END DO
#ifdef PARALL
      IF (MPIvar%ntor .GT. 1) THEN
         ntorloc = numer%ntor/MPIvar%ntor + 1
      ELSE
         ntorloc = numer%ntor
      END IF
#else
      ntorloc = numer%ntor
#endif

      nnodes_toroidal = ntorloc + 1 + (numer%ptor - 1)*ntorloc
      Mesh%Nnodes_toroidal = nnodes_toroidal
      ALLOCATE (Mesh%toroidal(nnodes_toroidal))

      DO itor = 1, ntorloc
#ifdef PARALL
         itorg = itor + (MPIvar%itor - 1)*numer%ntor/MPIvar%ntor
         IF (itorg == numer%ntor + 1) itorg = 1
#else
         itorg = itor
#endif
         tel = tdiv(itorg) + 0.5*(refElTor%coord1d + 1)*(tdiv(itorg + 1) - tdiv(itorg))
         ind = (itor - 1)*numer%ptor + (/(i, i=1, numer%ptor + 1)/)
         Mesh%toroidal(ind) = tel
      END DO
   END SUBROUTINE define_toroidal_discretization
#endif

   SUBROUTINE reduce_diffusion()
      ! Update the diffusion, the elemental matrices and the solution
      IF (MPIvar%glob_id .EQ. 0) THEN
         WRITE (6, *) "************************************************"
         WRITE (6, *) "Reducing diffusion: ", phys%diff_n*switch%diffred*simpar%refval_diffusion
         WRITE (6, *) "************************************************"
      END IF
      phys%diff_n = phys%diff_n*switch%diffred
      phys%diff_u = phys%diff_u*switch%diffred
#ifdef TEMPERATURE
      phys%diff_e = phys%diff_e*switch%diffred
      phys%diff_ee = phys%diff_ee*switch%diffred
#endif
#ifdef VORTICITY
      phys%diff_vort = phys%diff_vort*switch%diffred
      phys%diff_pot = phys%diff_pot*switch%diffred
#endif
   END SUBROUTINE reduce_diffusion

   SUBROUTINE error_checking()

      IF (lssolver%pctype .EQ. 'PCMG') THEN
         IF (lssolver%mglevels .LT. 2) THEN
            IF (MPIvar%glob_id .EQ. 0) THEN
               WRITE (*, *) "Warning: number of multigrid levels too low, set to ", lssolver%mglevels, ". It should be at least 2. stop"
            END IF
            CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         END IF
         IF (lssolver%mglevels .GT. refElPol%Ndeg/2) THEN
            IF (MPIvar%glob_id .EQ. 0) THEN
               WRITE (*, *) "Warning: number of multigrid levels too low, set to ", lssolver%mglevels, ". It should not be more than half of the polynomial degree ", refElPol%Ndeg, ". stop"
            END IF
            !call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
         END IF
      END IF

#ifndef TOR3D
      IF (switch%ripple) THEN
         IF (MPIvar%glob_id .EQ. 0) THEN
            WRITE (*, *) "Ripple is activated but the code is not compiled in TOR3D mode. Please recompile with -DTOR3D flag or deactivate ripple. stop"
         END IF
         CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
      END IF
#endif

      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

   END SUBROUTINE error_checking

   SUBROUTINE print_timing_infos()
      ! Code timing
      IF (MPIvar%glob_id .EQ. 0) THEN
         IF (utils%timing) THEN
            cputtot = 1e-8
            runttot = 1e-8
            cputtot = timing%cputpre + timing%cputjac + timing%cputbcd + timing%cputmap + timing%cputass + timing%cputglb + timing%cputsol + timing%cputadapt
            runttot = timing%runtpre + timing%runtjac + timing%runtbcd + timing%runtmap + timing%runtass + timing%runtglb + timing%runtsol + timing%runtadapt
#ifdef PARALL
            cputtot = cputtot + timing%cputcom
            runttot = runttot + timing%runtcom
#endif
            WRITE(6, *) " "
            WRITE(6, *) " "
            WRITE(6, *) " "
            WRITE(6, '(" *", 90("*"), "**")')
            WRITE(6, '(" *", 36X, "CODE TIMING ( Nthreads = ",i2,")", 26X, " *")') Nthreads
            WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
            WRITE(6, '(" *", 2X,  "Precal. matr     : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
            &timing%cputpre,timing%cputpre/cputtot*100,timing%runtpre,timing%runtpre/runttot*100,timing%cputpre/timing%runtpre/Nthreads
            WRITE(6, '(" *", 2X,  "Jacobian         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputjac,timing%cputjac/cputtot*100,timing%runtjac,timing%runtjac/runttot*100,timing%cputjac/timing%runtjac/Nthreads
            WRITE(6, '(" *", 2X,  "Mapping          : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputmap,timing%cputmap/cputtot*100,timing%runtmap,timing%runtmap/runttot*100,timing%cputmap/timing%runtmap/Nthreads
            WRITE(6, '(" *", 2X,  "Boundary cond.   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputbcd,timing%cputbcd/cputtot*100,timing%runtbcd,timing%runtbcd/runttot*100,timing%cputbcd/timing%runtbcd/Nthreads
            WRITE(6, '(" *", 2X,  "Assembly         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputass,timing%cputass/cputtot*100,timing%runtass,timing%runtass/runttot*100,timing%cputass/timing%runtass/Nthreads
            WRITE(6, '(" *", 2X,  "Solve glob. syst.: ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputglb,timing%cputglb/cputtot*100,timing%runtglb,timing%runtglb/runttot*100,timing%cputglb/timing%runtglb/Nthreads
            WRITE(6, '(" *", 2X,  "Element solution : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputsol,timing%cputsol/cputtot*100,timing%runtsol,timing%runtsol/runttot*100,timing%cputsol/timing%runtsol/Nthreads
            IF(adapt%adaptivity .EQV. .TRUE.) THEN
               WRITE(6, '(" *", 2X,  "Adaptivity       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
               &timing%cputadapt,timing%cputadapt/cputtot*100,timing%runtadapt,timing%runtadapt/runttot*100,timing%cputadapt/timing%runtadapt/Nthreads
            END IF
#ifdef PARALL
            WRITE(6, '(" *", 2X,  "Communications   : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")')   &
            &timing%cputcom,timing%cputcom/cputtot*100,timing%runtcom,timing%runtcom/runttot*100,timing%cputcom/timing%runtcom/Nthreads
#endif
            WRITE(6, '(" *", 2X, "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")')   &
               cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
            WRITE(6, '(" *", 90("*"), "**")')
            WRITE(6, *) " "
            WRITE(6, *) " "
            WRITE(6, *) " "
         END IF

         IF (lssolver%timing) THEN
            cputtot = timing%clstime1 + timing%clstime2 + timing%clstime3 + timing%clstime4 + timing%clstime5 + timing%clstime6
            runttot = timing%rlstime1 + timing%rlstime2 + timing%rlstime3 + timing%rlstime4 + timing%rlstime5 + timing%rlstime6
            WRITE(6, '(" *", 90("*"), "**")')
            IF (lssolver%sollib .EQ. 1) THEN
               WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PASTIX ( Nthreads = ",i2,")", 14X, " *")') Nthreads
            ELSE IF (lssolver%sollib .EQ. 2) THEN
               WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PSBLAS ( Nthreads = ",i2,")", 14X, " *")') Nthreads
            ELSE IF (lssolver%sollib .EQ. 3) THEN
               WRITE(6, '(" *", 24X, "LINEAR SYSTEM SOLVER TIMING: PETSc ( Nthreads = ",i2,")", 14X, " *")') Nthreads
            END IF
            WRITE(6, '(" *", 28X, "Cpu-time  (% tot)",6X,    "Run-time  (% tot)   Speedup/Nthreads ", 2X, " *")')
            IF (lssolver%sollib .EQ. 1) THEN
               WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
               &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
               WRITE(6, '(" *", 2X,  "Check mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
               WRITE(6, '(" *", 2X,  "Anal. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
               WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
               WRITE(6, '(" *", 2X,  "LU decomp.       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
               WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime6,timing%clstime6/cputtot*100,timing%rlstime6,timing%rlstime6/runttot*100,timing%clstime6/timing%rlstime6/Nthreads
               WRITE(6, '(" *", 2X,  "Total time      : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",7X,F5.1 , 10X, " *")') &
               &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
            ELSEIF (lssolver%sollib .EQ. 2) THEN
               WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
               &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
               WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
               WRITE(6, '(" *", 2X,  "Build prec       : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
               WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
               WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime5,timing%clstime5/cputtot*100,timing%rlstime5,timing%rlstime5/runttot*100,timing%clstime5/timing%rlstime5/Nthreads
               WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
               &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
            ELSEIF (lssolver%sollib .EQ. 3) THEN
               WRITE(6, '(" *", 2X,  "Init. mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,  F4.1 , 10X, " *")') &
               &timing%clstime1,timing%clstime1/cputtot*100,timing%rlstime1,timing%rlstime1/runttot*100,timing%clstime1/timing%rlstime1/Nthreads
               WRITE(6, '(" *", 2X,  "Build mat        : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime2,timing%clstime2/cputtot*100,timing%rlstime2,timing%rlstime2/runttot*100,timing%clstime2/timing%rlstime2/Nthreads
               WRITE(6, '(" *", 2X,  "Fill vec         : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime3,timing%clstime3/cputtot*100,timing%rlstime3,timing%rlstime3/runttot*100,timing%clstime3/timing%rlstime3/Nthreads
               WRITE(6, '(" *", 2X,  "Solve            : ", ES16.3," ("F4.1 "%)",1X,ES14.3," ("F4.1 "%)",8X,F4.1 , 10X, " *")') &
               &timing%clstime4,timing%clstime4/cputtot*100,timing%rlstime4,timing%rlstime4/runttot*100,timing%clstime4/timing%rlstime4/Nthreads
               WRITE(6, '(" *", 2X,  "Total time       : ", ES16.3," ("F5.1 "%)",1X,ES13.3," ("F5.1 "%)",6X,F5.1 , 10X, " *")') &
               &cputtot,cputtot/cputtot*100,runttot,runttot/runttot*100,cputtot/runttot/Nthreads
            END IF
            WRITE (6, '(" *", 90("*"), "**")')
            WRITE (6, *) " "
            WRITE (6, *) " "
            WRITE (6, *) " "
         END IF
      END IF
   END SUBROUTINE print_timing_infos

   SUBROUTINE free_main()
#ifdef WITH_PASTIX
      USE solve_pastix, ONLY: terminate_mat_PASTIX
#endif
#ifdef WITH_PETSC
      USE solve_petsc, ONLY: terminate_PETSC, FinalizePETSC
#endif
#ifdef WITH_PSBLAS
      USE solve_psblas, ONLY: terminate_PSBLAS
#endif

      IF(ALLOCATED(uiter)) DEALLOCATE (uiter)
      IF(ASSOCIATED(uiter_best)) DEALLOCATE(uiter_best)
      IF(ASSOCIATED(qiter_best)) DEALLOCATE(qiter_best)
      NULLIFY(uiter_best)
      NULLIFY(qiter_best)
      IF(ALLOCATED(mkelms)) DEALLOCATE (mkelms)
      IF(ALLOCATED(sol%u0)) DEALLOCATE (sol%u0)
      IF(ALLOCATED(oscillations)) DEALLOCATE(oscillations)

      CALL free_splines(splines)

      IF (lssolver%sollib .EQ. 1) THEN
#ifdef WITH_PASTIX
         CALL terminate_mat_PASTIX()
         ! MPI finalization
         CALL MPI_finalize(IERR)
#endif
      ELSEIF (lssolver%sollib .EQ. 2) THEN
#ifdef WITH_PSBLAS
         CALL terminate_PSBLAS()
#endif
      ELSEIF (lssolver%sollib .EQ. 3) THEN
#ifdef WITH_PETSC
         CALL terminate_PETSC()
         CALL FinalizePETSC()
         CALL MPI_finalize(IERR)
#endif
      END IF

      CALL free_all
      CALL free_mesh
#ifdef PARALL
      CALL free_mesh_loc(Mesh_glob)
#endif
      CALL free_mesh_loc(Mesh_prec)
      CALL free_reference_element
      CALL free_reference_element_pol(refElPol_prec)

   END SUBROUTINE free_main


END MODULE Main_utils
