MODULE solve_psblas
  USE matrices_types
  USE types
  USE globals
  USE psb_base_mod
  USE psb_prec_mod
  USE psb_krylov_mod
  USE psb_util_mod
  USE mpi_omp
#ifdef WITH_MLD2P4
  USE mld_prec_mod
#endif
  IMPLICIT NONE

  TYPE PSBLAS_STRUC ! A type to store structures for the PSBLAS library
     TYPE(psb_dspmat_type) :: mat          ! sparse matrix for PSBLAS
#ifndef WITH_MLD2P4
     TYPE(psb_dprec_type)  :: prec         ! preconditioner
#else
     TYPE(mld_dprec_type)  :: prec
#endif
     TYPE(psb_desc_type)   :: desc_a       ! descriptor
     TYPE(psb_d_vect_type) :: x, b         ! dense vector
     INTEGER(psb_ipk_)     :: ictxt, iam, np, info   ! parallel environment
  END TYPE PSBLAS_STRUC

  TYPE(PSBLAS_STRUC) :: matPSBLAS

CONTAINS

  !                real*8  :: tps, tpe
  !                integer :: cks,clock_rate,cke
  !                if (lssolver%timing) then
  !                                call cpu_time(tps)
  !                                call system_clock(cks,clock_rate)
  !                endif

  !                if (lssolver%timing) then
  !                        call cpu_time(tpe)
  !                        call system_clock(cke,clock_rate)
  !                        timing%rlstime1 = timing%rlstime1+(cke-cks)/real(clock_rate)
  !                        timing%clstime1 = timing%clstime1+tpe-tps
  !                end if

  !***********************************************
  ! Initialization of matrix instance of the
  ! PSBLAS library
  !***********************************************
  SUBROUTINE init_mat_PSBLAS(matPSBLAS, matK)
    TYPE(MAT_CSR_TYP) :: matK
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    INTEGER(psb_lpk_), ALLOCATABLE   :: ind(:)
    INTEGER(psb_lpk_), ALLOCATABLE   :: rows(:), cols(:)
    INTEGER(psb_ipk_)               :: info
    INTEGER(psb_ipk_)               :: iam, np, ncoef
    INTEGER                         :: irow, i
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    !******************************************
    !Initialize parallel environment for psblas
    !******************************************
    CALL psb_init(matPSBLAS%ictxt)

    ! get information about the PSBLAS parallel environment
    ! psb_info: in{icontxt} out{iam,np}
    ! iam is the current process
    ! np number of processes
    CALL psb_info(matPSBLAS%ictxt, iam, np)

    IF (iam < 0) THEN
       ! This should not happen, but just in case
       CALL psb_exit(matPSBLAS%ictxt)
       STOP
    ENDIF
    !  if(psb_get_errstatus() /= 0) goto 9999
    !  name='mld_s_pde3d'
    CALL psb_set_errverbosity(itwo)

    !******************************************
    !  Allocate a communication descriptor
    !******************************************
    !   call psb_cdall(matPSBLAS%ictxt,matPSBLAS%desc_a,info,nl=matK%n)

    ALLOCATE (rows(matK%n))
    rows = matK%loc2glob
    CALL psb_cdall(matPSBLAS%ictxt, matPSBLAS%desc_a, info, vl=rows)
    DEALLOCATE (rows)

    !******************************************
    !  Create the communication descriptor
    !******************************************
    DO irow = 1, matK%n
       ncoef = matK%rowptr(irow + 1) - matK%rowptr(irow)
       ALLOCATE (rows(ncoef), ind(ncoef), cols(ncoef))
       ind = (/(i, i=matK%rowptr(irow), matK%rowptr(irow + 1) - 1)/)
       rows = matK%loc2glob(irow)
       cols = matK%cols(ind)
       CALL psb_cdins(ncoef, rows, cols, matPSBLAS%desc_a, info)
       IF (info .NE. psb_success_) WRITE (6, *) "failed filling matrix!"
       DEALLOCATE (rows, ind, cols)
    END DO

    !******************************************
    !  Assembly the communication descriptor
    !******************************************
    CALL psb_cdasb(matPSBLAS%desc_a, info)

    !******************************************
    !  Allocate matrix,unknowns and rhs
    !******************************************
    IF (info == psb_success_) CALL psb_spall(matPSBLAS%mat, matPSBLAS%desc_a, info, nnz=matK%nnz)
    IF (info == psb_success_) CALL psb_geall(matPSBLAS%x, matPSBLAS%desc_a, info)
    IF (info == psb_success_) CALL psb_geall(matPSBLAS%b, matPSBLAS%desc_a, info)

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime1 = timing%rlstime1 + (cke - cks)/REAL(clock_rate)
       timing%clstime1 = timing%clstime1 + tpe - tps
    END IF
  END SUBROUTINE init_mat_PSBLAS

  !***********************************************
  ! Here the matrix is passed to the PSBLAS
  ! library
  !***********************************************
  SUBROUTINE build_mat_PSBLAS(matPSBLAS, matK)
    TYPE(MAT_CSR_TYP) :: matK
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    INTEGER(psb_lpk_), ALLOCATABLE   :: ind(:)
    INTEGER(psb_lpk_), ALLOCATABLE   :: rows(:), cols(:)
    REAL(psb_dpk_), ALLOCATABLE     :: vals(:)
    INTEGER(psb_ipk_)               :: info, irow, i, ncoef
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0


    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    !******************************************
    !  Initialize matrix to zero
    !******************************************
    CALL psb_sprn(matPSBLAS%mat, matPSBLAS%desc_a, info)

    !******************************************
    !  Fill matrix
    !******************************************
    DO irow = 1, matK%n
       ncoef = matK%rowptr(irow + 1) - matK%rowptr(irow)
       ALLOCATE (rows(ncoef), ind(ncoef), cols(ncoef), vals(ncoef))
       ind = (/(i, i=matK%rowptr(irow), matK%rowptr(irow + 1) - 1)/)
       rows = matK%loc2glob(irow)
       cols = matK%cols(ind)
       vals = matK%vals(ind)
       CALL psb_spins(ncoef, rows, cols, vals, matPSBLAS%mat, matPSBLAS%desc_a, info)
       IF (info .NE. psb_success_) WRITE (6, *) "failed filling matrix!"
       DEALLOCATE (rows, ind, cols, vals)
    END DO

    !******************************************
    !  Assembly matrix
    !******************************************
    CALL psb_spasb(matPSBLAS%mat, matPSBLAS%desc_a, info, dupl=psb_dupl_err_)
    IF (info .NE. psb_success_) THEN
       WRITE (6, *) "failed assemblying the matrix"
       STOP
    ENDIF

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime2 = timing%rlstime2 + (cke - cks)/REAL(clock_rate)
       timing%clstime2 = timing%clstime2 + tpe - tps
    END IF
  END SUBROUTINE build_mat_PSBLAS

  !***********************************************
  ! Build the preconditioner
  !
  !***********************************************
#ifndef WITH_MLD2P4
  SUBROUTINE build_prec_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    INTEGER(psb_ipk_)  :: info
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    CALL matPSBLAS%prec%init(matPSBLAS%ictxt, lssolver%ptype, info)
    CALL matPSBLAS%prec%build(matPSBLAS%mat, matPSBLAS%desc_a, info)

    IF (info .NE. psb_success_) THEN
       WRITE (6, *) "Something wrong in the preconditioner construction"
       STOP
    END IF
    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime3 = timing%rlstime3 + (cke - cks)/REAL(clock_rate)
       timing%clstime3 = timing%clstime3 + tpe - tps
    END IF
  END SUBROUTINE build_prec_PSBLAS
#else
  SUBROUTINE build_prec_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    INTEGER(psb_ipk_)  :: info
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke
    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    CALL matPSBLAS%prec%init(matPSBLAS%ictxt, lssolver%ptype, info)
    SELECT CASE (lssolver%ptype)
    CASE ('NONE')
       ! Do nothing, keep defaults

    CASE ('DIAG', 'JACOBI', 'GS', 'FBGS')
       ! 1-level sweeps from "outer_sweeps"
       CALL matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)

    CASE ('BJAC')
       CALL matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
       CALL matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
       CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
       CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)

    CASE ('AS')
       CALL matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)
       CALL matPSBLAS%prec%set('sub_ovr', lssolver%novr, info)
       CALL matPSBLAS%prec%set('sub_restr', lssolver%restr, info)
       CALL matPSBLAS%prec%set('sub_prol', lssolver%prol, info)
       CALL matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
       CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
       CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)

    CASE ('ML')
       ! multilevel preconditioner

       CALL matPSBLAS%prec%set('ml_cycle', lssolver%mlcycle, info)
       CALL matPSBLAS%prec%set('outer_sweeps', lssolver%outer_sweeps, info)
       IF (lssolver%csize > 0)&
            & CALL matPSBLAS%prec%set('min_coarse_size', lssolver%csize, info)
       IF (lssolver%mncrratio > 1)&
            & CALL matPSBLAS%prec%set('min_cr_ratio', lssolver%mncrratio, info)
       IF (lssolver%maxlevs > 0)&
            & CALL matPSBLAS%prec%set('max_levs', lssolver%maxlevs, info)
       IF (lssolver%athres >= dzero) &
            & CALL matPSBLAS%prec%set('aggr_thresh', lssolver%athres, info)
       !    if (lssolver%thrvsz>0) then
       !      do k=1,min(lssolver%thrvsz,size(prec%precv)-1)
       !        call matPSBLAS%prec%set('aggr_thresh',     lssolver%athresv(k),  info,ilev=(k+1))
       !      end do
       !    end if

       CALL matPSBLAS%prec%set('aggr_prol', lssolver%aggr_prol, info)
       CALL matPSBLAS%prec%set('par_aggr_alg', lssolver%par_aggr_alg, info)
       CALL matPSBLAS%prec%set('aggr_ord', lssolver%aggr_ord, info)
       CALL matPSBLAS%prec%set('aggr_filter', lssolver%aggr_filter, info)

       IF ((psb_toupper(lssolver%smther) /= 'ML').AND.(psb_toupper(lssolver%smther) /= 'NONE')) THEN
          CALL matPSBLAS%prec%set('smoother_type', lssolver%smther, info)
       END IF

       CALL matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps, info)

       SELECT CASE (psb_toupper(lssolver%smther))
       CASE ('NONE', 'DIAG', 'JACOBI', 'GS', 'FBGS', 'ML')
          ! do nothing
       CASE ('BJAC')
          CALL matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
          CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
          CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)
       CASE ('AS')
          CALL matPSBLAS%prec%set('sub_ovr', lssolver%novr, info)
          CALL matPSBLAS%prec%set('sub_restr', lssolver%restr, info)
          CALL matPSBLAS%prec%set('sub_prol', lssolver%prol, info)
          CALL matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
          CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
          CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)
       CASE default
          CALL matPSBLAS%prec%set('sub_ovr', lssolver%novr, info)
          CALL matPSBLAS%prec%set('sub_restr', lssolver%restr, info)
          CALL matPSBLAS%prec%set('sub_prol', lssolver%prol, info)
          CALL matPSBLAS%prec%set('sub_solve', lssolver%solve, info)
          CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill, info)
          CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr, info)
       END SELECT

       IF ((psb_toupper(lssolver%smther2) /= 'ML').AND.(psb_toupper(lssolver%smther2) /= 'NONE')) THEN
          CALL matPSBLAS%prec%set('smoother_type', lssolver%smther2, info, pos='post')
       END IF

       CALL matPSBLAS%prec%set('smoother_sweeps', lssolver%jsweeps2, info, pos='post')

       SELECT CASE (psb_toupper(lssolver%smther2))
       CASE ('NONE', 'DIAG', 'JACOBI', 'GS', 'FBGS','ML')
          ! do nothing
       CASE ('BJAC')
          CALL matPSBLAS%prec%set('sub_solve', lssolver%solve2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr2, info, pos='post')
       CASE ('AS')
          CALL matPSBLAS%prec%set('sub_ovr', lssolver%novr2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_restr', lssolver%restr2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_prol', lssolver%prol2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_solve', lssolver%solve2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr2, info, pos='post')
       CASE default
          CALL matPSBLAS%prec%set('sub_restr', lssolver%restr2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_prol', lssolver%prol2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_solve', lssolver%solve2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_fillin', lssolver%fill2, info, pos='post')
          CALL matPSBLAS%prec%set('sub_iluthrs', lssolver%thr2, info, pos='post')
       END SELECT

       IF (psb_toupper(lssolver%csolve) /= 'NONE') THEN
          CALL matPSBLAS%prec%set('coarse_solve', lssolver%csolve, info)
          IF (psb_toupper(lssolver%csolve) == 'BJAC') &
               &  CALL matPSBLAS%prec%set('coarse_subsolve', lssolver%csbsolve, info)
          CALL matPSBLAS%prec%set('coarse_mat', lssolver%cmat, info)
          CALL matPSBLAS%prec%set('coarse_fillin', lssolver%cfill, info)
          CALL matPSBLAS%prec%set('coarse_iluthrs', lssolver%cthres, info)
          CALL matPSBLAS%prec%set('coarse_sweeps', lssolver%cjswp, info)
       ENDIF

    END SELECT

    ! build the preconditioner
    CALL matPSBLAS%prec%hierarchy_build(matPSBLAS%mat, matPSBLAS%desc_a, info)
    CALL matPSBLAS%prec%smoothers_build(matPSBLAS%mat, matPSBLAS%desc_a, info)

    IF (info .NE. psb_success_) THEN
       WRITE (6, *) "Something wrong in the preconditioner construction"
       STOP
    END IF
    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime3 = timing%rlstime3 + (cke - cks)/REAL(clock_rate)
       timing%clstime3 = timing%clstime3 + tpe - tps
    END IF
  END SUBROUTINE build_prec_PSBLAS
#endif

  !***********************************************
  ! Fill RHS and initial guess
  !
  !***********************************************
  SUBROUTINE fill_vec_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    INTEGER(psb_lpk_), ALLOCATABLE   :: rows(:)
    INTEGER(psb_ipk_)  :: info, n
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0


    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF

    !     call matPSBLAS%b%set(rhs%vals)
    !     call matPSBLAS%x%set(sol%u_tilde)

    n = matPSBLAS%mat%get_nrows()
    CALL matPSBLAS%x%zero()
    CALL matPSBLAS%b%zero()
    ALLOCATE (rows(n))
    rows = matK%loc2glob
    CALL psb_geins(n, rows, rhs%vals, matPSBLAS%b, matPSBLAS%desc_a, info)
#ifdef PARALL
    CALL getFirstGuessInParallel(matPSBLAS, rows)
#else
    CALL psb_geins(n, rows, sol%u_tilde, matPSBLAS%x, matPSBLAS%desc_a, info)
#endif
    DEALLOCATE (rows)

    !     n = matPSBLAS%mat%get_nrows()
    !     call matPSBLAS%x%zero()
    !     call matPSBLAS%b%zero()
    !     allocate(ind(n))
    !     ind = (/(i,i=1,n)/)
    !     call psb_geins(n,ind,rhs%vals,matPSBLAS%b,matPSBLAS%desc_a,info)
    !     call psb_geins(n,ind,sol%u_tilde,matPSBLAS%x,matPSBLAS%desc_a,info)
    !     deallocate(ind)

    !******************************************
    !  Assembly rhs and initial guess
    !******************************************

    ! Init. guess
    CALL psb_geasb(matPSBLAS%x, matPSBLAS%desc_a, info)
    IF (info .NE. psb_success_) THEN
       WRITE (6, *) "Something wrong in the initial guess assembly"
       STOP
    END IF

    ! RHS
    CALL psb_geasb(matPSBLAS%b, matPSBLAS%desc_a, info)
    IF (info .NE. psb_success_) THEN
       WRITE (6, *) "Something wrong in the rhs assembly"
       STOP
    END IF

    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime4 = timing%rlstime4 + (cke - cks)/REAL(clock_rate)
       timing%clstime4 = timing%clstime4 + tpe - tps
    END IF
  END SUBROUTINE fill_vec_PSBLAS

  !***********************************************
  ! Solve call for PSBLAS, assign output to
  ! u_tilde
  !***********************************************
  SUBROUTINE solve_mat_PSBLAS(matPSBLAS)
    TYPE(PSBLAS_STRUC) :: matPSBLAS
    INTEGER(psb_ipk_)     :: info   ! parallel environment
    INTEGER(psb_ipk_)     :: iter
    REAL(psb_dpk_)        :: err
    REAL*8  :: tps, tpe
    INTEGER :: cks, clock_rate, cke

    cks = 0
    cke = 0
    clock_rate = 0


    IF (lssolver%timing) THEN
       CALL cpu_TIME(tps)
       CALL system_CLOCK(cks, clock_rate)
    ENDIF


    !write(6,*) "saving matrix"
    !!call mm_mat_write(matPSBLAS%mat,'sparse_matrix' ,ierr,1,'psblas_mat')
    !!call mm_array_write(matPSBLAS%b%v%v,'rhs',ierr,1,'psblas_b')
    !call mm_array_write(matPSBLAS%x%v%v,'solution_pastix',ierr,1,'psblas_solpastix')
    !!call mm_array_write(rhs%vals,ierr,1,'psblas_b')
    !!call mm_array_write(sol%u_tilde,ierr,1,'psblas_x0')
    !write(6,*) "done saving matrix"
    !stop
    !******************************************
    !  Solve
    !******************************************
    CALL psb_krylov(lssolver%kmethd, matPSBLAS%mat, matPSBLAS%prec, matPSBLAS%b, matPSBLAS%x, lssolver%tol, matPSBLAS%desc_a, info,&
         & itmax=lssolver%itmax, iter=iter, err=err, itrace=lssolver%itrace, istop=lssolver%istop, irst=lssolver%rest)

    !******************************************
    !  Store solution
    !******************************************
    !   sol%u_tilde = (1.-numer%dumpnr)*sol%u_tilde + numer%dumpnr*matPSBLAS%x%get_vect()
    IF (lssolver%timing) THEN
       CALL cpu_TIME(tpe)
       CALL system_CLOCK(cke, clock_rate)
       timing%rlstime5 = timing%rlstime5 + (cke - cks)/REAL(clock_rate)
       timing%clstime5 = timing%clstime5 + tpe - tps
    END IF
  END SUBROUTINE solve_mat_PSBLAS

  !******************************************
  !  Free memory
  !******************************************
  SUBROUTINE terminate_PSBLAS()
    INTEGER(psb_ipk_)     :: info

    CALL psb_gefree(matPSBLAS%b, matPSBLAS%desc_a, info)
    CALL psb_gefree(matPSBLAS%x, matPSBLAS%desc_a, info)
    CALL psb_spfree(matPSBLAS%mat, matPSBLAS%desc_a, info)
    CALL matPSBLAS%prec%free(info)
    CALL psb_cdfree(matPSBLAS%desc_a, info)
    CALL psb_exit(matPSBLAS%ictxt)
  END SUBROUTINE terminate_PSBLAS

#ifdef PARALL
  !********************************************
  !  Get the first guess in parallel
  !********************************************
  SUBROUTINE getFirstGuessInParallel(matPSBLAS, rows)
    TYPE(PSBLAS_STRUC)            :: matPSBLAS
    INTEGER(psb_lpk_), INTENT(in) :: rows(:)
    INTEGER(psb_ipk_)             :: info

    INTEGER*4                     :: i, j, n
#ifdef PARALL
    INTEGER*4                     :: Neq, Nfp
    INTEGER*4                     :: ct, indg(refElPol%Nfacenodes*phys%Neq), indl(refElPol%Nfacenodes*phys%Neq)
#ifdef TOR3D
    INTEGER*4                     :: Nfl, N2d, Np2d, itor, Nfaces, Nfdir, Nghostf, Nghoste, dd, ddl, ntorass
    INTEGER*4                     :: indgp(refElPol%Nnodes2D*phys%Neq),indlp(refElPol%Nnodes2D*phys%Neq),indgt(refElTor%Nfl*phys%Neq),indlt(refElTor%Nfl*phys%Neq)
#endif
#endif

#ifdef TOR3D
    ! ********* Parallel 3D ***********
    IF (MPIvar%ntor .GT. 1) THEN
       ntorass = numer%ntor/MPIvar%ntor
    ELSE
       ntorass = numer%ntor
    ENDIF

    Neq = phys%Neq
    Nfl = refElTor%Nfl
    N2d = Mesh%Nelems
    Np2d = refElPol%Nnodes2D
    Nfaces = Mesh%Nfaces
    Nfdir = Mesh%Ndir
    Nghostf = Mesh%nghostfaces
    Nghoste = Mesh%nghostelems
    IF (MPIvar%ntor .GT. 1 .AND. MPIvar%itor .EQ. MPIvar%ntor) THEN
       ct = 0
       dd = 1 + ntorass*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
       ddl = 1
       n = Np2D*Neq
       DO i = 1, N2D
          IF (Mesh%ghostelems(i) .EQ. 1) CYCLE
          ct = ct + 1
          indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
          indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
          CALL psb_geins(n, rows(indlp), sol%u_tilde(indgp), matPSBLAS%x, matPSBLAS%desc_a, info)

          !        aux_sol(indlp) = sol%u_tilde(indgp)
       END DO
    END IF
    DO itor = 1, ntorass
       ct = 0
       dd = 1 + (itor - 1)*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
       ddl = 1 + (itor - 1)*((N2D-Nghoste)*Np2D+(Nfaces - Nfdir - Nghostf)*Nfl)*Neq
       IF (MPIvar%ntor .GT. 1 .AND. MPIvar%itor .NE. MPIvar%ntor) THEN
          ddl = ddl - (N2D-Nghoste)*Np2d*Neq
       END IF
       n = Np2D*Neq
       DO i = 1, N2D
          IF (Mesh%ghostelems(i) .EQ. 1) CYCLE
          IF (MPIvar%ntor .GT. 1 .AND. itor == 1) CYCLE
          ct = ct + 1
          indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
          indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
          CALL psb_geins(n, rows(indlp), sol%u_tilde(indgp), matPSBLAS%x, matPSBLAS%desc_a, info)
          !        sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
       END DO
       ct = 0
       dd = dd + (N2D*Np2D)*Neq
       ddl = ddl + ((N2D-Nghoste)*Np2D)*Neq

       n = Neq*Nfl
       DO i = 1, Mesh%Nfaces
          IF (Mesh%ghostfaces(i) .EQ. 1) CYCLE
          ct = ct + 1
          indgt = dd + (i - 1)*Nfl*Neq + (/(j, j=0, Neq*Nfl - 1)/)
          indlt = ddl + (ct - 1)*Nfl*Neq + (/(j, j=0, Neq*Nfl - 1)/)
          CALL psb_geins(n, rows(indlt), sol%u_tilde(indgt), matPSBLAS%x, matPSBLAS%desc_a, info)
          !        sol%u_tilde(indgt) = (1.-numer%dumpnr)*sol%u_tilde(indgt) + numer%dumpnr*aux_sol(indlt)
       END DO
    END DO
    IF (MPIvar%ntor .GT. 1 .AND. MPIvar%itor .NE. MPIvar%ntor) THEN
       ct = 0
       dd = 1 + ntorass*(N2D*Np2D+(Nfaces - Nfdir)*Nfl)*Neq
       ddl = 1 + ntorass*((N2D-Nghoste)*Np2D+(Nfaces - Nfdir - Nghostf)*Nfl)*Neq
       ddl = ddl - (N2D-Nghoste)*Np2d*Neq
       n = Np2D*Neq
       DO i = 1, N2D
          IF (Mesh%ghostelems(i) .EQ. 1) CYCLE
          ct = ct + 1
          indgp = dd + (i - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
          indlp = ddl + (ct - 1)*Np2D*Neq + (/(j, j=0, Np2D*Neq - 1)/)
          CALL psb_geins(n, rows(indlp), sol%u_tilde(indgp), matPSBLAS%x, matPSBLAS%desc_a, info)
          !        sol%u_tilde(indgp) = (1.-numer%dumpnr)*sol%u_tilde(indgp) + numer%dumpnr*aux_sol(indlp)
       END DO
    END IF

#else
    ! ********* Parallel 2D ***********
    ct = 0
    Neq = phys%Neq
    Nfp = Mesh%Nnodesperface
    n = Neq*Nfp
    DO i = 1, Mesh%Nfaces
       IF (Mesh%ghostfaces(i) .EQ. 1) CYCLE
       ct = ct + 1
       indg = (i - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
       indl = (ct - 1)*Neq*Nfp + (/(j, j=1, Neq*Nfp)/)
       CALL psb_geins(n, rows(indl), sol%u_tilde(indg), matPSBLAS%x, matPSBLAS%desc_a, info)
       !     sol%u_tilde(indg) = (1.-numer%dumpnr)*sol%u_tilde(indg) + numer%dumpnr*aux_sol(indl)
    END DO
#endif

  END SUBROUTINE getFirstGuessInParallel
#endif

END MODULE solve_psblas
