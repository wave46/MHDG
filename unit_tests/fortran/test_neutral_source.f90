#if !defined(NEUTRAL) || !defined(TEMPERATURE)
#error This check requires a temperature model with neutrals.
#endif
PROGRAM test_neutral_source
  USE globals, ONLY: phys
  USE physics, ONLY: neutral_rt, compute_niz, compute_dniz_dU, compute_nrec, compute_dnrec_dU, &
       compute_fGammacx, compute_dfGammacx_dU, compute_fGammarec, compute_dfGammarec_dU, &
       compute_fEiiz, compute_dfEiiz_dU, compute_fEirec, compute_dfEirec_dU, &
       compute_fEicx, compute_dfEicx_dU
#ifdef NEUTRALGAMMA
  USE physics, ONLY: compute_fGammaN, compute_dfGammaN_dU, compute_fEiN, compute_dfEiN_dU
#endif
  IMPLICIT NONE
#ifdef NEUTRALGAMMA
  INTEGER, PARAMETER :: neq=6
#else
  INTEGER, PARAMETER :: neq=5
#endif
  REAL*8 :: U(neq), baseline(neq)
  INTEGER :: state_case

  phys%neq=neq
  phys%idx_rhon_eq=5
  phys%idx_gamman_eq=0
#ifdef NEUTRALGAMMA
  phys%idx_gamman_eq=6
#endif
  phys%ionization_ion_energy_fraction=0.3d0
  neutral_rt%recombination_energy=0.7d0
  baseline=0.d0
  baseline(1:5)=[0.8d0,0.2d0,1.2d0,1.7d0,0.4d0]
#ifdef NEUTRALGAMMA
  baseline(6)=0.05d0
#endif

  DO state_case=1,6
    U=baseline
    SELECT CASE(state_case)
    CASE(2)
      ! The negative-neutral Gauss state from the failed bootstrap.
      U(1:5)=[0.10611394602704305d0,-0.4138862937912107d0, &
           0.8311894963983084d0,4.145941890942085d0,-0.03230005490588605d0]
    CASE(3)
      U(5)=0.d0
    CASE(4)
      U(3)=-0.2d0
    CASE(5)
      U(1)=-0.2d0
      U(2)=0.d0
    CASE(6)
      U(5)=0.5d0*DBLE(1.e-20)
    END SELECT
    CALL check_prefactors(U,state_case==1)
    CALL check_assembled_source(U)
  ENDDO
  WRITE(*,'(A,I0,A)') 'Neutral sources: PASS; ',neq,' equations, six states'

CONTAINS

  SUBROUTINE evaluate_prefactor(kind,U,f,df)
    INTEGER, INTENT(IN) :: kind
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(OUT) :: f,df(:)
    SELECT CASE(kind)
    CASE(1)
      CALL compute_niz(U,f)
      CALL compute_dniz_dU(U,df)
    CASE(2)
      CALL compute_nrec(U,f)
      CALL compute_dnrec_dU(U,df)
    CASE(3)
      CALL compute_fGammacx(U,f)
      CALL compute_dfGammacx_dU(U,df)
    CASE(4)
      CALL compute_fGammarec(U,f)
      CALL compute_dfGammarec_dU(U,df)
    CASE(5)
      CALL compute_fEiiz(U,f)
      CALL compute_dfEiiz_dU(U,df)
    CASE(6)
      CALL compute_fEirec(U,f)
      CALL compute_dfEirec_dU(U,df)
    CASE(7)
      CALL compute_fEicx(U,f)
      CALL compute_dfEicx_dU(U,df)
    END SELECT
  END SUBROUTINE evaluate_prefactor

  REAL*8 FUNCTION difference_step(U,j) RESULT(h)
    REAL*8, INTENT(IN) :: U(:)
    INTEGER, INTENT(IN) :: j
    h=3.d-6*MAX(1.d0,ABS(U(j)))
    ! Keep the stencil on one side of a floor; no derivative exists at the join.
    IF (j==1 .OR. j==3 .OR. j==5) h=MIN(h,0.25d0*ABS(U(j)-DBLE(1.e-20)))
  END FUNCTION difference_step

  SUBROUTINE check_prefactors(U,check_homogeneity)
    REAL*8, INTENT(IN) :: U(:)
    LOGICAL, INTENT(IN) :: check_homogeneity
    REAL*8 :: f,df(neq),plus(neq),minus(neq),scratch(neq),fp,fm,h,fd,scale
    INTEGER :: kind,j
    DO kind=1,7
      CALL evaluate_prefactor(kind,U,f,df)
      IF (check_homogeneity) THEN
        IF (ABS(DOT_PRODUCT(df,U)-2.d0*f)>1.d-13*MAX(ABS(f),1.d-25)) &
             ERROR STOP 'positive-state prefactor homogeneity'
      ENDIF
      DO j=1,neq
        h=difference_step(U,j)
        plus=U; plus(j)=plus(j)+h
        minus=U; minus(j)=minus(j)-h
        CALL evaluate_prefactor(kind,plus,fp,scratch)
        CALL evaluate_prefactor(kind,minus,fm,scratch)
        fd=(fp-fm)/(2.d0*h)
        scale=MAX(ABS(df(j)),ABS(fd),ABS(f),1.d-25)
        IF (ABS(fd-df(j))>8.d-8*scale) THEN
          WRITE(*,*) 'prefactor, column, state:',kind,j,U
          WRITE(*,*) 'implemented, finite difference:',df(j),fd
          ERROR STOP 'atomic prefactor derivative'
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE check_prefactors

  SUBROUTINE evaluate_source(U,source,Sn,Sn0)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8, INTENT(OUT) :: source(:),Sn(:,:),Sn0(:)
    REAL*8 :: f(7),df(neq,7),rates(5),drates(neq,5),cooling,dcooling(neq)
    REAL*8 :: fGammaN,dfGammaN(neq),fEiN,dfEiN(neq)
    INTEGER :: kind
    DO kind=1,7
      CALL evaluate_prefactor(kind,U,f(kind),df(:,kind))
    ENDDO
    ! Nonhomogeneous rates exercise the affine formula independently of Euler identities.
    rates=[2.d0,0.3d0,1.1d0,0.7d0,0.2d0]*(1.d0+0.05d0*U(1)+0.02d0*U(4))
    drates=0.d0
    drates(1,:)=[2.d0,0.3d0,1.1d0,0.7d0,0.2d0]*0.05d0
    drates(4,:)=[2.d0,0.3d0,1.1d0,0.7d0,0.2d0]*0.02d0
    cooling=0.4d0+0.03d0*U(4)
    dcooling=0.d0; dcooling(4)=0.03d0
    fGammaN=0.d0; dfGammaN=0.d0
    fEiN=0.d0; dfEiN=0.d0
#ifdef NEUTRALGAMMA
    CALL compute_fGammaN(U,fGammaN)
    CALL compute_dfGammaN_dU(U,dfGammaN)
    ! This unregularized kinetic-energy closure requires nonzero neutral density.
    IF (U(5)>1.d-12) THEN
      CALL compute_fEiN(U,fEiN)
      CALL compute_dfEiN_dU(U,dfEiN)
    ENDIF
#endif
    CALL assemblyNeutral(U=U,niz=f(1),dniz_dU=df(:,1),nrec=f(2),dnrec_dU=df(:,2), &
         sigmaviz=rates(1),dsigmaviz_dU=drates(:,1),sigmavrec=rates(2),dsigmavrec_dU=drates(:,2), &
         fGammacx=f(3),dfGammacx_dU=df(:,3),fGammarec=f(4),dfGammarec_dU=df(:,4), &
#ifdef NEUTRALGAMMA
         fGammaN=fGammaN,dfGammaN_dU=dfGammaN, &
#endif
         sigmavcx=rates(3),dsigmavcx_dU=drates(:,3),fEiiz=f(5),dfEiiz_dU=df(:,5), &
         fEirec=f(6),dfEirec_dU=df(:,6),fEicx=f(7),dfEicx_dU=df(:,7), &
#ifdef NEUTRALGAMMA
         fEiN=fEiN,dfEiN_dU=dfEiN, &
#endif
         Sn=Sn,Sn0=Sn0,sigmavEiz=rates(4),dsigmavEiz_dU=drates(:,4), &
         sigmavErec=rates(5),dsigmavErec_dU=drates(:,5),cooling_factor=cooling,dcooling_factor_dU=dcooling)
    source=0.d0
    source(1)=f(1)*rates(1)-f(2)*rates(2)
    source(2)=-f(3)*rates(3)-f(4)*rates(2)+fGammaN*(rates(1)+rates(3))
    source(3)=phys%ionization_ion_energy_fraction*f(5)*rates(1)-f(6)*rates(2)-f(7)*rates(3) + &
         fEiN*(rates(1)+rates(3))
    source(4)=-f(1)*rates(4)-f(2)*(rates(5)+cooling-neutral_rt%recombination_energy*rates(2))
    source(5)=-source(1)
#ifdef NEUTRALGAMMA
    source(6)=-source(2)
#endif
  END SUBROUTINE evaluate_source

  SUBROUTINE check_assembled_source(U)
    REAL*8, INTENT(IN) :: U(:)
    REAL*8 :: source(neq),Sn(neq,neq),Sn0(neq),plus(neq),minus(neq),sp(neq),sm(neq)
    REAL*8 :: scratch(neq,neq),affine_scratch(neq),h,fd(neq)
    INTEGER :: j
    CALL evaluate_source(U,source,Sn,Sn0)
    IF (MAXVAL(ABS(-MATMUL(Sn,U)-Sn0-source))>2.d-14*MAX(1.d0,MAXVAL(ABS(source)))) &
         ERROR STOP 'assembled affine source value'
    IF (ANY(Sn(1,:)+Sn(5,:)/=0.d0) .OR. Sn0(1)+Sn0(5)/=0.d0) &
         ERROR STOP 'plasma-neutral particle cancellation'
    DO j=1,neq
      h=difference_step(U,j)
      plus=U; plus(j)=plus(j)+h
      minus=U; minus(j)=minus(j)-h
      CALL evaluate_source(plus,sp,scratch,affine_scratch)
      CALL evaluate_source(minus,sm,scratch,affine_scratch)
      fd=(sp-sm)/(2.d0*h)
      IF (ANY(ABS(fd+Sn(:,j))>8.d-8*MAX(1.d0,ABS(fd),ABS(Sn(:,j))))) &
           ERROR STOP 'assembled source Jacobian'
    ENDDO
  END SUBROUTINE check_assembled_source

#include "neutral_source_assembly_routine.F"

END PROGRAM test_neutral_source
