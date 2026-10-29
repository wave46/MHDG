! Filled by check_heat_flux_assembly.py with unmodified production heat blocks.
MODULE extracted_heat_assembly
  USE iso_c_binding
  USE globals
  USE physics
  IMPLICIT NONE
  INTEGER, PARAMETER :: Neq=5,Ndim=2,Npfl=1
CONTAINS
  SUBROUTINE probe(state,gradient,route,enabled,ju,jq,affine,phi) BIND(C)
    REAL(c_double), INTENT(IN) :: state(Neq),gradient(Ndim,Neq)
    INTEGER(c_int), VALUE :: route,enabled
    REAL(c_double), INTENT(OUT) :: ju(Neq,Neq),jq(Neq,2*Neq),affine(Neq),phi(2)
    REAL*8 :: ue(Neq),uf(Neq),ufg(Neq),qpr(Ndim,Neq),b(Ndim),bg(Ndim),n(Ndim),bn
    REAL*8 :: Vveci(Neq),Vvece(Neq),dVi(Neq,Neq),dVe(Neq,Neq),Taui(Ndim,Neq),Taue(Ndim,Neq)
    REAL*8 :: Alphai,Alphae,dAlpha_dUi(Neq),dAlpha_dUe(Neq),coefi,coefe,gmi,gme
    REAL*8 :: q_sh_i,q_sh_e,flux_limiter_i,flux_limiter_e,Abohm(2,Neq),kmult(1,1),kmultf(1)
    REAL*8 :: Auu(1,1,Neq*Neq),Auq(1,1,2*Neq*Neq),rhs(1,Neq),s,ds_dU(Neq),W,Zet(Ndim,Neq)
    REAL*8 :: QdW3(Ndim,Neq),QdW4(Ndim,Neq),dSohmic_dU(Neq),Sohmic,Jtor
    REAL*8 :: NNxy(1,1),NNi(1,1),NNbb(1),Ni(1),Nif(1),Nfbn(1),NNif(1,1),NiNi(1,1),NxyzNi(1,1,Ndim)
    LOGICAL :: isdir
    INTEGER :: i,j,k,z,iel,indi(1),indj(1),ind_if(1),ind_jf(1),ind_kf(1),ind_asf(1),ind_ash(1)
    INTEGER :: ind_ff(Neq),ind_fe(Neq),ind_fg(2*Neq)

    isdir=.FALSE.
    phys%Mref=1.3d0;phys%epn=2.5d0
    phys%diff_pari=.4d0;phys%diff_pare=.7d0
    phys%T_fluxlim_maxi=.6d0;phys%T_fluxlim_maxe=.9d0
    phys%c_fli=.2d0;phys%c_fle=.3d0
    simpar%refval_temperature=1.d0
    switch%testcase=50;switch%flux_limiter=enabled/=0;switch%ohmicsrc=.FALSE.
    b=[.2d0,.4d0];bg=b;n=[.6d0,.8d0];bn=DOT_PRODUCT(b,n)
    ue=state;uf=state;ufg=state;qpr=gradient
    CALL computeVi(state,Vveci);CALL computeVe(state,Vvece)
    CALL compute_dV_dUi(state,dVi);CALL compute_dV_dUe(state,dVe)
    Taui=MATMUL(qpr,dVi);Taue=MATMUL(qpr,dVe)
    gmi=DOT_PRODUCT(MATMUL(qpr,Vveci),b);gme=DOT_PRODUCT(MATMUL(qpr,Vvece),b)
    coefi=phys%diff_pari*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
    coefe=phys%diff_pare*(2.d0/(3.d0*phys%Mref))**(1.d0+phys%epn)
    IF (.NOT.ALLOCATED(elMat%Aul)) THEN
      ALLOCATE(elMat%Aul(Neq,Neq,1),elMat%Auq(Neq,2*Neq,1),elMat%S(Neq,1))
      ALLOCATE(elMat%All(Neq,Neq,1),elMat%Alq(Neq,2*Neq,1),elMat%fh(Neq,1),elMat%Aul_dir(Neq,1))
    ENDIF
    elMat%Aul=0;elMat%Auq=0;elMat%S=0;elMat%All=0;elMat%Alq=0;elMat%fh=0;elMat%Aul_dir=0
    Auu=0;Auq=0;rhs=0;Abohm=0
    s=0;ds_dU=0;W=0;Zet=0;QdW3=0;QdW4=0;dSohmic_dU=0;Sohmic=0;Jtor=0;NxyzNi=0
    iel=1;ind_asf=0;ind_ash=0;ind_ff=[(i,i=1,Neq)];ind_fe=ind_ff;ind_fg=[(i,i=1,2*Neq)]
    Nif=1;Ni=1;Nfbn=bn;NNif=1;NiNi=1;NNxy=1;NNi=1;NNbb=1
    SELECT CASE(route)
    CASE(1)
      DO i=3,4
! PRODUCTION_VOLUME
      ENDDO
      ju=RESHAPE(Auu,[Neq,Neq]);jq=RESHAPE(Auq,[Neq,2*Neq]);affine=rhs(1,:)
    CASE(2)
      DO i=3,4
        ind_if=i
! PRODUCTION_INTERIOR
      ENDDO
      ju=elMat%Aul(:,:,1);jq=elMat%Auq(:,:,1);affine=elMat%S(:,1)
      IF (ANY(ju/=elMat%All(:,:,1)).OR.ANY(jq/=elMat%Alq(:,:,1)).OR.ANY(affine/=elMat%fh(:,1))) &
           ERROR STOP 'interior element/trace heat differs'
    CASE(3)
      DO i=3,4
        ind_if=i
! PRODUCTION_EXTERIOR
      ENDDO
      ju=elMat%Aul(:,:,1);jq=elMat%Auq(:,:,1);affine=elMat%S(:,1)
    CASE(4)
! PRODUCTION_BOHM
      ju=elMat%All(:,:,1);jq=elMat%Alq(:,:,1);affine=elMat%fh(:,1)
    CASE(5)
! PRODUCTION_PRESCRIBED
      ju=0;jq=0;affine=elMat%Aul_dir(:,1)
    END SELECT
    phi=[flux_limiter_i,flux_limiter_e]
  END SUBROUTINE probe
END MODULE extracted_heat_assembly
