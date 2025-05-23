MODULE Solids_UPAS
  USE precision;
  USE Parallel_supplementary_Maths;
  CONTAINS
!-------------------------------
! Integrates the Element-Jacobian-Residual (Voight-notation)
! Total Lagrangian pressure-displacement-formulation (Active Strain)
!-------------------------------
SUBROUTINE STATIC_SOLIDUPAS_ELEMENT(Km, Rm, utemp, astrain, fibre, MATPROP   &
                                  , coord, Ni_p, dNi_u, Ni_u, weights, nprop &
                                  , ntots, ndofU, ndofP, ndim , nst  &
                                  , nip, nodU, nodp, MAP, UMAP)
  IMPLICIT NONE
  INTEGER                 :: Iel, Ig, i, j, k, l, p, q, m, n, o, s, t;
  INTEGER,  INTENT(IN)    :: nprop, ntots, ndim, nst;
  INTEGER,  INTENT(IN)    :: nip, nodU, nodp, ndofU, ndofP;
  INTEGER,  INTENT(IN)    :: MAP(ndim+1,nodU), UMAP(ndim,nodU)
  REAL(iwp),INTENT(IN)    :: astrain(nodU), fibre(ndim,ndim);
  REAL(iwp),INTENT(IN)    :: utemp(ntots), MATPROP(nprop);
  REAL(iwp),INTENT(IN)    :: coord(nodU,ndim), weights(nip)
  REAL(iwp),INTENT(IN)    :: dNi_u(ndim,nodU,nip), Ni_u(nodU,nip), Ni_p(nodP,nip)
  REAL(iwp),INTENT(INOUT) :: Km(ntots,ntots), Rm(ntots);
  REAL(iwp),PARAMETER     :: zero = 0._iwp, one = 1._iwp, two = 2._iwp, three=3._iwp;

  !Deformation related terms
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndim,ndim), det;
  REAL(iwp):: PK2, Ctang, S_bulk, K_bulk, S_kk, C_kkpp
  REAL(iwp):: auxm(ndim,nodU), pxm(nodP)
  REAL(iwp):: C_ijkl(nst,nst), S_ij(nst), CdefInv(ndim,ndim)
  REAL(iwp):: strain, gamma_f, gamma_n, gamma_s
  REAL(iwp):: Fdef(ndim,ndim), Fedef(ndim,ndim), F0def(ndim,ndim), press;
  REAL(iwp):: Cdef(ndim,ndim), F0inv(ndim,ndim), dE(nst,ndofU), d2E(nst,ndofU,ndofU);

  Km     = zero;  Rm      = zero;   dE    = zero; d2E    = zero;
  Cdef   = zero;  CdefInv = zero;   Fdef  = zero; C_kkpp = zero;
  C_ijkl = zero;  S_ij    = zero;   Fedef = zero; S_kk   = zero;
  Jac    = zero;  Jacinv  = zero;   det   = zero;

  !-----
  ! Calculate the stress,
  ! tangent-stiffness and
  ! deformation measures at
  ! the gauss-points
  !-----
  DO I=1,ndim
    auxm(I,:) = utemp(Map(I,:));
  ENDDO
  pxm = utemp(Map(ndim+1,1:ndofP));

  GAUSS_PTS: DO Ig = 1,nip
    ! Calculate the coordinate
    ! tranform from ref-to-physical
    Jac = MATMUL(dNi_u(:,:,Ig), coord(:,:))
    CALL Invert2(Jac, Jacinv, nDIM)
    det = DETERMINANT(Jac)

    ! Calculate the active strain and the
    ! displacement components of the def-tensors
    Fdef = zero;
    strain = DOT_PRODUCT(Ni_p(:,Ig),  astrain(1:nodP) )
    gamma_f = -strain;
    gamma_n =  4._iwp*gamma_f;
    gamma_s = (one/((one + gamma_f)*(one + gamma_n))) - one;
    DO i = 1,ndim
      DO j = 1,ndim
        Fdef(I,J) = Fdef(I,J) + DOT_PRODUCT(dNi_u(J,:,Ig),auxm(I,:));
        IF(I==J) Fdef(I,J) = Fdef(I,J) + one;
        F0def(i,j) = Kdelta(i,j) + gamma_f*fibre(i,1)*fibre(j,1)
        IF(ndim > 1) F0def(i,j) = F0def(i,j) + gamma_s*fibre(i,2)*fibre(j,2)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_n*fibre(i,3)*fibre(j,3)
      ENDDO
    ENDDO
    CALL INVERT2(F0def,F0inv,ndim)

    ! Calculate the stress, stiffness
    ! and strain measures/derivatives
    Fedef = MATMUL(Fdef,F0Inv)
    Cdef = MATMUL(TRANSPOSE(Fedef),Fedef)
    CALL GREENLAGRANGE_DERIVS_AS(Fedef,F0Inv,dNi_u(:,:,Ig),dE,d2E,UMAP,ndim,ndofU,nodU,nst)
    CALL INVERT2(Cdef,CdefInv(:,:),ndim);
    CALL MATMOD_NOH(C_ijkl,S_ij,MATPROP,Fedef,ndim,nst,nprop)

    ! Calculate the volumetric stress
    ! and the Bulk stiffness at the
    ! Gauss-Ptns
    S_kk    = zero;
    C_kkpp  = zero;
    DO m = 1,ndim
      S_kk = S_kk + S_ij(m);
      DO n = 1,ndim
        C_kkpp = C_kkpp + C_ijkl(m,n)
      ENDDO
    ENDDO

    !Calculating the pressure
    !at the Gauss points
    press = DOT_PRODUCT(Ni_p(:,Ig),pxm)

    !-----
    ! Integrate the residuals
    ! and the Jacobians
    !-----
    ! UU - Block
    ! Pure displacement block
    !
    DO M = 1,ndofU
      DO s = 1,nst; CALL VOIGHT_ITERATOR(s,I,J,nst)
        PK2 = S_ij(s) + press*CdefInv(I,J) - S_kk*Kdelta(I,J)/three;
        Rm(M) = Rm(M) + dE(s,M)*PK2*det*weights(Ig);
        Km(M,1:ndofU) = Km(M,1:ndofU) + PK2*d2E(s,:,M)*det*weights(Ig);
        DO t = 1,nst; CALL VOIGHT_ITERATOR(t,K,L,nst)
          Ctang = C_ijkl(s,t) - two*press*CdefInv(K,I)*CdefInv(J,L)
          Ctang = Ctang - C_kkpp*Kdelta(I,J)*Kdelta(K,L)/three;
          Km(M,1:ndofU) = Km(M,1:ndofU) + dE(s,M)*Ctang*dE(t,:)*det*weights(Ig);
        ENDDO
      ENDDO
    ENDDO

    ! UP-PU-PP-Block(s)
    ! Pressure-displacement mixed block
    !
    DO M = (ndofU+1),(ndofU+ndofP)
      S_bulk = S_kk/C_kkpp;
      K_bulk = C_kkpp/three;
      Rm(M) = Rm(M) + (S_bulk - press/K_bulk)*Ni_p(M,Ig)*det*weights(Ig);

      ! UP-PU Blocks
      DO s = 1,nst; CALL VOIGHT_ITERATOR(s,I,J,nst)
        Km(1:ndofU,M) = Km(1:ndofU,M) + CdefInv(I,J)*dE(s,:)*Ni_p(M,Ig)*det*weights(Ig);
        Km(M,1:ndofU) = Km(M,1:ndofU) + CdefInv(I,J)*dE(s,:)*Ni_p(M,Ig)*det*weights(Ig);
      ENDDO

      ! PP-Block
      p = (ndofU+1);
      q = (ndofU+ndofP);
      Km(p:q,M) = Km(p:q,M) - (one/K_bulk)*Ni_p(M,Ig)*Ni_p(:,Ig)*det*weights(Ig);
    ENDDO
  ENDDO GAUSS_PTS
  RETURN
END SUBROUTINE STATIC_SOLIDUPAS_ELEMENT

!-----------------
! Green-Lagrange-strain derivatives Active strain Variant
!-----------------
SUBROUTINE GREENLAGRANGE_DERIVS_AS(Fdef,F0Inv,dNi,dE,d2E,Map,ndim,ntots,nod,nst)
  IMPLICIT NONE
  INTEGER                :: i, j, k, l, p, q, r, m, n;
  INTEGER,   INTENT(IN)   :: nod, ntots, ndim, nst;
  INTEGER,   INTENT(IN)   :: Map(ndim,nod);
  REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), F0Inv(ndim,ndim), dNi(ndim,nod);
  REAL(iwp), INTENT(INOUT):: d2E(nst,ntots,ntots), dE(nst,ntots);
  REAL(iwp)               :: a, d2temp(ntots,ntots), d1temp(ntots);

  d2E = 0._iwp; dE = 0._iwp;
  DO m = 1,nst; CALL VOIGHT_ITERATOR(m, i, j, nst)
    DO n = 1,nst; CALL VOIGHT_ITERATOR(n, p, q, nst)
      d2temp = 0._iwp
      d1temp = 0._iwp
      a = 0.5_iwp*F0Inv(p,i)*F0Inv(q,j)
      IF(p/=q) a = F0Inv(p,i)*F0Inv(q,j)
      DO k = 1,ndim
        d1temp(Map(k,:)) = d1temp(Map(k,:)) + (dNi(p,:)*Fdef(k,q) + dNi(q,:)*Fdef(k,p))
        DO l = 1,ndim
          DO r = 1,nod
            d2temp(Map(l,r),Map(k,:)) = d2temp(Map(l,r),Map(k,:)) + (dNi(p,r)*dNi(q,:) + dNi(q,r)*dNi(p,:))
          ENDDO
        ENDDO
      ENDDO
      dE(m,:) = dE(m,:) + a*d1temp
      d2E(m,:,:) = d2E(m,:,:) + a*d2temp
    ENDDO
    IF(I/=J) dE(m,:) = 2._iwp*dE(m,:);
    IF(I/=J) d2E(m,:,:) = 2._iwp*d2E(m,:,:);
  ENDDO
  RETURN
ENDSUBROUTINE GREENLAGRANGE_DERIVS_AS

!-------------------------------
! Neo-Hookean
!-------------------------------
SUBROUTINE MATMOD_NOH(C_ijkl,S_ij,Matprops,Fdef,ndim,nst,nprop)
   IMPLICIT NONE
   INTEGER                 :: i, j, k, l, m, n;
   INTEGER,   INTENT(IN)   :: ndim, nst, nprop
   REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), Matprops(nprop)
   REAL(iwp), INTENT(INOUT):: C_ijkl(nst,nst), S_ij(nst);
   REAL(iwp)               :: LJ3, J3, mu0, mu, lmbda, Y, nu;
   REAL(iwp)               :: Cdef(ndim,ndim), C_inv(ndim,ndim);

   J3    = determinant(Fdef)
   LJ3   = DLOG(J3)
   Y     = 1.00_iwp;    !   Y     = Matprops(1); !1.00_iwp    !
   nu    = 0.49999_iwp; !   nu    = Matprops(2); !0.49999_iwp !
   mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
   lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))
   Cdef = MATMUL(TRANSPOSE(Fdef),Fdef)
   CALL INVERT2(Cdef,C_inv,ndim)
   DO m = 1,nst; CALL VOIGHT_ITERATOR(m, i, j, nst)
     S_ij(m) = mu*(Kdelta(i,j) - C_inv(i,j));
     DO n = 1,nst; CALL VOIGHT_ITERATOR(n, k, l, nst)
       C_ijkl(m,n) =  2._iwp*mu*C_inv(i,k)*C_inv(j,l);
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE MATMOD_NOH

!-------------------------------------
! SDF map
!-------------------------------------
SUBROUTINE SDF_MAP(Map, nod, ndim)
  IMPLICIT NONE
  INTEGER               :: I, J;
  INTEGER, INTENT(IN)   :: nod, ndim
  INTEGER, INTENT(INOUT):: Map(ndim,nod);

  DO I = 1,ndim
    DO J = 1,nod
      Map(I,J) = (I-1)*nod + J;
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE SDF_MAP
!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Solids_UPAS
