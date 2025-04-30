MODULE Solids_UPAS
  USE precision;
  USE new_library;
  USE Parallel_supplementary_Maths;
  CONTAINS

!-------------------------------
! Integrates the Element-Jacobian-Residual (Voight-notation)
! Total Lagrangian pressure-displacement-formulation (Active Strain)
!-------------------------------
SUBROUTINE STATIC_SOLIDUPAS_ELEMENT(Km, Rm, utemp, astrain, fibre, MATPROP   &
                                  , coord, dNi_p, Ni_p, dNi_u, Ni_u, weights &
                                  , nprop, ntots, ndofU, ndofP, ndim , nst   &
                                  , nip, nodU, nodp)
  IMPLICIT NONE
  INTEGER                :: Iel, Ig, i, j, k, l, p, q, m, n, o, s, t;
  INTEGER,  INTENT(IN)   :: nprop, ntots, ndim, nst;
  INTEGER,  INTENT(IN)   :: nip, nodU, nodp, ndofU, ndofP
  REAL(iwp),INTENT(IN)   :: astrain(nodU,nel_pp), fibre(ndim,ndim,nel_pp);
  REAL(iwp),INTENT(IN)   :: utemp(ntots,nel_pp), MATPROP(nprop);
  REAL(iwp),INTENT(IN)   :: coord(nodU,ndim,nel_pp), weights(nip)
  REAL(iwp),INTENT(IN)   :: dNi_p(ndim,nodP,nip), Ni_p(nodP,nip)
  REAL(iwp),INTENT(IN)   :: dNi_u(ndim,nodU,nip), Ni_u(nodU,nip)
  REAL(iwp),INTENT(INOUT):: Km(ntots,ntots,nel_pp), Rm(ntots,nel_pp);
  REAL(iwp),PARAMETER    :: zero = 0._iwp, one = 1._iwp, two = 2._iwp, three=3._iwp;

  !Deformation related terms
  REAL(iwp):: PK2, Ctang, S_bulk, K_bulk
  REAL(iwp):: S_kk(nip), C_kkpp(nip), auxm(ndim,nodU), pxm(nodP)
  REAL(iwp):: C_ijkl(nst,nst,nip), S_ij(nst,nip), CdefInv(ndim,ndim,nip)
  REAL(iwp):: strain, gamma_f, gamma_n, gamma_s, press(nip);
  REAL(iwp):: Fdef(ndim,ndim), Fedef(ndim,ndim), F0def(ndim,ndim);
  REAL(iwp):: Cdef(ndim,ndim), F0inv(ndim,ndim);
  REAL(iwp):: dE(nst,ndofU,nip), d2E(nst,ndofU,ndofU,nip);
  INTEGER  :: Map(ndim+1,nodU);

  !Shape function and derivatives
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndim,ndim), det(nip);

  Km     = zero;  Rm      = zero;
  C_ijkl = zero;  S_ij    = zero;   Fedef  = zero;
  dE     = zero;  d2E     = zero;
  Jac    = zero;  Jacinv  = zero;   det    = zero;
  Cdef   = zero;  CdefInv = zero;   Fdef   = zero;

  !===
  !
  ! Integrate all the elements
  !
  !====
  ELEMENTS:DO Iel = 1, nel_pp
    !-----
    ! Calculate the stress,
    ! tangent-stiffness and
    ! deformation measures at
    ! the gauss-points
    !-----
    CALL SDF_MAP(Map, nodU, ndim+1)
    DO I=1,ndim
      auxm(I,:) = utemp(Map(I,:),Iel);
    ENDDO
    pxm = utemp(Map(ndim+1,1:nodP),Iel);

    GAUSS_PTS1: DO Ig = 1,nip
      ! Calculate the coordinate
      ! tranform from ref-to-physical
      Jac = MATMUL(dNi(:,:,Ig), coord(:,:,Iel))
      CALL Invert2(Jac, Jacinv, nDIM)
      det(Ig) = DETERMINANT(Jac)

      ! Calculate the active strain
      ! components of the tensors
      strain = DOT_PRODUCT(Ni_p(:,Ig),  astrain(1:nodP,Iel) )
      gamma_f = -strain;
      gamma_n =  4._iwp*gamma_f;
      gamma_s = (one/((one + gamma_f)*(one + gamma_n))) - one;
      DO i = 1,ndim
        DO j = 1,ndim
          F0def(i,j) = Kdelta(i,j) + gamma_f*fibre(i,1,Iel)*fibre(j,1,Iel)
          IF(ndim > 1) F0def(i,j) = F0def(i,j) + gamma_s*fibre(i,2,Iel)*fibre(j,2,Iel)
          IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_n*fibre(i,3,Iel)*fibre(j,3,Iel)
        ENDDO
      ENDDO GAUSS_PTS1
      CALL INVERT2(F0def,F0inv,ndim)

      ! Calculate the standard
      ! displacement component
      Fdef = zero;
      DO I = 1,ndim; DO J = 1,ndim
        Fdef(I,J) = Fdef(I,J) + DOT_PRODUCT(dNi_u(J,:,Ig),auxm(I,:));
        IF(I==J) Fdef(I,J) = Fdef(I,J) + one;
      ENDDO; ENDDO

      ! Calculate the stress, stiffness
      ! and strain measures/derivatives
      Fedef = MATMUL(Fdef,F0Inv)
      CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,dNi,dE(:,:,Ig),d2E(:,:,:,Ig),Map(1:ndim,:),ndim,ndofU,nodU,nst)
      CALL INVERT2(Cdef,CdefInv(:,:,Ig),ndim); 
      CALL MATMOD_NOH(C_ijkl(:,:,Ig),S_ij(:,Ig),Matprops,Fdef,ndim,nst,nprop)

      ! Calculate the volumetric stress
      ! and the Bulk stiffness at the
      ! Gauss-Ptns
      DO m = 1,ndim
        S_kk(Ig) = S_kk(Ig) + S_ij(m,Ig);
        DO n = 1,ndim
          C_kkpp(Ig) = C_kkpp(Ig) + C_ijkl(m,n,Ig)
        ENDDO
      ENDDO

      !Calculating the pressure
      !at the Gauss points
      press(Ig) = DOT_PRODUCT(Ni_p(:,Ig),pxm)
    ENDDO

    !-----
    ! Integrate the residuals
    ! and the Jacobians
    !-----
    ! UU -Block
    ! Pure displacement block
    !
    DO M = 1,ndofU
      GAUSS_PTS2: DO Ig = 1,nip
        DO s = 1,nst; CALL VOIGHT_ITERATOR(s,I,J,nst)
          PK2 = S_ij(s,Ig) + press(Ig)*CdefInv(I,J,Ig) - S_kk(Ig)*Kdelta(I,J)/three;
          Rm(M,Iel) = Rm(M,Iel) + dE(s,M,Ig)*PK2*det(Ig)*weights(Ig);
          Km(M,1:ndofU,Iel) = Km(M,1:ndofU,Iel) + PK2*d2E(s,:,M,Ig)*det(Ig)*weights(Ig);
          DO t = 1,nst; CALL VOIGHT_ITERATOR(t,K,L,nst)
            Ctang = C_ijkl(s,t,Ig) - two* press(Ig)*CdefInv(K,I,Ig)*CdefInv(J,L,Ig)
            Ctang = Ctang - C_kkpp(Ig)*Kdelta(I,J)*Kdelta(K,L)/three;
            Km(M,1:ndofU,Iel) = Km(M,1:ndofU,Iel) + dE(s,M,Ig)*Ctang*dE(t,:,Ig)*det(Ig)*weights(Ig);
          ENDDO
        ENDDO
      END DO GAUSS_PTS2
    ENDDO

    ! UP-PU-PP-Block(s)
    ! Pressure-displacement mixed block
    !
    DO M = (ndofU+1),(ndofU+ndofP)
      GAUSS_PTS2: DO Ig = 1,nip
        S_bulk = S_kk(Ig)/C_kkpp(Ig);
        K_bulk = C_kkpp(Ig)/three;
        Rm(M,Iel) = Rm(M,Iel) + (S_bulk - press(Ig)/K_bulk)*Ni_p(M,Ig)*det(Ig)*weights(Ig);


        ! UP-PU Blocks
        DO s = 1,nst; CALL VOIGHT_ITERATOR(s,I,J,nst)
          Km(1:ndofU,M,Iel) = Km(1:ndofU,M,Iel) + CdefInv(I,J,Ig)*dE(s,:,Ig)*Ni_p(M,Ig)*det(Ig)*weights(Ig);
          Km(M,1:ndofU,Iel) = Km(M,1:ndofU,Iel) + CdefInv(I,J,Ig)*dE(s,:,Ig)*Ni_p(M,Ig)*det(Ig)*weights(Ig);
        ENDDO

        ! PP-Block
        p = (ndofU+1);
        q = (ndofU+ndofP);
        Km(p:q,M,Iel) = Km(p:q,M,Iel) - (one/K_bulk)*Ni_p(M,Ig)*Ni_p(:,Ig)*det(Ig)*weights(Ig);
      END DO GAUSS_PTS2
    ENDDO
  ENDDO ELEMENTS
  RETURN
END SUBROUTINE STATIC_SOLIDUPAS_ELEMENT

!-----------------
! Green-Lagrange-strain derivatives Active strain Variant
!-----------------
SUBROUTINE GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,dNi,dE,d2E,Map,ndim,ntots,nod,nst)
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
ENDSUBROUTINE GREENLAGRANGE_DERIVATIVES_AS


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
   Y     = Matprops(1); !1.00_iwp    !
   nu    = Matprops(2); !0.49999_iwp !

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
      Map(I,J) = (J-1)*ndim + I;
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE SDF_MAP
!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Solids_UPAS
