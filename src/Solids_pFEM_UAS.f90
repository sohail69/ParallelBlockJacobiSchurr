MODULE Static_SolidsUAS
  CONTAINS

!-------------------------------
! Integrates the Element-Jacobian-Residual (Voight-notation)
! Total Lagrangian displacement-formulation (Active Strain)
!-------------------------------
SUBROUTINE STATIC_SOLIDUAS_ELEMENT(Km, Rm, utemp, gama, fibre, coord   &
                                  , MATPROP, Ni, dNi, weights  &
                                  , nprop, ntots, ndim, nst, nip, nod)
  IMPLICIT NONE
  INTEGER                 :: Ig, I, J, K, M, N;
  INTEGER,   INTENT(IN)   :: nprop, ntots, ndim, nip, nod, nst;
  REAL(iwp), INTENT(IN)   :: utemp(ntots), MATPROP(nprop);
  REAL(iwp), INTENT(IN)   :: coord(nod,ndim), weights(nip);
  REAL(iwp), INTENT(INOUT):: Km(ntots,ntots), Rm(ntots);
  REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
  REAL(iwp):: C_ijkl(nst,nst), S_ij(nst), Fdef(ndim,ndim), Fedef(ndim,ndim);
  REAL(iwp):: dE(nst,ntots), d2E(nst,ntots,ntots), auxm(nod,ndim);
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndim,ndim), det(nip);
  
  Km     = zero;    Rm     = zero; 
  C_ijkl = zero;    S_ij   = zero;
  Fdef   = zero;    Fedef  = zero;    dE  = zero;    d2E = zero;    
  Jac    = zero;    Jacinv = zero;    det = zero;

  !-----
  ! Calculate the stress,
  ! tangent-stiffness and
  ! deformation measure at
  ! the gauss-ptns
  !-----
  GAUSS_PTS1: DO Ig=1,nip
    strain = DOT_PRODUCT(Ni(:,Ig),gama)
    gamma_f =  -strain;
    gamma_n =  4._iwp*gamma_f;
    gamma_s = (one/((one + gamma_f)*(one + gamma_n))) - one;
    DO i = 1,ndim
      DO j = 1,ndim
        F0def(i,j) = Kdelta(i,j) + gamma_f*fibre(i,1)*fibre(j,1)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_s*fibre(i,2)*fibre(j,2)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_n*fibre(i,3)*fibre(j,3)
      ENDDO
    ENDDO
    CALL INVERT2(F0def,F0inv,ndim)

    Fdef  = TRANSPOSE(MATMUL(deriv,auxm))
    DO i = 1,ndim; Fdef(i,i) = Fdef(i,i) + one; ENDDO
    Fedef = MATMUL(Fdef,F0Inv)
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,dNi,dE,d2E,ndim,ntots,nod,nst)
    CALL MATMOD_NOH(C_ijkl,S_ij,Matprops,Fedef,ndim,nst,nprop)
  ENDDO GAUSS_PTS1

  !-----
  ! Integrate the residuals
  ! and the Jacobians
  !-----
  DO M = 1,ntots
    GAUSS_PTS2: DO Ig = 1,nip
      DO I = 1,nst
        Rm(M) = Rm(M) + dE(I,M),S_ij(I)*det(Ig)*weights(Ig);
        Km(M,:) = Km(M,:) + S_ij(i)*d2E(i,:,M);
        DO J = 1,nst
          Km = Km + dE(I,M)*C_ijkl(I,J)*dE(J,:)*det(Ig)*weights(Ig);
        ENDDO
      ENDDO
    END DO GAUSS_PTS2
  ENDDO

  RETURN
END SUBROUTINE STATIC_SOLIDUAS_ELEMENT

!-----------------
! Green-Lagrange-strain derivatives Active strain Variant
!-----------------
SUBROUTINE GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dE,d2E,ndim,ntots &
                                        ,nod, nst)
  IMPLICIT NONE
   INTEGER                :: i, j, k, l, p, q, m, n;
  INTEGER,   INTENT(IN)   :: ndim, ntots, nod, nst;
  REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), F0Inv(ndim,ndim),SDF(ndim,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: d2E(nst,ntots,ntots), dE(nst,ntots);
  REAL(iwp)               :: a, d2temp(ntots,ntots), d1temp(ntots);

  d2E = 0._iwp; d2temp = 0._iwp;
  dE = 0._iwp;  d1temp = 0._iwp;
  DO m = 1,nst; CALL VOIGHT_ITERATOR(m, i, j, nst)
    DO n = 1,nst; CALL VOIGHT_ITERATOR(n, p, q, nst)
      DO k = 1,ndim
        d1temp = d1temp + (SDF(k,p,:)*Fdef(k,q) + SDF(k,q,:)*Fdef(k,p))

          d2temp(l,:) = d2temp(l,:) + (SDF(k,p,l)*SDF(k,q,:)+SDF(k,q,l)*SDF(k,p,:))
      ENDDO
      a = 0.5_iwp*F0Inv(p,i)*F0Inv(q,j)
      IF(p/=q) a = F0Inv(p,i)*F0Inv(q,j)
      dE(m,:) = dE(m,:) + a*d1temp
      d2E(m,:,:) = d2E(m,:,:) + a*d2temp
    ENDDO
    IF(i/=j) dE(m,:) = 2._iwp*dE(m,:);
    IF(i/=j) d2E(m,:,:) = 2._iwp*d2E(m,:,:);
  ENDDO
  RETURN
ENDSUBROUTINE GREENLAGRANGE_DERIVATIVES_AS

!-------------------------------
! Takes a
!-------------------------------
PURE FUNCTION SDF(dNi, I, J, nDIM, nod) RESULT(dNij)
   INTEGER,  INTENT(IN) :: I, J, nDIM, nod
  REAL(iwp), INTENT(IN) :: dNi(ndim,nod)
  REAL(iwp)             :: dNij

  IF(I /= J) dNij = 0._iwp
  IF(I == J) dNij = 
END FUNCTION SDF

!-------------------------------
! Neo-Hookean
!-------------------------------
SUBROUTINE MATMOD_NOH(C_ijkl,S_ij,Matprops,Fdef,ndim,nst,nprop)
   IMPLICIT NONE
   INTEGER                 :: i, j, k, l, m, n;
   INTEGER,   INTENT(IN)   :: ndim, nst, nprop
   REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), Matprops(nprop)
   REAL(iwp), INTENT(INOUT):: C_ijkl(nst,nst), S_ij(nst);
   REAL(iwp)               :: LJ3, J3, mu, lmbda, Y, nu;
   REAL(iwp)               :: Cdef(ndim,ndim), C_inv(ndim,ndim);

   J3    = determinant(Fdef)
   LJ3   = DLOG(J3)
   Y     = Matprops(1);  ! 1.00_iwp 
   nu    = Matprops(2); !0.49999_iwp !
   mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
   lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))
   Cdef  = MATMUL(TRANSPOSE(Fdef),Fdef)
   CALL INVERT2(Cdef,C_inv,ndim)
   DO m = 1,nst
     CALL VOIGHT_ITERATOR(m, i, j, nst)
     S_ij(m) = lmbda*LJ3*C_inv(i,j) + mu*(Kdelta(i,j) - C_inv(i,j));
     DO n = 1,nst
       CALL VOIGHT_ITERATOR(n, k, l, nst)
       C_ijkl(m,n) = 2._iwp*mu*C_inv(i,k)*C_inv(l,j) &
                   + 2._iwp*lmbda*( 0.5_iwp*C_inv(k,l)*C_inv(i,j) - LJ3*C_inv(i,k)*C_inv(l,j))
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE MATMOD_NOH
!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Static_SolidsUAS
