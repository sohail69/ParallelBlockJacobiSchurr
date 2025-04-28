MODULE SolidPFEM_UAS
  USE TensorElement;
  CONTAINS
!-------------------------------
! Integrates the Element-Jacobian-Residual (Voight-notation)
! Total Lagrangian displacement-formulation (Active Strain)
!-------------------------------
SUBROUTINE SOLID_PFEM_UAS_ELM(Km, Rm, utemp, gama, fibre, MATPROP   &
                            , coord, Ni, dNi, weights, nprop, ntots &
                            , ndim, nst, nip, nod)
  IMPLICIT NONE
  INTEGER                 :: Ig, I, J, K, M, N;
  INTEGER,   INTENT(IN)   :: nprop, ntots, ndim, nip, nod, nst;
  REAL(iwp), INTENT(IN)   :: gama(nod), fibre(ndim,ndim), utemp(ntots), MATPROP(nprop);
  REAL(iwp), INTENT(IN)   :: Ni(nod,nip), dNi(ndim,nod,nip);
  REAL(iwp), INTENT(IN)   :: coord(nod,ndim), weights(nip);
  REAL(iwp), INTENT(INOUT):: Km(ntots,ntots), Rm(ntots);
  REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
  REAL(iwp):: C_ijkl(nst,nst,nip), S_ij(nst,nip), Fdef(ndim,ndim), Fedef(ndim,ndim);
  REAL(iwp):: dE(nst,ntots,nip), d2E(nst,ntots,ntots,nip), auxm(ndim,nod);
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndim,ndim), dNi_tmp(ndim,nod), det(nip);
  REAL(iwp):: strain, gamma_f, gamma_n,  gamma_s, F0def(ndim,ndim),F0inv(ndim,ndim)
  INTEGER  :: Map(ndim,nod);
  Km     = zero;    Rm     = zero; 
  C_ijkl = zero;    S_ij   = zero;
  Fdef   = zero;    Fedef  = zero;    dE  = zero;    d2E = zero;    
  Jac    = zero;    Jacinv = zero;    det = zero;


  !-----
  ! Calculate the stress,
  ! tangent-stiffness and
  ! deformation measures at
  ! the gauss-points
  !-----
  CALL SDF_MAP(Map, nod, ndim)
  DO I=1,ndim
    auxm(I,:) = utemp(Map(I,:));
  ENDDO

  GAUSS_PTS1: DO Ig=1,nip
    ! Calculate the coordinate
    ! tranform from ref-to-physical
    Jac = MATMUL(dNi(:,:,Ig), coord(:,:))
    CALL Invert2(Jac, Jacinv, nDIM)
    dNi_tmp = MATMUL(Jacinv,dNi(:,:,Ig))
    det(Ig) = DETERMINANT(Jac)

    ! Calculate the active strain
    ! components of the tensors
    strain = DOT_PRODUCT(Ni(:,Ig),gama(:))
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

   ! Calculate the standard
   ! displacement component

    Fdef = zero;
    DO I = 1,ndim
      DO J = 1,ndim
        Fdef(I,J) = Fdef(I,J) + DOT_PRODUCT(dNi(J,:,Ig),auxm(I,:));
        IF(I==J) Fdef(I,J) = Fdef(I,J) + one;
      ENDDO
    ENDDO

    ! Calculate the stress, stiffness
    ! and strain measures/derivatives
    Fedef = MATMUL(Fdef,F0Inv)
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,dNi(:,:,Ig),dE(:,:,Ig),d2E(:,:,:,Ig),Map,ndim,ntots,nod,nst)
    CALL MATMOD_NOH(C_ijkl(:,:,Ig),S_ij(:,Ig),MATPROP,Fedef,ndim,nst,nprop)
  ENDDO GAUSS_PTS1

  !-----
  ! Integrate the residuals
  ! and the Jacobians
  !-----
  DO M = 1,ntots
    GAUSS_PTS2: DO Ig = 1,nip
      DO I = 1,nst
        Rm(M) = Rm(M) + dE(I,M,Ig)*S_ij(I,Ig)*det(Ig)*weights(Ig);
        Km(M,:) = Km(M,:) + S_ij(I,Ig)*d2E(I,:,M,Ig)*det(Ig)*weights(Ig);
        DO J = 1,nst
          Km(M,:) = Km(M,:) + dE(I,M,Ig)*C_ijkl(I,J,Ig)*dE(J,:,Ig)*det(Ig)*weights(Ig);
        ENDDO
      ENDDO
    END DO GAUSS_PTS2
  ENDDO
  RETURN
END SUBROUTINE SOLID_PFEM_UAS_ELM

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
   REAL(iwp)               :: LJ3, J3, mu, lmbda, Y, nu;
   REAL(iwp)               :: Cdef(ndim,ndim), C_inv(ndim,ndim);

   J3    = DETERMINANT(Fdef)
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
! Kroneckecker delta function
!-------------------------------
PURE FUNCTION Kdelta(i,j) RESULT(Kd)
   INTEGER, INTENT(IN) :: i,j;
   REAL(iwp)           :: Kd;
   IF(i == j)THEN; Kd = 1._iwp;
   ELSE;           Kd = 0._iwp;
   ENDIF
END FUNCTION Kdelta

!-------------------------------
! Voight iterator to tensor notation
!-------------------------------
SUBROUTINE VOIGHT_ITERATOR(k, i, j, nst)
   IMPLICIT NONE
   INTEGER, INTENT(IN)   :: k, nst
   INTEGER, INTENT(INOUT):: i, j
   INTEGER, ALLOCATABLE  :: i_iterator(:), j_iterator(:);

   IF((nst < 1).OR.(k > nst))THEN
      WRITE(*,*) "Wrong number of stresses"
   ELSE
      ALLOCATE(i_iterator(nst), j_iterator(nst))
      SELECT CASE(nst)
      CASE(1)      !1-Dimensional
         i_iterator(:) = (/1/);
         j_iterator(:) = (/1/);
         i = i_iterator(k);
         j = j_iterator(k);
      CASE(3)       !2-Dimensional
         i_iterator(:) = (/1,2,1/);
         j_iterator(:) = (/1,2,2/);
         i = i_iterator(k);
         j = j_iterator(k);
      CASE(6)       !3-Dimensional
         i_iterator(:) = (/1,2,3,2,1,1/);
         j_iterator(:) = (/1,2,3,3,3,2/);
         i = i_iterator(k);
         j = j_iterator(k);
      CASE DEFAULT
         WRITE(*,*) "Wrong number of stresses"
      END SELECT
      DEALLOCATE(i_iterator, j_iterator)
   ENDIF
   RETURN
ENDSUBROUTINE VOIGHT_ITERATOR

!-------------------------------------
!  Determinant for small Matrices
!-------------------------------------
PURE FUNCTION DETERMINANT(Amat) RESULT(det)
  IMPLICIT NONE
  REAL(iwp), INTENT(IN) :: Amat(:,:)
  REAL(iwp)             :: det
  INTEGER               :: I
  I   = UBOUND(Amat,1)
  det = Amat(1,1)
  IF(I == 2) det = DETERMINANT2(Amat) 
  IF(I == 3) det = DETERMINANT3(Amat)
ENDFUNCTION DETERMINANT


PURE FUNCTION DETERMINANT2(Amat) RESULT(det)
  IMPLICIT NONE
  REAL(iwp), INTENT(IN) :: Amat(2,2)
  REAL(iwp)             :: det

  det = 0._iwp
  det = Amat(1,1)*Amat(2,2) - Amat(1,2)*Amat(2,1)
ENDFUNCTION DETERMINANT2

PURE FUNCTION DETERMINANT3(Amat) RESULT(det)
  IMPLICIT NONE
  REAL(iwp), INTENT(IN) :: Amat(3,3)
  REAL(iwp)             :: det

  det = 0._iwp
  det = det + Amat(1,1)*DETERMINANT2(Amat( (/2,3/), (/2,3/) ))
  det = det - Amat(1,2)*DETERMINANT2(Amat( (/2,3/), (/1,3/) ))
  det = det + Amat(1,3)*DETERMINANT2(Amat( (/2,3/), (/1,2/) ))
ENDFUNCTION DETERMINANT3

!-------------------------------------
!  Inverts small matrices using Gaussian elimination
!-------------------------------------
SUBROUTINE Invert2(MatA, MatAinv, n)
  IMPLICIT NONE
  INTEGER                 :: i, j;
  INTEGER,   INTENT(IN)   :: n
  REAL(iwp), INTENT(IN)   :: MatA(n,n);
  REAL(iwp), INTENT(INOUT):: MatAinv(n,n);
  REAL(iwp)               :: con
  REAL(iwp)               :: det;
  REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;

  MatAinv = MatA;
  DO i = 1,n
    con = MatAinv(i,i);
    MatAinv(i,i) = one;
    MatAinv(i,:) = MatAinv(i,:)/con
    DO j = 1,n
      IF(i /= j)THEN
        con = MatAinv(j,i)
        MatAinv(j,i) = zero;
        MatAinv(j,:) =  MatAinv(j,:) - MatAinv(i,:)*con;
      ENDIF
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE Invert2

!-------------------------------------
!  SDF map
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
ENDMODULE SolidPFEM_UAS
