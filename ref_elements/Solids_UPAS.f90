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
  INTEGER                :: Ig, i, j, k, l, p, q, m, n, o, s, t;
  INTEGER,  INTENT(IN)   :: nprop, ntots, ndim, nst;
  INTEGER,  INTENT(IN)   :: nip, nodU, nodp, ndofU, ndofP
  REAL(iwp),INTENT(IN)   :: astrain(nodU,nel_pp), fibre(ndim,ndim,nel_pp);
  REAL(iwp),INTENT(IN)   :: utemp(ntots,nel_pp), MATPROP(nprop);
  REAL(iwp),INTENT(IN)   :: coord(nodU,ndim,nel_pp), weights(nip)
  REAL(iwp),INTENT(IN)   :: dNi_p(ndim,nodP,nip), Ni_p(nodP,nip)
  REAL(iwp),INTENT(IN)   :: dNi_u(ndim,nodU,nip), Ni_u(nodU,nip)
  REAL(iwp),INTENT(INOUT):: Km(ntots,ntots,nel_pp), Rm(ntots,nel_pp);
  REAL(iwp),PARAMETER    :: zero = 0._iwp, one = 1._iwp;

  !Deformation related terms
  REAL(iwp):: S_kk(nip), C_kkpp(nip)
  REAL(iwp):: C_ijkl(nst,nst,nip), S_ij(nst,nip), CdefInv(ndim,ndim,nip)
  REAL(iwp):: strain, gamma_f, gamma_n, gamma_s, press(nip);
  REAL(iwp):: Fdef(ndim,ndim), Fedef(ndim,ndim), F0def(ndim,ndim);
  REAL(iwp):: Cdef(ndim,ndim), F0inv(ndim,ndim);
  REAL(iwp):: dE(nst,ndofU,nip), d2E(nst,ndofU,ndofU,nip);


  !Shape function and derivatives
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndim,ndim), det(nip);
  REAL(iwp):: Y, nu, mu, lmbda, LJ3;

  !Residual and jacobian components
  REAL(iwp):: R1(ndofU), R2(ndofP), J3;
  REAL(iwp):: J11(ndofU,ndofU), J12(ndofU,ndofP);
  REAL(iwp):: J21(ndofP,ndofU), J22(ndofP,ndofP);


  Km     = zero;  Rm      = zero;
  C_ijkl = zero;  S_ij    = zero;   Fedef  = zero;
  dE     = zero;  d2E     = zero;
  Jac    = zero;  Jacinv  = zero;   det    = zero;
  Cdef   = zero;  CdefInv = zero;   Fdef   = zero;

!===
!
! Selective reduced integration
!
!====
ELEMENTS:DO Iel = 1, nel_pp
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

  GAUSS_PTS1: DO Ig = 1,nip
    ! Calculate the coordinate
    ! tranform from ref-to-physical
    Jac = MATMUL(dNi(:,:,Ig), coord(:,:))
    CALL Invert2(Jac, Jacinv, nDIM)
    det(Ig) = DETERMINANT(Jac)

    ! Calculate the active strain
    ! components of the tensors
    strain = DOT_PRODUCT(Ni_u(:,Ig),gama(:))
    gamma_f =  -strain;
    gamma_n =  4._iwp*gamma_f;
    gamma_s = (one/((one + gamma_f)*(one + gamma_n))) - one;
    DO i = 1,ndim
      DO j = 1,ndim
        F0def(i,j) = Kdelta(i,j) + gamma_f*fibre(i,1)*fibre(j,1)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_s*fibre(i,2)*fibre(j,2)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_n*fibre(i,3)*fibre(j,3)
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
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dE(:,:,Ig),d2E(:,:,:,Ig),ndim,ndofU,nodU,nst)
    CALL INVERT2(Cdef,CdefInv(:,:,Ig),ndim); 
    CALL MATMOD_NOH(C_ijkl(:,:,Ig),S_ij(:,Ig),Matprops,Fdef,ndim,nst,nprop)
    DO m = 1,ndim
      S_kk(Ig) = S_kk(Ig) + S_ij(m,Ig);
      DO n = 1,ndim
        C_kkpp(Ig) = C_kkpp(Ig) + C_ijkl(m,n,Ig)
      ENDDO
    ENDDO

    !Calculating the pressure
    !at the Gauss points
    press(Ig) = DOT_PRODUCT()
  ENDDO

!xCoord(ndim,nod,nel_pp)
!dNi_p(ndim,nodP,nip), Ni_p(nodP,nip)
!dNi_u(ndim,nodU,nip), Ni_u(nodU,nip)
!nip, nodU, nodp, ndofU, ndofP

  !-----
  ! Integrate the residuals
  ! and the Jacobians
  !-----
  DO M = 1,ntots
    GAUSS_PTS2: DO Ig = 1,nip

      ! UU-block
      DO s = 1,nst
        CALL VOIGHT_ITERATOR(s,I,J,nst)
        Rm(M) = Rm(M) + dE(s,M,Ig)*S_ij(s,Ig)*det(Ig)*weights(Ig);
                          CdefInv(I,J,Ig)


        Km(M,:) = Km(M,:) + S_ij(s,Ig)*d2E(s,:,M,Ig)*det(Ig)*weights(Ig);
        DO t = 1,nst
          Km(M,:) = Km(M,:) + dE(s,M,Ig)*C_ijkl(s,t,Ig)*dE(t,:,Ig)*det(Ig)*weights(Ig);
        ENDDO
      ENDDO


    END DO GAUSS_PTS2
  ENDDO


ENDDO ELEMENTS



    J11 = J11 + MATMUL(TRANSPOSE(dE),MATMUL(C_ijkl,dE))
    DO s =1,nst
      R1  = R1  + S_ij(s)*dE(s,:)
      J11 = J11 + S_ij(s)*d2E(s,:,:)
    ENDDO

    CALL SHAPE_FUN(funP,points2,Ig);
    funnyP(1,:) = funP(abaqustosg)
    Cdef = MATMUL(TRANSPOSE(Fedef),Fedef);
    CALL INVERT2(Cdef,CdefInv,ndim);
    J3 = determinant(Fedef)
    pressure = DOT_PRODUCT(funP(abaqustosg),ptemp)
    DO s =1,nst
      CALL VOIGHT_ITERATOR(s,i,j,nst);
      R1  = R1  + pressure*CdefInv(i,j)*dE(s,:)
      J11 = J11 + pressure*CdefInv(i,j)*d2E(s,:,:)
      DO m = 1,ndofU;
        J12(m,:)=J12(m,:) + CdefInv(i,j)*funnyP(1,:)*dE(s,m)
        DO n = 1,ndofU; DO t =1,nst
          CALL VOIGHT_ITERATOR(t,k,l,nst);
          J11(m,n)=J11(m,n)-2._iwp*pressure*CdefInv(i,k)*CdefInv(j,l) &
                                  *dE(s,m)*dE(t,n)
        ENDDO; ENDDO
      ENDDO
    ENDDO

  Y     = 1.00_iwp;
  nu    = 0.49999_iwp;
  mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
  lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))


    R2  = R2  + (DLOG(J3) - (pressure/lmbda))*funnyP(1,:)
    J22 = (-J3/lmbda)*MATMUL(TRANSPOSE(funnyP(:,:)),funnyP(:,:))


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


  !
  ! Residual and Jacobian eliminated Pressure DOFs
  !
  DO i = (ndofU+ndofP+1),ntots
    Rm(I)   = zero;
    Km(I,:) = zero;
    Km(I,I) = one;
  ENDDO

  RETURN
END SUBROUTINE STATIC_SOLIDUPAS_ELEMENT

!-----------------
! Green-Lagrange-strain derivatives Active strain Variant
!-----------------
SUBROUTINE GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dE,d2E,ndim,ntots,nod,nst)
  IMPLICIT NONE
   INTEGER                :: i, j, k, l, p, q, m, n;
  INTEGER,   INTENT(IN)   :: ndim, ntots, nod, nst;
  REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), F0Inv(ndim,ndim),SDF(ndim,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: d2E(nst,ntots,ntots), dE(nst,ntots);
  REAL(iwp)               :: a, d2temp(ntots,ntots), d1temp(ntots);

  d2E = 0._iwp; d2temp = 0._iwp;
  dE = 0._iwp;  d1temp = 0._iwp;
  DO m = 1,nst
    CALL VOIGHT_ITERATOR(m, i, j, nst)
    DO n = 1,nst
      CALL VOIGHT_ITERATOR(n, p, q, nst)
      DO k = 1,ndim
        d1temp = d1temp + (SDF(k,p,:)*Fdef(k,q) + SDF(k,q,:)*Fdef(k,p))
        DO l = 1,ntots
          d2temp(l,:) = d2temp(l,:) + (SDF(k,p,l)*SDF(k,q,:)+SDF(k,q,l)*SDF(k,p,:))
        ENDDO
      ENDDO
      a = 0.5_iwp*F0Inv(p,i)*F0Inv(q,j)
      dE(m,:) = dE(m,:) + a*d1temp
      IF(p/=q) dE(m,:) = dE(m,:) + a*d1temp
      d2E(m,:,:) = d2E(m,:,:) + a*d2temp
      IF(p/=q) d2E(m,:,:) = d2E(m,:,:) + a*d2temp
      d2temp = 0._iwp
      d1temp = 0._iwp
    ENDDO
    IF(i/=j) dE(m,:) = 2._iwp*dE(m,:);
    IF(i/=j) d2E(m,:,:) = 2._iwp*d2E(m,:,:);
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
   DO m = 1,nst
     CALL VOIGHT_ITERATOR(m, i, j, nst)
     S_ij(m) = mu*(Kdelta(i,j) - C_inv(i,j));
     DO n = 1,nst
       CALL VOIGHT_ITERATOR(n, k, l, nst)
       C_ijkl(m,n) =  2._iwp*mu*C_inv(i,k)*C_inv(j,l);
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE MATMOD_NOH

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
ENDMODULE Solids_UPAS
