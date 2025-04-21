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
  REAL(iwp),INTENT(IN)   :: astrain(nodU), fibre(ndim,ndim);
  REAL(iwp),INTENT(IN)   :: utemp(ntots), MATPROP(nprop);
  REAL(iwp),INTENT(IN)   :: coord(nodU,ndim), weights(nip)
  REAL(iwp),INTENT(IN)   :: dNi_p(ndim,nodP,nip), Ni_p(nodP,nip)
  REAL(iwp),INTENT(IN)   :: dNi_u(ndim,nodU,nip), Ni_u(nodU,nip)
  REAL(iwp),INTENT(INOUT):: Km(ntots,ntots), Rm(ntots);
  REAL(iwp),PARAMETER    :: zero = 0._iwp, one = 1._iwp;


  !Residual and jacobian components
  REAL(iwp):: R1(ndofU), R2(ndofP), J3;
  REAL(iwp):: J11(ndofU,ndofU), J12(ndofU,ndofP);
  REAL(iwp):: J21(ndofP,ndofU), J22(ndofP,ndofP);

  !Deformation related terms
  REAL(iwp):: C_ijkl(nst,nst,nip1), S_ij(nst)
  REAL(iwp):: pressure, strain, gamma_f, gamma_n, gamma_s;
  REAL(iwp):: Fdef(ndim,ndim), Fedef(ndim,ndim), F0def(ndim,ndim);
  REAL(iwp):: Cdef(ndim,ndim), CdefInv(ndim,ndim), F0inv(ndim,ndim);
  REAL(iwp):: dE(nst,ndofU,nip), d2E(nst,ndofU,ndofU,nip);

  !Shape function and derivatives
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndimdddddddd,ndim), det(nip);
  REAL(iwp):: Y, nu, mu, lmbda, LJ3;


xCoord(ndim,nod,nel_pp)
dNi_p(ndim,nodP,nip), Ni_p(nodP,nip)
dNi_u(ndim,nodU,nip), Ni_u(nodU,nip)


  Km     = zero;    Rm   = zero;
  C_ijkl = zero;    S_ij = zero;    Fedef  = zero;
  dE     = zero;    d2E  = zero;


  SDF      = zero;    Jac     = zero;    Jacinv = zero;    det   = zero;
  Cdef     = zero;    CdefInv = zero;    Fdef   = zero;

!===
!
! Selective reduced integration
!
!====
ELEMENTS:DO Iel = 1, nel_pp
  !Calculate the stiffness, stress
  !and strains measure at sample
  !points
  DO Ig = 1,nip


    ! Calculate the stress, stiffness
    ! and strain measures/derivatives
    Fedef = MATMUL(Fdef,F0Inv)
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dE(:,:,Ig),d2E(:,:,:,Ig),ndim,ndofU,nodU,nst)
    CALL MATMOD_NOH(C_ijkl(:,:,Ig),S_ij(:,Ig),Matprops,Fdef,ndim,nst,nprop)
  ENDDO

  !Integrate the residual and
  !Jacobian contributions
  DO Ig = 1,nip

  ENDDO
ENDDO ELEMENTS


  GAUSS_PTS1: DO Ig=1,nip
    !-----
    !Deformation Measure Calcs
    !-----
    strain  = DOT_PRODUCT(funU,astrain)
    gamma_f = -strain;
    gamma_n =  4._iwp*gamma_f;
    gamma_s = (one/((one + gamma_f)*(one + gamma_n))) - one;
    F0def = zero;
    DO i = 1,ndim
      DO j = 1,ndim
        F0def(i,j) = Kdelta(i,j) + gamma_f*fibre(i,1)*fibre(j,1)
        IF(ndim > 1) F0def(i,j) = F0def(i,j) + gamma_s*fibre(i,2)*fibre(j,2)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_n*fibre(i,3)*fibre(j,3)
      ENDDO
    ENDDO
    CALL INVERT2(F0def, F0inv, ndim);


    Fdef  = MATMUL(derivU,auxm)
    DO i = 1,ndim; Fdef(i,i) = Fdef(i,i) + one; ENDDO;
    Fedef = MATMUL(Fdef,F0inv)
    CALL SDF_Calc(SDF, derivU, SFMAP, ndim, nodU, ndofU)
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dE,d2E,ndim,ndofU,nodU,nst)

    !-----
    !Material Model calcs
    !-----
    CALL MATERIAL_MODEL_ISO(S_ij, C_ijkl, Fedef, MATPROP, nst, ndim, nprop, material)
    
    !-----
    !Residual and jacobian calcs
    !-----
    R1  = zero;
    J11 = zero;

    J11 = J11 + MATMUL(TRANSPOSE(dE),MATMUL(C_ijkl,dE))
    DO s =1,nst
      R1  = R1  + S_ij(s)*dE(s,:)
      J11 = J11 + S_ij(s)*d2E(s,:,:)
    ENDDO

    !
    ! Residuals and Jacobian components
    !
    Rtemp = zero;
    Jtemp = zero;

    Rtemp(1:ndofU) = R1
    Jtemp(1:ndofU,1:ndofU)   = J11

    Rm = Rm + Rtemp*det*weights(Ig);
    Km = Km + Jtemp*det*weights(Ig);
  END DO GAUSS_PTS1

!
! Reduced Integration
!
  Y     = 1.00_iwp;
  nu    = 0.49999_iwp;
  mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
  lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))
  GAUSS_PTS2: DO Ig=1,nip2
    !-----
    !Discrete field interpolation functions
    !-----
    IF(standardTest.OR.Hexa27test)THEN
      Jac = MATMUL(derU,coord)
      det = determinant(Jac);
      CALL INVERT2(Jac, Jacinv, ndim);
      derivU = MATMUL(Jacinv,derU);
    ENDIF

    !-----
    !Deformation Measure Calcs
    !-----
    strain  = DOT_PRODUCT(funU,astrain)
    gamma_f = -strain;
    gamma_n =  4._iwp*gamma_f;
    gamma_s = (one/((one + gamma_f)*(one + gamma_n))) - one;
    F0def = zero;
    DO i = 1,ndim
      DO j = 1,ndim
        F0def(i,j) = Kdelta(i,j) + gamma_f*fibre(i,1)*fibre(j,1)
        IF(ndim > 1) F0def(i,j) = F0def(i,j) + gamma_s*fibre(i,2)*fibre(j,2)
        IF(ndim > 2) F0def(i,j) = F0def(i,j) + gamma_n*fibre(i,3)*fibre(j,3)
      ENDDO
    ENDDO
    CALL INVERT2(F0def, F0inv, ndim);

    Fdef  = MATMUL(derivU,auxm)
    DO i = 1,ndim; Fdef(i,i) = Fdef(i,i) + one; ENDDO;
    Fedef = MATMUL(Fdef,F0inv)
    CALL SDF_Calc(SDF, derivU, SFMAP, ndim, nodU, ndofU)
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dE,d2E,ndim,ndofU,nodU,nst)


   !===========================================
   !===========================================
   !  Integrate the pressure-displacement
   !  and pure-pressure terms Peturbed Lagrangian
   !===========================================
   !==========================================
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
	
    R2  = R2  + (DLOG(J3) - (pressure/lmbda))*funnyP(1,:)
    J22 = (-J3/lmbda)*MATMUL(TRANSPOSE(funnyP(:,:)),funnyP(:,:))


    Rm = Rm + Rtemp*det*weights2(Ig);
    Km = Km + Jtemp*det*weights2(Ig);
  END DO GAUSS_PTS2

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
   Y     = 1.00_iwp !Matprops(1);
   nu    = 0.49999_iwp !Matprops(2);

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
