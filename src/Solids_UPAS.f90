MODULE Static_SolidsUPAS2
  USE precision;
  USE new_library;
  USE Parallel_supplementary_Maths;
  USE Parallel_BoundaryConditions;
  CONTAINS
!-------------------------------
!  Calculates ntots for a Solid mechanics displacement element
!-------------------------------
Function Static_SolidUPAS_Element_ntots(nod, ndim) RESULT(ntots)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nod, ndim
  INTEGER             :: ntots
  ntots = ndim*nod;
END Function Static_SolidUPAS_Element_ntots


Elm_Jac(:,:,:);
Elm_Res(:,:);

xCoord(ndim,nod,nel_pp)
Wi_detJ(nip, nel_pp)
dNi(ndim,nod,nip)
Ni(nod,nip)



!-------------------------------
! Integrates the Element-Jacobian-Residual (Voight-notation)
! Total Lagrangian pressure-displacement-formulation (Active Strain)
! Legal combinations of elements
! Hexa27U-Hexa20P
! Hexa27U-Hexa8P
! Tetra10U-Tetra4P
!-------------------------------
SUBROUTINE STATIC_SOLIDUPAS_ELEMENT(Km, Rm, auxm, ptemp, astrain, fibre, coord, SFMAP     &
                                  , material, MATPROP , points, weights, points2 &
                                  , weights2, nprop, ntots, ndofU, ndofP, ndim   &
                                  , nst, nip, nip2, nodU, nodp)
  IMPLICIT NONE
  INTEGER                :: igauss, i, j, k, l, p, q, m, n, o, s, t;
  INTEGER,  INTENT(IN)   :: nprop, ntots, ndim, nst;
  INTEGER,  INTENT(IN)   :: nip, nip2, nodU, nodp, ndofU, ndofP
  INTEGER,  INTENT(IN)   :: SFMAP(ndim,ndofU), material;
  REAL(iwp),INTENT(IN)   :: astrain(nodU), fibre(ndim,ndim);
  REAL(iwp),INTENT(IN)   :: auxm(nodU,ndim), ptemp(ndofP), MATPROP(nprop);
  REAL(iwp),INTENT(IN)   :: coord(nodU,ndim)
  REAL(iwp),INTENT(IN)   :: points(nip,ndim), weights(nip);
  REAL(iwp),INTENT(IN)   :: points2(nip2,ndim), weights2(nip2)
  REAL(iwp),INTENT(INOUT):: Km(ntots,ntots), Rm(ntots);
  REAL(iwp),PARAMETER    :: zero = 0._iwp, one = 1._iwp;

  !Residual and jacobian components
  REAL(iwp):: R1(ndofU), R2(ndofP), J3;
  REAL(iwp):: J11(ndofU,ndofU), J12(ndofU,ndofP);
  REAL(iwp):: J21(ndofP,ndofU), J22(ndofP,ndofP);
  REAL(iwp):: Jtemp(ntots,ntots), Rtemp(ntots);

  !Deformation related terms
  REAL(iwp):: C_tang(nst,nst), PK2(nst)
  REAL(iwp):: pressure, strain, gamma_f, gamma_n, gamma_s;
  REAL(iwp):: Fdef(ndim,ndim), Fedef(ndim,ndim), F0def(ndim,ndim);
  REAL(iwp):: Cdef(ndim,ndim), CdefInv(ndim,ndim), F0inv(ndim,ndim);
  REAL(iwp):: dEdef(nst,ndofU), d2Edef(nst,ndofU,ndofU);

  !Shape function and derivatives
  REAL(iwp):: derU(ndim,nodU), derivU(ndim,nodU), funP(nodP), funU(nodU);
  REAL(iwp):: funnyP(1,ndofP), SDF(ndim,ndim,ndofU);
  REAL(iwp):: Jac(ndim,ndim), Jacinv(ndim,ndim), det;
  REAL(iwp):: Y, nu, mu, lmbda, LJ3;

  LOGICAL  :: Hexa27test, Tetra10Test, standardTest;
  INTEGER  :: abaqustosg(nodP) !additional stuff

  Km       = zero;    Rm       = zero;    Jtemp  = zero;    Rtemp  = zero;
  C_tang   = zero;    PK2      = zero;    Fedef  = zero;
  dEdef    = zero;    d2Edef   = zero;

  funP     = zero;    funnyP   = zero;
  derU     = zero;    derivU   = zero;

  SDF      = zero;    Jac     = zero;    Jacinv = zero;    det   = zero;
  Cdef     = zero;    CdefInv = zero;    Fdef   = zero;


  Hexa27test   = ((nodU==27).AND.(ndim==3))
  Tetra10Test  = ((nodU==10).AND.(ndim==3))
  standardTest = (.NOT.(Hexa27test)).AND.(.NOT.(Tetra10Test))


!===================
! MAPPING Stuff
!===================
!IF(Hexa27test.AND.(nodp==8)) abaqustosg =(/1,5,6,2,4,8,7,3/)
IF(Hexa27test.AND.(nodp==8)) abaqustosg =(/1,2,3,4,5,6,7,8/)
IF(Hexa27test.AND.(nodp==20))abaqustosg =(/4,12,1,9,2,10,3,11,20,17,18,19,8,16,5,13,6,14,7,15/)
IF(Tetra10Test.AND.(nodp==4)) abaqustosg=(/1,2,3,4/)
!===================
! MAPPING Stuff
!===================

!===
!
! Selective reduced integration
!
!====
!
! Full Integration
!
  GAUSS_PTS1: DO igauss=1,nip
    !-----
    !Discrete field interpolation functions
    !-----
    !If using Hexa27 elements
    IF(Hexa27test) CALL SHAPE_FUN_HEX27(funU,points,igauss)
    IF(Hexa27test) CALL SHAPE_DER_HEX27(derU,points,igauss)

    !If using Tetra10 elements (special treatment of deriv and det)
    IF(Tetra10Test) CALL SHAPE_FUN_TET10(funU,points,igauss)
    IF(Tetra10Test) CALL JACOBIAN_DERIV_TET10(derivU, det, coord, points, igauss)

    !Other parafem standard elements
    IF(standardTest) CALL SHAPE_FUN(funU,points,igauss)
    IF(standardTest) CALL SHAPE_DER(derU,points,igauss);

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
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dEdef,d2Edef,ndim,ndofU,nodU,nst)

    !-----
    !Material Model calcs
    !-----
    CALL MATERIAL_MODEL_ISO(PK2, C_tang, Fedef, MATPROP, nst, ndim, nprop, material)
    
    !-----
    !Residual and jacobian calcs
    !-----
    R1  = zero;
    J11 = zero;

    J11 = J11 + MATMUL(TRANSPOSE(dEdef),MATMUL(C_tang,dEdef))
    DO s =1,nst
      R1  = R1  + PK2(s)*dEdef(s,:)
      J11 = J11 + PK2(s)*d2Edef(s,:,:)
    ENDDO

    !
    ! Residuals and Jacobian components
    !
    Rtemp = zero;
    Jtemp = zero;

    Rtemp(1:ndofU) = R1
    Jtemp(1:ndofU,1:ndofU)   = J11

    Rm = Rm + Rtemp*det*weights(igauss);
    Km = Km + Jtemp*det*weights(igauss);
  END DO GAUSS_PTS1

!
! Reduced Integration
!
  Y     = 1.00_iwp;
  nu    = 0.49999_iwp;
  mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
  lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))
  GAUSS_PTS2: DO igauss=1,nip2
    !-----
    !Discrete field interpolation functions
    !-----
    IF(Hexa27test) CALL SHAPE_FUN_HEX27(funU,points2,igauss)
    IF(Hexa27test) CALL SHAPE_DER_HEX27(derU,points2,igauss)

    !If using Tetra10 elements (special treatment of deriv and det)
    IF(Tetra10Test) CALL SHAPE_FUN_TET10(funU,points2,igauss)
    IF(Tetra10Test) CALL JACOBIAN_DERIV_TET10(derivU, det, coord, points2, igauss)

    !Other parafem standard elements
    IF(standardTest) CALL SHAPE_FUN(funU,points2,igauss)
    IF(standardTest) CALL SHAPE_DER(derU,points2,igauss);


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
    CALL GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dEdef,d2Edef,ndim,ndofU,nodU,nst)


    !-----
    !Residual and jacobian calcs
    !-----
    R1  = zero;    R2  = zero;
    J11 = zero;    J12 = zero;
    J22 = zero;


   !===========================================
   !===========================================
   !  Integrate the pressure-displacement
   !  and pure-pressure terms Peturbed Lagrangian
   !===========================================
   !==========================================
    CALL SHAPE_FUN(funP,points2,igauss);
    funnyP(1,:) = funP(abaqustosg)
    Cdef = MATMUL(TRANSPOSE(Fedef),Fedef);
    CALL INVERT2(Cdef,CdefInv,ndim);
    J3 = determinant(Fedef)
    pressure = DOT_PRODUCT(funP(abaqustosg),ptemp)
    DO s =1,nst
      CALL VOIGHT_ITERATOR(s,i,j,nst);
      R1  = R1  + pressure*CdefInv(i,j)*dEdef(s,:)
      J11 = J11 + pressure*CdefInv(i,j)*d2Edef(s,:,:)
      DO m = 1,ndofU;
        J12(m,:)=J12(m,:) + CdefInv(i,j)*funnyP(1,:)*dEdef(s,m)
        DO n = 1,ndofU; DO t =1,nst
          CALL VOIGHT_ITERATOR(t,k,l,nst);
          J11(m,n)=J11(m,n)-2._iwp*pressure*CdefInv(i,k)*CdefInv(j,l) &
                                  *dEdef(s,m)*dEdef(t,n)
        ENDDO; ENDDO
      ENDDO
    ENDDO
	
    R2  = R2  + (DLOG(J3) - (pressure/lmbda))*funnyP(1,:)
    J22 = (-J3/lmbda)*MATMUL(TRANSPOSE(funnyP(:,:)),funnyP(:,:))

    !
    ! Residuals and Jacobian components
    !
    Rtemp = zero;
    Jtemp = zero;

    Rtemp(1:ndofU) = R1
    Rtemp((ndofU+1):(ndofU+ndofP)) = R2

    Jtemp(1:ndofU,1:ndofU)                  = J11;
    Jtemp(1:ndofU,(ndofU+1):(ndofU+ndofP))  = J12;
    Jtemp((ndofU+1):(ndofU+ndofP),1:ndofU)  = TRANSPOSE(J12);
    Jtemp((ndofU+1):(ndofU+ndofP),(ndofU+1):(ndofU+ndofP)) = J22;
    Rm = Rm + Rtemp*det*weights2(igauss);
    Km = Km + Jtemp*det*weights2(igauss);
  END DO GAUSS_PTS2


  !
  ! Residual and Jacobian eliminated Pressure DOFs
  !
  DO i = (ndofU+ndofP+1),ntots
    Rm(I) = zero;
    Km(:,I) = zero;
    Km(I,:) = zero;
    Km(I,I) = one;
  ENDDO

  RETURN
END SUBROUTINE STATIC_SOLIDUPAS_ELEMENT


!-------------------------------
! Deformation Measures
!-------------------------------
SUBROUTINE DEFORMATION_MEASURES(Fdef, Edef, Cdef, ndim)
   IMPLICIT NONE
   INTEGER                 :: i;
   INTEGER,   INTENT(IN)   :: ndim
   REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim)
   REAL(iwp), INTENT(INOUT):: Edef(ndim,ndim), Cdef(ndim,ndim)
   REAl(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
   Edef = zero;
   Cdef = zero;
   Cdef = MATMUL(TRANSPOSE(Fdef),Fdef)
   Edef = 0.5_iwp*Cdef;
   DO i = 1,ndim
     Edef(i,i) = Edef(i,i) - 0.5_iwp*one;
   ENDDO
   RETURN
ENDSUBROUTINE DEFORMATION_MEASURES

!-------------------------------
! Deformation Invariant calculator
!-------------------------------
SUBROUTINE DEFORMATION_INVARIANTS(InvarE, Edef, ndim)   
   IMPLICIT NONE
   INTEGER                 :: i;
   INTEGER,   INTENT(IN)   :: ndim
   REAL(iwp), INTENT(IN)   :: Edef(ndim,ndim)
   REAL(iwp), INTENT(INOUT):: InvarE(3);
   REAl(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
   REAL(iwp)               :: Edef2(ndim,ndim);
   REAL(iwp)               :: I1, I2, I3, I1_2;
   I1 = zero;    I2 = zero;
   I3 = zero;    I1_2 = zero;
   Edef2 = MATMUL(Edef,Edef)
   DO i = 1,ndim
     I1   = I1 + Edef(i,i)
     I1_2 = I1_2 + Edef2(i,i)
   ENDDO
   I2 = 0.5*(I1**2 - I1_2)
   I3 = determinant(Edef)
   InvarE(1) = I1;
   InvarE(2) = I2;
   InvarE(3) = I3;
   RETURN
ENDSUBROUTINE DEFORMATION_INVARIANTS

!-----------------
! Green-Lagrange-strain derivatives Active strain Variant
!-----------------
SUBROUTINE GREENLAGRANGE_DERIVATIVES_AS(Fdef,F0Inv,SDF,dEdef,d2Edef,ndim,ntots &
                                        ,nod, nst)
  IMPLICIT NONE
   INTEGER                :: i, j, k, l, p, q, m, n;
  INTEGER,   INTENT(IN)   :: ndim, ntots, nod, nst;
  REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), F0Inv(ndim,ndim),SDF(ndim,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: d2Edef(nst,ntots,ntots), dEdef(nst,ntots);
  REAL(iwp)               :: a, d2temp(ntots,ntots), d1temp(ntots);

  d2Edef = 0._iwp; d2temp = 0._iwp;
  dEdef = 0._iwp;  d1temp = 0._iwp;
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
      dEdef(m,:) = dEdef(m,:) + a*d1temp
      IF(p/=q) dEdef(m,:) = dEdef(m,:) + a*d1temp
      d2Edef(m,:,:) = d2Edef(m,:,:) + a*d2temp
      IF(p/=q) d2Edef(m,:,:) = d2Edef(m,:,:) + a*d2temp
      d2temp = 0._iwp
      d1temp = 0._iwp
    ENDDO
    IF(i/=j) dEdef(m,:) = 2._iwp*dEdef(m,:);
    IF(i/=j) d2Edef(m,:,:) = 2._iwp*d2Edef(m,:,:);
  ENDDO
  RETURN
ENDSUBROUTINE GREENLAGRANGE_DERIVATIVES_AS


!-------------------------------
! Utemp to auxm
!-------------------------------
SUBROUTINE Utemp_to_auxm(auxm, utemp, ndim, nod, ntots)
   IMPLICIT NONE
   INTEGER                :: i, j, k, l;
   INTEGER,  INTENT(IN)   :: ndim, nod, ntots;
   REAL(iwp),INTENT(IN)   :: utemp(ntots);
   REAL(iwp),INTENT(INOUT):: auxm(nod,ndim);
   l = 1;
   DO i = 1,nod
      DO  j = 1,ndim
         auxm(i,j) = utemp(l)
         l = l + 1;
      ENDDO
   ENDDO
   RETURN
END SUBROUTINE Utemp_to_auxm


!-------------------------------
! Element pressure node map size (mixed U-P elements)
!-------------------------------
SUBROUTINE Pressure_ELEMENT_NODE_MAP_SIZE(nodP, nodMax, ndim)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nodMax, ndim;
  INTEGER,   INTENT(INOUT):: nodP;
  nodP = 1;

  !---
  !3D elements
  !---
  IF((nodMax==4 ).AND.(ndim==3)) nodP = 1;  !Tetrahedra  4-node
  IF((nodMax==10).AND.(ndim==3)) nodP = 4;  !Tetrahedra 10-node
  IF((nodMax==8 ).AND.(ndim==3)) nodP = 1;  !Hexahedra   8-node
  IF((nodMax==14).AND.(ndim==3)) nodP = 6;  !Hexahedra  14-node
  IF((nodMax==20).AND.(ndim==3)) nodP = 8;  !Hexahedra  20-node
  IF((nodMax==27).AND.(ndim==3)) nodP = 8;  !Hexahedra  27-node

  !---
  !2D elements
  !---
  IF((nodMax==3 ).AND.(ndim==2)) nodP = 1; !Triangle  3-node
  IF((nodMax==6 ).AND.(ndim==2)) nodP = 3; !Triangle  6-node
  IF((nodMax==15).AND.(ndim==2)) nodP = 6; !Triangle 15-node
  IF((nodMax==4 ).AND.(ndim==2)) nodP = 1; !Quadrilateral 4-node
  IF((nodMax==8 ).AND.(ndim==2)) nodP = 4; !Quadrilateral 8-node
  IF((nodMax==9 ).AND.(ndim==2)) nodP = 4; !Quadrilateral 9-node

  !---
  !1D elements
  !---
  IF((nodMax==2 ).AND.(ndim==1)) nodP = 1; !Line  2-node
  IF((nodMax==3 ).AND.(ndim==1)) nodP = 2; !Line  3-node
  IF((nodMax==4 ).AND.(ndim==1)) nodP = 2; !Line  4-node
  IF((nodMax==5 ).AND.(ndim==1)) nodP = 3; !Line  5-node

  RETURN
END SUBROUTINE Pressure_ELEMENT_NODE_MAP_SIZE


!-------------------------------
! Element pressure node map (mixed U-P elements)
!-------------------------------
SUBROUTINE Pressure_ELEMENT_NODE_MAP(pressure_map, nodP, nodMax, ndim)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nodP, nodMax, ndim
  INTEGER,   INTENT(INOUT):: pressure_map(nodP);
  pressure_map = 1;

  !---
  !3D elements
  !---
  !Point source 1-node
  IF((nodP==1).AND.(ndim==3)) pressure_map = 1;
  !Tetrahedra 4-node
  IF((nodP==1).AND.(ndim==3)) pressure_map = (/1,2,3,4/);
  !Hexahedra 8-node coupled to 14-node hexahedra
  IF((nodP==8).AND.(nodMax==14).AND.(ndim==3)) pressure_map=(/1,2,3,4,10,11,12,13/);
  !Hexahedra 8-node coupled to  20-node or 27-Node Hexahedra
  IF((nodP==8).AND.(ndim==3)) pressure_map=(/1,2,3,4,5,6,7,8/);
  !Hexahedra 20-node
  IF((nodP==20).AND.(ndim==3)) pressure_map=(/1,2,3,4,5,6,7,8,9,10,11,12 & 
                                               ,13,14,15,16,17,18,19,20/);
  !---
  !2D elements
  !---
  !Triangle  3-node
  IF((nodMax==3 ).AND.(ndim==2)) pressure_map = 1;
  !Triangle  6-node
  IF((nodMax==6 ).AND.(ndim==2)) pressure_map = (/1,3,5/);
  !Triangle 15-node
  IF((nodMax==15).AND.(ndim==2)) pressure_map = (/1,3,5,7,9,11/);
  !Quadrilateral 4-node
  IF((nodMax==4 ).AND.(ndim==2)) pressure_map = 1;
  !Quadrilateral 8-node
  IF((nodMax==8 ).AND.(ndim==2)) pressure_map = (/1,3,5,7/);
  !Quadrilateral 9-node
  IF((nodMax==9 ).AND.(ndim==2)) pressure_map = (/1,3,5,7/);

  !---
  !1D elements
  !---
  IF((nodMax==2 ).AND.(ndim==1)) pressure_map = 1;       !Line  2-node
  IF((nodMax==3 ).AND.(ndim==1)) pressure_map=(/1,3/);   !Line  3-node
  IF((nodMax==4 ).AND.(ndim==1)) pressure_map=(/1,4/);   !Line  4-node
  IF((nodMax==5 ).AND.(ndim==1)) pressure_map=(/1,3,5/); !Line  5-node
  RETURN
END SUBROUTINE Pressure_ELEMENT_NODE_MAP

!-------------------------------
! Shape function map generation
!-------------------------------
SUBROUTINE SHAPE_FUNCTION_MAP(SFMAP, ntots, nod, nodof)
  IMPLICIT NONE
  INTEGER               :: i, j, k, l, m, n;
  INTEGER, INTENT(IN)   :: ntots, nod, nodof;
  INTEGER, INTENT(INOUT):: SFMAP(nodof,ntots);
  INTEGER               :: SF_MAP(ntots,nodof);

  SF_MAP = 0;
  DO I = 1,nodof; L = 0;
    DO J = 1,nodof
      DO K = 1,nod
        L = L + 1;
        IF(I==J) SF_MAP(L,I) = K;
      ENDDO
    ENDDO
  ENDDO
  SF_MAP = SF_MAP + 1;
  SFMAP  = TRANSPOSE(SF_MAP)
  RETURN
ENDSUBROUTINE SHAPE_FUNCTION_MAP

!-------------------------------
! Face Shape function map generation
!-------------------------------
SUBROUTINE FACE_SHAPE_FUNCTION_MAP(SFMAP, NodalFaceMask, nFace, nodFace, nod, nodof)
  IMPLICIT NONE
  INTEGER               :: i, j, k, l, m, n;
  INTEGER, INTENT(IN)   :: nod, nodof, nFace, nodFace;
  INTEGER, INTENT(IN)   :: NodalFaceMask(nFace,nodFace);
  INTEGER, INTENT(INOUT):: SFMAP(nFace,nodof,nodof*nod);

  SFMAP = 0;
  DO i = 1,nFace
    DO j = 1,nodof
      DO k = 1,nodFace
        DO l = 1,nodof
           m = (k-1)*nodof + l;
           n = (NodalFaceMask(i,k) - 1)*nodof+l;
           SFMAP(i,j,n) = m;
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE FACE_SHAPE_FUNCTION_MAP


!-------------------------------
! ShapeFunction derivative and funny
!-------------------------------
SUBROUTINE SDF_Calc(SDF, deriv, SFMAP, ndim, nod, ntots)
   IMPLICIT NONE
   INTEGER                 :: i, j, k, l;
   INTEGER,   INTENT(IN)   :: ndim, nod, ntots
   INTEGER,   INTENT(IN)   :: SFMAP(ndim,ntots)
   REAL(iwp), INTENT(IN)   :: deriv(ndim,nod);
   REAL(iwp), INTENT(INOUT):: SDF(ndim,ndim,ntots);
   REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
   REAL(iwp)               :: deriv2(ndim,nod+1)

   deriv2 = zero; deriv2(:,2:) = deriv;
   DO i = 1,ndim
     DO j = 1,ntots
       SDF(i,:,j) = deriv2(:,SFMAP(i,j));
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE SDF_Calc


!-------------------------------
! Isotropic Material Model-Selector
!-------------------------------
SUBROUTINE MATERIAL_MODEL_ISO(PK2, C_tang, Fdef, MATPROP, nst, ndim, nprop, material)
   IMPLICIT NONE
   INTEGER,   INTENT(IN)   :: nst, ndim, nprop, material;
   REAL(iwp), INTENT(IN)   :: Fdef(ndim,ndim), MATPROP(nprop);
   REAL(iwp), INTENT(INOUT):: PK2(nst), C_tang(nst,nst);
   REAL(iwp)               :: Edef(ndim,ndim), Cdef(ndim,ndim)
   REAL(iwp)               :: InvarE(3), InvarC(3), InvarF(3)

   CALL DEFORMATION_MEASURES(Fdef, Edef, Cdef, ndim)
   CALL DEFORMATION_INVARIANTS(InvarE, Edef, ndim)   
   CALL DEFORMATION_INVARIANTS(InvarC, Cdef, ndim)   
   CALL DEFORMATION_INVARIANTS(InvarF, Fdef, ndim)   

   SELECT CASE(material)
     CASE(1)
       !Saint-Venant-Kirchoff
        CALL MATMOD_SVN(C_tang,PK2,MATPROP,InvarE,Edef,ndim,nst,nprop)
     CASE(2)
       !Mooney Rivlin
        CALL MATMOD_MRN(C_tang,PK2,MATPROP,InvarE,Edef,ndim,nst,nprop)
     CASE(3)
       !Neo-hookean
        CALL MATMOD_NOH(C_tang,PK2,MATPROP,InvarF,Cdef,ndim,nst,nprop)
     CASE DEFAULT
       !Neo-hookean
       CALL MATMOD_NOH(C_tang,PK2,MATPROP,InvarF,Cdef,ndim,nst,nprop)
       WRITE(*,*) "ERROR INVALID MATERIAL CHOICE"
       WRITE(*,*) "Default to neo-hookean"
   ENDSELECT
   RETURN
ENDSUBROUTINE MATERIAL_MODEL_ISO


!-------------------------------
! Saint-Venant-Kirchoff
!-------------------------------
SUBROUTINE MATMOD_SVN(C_tang,PK2,Matprops,InvarE,Edef,ndim,nst,nprop)
   IMPLICIT NONE
   INTEGER                 :: i, j, k, l, m, n;
   INTEGER,   INTENT(IN)   :: ndim, nst, nprop
   REAL(iwp), INTENT(IN)   :: Edef(ndim,ndim), InvarE(3), Matprops(nprop)
   REAL(iwp), INTENT(INOUT):: C_tang(nst,nst), PK2(nst);
   REAL(iwp)               :: I1, mu, lmbda, Y, nu;
   I1    = InvarE(1);
   Y     = Matprops(1);
   nu    = Matprops(2);
   mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
   lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))

   DO m = 1,nst
     CALL VOIGHT_ITERATOR(m, i, j, nst)
     PK2(m) = lmbda*I1*Kdelta(i,j) + 2*mu*Edef(i,j);
     DO n = 1,nst
       CALL VOIGHT_ITERATOR(n, k, l, nst)
       C_tang(m,n) = lmbda*Kdelta(i,j)*Kdelta(k,l) + 2*mu*Kdelta(i,k)*Kdelta(j,l);
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE MATMOD_SVN


!-------------------------------
! Mooney-Rivlin
!-------------------------------
SUBROUTINE MATMOD_MRN(C_tang,PK2,Matprops,InvarE,Edef,ndim,nst,nprop)
   IMPLICIT NONE
   INTEGER                 :: i, j, k, l, m, n;
   INTEGER,   INTENT(IN)   :: ndim, nst, nprop
   REAL(iwp), INTENT(IN)   :: Edef(ndim,ndim), InvarE(3), Matprops(nprop)
   REAL(iwp), INTENT(INOUT):: C_tang(nst,nst), PK2(nst);
   REAL(iwp)               :: I1, mu, lmbda, Y, nu;
   I1    = InvarE(1);
   Y     = Matprops(1);
   nu    = Matprops(2);
   mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
   lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))

   DO m = 1,nst
     CALL VOIGHT_ITERATOR(m, i, j, nst)
     PK2(m) = lmbda*I1*Kdelta(i,j) + 2*mu*Edef(i,j);
     DO n = 1,nst
       CALL VOIGHT_ITERATOR(n, k, l, nst)
       C_tang(m,n) = lmbda*Kdelta(i,j)*Kdelta(k,l) + 2*mu*Kdelta(i,k)*Kdelta(j,l);
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE MATMOD_MRN


!-------------------------------
! Neo-Hookean
!-------------------------------
SUBROUTINE MATMOD_NOH(C_tang,PK2,Matprops,InvarF,Cdef,ndim,nst,nprop)
   IMPLICIT NONE
   INTEGER                 :: i, j, k, l, m, n;
   INTEGER,   INTENT(IN)   :: ndim, nst, nprop
   REAL(iwp), INTENT(IN)   :: Cdef(ndim,ndim), InvarF(3), Matprops(nprop)
   REAL(iwp), INTENT(INOUT):: C_tang(nst,nst), PK2(nst);
   REAL(iwp)               :: LJ3, J3, mu0, mu, lmbda, C_inv(ndim,ndim);
   REAL(iwp)               :: Y, nu;
   J3    = InvarF(3);
   LJ3   = DLOG(J3)
   Y     = 1.00_iwp !Matprops(1);
   nu    = 0.49999_iwp !Matprops(2);

   mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
   lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))
   CALL INVERT2(Cdef,C_inv,ndim)
   DO m = 1,nst
     CALL VOIGHT_ITERATOR(m, i, j, nst)
     PK2(m) = mu*(Kdelta(i,j) - C_inv(i,j));
     DO n = 1,nst
       CALL VOIGHT_ITERATOR(n, k, l, nst)
       C_tang(m,n) =  2._iwp*mu*C_inv(i,k)*C_inv(j,l);
     ENDDO
   ENDDO
   RETURN
ENDSUBROUTINE MATMOD_NOH

!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Static_SolidsUPAS2
