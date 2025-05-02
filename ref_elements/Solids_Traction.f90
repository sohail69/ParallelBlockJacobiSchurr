MODULE Solids_Traction
  USE precision;
  USE new_library;
  USE Parallel_supplementary_Maths;
  USE Parallel_BoundaryConditions;
  CONTAINS
!-------------------------------
!  Integrate Solid-mechanics Traction boundary Element (displacement Formulation)
!-------------------------------
SUBROUTINE SOLID_Traction_element(Residual, StoreKE, Element, gg_coord, Stress_pp, utemp  &
                                , disp_MAP, gg_Face, nloadedFace, ndim, nod, ntots, ndofU &
								, nodU, nst, nFace, nodFace, nipFace, nip, nel_pp)
  IMPLICIT NONE
  INTEGER                       :: Iel, Igauss, i, j, k, l, m, n, IFace;
  CHARACTER(LEN=15)             :: FaceElement;
  CHARACTER(LEN=15), INTENT(IN) :: Element;
  INTEGER  , INTENT(IN)         :: nloadedFace, ndim, nod, ntots, nst, nip, nel_pp;
  INTEGER  , INTENT(IN)         :: nFace, nodFace, nipFace, ndofU, nodU;
  INTEGER  , INTENT(IN)         :: gg_Face(nFace+1,nel_pp), disp_MAP(ndim,nodU);
  REAL(iwp), INTENT(IN)         :: gg_coord(nod,ndim,nel_pp);
  REAL(iwp), INTENT(IN)         :: Stress_pp(nst*nodFace,nloadedFace), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT)      :: Residual(ntots,nel_pp), StoreKE(ntots,ntots,nel_pp);
  REAL(iwp), PARAMETER   :: zero = 0._iwp, one = 1._iwp;

  REAL(iwp), ALLOCATABLE :: Ni(:,:), dNdxi(:,:,:), auxm(:,:);
  REAL(iwp), ALLOCATABLE :: der_F(:,:), der_F2(:,:);
  REAL(iwp), ALLOCATABLE :: funf(:), funf2(:), Jac_F(:,:);
  REAL(iwp), ALLOCATABLE :: Gcoords(:,:), Fcoords(:,:);
  REAL(iwp), ALLOCATABLE :: NormalDer(:,:), normal(:);
  REAL(iwp), ALLOCATABLE :: Rtemp(:), Ktemp(:,:);
  REAL(iwp), ALLOCATABLE :: pointsF(:,:), weightsF(:), Stress(:);

  INTEGER,   ALLOCATABLE :: NodalFaceMask(:,:), StressMap(:,:), SHAPEFUNDISPMASK(:,:,:)

  INTEGER                :: nvectors, nipFace2;
  LOGICAL                :: Hexa27test, Tetra10Test, standardTest;

  REAL(iwp):: x,y,z


  !-----
  ! Sample element Gauss integration points and weights
  !-----
  !Note: for boundary elements ntots=ndim*nod
  IF(ndim == 1) nvectors = 1;         !Special case
  IF(ndim  > 1) nvectors = ndim - 1;  !General case
  nipFace2 = nipFace;

  ALLOCATE(Ni(ndim,ndofU), dNdxi(ndim-1,ndim,ndofU), auxm(nodU,ndim) )
  ALLOCATE(der_F(ndim-1,nodFace), der_F2(ndim-1,nodFace+1) )
  ALLOCATE(funf(nodface), funf2(nodface+1), Jac_F(ndim-1,ndim) )
  ALLOCATE(Gcoords(nod,ndim), Fcoords(nodFace,ndim) )
  ALLOCATE(NormalDer(ndim,ndofU), normal(ndim) )
  ALLOCATE(Rtemp(ndofU), Ktemp(ndofU,ndofU) )
  ALLOCATE(pointsF(nipFace2,nvectors), weightsF(nipFace2), Stress(nst) )

  ALLOCATE(NodalFaceMask(nFace,nodFace), StressMap(nst,nodFace))
  ALLOCATE(SHAPEFUNDISPMASK(nFace,ndim,ndofU))

  CALL Element_FaceMASK(NodalFaceMask,element,FaceElement,nFace,nodFace,ndim);
  CALL SAMPLE2(FaceElement,pointsF,weightsF);

  Hexa27test   = ((nodU==27).AND.(ndim==3))
  Tetra10Test  = ((nodU==10).AND.(ndim==3))
  standardTest = (.NOT.(Hexa27test)).AND.(.NOT.(Tetra10Test))


  SHAPEFUNDISPMASK = 1;
  FACES:DO IFace=1,nFace
    ELM_NODES:DO I=1,nod
      F_NODES:DO J=1,nodFace
        IF(I==NodalFaceMask(IFace,J))THEN
          DIMS:DO K=1,ndim
            SHAPEFUNDISPMASK(IFace,K,K*I) = J+1;
          ENDDO DIMS
          EXIT F_NODES;
        ENDIF
      ENDDO F_NODES
    ENDDO ELM_NODES
  ENDDO FACES

  StressMap = 0;
  DO J = 1,nodFace
    DO K = 1,nst
      L = (J-1)*nst + K;
      StressMap(K,J) = L;
    ENDDO
  ENDDO 


  IF(nloadedFace /= 0)THEN
    !-----
    ! Integrate all the surface elements on reference element
    !-----
    L = 0;
    ELEMENTS:DO IEL = 1,nel_pp
      !--
      !Split the pressure and displacement for separate integration
      !--
      auxm(:,1) = utemp(disp_MAP(1,:),iel)
      IF(ndim >= 2) auxm(:,2) = utemp(disp_MAP(2,:),iel)
      IF(ndim >= 3) auxm(:,3) = utemp(disp_MAP(3,:),iel)

      Gcoords = gg_coord(:,:,Iel) + auxm;

      IF((gg_Face(1,iel)/=0).AND.(L/=nloadedFace))THEN
        ElementFaces:DO IFace = 1,nFace
          IF((gg_Face(IFace+1,Iel) /= 0).AND.(L<nloadedFace))THEN
            L = L + 1;
            Fcoords = Gcoords(NodalFaceMask(IFace,:),:);
            GAUSS_PTS:DO igauss = 1,nipFace2
              !-----
              !Face element
              !-----
              Rtemp = zero;  Ktemp = zero;
              funf  = zero;  der_F = zero;
              CALL SHAPE_FUN(funf,pointsF,igauss)
              CALL SHAPE_DER(der_F,pointsF,igauss)
              Jac_F  = MATMUL(der_F,Fcoords)
              funf2  = zero;
              der_F2 = zero;
              funf2(2:(nodface+1)) = funf;
              der_F2(:,2:(nodface+1)) = der_F;
              Ni = zero;
              dNdxi = zero;


              DO J=1,ndim
                Ni(J,:) = funf2(SHAPEFUNDISPMASK(IFace,J,:));
                dNdxi(:,J,:) = der_f2(:,SHAPEFUNDISPMASK(IFace,J,:));
              ENDDO


              Stress = zero;
              DO J = 1,nst !error in stress calc
                DO N=1,nodFace
                  Stress(J) = Stress(J) + funf(N)*Stress_pp(StressMap(J,N),L);
                ENDDO
              ENDDO


              CALL FACE_NORMAL_VECTOR(normal,Gcoords,Fcoords,Jac_F,funf &
                                     ,nodU,nodFace,ndim,nvectors);

              CALL FACE_NORMAL_VECTOR_DER(NormalDer,Gcoords,Fcoords,dNdxi,Jac_F,funf &
                                         ,nodU,nodFace,ndofU,ndim,nvectors);

              CALL SOLID_Traction_Residual(Rtemp,Stress,normal,Ni,ndim,nod,nst);
              CALL SOLID_Traction_Jacobian(Ktemp,Stress,NormalDer,Ni,ndim,nodU,nst);

              Residual(1:ndofU,IEL) = Residual(1:ndofU,IEL) &
                                    + Rtemp(1:ndofU)*weightsF(igauss);
              StoreKE(1:ndofU,1:ndofU,IEL) = StoreKE(1:ndofU,1:ndofU,IEL) &
                                           + Ktemp(1:ndofU,1:ndofU)*weightsF(igauss);
            ENDDO GAUSS_PTS
          ENDIF
        ENDDO ElementFaces
      ENDIF
      IF(L == nloadedFace) EXIT ELEMENTS;
    ENDDO ELEMENTS
  ENDIF

  DEALLOCATE(Ni, dNdxi, auxm, der_F, der_F2)
  DEALLOCATE(funf, funf2, Jac_F, Gcoords, Fcoords)
  DEALLOCATE(NormalDer, normal, Rtemp, Ktemp)
  DEALLOCATE(pointsF, weightsF, Stress)
  DEALLOCATE(NodalFaceMask, StressMap, SHAPEFUNDISPMASK)
  RETURN
ENDSUBROUTINE SOLID_Traction_element


SUBROUTINE SOLID_Traction_JAC_LORHEX27(StoreKE, Element, gg_coord, Stress_pp, utemp   &
                            , disp_MAP, gg_Face, nloadedFace, ndim, nod, ntots, ndofU &
							, nodU, nst, nFace, nodFace, nipFace, nip, nel_pp)
  IMPLICIT NONE
  INTEGER                       :: Iel, Igauss, i, j, k, l, m, n, IFace, ISubFace;
  CHARACTER(LEN=15)             :: FaceElement;
  CHARACTER(LEN=15), INTENT(IN) :: Element;
  INTEGER  , INTENT(IN)         :: nloadedFace, ndim, nod, ntots, nst, nip, nel_pp;
  INTEGER  , INTENT(IN)         :: nFace, nodFace, nipFace, ndofU, nodU;
  INTEGER  , INTENT(IN)         :: gg_Face(nFace+1,nel_pp), disp_MAP(ndim,nodU);
  REAL(iwp), INTENT(IN)         :: gg_coord(nod,ndim,nel_pp);
  REAL(iwp), INTENT(IN)         :: Stress_pp(nst*nodFace,nloadedFace), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT)      :: StoreKE(ntots,ntots,nel_pp);
  REAL(iwp), PARAMETER   :: zero = 0._iwp, one = 1._iwp;

  REAL(iwp), ALLOCATABLE :: Ni(:,:), dNdxi(:,:,:), auxm(:,:);
  REAL(iwp), ALLOCATABLE :: der_F(:,:), der_F2(:,:);
  REAL(iwp), ALLOCATABLE :: funf(:), funf2(:), Jac_F(:,:);
  REAL(iwp), ALLOCATABLE :: Gcoords(:,:), Fcoords(:,:);
  REAL(iwp), ALLOCATABLE :: NormalDer(:,:), normal(:);
  REAL(iwp), ALLOCATABLE :: Rtemp(:), Ktemp(:,:);
  REAL(iwp), ALLOCATABLE :: pointsF(:,:), weightsF(:), Stress(:);

  INTEGER,   ALLOCATABLE :: NodalFaceMask(:,:), StressMap(:,:), SHAPEFUNDISPMASK(:,:,:,:)
  INTEGER                :: subfaceMap(4,4);
  INTEGER                :: nvectors, nipFace2;
  LOGICAL                :: Hexa27test, Tetra10Test, standardTest;

  REAL(iwp):: x,y,z


  !-----
  ! Sample element Gauss integration points and weights
  !-----
  !Note: for boundary elements ntots=ndim*nod
  IF(ndim == 1) nvectors = 1;         !Special case
  IF(ndim  > 1) nvectors = ndim - 1;  !General case
  nipFace2 = 4;


  ALLOCATE(Ni(ndim,ndofU), dNdxi(ndim-1,ndim,ndofU), auxm(nodU,ndim) )
  ALLOCATE(der_F(ndim-1,4), der_F2(ndim-1,4+1) )
  ALLOCATE(funf(4), funf2(4+1), Jac_F(ndim-1,ndim) )
  ALLOCATE(Gcoords(nod,ndim), Fcoords(4,ndim) )
  ALLOCATE(NormalDer(ndim,ndofU), normal(ndim) )
  ALLOCATE(Rtemp(ndofU), Ktemp(ndofU,ndofU) )
  ALLOCATE(pointsF(nipFace2,nvectors), weightsF(nipFace2), Stress(nst) )

  ALLOCATE(NodalFaceMask(nFace,nodFace), StressMap(nst,nodFace))
  ALLOCATE(SHAPEFUNDISPMASK(nFace,ndim,4,ndofU))

  CALL Element_FaceMASK(NodalFaceMask,element,FaceElement,nFace,nodFace,ndim);
  CALL SAMPLE2(FaceElement,pointsF,weightsF);

  Hexa27test   = ((nodU==27).AND.(ndim==3))
  Tetra10Test  = ((nodU==10).AND.(ndim==3))
  standardTest = (.NOT.(Hexa27test)).AND.(.NOT.(Tetra10Test))


  !
  ! Lower Order subface map
  !
  subfaceMap(1,:) = (/1,2,9,8/)
  subfaceMap(2,:) = (/2,3,4,9/)
  subfaceMap(3,:) = (/9,4,5,6/)
  subfaceMap(4,:) = (/8,9,6,7/)


  SHAPEFUNDISPMASK = 1;
  FACES:DO IFace=1,nFace
    SUBFACES:DO ISubFace=1,4
      ELM_NODES:DO I=1,nod
        F_NODES:DO J=1,4 
          IF( I==NodalFaceMask(IFace,subfaceMap(ISubFace,J)) )THEN
            DIMS:DO K=1,ndim
              SHAPEFUNDISPMASK(IFace,K,ISubFace,K*I) = J+1;
            ENDDO DIMS
            EXIT F_NODES;
          ENDIF
        ENDDO F_NODES
      ENDDO ELM_NODES
    ENDDO SUBFACES
  ENDDO FACES

  StressMap = 0;
  DO J = 1,nodFace
    DO K = 1,nst
      L = (J-1)*nst + K;
      StressMap(K,J) = L;
    ENDDO
  ENDDO 


  IF(nloadedFace /= 0)THEN
    !-----
    ! Integrate all the surface elements on reference element
    !-----
    L = 0;
    ELEMENTS:DO IEL = 1,nel_pp
      !--
      !Split the pressure and displacement for separate integration
      !--
      auxm(:,1) = utemp(disp_MAP(1,:),iel)
      IF(ndim >= 2) auxm(:,2) = utemp(disp_MAP(2,:),iel)
      IF(ndim >= 3) auxm(:,3) = utemp(disp_MAP(3,:),iel)

      Gcoords = gg_coord(:,:,Iel) + auxm;

      IF((gg_Face(1,iel)/=0).AND.(L/=nloadedFace))THEN
        ElementFaces:DO IFace = 1,nFace
          IF((gg_Face(IFace+1,Iel) /= 0).AND.(L<nloadedFace))THEN
            L = L + 1;
            SUB_FACES:DO ISubFace = 1,4
              Fcoords = Gcoords(NodalFaceMask(IFace,subfaceMap(ISubFace,:)),:);
              GAUSS_PTS:DO igauss = 1,nipFace2
                !-----
                !Face element
                !-----
                Rtemp = zero;  Ktemp = zero;
                funf  = zero;  der_F = zero;
                CALL SHAPE_FUN(funf,pointsF,igauss)
                CALL SHAPE_DER(der_F,pointsF,igauss)
                Jac_F  = MATMUL(der_F,Fcoords)
                funf2  = zero;
                der_F2 = zero;
                funf2(2:(4+1)) = funf;
                der_F2(:,2:(4+1)) = der_F;
                Ni = zero;
                dNdxi = zero;


                DO J=1,ndim
                  Ni(J,:) = funf2(SHAPEFUNDISPMASK(IFace,J,ISubFace,:));
                  dNdxi(:,J,:) = der_f2(:,SHAPEFUNDISPMASK(IFace,J,ISubFace,:));
                ENDDO


                Stress = zero;
                DO J = 1,nst !error in stress calc
                  DO N=1,nodFace
                    Stress(J) = Stress(J) + funf(N)*Stress_pp(StressMap(J,N),L);
                  ENDDO
                ENDDO

                CALL FACE_NORMAL_VECTOR_DER(NormalDer,Gcoords,Fcoords,dNdxi,Jac_F,funf &
                                         ,nodU,nodFace,ndofU,ndim,nvectors);

                CALL SOLID_Traction_Jacobian(Ktemp,Stress,NormalDer,Ni,ndim,nodU,nst);
                StoreKE(1:ndofU,1:ndofU,IEL) = StoreKE(1:ndofU,1:ndofU,IEL) &
                                             + Ktemp(1:ndofU,1:ndofU)*weightsF(igauss);
              ENDDO GAUSS_PTS
            ENDDO SUB_FACES
          ENDIF
        ENDDO ElementFaces
      ENDIF
      IF(L == nloadedFace) EXIT ELEMENTS;
    ENDDO ELEMENTS
  ENDIF

  DEALLOCATE(Ni, dNdxi, auxm, der_F, der_F2)
  DEALLOCATE(funf, funf2, Jac_F, Gcoords, Fcoords)
  DEALLOCATE(NormalDer, normal, Rtemp, Ktemp)
  DEALLOCATE(pointsF, weightsF, Stress)
  DEALLOCATE(NodalFaceMask, StressMap, SHAPEFUNDISPMASK)
  RETURN
ENDSUBROUTINE SOLID_Traction_JAC_LORHEX27

!-------------------------------
!  Integrate Solid-mechanics Traction residual
!-------------------------------
SUBROUTINE  SOLID_Traction_Residual(Residual,Stress,normal,Ni,ndim,nod,nst)
  IMPLICIT NONE
  INTEGER                  :: I, J, K;
  INTEGER  , INTENT(IN)    :: ndim, nod, nst;
  REAL(iwp), INTENT(IN)    :: Stress(nst), normal(ndim), Ni(ndim,nod*ndim);
  REAL(iwp), INTENT(INOUT) :: Residual(nod*ndim);

  Residual = 0._iwp;
  DO K=1,nst
    CALL VOIGHT_ITERATOR(K,I,J,nst);
    Residual = Residual + Stress(K)*normal(I)*Ni(J,:);
    IF(I/=J) Residual = Residual + Stress(K)*normal(I)*Ni(J,:);
  ENDDO
ENDSUBROUTINE SOLID_Traction_Residual

!-------------------------------
!  Integrate Solid-mechanics Traction Jacobian
!-------------------------------
SUBROUTINE SOLID_Traction_Jacobian(Jacobian,Stress,NormalDer,Ni,ndim,nod,nst);
  IMPLICIT NONE
  INTEGER                  :: i, j, k, l, m, n, p, s, t;
  INTEGER  , INTENT(IN)    :: ndim, nod, nst;
  REAL(iwp), INTENT(IN)    :: Stress(nst),Ni(ndim,ndim*nod)
  REAL(iwp), INTENT(IN)    :: NormalDer(ndim,ndim*nod);
  REAL(iwp), INTENT(INOUT) :: Jacobian(nod*ndim,nod*ndim);

  Jacobian = 0._iwp;
  DO M=1,ndim*nod
    DO K=1,nst
      CALL VOIGHT_ITERATOR(K,I,J,nst);
      Jacobian(M,:) = Jacobian(M,:) + Stress(K)*NormalDer(J,:)*Ni(I,M);
      IF(I/=J) Jacobian(M,:) = Jacobian(M,:) + Stress(K)*NormalDer(J,:)*Ni(I,M);
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE SOLID_Traction_Jacobian

!-------------------------------
! Shape function map generation
!-------------------------------
SUBROUTINE SHAPE_FUNCTION_MAP1(SFMAP, ntots, nod, nodof)
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
  SFMAP = TRANSPOSE(SF_MAP)
  RETURN
END SUBROUTINE SHAPE_FUNCTION_MAP1

ENDMODULE Solids_Traction
