MODULE Parallel_BoundaryConditions
  USE precision;      USE global_variables;
  USE maths;          USE new_library;
  USE MP_INTERFACE;   USE gather_scatter;
  CONTAINS

!-------------------------------
! Apply 3-D No-rotation BC for XY-plane with radius smoothing
!-------------------------------
SUBROUTINE SMOOTHED_NO_ROTATION_BC1(Km_mat, Rm_vec, utemp, gg_coord, gg_Constr &
                                     , centre, disp_MAP, ntots, ndim, nod, nel_pp)
  IMPLICIT NONE
  INTEGER                  :: Iel, Inode, I, J, K, L, M;
  INTEGER,   INTENT(IN)    :: ntots, ndim, nod, nel_pp;
  INTEGER,   INTENT(IN)    :: disp_MAP(nod*ndim), gg_Constr(nod*ndim,nel_pp);
  REAL(iwp), INTENT(IN)    :: centre(ndim);
  REAL(iwp), INTENT(IN)    :: gg_coord(nod,ndim,nel_pp), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT) :: Km_mat(ntots,ntots,nel_pp), Rm_vec(ntots,nel_pp);
  REAL(iwp), PARAMETER     :: zero=0._iwp, one=1._iwp, two = 2._iwp;
  REAL(iwp), PARAMETER     :: tol=1.0E-02_iwp;
  REAL(iwp)                :: cx, cy, ux, uy, penalty;
  REAL(iwp)                :: alpha, ri, r_avg, r2_accum, N_nodes;
  REAL(iwp)                :: rx, ry, Jxx, Jxy, Jyx, Jyy;
  REAL(iwp)                :: rx0, ry0, Jxx0, Jxy0, Jyx0, Jyy0;
  REAL(iwp)                :: rs_x, rs_y, rrot_x, rrot_y;
  REAL(iwp)                :: Js_xx, Js_xy, Js_yx, Js_yy;
  REAL(iwp)                :: Jrot_xx, Jrot_xy, Jrot_yx, Jrot_yy;
  REAL(iwp)                :: pi_rot, rx_rot, ry_rot;
  REAL(iwp)                :: pi_s, rx_s, ry_s;

  !
  ! Calculate the average radius
  !
  r2_accum = zero;
  N_nodes  = zero;
  M = 0;
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_Constr(L,iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        K  = ndim*(Inode-1);
        I  = disp_MAP(K+1); 
        J  = disp_MAP(K+2);
        ux = utemp(I,Iel)
        uy = utemp(J,Iel)

        r2_accum =  r2_accum + DSQRT((cx+ux)*(cx+ux) + (cy+uy)*(cy+uy));
        N_nodes  = N_nodes + one;
      ENDIF
    ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(N_nodes,N_nodes,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  CALL MPI_ALLREDUCE(r2_accum,r2_accum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  r_avg = DSQRT(r2_accum/N_nodes);

  !
  ! Apply the Boundary conditions
  !
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_Constr(L,Iel) == 0)THEN
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);

        K  = ndim*(Inode-1);
        I  = disp_MAP(K+1); 
        J  = disp_MAP(K+2);

        ux   = utemp(I,Iel);
        uy   = utemp(J,Iel);
        rx0  = Rm_vec(I,Iel);
        ry0  = Rm_vec(J,Iel);
        Jxx0 = Km_mat(I,I,Iel);
        Jxy0 = Km_mat(I,J,Iel);
        Jyx0 = Km_mat(J,I,Iel);
        Jyy0 = Km_mat(J,J,Iel);

        !
        ! No-rotation boundary condition
        !
        pi_rot  =  (cx + ux)*ux - (cy + uy)*uy;
		rx_rot  =  (cx + two*ux);
		ry_rot  = -(cy + two*uy);

        rrot_x  = rx_rot*pi_rot;
        rrot_y  = ry_rot*pi_rot;

        Jrot_xx = two*rx_rot*rx_rot + two*two*pi_rot;
        Jrot_xy = rx_rot*ry_rot;
        Jrot_yx = ry_rot*rx_rot;
        Jrot_yy = two*ry_rot*ry_rot - two*two*pi_rot;

        !
        ! Radial averaging boundary condition
        !
        ri   = DSQRT((cx+ux)*(cx+ux) + (cy+uy)*(cy+uy));
        pi_s = (r_avg - ri);
        rx_s = (cx+ux)/(r_avg*N_nodes) - (cx+ux)/(ri);
        ry_s = (cy+uy)/(r_avg*N_nodes) - (cy+uy)/(ri);

        rs_x = zero; !two*pi_s*rx_s;
        rs_y = zero; !two*pi_s*ry_s;

        Js_xx = two*pi_s*() + two*rx_s*rx_s;
        Js_xy = zero;
        Js_yx = zero;
        Js_yy = two*pi_s*() + two*rx_s*rx_s;;

        !
        ! Add on the BC penalty parameters
        !
        penalty = 7.00_iwp;
        rx  = rx0 + penalty*(rs_x + rrot_x);
        ry  = ry0 + penalty*(rs_y + rrot_y);

        Jxx = Jxx0 + penalty*(Js_xx + Jrot_xx);
        Jxy = Jxy0 + penalty*(Js_xy + Jrot_xy);
        Jyx = Jyx0 + penalty*(Js_yx + Jrot_yx);
        Jyy = Jyy0 + penalty*(Js_yy + Jrot_yy);

        Rm_vec(I,Iel)   = rx;
        Rm_vec(J,Iel)   = ry;
        Km_mat(I,I,Iel) = Jxx;
        Km_mat(I,J,Iel) = Jxy;
        Km_mat(J,I,Iel) = Jyx;
        Km_mat(J,J,Iel) = Jyy;
      ENDIF
    ENDDO
  ENDDO

  RETURN
ENDSUBROUTINE SMOOTHED_NO_ROTATION_BC1

!-------------------------------
! Detects local boundary Face elements and returns Face steering array
!-------------------------------
!Warning very slow algorithm
SUBROUTINE BoundaryFaceElement_detection(gg_Face, nloadedFace, gnum_pp, boundary_N &
                                       , boundaryID, element, ndim, nod, nFace     &
                                       , nodFace, nbnd, nel_pp, npess, numpes)
  IMPLICIT NONE
  INTEGER                      :: iel, i, j, k, l;
  INTEGER                      :: BoundaryNods, BoundaryElements;
  CHARACTER(len=15)            :: face_element
  CHARACTER(len=15), INTENT(IN):: element
  INTEGER,   INTENT(IN)        :: ndim, nod, nFace, nodFace, nbnd, nel_pp;
  INTEGER,   INTENT(IN)        :: gnum_pp(nod,nel_pp), npess, numpes;
  INTEGER,   INTENT(IN)        :: boundary_N(nbnd), boundaryID;
  INTEGER,   INTENT(INOUT)     :: gg_Face(nFace+1,nel_pp), nloadedFace;
  INTEGER, ALLOCATABLE         :: BoundaryElementPos(:), FaceNodeMap(:,:);
  INTEGER, ALLOCATABLE         :: boundary_N2(:,:), gg_pp(:,:);
  INTEGER, ALLOCATABLE         :: gFace_temp(:,:), gnum_temp(:,:);
  LOGICAL                      :: FaceTest1 = .TRUE., FaceTest2 = .TRUE.;

  !-----
  ! Initialisation
  !-----
  ALLOCATE(FaceNodeMap(nFace,nodFace), gg_pp(nod,nel_pp))
  CALL Element_FaceMASK(FaceNodeMap, element, face_element, nface, nodFace, ndim)
  gg_pp = 0;
  
  !-----
  ! Flag all local boundary nodes
  !-----
  ALLOCATE(boundary_N2(nbnd,2))
  boundary_N2(:,1) = boundary_N(:);
  boundary_N2(:,2) = 0;
  ELEMENTS: DO iel = 1,nel_pp
    CALL FIND_G3(gnum_pp(:,iel),gg_pp(:,iel),boundary_N2)
  ENDDO ELEMENTS
  DEALLOCATE(boundary_N2)

  !
  ! First, second and third Pass are reduce the search space
  !

  !-----
  ! First-Pass approximating number of Elements with boundary faces
  !-----
  BoundaryElements = 0;
  DO iel = 1,nel_pp
    BoundaryNods = 0;
    DO i = 1,nod
      IF(gg_pp(i,iel)==0) BoundaryNods = BoundaryNods + 1;
    ENDDO
    IF(BoundaryNods >= nodFace) BoundaryElements = BoundaryElements + 1;
  ENDDO


  !-----
  ! Second pass Storing approximate Elements with boundary faces
  !-----
  IF(BoundaryElements > 0)THEN
    ALLOCATE(BoundaryElementPos(BoundaryElements));
    BoundaryElementPos = 0;
    k = 0;
    DO iel = 1,nel_pp
      BoundaryNods = 0;
      DO i = 1,nod
        IF(gg_pp(i,iel)==0) BoundaryNods = BoundaryNods + 1;
      ENDDO
      IF(BoundaryNods >= nodFace) THEN
        k = k + 1;
        BoundaryElementPos(k) = iel;
      ENDIF
    ENDDO
  ENDIF

  !-----
  ! Third Pass finding the element faces and setting them 
  !-----
  gg_Face = 0;
  IF(BoundaryElements > 0)THEN
    DO i = 1,BoundaryElements
      iel = BoundaryElementPos(i);
      DO j = 1,nFace
        FaceTest1 = .TRUE.
        DO k = 1,nodFace
          FaceTest2 = (gg_pp(FaceNodeMap(j,k),iel)==0);
          FaceTest1 = FaceTest2 .AND. FaceTest1;
        ENDDO
        IF(FaceTest1)THEN
          gg_Face(j+1,iel) = boundaryID;
          gg_Face(1,iel) = gg_Face(1,iel) + 1;
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  !-----
  ! Fourth pass calculate the number of loaded faces
  !-----
  nloadedFace = 0;
  DO iel = 1,nels_pp
    nloadedFace = nloadedFace + gg_Face(1,iel);
  ENDDO
  IF(BoundaryElements > 0) DEALLOCATE(BoundaryElementPos);
  DEALLOCATE(FaceNodeMap, gg_pp)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE BoundaryFaceElement_detection

!-------------------------------
! Calculates normal vector from tangents
!-------------------------------
SUBROUTINE NORMAL_VECTOR(normal, tangents, ndim, nvectors)
  IMPLICIT NONE
  INTEGER                 :: i, j, k, l;
  INTEGER,   INTENT(IN)   :: ndim, nvectors;
  REAL(iwp), INTENT(IN)   :: tangents(nvectors,ndim);
  REAL(iwp), INTENT(INOUT):: normal(ndim);

  !Calculate the normal vector
  normal = 0._iwp;
  SELECT CASE(ndim)
    CASE(1)
      normal(1) =  tangents(1,1);
    CASE(2)
      normal(1) =  tangents(1,2);
      normal(2) = -tangents(1,1);
    CASE(3)
      normal(1) =  tangents(1,2)*tangents(2,3) - tangents(1,3)*tangents(2,2);
      normal(2) =  tangents(1,3)*tangents(2,1) - tangents(1,1)*tangents(2,3);
      normal(3) =  tangents(1,1)*tangents(2,2) - tangents(1,2)*tangents(2,1);
    CASE DEFAULT
      WRITE(*,*) "ERROR IN NUMBER OF DIMENSIONS"
      WRITE(*,*) "only 1 to 3-D acceptable"
  END SELECT
  RETURN
END SUBROUTINE NORMAL_VECTOR

!-------------------------------
! Calculates normal vector derivative from tangents and SF-derivative
!-------------------------------
SUBROUTINE NORMAL_VECTOR_DER(NormalDer,tangents,Der,ndim,ntots,nvectors)
  IMPLICIT NONE
  INTEGER                 :: i, j, k, l;
  INTEGER,   INTENT(IN)   :: ndim, ntots, nvectors;
  REAL(iwp), INTENT(IN)   :: tangents(nvectors,ndim);
  REAL(iwp), INTENT(IN)   :: Der(nvectors,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: NormalDer(ndim,ntots);

  !Calculate the normal vector
  NormalDer = 0._iwp;
  SELECT CASE(ndim)
    CASE(1)
      NormalDer(1,:) =  tangents(1,1)*Der(1,1,:);
    CASE(2)
      NormalDer(1,:) =  tangents(1,2)*Der(1,1,:);
      NormalDer(2,:) = -tangents(1,1)*Der(1,2,:);
    CASE(3)
      NormalDer(1,:)  = NormalDer(1,:) + tangents(1,2)*Der(2,3,:);
      NormalDer(1,:)  = NormalDer(1,:) + Der(1,2,:)*tangents(2,3);
      NormalDer(1,:)  = NormalDer(1,:) - tangents(1,3)*Der(2,2,:);
      NormalDer(1,:)  = NormalDer(1,:) - Der(1,3,:)*tangents(2,2);

      NormalDer(2,:) =  NormalDer(2,:) + tangents(1,3)*Der(2,1,:);
      NormalDer(2,:) =  NormalDer(2,:) + Der(1,3,:)*tangents(2,1);
      NormalDer(2,:) =  NormalDer(2,:) - tangents(1,1)*Der(2,3,:);
      NormalDer(2,:) =  NormalDer(2,:) - Der(1,1,:)*tangents(2,3);

      NormalDer(3,:) =  NormalDer(3,:) + tangents(1,1)*Der(2,2,:);
      NormalDer(3,:) =  NormalDer(3,:) + Der(1,1,:)*tangents(2,2);
      NormalDer(3,:) =  NormalDer(3,:) - tangents(1,2)*Der(2,1,:);
      NormalDer(3,:) =  NormalDer(3,:) - Der(1,2,:)*tangents(2,1);
    CASE DEFAULT
      WRITE(*,*) "ERROR IN NUMBER OF DIMENSIONS"
      WRITE(*,*) "only 1 to 3-D acceptable"
  END SELECT
  RETURN
END SUBROUTINE NORMAL_VECTOR_DER

!-------------------------------
! Calculates normal vector to face
!-------------------------------
SUBROUTINE FACE_NORMAL_VECTOR(normal,Vcoords,Fcoords,tangents,Facefun &
                             ,nod,nodFace,ndim,nvectors)
  IMPLICIT NONE
  INTEGER                  :: i;
  INTEGER,   INTENT(IN)    :: ndim, nvectors, nod, nodFace;
  REAL(iwp), INTENT(IN)    :: tangents(nvectors,ndim), Facefun(nod);
  REAL(iwp), INTENT(IN)    :: Fcoords(nodFace,ndim), Vcoords(nod,ndim);
  REAL(iwp), INTENT(INOUT) :: normal(ndim);
  REAL(iwp), PARAMETER     :: one = 1._iwp, zero = 0._iwp;
  REAL(iwp)                :: VolCentroid(ndim),FaceCentroid(ndim)
  REAL(iwp)                :: FaceNormal(ndim), subprod

  DO I = 1,ndim
    VolCentroid(I)  = SUM(Vcoords(:,I))/REAL(nod,iwp)
    FaceCentroid(I) = DOT_PRODUCT(Facefun,Fcoords(:,I))
  ENDDO
  CALL NORMAL_VECTOR(normal, tangents, ndim, nvectors)
  FaceNormal = FaceCentroid-VolCentroid;
  subprod = DOT_PRODUCT(FaceNormal,normal)
  IF(subprod < zero) normal = -normal
  RETURN
END SUBROUTINE FACE_NORMAL_VECTOR

!-------------------------------
! Calculates normal vector to face
!-------------------------------
SUBROUTINE FACE_NORMAL_VECTOR_DER(NormalDer,Vcoords,Fcoords,Der,tangents,Facefun &
                             ,nod,nodFace,ntots,ndim,nvectors)
  IMPLICIT NONE
  INTEGER                 :: i;
  INTEGER,   INTENT(IN)   :: ndim, nvectors, nod, nodFace, ntots;
  REAL(iwp), INTENT(IN)   :: tangents(nvectors,ndim), Facefun(nod);
  REAL(iwp), INTENT(IN)   :: Fcoords(nodFace,ndim), Vcoords(nod,ndim);
  REAL(iwp), INTENT(IN)   :: Der(nvectors,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: NormalDer(ndim,ntots);
  REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;
  REAL(iwp)               :: VolCentroid(ndim),FaceCentroid(ndim)
  REAL(iwp)               :: FaceNormal(ndim), normal(ndim), subprod

  DO I = 1,ndim
    VolCentroid(I)  = SUM(Vcoords(:,I))/REAL(nod,iwp)
    FaceCentroid(I) = DOT_PRODUCT(Facefun,Fcoords(:,I))
  ENDDO
  CALL NORMAL_VECTOR_DER(NormalDer,tangents,Der,ndim,ntots,nvectors)
  CALL NORMAL_VECTOR(normal, tangents, ndim, nvectors)
  FaceNormal = FaceCentroid-VolCentroid;
  subprod = DOT_PRODUCT(FaceNormal,normal)
  IF(subprod < zero) NormalDer = -NormalDer;
  RETURN
END SUBROUTINE FACE_NORMAL_VECTOR_DER

!-------------------------------
! Returns the Element face/boundary element mapping
!-------------------------------
SUBROUTINE Element_FaceMASK(EF_MAP, element, face_element, nface, nod_face, ndim)
  IMPLICIT NONE
  INTEGER                          :: i, j;
  INTEGER,           INTENT(IN)    :: nface, nod_face, ndim
  CHARACTER(len=15), INTENT(IN)    :: element
  CHARACTER(len=15), INTENT(INOUT) :: face_element
  INTEGER,           INTENT(OUT)   :: EF_MAP(nface,nod_face)

  SELECT CASE(ndim)
    CASE(2)
      face_element = 'line';
      SELECT CASE(element)
        CASE('quadrilateral')
          SELECT CASE(nod_face)
            CASE(2)
              EF_MAP(1,:) = (/1,2/);
              EF_MAP(2,:) = (/2,3/);
              EF_MAP(3,:) = (/3,4/);
              EF_MAP(4,:) = (/4,1/);
            CASE(3)
              EF_MAP(1,:) = (/1,2,3/);
              EF_MAP(2,:) = (/3,4,5/);
              EF_MAP(3,:) = (/5,6,7/);
              EF_MAP(4,:) = (/7,8,1/);
          END SELECT
        CASE('triangle')
          SELECT CASE(nod_face)
            CASE(2)
              EF_MAP(1,:) = (/1,2/);
              EF_MAP(2,:) = (/2,3/);
              EF_MAP(3,:) = (/3,1/);
            CASE(3)
              EF_MAP(1,:) = (/1,2,3/);
              EF_MAP(2,:) = (/3,4,5/);
              EF_MAP(3,:) = (/5,6,1/);
            CASE(4)
              EF_MAP(1,:) = (/1,2,3,4/);
              EF_MAP(2,:) = (/4,5,6,7/);
              EF_MAP(3,:) = (/7,8,9,1/);
            CASE(5)
              EF_MAP(1,:) = (/1,2,3,4,5/);
              EF_MAP(2,:) = (/5,6,7,8,9/);
              EF_MAP(3,:) = (/9,10,11,12,1/);
          END SELECT
      ENDSELECT
    CASE(3)
      SELECT CASE(element)
        CASE('hexahedron')
          face_element = 'quadrilateral';
          SELECT CASE(nod_face)
            CASE(4)
              EF_MAP(1,:) = (/1,5,6,2/);
              EF_MAP(2,:) = (/2,6,7,3/);
              EF_MAP(3,:) = (/3,7,8,4/);
              EF_MAP(4,:) = (/4,8,5,1/);
              EF_MAP(5,:) = (/1,2,3,4/);
              EF_MAP(6,:) = (/5,8,7,6/);
            CASE(8)
              EF_MAP(1,:) = (/1,17,5,13,6,18,2,9/);
              EF_MAP(2,:) = (/2,18,6,14,7,19,3,10/);
              EF_MAP(3,:) = (/3,19,7,15,8,20,4,11/);
              EF_MAP(4,:) = (/4,20,8,16,5,17,1,12/);
              EF_MAP(5,:) = (/1,9,2,10,3,11,4,12/);
              EF_MAP(6,:) = (/5,16,8,15,7,14,6,13/);
            CASE(9)
              EF_MAP(1,:) = (/1,17,5,13,6,18,2,9,23/);
              EF_MAP(2,:) = (/2,18,6,14,7,19,3,10,22/);
              EF_MAP(3,:) = (/3,19,7,15,8,20,4,11,24/);
              EF_MAP(4,:) = (/4,20,8,16,5,17,1,12,21/);
              EF_MAP(5,:) = (/1,9,2,10,3,11,4,12,25/);
              EF_MAP(6,:) = (/5,16,8,15,7,14,6,13,26/);
          END SELECT
        CASE('tetrahedron')
          face_element = 'triangle';
          SELECT CASE(nod_face)
            CASE(3)
              EF_MAP(1,:)  = (/1,2,3/);
              EF_MAP(2,:)  = (/1,4,2/);
              EF_MAP(3,:)  = (/2,4,3/);
              EF_MAP(4,:)  = (/3,4,1/);
            CASE(6)
              EF_MAP(1,:)  = (/1,5,2,6,3,7/);
              EF_MAP(2,:)  = (/1,8,4,9,2,5/);
              EF_MAP(3,:)  = (/2,9,4,10,3,6/);
              EF_MAP(4,:)  = (/3,10,4,8,1,7/);
          END SELECT
        CASE DEFAULT
          WRITE(*,*) "ERROR Invalid elements"
          EF_MAP = 1;
      ENDSELECT
  END SELECT
  RETURN
END SUBROUTINE Element_FaceMASK

!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Parallel_BoundaryConditions
