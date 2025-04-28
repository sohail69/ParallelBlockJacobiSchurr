MODULE meshEnrichment
USE TensorElement;
INTEGER, PARAMETER :: iwp = SELECTED_REAL_KIND(15,300)
CONTAINS

!/***************************************\
! Node numberings for tensor
! element (Line, quad, hex etc..)
!
! The Nodal ordering in terms of entities
! is:
! Node_no = (Vertices, Edges, Faces, Volms)
!
!\***************************************/
SUBROUTINE NODE_NUMBERINGS(Node_no, pOrder, nod, ndim)
  IMPLICIT NONE
  INTEGER               :: I, J, K, L;
  INTEGER, INTENT(IN)   :: ndim, pOrder, nod;
  INTEGER, INTENT(INOUT):: Node_no(ndim,nod) !gives the map of 1D-SF to ND elm
  INTEGER, ALLOCATABLE  :: Verts(:,:), Edges(:,:), Faces(:,:), Volms(:,:);
  INTEGER               :: nVertsQ, nEdgesQ, nFacesQ, nVolmsQ;
  INTEGER               :: nVerts, nEdges, nFaces, nVolms;

  ! Calculate the number of entities
  ! for a reference quadratic tensor
  ! element upto 3-Dimensions
  nVertsQ = 2**ndim
  nEdgesQ = ndim*2**(ndim-1)
  nFacesQ = (ndim-1)*ndim*2**(ndim-3)
  nVolmsQ = ((ndim-2)*(ndim-1)*ndim*2**(ndim-4))/3;

  nVerts = nVertsQ;
  nEdges = nEdgesQ*(pOrder-1);
  nFaces = nFacesQ*(pOrder-1)*(pOrder-1);
  nVolms = nVolmsQ*(pOrder-1)*(pOrder-1)*(pOrder-1);

  ! Calculate the vertex numbering
  IF(nVerts /= 0)THEN
    ALLOCATE(Verts(ndim,nVerts))
    CALL VertexModes(Verts, nVerts, ndim)
    I = 1;
    J = nVerts;
    Node_no(:,I:J) = Verts(:,:)
    DEALLOCATE(Verts)
  ENDIF

  ! Calculate the edge numbering
  IF(nEdges /= 0)THEN
    ALLOCATE(Edges(ndim,nEdges))

    I = nVerts + 1;
    J = nVerts + nEdges;
    Node_no(:,I:J) = Edges(:,:)
    DEALLOCATE(Edges)
  ENDIF

  ! Calculate the face numbering
  IF(nFaces /= 0)THEN
    ALLOCATE(Faces(ndim,nFaces))

    I = nVerts + nEdges + 1;
    J = nVerts + nEdges + nFaces;
    Node_no(:,I:J) = Faces(:,:)
    DEALLOCATE(Faces)
  ENDIF

  ! Calculate the vertex numbering
  IF(nVolms /= 0)THEN
    ALLOCATE(Volms(ndim,nVolms))
    CALL VolmModes(Volms, nVolms, pOrder, ndim)
    I = nVerts + nEdges + nFaces + 1;
    J = nVerts + nEdges + nFaces + nVolms;
    Node_no(:,I:J) = Volms(:,:)
    DEALLOCATE(Volms)
  ENDIF
  RETURN
END SUBROUTINE NODE_NUMBERINGS


!/***************************************\
!Vertex mode numberings
!
!\***************************************/
SUBROUTINE VertexModes(Verts, nVerts, ndim)
  INTEGER, INTENT(IN)   :: nVerts, ndim
  INTEGER, INTENT(INOUT):: Verts(ndim,nVerts);
  INTEGER               :: Ivert, I, J, K;
  REAL(iwp)             :: IvertR;
  INTEGER, PARAMETER    :: modes(4) = (/1,2,2,1/);

  DO Ivert = 1,nVerts
    DO I = 1,ndim
      IvertR = REAL(Ivert,iwp)/(2._iwp**(I-1))
      K = CEILING(IvertR);
      J = MOD(K - 1, 4) + 1;
      Verts(I,Ivert) = modes(J);
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE VertexModes

!/***************************************\
!Edge mode numberings
!
!\***************************************/
SUBROUTINE EdgeModes(Edges, nEdgesQ, nEdges, pOrder, ndim)
  INTEGER, INTENT(IN)   :: nEdges, nEdgesQ, pOrder, ndim
  INTEGER, INTENT(INOUT):: Edges(ndim,nEdges);
  INTEGER               :: Iedge, Idir, IedgeQ, I, J, K;
  REAL(iwp)             :: IedgeR;
  INTEGER               :: cNodesI(ndim,nEdgesQ);

  !Can be 1, 2 or zero
  SELECT CASE(ndim)
    CASE(1)
      cNodesI(1,:) = (/0/);
    CASE(2)
      cNodesI(1,:) = (/1,0,1,0/);
      cNodesI(2,:) = (/0,1,0,1/);
    CASE(3)
      cNodesI(1,:) = (/1,2,1,2,3,3,3,3,1,2,1,2/);
      cNodesI(2,:) = (/1,2,1,2,3,3,3,3,1,2,1,2/);
      cNodesI(3,:) = (/1,2,1,2,3,3,3,3,1,2,1,2/);
    CASE DEFAULT
      WRITE(*,*) "Error higher than 3D edges not supported"
  END SELECT

  DO Iedge = 1,nEdgesQ
    K = 
    L = (Iedge-1)*(pOrder-1)
    DO I = 1,(ndim-1)
      Edges()
    ENDDO
    DO J = 1,(pOrder-1)
      Edges(K,L+J) = J + 2;
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE EdgeModes

!/***************************************\
!Face mode numberings
!
!\***************************************/
SUBROUTINE FaceModes(Faces, nFaces, pOrder, ndim)
  INTEGER, INTENT(IN)   :: nFaces, pOrder, ndim
  INTEGER, INTENT(INOUT):: Faces(ndim,nFaces);
  INTEGER               :: Iface, I, J, K;
  REAL(iwp)             :: IfaceR;
  INTEGER, PARAMETER    :: modes(4) = (/1,2,2,1/);

  nFacesQ = (ndim-1)*ndim*2**(ndim-3)
  nFaces = nFacesQ*(pOrder-1)*(pOrder-1);
  DO
   Iface = 1,nFaces
    DO I = 1,ndim-1
      DO


    ENDDO
  ENDDO
  RETURN
END SUBROUTINE FaceModes

!/***************************************\
!Volume mode numberings
!
!\***************************************/
SUBROUTINE VolmModes(Volms, nVolms, pOrder, ndim)
  INTEGER, INTENT(IN)   :: nVolms, pOrder, ndim
  INTEGER, INTENT(INOUT):: Volms(ndim,nVolms);
  INTEGER               :: Ivolm, Iters(ndim);
  INTEGER, PARAMETER    :: modes(4) = (/1,2,2,1/);

  DO Ivolm = 1,nVolms
    CALL InverseIterator(Ivolm,Iters,pOrder-1,nDIM)
    Volms(:,Ivolm) = Iters;
  ENDDO
  RETURN
END SUBROUTINE VolmModes
!/***************************************\
!\***************************************/
END MODULE meshEnrichment


























