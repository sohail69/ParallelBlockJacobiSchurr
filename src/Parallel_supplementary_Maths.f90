MODULE Parallel_supplementary_Maths
  USE precision;      USE global_variables;
  USE maths;          USE gather_scatter;
  USE new_library;    USE MP_INTERFACE;
  USE OMP_LIB;
  IMPLICIT NONE
  CONTAINS
!*> Additional supplementary routines for calculation purposes
!-------------------------------------------------------------------!
!           Gather and scatter functions for multi-fields           !
!-------------------------------------------------------------------!
!Multi-field gather and scatter functions treat each DOF individually
!by assuming a 1-DOF reference field has nod*nels_pp nodal element
!freedoms and nn_pp global nodal freedoms. Dirchelet boundary conditions
!and special elements are handled by matrix elimination 

!-------------------------------
! Parafem Thread Partitioning
!-------------------------------
!
! Smallest partitionSize
!
PURE FUNCTION PARTITION_SIZE(threadID, nthreads, nsize) RESULT(psize)
  INTEGER,   INTENT(IN) :: threadID, nthreads, nsize
  INTEGER               :: psize

  psize = nsize/nthreads
END FUNCTION PARTITION_SIZE

!
! calculate remainder
!
PURE FUNCTION REMAINDER_CALC(threadID, nthreads, nsize) RESULT(tremainder)
  INTEGER,   INTENT(IN) :: threadID, nthreads, nsize
  INTEGER               :: tremainder

  tremainder = nsize - nthreads*PARTITION_SIZE(threadID, nthreads, nsize);
END FUNCTION REMAINDER_CALC

!
! First iterator
!
PURE FUNCTION ITERATOR_START(threadID, nthreads, nsize) RESULT(tnstart)
  INTEGER,   INTENT(IN) :: threadID, nthreads, nsize
  INTEGER               :: tnstart

  IF(REMAINDER_CALC(threadID, nthreads, nsize) /= 0)THEN
    IF(threadID < REMAINDER_CALC(threadID, nthreads, nsize) )THEN
      tnstart = threadID*(PARTITION_SIZE(threadID, nthreads, nsize)+1) + 1;
    ELSE
      tnstart = threadID*PARTITION_SIZE(threadID, nthreads, nsize) + REMAINDER_CALC(threadID, nthreads, nsize) + 1;
    ENDIF
  ELSE
    tnstart = threadID*PARTITION_SIZE(threadID, nthreads, nsize) + 1;
  ENDIF
END FUNCTION ITERATOR_START

!
! Last iterator
!
PURE FUNCTION ITERATOR_END(threadID, nthreads, nsize) RESULT(tnend)
  INTEGER,   INTENT(IN) :: threadID, nthreads, nsize
  INTEGER               :: tnend

  IF(REMAINDER_CALC(threadID, nthreads, nsize) /= 0)THEN
    IF(threadID < REMAINDER_CALC(threadID, nthreads, nsize) )THEN
      tnend = (threadID+1)*(PARTITION_SIZE(threadID, nthreads, nsize)+1);
    ELSE
      tnend = (threadID+1)*PARTITION_SIZE(threadID, nthreads, nsize) + REMAINDER_CALC(threadID, nthreads, nsize);
    ENDIF
  ELSE
    tnend = (threadID+1)*PARTITION_SIZE(threadID, nthreads, nsize);
  ENDIF
END FUNCTION ITERATOR_END


!-------------------------------
! Parafem vector parallel dot product
!-------------------------------
FUNCTION DOT_PRODUCT_P2(vectora, vectorb, neqs_pp)  RESULT(dot_product_p)
  REAL(iwp), INTENT(IN) :: vectora(:), vectorb(:)
  INTEGER,   INTENT(IN) :: neqs_pp
  INTEGER               :: I;                 !OpenMP stuff
  INTEGER               :: bufsize, ier       !MPI Stuff
  REAL(iwp)             :: l_product, tl_product, dot_product_p
  INTEGER               :: nthreads, threadID, tnstart, tnend; !OpenMP stuff

   l_product = 0._iwp;
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, threadID, nthreads, tnstart, tnend, tl_product)
   threadID = OMP_GET_THREAD_NUM();
   nthreads = OMP_GET_MAX_THREADS();
   tnstart  = ITERATOR_START(threadID, nthreads, neqs_pp)
   tnend    =   ITERATOR_END(threadID, nthreads, neqs_pp) 

    tl_product = 0._iwp
    DO I = tnstart,tnend
      tl_product = tl_product + vectora(I)*vectorb(I)
    ENDDO

   !$OMP CRITICAL
    l_product = l_product + tl_product;
   !$OMP END CRITICAL

   !$OMP BARRIER
   !$OMP END PARALLEL

  bufsize = 1
  dot_product_p = 0._iwp;
  CALL MPI_ALLREDUCE(l_product,dot_product_p,bufsize,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
END FUNCTION DOT_PRODUCT_P2


!-------------------------------
! Parafem element by element matrix transpose vector product (Multi-field)
!-------------------------------
SUBROUTINE PARAMATVECT(storA,x,b,pmul,qmul,MASK,nel,nn_pp,neqs_pp,ntots,nod,nodof)
   IMPLICIT NONE
   INTEGER                 :: iel, i;
   INTEGER,   INTENT(IN)   :: nel,ntots,nod,nodof,nn_pp,neqs_pp
   INTEGER,   INTENT(IN)   :: MASK(nodof,nod)
   REAL(iwp), INTENT(IN)   :: storA(ntots,ntots,nel), x(neqs_pp)
   REAL(iwp), INTENT(INOUT):: pmul(ntots,nel), qmul(ntots,nel);
   REAL(iwp), INTENT(INOUT):: b(neqs_pp)
   REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
   INTEGER                 :: nthreads, threadID, tnstart, tnend; !OpenMP stuff

   pmul = zero;
   qmul = zero;
   CALL GATHERM(x, pmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)

   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Iel, threadID, nthreads, tnstart, tnend)
   threadID = OMP_GET_THREAD_NUM();
   nthreads = OMP_GET_MAX_THREADS();
   tnstart = ITERATOR_START(threadID, nthreads, nel)
   tnend   =   ITERATOR_END(threadID, nthreads, nel) 

   DO iel = tnstart, tnend
      qmul(:,iel) = MATMUL(TRANSPOSE(storA(:,:,iel)),pmul(:,iel))
   ENDDO
   !$OMP BARRIER
   !$OMP END PARALLEL

   b = zero;
   CALL SCATTERM(b, qmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)
   RETURN
END SUBROUTINE PARAMATVECT

!-------------------------------
! Parafem element by element matrix vector product (Multi-field)
!-------------------------------
SUBROUTINE PARAMATVEC(storA,x,b,pmul,qmul,MASK,nel,nn_pp,neqs_pp,ntots,nod,nodof)
   IMPLICIT NONE
   INTEGER                 :: iel, i;
   INTEGER,   INTENT(IN)   :: nel,ntots,nod,nodof,nn_pp,neqs_pp
   INTEGER,   INTENT(IN)   :: MASK(nodof,nod)
   REAL(iwp), INTENT(IN)   :: storA(ntots,ntots,nel), x(neqs_pp)
   REAL(iwp), INTENT(INOUT):: pmul(ntots,nel), qmul(ntots,nel);
   REAL(iwp), INTENT(INOUT):: b(neqs_pp)
   REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
   INTEGER                 :: nthreads, threadID, tnstart, tnend; !OpenMP stuff

   pmul = zero;
   qmul = zero;
   CALL GATHERM(x, pmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp) 
   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Iel, threadID, nthreads, tnstart, tnend)
   threadID = OMP_GET_THREAD_NUM();
   nthreads = OMP_GET_MAX_THREADS();
   tnstart = ITERATOR_START(threadID, nthreads, nel)
   tnend   = ITERATOR_END(threadID, nthreads, nel) 

   DO iel = tnstart, tnend
      qmul(:,iel) = MATMUL(storA(:,:,iel),pmul(:,iel))
   ENDDO
   !$OMP BARRIER
   !$OMP END PARALLEL
   b = zero;
   CALL SCATTERM(b, qmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)
   RETURN
END SUBROUTINE PARAMATVEC

!-------------------------------
! Multi-field Gather operation
!-------------------------------
SUBROUTINE GATHERM(x, pmul, MASK1, ntots, nodof, nod, nel, neqs_pp, nn_pp)
  IMPLICIT NONE
  INTEGER                 :: iel, i, j, k, ndof_extra, neq, ndof;
  INTEGER,   INTENT(IN)   :: ntots, nodof, nod, nel, neqs_pp, nn_pp;
  INTEGER,   INTENT(IN)   :: MASK1(nodof,nod)
  REAL(iwp), INTENT(IN)   :: x(neqs_pp)
  REAL(iwp), INTENT(INOUT):: pmul(ntots,nel)
  REAL(iwp)               :: pmul2(nod,nel)
  pmul = 0._iwp;
  DO i = 1,nodof
    j = (i - 1)*nn_pp + 1;
    k = i*nn_pp;
    pmul2 = 0._iwp;
    CALL GATHER(x(j:k),pmul2)
    pmul(MASK1(i,:),:) = pmul2;
  ENDDO

   !Discontinuous galerkin additional terms (e.g. Bubble-functions)
   IF(neqs_pp > (nodof*nn_pp))THEN
      ndof       = nod*nodof
      ndof_extra = ntots - ndof
      DO iel = 1,nel
         j = (iel - 1)*ndof_extra + 1 + nodof*nn_pp;
         k = iel*ndof_extra + nodof*nn_pp;
         pmul((ndof+1):,iel) = x(j:k)
      ENDDO
   ENDIF
   RETURN
ENDSUBROUTINE GATHERM

!-------------------------------
! Multi-field Gather operation
!-------------------------------
SUBROUTINE SCATTERM(x, pmul, MASK1, ntots, nodof, nod, nel, neqs_pp, nn_pp)
   IMPLICIT NONE
   INTEGER                 :: iel, i, j, k, ndof_extra, neq, ndof;
   INTEGER,   INTENT(IN)   :: ntots, nodof, nod, nel, neqs_pp, nn_pp;
   INTEGER,   INTENT(IN)   :: MASK1(nodof,nod)
   REAL(iwp), INTENT(INOUT):: x(neqs_pp)
   REAL(iwp), INTENT(IN)   :: pmul(ntots,nel)
   REAL(iwp)               :: x2(nn_pp)
   x = 0._iwp;
   DO i = 1,nodof
      j = (i - 1)*nn_pp + 1;
      k = j + nn_pp - 1;
      x2 = 0._iwp;
      CALL SCATTER(x2,pmul(MASK1(i,:),:))
      x(j:k) = x2;
   ENDDO

   !Discontinuous galerkin additional terms (e.g. Bubble-functions)
   IF(neqs_pp > (nodof*nn_pp))THEN
      ndof = nod*nodof
      ndof_extra = ntots - ndof
      DO iel = 1,nel
         j = (iel - 1)*ndof_extra + 1 + nodof*nn_pp;
         k = iel*ndof_extra + nodof*nn_pp;
         x(j:k) = pmul((ndof+1):,iel)
      ENDDO
   ENDIF
   RETURN
ENDSUBROUTINE SCATTERM

!-------------------------------
! Vector Outer product of two small vectors
!-------------------------------
SUBROUTINE VEC_OUTERPRODUCT(Aij, xi, xj, m, n)
   IMPLICIT NONE
   INTEGER                 :: i, j;
   INTEGER,   INTENT(IN)   :: m, n;
   REAL(iwp), INTENT(IN)   :: xi(m), xj(n)
   REAL(iwp), INTENT(INOUT):: Aij(m,n)
   REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
   Aij = zero;
   DO i = 1,m
      DO j = 1,n
         Aij(i,j) = xi(i)*xj(j)
      ENDDO
   ENDDO
   RETURN
END SUBROUTINE VEC_OUTERPRODUCT

!-------------------------------
! Kroneckecker delta function
!-------------------------------
PURE FUNCTION Kdelta(i,j) RESULT(Kd)
   INTEGER, INTENT(IN) :: i,j;
   REAL(iwp)           :: Kd;
   IF(i == j)THEN; Kd = 1._iwp;
   ELSE;           Kd = 0._iwp;
   ENDIF
END

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


!-------------------------------
! Sample integration points in iso parametric reference frame
!-------------------------------
SUBROUTINE SAMPLE2(element, points, weights)
  IMPLICIT NONE
  CHARACTER(LEN=15), INTENT(IN):: element
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp),          PARAMETER:: zero = 0._iwp, one = 1._iwp;

  SELECT CASE(element)
    CASE('line');          CALL SAMPLE_LINE(points,weights);
    CASE('quadrilateral'); CALL SAMPLE_QUADRILATERAL(points,weights);
    CASE('tetrahedron');   CALL SAMPLE_TETAHEDRON(points,weights);
    CASE('hexahedron');    CALL SAMPLE_HEXAHEDRON(points,weights);
    CASE DEFAULT
     WRITE(*,*) "ERROR unavailable element  choice"
  ENDSELECT
  RETURN
ENDSUBROUTINE SAMPLE2

!
! Line sample points
!
SUBROUTINE SAMPLE_LINE(points,weights);
  IMPLICIT NONE
  INTEGER                      :: i, j, k, l;
  INTEGER                      :: nip, ndim;
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp),          PARAMETER:: zero = 0._iwp, one = 1._iwp;
  REAl(iwp)                    :: pa, pb, pc, wa, wb, wc; 

  nip = UBOUND(points,1)
  ndim = UBOUND(points,2)
  IF(ndim /= 1) WRITE(*,*) "ERROR Wrong dimension for line"
  IF(ndim /= 1) RETURN

  SELECT CASE (nip)
    CASE(1)
      pa = 0.0_iwp;
      wb = 2._iwp;
      points(:,1) = (/pa/);
      weights(:)  = (/wa/);
    CASE(2)
      pa = 1.0_iwp/dsqrt(3.0_iwp);
      wa = 1.0_iwp;
      points(:,1) = (/pa,-pa/);
      weights(:)  = (/wa,wa/);
    CASE(3)
      pa = 0.7745966692414834_iwp;
      pb = 0.0_iwp;
      
      wa = 0.5555555555555556_iwp;
      wb = 0.8888888888888888_iwp;
      points(:,1) = (/pa,pb,-pa/);
      weights(:)  = (/wa,wb,wa/);
    CASE(4)
      pa = 0.861136311594053_iwp;
      pb = 0.339981043584856_iwp;
      wa = 0.347854845137454_iwp;
      wb = 0.652145154862546_iwp;
      points(:,1) = (/pa,pb,-pb,-pa/);
      weights(:)  = (/wa,wb,wb,wa/);
    CASE(5)
      pa = 0.538469_iwp;
      pb = 0.90618_iwp;
      pc = zero;
      wa = 0.478629_iwp;
      wb = 0.236927_iwp;
      wc = 0.568889_iwp;
      points(:,1) = (/pa,pb,pc,-pb,-pa/);
      weights(:)  = (/wa,wb,wc,wb,wa/);
    CASE DEFAULT
      WRITE(*,*) "ERROR wrong number of integration points for Line"
      RETURN
  ENDSELECT
  RETURN
ENDSUBROUTINE SAMPLE_LINE

!
! Chebychev Line sample points
!
SUBROUTINE SAMPLE_CHEBY_LINE(points,weights);
  IMPLICIT NONE
  INTEGER                      :: i, j, k, l;
  INTEGER                      :: nip, ndim;
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp),          PARAMETER:: zero = 0._iwp, one = 1._iwp;
  REAl(iwp)                    :: pa, pb, pc, wa, wb, wc; 

  nip = UBOUND(points,1)
  ndim = UBOUND(points,2)
  IF(ndim /= 1) WRITE(*,*) "ERROR Wrong dimension for line"
  IF(ndim /= 1) RETURN

  SELECT CASE (nip)
    CASE(1)
      pa = 0.0_iwp;
      wb = 2._iwp;
      points(:,1) = (/pa/);
      weights(:)  = (/wa/);
    CASE(2)
      pa = 1.0_iwp/DSQRT(3.0_iwp);
      wa = 1.0_iwp;
      points(:,1) = (/pa,-pa/);
      weights(:)  = (/wa,wa/);
    CASE(3)
      pa = 0.7745966692414834_iwp;
      pb = 0.0_iwp;
      
      wa = 0.5555555555555556_iwp;
      wb = 0.8888888888888888_iwp;
      points(:,1) = (/pa,pb,-pa/);
      weights(:)  = (/wa,wb,wa/);
    CASE(4)
      pa = 0.861136311594053_iwp;
      pb = 0.339981043584856_iwp;
      wa = 0.347854845137454_iwp;
      wb = 0.652145154862546_iwp;
      points(:,1) = (/pa,pb,-pb,-pa/);
      weights(:)  = (/wa,wb,wb,wa/);
    CASE(5)
      pa = 0.538469_iwp;
      pb = 0.90618_iwp;
      pc = zero;
      wa = 0.478629_iwp;
      wb = 0.236927_iwp;
      wc = 0.568889_iwp;
      points(:,1) = (/pa,pb,pc,-pb,-pa/);
      weights(:)  = (/wa,wb,wc,wb,wa/);
    CASE DEFAULT
      WRITE(*,*) "ERROR wrong number of integration points for Line"
      RETURN
  ENDSELECT
  RETURN
ENDSUBROUTINE SAMPLE_CHEBY_LINE

!
! Quadrilateral sample points
!
SUBROUTINE SAMPLE_QUADRILATERAL(points,weights);
  IMPLICIT NONE
  INTEGER                      :: i, j, l, nip1D;
  INTEGER                      :: ndim, nip;
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp),          PARAMETER:: zero = 0._iwp, one = 1._iwp;
  REAl(iwp),        ALLOCATABLE:: points1D(:,:), weights1D(:);
  REAl(iwp)                    :: pa, pb, wa, wb; 

  nip = UBOUND(points,1)
  ndim = UBOUND(points,2)
  IF(ndim /= 2) WRITE(*,*) "ERROR Wrong dimension for Quadrilateral"
  IF(ndim /= 2) RETURN
  SELECT CASE (nip)
    CASE(1);  nip1D = 1;
    CASE(4);  nip1D = 2;
    CASE(9);  nip1D = 3;
    CASE(16); nip1D = 4;
    CASE DEFAULT
      WRITE(*,*) "ERROR wrong number of integration points for Quadrilateral"
      RETURN
  ENDSELECT


  ALLOCATE(points1D(nip1D,1), weights1D(nip1D));
  CALL SAMPLE_LINE(points1D,weights1D);
  l = 0;
  DO i = 1,nip1D
    DO j = 1,nip1D
      l = l + 1;
      weights(l)  = weights1D(i)*weights1D(j)
      points(l,1) = points1D(i,1)
      points(l,2) = points1D(j,1)
    ENDDO
  ENDDO
  DEALLOCATE(points1D, weights1D);
  RETURN
ENDSUBROUTINE SAMPLE_QUADRILATERAL

!
! Hexhedron sample points
!
SUBROUTINE SAMPLE_HEXAHEDRON(points,weights);
  IMPLICIT NONE
  INTEGER                      :: i, j, k, l, nip1D;
  INTEGER                      :: ndim, nip;
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp),          PARAMETER:: zero = 0._iwp, one = 1._iwp;
  REAl(iwp),        ALLOCATABLE:: points1D(:,:), weights1D(:);
  REAl(iwp)                    :: pa, pb, wa, wb, b, c; 

  nip = UBOUND(points,1)
  ndim = UBOUND(points,2)
  points = zero;
  weights = zero
  IF(ndim /= 3) WRITE(*,*) "ERROR Wrong dimension for Hexhedron"
  IF(ndim /= 3) RETURN

  SELECT CASE (nip)
    CASE(1);   nip1D = 1;
    CASE(8);   nip1D = 2;
    CASE(27);  nip1D = 3;
    CASE(64);  nip1D = 4;
    CASE(125); nip1D = 5;
    CASE DEFAULT
      WRITE(*,*) "ERROR wrong number of integration points for Hexhedron"
      RETURN
  ENDSELECT

  ALLOCATE(points1D(nip1D,1), weights1D(nip1D));
  CALL SAMPLE_LINE(points1D,weights1D);
  l = 0;
  DO i = 1,nip1D
    DO j = 1,nip1D
      DO k = 1,nip1D
        l = l + 1;
        weights(l)  = weights1D(i)*weights1D(j)*weights1D(k)
        points(l,1) = points1D(i,1)
        points(l,2) = points1D(j,1)
        points(l,3) = points1D(k,1)
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE(points1D, weights1D);
  RETURN
ENDSUBROUTINE SAMPLE_HEXAHEDRON

!Uses and alternativie quadrature rule for hexahedra
SUBROUTINE SAMPLE_HEXA_ALT(points,weights)
  IMPLICIT NONE
  INTEGER                      :: i, j, k, l;
  INTEGER                      :: ndim, nip;
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp),          PARAMETER:: zero = 0._iwp, one = 1._iwp;
  REAl(iwp)                    :: a, b;

!  a = 0.8_iwp;
!  b = 3.2_iwp;
a = 4.0_iwp/3.0_iwp;
b = 0.0_iwp;

!7-point Rule
  points(1,:) =  (/one,zero,zero/)
  points(2,:) =  (/zero,one,zero/)
  points(3,:) =  (/zero,zero,one/)
  points(4,:) = -(/one,zero,zero/)
  points(5,:) = -(/zero,one,zero/)
  points(6,:) = -(/zero,zero,one/)
  points(7,:) =  (/zero,zero,zero/)

  weights(1) = a
  weights(2) = a
  weights(3) = a
  weights(4) = a
  weights(5) = a
  weights(6) = a
  weights(7) = b
  RETURN
ENDSUBROUTINE SAMPLE_HEXA_ALT


SUBROUTINE SAMPLE_TETAHEDRON(points,weights);
  IMPLICIT NONE
  INTEGER                      :: i, j, k, l;
  INTEGER                      :: ndim, nod, nip;
  REAL(iwp),      INTENT(INOUT):: points(:,:), weights(:)
  REAl(iwp)                    :: pa, pb, pc, pd, pe, pf, pg, wa, wb, wc, wd; 

  nip = UBOUND(points,1)

  SELECT CASE(nip)
    CASE(1) !1-Point integration
      pa = 0.25_iwp 
      wa = 1.0_iwp/6.0_iwp

      points(1,:) = (/pa,pa,pa/);
      weights(:) = (/wa/);
    CASE(4) !4-Point Gauss integration
      pa = 0.58541020_iwp;
      pb = 0.13819660_iwp;
      wa = 0.041666667_iwp ;

      points(1,:) = (/pa,pb,pb/);
      points(2,:) = (/pb,pa,pb/);
      points(3,:) = (/pb,pb,pa/);
      points(4,:) = (/pb,pb,pb/);
      weights(:) = (/wa,wa,wa,wa/);

    CASE(5) !4-Point Gauss integration
      pa = 0.50000000_iwp;
      pb = 0.16666667_iwp;
      pc = 0.25000000_iwp;
      wa = 0.07500000_iwp;
      wb =-0.13333334_iwp;


      points(1,:) = (/pa,pb,pb/);
      points(2,:) = (/pb,pa,pb/);
      points(3,:) = (/pb,pb,pa/);
      points(4,:) = (/pb,pb,pb/);
      points(5,:) = (/pc,pc,pc/);
      weights(:) = (/wa,wa,wa,wa,wb/);


    CASE(8) !8-Point Gauss integration
      pa = .01583591_iwp;
      pb = .328054697_iwp;
      pc = .679143178_iwp;
      pd = .106952274_iwp;
      wa = .023087995_iwp;
      wb = .018578672_iwp;

      points(1,:) = (/pa,pb,pb/);
      points(2,:) = (/pb,pa,pb/);
      points(3,:) = (/pb,pb,pa/);
      points(4,:) = (/pb,pb,pb/);
      points(5,:) = (/pc,pd,pd/);
      points(6,:) = (/pd,pc,pd/);
      points(7,:) = (/pd,pd,pc/);
      points(8,:) = (/pd,pd,pd/);
      weights(:) = (1._iwp/6._iwp)*(/wa,wa,wa,wa,wb,wb,wb,wb/);

    CASE(15) !15-Point Gauss integration
      pa =  0.25_iwp
      pb =  0.333333333333333_iwp
      pc =  0.000000000000000_iwp
      pd =  0.090909090909091_iwp
      pe =  0.727272727272727_iwp
      pf =  0.433449846426336_iwp
      pg =  0.066550153573664_iwp

      wa =  0.030283678097089_iwp
      wb =  0.006026785714286_iwp
      wc =  0.011645249086029_iwp
      wd =  0.010949141561386_iwp

      points(1,:)  = (/pa,pa,pa/);
      points(2,:)  = (/pb,pb,pb/);
      points(3,:)  = (/pc,pb,pb/);
      points(4,:)  = (/pb,pc,pb/);
      points(5,:)  = (/pb,pb,pc/);
      points(6,:)  = (/pd,pd,pd/);
      points(7,:)  = (/pe,pd,pd/);
      points(8,:)  = (/pd,pe,pd/);
      points(9,:)  = (/pd,pd,pe/);
      points(10,:) = (/pf,pg,pg/);
      points(11,:) = (/pg,pf,pg/);
      points(12,:) = (/pg,pg,pf/);
      points(13,:) = (/pg,pf,pf/);
      points(14,:) = (/pf,pg,pf/);
      points(15,:) = (/pf,pf,pg/);
      weights(:) = (/wa,wb,wb,wb,wb,wc,wc,wc,wc,wd,wd,wd,wd,wd,wd/);
    CASE DEFAULT
      WRITE(*,*) "INCORRECT NUMBER OF INTEGRATION POITNS FOR TETRAHEDRA"
  ENDSELECT
  RETURN
ENDSUBROUTINE SAMPLE_TETAHEDRON


!-------------------------------
! Hexahedra8 Iso 27 shape-function and derivative
!-------------------------------
SUBROUTINE SHAPE_FUN_HEX8ISO27(fun,points,igauss)
  IMPLICIT NONE
  INTEGER                 :: i, j, k;
  INTEGER                 :: ndim, nod;
  INTEGER,   INTENT(IN)   :: igauss
  REAL(iwp), INTENT(IN)   :: points(:,:);
  REAL(iwp), INTENT(INOUT):: fun(:);
  REAl(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
  REAl(iwp), PARAMETER    :: half=0.5_iwp, eigth=0.125_iwp;
  REAl(iwp)               :: zeta, eta, xi;
  REAl(iwp)               :: etam, xim, zetam, etap, xip, zetap;


  ndim = UBOUND(points,2)
  nod  = UBOUND(fun,1)  
  fun  = zero;

  xi    = points(igauss,1)
  eta   = points(igauss,2)
  zeta  = points(igauss,3)

  etam  = one  - eta 
  xim   = one  - xi  
  zetam = one  - zeta
  etap  = eta  + one 
  xip   = xi   + one 
  zetap = zeta + one
         
  fun(1) = eigth*xim*etam*zetam
  fun(2) = eigth*xim*etam*zetap
  fun(3) = eigth*xip*etam*zetap
  fun(4) = eigth*xip*etam*zetam
  fun(5) = eigth*xim*etap*zetam
  fun(6) = eigth*xim*etap*zetap
  fun(7) = eigth*xip*etap*zetap
  fun(8) = eigth*xip*etap*zetam
  RETURN
ENDSUBROUTINE SHAPE_FUN_HEX8ISO27

SUBROUTINE APPLY_LOWER_ORDER_LINEAR_CONSTRAINTS_HEX27(LOR_MAT, LOR_VEC)
  IMPLICIT NONE
  INTEGER                 :: i;
  REAL(iwp), INTENT(INOUT):: LOR_MAT(27,27), LOR_VEC(27);
  REAl(iwp), PARAMETER    :: zero=0._iwp, one=1._iwp;
  REAl(iwp), PARAMETER    :: half=0.5_iwp, quarter=0.25_iwp, eigth=0.125_iwp;
  REAL(iwp)               :: UniqueMAT(8,8);

  UniqueMAT(:,:) = LOR_MAT(1:8,1:8)
  LOR_MAT = zero;

  !--
  !Implement linear constraints
  !--

  !Edge Constraints
  LOR_MAT(1,9)  = -half;   LOR_MAT(2,9)  = -half;   LOR_MAT(9,9)   = one; 
  LOR_MAT(2,10) = -half;   LOR_MAT(3,10) = -half;   LOR_MAT(10,10) = one; 
  LOR_MAT(3,11) = -half;   LOR_MAT(4,11) = -half;   LOR_MAT(11,11) = one; 
  LOR_MAT(4,12) = -half;   LOR_MAT(1,12) = -half;   LOR_MAT(12,12) = one; 

  LOR_MAT(5,13) = -half;   LOR_MAT(6,13) = -half;   LOR_MAT(13,13) = one; 
  LOR_MAT(6,14) = -half;   LOR_MAT(7,14) = -half;   LOR_MAT(14,14) = one; 
  LOR_MAT(7,15) = -half;   LOR_MAT(8,15) = -half;   LOR_MAT(15,15) = one; 
  LOR_MAT(8,16) = -half;   LOR_MAT(5,16) = -half;   LOR_MAT(16,16) = one; 

  LOR_MAT(1,17) = -half;   LOR_MAT(5,17) = -half;   LOR_MAT(17,17) = one; 
  LOR_MAT(2,18) = -half;   LOR_MAT(6,18) = -half;   LOR_MAT(18,18) = one; 
  LOR_MAT(3,19) = -half;   LOR_MAT(7,19) = -half;   LOR_MAT(19,19) = one; 
  LOR_MAT(4,20) = -half;   LOR_MAT(8,20) = -half;   LOR_MAT(20,20) = one; 


  !Face Constraints
  LOR_MAT(1,21) = -quarter;   LOR_MAT(5,21) = -quarter;   LOR_MAT(4,21) = -quarter; 
  LOR_MAT(8,21) = -quarter;   LOR_MAT(21,21) = one;

  LOR_MAT(2,22) = -quarter;   LOR_MAT(3,22) = -quarter;   LOR_MAT(6,22) = -quarter; 
  LOR_MAT(7,22) = -quarter;   LOR_MAT(22,22) = one;

  LOR_MAT(1,23) = -quarter;   LOR_MAT(2,23) = -quarter;   LOR_MAT(5,23) = -quarter; 
  LOR_MAT(6,23) = -quarter;   LOR_MAT(23,23) = one;

  LOR_MAT(3,24) = -quarter;   LOR_MAT(4,24) = -quarter;   LOR_MAT(7,24) = -quarter; 
  LOR_MAT(8,24) = -quarter;   LOR_MAT(24,24) = one;

  LOR_MAT(1,25) = -quarter;   LOR_MAT(2,25) = -quarter;   LOR_MAT(3,25) = -quarter; 
  LOR_MAT(4,25) = -quarter;   LOR_MAT(25,25) = one;

  LOR_MAT(5,26) = -quarter;   LOR_MAT(6,26) = -quarter;   LOR_MAT(7,26) = -quarter; 
  LOR_MAT(8,26) = -quarter;   LOR_MAT(26,26) = one;

  !Volume Constraints
  LOR_MAT(1,27) = -eigth;   LOR_MAT(2,27) = -eigth;   LOR_MAT(3,27)  = -eigth; 
  LOR_MAT(4,27) = -eigth;   LOR_MAT(5,27) = -eigth;   LOR_MAT(6,27)  = -eigth; 
  LOR_MAT(7,27) = -eigth;   LOR_MAT(8,27) = -eigth;   LOR_MAT(27,27) = one; 

  !Apply symmetry
  LOR_MAT = LOR_MAT + TRANSPOSE(LOR_MAT)
 
 DO I=1,27; LOR_MAT(I,I) = half*LOR_MAT(I,I); ENDDO;

  !Re-add unique components
  LOR_MAT(1:8,1:8) = LOR_MAT(1:8,1:8) + UniqueMAT
  RETURN
ENDSUBROUTINE APPLY_LOWER_ORDER_LINEAR_CONSTRAINTS_HEX27

!-------------------------------
! Hexahedra27 shape-function and derivative
!-------------------------------
SUBROUTINE SHAPE_FUN_HEX27(fun,points,igauss)
  IMPLICIT NONE
  INTEGER                 :: i, j, k;
  INTEGER                 :: ndim, nod;
  INTEGER,   INTENT(IN)   :: igauss
  REAL(iwp), INTENT(IN)   :: points(:,:);
  REAL(iwp), INTENT(INOUT):: fun(:);
  REAl(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp, half = 0.5_iwp;
  REAl(iwp)               :: zeta, eta, gama;
  REAl(iwp)               :: zetaposP1, zetapos0, zetaposM1;
  REAl(iwp)               :: etaposP1,  etapos0,  etaposM1;
  REAl(iwp)               :: gamaposP1, gamapos0, gamaposM1;

  ndim = UBOUND(points,2)
  nod  = UBOUND(fun,1)  

  eta  = points(igauss,1)
  gama = points(igauss,2)
  zeta = points(igauss,3)


  zetaposM1 =  half*zeta*(zeta - one);
  zetapos0  = (one - zeta*zeta);
  zetaposP1 =  half*zeta*(zeta + one);

  etaposM1  =  half*eta*(eta - one);
  etapos0   = (one - eta*eta);
  etaposP1  =  half*eta*(eta + one);

  gamaposM1 =  half*gama*(gama - one);
  gamapos0  = (one - gama*gama);
  gamaposP1 =  half*gama*(gama + one);


  fun(1) = zetaposM1*etaposM1*gamaposM1
  fun(2) = zetaposP1*etaposM1*gamaposM1
  fun(3) = zetaposP1*etaposP1*gamaposM1
  fun(4) = zetaposM1*etaposP1*gamaposM1

  fun(5) = zetaposM1*etaposM1*gamaposP1
  fun(6) = zetaposP1*etaposM1*gamaposP1
  fun(7) = zetaposP1*etaposP1*gamaposP1
  fun(8) = zetaposM1*etaposP1*gamaposP1

  fun(9)  = zetapos0*etaposM1*gamaposM1
  fun(10) = zetaposP1*etapos0*gamaposM1
  fun(11) = zetapos0*etaposP1*gamaposM1
  fun(12) = zetaposM1*etapos0*gamaposM1

  fun(13) = zetapos0*etaposM1*gamaposP1
  fun(14) = zetaposP1*etapos0*gamaposP1
  fun(15) = zetapos0*etaposP1*gamaposP1
  fun(16) = zetaposM1*etapos0*gamaposP1 

  fun(17) = zetaposM1*etaposM1*gamapos0
  fun(18) = zetaposP1*etaposM1*gamapos0
  fun(19) = zetaposP1*etaposP1*gamapos0
  fun(20) = zetaposM1*etaposP1*gamapos0


  fun(21) = zetaposM1*etapos0*gamapos0
  fun(22) = zetaposP1*etapos0*gamapos0
  fun(23) = zetapos0*etaposM1*gamapos0
  fun(24) = zetapos0*etaposP1*gamapos0
  fun(25) = zetapos0*etapos0*gamaposM1
  fun(26) = zetapos0*etapos0*gamaposP1

  fun(27) = zetapos0*etapos0*gamapos0
  RETURN
ENDSUBROUTINE SHAPE_FUN_HEX27


SUBROUTINE SHAPE_DER_HEX27(der,points,igauss)
  IMPLICIT NONE
  INTEGER                 :: i, j, k;
  INTEGER                 :: ndim, nod;
  INTEGER,   INTENT(IN)   :: igauss;
  REAL(iwp), INTENT(IN)   :: points(:,:);
  REAL(iwp), INTENT(INOUT):: der(:,:);
  REAl(iwp), PARAMETER    :: one = 1._iwp, two = 2._iwp, half = 0.5_iwp;
  REAl(iwp)               :: zeta, eta, gama;
  REAl(iwp)               :: zetaposP1, zetapos0, zetaposM1;
  REAl(iwp)               :: etaposP1,  etapos0,  etaposM1;
  REAl(iwp)               :: gamaposP1, gamapos0, gamaposM1;

  REAl(iwp)               :: der_zetaposP1, der_zetapos0, der_zetaposM1;
  REAl(iwp)               :: der_etaposP1,  der_etapos0,  der_etaposM1;
  REAl(iwp)               :: der_gamaposP1, der_gamapos0, der_gamaposM1;

  ndim = UBOUND(der , 1)
  nod  = UBOUND(der , 2) 

  eta  = points(igauss,1)
  gama = points(igauss,2)
  zeta = points(igauss,3)


  zetaposM1 =  half*zeta*(zeta - one);
  zetapos0  = (one - zeta*zeta);
  zetaposP1 =  half*zeta*(zeta + one);

  etaposM1  =  half*eta*(eta - one);
  etapos0   = (one - eta*eta);
  etaposP1  =  half*eta*(eta + one);

  gamaposM1 =  half*gama*(gama - one);
  gamapos0  = (one - gama*gama);
  gamaposP1 =  half*gama*(gama + one);

  der_zetaposM1 =  zeta-half;
  der_zetapos0  = -two*zeta;
  der_zetaposP1 =  zeta+half;

  der_etaposM1  =  eta-half;
  der_etapos0   = -two*eta;
  der_etaposP1  =  eta+half;

  der_gamaposM1 =  gama-half;
  der_gamapos0  = -two*gama;
  der_gamaposP1 =  gama+half;

  !-------------
  !Derivative in first direction
  !-------------
  der(1,1) = der_zetaposM1*etaposM1*gamaposM1
  der(1,2) = der_zetaposP1*etaposM1*gamaposM1
  der(1,3) = der_zetaposP1*etaposP1*gamaposM1
  der(1,4) = der_zetaposM1*etaposP1*gamaposM1

  der(1,5) = der_zetaposM1*etaposM1*gamaposP1
  der(1,6) = der_zetaposP1*etaposM1*gamaposP1
  der(1,7) = der_zetaposP1*etaposP1*gamaposP1
  der(1,8) = der_zetaposM1*etaposP1*gamaposP1

  der(1,9)  = der_zetapos0*etaposM1*gamaposM1
  der(1,10) = der_zetaposP1*etapos0*gamaposM1
  der(1,11) = der_zetapos0*etaposP1*gamaposM1
  der(1,12) = der_zetaposM1*etapos0*gamaposM1

  der(1,13) = der_zetapos0*etaposM1*gamaposP1
  der(1,14) = der_zetaposP1*etapos0*gamaposP1
  der(1,15) = der_zetapos0*etaposP1*gamaposP1
  der(1,16) = der_zetaposM1*etapos0*gamaposP1 

  der(1,17) = der_zetaposM1*etaposM1*gamapos0
  der(1,18) = der_zetaposP1*etaposM1*gamapos0
  der(1,19) = der_zetaposP1*etaposP1*gamapos0
  der(1,20) = der_zetaposM1*etaposP1*gamapos0


  der(1,21) = der_zetaposM1*etapos0*gamapos0
  der(1,22) = der_zetaposP1*etapos0*gamapos0
  der(1,23) = der_zetapos0*etaposM1*gamapos0
  der(1,24) = der_zetapos0*etaposP1*gamapos0
  der(1,25) = der_zetapos0*etapos0*gamaposM1
  der(1,26) = der_zetapos0*etapos0*gamaposP1

  der(1,27) = der_zetapos0*etapos0*gamapos0


  !-------------
  !Derivative in second direction
  !-------------
  der(2,1) = zetaposM1*der_etaposM1*gamaposM1
  der(2,2) = zetaposP1*der_etaposM1*gamaposM1
  der(2,3) = zetaposP1*der_etaposP1*gamaposM1
  der(2,4) = zetaposM1*der_etaposP1*gamaposM1

  der(2,5) = zetaposM1*der_etaposM1*gamaposP1
  der(2,6) = zetaposP1*der_etaposM1*gamaposP1
  der(2,7) = zetaposP1*der_etaposP1*gamaposP1
  der(2,8) = zetaposM1*der_etaposP1*gamaposP1

  der(2,9)  = zetapos0*der_etaposM1*gamaposM1
  der(2,10) = zetaposP1*der_etapos0*gamaposM1
  der(2,11) = zetapos0*der_etaposP1*gamaposM1
  der(2,12) = zetaposM1*der_etapos0*gamaposM1

  der(2,13) = zetapos0*der_etaposM1*gamaposP1
  der(2,14) = zetaposP1*der_etapos0*gamaposP1
  der(2,15) = zetapos0*der_etaposP1*gamaposP1
  der(2,16) = zetaposM1*der_etapos0*gamaposP1 

  der(2,17) = zetaposM1*der_etaposM1*gamapos0
  der(2,18) = zetaposP1*der_etaposM1*gamapos0
  der(2,19) = zetaposP1*der_etaposP1*gamapos0
  der(2,20) = zetaposM1*der_etaposP1*gamapos0


  der(2,21) = zetaposM1*der_etapos0*gamapos0
  der(2,22) = zetaposP1*der_etapos0*gamapos0
  der(2,23) = zetapos0*der_etaposM1*gamapos0
  der(2,24) = zetapos0*der_etaposP1*gamapos0
  der(2,25) = zetapos0*der_etapos0*gamaposM1
  der(2,26) = zetapos0*der_etapos0*gamaposP1

  der(2,27) = zetapos0*der_etapos0*gamapos0


  !-------------
  !Derivative in third direction
  !-------------
  der(3,1) = zetaposM1*etaposM1*der_gamaposM1
  der(3,2) = zetaposP1*etaposM1*der_gamaposM1
  der(3,3) = zetaposP1*etaposP1*der_gamaposM1
  der(3,4) = zetaposM1*etaposP1*der_gamaposM1

  der(3,5) = zetaposM1*etaposM1*der_gamaposP1
  der(3,6) = zetaposP1*etaposM1*der_gamaposP1
  der(3,7) = zetaposP1*etaposP1*der_gamaposP1
  der(3,8) = zetaposM1*etaposP1*der_gamaposP1

  der(3,9)  = zetapos0*etaposM1*der_gamaposM1
  der(3,10) = zetaposP1*etapos0*der_gamaposM1
  der(3,11) = zetapos0*etaposP1*der_gamaposM1
  der(3,12) = zetaposM1*etapos0*der_gamaposM1

  der(3,13) = zetapos0*etaposM1*der_gamaposP1
  der(3,14) = zetaposP1*etapos0*der_gamaposP1
  der(3,15) = zetapos0*etaposP1*der_gamaposP1
  der(3,16) = zetaposM1*etapos0*der_gamaposP1 

  der(3,17) = zetaposM1*etaposM1*der_gamapos0
  der(3,18) = zetaposP1*etaposM1*der_gamapos0
  der(3,19) = zetaposP1*etaposP1*der_gamapos0
  der(3,20) = zetaposM1*etaposP1*der_gamapos0

  der(3,21) = zetaposM1*etapos0*der_gamapos0
  der(3,22) = zetaposP1*etapos0*der_gamapos0
  der(3,23) = zetapos0*etaposM1*der_gamapos0
  der(3,24) = zetapos0*etaposP1*der_gamapos0
  der(3,25) = zetapos0*etapos0*der_gamaposM1
  der(3,26) = zetapos0*etapos0*der_gamaposP1

  der(3,27) = zetapos0*etapos0*der_gamapos0
  RETURN
ENDSUBROUTINE SHAPE_DER_HEX27


!=====================================================
! Shape function
!=====================================================
SUBROUTINE SHAPE_FUN_TET10(fun,points,i)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: i
  REAL(iwp), INTENT(IN)   :: points(:,:)
  REAL(iwp), INTENT(INOUT):: fun(:)
  REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
  REAL(iwp), PARAMETER    :: four = 4._iwp, two = 2._iwp;
  REAL(iwp)               :: r, s, t, t1, t2, t3, t4

  r = points(i,1)
  s = points(i,2)
  t = points(i,3)

  t1 = one - r - s -t
  t2 = r
  t3 = s
  t4 = t

  fun(1)  = t1*(two*t1 - one)
  fun(2)  = t2*(two*t2 - one)
  fun(3)  = t3*(two*t3 - one)
  fun(4)  = t4*(two*t4 - one)

  fun(5)  = four*t1*t2
  fun(6)  = four*t2*t3
  fun(7)  = four*t3*t1
  fun(8)  = four*t1*t4
  fun(9)  = four*t2*t4
  fun(10) = four*t3*t4
  RETURN
ENDSUBROUTINE SHAPE_FUN_TET10
  
!=====================================================
! Shape function first derivative iin Iso-P frame
!=====================================================
SUBROUTINE SHAPE_DER_TET10(der,points,i)
  INTEGER,   INTENT(IN)   :: i
  REAL(iwp), INTENT(IN)   :: points(:,:)
  REAL(iwp), INTENT(INOUT):: der(4,10)
  REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
  REAL(iwp), PARAMETER    :: four = 4._iwp, two = 2._iwp;
  REAL(iwp)               :: r, s, t, t1, t2, t3, t4

  r = points(i,1)
  s = points(i,2)
  t = points(i,3)

  t1 = one - r - s -t
  t2 = r
  t3 = s
  t4 = t

  der(1,1)  = four*t1 - one
  der(1,2)  = zero
  der(1,3)  = zero
  der(1,4)  = zero
  der(1,5)  = four*t2
  der(1,6)  = zero
  der(1,7)  = four*t3
  der(1,8)  = four*t4
  der(1,9)  = zero
  der(1,10) = zero

  der(2,1)  = zero
  der(2,2)  = four*t2 - one
  der(2,3)  = zero
  der(2,4)  = zero
  der(2,5)  = four*t1
  der(2,6)  = four*t3
  der(2,7)  = zero
  der(2,8)  = zero
  der(2,9)  = four*t4
  der(2,10) = zero

  der(3,1)  = zero
  der(3,2)  = zero
  der(3,3)  = four*t3 - one
  der(3,4)  = zero
  der(3,5)  = zero
  der(3,6)  = four*t2
  der(3,7)  = four*t1
  der(3,8)  = zero
  der(3,9)  = zero
  der(3,10) = four*t4

  der(4,1)  = zero
  der(4,2)  = zero
  der(4,3)  = zero
  der(4,4)  = four*t4 - one
  der(4,5)  = zero
  der(4,6)  = zero
  der(4,7)  = zero
  der(4,8)  = four*t1
  der(4,9)  = four*t2
  der(4,10) = four*t3
  RETURN
ENDSUBROUTINE SHAPE_DER_TET10

!=====================================================
! Jacobian-determinant and shapefunction derivative
! in global frame
!=====================================================
SUBROUTINE JACOBIAN_DERIV_TET10(deriv, det, coord, points, igauss)
  INTEGER                 :: i, j, m, n;
  INTEGER                 :: iterator1(6), iterator2(6);
  INTEGER,   INTENT(IN)   :: igauss
  REAL(iwp), INTENT(IN)   :: points(:,:), coord(10,3);
  REAL(iwp), INTENT(INOUT):: deriv(3,10), det
  REAL(iwp)               :: der(4,10)
  REAL(iwp), PARAMETER    :: zero = 0._iwp, one = 1._iwp;
  REAL(iwp), PARAMETER    :: four = 4._iwp, two = 2._iwp;
  REAL(iwp)               :: r, s, t, t_affine(4)
  REAL(iwp)               :: Jacobian(4,4), JacInv(4,4), Pmat(4,3)

  r = points(igauss,1)
  s = points(igauss,2)
  t = points(igauss,3)

  t_affine(1) = one - r - s -t
  t_affine(2) = r
  t_affine(3) = s
  t_affine(4) = t

  CALL SHAPE_DER_TET10(der,points,igauss)
  DO i = 1,4
    Jacobian(1,i) = one
    DO j = 1,3
      Jacobian(j+1,i) = DOT_PRODUCT(der(i,:), coord(:,j))
    ENDDO
  ENDDO
  det = DETERMINANT4(Jacobian)
  CALL INVERT2(Jacobian,JacInv,4)
  Pmat = TRANSPOSE(JacInv(2:4,1:4))

  iterator1 = (/1,2,3,1,2,3/)
  iterator2 = (/2,3,1,4,4,4/)
  deriv = zero;
  DO n = 1,4
    deriv(1,n) = (four*t_affine(n) - one)*Pmat(n,1)
    deriv(2,n) = (four*t_affine(n) - one)*Pmat(n,2)
    deriv(3,n) = (four*t_affine(n) - one)*Pmat(n,3)
  ENDDO

  DO m = 1,6 
    i = iterator1(m)
    j = iterator2(m)
    n = m + 4

    deriv(1,n) = four*(t_affine(i)*Pmat(j,1) + t_affine(j)*Pmat(i,1))
    deriv(2,n) = four*(t_affine(i)*Pmat(j,2) + t_affine(j)*Pmat(i,2))
    deriv(3,n) = four*(t_affine(i)*Pmat(j,3) + t_affine(j)*Pmat(i,3))
  ENDDO
  RETURN
ENDSUBROUTINE JACOBIAN_DERIV_TET10


!-------------------------------------
!  Determinant for 2by2 Matrix
!-------------------------------------
PURE FUNCTION DETERMINANT2(Amat) RESULT(det)
  IMPLICIT NONE
  REAL(iwp), INTENT(IN) :: Amat(2,2)
  REAL(iwp)             :: det

  det = 0._iwp
  det = Amat(1,1)*Amat(2,2) - Amat(1,2)*Amat(2,1)
ENDFUNCTION DETERMINANT2

!-------------------------------------
!  Determinant for 3by3 Matrix
!-------------------------------------
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
!  Determinant for 4by4 Matrix
!-------------------------------------
PURE FUNCTION DETERMINANT4(Amat) RESULT(det)
  IMPLICIT NONE
  REAL(iwp), INTENT(IN) :: Amat(4,4)
  REAL(iwp)             :: det

  det = 0._iwp
  det = det + Amat(1,1)*DETERMINANT3(Amat( (/2,3,4/), (/2,3,4/) ))
  det = det - Amat(1,2)*DETERMINANT3(Amat( (/2,3,4/), (/1,3,4/) ))
  det = det + Amat(1,3)*DETERMINANT3(Amat( (/2,3,4/), (/1,2,4/) ))
  det = det - Amat(1,4)*DETERMINANT3(Amat( (/2,3,4/), (/1,2,3/) ))
ENDFUNCTION DETERMINANT4

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
  REAL(iwp)               :: det,j11,j12,j13,j21,j22,j23,j31,j32,j33;
  REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;

    IF(n==4)THEN
      det = DETERMINANT4(MatA)
      MatAinv(1,1) =  DETERMINANT3(MatA( (/2,3,4/), (/2,3,4/) ))
      MatAinv(1,2) = -DETERMINANT3(MatA( (/2,3,4/), (/1,3,4/) ))
      MatAinv(1,3) =  DETERMINANT3(MatA( (/2,3,4/), (/1,2,4/) ))
      MatAinv(1,4) = -DETERMINANT3(MatA( (/2,3,4/), (/1,2,3/) ))

      MatAinv(2,1) = -DETERMINANT3(MatA( (/1,3,4/), (/2,3,4/) ))
      MatAinv(2,2) =  DETERMINANT3(MatA( (/1,3,4/), (/1,3,4/) ))
      MatAinv(2,3) = -DETERMINANT3(MatA( (/1,3,4/), (/1,2,4/) ))
      MatAinv(2,4) =  DETERMINANT3(MatA( (/1,3,4/), (/1,2,3/) ))

      MatAinv(3,1) =  DETERMINANT3(MatA( (/1,2,4/), (/2,3,4/) ))
      MatAinv(3,2) = -DETERMINANT3(MatA( (/1,2,4/), (/1,3,4/) ))
      MatAinv(3,3) =  DETERMINANT3(MatA( (/1,2,4/), (/1,2,4/) ))
      MatAinv(3,4) = -DETERMINANT3(MatA( (/1,2,4/), (/1,2,3/) ))

      MatAinv(4,1) = -DETERMINANT3(MatA( (/1,2,3/), (/2,3,4/) ))
      MatAinv(4,2) =  DETERMINANT3(MatA( (/1,2,3/), (/1,3,4/) ))
      MatAinv(4,3) = -DETERMINANT3(MatA( (/1,2,3/), (/1,2,4/) ))
      MatAinv(4,4) =  DETERMINANT3(MatA( (/1,2,3/), (/1,2,3/) ))
      MatAinv      =  MatAinv/det

    ELSE IF(n==3)THEN
      det =  MatA(1,1)*(MatA(2,2)*MatA(3,3)-MatA(3,2)*MatA(2,3))
      det =  det-MatA(1,2)*(MatA(2,1)*MatA(3,3)-MatA(3,1)*MatA(2,3))
      det =  det+MatA(1,3)*(MatA(2,1)*MatA(3,2)-MatA(3,1)*MatA(2,2))
      j11 =  MatA(2,2)*MatA(3,3)-MatA(3,2)*MatA(2,3)
      j21 = -MatA(2,1)*MatA(3,3)+MatA(3,1)*MatA(2,3)
      j31 =  MatA(2,1)*MatA(3,2)-MatA(3,1)*MatA(2,2)
      j12 = -MatA(1,2)*MatA(3,3)+MatA(3,2)*MatA(1,3)
      j22 =  MatA(1,1)*MatA(3,3)-MatA(3,1)*MatA(1,3)
      j32 = -MatA(1,1)*MatA(3,2)+MatA(3,1)*MatA(1,2)
      j13 =  MatA(1,2)*MatA(2,3)-MatA(2,2)*MatA(1,3)
      j23 = -MatA(1,1)*MatA(2,3)+MatA(2,1)*MatA(1,3)
      j33 =  MatA(1,1)*MatA(2,2)-MatA(2,1)*MatA(1,2)
      MatAinv(1,1) = j11
      MatAinv(1,2) = j12
      MatAinv(1,3) = j13
      MatAinv(2,1) = j21
      MatAinv(2,2) = j22
      MatAinv(2,3) = j23
      MatAinv(3,1) = j31
      MatAinv(3,2) = j32
      MatAinv(3,3) = j33
      MatAinv      = MatAinv/det
    ELSE IF(n==2)THEN
      det          =  MatA(1,1)*MatA(2,2)-MatA(1,2)*MatA(2,1)
      j11          =  MatA(1,1)
      MatAinv(1,1) =  MatA(2,2)
      MatAinv(2,2) =  j11
      MatAinv(1,2) = -MatA(1,2)
      MatAinv(2,1) = -MatA(2,1)
      MatAinv      =  MatA/det
    ELSE
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
    ENDIF
  RETURN
END SUBROUTINE Invert2

!-------------------------------------
!  Cholesky-factorisation
!-------------------------------------
SUBROUTINE CHOL_FACTORISATION(lower, Amat, n)
  IMPLICIT NONE
  INTEGER                   :: i,j,k; 
  INTEGER,     INTENT(IN)   :: n; 
  REAL(iwp),   INTENT(IN)   :: Amat(n,n);
  REAL(iwp),   INTENT(INOUT):: Lower(n,n);
  REAL(iwp),   PARAMETER    :: one = 1.0_iwp, zero = 0.0_iwp;
  REAL(iwp)                 :: partial;

  Lower = zero;
  DO I = 1,N
    DO J = 1,I
      partial = zero;
      DO K = 1,J
        partial = partial + Lower(I,K)*Lower(J,K);
      ENDDO
      IF(I == J) Lower(I,J) = DSQRT(Amat(I,J) - partial);
      IF(I /= J) Lower(I,J) = (one/Lower(J,J))*(Amat(I,J) - partial)
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE CHOL_FACTORISATION



!-------------------------------------
!  Incomplete Cholesky-factorisation
!-------------------------------------
SUBROUTINE ICHOL_FACTORISATION(lower, Amat, n, tol)
  IMPLICIT NONE
  INTEGER                   :: i,j,k; 
  INTEGER,     INTENT(IN)   :: n; 
  REAL(iwp),   INTENT(IN)   :: Amat(n,n), tol;
  REAL(iwp),   INTENT(INOUT):: lower(n,n);
  REAL(iwp),   PARAMETER    :: one = 1.0_iwp, zero = 0.0_iwp
  REAL(iwp)                 :: partial;

  Lower = Amat
  DO  k=1,n
    Lower(k,k) = DSQRT(Lower(k,k));
    DO i=(k+1),n
      IF(Lower(i,k)>tol)THEN
        Lower(i,k) = Lower(i,k)/Lower(k,k);            
      ENDIF
    ENDDO
    DO j=(k+1),n
      DO i=j,n
        IF(Lower(i,j)>tol)THEN
          Lower(i,j) = Lower(i,j)-Lower(i,k)*Lower(j,k);  
        ENDIF
      ENDDO
    ENDDO
  ENDDO

  DO i=1,n
    DO j=(i+1),n
      Lower(i,j) = 0;
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE ICHOL_FACTORISATION


!-------------------------------------
!  LU-factorisation
!-------------------------------------
SUBROUTINE LU_FACTORISATION(lower,upper,a,n)
 USE precision; IMPLICIT NONE
 INTEGER,  INTENT(IN)  :: n; 
 REAL(iwp),INTENT(IN)  :: a(n,n); 
 REAL(iwp),INTENT(OUT) :: lower(n,n),upper(n,n)
 INTEGER               :: i,j,k,l; 
 REAL(iwp)             :: total,zero=.0_iwp;

 upper=zero;
 lower=zero; 
 upper(1,:)=a(1,:)
 DO i=1,n; lower(i,i)=1.0_iwp; end do
 DO k=1,n-1
   IF(ABS(upper(k,k))>1.e-10_iwp)THEN
     DO i=k+1,n
!---Lower Triangular Factors---
       DO j=1,i-1
         total=zero
         DO l=1,j-1
           total= total-lower(i,l)*upper(l,j)
         END DO
         lower(i,j)=(a(i,j)+total)/upper(j,j)
       END DO
!---Upper Triangular Factors---
       DO j=1,n
         total=zero
         DO l=1,i-1
           total=total-lower(i,l)*upper(l,j)
         END DO
         upper(i,j)=a(i,j)+total
       END DO
     END DO
   ELSE
     WRITE(11,*)"Zero pivot found in row", k; EXIT
   END IF
 END DO
RETURN
END SUBROUTINE LU_FACTORISATION

!-------------------------------------
!  Forms (S)-vector for BICGSTAB(L) linear solver
!-------------------------------------
SUBROUTINE FORM_S2(GG, ell, kappa, omega, GAMMA, S)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)  :: ell;
  REAL(iwp), INTENT(IN)  :: GG(ell+1,ell+1), kappa;
  REAL(iwp), INTENT(OUT) :: omega, GAMMA(ell+1), S(ell+1);
  REAL(iwp)              :: p(ell-1), q(ell-1);
  REAL(iwp)              :: HH(ell-1,ell-1), GAMMA0(ell+1), GAMMA1(ell+1);
  REAL(iwp)              :: NGAMMA0, NGAMMA1, cosine;
  REAL(iwp), PARAMETER   :: one = 1._iwp, zero = 0._iwp;

  CALL INVERT2(-GG(2:ell,2:ell), HH, ell-1);
  p = MATMUL(HH, GG(2:ell, 1));
  q = MATMUL(HH, GG(2:ell, ell+1));

  GAMMA0(1)     = one;
  GAMMA0(ell+1) = zero;
  GAMMA0(2:ell) = p;

  GAMMA1(1)     = zero;
  GAMMA1(ell+1) = one;
  GAMMA1(2:ell) = q;

  NGAMMA0 = DOT_PRODUCT(GAMMA0,MATMUL(GG,GAMMA0))
  NGAMMA1 = DOT_PRODUCT(GAMMA1,MATMUL(GG,GAMMA1))
  omega   = DOT_PRODUCT(GAMMA0,MATMUL(GG,GAMMA1))
  cosine  = DABS(omega)/DSQRT(DABS(NGAMMA0*NGAMMA1))
  omega   = omega/NGAMMA1
  RETURN
END SUBROUTINE FORM_S2
!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Parallel_supplementary_Maths











