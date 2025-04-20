MODULE TensorElement
INTEGER, PARAMETER :: iwp = SELECTED_REAL_KIND(15,300)
CONTAINS

!/***************************************\
!|***************************************|
! Calculations relating to construction
! of the Gauss-Legendre integration
! rules in 1-D using the Golub-Welsch
! algorithm
! This include the sample points 
! and weights
! Inputs:
!   nip = The integration order
! Outputs:
!   weights = The Gauss-Legendre
!             integration weights
!   points  = The Gauss-Legendre
!             sample points
!|***************************************|
!\***************************************/
SUBROUTINE GaussLengendre1D(weights,points,nip)
  IMPLICIT NONE
  INTEGER                 :: I, J, K;
  INTEGER,   INTENT(IN)   :: nip;
  REAL(iwp), INTENT(INOUT):: weights(nip),points(nip);
  REAL(iwp), PARAMETER    :: one=1._iwp, two=2._iwp, four=4._iwp;
  REAL(iwp), PARAMETER    :: zero=0._iwp, eps=1.0E-15_iwp;
  REAL(iwp)               :: JM_e(nip-1), OffDiag, I_d;
  REAL(iwp)               :: eigVecs(nip,nip);
  REAL(iwp), ALLOCATABLE  :: work(:)
  INTEGER(8)              :: info, lwork, nip2

  nip2 = nip;
  lwork = MAX(1,2*nip2-2)
  ALLOCATE(work(lwork))

  !
  ! Construct the Jacobi
  ! matrix
  weights=0._iwp;
  points=0._iwp;
  DO I=1,nip-1
    I_d = REAL(I,iwp);
    OffDiag = I_d/DSQRT(four*I_d*I_d - one);
    JM_e(I) = OffDiag;
  ENDDO

  !
  ! Solve the Eigen
  ! value problem
  CALL DSTEV('V', nip2, points, JM_e, eigVecs, nip2, work, info)
  IF (info /= 0) THEN
    PRINT*, "Error: DSTEV did not converge, INFO =", info
    STOP
  ENDIF

  !
  ! Get the integration
  ! points and weights
  DO I = 1, nip
    weights(I) = two * ( eigVecs(1,I)*eigVecs(1,I))
    IF(weights(I)       < eps) weights(I) = zero;
    IF(DABS(points(I))  < eps) points(I)  = zero;
  ENDDO

  DEALLOCATE(work)
  RETURN
END SUBROUTINE

!/***************************************\
!|***************************************|
! Calculations relating to
! recursive formulas of legendre
! polynomials, their derivatives
! shapefunctions and their 
! derivatives
!|***************************************|
!\***************************************/

!/***************************************\
! Calculates the i'th legendre
! polynomial weight at a
! integration point
! Inputs:
!   I   = the term order
!   r   = position of sampled
!         integration point
!   Ln1 = The legendre poly-function
!         weight of order I-1
!         sampled at r
!   Ln2 = The legendre poly-function
!         weight of order I-2
!         sampled at r
! Outputs:
!   Ln0 = The legendre poly-function
!         weight of order I
!         sampled at r
!\***************************************/
FUNCTION LEGENDRE_PolyFunc(I, r, Ln1, Ln2) RESULT(Ln0)
  INTEGER,   INTENT(IN):: I
  REAL(iwp), INTENT(IN):: r, Ln1, Ln2
  REAL(iwp)            :: Ln0, N
  REAL(iwp), PARAMETER :: one=1._iwp, two=2._iwp;
  N = REAL(I,iwp)
  
  IF(I==1)              Ln0 = one
  IF(I==2)              Ln0 = r
  IF((I/=1).AND.(I/=2)) Ln0 = (one/N)*( (two*N - one)*r*Ln1 - (N- one)*Ln2);
ENDFUNCTION LEGENDRE_PolyFunc

!/***************************************\
! Calculates the i'th legendre
! polynomial weight at a
! integration point
! Inputs:
!   I    = the term order
!   r    = position of sampled
!         integration point
!   Ln1  = The legendre poly-function
!          weight of order I-1
!          sampled at r
!   dLn1 = The legendre poly-function
!          derivative weight of order I-1
!          sampled at r
!   dLn2 = The legendre poly-function
!          derivative weight of order I-2
!          sampled at r
! Outputs:
!   dLn0 = The legendre poly-function
!          derivative weight of order I
!          sampled at r
!\***************************************/
FUNCTION LEGENDRE_PolyDerFunc(I, r, Ln1, dLn1, dLn2) RESULT(dLn0)
  INTEGER,   INTENT(IN):: I
  REAL(iwp), INTENT(IN):: r, Ln1, dLn1, dLn2
  REAL(iwp)            :: dLn0, Ni, N
  REAL(iwp),  PARAMETER:: one=1._iwp, two=2._iwp;
  N = REAL(I,iwp)
  
  IF(I==1)              dLn0 = 0._iwp
  IF(I==2)              dLn0 = one
  IF((I/=1).AND.(I/=2)) dLn0 = (one/N)*( (two*N - one)*(Ln1 + r*dLn1) - (N- one)*dLn2 );
ENDFUNCTION LEGENDRE_PolyDerFunc

!/***************************************\
! Calculates the i'th legendre
! shape function weight at a
! integration point
! Inputs:
!   I    = the term order
!   r    = position of sampled
!         integration point
!   Ln0  = The legendre poly-function
!          weight of order I
!          sampled at r
!   Ln2 = The legendre poly-function
!         weight of order I-2
!         sampled at r
! Outputs:
!   Ni = The legendre polynomial 
!        shape function sampled at
!        point r
!\***************************************/
FUNCTION LEGENDRE_ShapeFunc(I, r, Ln0, Ln2) RESULT(Ni)
  INTEGER,   INTENT(IN):: I
  REAL(iwp), INTENT(IN):: r, Ln0, Ln2
  REAL(iwp)            :: Ni, N
  REAL(iwp),  PARAMETER:: half=0.5_iwp, one=1._iwp, two=2._iwp, four=4._iwp;
  N = REAL(I,iwp)

  IF(I==1)              Ni = half*(one + r);
  IF(I==2)              Ni = half*(one - r);
  IF((I/=1).AND.(I/=2)) Ni = (Ln0 - Ln2)/DSQRT(four*N - two);
ENDFUNCTION LEGENDRE_ShapeFunc

!/***************************************\
! Calculates the derivative of the
! i'th legendre shape function 
! weight at a integration point
! Inputs:
!   I    = the term order
!   dLn0 = The legendre poly-function
!          derivative weight of order I
!          sampled at r
!   dLn2 = The legendre poly-function
!          derivative weight of order I-2
!          sampled at r
! Outputs:
!   dNi = The legendre polynomial 
!         shape function derivative
!         sampled at point r
!\***************************************/
FUNCTION LEGENDRE_ShapeDerFunc(I, dLn0, dLn2) RESULT(dNi)
  INTEGER,   INTENT(IN):: I
  REAL(iwp), INTENT(IN):: dLn0, dLn2
  REAL(iwp)            :: dNi, N
  REAL(iwp),  PARAMETER:: half=0.5_iwp, one=1._iwp, two=2._iwp, four=4._iwp;
  N = REAL(I,iwp)
  
  IF(I==1)              dNi =  half
  IF(I==2)              dNi = -half;
  IF((I/=1).AND.(I/=2)) dNi = (dLn0 - dLn2)/DSQRT(four*N - two);
ENDFUNCTION LEGENDRE_ShapeDerFunc


!/***************************************\
! Using the recursion formula, form the 
! lengedre polynomial weight series for
! a polynomial of order-pOrder and use
! that to form the shape functions and
! derivatives for a 1-D hierachical
! Legendre-line-element
! Inputs:
!   points = Integration sample points
!   pOrder = The polynomial order
!   nip    = The number of integration
!            sample points
! Outputs:
!   Legendre1D_Ni  = The 1D legendre
!                    polynomial shape
!                    function weights
!   Legendre1D_dNi = The 1D legendre
!                    polynomial shape
!                    function derivative
!                    weights
!\***************************************/
SUBROUTINE TENSOR_ELEMENT_1DPoly(Legendre1D_Ni, Legendre1D_dNi, points, pOrder, nip)
  INTEGER  :: Igauss, I;
  INTEGER  :: pOrder, nip;
  REAL(iwp):: Ln0, Ln1, Ln2, dLn0, dLn1, dLn2, r;
  REAL(iwp):: Legendre_P(pOrder+1,nip), Legendre_dP(pOrder+1,nip)
  REAL(iwp):: Legendre1D_Ni(pOrder+1,nip), Legendre1D_dNi(pOrder+1,nip)
  REAL(iwp):: points(nip)
  REAL(iwp),PARAMETER :: zero=0._iwp

  DO Igauss=1,nip
    r = points(Igauss)
    DO I = 1,pOrder+1
	  Ln1=zero; Ln2=zero; dLn1=zero; dLn2=zero;
      IF(I>2)THEN
        Ln1  = Legendre_P(I-1,Igauss)
        Ln2  = Legendre_P(I-2,Igauss)
        dLn1 = Legendre_dP(I-1,Igauss)
        dLn2 = Legendre_dP(I-2,Igauss)
      ENDIF
      Legendre_P(I,Igauss)  = LEGENDRE_PolyFunc(I, r, Ln1, Ln2);
	  Legendre_dP(I,Igauss) = LEGENDRE_PolyDerFunc(I, r, Ln1, dLn1, dLn2);
    ENDDO

    DO I = 1,pOrder+1
	  Ln0=zero; Ln2=zero; dLn0=zero; dLn2=zero;
      IF(I>2)THEN
        Ln0  = Legendre_P(I,Igauss)
        Ln2  = Legendre_P(I-2,Igauss)
        dLn0 = Legendre_dP(I,Igauss)
        dLn2 = Legendre_dP(I-2,Igauss)
      ENDIF
      Legendre1D_Ni(I,Igauss)  = LEGENDRE_ShapeFunc(I, r, Ln0, Ln2);
	  Legendre1D_dNi(I,Igauss) = LEGENDRE_ShapeDerFunc(I, dLn0, dLn2);
    ENDDO
  ENDDO
  RETURN 
ENDSUBROUTINE TENSOR_ELEMENT_1DPoly


!/***************************************\
! Using the recursion formula, form the 
! lengedre polynomial weight series for
! a polynomial of order-pOrder and use
! that to form the shape functions and
! derivatives for a N-D hierachical
! Legendre-element (1-Line, 2-Quad,
! 3-Hex..)
! Inputs:
!   points   = Integration sample points
!              for 1D element
!   pOrder1D = The polynomial order
!              for 1D element
!   nip1D    = The number of integration
!              sample points  in 1D-element
!   nip      = The number of integration
!              sample points in ND-element
!   nod      = The number of DOFs/nodes on
!              the ND-element
!   nDIM     = The dimension of the element
! Outputs:
!   Legendre_Ni  = The ND legendre
!                  polynomial shape
!                  function weights
!   Legendre_dNi = The ND legendre
!                  polynomial shape
!                  function derivative
!                  weights
!\***************************************/
SUBROUTINE TENSOR_ELEMENT_NDPoly(Legendre_Ni, Legendre_dNi, points, pOrder1D, nip1D, nod, nip, nDIM)
  INTEGER                :: Lgauss, Jgauss, L, P, Q, K;
  INTEGER                :: Igauss(ndim), I(ndim);
  INTEGER,  INTENT(IN)   :: nDIM, pOrder1D, nip1D, nod, nip;
  REAL(iwp),INTENT(IN)   :: points(nip);
  REAL(iwp),INTENT(INOUT):: Legendre_Ni(nod,nip), Legendre_dNi(nDIM,nod,nip)
  REAL(iwp)              :: Legendre_1DNi(pOrder1D+1,nip1D), Legendre_1DdNi(pOrder1D+1,nip1D)
  REAL(iwp),  PARAMETER  :: zero=0._iwp, one=1._iwp, two=2._iwp;

  CALL TENSOR_ELEMENT_1DPoly(Legendre_1DNi, Legendre_1DdNi, points, pOrder1D, nip1D)
  
  DO Lgauss = 1,nip
    CALL InverseIterator(Lgauss,Igauss,nip1D,nDIM);
    DO L = 1,nod
      CALL InverseIterator(L,I,pOrder1D+1,nDIM);
      Legendre_Ni(L,Lgauss) = one;
      DO P = 1,nDIM
        Jgauss = Igauss(P)
        K = I(P)
        Legendre_Ni(L,Lgauss) = Legendre_Ni(L,Lgauss)*Legendre_1DNi(K,Jgauss);
        Legendre_dNi(P,K,Lgauss) = one;
        DO Q = 1,nDIM
          IF(Q==P) Legendre_dNi(P,K,Lgauss) = Legendre_dNi(P,K,Lgauss)*Legendre_1DdNi(K,Jgauss);
          IF(Q/=P) Legendre_dNi(P,K,Lgauss) = Legendre_dNi(P,K,Lgauss)*Legendre_1DNi(K,Jgauss);
        ENDDO
      ENDDO 
    ENDDO
  ENDDO
  RETURN 
ENDSUBROUTINE TENSOR_ELEMENT_NDPoly


!/***************************************\
! Calculates the gaussian integration rules
! weights for an N-dimensional Tensor
! element
! Inputs:\
!   M           = Flattened ND iterator
!   nDIM        = Dimension of the problem
!   nSize1D     = Size of each iterator
!                 in reference 1D
! Outputs:
!   Iters(nDIM) = Unpacked ND Iterators
!                 in each dimension
!\***************************************/
SUBROUTINE CalculateNDWeights(weightsND,weights1D,nip,nip1D,nDIM)
  INTEGER, INTENT(IN)    :: nip, nip1D, nDIM;
  REAL(iwp),INTENT(IN)   :: weights1D(nip1D);
  REAL(iwp),INTENT(INOUT):: weightsND(nip)
  INTEGER                :: I, J, K, Iters(nDIM);

  DO I=1,nip
    CALL InverseIterator(I,Iters,nip1D,nDIM);
    weightsND(I) = 1._iwp;
    DO J=1,nDIM
      weightsND(I) = weightsND(I)*weights1D(Iters(J));
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE CalculateNDWeights


!/***************************************\
! Inverse iterator, takes a 1D array input
! and converts it to n-D input used to help
! construct arbitrary dimensional tensor 
! shape functions from a single dimensional 
! line shape functions
! Inputs:
!   M           = Flattened ND iterator
!   nDIM        = Dimension of the problem
!   nSize1D     = Size of each iterator
!                 in reference 1D
! Outputs:
!   Iters(nDIM) = Unpacked ND Iterators
!                 in each dimension
!\***************************************/
SUBROUTINE InverseIterator(M,Iters,nSize1D,nDIM)
  INTEGER, INTENT(IN)   :: M, nSize1D, nDIM
  INTEGER, INTENT(INOUT):: Iters(nDIM);
  INTEGER               :: p, q;
  INTEGER               :: I, J, K;
  Iters = 0
  DO P = 1,nDIM
    I = P-1;
    K = M-1;
    DO Q = 1,P-1
      K = (K - (Iters(Q)) );
      IF(P/=1) K = K/nSize1D
    ENDDO
    Iters(P) = MOD(K,nSize1D);
  ENDDO
  Iters = Iters + 1;
  RETURN
ENDSUBROUTINE InverseIterator


!/***************************************\
!\***************************************/
ENDMODULE TensorElement
