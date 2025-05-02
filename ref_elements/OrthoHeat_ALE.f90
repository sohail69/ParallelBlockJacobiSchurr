MODULE Orthotropic_Heat_ALE
  USE new_library;
  USE Parallel_supplementary_Maths;
  CONTAINS
!-------------------------------
!  Integrate ALE Orthotropic heat equation
!-------------------------------
SUBROUTINE Heat_Orthotropic_ALE_Element(StorKA, StorKB, coord, utemp, diff, fibre &
                                      , theta, dtim , points, weights, ndim, nod  &
                                      , ntots, nip)
  IMPLICIT NONE
  INTEGER                  :: igauss, i, j, k, l
  INTEGER  , INTENT(IN)    :: ndim, nod, ntots, nip;
  REAL(iwp), INTENT(IN)    :: points(nip,ndim), weights(nip), theta, dtim;
  REAL(iwp), INTENT(IN)    :: coord(nod,ndim), utemp(ndim*nod);
  REAL(iwp), INTENT(IN)    :: fibre(ndim,ndim), diff(ndim);
  REAL(iwp), INTENT(INOUT) :: StorKA(ntots,ntots), StorKB(ntots,ntots);
  REAL(iwp), PARAMETER     :: one = 1._iwp, zero = 0._iwp;
  REAL(iwp)                :: der(ndim,nod), deriv(ndim,nod);
  REAL(iwp)                :: fun(nod), funny(1,nod);
  REAL(iwp)                :: Jac(ndim,ndim), Jacinv(ndim,ndim), det;
  REAL(iwp)                :: pm(ntots,ntots), kc(ntots,ntots), auxm(nod,ndim);
  REAL(iwp)                :: Kay(ndim,ndim), kay_t(ndim,ndim), kay2(ndim,ndim);
  REAL(iwp)                :: Fdef(ndim,ndim), Cdef(ndim,ndim), CdefInv(ndim,ndim)
  LOGICAL                  :: Hexa27test, Tetra10Test, standardTest;

  !-----
  !Calculate material properties
  !-----
  der = zero;  deriv  = zero;
  fun = zero;  funny  = zero;
  Jac = zero;  Jacinv = zero;  det    = zero;
  Kay = zero;  kay_t  = zero;  StorKA = zero;
  pm  = zero;  kc     = zero;  StorKB = zero;


  DIMENSIONS:DO j = 1,ndim
    kay_t = zero;
    CALL VEC_OUTERPRODUCT(kay_t, fibre(:,j), fibre(:,j), ndim, ndim);
    kay = kay + diff(j)*kay_t
  ENDDO DIMENSIONS


  !-----
  !Update the nodal coordinates element
  !-----
  l = 1;
  DO i = 1,nod
    DO j = 1,ndim
      auxm(i,j) = utemp(l);
      l = l + 1;
    ENDDO
  ENDDO

  Hexa27test   = ((nod==27).AND.(ndim==3))
  Tetra10Test  = ((nod==10).AND.(ndim==3))
  standardTest = (.NOT.(Hexa27test)).AND.(.NOT.(Tetra10Test))

  !-----
  !Integrate elements in current ref-frame
  !-----
  kay2 = zero;
  GAUSS_PTS2:DO igauss = 1,nip
    !If using Hexa27 elements
    IF(Hexa27test)   CALL SHAPE_FUN_HEX27(fun,points,igauss)
    IF(Hexa27test)   CALL SHAPE_DER_HEX27(der,points,igauss)

    !If using Tetra10 elements (special treatment of deriv and det)
    IF(Tetra10Test) CALL SHAPE_FUN_TET10(fun,points,igauss);
    IF(Tetra10Test) CALL JACOBIAN_DERIV_TET10(deriv, det, coord, points, igauss)

    !Other parafem standard elements
    IF(standardTest) CALL SHAPE_FUN(fun,points,igauss);
    IF(standardTest) CALL SHAPE_DER(der,points,igauss);

    IF(standardTest.OR.Hexa27test)THEN
      Jac = MATMUL(der,coord)
      det = determinant(Jac);
      CALL INVERT2(Jac, Jacinv, ndim);
      deriv = MATMUL(Jacinv,der);
    ENDIF

    funny(1,:) = fun(:);
    Fdef = MATMUL(deriv,auxm)
    DO i = 1,ndim; Fdef(i,i) = Fdef(i,i) + one; ENDDO;
    Cdef = MATMUL(TRANSPOSE(Fdef),Fdef)
    CALL INVERT2(Cdef, CdefInv, ndim)
    kay2 = MATMUL(kay,CdefInv)
    kc=kc+MATMUL(MATMUL(TRANSPOSE(deriv),kay2),deriv)*det*weights(igauss);
    pm=pm+MATMUL(TRANSPOSE(funny),funny)*det*weights(igauss);
  ENDDO GAUSS_PTS2
  StorKA = pm+kc*theta*dtim;
  StorKB = pm-kc*(one - theta)*dtim;
  RETURN
ENDSUBROUTINE Heat_Orthotropic_ALE_Element

!-------------------------------
!-------------------------------
!-------------------------------
END MODULE Orthotropic_Heat_ALE















