PROGRAM MAIN
  USE SolidPFEM_UAS;
  USE TensorElement
  IMPLICIT NONE
  INTEGER :: I, J;
  INTEGER,   PARAMETER :: nprop=2, ndim=1, nip1D=2, nod=5, pOrder=4;
  INTEGER,   PARAMETER :: ntots=ndim*nod;
  INTEGER,   PARAMETER :: nst=(((ndim+1)*ndim)/2)
  INTEGER,   PARAMETER :: nip = nip1D**ndim;
  REAL(iwp), PARAMETER :: zero = 0._iwp, one = 1._iwp;
  REAL(iwp)            :: gama(nod), fibre(ndim,ndim), utemp(ntots), MATPROP(nprop);
  REAL(iwp)            :: Ni1D(nod,nip1D), dNi1D(nod,nip1D), points1D(nip1D),  weights1D(nip1D);  
  REAL(iwp)            :: Ni(nod,nip), dNi(ndim,nod,nip);
  REAL(iwp)            :: coord(nod,ndim), weights(nip);
  REAL(iwp)            :: Km(ntots,ntots), Rm(ntots);

  CALL GaussLengendre1D(weights1D,points1D,nip1D)
  CALL CalculateNDWeights(weights,weights1D,nip,nip,nDIM)
  CALL TENSOR_ELEMENT_NDPoly(Ni, dNi, points1D, pOrder, nip1D, nod, nip, nDIM)

   MATPROP(1) = 1.00_iwp;
   MATPROP(2) = 0.49999_iwp;
   gama  = 0.01;
   fibre = 1.0_iwp
   utemp = 0.0001_iwp;

   coord(1,1) = 0._iwp
   coord(2,1) = 1._iwp

  CALL SOLID_PFEM_UAS_ELM(Km, Rm, utemp, gama, fibre, MATPROP   &
                        , coord, Ni, dNi, weights, nprop, ntots &
                        , ndim, nst, nip, nod)



  OPEN(UNIT=(12),FILE = "testResult/HierTest.txt")
  DO I = 1,ntots
    WRITE(12,*) Km(I,:), " | ", Rm(I)
  ENDDO
  CLOSE(12)



END PROGRAM MAIN
