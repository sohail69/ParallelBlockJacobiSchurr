!------------------------------------------------------------------------------
! This is a file that interfaces the Fortran90/95 functions
! in the PARAFEM, src and Ref_element libraries to be made available
! to InterfaceC.h file and C/C++ functions
! This acts as a single interface for the sake of brevity though it
! could be easily split up into various components e.g IO, mesh, linear
! solvers, elements etc.
!
! Author: Sohail Rathore
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!
!                     Parallel Reading/Input functions
!
!------------------------------------------------------------------------------
!-------------------------------
! Initialise parafem and read in Job
!-------------------------------
SUBROUTINE READ_JOBInitialise(job_name, nlen, element, partitioner, npri, numpes   &
                            , mesh, npess, nod, nip, nbnd, nFace, nodFace, nipFace &
                            , iels_start, nel_pp, nels, nn, np_types)
  USE precision;      USE global_variables;  USE MP_INTERFACE;
  USE maths;          USE gather_scatter;    USE new_library;
  USE Parallel_IO;
  IMPLICIT NONE
  CHARACTER(LEN=50), INTENT(IN)   :: job_name
  CHARACTER(LEN=15), INTENT(INOUT):: element;
  INTEGER,           INTENT(IN)   :: nlen;
  INTEGER,           INTENT(INOUT):: partitioner, npri;
  INTEGER,           INTENT(INOUT):: numpes, mesh, npess;
  INTEGER,           INTENT(INOUT):: nod, nip, nbnd, nFace, nodFace, nipFace;
  INTEGER,           INTENT(INOUT):: iels_start, nel_pp, nels, nn, np_types;

  CALL find_pe_procs(numpes,npess);
  CALL READ_JOB_DATA(job_name, nlen, numpes, element, mesh, partitioner   &
                   , np_types, nels, nn, nod, nip, nbnd, nodFace, nipFace &
                   , nFace, npri);
  CALL calc_nels_pp(job_name,nels,npess,numpes,partitioner,nels_pp)
  iels_start = iel_start;
  ntot       = nod;
  nel_pp     = nels_pp;
  RETURN
END SUBROUTINE READ_JOBInitialise


!-------------------------------
! Read in the Mesh
!-------------------------------
SUBROUTINE READ_MESH(JOB_NAME, element, coord_pp, gnum_pp, etype_pp         &
                   , nod, ndim, nn, nel_pp, neqs_pp, iels_start, ieqs_start &
                   , meshgen, numpes, npess, npes_pp)
  USE precision;      USE global_variables;  USE MP_INTERFACE;
  USE maths;          USE gather_scatter;    USE new_library;
  USE Input;
  IMPLICIT NONE
  CHARACTER(LEN=50), INTENT(IN):: JOB_NAME;
  CHARACTER(LEN=15), INTENT(IN):: element;
  INTEGER                      :: iel, i;
  INTEGER,   INTENT(IN)        :: nel_pp, iels_start;
  INTEGER,   INTENT(IN)        :: nn, npess, numpes;
  INTEGER,   INTENT(IN)        :: meshgen, nod, ndim;
  INTEGER,   INTENT(INOUT)     :: neqs_pp, npes_pp, ieqs_start;
  INTEGER,   INTENT(INOUT)     :: etype_pp(nel_pp), gnum_pp(nod,nel_pp)
  REAL(iwp), INTENT(INOUT)     :: coord_pp(nod,ndim,nel_pp);

  iel_start = iels_start;
  nels_pp   = nel_pp;
  npes      = npess
  numpe     = numpes;
  ntot      = nod;
  CALL READ_ELEMENTS(JOB_NAME,iel_start,nn,npes,numpe,etype_pp,gnum_pp)
  IF((meshgen==2).AND.(nod/=27)) CALL abaqus2sg(element,gnum_pp)
  CALL READ_G_COORD_PP(JOB_NAME,gnum_pp,nn,npes,numpe,coord_pp)
  CALL CALC_NEQ_PP2(nn, neqs_pp, ieqs_start)
!  CALL CALC_NPES_PP(npes,npes_pp);
  npes_pp = npes;
  neq_pp = neqs_pp;
  ieq_start = ieqs_start;
  CALL MAKE_GGL3(npes,npes,gnum_pp);
  RETURN
END SUBROUTINE READ_MESH


!-------------------------------
! Read Boundary
!-------------------------------
SUBROUTINE READ_BOUNDARYNODES(job_name, boundary_N, nbnd, nodof, numpes)
  USE new_library;
  USE Parallel_IO;
  IMPLICIT NONE
  CHARACTER(LEN=50), INTENT(IN) :: job_name;
  INTEGER,INTENT(IN)            :: numpes, nodof, nbnd;
  INTEGER,INTENT(INOUT)         :: boundary_N(nbnd,nodof+1)

  IF(nbnd /= 0)THEN
    CALL READ_BOUNDARY(job_name,numpes,boundary_N,nodof,nbnd);
    CALL REARRANGE(boundary_N);
  ENDIF
  RETURN
END SUBROUTINE READ_BOUNDARYNODES

!-------------------------------
! Read material properties
!-------------------------------
SUBROUTINE READ_MATERIALPROPS(job_name, numpes, matprops, nmat, np_types)
  USE precision;
  USE new_library;
  USE Parallel_IO;
  IMPLICIT NONE
  CHARACTER(LEN=50),INTENT(IN)   :: job_name;
  INTEGER,          INTENT(IN)   :: numpes, nmat, np_types;
  REAL(iwp),        INTENT(INOUT):: matprops(nmat,np_types);

  CALL READ_MATERIAL_DATA(job_name, numpes, matprops, nmat, np_types)
  RETURN
END SUBROUTINE READ_MATERIALPROPS

!-------------------------------
! Read partitioned Real data
!-------------------------------
SUBROUTINE READ_DATAR(JOB_NAME, numpes, npess, data_pp,  n_pp, dof)
  USE precision;
  USE new_library;
  USE Parallel_IO;
  IMPLICIT NONE
  CHARACTER(LEN=50),INTENT(IN)   :: JOB_NAME;
  INTEGER,          INTENT(IN)   :: numpes, npess, n_pp, dof;
  REAL(iwp),        INTENT(INOUT):: data_pp(dof,n_pp);

  CALL READ_DATA_REAL(JOB_NAME,numpes,npess,data_pp,n_pp,dof)
  RETURN
END SUBROUTINE READ_DATAR

!-------------------------------
! Read partitioned Integer data
!-------------------------------
SUBROUTINE READ_DATAI(JOB_NAME, numpes, npess, data_pp,  n_pp, dof)
  USE precision;
  USE new_library;
  USE Parallel_IO;
  IMPLICIT NONE
  CHARACTER(LEN=50),INTENT(IN)   :: JOB_NAME;
  INTEGER,          INTENT(IN)   :: numpes, npess, n_pp, dof;
  INTEGER,          INTENT(INOUT):: data_pp(dof,n_pp);

  CALL READ_DATA_INTEGER(JOB_NAME,numpes,npess,data_pp,n_pp,dof)
  RETURN
END SUBROUTINE READ_DATAI

!------------------------------------------------------------------------------
!
!                     Parallel calculation functions
!
!------------------------------------------------------------------------------
!-------------------------------
! Set nodal freedoms from boundary nodes
!-------------------------------
SUBROUTINE BoundarySetNF(gg_pp, gnum_pp, boundary_N, ndim, nodof &
                       , ntots, nod, nbnd, nel_pp, npess, numpes)
  USE precision;
  USE new_library;
  USE MP_INTERFACE;
  IMPLICIT NONE
  INTEGER                    :: IEL, I, M, N;
  INTEGER,   INTENT(IN)      :: ndim,nodof,nod,ntots,nbnd,nel_pp,npess,numpes;
  INTEGER,   INTENT(IN)      :: gnum_pp(nod,nel_pp), boundary_N(nbnd,nodof+1);
  INTEGER,   INTENT(INOUT)   :: gg_pp(nod*nodof,nel_pp);

  DO IEL = 1,nel_pp
    DO I = 1,nodof
      M=(I-1)*nod+1
      N=I*nod;
      CALL FIND_G3(gnum_pp(:,IEL),gg_pp(M:N,IEL),boundary_N(:,(/1,I+1/)))
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE BoundarySetNF

!-------------------------------
! Find faces on the Boundary based on given boundary nodes
! and makes face steering array
!-------------------------------
SUBROUTINE BoundaryFaceDetection(gg_Face, nloadedFace, gnum_pp, boundary_N    &
                               , boundaryID, element, ndim, nodof, nod, nFace &
                               , nodFace, nbnd, nel_pp, npess, numpes)
  USE precision;
  USE new_library;
  USE Parallel_BoundaryConditions;
  USE MP_INTERFACE;
  IMPLICIT NONE
  CHARACTER(len=15), INTENT(IN):: element
  INTEGER,   INTENT(IN)        :: ndim, nodof, nod, nFace, nodFace, nbnd, nel_pp;
  INTEGER,   INTENT(IN)        :: gnum_pp(nod,nel_pp), npess, numpes;
  INTEGER,   INTENT(IN)        :: boundary_N(nbnd,nodof+1), boundaryID;
  INTEGER,   INTENT(INOUT)     :: gg_Face(nFace+1,nel_pp), nloadedFace;

  gg_Face = 0;
  CALL BoundaryFaceElement_detection(gg_Face, nloadedFace, gnum_pp, boundary_N(:,1) &
                                   , boundaryID, element, ndim, nod, nFace, nodFace &
                                   , nbnd, nel_pp, npess, numpes)
  RETURN
END SUBROUTINE BoundaryFaceDetection

!-------------------------------
! Preconditioner Graph colouring
!-------------------------------
SUBROUTINE PRECONDITIONER_COLOURING(gg_colour,ncolour,gg_pp,nod,nn_pp,nel_pp &
                                   ,nn,nels,npess,numpes)
  USE precision;      USE global_variables;  USE MP_INTERFACE;
  USE maths;          USE gather_scatter;    USE new_library;
  USE Parallel_supplementary_Maths;
  USE GRAPH_COLOURING;
  IMPLICIT NONE
  INTEGER                :: ncol0
  INTEGER,  INTENT(IN)   :: nod, nn_pp, nel_pp, nn, nels, npess, numpes;
  INTEGER,  INTENT(IN)   :: gg_pp(nod,nel_pp);
  INTEGER,  INTENT(INOUT):: gg_colour(nel_pp), ncolour;

  ncolour = 0;
  CALL ELEMENT_COLOURING(gg_colour,gg_pp,nod,nn_pp,nel_pp,nn,nels,npess,numpes)
  ncol0 = MAXVAL(gg_colour)
  CALL MPI_ALLREDUCE(ncol0,ncolour,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier);
  IF(numpes==1) WRITE(*,*) "ELEMENT CHROMATIC NUMBER: ", ncolour
  RETURN
END SUBROUTINE PRECONDITIONER_COLOURING

!-------------------------------
! Linear Solvers
!-------------------------------
SUBROUTINE LinearSolve(A_mat, M_mat, x_vec, b_vec, NodalMask  &
                     , gg_colour, ncolours, ntots, nod, nodof, nel_pp, neqs_pp &
                     , nn_pp, ltol, limit, iters, ell, error, solver, precon)
  USE precision;      USE global_variables;  USE MP_INTERFACE;
  USE maths;          USE gather_scatter;    USE new_library;
  USE Parallel_supplementary_Maths;
  USE Parallel_FEA_LinearSolvers;
  USE PRECONDITIONERS;
  IMPLICIT NONE
  INTEGER,   INTENT(INOUT):: iters;
  INTEGER,   INTENT(IN)   :: nodof, nod, nn_pp, ncolours, solver, precon;
  INTEGER,   INTENT(IN)   :: NodalMask(nodof,nod), gg_colour(nel_pp);
  INTEGER,   INTENT(IN)   :: ntots, nel_pp, neqs_pp, limit, ell;
  REAL(iwp), INTENT(IN)   :: M_mat(ntots,ntots,nel_pp); !preconditioner Matrix
  REAL(iwp), INTENT(IN)   :: A_mat(ntots,ntots,nel_pp), b_vec(neqs_pp), ltol;
  REAL(iwp), INTENT(INOUT):: x_vec(neqs_pp), error;

  SELECT CASE(solver)
    CASE(1) !Conjugate gradient method
      CALL SolveLinearSystem_CG(A_mat, x_vec, b_vec &
                              , NodalMask, gg_colour, ntots, nod, ncolours  &
                              , nodof, nel_pp, neqs_pp, nn_pp, ltol, limit  &
                              , iters, error, precon)

    CASE(2) !Stabilised Hybrid Bi-Conjugate polynomial gradient(l) method
      CALL SolveLinearSystem_BICGSTABL(A_mat, x_vec, b_vec, NodalMask, ntots, nod &
                                     , nodof, nel_pp, neqs_pp, nn_pp, ltol, limit &
                                     , iters, ell, error)

    CASE(3) !Stabilised Bi-Conjugate method
      CALL SolveLinearSystem_BICGSTAB(A_mat, x_vec, b_vec  &
                                   , NodalMask, gg_colour, ncolours, ntots, nod  &
                                   , nodof, nel_pp, neqs_pp, nn_pp, ltol, limit  &
                                   , iters, ell, error, precon)

    CASE(4) !Restarted Generalised minimum residual method Variant 1 Full GMRES
      CALL SolveLinearSystem_GMRESR(A_mat, M_mat, x_vec, b_vec  &
                                   , NodalMask, gg_colour, ncolours, ntots, nod  &
                                   , nodof, nel_pp, neqs_pp, nn_pp, ltol, limit  &
                                   , iters, ell, error, precon, 0)

    CASE(5) !Restarted Generalised minimum residual method Variant 2 DQGMRES
      CALL SolveLinearSystem_GMRESR(A_mat, M_mat, x_vec, b_vec  &
                                   , NodalMask, gg_colour, ncolours, ntots, nod  &
                                   , nodof, nel_pp, neqs_pp, nn_pp, ltol, limit  &
                                   , iters, ell, error, precon, 1)

    CASE DEFAULT !Default to Conjugate gradient method
      CALL SolveLinearSystem_CG(A_mat, x_vec, b_vec &
                              , NodalMask, gg_colour, ntots, nod, ncolours  &
                              , nodof, nel_pp, neqs_pp, nn_pp, ltol, limit  &
                              , iters, error, precon)
  END SELECT
  RETURN
END SUBROUTINE LinearSolve

!-------------------------------
! Multifield Gather operation
!-------------------------------
SUBROUTINE GATHER_M(x, pmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)
  USE precision;
  USE global_variables;
  USE gather_scatter;
  USE MP_INTERFACE;
  USE Parallel_supplementary_Maths;
  IMPLICIT NONE
  INTEGER,    INTENT(IN)   :: ntots, nodof, nod, nel, neqs_pp, nn_pp
  INTEGER,    INTENT(IN)   :: MASK(nodof,nod)
  REAL(iwp),  INTENT(IN)   :: x(neqs_pp)
  REAL(iwp),  INTENT(INOUT):: pmul(ntots,nel)

  pmul = 0._iwp;
  CALL GATHERM(x, pmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)
  RETURN
END SUBROUTINE GATHER_M

!-------------------------------
! Multifield Gather operation
!-------------------------------
SUBROUTINE SCATTER_M(x, pmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)
  USE precision;
  USE global_variables;
  USE gather_scatter;
  USE MP_INTERFACE;
  USE Parallel_supplementary_Maths;
  IMPLICIT NONE
  INTEGER,    INTENT(IN)   :: ntots, nodof, nod, nel, neqs_pp, nn_pp
  INTEGER,    INTENT(IN)   :: MASK(nodof,nod)
  REAL(iwp),  INTENT(IN)   :: pmul(ntots,nel)
  REAL(iwp),  INTENT(INOUT):: x(neqs_pp);

  x = 0._iwp;
  CALL SCATTERM(x, pmul, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp)
  RETURN
END SUBROUTINE SCATTER_M

!-------------------------------
! Multifield matrix-vector multiplication
!-------------------------------
SUBROUTINE PARA_MATVEC(storA,x,b,MASK,nel,nn_pp,neqs_pp,ntots,nod,nodof)
  USE precision;
  USE global_variables;
  USE gather_scatter;
  USE MP_INTERFACE;
  USE Parallel_supplementary_Maths;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nel,ntots,nod,nodof,nn_pp,neqs_pp
  INTEGER,   INTENT(IN)   :: MASK(nodof,nod)
  REAL(iwp), INTENT(IN)   :: storA(ntots,ntots,nel), x(neqs_pp)
  REAL(iwp), INTENT(INOUT):: b(neqs_pp)
  REAL(iwp), ALLOCATABLE  :: pmul(:,:), qmul(:,:);

  ALLOCATE(pmul(ntots,nel), qmul(ntots,nel))
  pmul = 0._iwp;
  qmul = 0._iwp;
  b    = 0._iwp;
  CALL PARAMATVEC(storA,x,b,pmul,qmul,MASK,nel,nn_pp,neqs_pp,ntots,nod,nodof)
  DEALLOCATE(pmul, qmul)
  RETURN
END SUBROUTINE PARA_MATVEC

!-------------------------------
! Increment vector
!-------------------------------
SUBROUTINE INCREMENT(unew, uold, du, a, b, n_pp)
  USE precision;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: n_pp
  REAL(iwp), INTENT(IN)   :: a, b
  REAL(iwp), INTENT(IN)   :: uold(n_pp), du(n_pp)
  REAL(iwp), INTENT(INOUT):: unew(n_pp)

  unew = a*uold + b*du;
  RETURN
END SUBROUTINE INCREMENT

!-------------------------------
! Dot_product vector
!-------------------------------
SUBROUTINE DOT_PRODUCT_pp(dproduct, u_pp, v_pp, neqs_pp)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: neqs_pp;
  REAL(iwp), INTENT(IN)   :: u_pp(neqs_pp), v_pp(neqs_pp)
  REAL(iwp), INTENT(INOUT):: dproduct

  dproduct = DOT_PRODUCT_P(u_pp, v_pp)
  RETURN
END SUBROUTINE DOT_PRODUCT_pp

!-------------------------------
! Vector divided by scalar
!-------------------------------
SUBROUTINE Scalar_Vector_Product(unew, a, uold, n_pp)
  USE precision;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: n_pp;
  REAL(iwp), INTENT(IN)   :: uold(n_pp), a
  REAL(iwp), INTENT(INOUT):: unew(n_pp)

  unew = a*uold;
  RETURN
END SUBROUTINE Scalar_Vector_Product

!-------------------------------
! Copies a vector to another vector
!-------------------------------
SUBROUTINE Copy_Vector(unew, uold, n_pp)
  USE precision;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: n_pp;
  REAL(iwp), INTENT(IN)   :: uold(n_pp);
  REAL(iwp), INTENT(INOUT):: unew(n_pp);

  unew = uold
  RETURN
END SUBROUTINE Copy_Vector

!-------------------------------
! Parallel vector 2-norm
!-------------------------------
SUBROUTINE norm_pp(twonorm, u_pp, neqs_pp)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: neqs_pp;
  REAL(iwp), INTENT(IN)   :: u_pp(neqs_pp)
  REAL(iwp), INTENT(INOUT):: twonorm

  twonorm = norm_p(u_pp)
  RETURN
END SUBROUTINE norm_pp

!-------------------------------
! Parallel scalar global sum
!-------------------------------
SUBROUTINE SUM_PP(Total, partial)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  REAL(iwp), INTENT(IN)   :: partial
  REAL(iwp), INTENT(INOUT):: Total
  CALL MPI_ALLREDUCE(partial,Total,1,MPI_REAL8,MPI_SUM &
                   , MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE SUM_PP

!-------------------------------
! Parallel Stop condition
!-------------------------------
SUBROUTINE STOP_COND_PP(IsConverged, error, rtol)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER,   INTENT(INOUT):: IsConverged
  REAL(iwp), INTENT(IN)   :: error, rtol

  IsConverged = 0;
  IF(error < rtol) IsConverged = 1;
  CALL MPI_ALLREDUCE(IsConverged,IsConverged,1,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE STOP_COND_PP

!-------------------------------
! Parallel Integer maxval
!-------------------------------
SUBROUTINE INTEGER_MAXVAL_PP(globalMax, localMax)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER, INTENT(IN)   :: localMax
  INTEGER, INTENT(INOUT):: globalMax
  CALL MPI_ALLREDUCE(localMax,globalMax,1,MPI_INTEGER,MPI_MAX &
                   , MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE INTEGER_MAXVAL_PP

!-------------------------------
! Parallel REAL8/double maxval
!-------------------------------
SUBROUTINE REAL8_MAXVAL_PP(globalMax, localMax)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  REAL(iwp), INTENT(IN)   :: localMax
  REAL(iwp), INTENT(INOUT):: globalMax
  CALL MPI_ALLREDUCE(localMax,globalMax,1,MPI_REAL8,MPI_MAX &
                   , MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE REAL8_MAXVAL_PP

!-------------------------------
! Parallel REAL8/double maxval of a vector
!-------------------------------
SUBROUTINE REAL8_MAXVALVEC_PP(globalMax, Realvector, NEQs_pp)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: NEQs_pp;
  REAL(iwp), INTENT(IN)   :: Realvector(NEQs_pp);
  REAL(iwp), INTENT(INOUT):: globalMax;
  REAL(iwp)               :: localMax;

  localMax = MAXVAL(Realvector)
  CALL MPI_ALLREDUCE(localMax,globalMax,1,MPI_REAL8,MPI_MAX &
                   , MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE REAL8_MAXVALVEC_PP

!-------------------------------
! Parallel REAL8/double Absolute maxval of a vector
!-------------------------------
SUBROUTINE REAL8_ABSMAXVALVEC_PP(globalMax, Realvector, NEQs_pp)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: NEQs_pp;
  REAL(iwp), INTENT(IN)   :: Realvector(NEQs_pp);
  REAL(iwp), INTENT(INOUT):: globalMax;
  REAL(iwp)               :: localMax;

  localMax = MAXVAL(DABS(Realvector(:)))
  CALL MPI_ALLREDUCE(localMax,globalMax,1,MPI_REAL8,MPI_MAX &
                   , MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE REAL8_ABSMAXVALVEC_PP

!------------------------------------------------------------------------------
!
!                     Parallel Writing/Output functions
!
!------------------------------------------------------------------------------
!-------------------------------
! Alya ENSI Gold geometry output functions
!-------------------------------
SUBROUTINE ENSI_GEO_OUTPUT(argv, element, nlen, gcoord_pp, gnum_pp, numpes &
                         , npess, ndim, nod, nn, nel_pp, nn_pp)
  USE Parallel_IO;   USE precision;
  USE new_library;   USE MP_INTERFACE;
  USE gather_scatter;
  IMPLICIT NONE
  CHARACTER(LEN=15), INTENT(IN)   :: element
  CHARACTER(LEN=50), INTENT(IN)   :: argv;
  INTEGER,           INTENT(IN)   :: numpes,npess,nlen;
  INTEGER,           INTENT(IN)   :: ndim, nod, nn, nel_pp, nn_pp;
  INTEGER,           INTENT(IN)   :: gnum_pp(nod,nel_pp);
  REAL(iwp),         INTENT(IN)   :: gcoord_pp(nod,ndim,nel_pp)

  CALL ENSI_GEOMETRY_OUTPUT(argv, element, nlen, gcoord_pp, gnum_pp, numpes &
                          , npess, ndim, nod, nn, nel_pp, nn_pp);
  RETURN
END SUBROUTINE ENSI_GEO_OUTPUT

!-------------------------------
! Alya ENSI Gold Output element colour
!-------------------------------
SUBROUTINE ENSI_ELMCOLOUR_OUTPUT(argv,gg_colour,numpes,npess,nlen,nel_pp &
                                ,ndim,nod,element)
  USE Parallel_IO;   USE precision;
  USE new_library;   USE MP_INTERFACE;
  USE gather_scatter;
  IMPLICIT NONE
  CHARACTER(LEN=15), INTENT(IN)   :: element;   
  CHARACTER(LEN=50), INTENT(IN)   :: argv;
  INTEGER,           INTENT(IN)   :: numpes,npess,nlen
  INTEGER,           INTENT(IN)   :: nel_pp,ndim,nod;
  INTEGER,           INTENT(IN)   :: gg_colour(nel_pp);
  INTEGER,           PARAMETER    :: j=0,var=2;
  REAL(iwp),         ALLOCATABLE  :: gg_colourR(:,:);
  INTEGER                         :: Iel, I, K, L, nel2;

  IF(nod/=27) L=1;
  IF(nod==27) L=8; 
  nel2 = L*nel_pp
  ALLOCATE(gg_colourR(1,nel2))
  DO Iel=1,nel_pp
    DO I=1,L
      K = L*(Iel-1) + I;
      gg_colourR(1,K) = REAL(gg_colour(Iel),iwp)
    ENDDO
  ENDDO
  CALL ENSI_FILE_OUTPUT(argv,gg_colourR,numpes,npess,nlen,nel2,nel2 &
                       ,1,ndim,j,var,nod,element)
  DEALLOCATE(gg_colourR)
  RETURN
END SUBROUTINE ENSI_ELMCOLOUR_OUTPUT

!-------------------------------
! Alya ENSI Gold Real data output functions
!-------------------------------
SUBROUTINE ENSI_DATA_OUTPUT(argv,xnew_pp,numpes,npess,nlen,neq_p &
                           ,nn_pp,nodof,ndim,j,var,nod,element)
  USE Parallel_IO;   USE precision;
  USE new_library;   USE MP_INTERFACE;
  USE gather_scatter;
  IMPLICIT NONE
   CHARACTER(LEN=15), INTENT(IN)   :: element;
   CHARACTER(LEN=50), INTENT(IN)   :: argv;
   INTEGER,           INTENT(IN)   :: numpes,npess,nlen,j,var;
   INTEGER,           INTENT(IN)   :: neq_p,nodof,ndim,nn_pp,nod;
   REAL(iwp),         INTENT(IN)   :: xnew_pp(neq_p);


  CALL ENSI_FILE_OUTPUT(argv,xnew_pp,numpes,npess,nlen,neq_p,nn_pp &
                       ,nodof,ndim,j,var,nod,element)
  RETURN
END SUBROUTINE ENSI_DATA_OUTPUT

!-------------------------------
! Alya ENSI Gold Output partition
!-------------------------------
SUBROUTINE ENSI_PARTITION_DATA(argv,numpes,npess,nlen &
                              ,nn_pp,ndim)
  USE Parallel_IO;   USE precision;
  USE new_library;   USE MP_INTERFACE;
  USE gather_scatter;
  IMPLICIT NONE
   CHARACTER(LEN=15)               :: element = " ";
   CHARACTER(LEN=50), INTENT(IN)   :: argv;
   INTEGER,           INTENT(IN)   :: numpes,npess,nlen;
   INTEGER,           INTENT(IN)   :: ndim,nn_pp;
   INTEGER,           PARAMETER    :: nodof=1, j=1, var=1, nod=1;
   REAL(iwp),        ALLOCATABLE   :: xnew_pp(:);

  ALLOCATE(xnew_pp(nn_pp))
  xnew_pp = REAL(numpes,iwp)
  CALL ENSI_FILE_OUTPUT(argv,xnew_pp,numpes,npess,nlen,nn_pp,nn_pp &
                       ,nodof,ndim,j,var,nod,element)
  DEALLOCATE(xnew_pp)
  RETURN
END SUBROUTINE ENSI_PARTITION_DATA

!-------------------------------
! Alya ENSI Gold Real Traction
!-------------------------------
SUBROUTINE ENSI_TRACTION_DATA(element, gg_Face, numpes, npess, ndim  &
                            , nFace, nodFace, nod, nel_pp, nn_pp)
  USE Parallel_IO;   USE precision;
  USE new_library;   USE MP_INTERFACE;
  USE gather_scatter;
  USE Parallel_BoundaryConditions;
  IMPLICIT NONE
   INTEGER                         :: Iel, I, K, L;
   INTEGER                         :: nlen=23;
   CHARACTER(LEN=50)               :: argv='Results/TractionSurface';
   CHARACTER(LEN=15)               :: face_element;
   CHARACTER(LEN=15), INTENT(IN)   :: element;
   INTEGER,           INTENT(IN)   :: numpes,npess;
   INTEGER,           INTENT(IN)   :: ndim, nFace, nod, nodFace, nel_pp, nn_pp;
   INTEGER,           INTENT(IN)   :: gg_Face(nFace+1,nel_pp);
   INTEGER,           PARAMETER    :: nodof=1, j=1, var=1;
   INTEGER,          ALLOCATABLE   :: EF_MAP(:,:);
   REAL(iwp),        ALLOCATABLE   :: utemp_pp(:,:), xnew_pp(:);
   REAL(iwp),         PARAMETER    :: zero=0._iwp, one=1._iwp,tol=1.0E-4_iwp;

  ALLOCATE(xnew_pp(nn_pp), utemp_pp(nod,nel_pp),EF_MAP(nFace,nodFace))
  CALL Element_FaceMASK(EF_MAP, element, face_element, nface, nodFace, ndim)
  utemp_pp = zero; xnew_pp = zero;
  DO Iel=1,nel_pp
    DO I=1,nFace
      IF(gg_Face(I+1,Iel) /= 0)THEN
        DO L=1,nodFace
          K = EF_MAP(I,L);
          utemp_pp(K,Iel) = one;
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  CALL SCATTER(xnew_pp, utemp_pp)
  DO I=1,nn_pp
    IF(xnew_pp(I) < tol) xnew_pp(I) = zero;
    IF(xnew_pp(I) > tol) xnew_pp(I) = one;
  ENDDO
  CALL ENSI_FILE_OUTPUT(argv,xnew_pp,numpes,npess,nlen,nn_pp,nn_pp &
                       ,nodof,ndim,j,var,nod,element)
  DEALLOCATE(xnew_pp, utemp_pp, EF_MAP)
  RETURN
END SUBROUTINE ENSI_TRACTION_DATA

!-------------------------------
! finalise parafem mpi processes
!-------------------------------
SUBROUTINE finalise()
   USE MP_Interface;
   IMPLICIT NONE
   CALL SHUTDOWN()
   RETURN
END SUBROUTINE finalise

!------------------------------------------------------------------------------
!
!                     Parallel PROBLEM Element integration functions
!
!------------------------------------------------------------------------------
!-------------------------------
! Heat equation Element nodal mask
!-------------------------------
SUBROUTINE CalcHeatMasks(HeatMask, nod, nodof)
  IMPLICIT NONE
  INTEGER               :: i, j, k;
  INTEGER, INTENT(IN)   :: nod, nodof;
  INTEGER, INTENT(INOUT):: HeatMask(nodof,nod)

  DO i = 1,nod
    HeatMask(1,i) = i;
  ENDDO
  RETURN
END SUBROUTINE CalcHeatMasks


!-------------------------------
! Solid mechanics Element nodal mask
!-------------------------------
SUBROUTINE CalcSolidsMasks(SolidsMask, nod, nodof)
  IMPLICIT NONE
  INTEGER               :: i, j, k, l;
  INTEGER, INTENT(IN)   :: nod, nodof;
  INTEGER, INTENT(INOUT):: SolidsMask(nodof,nod);
  l=0;
  DO I = 1,nodof
    DO J = 1,nod
      l = l + 1;
      SolidsMask(I,J) = l;
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE CalcSolidsMasks

!-------------------------------
! Heat element Integration
!-------------------------------
SUBROUTINE HEAT_Integration(StorKA, StorKB, element, coord, utemp, diff, fibre &
                          , theta, dtim, nel_pp, ndim, nod, ntots, nip)
  USE precision;
  USE new_library;
  USE Orthotropic_Heat_ALE;
  IMPLICIT NONE
  INTEGER                     :: iel, i, j, k, l;
  CHARACTER(LEN=15),INTENT(IN):: element
  INTEGER,          INTENT(IN):: nel_pp, ndim, nod, ntots, nip;
  REAL(iwp),        INTENT(IN):: theta, dtim;
  REAL(iwp),        INTENT(IN):: coord(nod,ndim,nel_pp), utemp(ndim*nod,nels_pp);
  REAL(iwp),        INTENT(IN):: diff(ndim,nel_pp), fibre(ndim,ndim,nel_pp);
  REAL(iwp),     INTENT(INOUT):: StorKA(ntots,ntots,nel_pp)
  REAL(iwp),     INTENT(INOUT):: StorKB(ntots,ntots,nel_pp)
  REAL(iwp),     ALLOCATABLE  :: points(:,:), weights(:)
  INTEGER                     :: nip2;


  nip2 = nip;
  ALLOCATE(points(nip2,ndim), weights(nip2))
  StorKA = 0._iwp;  points  = 0._iwp;
  StorKB = 0._iwp;  weights = 0._iwp;
  CALL SAMPLE2(element, points, weights)
  ELEMENTS:DO iel = 1,nel_pp
      CALL Heat_Orthotropic_ALE_Element(StorKA(:,:,iel), StorKB(:,:,iel)    &
                             , coord(:,:,iel), utemp(:,iel), diff(:,iel)    &
                             , fibre(:,:,iel), theta, dtim, points, weights &
                             , ndim, nod, ntots, nip2)
  ENDDO ELEMENTS
  DEALLOCATE(points, weights)
  RETURN
END SUBROUTINE HEAT_Integration

!-------------------------------
! Solid mechanics element Residual-Jacobian Integration
! for Displacement-Pressure element
!-------------------------------
SUBROUTINE SOLID_Integration_UP(Residual, StoreKE, utemp     &
                           , coord, gg_pp, gg_Face, val_pp, Stress &
                           , MATPROP, nel_pp, ntots, ndim, nst, nip,nodMax  &
                           , nodof, nFace, nodFace, nipFace, nloadedFace    &
                           , nr, nprop, material, element)
  USE precision;
  USE new_library;
  USE Solids_UPAS;
  USE Static_SolidsTraction;
  USE Parallel_BoundaryConditions;
  USE Parallel_supplementary_Maths;
  USE PRECONDITIONERS;
  IMPLICIT NONE
  INTEGER                     :: iel, i, j, k, l, m, n, igauss, nip2, nip3;
  CHARACTER(LEN=15),INTENT(IN):: element;
  INTEGER,          INTENT(IN):: nprop, material, nodMax;
  INTEGER,          INTENT(IN):: nel_pp, ntots, ndim, nst, nip, nodof;
  INTEGER,          INTENT(IN):: nFace, nodFace, nipFace, nloadedFace, nr;
  INTEGER,          INTENT(IN):: gg_pp(ndim*nodMax,nel_pp)
  INTEGER,          INTENT(IN):: gg_Face(nFace+1,nel_pp);
  REAL(iwp),        INTENT(IN):: MATPROP(nprop), utemp(ntots,nel_pp)
  REAL(iwp),        INTENT(IN):: val_pp(ndim*nodMax,nel_pp), coord(nodMax,ndim,nel_pp);
  REAL(iwp),        INTENT(IN):: Stress(nst*nodFace,nloadedFace);
  REAL(iwp),     INTENT(INOUT):: Residual(ntots,nel_pp)
  REAL(iwp),     INTENT(INOUT):: StoreKE(ntots,ntots,nel_pp);
  REAL(iwp),         PARAMETER:: one = 1._iwp, zero = 0._iwp;
  REAL(iwp),       ALLOCATABLE:: centre(:), points(:,:), weights(:);   !Full integration
  INTEGER                     :: ndofU, ndofP, nodU, nodP, nodofP, nodofU;
  REAL(iwp)                   :: cx, cy, cr, ctol, rad0, rad1; !Rotation test


  ! Though the maps in this first component can technically be stored
  ! and reused in the context again later of the shear number of elements
  ! being integrated this cost of recalculating them only once per every nel_pp
  ! elements makes the cost irrelevant, more time can be saved in optimisations
  ! on the individual element level

  !---
  ! Element DOF and integration sizing routines
  !---
  nodofU = ndim;
  nodofP = 1;
  nodU = nodMax;
  ndofU=nodofU*nodU
  ndofP=nodofP*nodP


  !---
  !Array allocations
  !---
  ALLOCATE(centre(ndim), points(nip,ndim), weights(nip))


  !---
  !Numerical Integration and element mapping routines
  !---
  CALL SAMPLE2(element, points, weights)

  !---
  ! Element integration routines
  !---
  CALL STATIC_SOLIDUPAS_ELEMENT(Km, Rm, utemp, astrain, fibre, MATPROP   &
                              , coord, Ni_p, dNi_u, Ni_u, weights, nprop &
                              , ntots, ndofU, ndofP, ndim , nst, nip, nodU, nodp)

  !--
  !Apply Traction boundary conditions using surface elements
  !--
  IF(nloadedFace /= 0)THEN
  !  CALL SOLID_Traction_element(Residual, StoreKE, Element, coord, Stress, utemp &
  !                            , disp_MAP, gg_Face, nloadedFace, ndim, nodMax     &
  !                            , ntots, ndofU, nodU, nst, nFace, nodFace, nipFace &
  !                            , nip2, nel_pp)
  ENDIF

  !--
  !Apply Dirchelet boundary conditions using elimination scheme
  !--
  DO IEL = 1,nel_pp
    DO J = 1,nodU
      DO I = 1,ndim
        IF(gg_pp(M,iel)==0)THEN
          Residual(M,IEL)  = zero;
          StoreKE(M,:,IEL) = zero;
          StoreKE(M,M,IEL) = one;
        ENDIF
      ENDDO
    ENDDO
  ENDDO


  DEALLOCATE(points, weights, centre)
  RETURN
END SUBROUTINE SOLID_Integration_UP


!------------------------------------------------------------------------------
!
!                Data separation, generation and handling functions
!
!------------------------------------------------------------------------------
!-------------------------------
! Quick Pressure Boundary
!-------------------------------
SUBROUTINE Quick_Pressure(Stress, pressure, ntotsStress, ndim, nloadedFace)
  USE precision;
  INTEGER                  :: i, j, k, l, nst, nodFace;
  INTEGER,   INTENT(IN)    :: ntotsStress, ndim, nloadedFace;
  REAL(iwp), INTENT(IN)    :: pressure;
  REAL(iwp), INTENT(INOUT) :: Stress(ntotsStress,nloadedFace);

  nst     = ((ndim+1)*ndim)/2;
  nodFace = ntotsStress/nst;
  IF(nloadedFace > 0)THEN
    Stress = 0._iwp
    DO i = 1,nloadedFace
      DO k = 1,nodFace
        DO j = 1,nst
          l = (k-1)*nst + j;
          IF(j <= ndim) Stress(l,i) = pressure;
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE Quick_Pressure

!-------------------------------
! Split up fibre and diffusion 
!-------------------------------
SUBROUTINE Fibre_Diffusion(Fibre, Diff, data_pp, dof, ndim, nels_pp)
  USE precision;
  IMPLICIT NONE
  INTEGER                  :: iel, i, j, k, l;
  INTEGER,   INTENT(IN)    :: dof, ndim, nels_pp;
  REAL(iwp), INTENT(IN)    :: data_pp(dof,nels_pp);
  REAL(iwp), INTENT(INOUT) :: Fibre(ndim,ndim,nels_pp), Diff(ndim,nels_pp);

  DO iel = 1,nels_pp
    k = 1;
    DO i = 1,ndim
      DO j = 1,ndim
        fibre(j,i,iel) = data_pp(k,iel);
        k = k + 1;
      ENDDO
    ENDDO
    DO i = 1,ndim
      Diff(i,iel) = data_pp(k,iel);
      k = k + 1;
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE Fibre_Diffusion


!-------------------------------
! Split up displacement and pressure Global nodal
!-------------------------------
SUBROUTINE SPLIT_U_P(u_pp, pressure, data_pp, nodof, nn_pp, neqs_pp)
  USE precision;
  IMPLICIT NONE
  INTEGER                  :: iel, i, j, k, l;
  INTEGER,   INTENT(IN)    :: nodof, nn_pp, neqs_pp;
  REAL(iwp), INTENT(IN)    :: data_pp(neqs_pp);
  REAL(iwp), INTENT(INOUT) :: u_pp(nodof*nn_pp), pressure(nn_pp);

  u_pp(:)     = data_pp(1:(nodof*nn_pp));
  pressure(:) = data_pp((nodof*nn_pp + 1):);
  RETURN
END SUBROUTINE SPLIT_U_P

!-------------------------------
! Split up CellIDs and stimulus definitions 
!-------------------------------
SUBROUTINE CISD(CellID, StimDef, data_pp, dof, ndim, nn_pp)
  IMPLICIT NONE
  INTEGER                :: iel, i, j, k, l;
  INTEGER, INTENT(IN)    :: dof, ndim, nn_pp;
  INTEGER, INTENT(IN)    :: data_pp(dof,nn_pp);
  INTEGER, INTENT(INOUT) :: CellID(nn_pp), StimDef(nn_pp);

  CellID(:)  = data_pp(1,:);
  StimDef(:) = data_pp(2,:);
  RETURN
END SUBROUTINE CISD
