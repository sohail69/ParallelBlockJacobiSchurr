MODULE Parallel_IO
  USE precision;      USE global_variables;
  USE maths;          USE gather_scatter;
  USE new_library;    USE MP_INTERFACE;
!
! My own libraries
!
  USE Parallel_supplementary_Maths;
  CONTAINS
!------------------------------------------------------------------------------
!
!                     Parallel Reading/Input functions
!
!------------------------------------------------------------------------------

!-------------------------------
! Read in Job data
!-------------------------------
  SUBROUTINE READ_JOB_DATA(job_name, nlen, numpe, element, mesh, partitioner    &
                         , np_types, nels, nn, nod, nip, nbnd, nodFace, nipFace &
                         , nFace, npri)
    IMPLICIT NONE
    INTEGER,   INTENT(IN)           :: numpe, nlen;
    CHARACTER(LEN=50), INTENT(IN)   :: job_name
    CHARACTER(LEN=15), INTENT(INOUT):: element;
    INTEGER,   INTENT(INOUT)        :: mesh, partitioner, np_types;
    INTEGER,   INTENT(INOUT)        :: nels, nn, nod, nip;
    INTEGER,   INTENT(INOUT)        :: nbnd, nFace, nipFace, nodFace, npri;
    INTEGER                         :: bufsize,ier,integer_store(12)
    CHARACTER(LEN=50)               :: fname
  
!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------
    integer_store      = 0       !Initialising store arrays
    IF (numpe==1) THEN
      fname = job_name(1:INDEX(job_name, " ") -1) // ".dat"
      OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
      READ(10,*) element, mesh, partitioner, np_types, nels, nn, nod, nip  &
               , nbnd, nFace, nodFace, nipFace, npri
      CLOSE(10)
      integer_store(1)   = mesh
      integer_store(2)   = partitioner
      integer_store(3)   = np_types
      integer_store(4)   = nels
      integer_store(5)   = nn
      integer_store(6)   = nod
      integer_store(7)   = nip
      integer_store(8)   = nbnd
      integer_store(9)   = nFace
      integer_store(10)  = nodFace
      integer_store(11)  = nipFace
      integer_store(12)  = npri
    END IF

  bufsize = 12;
  CALL MPI_BCAST(integer_store,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier);

  bufsize = 15;
  CALL MPI_BCAST(element, bufsize, MPI_CHARACTER, 0, MPI_COMM_WORLD, ier);
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);

!------------------------------------------------------------------------------
! 4. Slave processors extract the variables from the temporary arrays
!------------------------------------------------------------------------------
  IF (numpe/=1) THEN
    mesh        = integer_store(1);
    partitioner = integer_store(2);
    np_types    = integer_store(3);
    nels        = integer_store(4);
    nn          = integer_store(5);
    nod         = integer_store(6);
    nip         = integer_store(7);
    nbnd        = integer_store(8);
    nFace       = integer_store(9)
    nodFace     = integer_store(10);
    nipFace     = integer_store(11);
    npri        = integer_store(12);
  END IF

  RETURN
  END SUBROUTINE READ_JOB_DATA
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!-------------------------------
! Read in material data
!-------------------------------
  SUBROUTINE READ_MATERIAL_DATA(job_name, numpe, matprops, nmat, np_types) 
  IMPLICIT NONE
  INTEGER,           INTENT(IN)   :: numpe, nmat, np_types
  CHARACTER(LEN=50), INTENT(IN)   :: job_name
  REAL(iwp),         INTENT(INOUT):: matprops(nmat,np_types)
!------------------------------------------------------------------------------
! 1. Local variables
!------------------------------------------------------------------------------
  INTEGER                         :: bufsize,ier, i, k;
  REAL(iwp)                       :: real_store(nmat)
  CHARACTER(LEN=50)               :: fname
  CHARACTER(LEN=50)               :: program_name

!------------------------------------------------------------------------------
! 2. Master processor reads the data and copies it into temporary arrays
!------------------------------------------------------------------------------
  real_store         = 0.0_iwp
  matprops           = 0.0_iwp;

  IF (numpe==1) THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".mat"
    OPEN(10,FILE=fname,STATUS='OLD',ACTION='READ')
    READ(10,*) !Headers
    READ(10,*) !Headers
    DO i = 1,np_types
       READ(10,*) k, matprops(:,i)
    ENDDO
    CLOSE(10)
  END IF

!------------------------------------------------------------------------------
! 3. Master processor broadcasts the temporary arrays to the slave processors
!------------------------------------------------------------------------------
  bufsize = nmat
  DO i = 1,np_types
    IF (numpe==1) THEN
      real_store = matprops(:,i)
    END IF
    bufsize = nmat
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
    CALL MPI_BCAST(real_store,bufsize,MPI_REAL8,0,MPI_COMM_WORLD,ier)
    IF (numpe /= 1) THEN
      matprops(:,i) = real_store;
    END IF
  ENDDO
  RETURN
  END SUBROUTINE READ_MATERIAL_DATA
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_BOUNDARY(JOB_NAME,NUMPE1,boundary_N,nodof,nbnd)
  IMPLICIT NONE  
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  CHARACTER(LEN=50)             :: fname
  INTEGER                       :: i,ier
  INTEGER,INTENT(IN)            :: numpe1, nodof, nbnd;
  INTEGER,INTENT(INOUT)         :: boundary_N(nbnd,nodof+1)

  boundary_N = 1;
  IF(numpe1==1)THEN
    fname = job_name(1:INDEX(job_name, " ")-1) // ".bnd"
    OPEN(23, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,nbnd
      READ(23,*) boundary_N(i,:)  !First number is node Number
    END DO
    CLOSE(23)
  END IF
  IF(npes > 1)THEN
    bufsize = nbnd*(nodof+1)
    CALL MPI_BCAST(boundary_N,bufsize,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  ENDIF
  RETURN
  END SUBROUTINE READ_BOUNDARY
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

  SUBROUTINE READ_DATA_REAL(JOB_NAME,NUMPE1,npes1,data_pp,n_pp,dof)
  IMPLICIT NONE
  CHARACTER(LEN=50)             :: fname
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  REAL(iwp), INTENT(INOUT)      :: data_pp(dof,n_pp);
  INTEGER,   INTENT(IN)         :: npes1, numpe1, n_pp, dof;
  INTEGER                       :: i,j,bufsize,ier
  INTEGER                       :: STATUS(MPI_STATUS_SIZE);
  INTEGER,   ALLOCATABLE        :: n_p(:), n_p1(:);
  REAL(iwp), ALLOCATABLE        :: readdat(:,:)

  ALLOCATE(n_p(npes1), n_p1(npes1))
  n_p     = 0;
  n_p1    = 0;
  data_pp = 0._iwp;
  n_p(numpe1)  = n_pp;
  n_p1(numpe1) = n_pp;
  bufsize = npes1;
  CALL MPI_ALLREDUCE(n_p1,n_p,bufsize,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  DEALLOCATE(n_p1)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);

  IF(numpe1 == 1)THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat";
    OPEN(11, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,npes
      ALLOCATE(readdat(dof,n_p(i)))
      DO j = 1,n_p(i)
        READ(11,*) readdat(:,j)
      ENDDO
      IF(numpe1 == i)THEN
        data_pp = readdat;
      ELSE
        bufsize = dof*n_p(i)
        CALL MPI_SEND(readdat,bufsize,MPI_REAL8,i-1,i &
                     ,MPI_COMM_WORLD,status,ier)
      ENDIF
      DEALLOCATE(readdat)
    ENDDO
    CLOSE(11)
    DEALLOCATE(n_p)
  ENDIF

  IF(numpe1 /= 1)THEN
    DEALLOCATE(n_p)
    bufsize = dof*n_pp
    CALL MPI_RECV(data_pp,bufsize,MPI_REAL8,0,numpe1 &
                , MPI_COMM_WORLD,status,ier)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  RETURN
  END SUBROUTINE READ_DATA_REAL
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!


  SUBROUTINE READ_DATA_INTEGER(JOB_NAME,NUMPE1,npes1,data_pp,n_pp,dof)
  IMPLICIT NONE
  CHARACTER(LEN=50)             :: fname
  CHARACTER(LEN=50), INTENT(IN) :: job_name
  INTEGER,   INTENT(INOUT)      :: data_pp(dof,n_pp);
  INTEGER,   INTENT(IN)         :: npes1, numpe1, n_pp, dof;
  INTEGER                       :: i,j,bufsize,ier
  INTEGER                       :: STATUS(MPI_STATUS_SIZE);
  INTEGER,   ALLOCATABLE        :: n_p(:), n_p1(:);
  INTEGER,   ALLOCATABLE        :: readdat(:,:)

  ALLOCATE(n_p(npes1), n_p1(npes1))
  n_p     = 0;
  n_p1    = 0;
  data_pp = 0._iwp;
  n_p(numpe1)  = n_pp;
  n_p1(numpe1) = n_pp;
  bufsize = npes1;
  CALL MPI_ALLREDUCE(n_p1,n_p,bufsize,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,ier)
  DEALLOCATE(n_p1)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);

  IF(numpe1 == 1)THEN
    fname = job_name(1:INDEX(job_name, " ") -1) // ".dat";
    OPEN(11, FILE=fname, STATUS='OLD', ACTION='READ')
    DO i = 1,npes
      ALLOCATE(readdat(dof,n_p(i)))
      DO j = 1,n_p(i)
        READ(11,*) readdat(:,j)
      ENDDO
      IF(numpe1 == i)THEN
        data_pp = readdat;
      ELSE
        bufsize = dof*n_p(i)
        CALL MPI_SEND(readdat,bufsize,MPI_INTEGER,i-1,i &
                     ,MPI_COMM_WORLD,status,ier)
      ENDIF
      DEALLOCATE(readdat)
    ENDDO
    CLOSE(11)
    DEALLOCATE(n_p)
  ENDIF

  IF(numpe1 /= 1)THEN
    DEALLOCATE(n_p)
    bufsize = dof*n_pp
    CALL MPI_RECV(data_pp,bufsize,MPI_INTEGER,0,numpe1 &
                , MPI_COMM_WORLD,status,ier)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  RETURN
  END SUBROUTINE READ_DATA_INTEGER


!------------------------------------------------------------------------------
!
!                     Parallel Writing/Output functions
!
!------------------------------------------------------------------------------

!-------------------------------------------------------------------!
!                      Output ENSI Geometry-File                    !
!-------------------------------------------------------------------!
SUBROUTINE ENSI_GEOMETRY_OUTPUT(argv, element, nlen, gcoord_pp, gnum_pp, numpes &
                              , npess, ndim, nod, nn, nel_pp, nn_pp);

   IMPLICIT NONE
   INTEGER                         :: iel, i, j, k, l, m, nels, nels2, nod2;
   INTEGER                         :: bufsize1, bufsize2, request;
   INTEGER                         :: statu(MPI_STATUS_SIZE)
   CHARACTER(LEN=15)               :: EnsiElement;
   CHARACTER(LEN=15), INTENT(IN)   :: element;
   CHARACTER(LEN=50), INTENT(IN)   :: argv;
   INTEGER,           INTENT(IN)   :: numpes,npess,nlen;
   INTEGER,           INTENT(IN)   :: ndim, nod, nn, nel_pp, nn_pp;
   INTEGER,           INTENT(IN)   :: gnum_pp(nod,nel_pp);
   REAL(iwp),         INTENT(IN)   :: gcoord_pp(nod,ndim,nel_pp)
   REAL(iwp),         PARAMETER    :: one = 1._iwp, zero = 0._iwp;
   REAL(iwp),         ALLOCATABLE  :: x_pp(:), gcoordn_pp(:,:);
   REAL(iwp),         ALLOCATABLE  :: element_constant(:,:), nodal_constant(:);
   INTEGER,           ALLOCATABLE  :: nels_array(:), nels_array2(:), nn_array(:); ! Sizing arrays
   INTEGER,           ALLOCATABLE  :: size_data(:), gtemp(:,:);   ! Data packets
   INTEGER,           ALLOCATABLE  :: ENSI_NODAL_MASK(:), gtemp2(:,:)


   CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
   !-----------
   ! Sizing array calcs and allocations
   !-----------
   ALLOCATE(size_data(2))
   IF(numpes==1)THEN
     ALLOCATE(nels_array(npess), nn_array(npess), nels_array2(npess));
     nn_array(1)   = nn_pp;
     nels_array(1) = nel_pp;
   ENDIF

  IF(npess > 1)THEN
    request = 0;
    IF(numpes/=1)THEN
      size_data(1) = nn_pp;
      size_data(2) = nel_pp;
      bufsize1 = 2;
      CALL MPI_ISEND(size_data,bufsize1,MPI_INTEGER,0,numpes &
                    ,MPI_COMM_WORLD,request,ier)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
    IF(numpes==1)THEN
      PROCESSES1:DO i = 2,npess
        bufsize1 = 2;
        size_data = 0;
        CALL MPI_RECV(size_data,bufsize1,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
        nn_array(i)   = size_data(1);
        nels_array(i) = size_data(2);
       ENDDO PROCESSES1
     ENDIF
     IF(numpes/=1) CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
   ENDIF
   DEALLOCATE(size_data);

   !-----------
   ! Writeout File headers
   !-----------
  IF(numpes==1)THEN
    OPEN(12,FILE=argv(1:nlen)//".ensi.geo",STATUS='REPLACE',ACTION='WRITE')
    WRITE(12,'(/2A)')   "Problem name: ", argv(1:nlen)
    WRITE(12,'(A/A/A)') "Geometry files","node id given","element id given"
    WRITE(12,'(A/A)')   "part","      1"
    WRITE(12,'(A)') "Volume Mesh"
    WRITE(12,'(A)') "coordinates"
    WRITE(12,'(I10)') nn
  ENDIF

  !-----------
  ! Element to nodal Coordinates
  !-----------
  ALLOCATE(element_constant(nod,nel_pp),nodal_constant(nn_pp))
  element_constant = one;
  nodal_constant   = zero;
  CALL SCATTER(nodal_constant,element_constant)
  DEALLOCATE(element_constant)


  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  ALLOCATE(gcoordn_pp(ndim,nn_pp), x_pp(nn_pp));
  gcoordn_pp(:,:) = zero;
  DIMENSION1:DO i = 1,ndim
    x_pp = zero;
    CALL SCATTER(x_pp,gcoord_pp(:,i,:))
    gcoordn_pp(i,:) = x_pp/nodal_constant;
  ENDDO DIMENSION1
  DEALLOCATE(x_pp, nodal_constant)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);


  !-----------
  ! Send and writeout Nodal Coordinates
  !-----------
  DIMENSIONS:DO i = 1,ndim
    IF(numpes==1)THEN
      DO k =1,nn_pp
        WRITE(12,'(E12.5)') gcoordn_pp(i,k)
      ENDDO
    ENDIF

    IF(npess > 1)THEN
      request = 0;
      IF(numpes/=1)THEN
        ALLOCATE(x_pp(nn_pp))
        bufsize1 = nn_pp;
        x_pp = gcoordn_pp(i,:);
        CALL MPI_ISEND(x_pp,bufsize1,MPI_REAL8,0,numpes,MPI_COMM_WORLD,request,ier)
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
      IF(numpes==1)THEN
        DO j = 2,npess
          ALLOCATE(x_pp(nn_array(j)))
          x_pp = zero;
          bufsize1 = nn_array(j)
          CALL MPI_RECV(x_pp,bufsize1,MPI_REAL8,j-1,j,MPI_COMM_WORLD,statu,ier)
          DO k =1,nn_array(j)
            WRITE(12,'(E12.5)') x_pp(k)
          ENDDO
          DEALLOCATE(x_pp)
        ENDDO
      ENDIF
      IF(numpes/=1) CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
      IF(numpes/=1) DEALLOCATE(x_pp)
    ENDIF
  ENDDO DIMENSIONS


  !Additional zeros as ENSI requires 3-D coordinates
  IF(ndim < 3)THEN
    IF(numpes == 1)THEN
      DO i = 1,(3-ndim)
        DO j = 1,nn
          WRITE(12,'(E12.5)') zero;
        ENDDO
      ENDDO
    ENDIF
  ENDIF

   !-----------
   ! Send and writeout Element connectivity/steering
   !-----------
   IF(nod == 27)THEN
     nod2                      = 8
     nels2                     = 8*nel_pp
     IF(numpes==1) nels_array2 = 8*nels_array
     ALLOCATE(gtemp2(nod2,nels2))
     DO iel = 1,nel_pp
       k = 8*(iel-1);
       gtemp2(:,k+1) = gnum_pp((/6,14,26,13,18,22,27,23/),iel)
       gtemp2(:,k+2) = gnum_pp((/14,7,15,26,22,19,24,27/),iel)
       gtemp2(:,k+3) = gnum_pp((/18,22,27,23,2,10,25,9/),iel)
       gtemp2(:,k+4) = gnum_pp((/22,19,24,27,10,3,11,25/),iel)
       gtemp2(:,k+5) = gnum_pp((/13,26,16,5,23,27,21,17/),iel)
       gtemp2(:,k+6) = gnum_pp((/26,15,8,16,27,24,20,21/),iel)
       gtemp2(:,k+7) = gnum_pp((/23,27,21,17,9,25,12,1/),iel)
       gtemp2(:,k+8) = gnum_pp((/27,24,20,21,25,11,4,12/),iel)
     ENDDO
   ELSE IF(nod == 9)THEN
     nod2                      = 4
     nels2                     = 4*nel_pp
     IF(numpes==1) nels_array2 = 4*nels_array
     ALLOCATE(gtemp2(nod2,nels2))
     DO iel = 1,nel_pp
       k = 4*(iel-1);
       gtemp2(:,k+1) = gnum_pp((/1,2,9,8/),iel)
       gtemp2(:,k+2) = gnum_pp((/2,3,4,9/),iel)
       gtemp2(:,k+3) = gnum_pp((/8,9,6,7/),iel)
       gtemp2(:,k+4) = gnum_pp((/9,4,5,6/),iel)
     ENDDO
   ELSE
     nod2  = nod
     nels2 = nel_pp
     IF(numpes==1) nels_array2 = nels_array
     ALLOCATE(gtemp2(nod2,nels2))
     gtemp2 = gnum_pp
   ENDIF


   ALLOCATE(ENSI_NODAL_MASK(nod2))
   CALL ENSI_ELEMENT_NODAL_NUMBERING(ENSI_NODAL_MASK, element, nod2)
   IF(numpes==1)THEN
     CALL Element_to_ENSIName(EnsiElement, element, nod2)
     nels = SUM(nels_array2(:));
     WRITE(12,'(A/I10)') EnsiElement(1:INDEX(EnsiElement, " ") -1), nels
     ELEMENTS1:DO j =1,nels_array2(1)
       WRITE(12,*) gtemp2(ENSI_NODAL_MASK(:),j)
     ENDDO ELEMENTS1
   ENDIF


  
  IF(npess > 1)THEN
    request = 0;
    IF(numpes/=1)THEN
      bufsize1 = nod2*nels2
      CALL MPI_ISEND(gtemp2,bufsize1,MPI_INTEGER,0,numpes,MPI_COMM_WORLD,request,ier)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
    IF(numpes==1)THEN
      PROCESSES3:DO i = 2,npess
        ALLOCATE( gtemp(nod2,nels_array2(i)) )
        bufsize2 = nod2*nels_array2(i)
        CALL MPI_RECV(gtemp,bufsize2,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
        ELEMENTS2:DO j =1,nels_array2(i)
          WRITE(12,*) gtemp(ENSI_NODAL_MASK(:),j)
        ENDDO ELEMENTS2
        DEALLOCATE(gtemp)
      ENDDO PROCESSES3
    ENDIF
    IF(numpes/=1) CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
  ENDIF
  

  DEALLOCATE(ENSI_NODAL_MASK, gtemp2)
  IF(numpes==1) CLOSE(12)
  IF(numpes==1) DEALLOCATE(nels_array, nn_array, nels_array2);
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE ENSI_GEOMETRY_OUTPUT


!-------------------------------------------------------------------!
!                   Output nodal Real ENSI file data                !
!-------------------------------------------------------------------!
SUBROUTINE ENSI_FILE_OUTPUT(argv,xnew_pp,numpes,npess,nlen,neq_p,nn_pp &
                           ,nodof,ndim,j,var,nod,element)
   USE OUTPUT;
   IMPLICIT NONE
   INTEGER                         :: i, m, n;
   CHARACTER(LEN=6)                :: ch;
   CHARACTER(LEN=15)               :: variableType, varRef, rFrame;
   CHARACTER(LEN=15), INTENT(IN)   :: element;   
   CHARACTER(LEN=50), INTENT(IN)   :: argv;
   INTEGER,           INTENT(IN)   :: numpes,npess,nlen,j,var;
   INTEGER,           INTENT(IN)   :: neq_p,nodof,ndim,nn_pp,nod;
   REAL(iwp),         INTENT(IN)   :: xnew_pp(neq_p);
   REAL(iwp),         ALLOCATABLE  :: ttr_pp(:);
   REAL(iwp),         PARAMETER    :: one = 1._iwp, zero = 0._iwp;

   CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
   !-----------
   ! Element or nodal reference frame
   !-----------
   IF(var==1) THEN 
     rFrame    = "coordinates"
     varRef = "per-node";
   ENDIF

   IF(var==2) THEN
     IF(nod/=27) CALL Element_to_ENSIName(rFrame, element, nod)
     IF(nod==27) CALL Element_to_ENSIName(rFrame, element, 8)
     varRef = "per-element";
   ENDIF

   IF((var/=2).AND.(var/=1)) THEN
     rFrame    = "coordinates"
     varRef = "per-node";
   ENDIF

   !-----------
   ! Variable type
   !-----------
   IF(nodof==1)                   variableType = "Scalar";
   IF(nodof==ndim)                variableType = "Vector";
   IF(nodof==((ndim*(ndim+1))/2)) variableType = "Tensor";


   ALLOCATE(ttr_pp(nn_pp))
   IF(numpes==1)THEN
     WRITE(ch,'(I6.6)') j
     OPEN(12,FILE=argv(1:nlen)//".ensi.DISPL-"//ch,STATUS='REPLACE',ACTION='WRITE')
   ENDIF

   IF(numpes==1) WRITE(12,'(A)')"Alya Ensight Gold ---"// &
                                variableType(1:INDEX(variableType," ")-1)// &
                                " "// varRef(1:INDEX(varRef," ")-1) // &
                                " variable file"

   IF(numpes==1) WRITE(12,'(A/A/A)') "part", "    1", rFrame(1:INDEX(rFrame," ")-1)

   !-----------
   !Scalar Case
   !-----------
   IF(nodof==1)THEN
     ttr_pp=xnew_pp;
     CALL dismsh_ensi_p(12,1,nn_pp,npess,numpes,1,ttr_pp)
   ENDIF

   !-----------
   !Vector case
   !-----------
   IF(nodof==ndim)THEN
     DOFS1:DO i = 1,nodof
       m = (i - 1)*nn_pp + 1;
       n = i*nn_pp;
       ttr_pp=xnew_pp(m:n);
       CALL dismsh_ensi_p(12,1,nn_pp,npess,numpes,1,ttr_pp)
     ENDDO DOFS1
     IF(ndim < 3)THEN
       DOFS2:DO i = 1,(3-nodof)
         ttr_pp=zero;
         CALL dismsh_ensi_p(12,1,nn_pp,npess,numpes,1,ttr_pp)
       ENDDO DOFS2
     ENDIF
   ENDIF

   !-----------
   !Tensor Case
   !-----------
   IF(nodof==((ndim*(ndim+1))/2))THEN
     DOFS3:DO i = 1,nodof
       m = (i - 1)*nn_pp + 1;
       n = i*nn_pp;
       ttr_pp=xnew_pp(m:n);
       CALL dismsh_ensi_p(12,1,nn_pp,npess,numpes,1,ttr_pp)
     ENDDO DOFS3
     IF(ndim < 3)THEN
       DOFS4:DO i = 1,(6-nodof)
         ttr_pp=zero;
         CALL dismsh_ensi_p(12,1,nn_pp,npess,numpes,1,ttr_pp)
       ENDDO DOFS4
     ENDIF
   ENDIF
   DEALLOCATE(ttr_pp)
   IF(numpes==1) CLOSE(12)
   CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
   RETURN
END SUBROUTINE ENSI_FILE_OUTPUT


!-------------------------------------------------------------------!
!                   parafem to ENSI element type                    !
!-------------------------------------------------------------------!
SUBROUTINE Element_to_ENSIName(EnsiElement, element, nod)
  IMPLICIT NONE
  INTEGER                         :: i, j, k, l, m;
  CHARACTER(LEN=6)                :: ch;
  CHARACTER(LEN=15), INTENT(INOUT):: EnsiElement;
  CHARACTER(LEN=15), INTENT(IN)   :: element;
  INTEGER,           INTENT(IN)   :: nod;


  WRITE(ch,'(I6.1)') nod
  j = 0;
  l = nod;
  DO k = 1,6
    l = l/10;
    j = j + 1;
    IF(l == 0)EXIT
  ENDDO
  m = 10**j
  IF(MODULO(j,m) == 0) i = j;
  IF(MODULO(j,m) /= 0) i = j - 1;
  i = 6 - i
  SELECT CASE(element)
    CASE('line')
      EnsiElement = "bar"  // ch(i:6);
    CASE('triangle')
      EnsiElement = "tria"  // ch(i:6);
    CASE('quadrilateral')
      EnsiElement = "quad"  // ch(i:6);
    CASE('tetrahedron')
      EnsiElement = "tetra" // ch(i:6);
    CASE('hexahedron')
      EnsiElement = "hexa"  // ch(i:6);
  END SELECT
  RETURN
END SUBROUTINE Element_to_ENSIName

!-------------------------------------------------------------------!
!           parafem to ENSI element nodal numbering                 !
!-------------------------------------------------------------------!
SUBROUTINE ENSI_ELEMENT_NODAL_NUMBERING(ENSI_node_numbering, element, nod)
  IMPLICIT NONE
  INTEGER                         :: i, j, k, l, m;
  CHARACTER(LEN=15), INTENT(IN)   :: element;
  INTEGER,           INTENT(IN)   :: nod;
  INTEGER,           INTENT(INOUT):: ENSI_node_numbering(nod);

  ENSI_node_numbering = 1;
  SELECT CASE(element)
    CASE('line')
      DO i = 1,nod
        ENSI_node_numbering(i) = i;
      ENDDO
    CASE('triangle')
      IF(nod == 3)  ENSI_node_numbering(:) = (/3, 2, 1/);
    CASE('quadrilateral')
      IF(nod == 4)  ENSI_node_numbering(:) = (/1, 4, 3, 2/);
      IF(nod == 8)  ENSI_node_numbering(:) = (/1, 7, 5, 3, 8, 6, 4, 2/);
      IF(nod == 9)  ENSI_node_numbering(:) = (/1, 7, 5, 3, 8, 6, 4, 2, 9/);
    CASE('tetrahedron')
      IF(nod == 4)  ENSI_node_numbering(:)=(/1,3,2,4/);
   !  IF(nod == 10) ENSI_node_numbering(:)=(/1,3,2,4,5,6,7,8,9,10/);
      IF(nod == 10) ENSI_node_numbering(:)=(/1,2,3,4,5,6,7,8,9,10/);
    CASE('hexahedron')
      IF(nod == 8)  ENSI_node_numbering(:)=(/1,4,8,5,2,3,7,6/);
      IF(nod == 20) ENSI_node_numbering(:)=(/4,3,7,8,1,2,6,5,11,19,15 &
                                            ,20,9,18,13,17,12,10,14,16/);

      !IF(nod == 27) ENSI_node_numbering(1:20)=(/4,3,7,8,1,2,6,5,11,19,15 &
      !                                         ,20,9,18,13,17,12,10,14,16/);
      IF(nod == 27) ENSI_node_numbering(1:20)=(/1,2,3,4,5,6,7,8,9,10,11,12 &
                                               ,13,14,15,16,17,18,19,20/);

    CASE DEFAULT
      WRITE(*,*) "ERROR undefined element"
      ENSI_node_numbering(:) = 1;
  END SELECT
  RETURN
END SUBROUTINE ENSI_ELEMENT_NODAL_NUMBERING
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!

!-------------------------------------------------------------------!
!                   Outputs a nodal array for data                  !
!                    1-D time series simulations                    !
!-------------------------------------------------------------------!
SUBROUTINE GLOBAL_1D_ARRAY(u_global, u_pp, numpes, npess, nn, nn_pp);
   IMPLICIT NONE
   INTEGER                         :: iel, i, j, k, l, m;
   INTEGER                         :: bufsize1, request;
   INTEGER                         :: statu(MPI_STATUS_SIZE)
   INTEGER,           INTENT(IN)   :: numpes,npess;
   INTEGER,           INTENT(IN)   :: nn, nn_pp;
   REAL(iwp),         INTENT(IN)   :: u_pp(nn_pp);
   REAL(iwp),         INTENT(INOUT):: u_global(nn)
   REAL(iwp),         PARAMETER    :: one = 1._iwp, zero = 0._iwp;
   REAL(iwp),         ALLOCATABLE  :: x_pp(:);
   INTEGER,           ALLOCATABLE  :: nn_array(:);   ! Sizing arrays
   INTEGER                         :: temp

   CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
   !-----------
   ! Sizing array calcs and allocations
   !-----------
   u_global = zero;
   IF(numpes==1)THEN
     ALLOCATE(nn_array(npess));
     nn_array(1)   = nn_pp;
   ENDIF

   IF(npess > 1)THEN
     request = 0;
     IF(numpes/=1)THEN
       bufsize1 = 1;
       temp = nn_pp;
       CALL MPI_ISEND(temp,bufsize1,MPI_INTEGER,0,numpes,MPI_COMM_WORLD,request,ier)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
     IF(numpes==1)THEN
       PROCESSES1:DO i = 2,npess
         bufsize1 = 1;
         CALL MPI_RECV(temp,bufsize1,MPI_INTEGER,i-1,i,MPI_COMM_WORLD,statu,ier)
         nn_array(i) = temp
        ENDDO PROCESSES1
     ENDIF
     IF(numpes/=1) CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
   ENDIF

  !-----------
  ! Send and accumulate nodal data on numpe=1
  !-----------
  IF(numpes==1)THEN
    k = 1;
    l = nn_pp
    u_global(k:l) = u_pp
  ENDIF
      
  IF(npess > 1)THEN
    request = 0;
    IF(numpes/=1)THEN
      bufsize1 = nn_pp;
      CALL MPI_ISEND(u_pp,bufsize1,MPI_REAL8,0,numpes,MPI_COMM_WORLD,request,ier)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
    IF(numpes==1)THEN
      k = 1;
      DO j = 2,npess
        ALLOCATE(x_pp(nn_array(j)))
        x_pp = zero;
        bufsize1 = nn_array(j)
        CALL MPI_RECV(x_pp,bufsize1,MPI_REAL8,j-1,j,MPI_COMM_WORLD,statu,ier)
        k = k + nn_array(j-1)
        l = k + nn_array(j)
        u_global(k:l) = x_pp
        DEALLOCATE(x_pp)
      ENDDO
    ENDIF
    IF(numpes/=1) CALL MPI_WAIT(request, MPI_STATUS_IGNORE, ier)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  RETURN
ENDSUBROUTINE GLOBAL_1D_ARRAY
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
!-------------------------------------------------------------------!
END MODULE Parallel_IO
