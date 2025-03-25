MODULE BlockJacobiPrecon
  INTEGER, PARAMETER :: iwp = SELECTED_REAL_KIND(15,300)
  INCLUDE "mpif.h"
  CONTAINS
!================================================
!================================================
!================================================

!/***************************\
! Parafem Thread Partitioning
!
! Starts numbering from 0
!
!\***************************/
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


!/***************************\
! PARAFEMS classic Partitioner
!
! Starts numbering from 1
!
!\***************************/
SUBROUTINE PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nprocs, procID)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nn, nprocs, procID
  INTEGER,   INTENT(INOUT):: nn_pp1, nn_pp2, nn_pp, num

  nn_pp2 = 0;
  nn_pp1 = 0;
  nn_pp  = 0;
  num    = 0;

  nn_pp2 = nn/nprocs
  num    = nn - nn_pp2*nprocs
  nn_pp1 = nn_pp2
  IF(num/=0) nn_pp1 = nn_pp1 + 1
  IF((procID <= num).OR.(num == 0))THEN
    nn_pp = nn_pp1;
  ELSE
    nn_pp = nn_pp2;
  ENDIF
  RETURN
ENDSUBROUTINE PF_PARTITIONER1

!/***************************\
! Calculate the process that 
! owns a particular node/element
!
! Starts numbering from 1
!\***************************/
PURE FUNCTION FIND_NODE_PROC1(GnodeID,nn_pp1,nn_pp2,num,nprocs) RESULT(procID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: GnodeID, nprocs, num, nn_pp1, nn_pp2;
  INTEGER            :: PID1, A, B;
  INTEGER            :: procID

  A = GnodeID - num*nn_pp1
  IF(A == 0) procID = num;
  IF(A <  0)THEN
    PID1 = GnodeID/nn_pp1
    B = GnodeID - PID1*nn_pp1;
    IF(B == 0) procID = PID1
    IF(B /= 0) procID = PID1 + 1;
  ENDIF
  IF(A >  0)THEN
    PID1 = A/nn_pp2;
    B = A - PID1*nn_pp2
    IF(B == 0) procID = PID1 + num
    IF(B /= 0) procID = PID1 + num + 1;
  ENDIF
ENDFUNCTION FIND_NODE_PROC1

!/***************************\
! Calculates the local 
! numberingID from global
! numberingID of a particular
! node/element
!\***************************/

PURE FUNCTION GLOBAL_TO_LOCAL1(G_nodeID,nn_pp1,nn_pp2,num,nprocs) RESULT(L_nodeID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: G_nodeID, num, nn_pp1, nn_pp2, nprocs;
  INTEGER            :: L_nodeID, procID, threshhold;

  threshhold = num*nn_pp1;
  procID = FIND_NODE_PROC1(G_nodeID,nn_pp1,nn_pp2,num,nprocs)
  procID = procID - 1;
  L_nodeID = G_nodeID - procID*nn_pp1
  IF(num /= 0)THEN
    IF(G_nodeID > threshhold)THEN
      L_nodeID = G_nodeID - num*nn_pp1 - (procID - num)*nn_pp2
    ENDIF
  ENDIF
ENDFUNCTION GLOBAL_TO_LOCAL1

!/***************************\
!  Find the first position
!  instance of an integer in
!  an integer array
!\***************************/
FUNCTION FIND_POS1(a, array, nsize) RESULT(I)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: nsize
  INTEGER, INTENT(IN):: a, array(nsize)
  INTEGER            :: I, J

  I = -1;
  DO J = 1,nsize
    IF( a == array(J) ) I = J;
  ENDDO
ENDFUNCTION FIND_POS1


!/***************************\
!  Size up the message tables
!  for local processor using
!  the element steering array
!\***************************/
SUBROUTINE SIZE_MESSAGE_TABLES_FE(nSends, nRecvs, gg_pp, LprocID &
                                , nProcs, nn, nel, nod, nel_pp)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L, procID, errMPI;
  INTEGER,   INTENT(IN)   :: LprocID, nProcs;
  INTEGER,   INTENT(IN)   :: nn, nel, nod, nel_pp;
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel_pp);
  INTEGER,   INTENT(INOUT):: nSends, nRecvs
  INTEGER,   ALLOCATABLE  :: procTable(:), procTable_tmp(:);
  INTEGER                 :: nn_pp1, nn_pp2, nn_pp, num;

  !
  ! Do a paritioning of
  ! the problem
  !
  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nProcs, LprocID+1)

  !
  ! Allocate and initialise
  !
  ALLOCATE(procTable(nProcs), procTable_tmp(nProcs))
  nSends = 0;
  nRecvs = 0;
  procTable = 0;
  procTable_tmp = 0;

  !
  ! Calculate the number
  ! of sends by going through
  ! all of the element nodes
  !
  DO Iel = 1,nel_pp
    DO I = 1,nod
      procID = FIND_NODE_PROC1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nProcs);
      IF(procID /= (LprocID+1)) procTable_tmp(procID) = 1;
    ENDDO
  ENDDO

  DO J = 1,nProcs
    IF(procTable_tmp(J) /= 0) nSends = nSends + 1;
  ENDDO

  !
  ! Calculate the number of
  ! recieves using MPI_ALLREDUCE
  ! to find the locally owned
  ! remote element nodes
  !
  CALL MPI_ALLREDUCE(procTable_tmp,procTable,nProcs,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,errMPI);
  nRecvs = procTable(LprocID+1);
  DEALLOCATE(procTable, procTable_tmp)
  RETURN
ENDSUBROUTINE SIZE_MESSAGE_TABLES_FE


!/***************************\
!  Forms the message tables used
!  in send-recieve to exchange
!  nodes data between the processes
!  repeats a lot of calculation
!  from the table sizing
!\***************************/
SUBROUTINE FORM_MESSAGE_TABLES(sendTables, recvTables, sendNodes &
                             , recvNodes, gg_pp, nRecvs, nSends  &
                             , nn , nod, nel_pp, nn_ppMax, nProcs&
                             , LprocID)
  IMPLICIT NONE
  INTEGER               :: Iel, I, J, K, L, procID, LNodeID, errMPI;
  INTEGER, INTENT(IN)   :: nRecvs, nSends, nProcs, LprocID
  INTEGER, INTENT(IN)   :: nn, nod, nel_pp, nn_ppMax;
  INTEGER, INTENT(IN)   :: gg_pp(nod,nel_pp);
  INTEGER, INTENT(INOUT):: sendTables(2,nSends), recvTables(2,nRecvs)
  INTEGER, INTENT(INOUT):: sendNodes(nSends,nn_ppMax), recvNodes(nRecvs,nn_ppMax)
  INTEGER               :: nn_pp1, nn_pp2, nn_pp, num;
  INTEGER               :: nSends_tmp;

  !
  ! Do a paritioning of
  ! the problem
  !
  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nProcs, LprocID)

  !
  ! Mark all the node that
  ! belong to a non-local 
  ! process and all the
  ! procIDs of neighbouring
  ! processes
  !
  nSends_tmp = 0;
  sendTables = -1;
  DO Iel = 1,nel_pp
    DO I = 1,nod
      procID = FIND_NODE_PROC1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nProcs)-1;
      IF(procID /= LprocID)THEN
        LNodeID = GLOBAL_TO_LOCAL1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nprocs)
        L = -1;
	PROC_SEARCH:DO J = 1,nSends_tmp
          IF(sendTables(1,J) == procID) sendNodes(J,LNodeID) = 1;
          IF(sendTables(1,J) == procID) L=1;
          IF(sendTables(1,J) == procID) EXIT PROC_SEARCH;
        ENDDO PROC_SEARCH
        IF(L == -1) nSends_tmp = nSends_tmp + 1;
        IF(L == -1) sendTables(1,nSends_tmp) = procID;
        IF(L == -1) sendNodes(nSends_tmp,LNodeID) = 1;
      ENDIF
    ENDDO
  ENDDO

  !
  ! The number of messages
  ! is known however the recv
  ! processes are unknown
  ! so send and recv the procIDs
  !
  CALL SEND_RECV_DATAI(recvNodes, recvTables   &
                     , sendNodes, sendTables   &
                     , LprocID, nRecvs, nSends &
                     , nn_ppMax)

  !
  ! Number the non-zero recv 
  ! nodes chronologically to
  ! simplify matrix assembly
  !
  DO I = 1,nRecvs
    K=0;
    DO J = 1,nn_ppMax
      IF(recvNodes(I,J) /= 0) K = K + 1;
      IF(recvNodes(I,J) /= 0) recvNodes(I,J) = K;
    ENDDO
  ENDDO

  RETURN
ENDSUBROUTINE FORM_MESSAGE_TABLES

!/***************************\
!  Send and Recieve messages
!  of integer data arrays using
!  message tables
!\***************************/
SUBROUTINE SEND_RECV_DATAI(recv_array, recvTable   &
                         , send_array, sendTable   &
                         , LprocID, nRecvs, nSends &
                         , maxSize)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L;
  INTEGER                 :: pID, psize, maxMessages;
  INTEGER                 :: RecvTest, errMPI;
  INTEGER,   INTENT(IN)   :: LprocID, nRecvs, nSends, maxSize;
  INTEGER,   INTENT(IN)   :: sendTable(2,nSends)
  INTEGER,   INTENT(IN)   :: send_array(nSends,maxSize);
  INTEGER,   INTENT(INOUT):: recv_array(nRecvs,maxSize);
  INTEGER,   INTENT(INOUT):: recvTable(2,nRecvs);
  INTEGER,   ALLOCATABLE  :: vrequest(:), vstatus(:,:)
  LOGICAL                 :: flagger;

  maxMessages = MAX(nSends,nRecvs);
  ALLOCATE(vrequest(maxMessages), vstatus(MPI_STATUS_SIZE,maxMessages))
  vrequest = MPI_REQUEST_NULL;


  !
  ! Non-blocking sends
  ! the non-locally Owned 
  ! data
  !
  IF(nSends>0)THEN
    DO I = 1,nSends
      pID   = sendTable(1,I);
      psize = sendTable(2,I);
!! Writes out a messages for analysis message table analysis
!! WRITE(*,*) "THIS IS :", LprocID, "SENDING TO : ", pID
      CALL MPI_ISEND(send_array(I,:),maxSize,MPI_INTEGER,pID  &
                   , LprocID, MPI_COMM_WORLD,vrequest(I),errMPI)

    ENDDO
  ENDIF

  !
  ! Recieve the Messages
  !
  recvTable = -1; !initialise recv table
  IF(nRecvs > 0)THEN
    DO I = 1,nRecvs
      !
      ! Message probing loop
      ! finds unique messages
      !
      flagger = .TRUE.
      DO WHILE(flagger)
        CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD &
                      ,vstatus(1,I),errMPI)
        pID = vstatus(MPI_tag,I);
        RecvTest = 0;
        DO k = 1,I-1
          IF(recvTable(1,K) == pID) RecvTest = RecvTest + 1;
        ENDDO
        IF(RecvTest == 0) flagger = .FALSE.;
      ENDDO

      !
      ! Get expected message
      ! size
      !
      CALL MPI_GET_COUNT(vstatus(1,I),MPI_INTEGER,psize,errMPI)

      !
      ! Recieve the unique
      ! message and update
      ! the recieve log
      !
      recvTable(1,I) = pID;
      recvTable(2,I) = maxSize;
      CALL MPI_RECV(recv_array(I,:),maxSize,MPI_INTEGER  &    
                  , pID,MPI_ANY_TAG, MPI_COMM_WORLD,vstatus(1,I),errMPI)
    ENDDO
  ENDIF
  CALL MPI_WAITALL(maxMessages,vrequest,vstatus,errMPI)
  DEALLOCATE(vrequest, vstatus)
  RETURN
ENDSUBROUTINE SEND_RECV_DATAI


!/***************************\
!  Send and Recieve messages
!  of real data arrays using
!  message tables
!\***************************/
SUBROUTINE SEND_RECV_DATAR(recv_array, recvTable &
                         , send_array, sendTable &
                         , nRecvs, nSends, maxSize)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L;
  INTEGER                 :: pID, procID, psize, maxMessages;
  INTEGER                 :: RecvTest, errMPI;
  INTEGER,   INTENT(IN)   :: nRecvs, nSends, maxSize;
  INTEGER,   INTENT(IN)   :: sendTable(2,nSends)
  REAL(iwp), INTENT(IN)   :: send_array(nSends,maxSize);
  REAL(iwp), INTENT(INOUT):: recv_array(nRecvs,maxSize);
  INTEGER,   INTENT(INOUT):: recvTable(2,nRecvs);
  INTEGER,   ALLOCATABLE  :: vrequest(:), vstatus(:,:)
  LOGICAL                 :: flagger;

  maxMessages = MAX(nSends,nRecvs);
  ALLOCATE(vrequest(maxMessages), vstatus(MPI_STATUS_SIZE,maxMessages))
  vrequest = MPI_REQUEST_NULL;

  !
  ! Non-blocking sends
  ! the non-locally Owned 
  ! data
  !
  IF(nSends>0)THEN
    DO I = 1,nSends
      pID   = sendTable(1,I);
      psize = sendTable(2,I);
!! Writes out a messages for analysis message table analysis
!! WRITE(*,*) "THIS IS :", LprocID, "SENDING TO : ", pID
      CALL MPI_ISEND(send_array(I,1:psize),psize,MPI_REAL8,pID  &
                   , procID,MPI_COMM_WORLD,vrequest(I),errMPI)
    ENDDO
  ENDIF

  !
  ! Recieve the Messages
  !
  recvTable = -1; !initialise recv table
  IF(nRecvs > 0)THEN
    DO I = 1,nRecvs 
      !
      ! Message probing loop
      ! finds unique messages
      !
      flagger = .TRUE.
      DO WHILE(flagger)
        CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD &
                      ,vstatus(1,I),errMPI)
        pID = vstatus(MPI_tag,I);
        RecvTest = 0;
        DO k = 1,I-1
          IF(recvTable(1,K) == pID) RecvTest = RecvTest + 1;
        ENDDO
        IF(RecvTest == 0) flagger = .FALSE.;
      ENDDO

      !
      ! Get expected message
      ! size from local table;
      !
      CALL MPI_GET_COUNT(vstatus(1,I),MPI_REAL8,psize,errMPI)

      !
      ! Recieve the unique
      ! message and update
      ! the recieve log
      !
      recvTable(1,I) = pID;
      recvTable(2,I) = psize;
      CALL MPI_RECV(recv_array(I,1:psize),psize,MPI_REAL8 &    
                  , pID,MPI_ANY_TAG, MPI_COMM_WORLD,vstatus(1,I),errMPI)
    ENDDO
  ENDIF
  CALL MPI_WAITALL(maxMessages,vrequest,vstatus,errMPI)
  DEALLOCATE(vrequest, vstatus)
  RETURN
ENDSUBROUTINE SEND_RECV_DATAR

!/***************************\
!  Local dist assembled matrix
!  size, finds the size of
!  the data structures wh
!\***************************/
SUBROUTINE LOCAL_MAT_SIZE(Offsets, gg_pp, nod, nodof, nn, nel_pp)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L, P, Q;
  INTEGER                 :: pID, procID1, procID2, psize;
  INTEGER,   INTENT(IN)   :: nn, nodof, nel_pp, nod;
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel_pp), Offsets(maxMessages);



  RETURN
ENDSUBROUTINE LOCAL_MAT_SIZE

!/***************************\
!  Assemble the element
!  matrices into a global
!  extended distributed matrix
!\***************************/
SUBROUTINE ASSEMBLE_DIST_MATRIX_EXT(Amat_dist, AmatElm_pp, Offsets, gg_pp    &
                                  , Nei_procs, nRecvs, nSends, maxMessages   &
                                  , maxSize, matSize, LprocID, nProcs, nodof &
                                  , nel_pp, nod, ntot, nn)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L, P, Q;
  INTEGER                 :: pID, procID1, procID2, psize;
  INTEGER                 :: RecvTest, errMPI;
  INTEGER,   INTENT(IN)   :: nRecvs, nSends, maxMessages, maxSize, matSize, LprocID, nProcs;
  INTEGER,   INTENT(IN)   :: nn, nodof, nel_pp, ntot, nod;
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel_pp), Offsets(maxMessages), Nei_procs(maxMessages);
  REAL(iwp), INTENT(IN)   :: AmatElm_pp(ntot,ntot,nel_pp);
  REAL(iwp), INTENT(INOUT):: Amat_dist(matSize,matSize);
  INTEGER                 :: nn_pp1, nn_pp2, nn_pp, num;


  !
  ! Do a paritioning of
  ! the problem
  !
  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nProcs, LprocID)


  !
  ! Assemble the local 
  ! components of the matrix
  !
  DO Iel = 1,nel_pp
    DO I = 1,nod
      procID2 = FIND_NODE_PROC1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nProcs);
      P = I
      IF(procID1 /= LprocID) P = FIND_POS1(procID1, Nei_procs, maxMessages);
      DO J = 1,nod
        procID2 = FIND_NODE_PROC1(gg_pp(J,Iel),nn_pp1,nn_pp2,num,nProcs);
        Q = I
        IF(procID2 /= LprocID) Q = FIND_POS1(procID2, Nei_procs, maxMessages);
        DO K = 1,nodof
          DO L = 1,nodof
!Amat_dist
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO


!Amat_dist() = Amat_dist() + 
!IF(procID /= LprocID)THEN
!LNodeID = GLOBAL_TO_LOCAL1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nprocs)
  RETURN
ENDSUBROUTINE ASSEMBLE_DIST_MATRIX_EXT   


!/***************************\
!  Generate the Schurr
!  complement block jacobi
!  matrix inverse used for
!  preconditioning
!\***************************/
SUBROUTINE SCHURR_COMPLEMENT(Pmat,Amat,neqsA,neqsB)
  IMPLICIT NONE
  INTEGER                 :: I, J, K, L
  INTEGER,   INTENT(IN)   :: neqsA, neqsB
  REAL(iwp), INTENT(IN)   :: Amat(neqsA+neqsB,neqsA+neqsB)
  REAL(iwp), INTENT(INOUT):: Pmat(neqsA,neqsA)
  REAL(iwp), ALLOCATABLE  :: BlockMM(:,:), BlockMN(:,:), BlockNM(:,:)
  REAL(iwp), ALLOCATABLE  :: BlockNN(:,:), BlockNNInv(:,:)

  ALLOCATE(BlockMM(neqsA,neqsA), BlockMN(neqsA,neqsB), BlockNM(neqsB,neqsA))
  ALLOCATE(BlockNN(neqsB,neqsB), BlockNNInv(neqsB,neqsB))

  I = 1
  J = neqsA
  K = neqsA + 1
  L = neqsA + neqsB

  BlockMM = Amat(I:J,I:J)
  BlockMN = Amat(I:J,K:L)
  BlockNM = Amat(K:L,I:J)
  BlockNN = Amat(K:L,K:L)
!  CALL INVERT2(BlockNN,BlockNNInv,neqsB)
  Pmat = BlockMM - MATMUL(BlockMN,MATMUL(BlockNNInv,BlockNM))

  DEALLOCATE(BlockMM, BlockMN, BlockNM)
  DEALLOCATE(BlockNN, BlockNNInv)
  RETURN
END SUBROUTINE SCHURR_COMPLEMENT

!================================================
!================================================
!================================================
ENDMODULE BlockJacobiPrecon
