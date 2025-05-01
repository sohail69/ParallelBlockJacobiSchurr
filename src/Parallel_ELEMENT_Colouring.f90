MODULE GRAPH_COLOURING
  USE PRECISION;      USE global_variables;  USE MP_INTERFACE;
  USE maths;          USE gather_scatter;    USE new_library;
  USE PARALLEL_SUPPLEMENTARY_MATHS;
  CONTAINS


!---
! PARAFEMS classic Partitioner
!---
SUBROUTINE PF_PARTITIONER(nn_pp1, nn_pp2, nn_pp, num, nn, npess, numpes)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nn, npess, numpes
  INTEGER,   INTENT(INOUT):: nn_pp1, nn_pp2, nn_pp, num

  nn_pp2 = 0;
  nn_pp1 = 0;
  nn_pp  = 0;
  num    = 0;

  nn_pp2 = nn/npess
  num    = nn - nn_pp2*npess
  nn_pp1 = nn_pp2
  IF(num/=0) nn_pp1 = nn_pp1 + 1
  IF((numpes <= num).OR.(num == 0))THEN
    nn_pp = nn_pp1;
  ELSE
    nn_pp = nn_pp2;
  ENDIF
  RETURN
ENDSUBROUTINE PF_PARTITIONER


!---
! Calculate the process that owns
! a particular node/element
!---
PURE FUNCTION FIND_NODE_PROC(nodeID,nn_pp1,nn_pp2,num,npess) RESULT(procID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: nodeID, npess, num, nn_pp1, nn_pp2;
  INTEGER            :: PID1, A, B;
  INTEGER            :: procID

  A = nodeID - num*nn_pp1
  IF(A == 0) procID = num;
  IF(A <  0)THEN
    PID1 = nodeID/nn_pp1
    B = nodeID - PID1*nn_pp1;
    IF(B == 0) procID = PID1
    IF(B /= 0) procID = PID1 + 1;
  ENDIF
  IF(A >  0)THEN
    PID1 = A/nn_pp2;
    B = A - PID1*nn_pp2
    IF(B == 0) procID = PID1 + num
    IF(B /= 0) procID = PID1 + num + 1;
  ENDIF
ENDFUNCTION FIND_NODE_PROC


!---
! Calculates the local numberingID from global
! numberingID of a particular node/element
!---
PURE FUNCTION GLOBAL_TO_LOCAL(G_nodeID,nn_pp1,nn_pp2,num,npess) RESULT(L_nodeID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: G_nodeID, num, nn_pp1, nn_pp2, npess;
  INTEGER            :: L_nodeID, procID, threshhold;

  threshhold = num*nn_pp1;
  procID = FIND_NODE_PROC(G_nodeID,nn_pp1,nn_pp2,num,npess)
  procID = procID - 1;
  L_nodeID = G_nodeID - procID*nn_pp1
  IF(num /= 0)THEN
    IF(G_nodeID > threshhold)THEN
      L_nodeID = G_nodeID - num*nn_pp1 - (procID - num)*nn_pp2
    ENDIF
  ENDIF
ENDFUNCTION GLOBAL_TO_LOCAL


!---
! Calculates the  global numberingID from local
! numberingID of a particular node/element
!---
PURE FUNCTION LOCAL_TO_GLOBAL(L_nodeID,nn_pp1,nn_pp2,num,numpes) RESULT(G_nodeID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: L_nodeID, numpes, num, nn_pp1, nn_pp2;
  INTEGER            :: G_nodeID, procID;
 
  procID = numpes - 1;
  G_nodeID = L_nodeID + procID*nn_pp1;
  IF(num/=0)THEN
    IF(num < procID)THEN
      G_nodeID = L_nodeID + num*nn_pp1 + (procID-num)*nn_pp2;
    ENDIF
  ENDIF
END FUNCTION LOCAL_TO_GLOBAL

!---
! Calculate the maximal number
! of elements connected to a node
!---
SUBROUTINE MAXIMAL_NODAL_CONNECTIVITY(n_nel, nod, nn_pp, nel_pp, npess, numpes)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nod, nn_pp, nel_pp, npess, numpes
  INTEGER,   INTENT(INOUT):: n_nel
  REAL(iwp), ALLOCATABLE  :: pmul(:,:), p_pp(:)
  REAL(iwp), PARAMETER    :: one=1._iwp, zero=0._iwp;
  REAL(iwp)               :: x, x1;

  ALLOCATE(pmul(nod,nel_pp), p_pp(nn_pp))
  pmul = one;
  p_pp = zero;
  CALL SCATTER(p_pp,pmul)
  x = MAXVAL(p_pp)
  CALL MPI_ALLREDUCE(x,x1,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier);
  x1 = x1 + 0.3_iwp !Prevents excessive rounding down
  n_nel = INT(x1)
  DEALLOCATE(pmul, p_pp)
  RETURN
ENDSUBROUTINE MAXIMAL_NODAL_CONNECTIVITY

!---
! Calculates the approximate sends and recieves of a
! 2-D array partitioned by ParaFEM partitioning
! Only works for non-negative array Indexes 
! Must also be contiguous
!---
SUBROUTINE APPROX_PARTITIONED_COMM(Aproxsend2, sends, recvs, recvs2, gg_pp &
                                 , nelms, n, n1, n2, num, numpes, npess)
  IMPLICIT NONE
  INTEGER                :: Iel, I, J, K, L;
  INTEGER, INTENT(IN)    :: nelms, n, n1, n2, num, numpes, npess;
  INTEGER, INTENT(IN)    :: gg_pp(nelms,n)
  INTEGER, INTENT(INOUT) :: Aproxsend2(npess), sends, recvs, recvs2
  INTEGER, ALLOCATABLE   :: AproxRecv(:), Aproxsend(:), Element_procs(:);
  INTEGER                :: procID, nodeID


  ALLOCATE(AproxRecv(npess), Aproxsend(npess), Element_procs(nelms))
  !Size-up the send buffers
  Aproxsend = 0;
  DO Iel = 1,n
    Element_procs = 0;
    DO I = 1,nelms
      nodeID  = gg_pp(I,Iel)
      IF(nodeID/=0)THEN
        procID  = FIND_NODE_PROC(nodeID, n1, n2, num,npess)
        !New-element-proc check
        K = 0;
        DO J = 1,I
          IF(Element_procs(J)==procID) K = K + 1;
        ENDDO
        IF(K==0) Aproxsend(procID)=Aproxsend(procID)+1; !new proc for element
        Element_procs(I) = procID;
      ENDIF
    ENDDO
  ENDDO

  !Size-up the send and receive tables for the processes
  sends=0;
  recvs=0;
  AproxRecv  = 0;
  Aproxsend2 = 0;
  DO I = 1,npess
    IF(Aproxsend(I)/=0) Aproxsend2(I) = 1;
    IF(Aproxsend(I)/=0) sends = sends + 1; !messages to be sent
  ENDDO

  CALL MPI_ALLREDUCE(Aproxsend2,AproxRecv,npess,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  recvs = AproxRecv(numpes) !expected number of messages to be recieved
  recvs2 = recvs;

  !Don't send yourself a message thats just sad
  IF(Aproxsend2(numpes) /= 0)THEN
    recvs2 = recvs2 - 1;
    sends  = sends - 1;
  ENDIF
  DEALLOCATE(AproxRecv, Aproxsend, Element_procs)
  RETURN
END SUBROUTINE APPROX_PARTITIONED_COMM


!---
! Construct the node-element adjacency graph
!---
SUBROUTINE NODE_ELEMENT_ADJACENCYGRAPH(NODE_ADJ, gg_pp, n_nel, nod, nel_pp, nn_pp &
                                     , nels, nn, numpes, npess)
  IMPLICIT NONE
  INTEGER                 :: Iel, Iel1, I, J, K, L;
  INTEGER,   INTENT(IN)   :: n_nel, nod, nn_pp, nel_pp, nn, nels, numpes, npess
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel_pp);
  INTEGER,   INTENT(INOUT):: NODE_ADJ(n_nel,nn_pp); !Adjacency graph-array
  INTEGER                 :: Iterators(nn_pp); !Adjacency graph-Iterators
  INTEGER                 :: procID,procID2,nodeID,LnodeID
  INTEGER                 :: nel_pp1,nel_pp2,nume,nn_pp1,nn_pp2,numn
  INTEGER                 :: sends, recvs, recvs2, bufsiz;
  INTEGER                 :: RecvTest,lflag
  INTEGER                 :: vrequest2(npess),vstatus2(MPI_STATUS_SIZE,npess);
  INTEGER                 :: Element_procs(nod);
  INTEGER,   ALLOCATABLE  :: Sendbuf(:), Recvbuf(:), gg2_pp(:,:,:);
  INTEGER,   ALLOCATABLE  :: Aproxsend2(:);

  CALL PF_PARTITIONER(nn_pp1,  nn_pp2, K,  numn, nn,   npess, numpes) !Nodes
  CALL PF_PARTITIONER(nel_pp1, nel_pp2, K, nume, nels, npess, numpes) !Elements


  !Size-up the send buffers
  ALLOCATE(Aproxsend2(npess))
  CALL APPROX_PARTITIONED_COMM(Aproxsend2, sends, recvs, recvs2, gg_pp &
                             , nod, nel_pp, nn_pp1, nn_pp2, numn, numpes, npess)

  !List the procs to send and recieve data
  ALLOCATE(Sendbuf(sends), Recvbuf(recvs), gg2_pp(recvs,nod,nel_pp1))
  Recvbuf = -1
  gg2_pp = 0;
  IF(Aproxsend2(numpes) /= 0)THEN
    gg2_pp(recvs,:,1:nel_pp) = gg_pp
    Recvbuf(recvs) = numpes
  ENDIF
  K=0;
  DO I = 1,npess
    IF(I/=numpes)THEN
      IF(Aproxsend2(I)/=0) K=K+1;
      IF(Aproxsend2(I)/=0) Sendbuf(K) = I-1;
    ENDIF
  ENDDO


  !---
  !Send and receive the elements from the neighboring procs
  !---
  !Send and receive data
  bufsiz = nod*nel_pp
  DO I = 1,sends
    CALL MPI_ISEND(gg_pp,bufsiz,MPI_INTEGER,Sendbuf(I),numpes &
                  ,MPI_COMM_WORLD,vrequest2(I),ier)
    IF(ier/=MPI_SUCCESS) CALL MPERROR('Error in (A4) isend',ier)
    CALL MPI_TEST(vrequest2(I),lflag,vstatus2(1,I),ier)
  ENDDO

  DO I = 1,recvs2
    RecvTest = 1; !Initialise the receive test
    DO WHILE(RecvTest/=0)
      CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus2(1,I),ier)
      IF(ier /= MPI_SUCCESS) CALL MPERROR('Error in (A5) probe',ier)
      CALL MPI_GET_COUNT(vstatus2(1,I),MPI_INTEGER,bufsiz,ier)
      IF(ier /= MPI_SUCCESS) CALL MPERROR('Error in (A5) get_count',ier)
      procID = vstatus2(MPI_tag,I)
      RecvTest = 0;
      DO J = 1,I-1
        IF(Recvbuf(J)==procID) RecvTest=RecvTest+1;
      ENDDO
    ENDDO
    L = bufsiz/nod;
    Recvbuf(I) = procID;
    CALL MPI_RECV(gg2_pp(I,:,1:L),bufsiz,MPI_INTEGER,procID-1,MPI_ANY_TAG &
                 ,MPI_COMM_WORLD,vstatus2(1,I),ier)
  ENDDO
  CALL MPI_WAITALL(sends,vrequest2,vstatus2,ier)
  IF (ier /= MPI_SUCCESS) CALL MPERROR ('Error in MPI_WAITALL', ier)


  !---
  !Form node-element Adjacency array
  !---
  NODE_ADJ  = 0;
  Iterators = 0;
  DO Iel = 1,nel_pp1
    DO I = 1,nod
      DO J = 1,recvs
        K = gg2_pp(J,I,Iel)
        nodeID = GLOBAL_TO_LOCAL(K,nn_pp1,nn_pp2,numn,npess)
        procID = FIND_NODE_PROC(K,nn_pp1,nn_pp2,numn,npess)
        !Adjacency graph-array
        IF(procID==numpes)THEN
          Iterators(nodeID) = Iterators(nodeID) + 1;
          L = Iterators(nodeID)
          procID2 = Recvbuf(J)
          Iel1 = LOCAL_TO_GLOBAL(Iel,nel_pp1,nel_pp2,nume,procID2)
          NODE_ADJ(L,nodeID) = Iel1;
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE(gg2_pp, Sendbuf, Recvbuf)
  RETURN
ENDSUBROUTINE NODE_ELEMENT_ADJACENCYGRAPH


!---
! Construct the element-element adjacency graph
!---
SUBROUTINE ELEMENT_ELEMENT_ADJACENCYGRAPH(ELM_ADJ,NODE_ADJ, gg_pp, n_nel, nod &
                                        , nel_pp, nn_pp, nels, nn, numpes, npess)
  IMPLICIT NONE
  INTEGER                 :: Iel, Iel1, I, J, K, L, M, N;
  INTEGER,   INTENT(IN)   :: n_nel, nod, nn_pp, nel_pp, nn, nels, numpes, npess
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel_pp),  NODE_ADJ(n_nel,nn_pp);
  INTEGER,   INTENT(INOUT):: ELM_ADJ(nod*n_nel,nn_pp); 
  INTEGER                 :: Iterators(nn_pp); !Adjacency graph-Iterators
  INTEGER                 :: procID,nodeID, ElmID
  INTEGER                 :: nel_pp1,nel_pp2,nume,nn_pp1,nn_pp2,numn
  INTEGER                 :: sends, recvs, recvs2, bufsiz;
  INTEGER                 :: RecvTest,lflag,vrequest2(npess)
  INTEGER                 :: vstatus2(MPI_STATUS_SIZE,npess);
  INTEGER                 :: Element_procs(nod);
  INTEGER,   ALLOCATABLE  :: Sendbuf(:), Recvbuf(:), NODE_ADJ_tmp(:,:,:);
  INTEGER,   ALLOCATABLE  :: Aproxsend2(:);

  CALL PF_PARTITIONER(nn_pp1,  nn_pp2, K,  numn, nn,   npess, numpes) !Nodes
  CALL PF_PARTITIONER(nel_pp1, nel_pp2, K, nume, nels, npess, numpes) !Elements


  !Size-up the send buffers
  ALLOCATE(Aproxsend2(npess))
  CALL APPROX_PARTITIONED_COMM(Aproxsend2, sends, recvs, recvs2, NODE_ADJ, n_nel &
                             , nn_pp, nel_pp1, nel_pp2, nume, numpes, npess)

  !List the procs to send and recieve data
  ALLOCATE(Sendbuf(sends), Recvbuf(recvs), NODE_ADJ_tmp(recvs,n_nel,nn_pp1))
  Recvbuf = -1
  IF(Aproxsend2(numpes) /= 0)THEN
    Recvbuf(recvs) = numpes
    NODE_ADJ_tmp(recvs,:,1:nn_pp) = NODE_ADJ
  ENDIF
  K=0;
  DO I = 1,npess
    IF(I/=numpes)THEN
      IF(Aproxsend2(I)/=0) K=K+1;
      IF(Aproxsend2(I)/=0) Sendbuf(K) = I-1;
    ENDIF
  ENDDO

  !---
  !Send and receive the elements from the neighboring procs
  !---
  bufsiz = n_nel*nn_pp
  DO I = 1,sends
    CALL MPI_ISEND(NODE_ADJ,bufsiz,MPI_INTEGER,Sendbuf(I),numpes &
                  ,MPI_COMM_WORLD,vrequest2(I),ier)
    IF(ier/=MPI_SUCCESS) CALL MPERROR('Error in (A4) isend',ier)
    CALL MPI_TEST(vrequest2(I),lflag,vstatus2(1,I),ier)
  ENDDO

  DO I = 1,recvs2
    RecvTest = 1; !Initialise the receive test
    DO WHILE(RecvTest/=0)
      CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus2(1,I),ier)
      IF(ier /= MPI_SUCCESS) CALL MPERROR('Error in (A5) probe',ier)
      CALL MPI_GET_COUNT(vstatus2(1,I),MPI_INTEGER,bufsiz,ier)
      IF(ier /= MPI_SUCCESS) CALL MPERROR('Error in (A5) get_count',ier)
      procID = vstatus2(MPI_tag,I)
      RecvTest = 0;
      DO J = 1,I-1
        IF(Recvbuf(J)==procID) RecvTest=RecvTest+1;
      ENDDO
    ENDDO
    L = bufsiz/n_nel;
    Recvbuf(I) = procID;
    CALL MPI_RECV(NODE_ADJ_tmp(I,:,1:L),bufsiz,MPI_INTEGER,procID-1,MPI_ANY_TAG &
                 ,MPI_COMM_WORLD,vstatus2(1,I),ier)
  ENDDO
  CALL MPI_WAITALL(sends,vrequest2,vstatus2,ier)
  IF (ier /= MPI_SUCCESS) CALL MPERROR('Error in MPI_WAITALL', ier)


  !---
  !Form element-element Adjacency array
  !---
  DO Iel = 1,nel_pp
    DO I = 1,nod
      M = gg_pp(I,Iel)
      procID = FIND_NODE_PROC(M,nn_pp1,nn_pp2,numn,npess)
      PROC_CHECK:DO J = 1,recvs
        N = J;
        IF(Recvbuf(J) == procID) EXIT PROC_CHECK;
      ENDDO PROC_CHECK
      K = (I-1)*n_nel + 1;
      L = I*n_nel;
      nodeID = GLOBAL_TO_LOCAL(M,nn_pp1,nn_pp2,numn,npess)
      ELM_ADJ(K:L,Iel) = NODE_ADJ_tmp(N,:,nodeID)
    ENDDO
  ENDDO
  DEALLOCATE(Sendbuf, Recvbuf)
  RETURN
ENDSUBROUTINE ELEMENT_ELEMENT_ADJACENCYGRAPH


!---
! Tenatively colour the graph initially
!---
SUBROUTINE TENATIVE_COLOURING(gg_colour, ELM_ADJ, n_nel, nod &
                            , ncol0, nel_pp, nels, numpes, npess)
  IMPLICIT NONE
  INTEGER                 :: Iel, Iel1, I, J, K, L;
  INTEGER,   INTENT(IN)   :: n_nel, nod, nel_pp, nels, numpes, npess, ncol0
  INTEGER,   INTENT(IN)   :: ELM_ADJ(nod*n_nel,nel_pp);
  INTEGER,   INTENT(INOUT):: gg_colour(nel_pp)
  INTEGER                 :: procID,ElmID,L_ElmID
  INTEGER                 :: nel_pp1, nel_pp2, nume
  INTEGER                 :: Pallette(ncol0)

  CALL PF_PARTITIONER(nel_pp1, nel_pp2, K, nume, nels, npess, numpes) !Elements
  gg_colour = 0;
  DO Iel = 1,nel_pp
    !Delete colours on pallette according to local-neighbors
    Pallette = 0;
    DO I = 1,(nod*n_nel)
      ElmID = ELM_ADJ(I,Iel)
      IF(ElmID/=0)THEN
        procID = FIND_NODE_PROC(ElmID,nel_pp1,nel_pp2,nume,npess)
        IF(procID==numpes)THEN
          L_ElmID = GLOBAL_TO_LOCAL(ElmID,nel_pp1,nel_pp2,nume,npess)
          IF(gg_colour(L_ElmID)/=0)THEN
            Pallette(gg_colour(L_ElmID)) = 1; !Colour occupied
          ENDIF
        ENDIF
      ENDIF
    ENDDO


    !Colour according to minimal local-order
    COLOURING:DO I = 1,ncol0
      IF(Pallette(I)==0)THEN
        gg_colour(Iel) = I;
        EXIT COLOURING;
      ENDIF
    ENDDO COLOURING
  ENDDO
  RETURN
END SUBROUTINE TENATIVE_COLOURING


!---
! colour the graph boundary nodes and recompare
!---
SUBROUTINE BOUNDARY_COLOURING(gg_colour, ncolours, ELM_ADJ, n_nel, nod &
                            , ncol0, nel_pp, nels, numpes, npess)
  IMPLICIT NONE
  INTEGER                 :: Iproc, Iel, Iel1, I, J, K, L, Iters;
  INTEGER,   INTENT(IN)   :: n_nel, nod, nel_pp, nels;
  INTEGER,   INTENT(IN)   :: numpes, npess, ncol0;
  INTEGER,   INTENT(IN)   :: ELM_ADJ(nod*n_nel,nel_pp);
  INTEGER,   INTENT(INOUT):: gg_colour(nel_pp), ncolours;
  INTEGER                 :: procID,ElmID,L_ElmID,lflag,bufsiz;
  INTEGER                 :: nel_pp1, nel_pp2, nume
  INTEGER                 :: Pallette(ncol0);
  INTEGER                 :: vrequest2(npess), vstatus2(MPI_STATUS_SIZE,npess);
  INTEGER,   ALLOCATABLE  :: LprocList(:), GprocList(:), ProcElmAdJ_colours(:,:)
  INTEGER,   ALLOCATABLE  :: Recv_colours(:,:),  Send_colour(:);
  INTEGER                 :: nLprocs, overlap, LprocPos, procPos;
  INTEGER                 :: vstatus(MPI_STATUS_SIZE,npess)
  INTEGER                 :: procTickets(npess), RecvTest;
  LOGICAL                 :: flagger;

  CALL PF_PARTITIONER(nel_pp1, nel_pp2, K, nume, nels, npess, numpes) !Elements

  !=======================
  ! The process connectivity
  !=======================
  !Find all procs that are connected to local process
  ALLOCATE(GprocList(npess)); GprocList=0;
  DO IEL=1,nel_pp
    DO I=1,(nod*n_nel)
      ElmID = ELM_ADJ(I,IEL);
      IF(ElmID/=0)THEN
        procID = FIND_NODE_PROC(ElmID,nel_pp1,nel_pp2,nume,npess)
        GprocList(procID)=1;
      ENDIF
    ENDDO
  ENDDO

  !Calculate the number of local procs
  nLprocs=0;
  DO I=1,npess
    IF(GprocList(I)==1) nLprocs=nLprocs+1;
  ENDDO

  !Form the list of connected procs
  ALLOCATE(LprocList(nLprocs)); LprocList=-1;
  K=0;
  DO I=1,npess
    IF(GprocList(I)==1)THEN
      K=K+1;
      LprocList(K)=I;
    ENDIF
  ENDDO
  DEALLOCATE(GprocList);

  !Find the local proc number in the list
  PROC_POS_FIND:DO I=1,nLprocs
    LprocPos = I;
    IF(LprocList(I)==numpes) EXIT PROC_POS_FIND;
  ENDDO PROC_POS_FIND

  !=======================
  ! Locate boundary elements
  !=======================
  ALLOCATE(ProcElmAdJ_colours(nLprocs,nel_pp))
  ProcElmAdJ_colours=0;
  DO IEL=1,nel_pp
    DO I = 1,(nod*n_nel)
      ElmID = ELM_ADJ(I,IEL);
      IF(ElmID/=0)THEN
        procID = FIND_NODE_PROC(ElmID,nel_pp1,nel_pp2,nume,npess)
        PROCS:DO J=1,nLprocs
          K = J;
          IF(LprocList(J)==procID) EXIT PROCS;
        ENDDO PROCS
        ProcElmAdJ_colours(K,IEL) = 1;
      ENDIF
    ENDDO
  ENDDO

  !LprocList(:), GprocList(:);
  !Go through all colours of the process

  !=======================
  ! Iterative colouring of boundaries
  !=======================
  ALLOCATE(Recv_colours(nLprocs,nel_pp1),  Send_colour(nel_pp1))
  Recv_colours = 0;
  Send_colour  = 0;
  COLOUR_ITERS:DO Iters=1,npess
    !=======================
    ! Colour local colours assuming proc Hierachy
    !=======================
    !colour the boundary elements
    DO IEL=1,nel_pp
      K=0;
      !Checks if element is on Boundary
      DO J = 1,nLprocs
        IF(ProcElmAdJ_colours(J,IEL)/=0) K=K+1;
      ENDDO

      !Only colour if is on processor boundary
      IF(K /= 0)THEN
        !Eliminate lower hierachy pallette colours
        Pallette = 0;
        DO I = 1,(nod*n_nel)
          IF(ELM_ADJ(I,IEL) /= 0)THEN
            ElmID = ELM_ADJ(I,IEL);
            procID = FIND_NODE_PROC(ElmID,nel_pp1,nel_pp2,nume,npess)
            DO J=1,nLprocs
              IF(LprocList(J)==procID) L=J;
            ENDDO
            L_ElmID = GLOBAL_TO_LOCAL(ElmID,nel_pp1,nel_pp2,nume,npess)

            IF( .NOT.((L_ElmID==IEL).AND.(procID==numpes)) )THEN
              IF(Recv_colours(L,L_ElmID)/=0)THEN
                IF(procID<=numpes) Pallette(Recv_colours(L,L_ElmID))=1;
              ENDIF;
            ENDIF

          ENDIF;
        ENDDO

        !Colour Local pallette colour
        COLOUR_SEARCH:DO I=1,ncol0
          IF(Pallette(I)==0)THEN
            Recv_colours(LprocPos,IEL) = I;
            EXIT COLOUR_SEARCH;
          ENDIF
        ENDDO COLOUR_SEARCH
      ENDIF
    ENDDO

    !=======================
    ! Send and recieve non-local colours
    !=======================
    !Barrier to make sure Sends are ready before they start
    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    J=0; bufsiz = nel_pp1;
    Send_colour = Recv_colours(LprocPos,:)
    DO I = 1,nLprocs
      procID = LprocList(I);
      IF(procID /=  numpes)THEN
        J=J+1;
        CALL MPI_ISEND(Send_colour,bufsiz,MPI_INTEGER,(procID-1),numpes &
                      ,MPI_COMM_WORLD,vrequest2(J),ier)
        IF(ier/=MPI_SUCCESS) CALL MPERROR('Error in (A4) isend',ier)
        CALL MPI_TEST(vrequest2(J),lflag,vstatus2(1,J),ier)
      ENDIF
    ENDDO

    bufsiz = nel_pp1;
    Recv_colours = 0;
    procTickets  = -1;
    DO I = 1,(nLprocs-1); flagger = .TRUE.;
      DO WHILE(flagger)
        CALL MPI_PROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,vstatus2(1,I),ier)
        IF(ier /= MPI_SUCCESS) CALL MPERROR('Error in (A5) probe',ier)
        procID = vstatus2(MPI_tag,I)
        RecvTest = 0;
        IF(I /= 1)THEN; DO J = 1,I-1
          IF(procTickets(J) == procID) RecvTest = RecvTest + 1;
        ENDDO; ENDIF
        IF(RecvTest == 0) flagger = .FALSE.;
      ENDDO
      procTickets(I) = procID
      DO J=1,nLprocs
        IF(LprocList(J) == procID) procPos=J;
      ENDDO
      CALL MPI_RECV(Recv_colours(procPos,1:nel_pp1),bufsiz,MPI_INTEGER &
                    , vstatus2(MPI_SOURCE,i),vstatus2(MPI_TAG,i)     &
					,MPI_COMM_WORLD,vstatus2(1,I),ier)
      IF(ier/=MPI_SUCCESS) CALL MPERROR('Error in Recv',ier)
    ENDDO

    CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
    CALL MPI_WAITALL((nLprocs-1),vrequest2,vstatus2,ier)
    IF (ier /= MPI_SUCCESS) CALL MPERROR ('Error in MPI_WAITALL', ier)

    !=======================
    ! Check for convergence
    !=======================
    !check for colour overlap
    overlap=0;
    Recv_colours(LprocPos,:) = Send_colour;
    DO IEL=1,nel_pp
      DO I = 1,(nod*n_nel)
        ElmID = ELM_ADJ(I,Iel)
        IF(ElmID/=0)THEN

          procID  = FIND_NODE_PROC(ElmID,nel_pp1,nel_pp2,nume,npess)
          L_ElmID = GLOBAL_TO_LOCAL(ElmID,nel_pp1,nel_pp2,nume,npess)
          L=-1;
          DO J=1,nLprocs !Boundary node check
            IF(LprocList(J)==procID) L=J;
          ENDDO
          IF(Send_colour(IEL) /= 0)THEN;
            IF( .NOT.((L_ElmID==IEL).AND.(procID==numpes)) )THEN
              IF(Send_colour(IEL)==Recv_colours(L,L_ElmID)) overlap=overlap+1;
            ENDIF
          ENDIF;


        ENDIF
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(overlap,overlap,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier);

    IF(numpes==1) WRITE(*,*) Iters, overlap;
    IF(overlap==0) EXIT COLOUR_ITERS;
  ENDDO COLOUR_ITERS
  IF(numpes==1) WRITE(*,*) "BOUNDARY-ELEMENT COLOURING COMPLETED"

  !=======================
  ! Colour local
  !=======================
  gg_colour(:) = 0;
  gg_colour(:) = Send_colour(1:nel_pp)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier)
  DO Iel = 1,nel_pp
    IF(gg_colour(Iel)==0)THEN
      !Delete colours on pallette according to local-neighbors
      Pallette = 0;
      DO I = 1,(nod*n_nel)
        ElmID = ELM_ADJ(I,Iel)
        IF(ElmID/=0)THEN
          procID = FIND_NODE_PROC(ElmID,nel_pp1,nel_pp2,nume,npess)
          IF(procID==numpes)THEN
            L_ElmID = GLOBAL_TO_LOCAL(ElmID,nel_pp1,nel_pp2,nume,npess)
            IF(gg_colour(L_ElmID)/=0)THEN
              Pallette(gg_colour(L_ElmID)) = 1; !Colour occupied
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      !Colour according to minimal local-order
      COLOURING:DO I = 1,ncol0
        IF(Pallette(I)==0)THEN
          gg_colour(Iel) = I;
          EXIT COLOURING;
        ENDIF
      ENDDO COLOURING
    ENDIF
  ENDDO
  IF(numpes==1) WRITE(*,*) "ELEMENT-GRAPH GRAPH COLOURING COMPLETED"
  ncolours = MAXVAL(gg_colour)
  CALL MPI_ALLREDUCE(ncolours,ncolours,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE BOUNDARY_COLOURING

!---------
! Colour the finite-element (Locally Optimal algorithm)
!--------
SUBROUTINE ELEMENT_COLOURING(gg_colour,gg_pp,nod,nn_pp,nel_pp,nn,nels,npess,numpes)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L;
  INTEGER                 :: n_nel, nproc_pp, ncolours;
  INTEGER,   INTENT(IN)   :: nod, nn_pp, nel_pp, nn, nels, npess, numpes
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel_pp)
  INTEGER,   INTENT(INOUT):: gg_colour(nel_pp)
  INTEGER,   PARAMETER    :: ncol0 = 128;   !Arbitrary Large Pallette Size
  INTEGER,   ALLOCATABLE  :: NODE_ADJ(:,:), ELM_ADJ(:,:)
  INTEGER,   ALLOCATABLE  :: ProcPalletteSize(:), MaximalProcPalletteSize(:)
  INTEGER,   ALLOCATABLE  :: proc_colour(:), Total_Pallette(:), local_pallette(:);

  !---
  ! Calculate the maximal number
  ! of elements connected to a node
  !---
  CALL MAXIMAL_NODAL_CONNECTIVITY(n_nel,nod,nn_pp,nel_pp,npess,numpes)
  IF(numpes==1) WRITE(*,*) "MAXMIMAL NODAL CONNECTIVITY :", n_nel 

  !---
  ! Calculate the Nodal-element
  ! Adjacency-graph/array
  !---
  ALLOCATE(NODE_ADJ(n_nel,nn_pp))
  NODE_ADJ = 0;
  CALL NODE_ELEMENT_ADJACENCYGRAPH(NODE_ADJ,gg_pp,n_nel,nod,nel_pp,nn_pp &
                                  ,nels,nn,numpes,npess)
  IF(numpes==1) WRITE(*,*) "NODE-ELEMENT-ADJACENCY GRAPH FORMED"

  !---
  ! Calculate the element-element
  ! Adjacency-graph/array
  !---
  ALLOCATE(ELM_ADJ(nod*n_nel,nel_pp))
  ELM_ADJ = 0;
  CALL ELEMENT_ELEMENT_ADJACENCYGRAPH(ELM_ADJ,NODE_ADJ,gg_pp,n_nel,nod,nel_pp &
                                     ,nn_pp, nels,nn,numpes,npess)
  DEALLOCATE(NODE_ADJ)
  IF(numpes==1) WRITE(*,*) "ELEMENT-ELEMENT-ADJACENCY GRAPH FORMED"

  !---
  ! Calculate the element-element
  ! Adjacency-graph/array
  !---
  CALL BOUNDARY_COLOURING(gg_colour, ncolours, ELM_ADJ, n_nel, nod &
                        , ncol0, nel_pp, nels, numpes, npess)
  RETURN
ENDSUBROUTINE ELEMENT_COLOURING
!--------------------------------------
!--------------------------------------
!--------------------------------------


! module load gcc7/7.2.0
! module load openmpi-gcc7/4.0.4-cuda11.2


ENDMODULE GRAPH_COLOURING
