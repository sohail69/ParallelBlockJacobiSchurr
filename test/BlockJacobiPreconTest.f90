!
! Patch test for the
! block-Jacobi precon on
! finite elements
!


! mpif90 -o main  BlockJacobiPreconTest.f90 BlockJacobiPrecon.f90 -fallow-argument-mismatch
!
!

PROGRAM MAIN
  USE BlockJacobiPrecon;
  IMPLICIT NONE
  INTEGER, PARAMETER   :: nod=4, nel=5, nn=10, ndim=2, nodof=2;
  INTEGER              :: gg_global(nod,nel);
  INTEGER              :: nn_pp, nn_pp1, nn_pp2, num_nn
  INTEGER              :: nel_pp, nel_pp1, nel_pp2, num_nel
  INTEGER              :: MPIErr, nprocs, procID, tnstart, tnend;
  INTEGER              :: IEL, I, J, K, L, P, Q, S, T;
  INTEGER              :: ntot, neq_pp, neq;
  INTEGER,  ALLOCATABLE:: gg_pp(:,:);
  REAL(iwp),ALLOCATABLE:: stork_pp(:,:,:), Amat(:,:);
  character(len=1024)  :: procID_string
  INTEGER              :: nSends, nRecvs, totSends, totRecvs, errMPI

  INTEGER,  ALLOCATABLE:: sendTables(:,:), recvTables(:,:);
  INTEGER,  ALLOCATABLE:: sendNodes(:,:), recvNodes(:,:)

!===============================
!
!   Initialise the system and
!   use a sample mesh and
!   partition the problem
!
!===============================
  CALL MPI_INIT(MPIErr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, MPIErr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, procID, MPIErr)

  !The element topology/steering array
  !for patch test mesh
  gg_global(:,1) = (/ 1, 2, 5, 4/);
  gg_global(:,2) = (/ 2, 3, 6, 5/);
  gg_global(:,3) = (/ 4, 5, 6, 8/);
  gg_global(:,4) = (/ 4, 8, 9, 7/);
  gg_global(:,5) = (/ 8, 6,10, 9/);


  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num_nn, nn, nprocs, procID+1)
  CALL PF_PARTITIONER1(nel_pp1, nel_pp2, nel_pp, num_nel, nel, nprocs, procID+1)

  ntot = nod*nodof;
  neq_pp = nn_pp*nodof;
  neq = nn*nodof;


!===============================
!
!  Assemble the local system
!  (using a global size matrix
!  its a small problem)
!
!===============================
  ALLOCATE(Amat(neq,neq))
  IF(nel_pp /= 0) ALLOCATE(gg_pp(nod,nel_pp))
  IF(nel_pp /= 0) ALLOCATE(stork_pp(ntot,ntot,nel_pp))
  stork_pp = 1._iwp;

  tnstart = ITERATOR_START(procID, nprocs, nel)
  tnend   = ITERATOR_END(procID, nprocs, nel) 
  DO IEL = 1,nel_pp
    gg_pp(:,IEL) = gg_global(:,tnstart+IEL-1)
  ENDDO

!===============================
!
! Constructing the message tables
! and making certain there is
! message parity
!
!===============================
  CALL SIZE_MESSAGE_TABLES(nSends, nRecvs, gg_pp, procID &
                         , nProcs, nn, nel, nod, nel_pp)

  ALLOCATE(sendTables(2,nSends), recvTables(2,nRecvs))
  ALLOCATE(sendNodes(nSends,nn_pp1), recvNodes(nRecvs,nn_pp1))


totSends=0;
totRecvs=0;
CALL MPI_ALLREDUCE(nRecvs,totRecvs,1,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,errMPI);
CALL MPI_ALLREDUCE(nSends,totSends,1,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,errMPI);
				   
WRITE(*,*) "TOTALS", totSends, totRecvs
IF(procID==0) WRITE(*,*) "/*-----------------------------------------*\"

  CALL FORM_MESSAGE_TABLES(sendTables, recvTables, sendNodes &
                         , recvNodes, gg_pp, nRecvs, nSends  &
                         , nn , nod, nel_pp, nn_pp1, nprocs  &
                         , procID)


!===============================
!
!  Output the data
!
!===============================
  WRITE(procID_string, "(I2)") (procID+1) 
  OPEN(UNIT=(procID+1),FILE = "name"//trim(procID_string)//".txt")
  DO I = 1,neq
    WRITE( (procID+1),*) Amat(I,:)
  ENDDO

  DO IEL = 1,nel_pp
    WRITE( (procID+1),*) gg_pp(:,IEL)
  ENDDO
  CLOSE(procID+1)


  IF(nel_pp /= 0) DEALLOCATE(gg_pp, stork_pp)
  DEALLOCATE(Amat)
  DEALLOCATE(sendTables, recvTables)
  DEALLOCATE(sendNodes, recvNodes)

  CALL MPI_FINALIZE(MPIErr)
ENDPROGRAM MAIN