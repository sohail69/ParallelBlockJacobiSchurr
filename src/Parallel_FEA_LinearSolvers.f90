MODULE Parallel_FEA_LinearSolvers
  USE precision;      USE global_variables;
  USE maths;          USE gather_scatter;
  USE new_library;    USE MP_INTERFACE;
!
! My own libraries
!
  USE Parallel_supplementary_Maths;
  CONTAINS

SUBROUTINE STOP_COND(IsConverged, error, rtol)
  USE precision;       USE global_variables;
  USE gather_scatter;  USE MP_INTERFACE;
  USE maths;
  IMPLICIT NONE
  INTEGER                 :: IsConverged_tmp
  INTEGER,   INTENT(INOUT):: IsConverged
  REAL(iwp), INTENT(IN)   :: error, rtol

  IsConverged = 0;
  IF(error < rtol) IsConverged = 1;
  CALL MPI_ALLREDUCE(IsConverged,IsConverged_tmp,1,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,ier);
  IsConverged = IsConverged_tmp;
  RETURN
END SUBROUTINE STOP_COND

!-------------------------------
! Conjugate gradient method
!-------------------------------
  SUBROUTINE SolveLinearSystem_CG(A_mat, x_vec, b_vec, NodalMask, gg_colour    &
                                , ntots, nod, ncolours, nodof, nel_pp, neqs_pp &
                                , nn_pp, ltol, limit, iters, error, precon)
    IMPLICIT NONE
    INTEGER                 :: Iel, i, j, IsConverged;
    INTEGER,   INTENT(INOUT):: iters;
    INTEGER,   INTENT(IN)   :: nodof, nod, nn_pp, precon, ncolours;
    INTEGER,   INTENT(IN)   :: NodalMask(nodof,nod), gg_colour(nel_pp);
    INTEGER,   INTENT(IN)   :: ntots, nel_pp, neqs_pp, limit;
    REAL(iwp), INTENT(IN)   :: A_mat(ntots,ntots,nel_pp), b_vec(neqs_pp), ltol;
    REAL(iwp), INTENT(INOUT):: x_vec(neqs_pp), error;
    REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;
    REAL(iwp), ALLOCATABLE  :: r_pp(:), p_pp(:), ap_pp(:), d_pp(:)
    REAL(iwp), ALLOCATABLE  :: pmul(:,:), qmul(:,:)
    REAL(iwp)               :: alpha, beta, rho0, b_norm, rpap_norm;

    !-----
    ! Allocate and initialise arrays
    !-----
    ALLOCATE(pmul(ntots,nel_pp), qmul(ntots,nel_pp));
    ALLOCATE(r_pp(neqs_pp), p_pp(neqs_pp), ap_pp(neqs_pp),d_pp(neqs_pp))
    pmul  = zero;
    qmul  = zero;
    r_pp  = zero;
    p_pp  = zero;
    ap_pp = zero;
    x_vec = zero;
    d_pp  = zero;


    !-----
    ! Initialise linear system
    !-----
    CALL PARAMATVEC(A_mat,x_vec,r_pp,pmul,qmul,NodalMask &
                    ,nel_pp,nn_pp,neqs_pp,ntots,nod,nodof)
    !PARAMATVEC(A_mat, x_vec, r_pp)
    b_norm = DOT_PRODUCT_P(b_vec,b_vec)
    b_norm = DSQRT(b_norm)
    r_pp  = b_vec - r_pp;
    error = norm_p(r_pp)
    error = error/b_norm
    CALL STOP_COND(IsConverged, error, ltol)
    IF(IsConverged /= 0)RETURN


   !!!d_pp = M_Inv * r_pp;
    d_pp = r_pp
    p_pp  = d_pp;

IF(numpe==1) WRITE(*,*) "CG Solver started"

    !-----
    ! Solve Linear system using CG
    !-----
    iters = 0;
    ITERATIONS:DO i = 1,limit
      iters = iters + 1;
      CALL PARAMATVEC(A_mat,p_pp,ap_pp,pmul,qmul,NodalMask &
                      ,nel_pp,nn_pp,neqs_pp,ntots,nod,nodof)
      !PARAMATVEC(A_mat, p_pp, ap_pp)
      rho0  = DOT_PRODUCT_P(d_pp,r_pp)
      rpap_norm = DOT_PRODUCT_P(p_pp,ap_pp)
      alpha = rho0/rpap_norm;
      x_vec = x_vec + alpha*p_pp
      r_pp  = r_pp - alpha*ap_pp
      error = norm_p(r_pp)
      error = error/b_norm
      !!!d_pp = M_Inv * r_pp;
      d_pp = r_pp
      CALL STOP_COND(IsConverged, error, ltol)
      IF(IsConverged /= 0) EXIT ITERATIONS
      beta  = DOT_PRODUCT_P(r_pp,d_pp)
      beta  = beta/rho0
      beta  = MAX(zero,beta)
      p_pp  = d_pp + beta*p_pp     
    ENDDO ITERATIONS
    
    !-----
    ! Deallocate arrays and exit
    !-----
    DEALLOCATE(pmul, qmul);
    DEALLOCATE(r_pp, p_pp, ap_pp, d_pp);
    RETURN
  END SUBROUTINE SolveLinearSystem_CG

!-------------------------------
! Generalised minimum residual method (GMRES) and Direct-Quasi GMRES
!-------------------------------
  SUBROUTINE SolveLinearSystem_GMRESR(A_mat, M_mat, x_vec, b_vec, NodalMask, gg_colour  &
                                   ,ncolours, ntots, nod, nodof, nel_pp, neqs_pp &
                                   ,nn_pp,ltol,limit,iters,ell,error,precon,option)
    IMPLICIT NONE
    INTEGER                 :: Iel, i, j, k, riters, m, n, IsConverged;
    REAL(iwp)               :: BiCGerr;
    INTEGER,   INTENT(INOUT):: iters
    INTEGER,   INTENT(IN)   :: nodof, nod, nn_pp, ncolours, option, precon;
    INTEGER,   INTENT(IN)   :: NodalMask(nodof,nod), gg_colour(nel_pp);
    INTEGER,   INTENT(IN)   :: ntots, nel_pp, neqs_pp, limit, ell;
    REAL(iwp), INTENT(IN)   :: A_mat(ntots,ntots,nel_pp), b_vec(neqs_pp), ltol;
    REAL(iwp), INTENT(IN)   :: M_mat(ntots,ntots,nel_pp);
    REAL(iwp), INTENT(INOUT):: x_vec(neqs_pp), error;
    REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;
    REAL(iwp), ALLOCATABLE  :: r0_pp(:), w_pp(:), v_pp(:,:), z_pp(:,:);
    REAL(iwp), ALLOCATABLE  :: Hess(:,:), HessInv(:,:), sn(:), cs(:)
    REAL(iwp), ALLOCATABLE  :: beta(:), ytemp(:);
    REAL(iwp), ALLOCATABLE  :: pmul(:,:), qmul(:,:), x_pp(:);
    REAL(iwp)               :: temp1, temp2, cs_k, sn_k, b_norm, errnorm0, errnorm1;

character(len=10) :: file_id

    !-----
    ! Allocate and initialise Solver arrays
    !-----
    ALLOCATE(pmul(ntots,nel_pp), qmul(ntots,nel_pp), x_pp(neqs_pp));
    ALLOCATE(r0_pp(neqs_pp), Hess(ell+1,ell+1), sn(ell+1), cs(ell+1))
    ALLOCATE(w_pp(neqs_pp), v_pp(neqs_pp,ell+1), z_pp(neqs_pp,ell+1), beta(ell+1))
    Hess  = zero;    sn    = zero;
    cs    = zero;    beta  = zero;
    pmul  = zero;    qmul  = zero;
    x_vec = zero;    r0_pp = zero;
    v_pp  = zero;    w_pp  = zero;
    x_pp  = zero;
    b_norm = norm_p(b_vec)

write(file_id, '(i0)') precon
 IF(numpe==1)OPEN(43,FILE="Results2/LinearSolver"//trim(adjustl(file_id))//".dat",STATUS='REPLACE',ACTION='WRITE')

    iters = 0;
    RESTART_ITERATIONS:DO riters = 1,limit
      !-----
      ! Initialise linear system
      !-----
      CALL PARAMATVEC(A_mat,x_pp,r0_pp,pmul,qmul,NodalMask &
                     ,nel_pp,nn_pp,neqs_pp,ntots,nod,nodof)
      !PARAMATVEC(A_mat, x_vec, r0_pp)
      r0_pp = b_vec - r0_pp;
      beta(1)   = norm_p(r0_pp);
      v_pp(:,1) = r0_pp/beta(1)
      errnorm0  = beta(1)

      !-----
      ! Solve Linear system using GMRES(l)
      !-----
      GMRES_ITERATIONS:DO i = 1,ell
        iters = iters +1;
        !-----
        ! Flexible Preconditioned Arnoldi iterations
        !-----
!==================================================================
!---------------------------Preconditioner-------------------------
!==================================================================
           z_pp(:,i) = v_pp(:,i);

           CALL PARAMATVEC(A_mat,z_pp(:,i),w_pp,pmul,qmul,NodalMask &
                          ,nel_pp,nn_pp,neqs_pp,ntots,nod,nodof)
!==================================================================
!---------------------------Preconditioner-------------------------
!==================================================================
        m = i - riters + 1;
        IF(option == 0) n = 1          !Fully orthogonalized GMRESR
        IF(option == 1) n = MAX(1,m)   !Direct Incomplete-Orthogonalized GMRESR
        ARNOLDI_ITERATION:DO j = n,i
          Hess(j,i) = DOT_PRODUCT_P(v_pp(:,j), w_pp)
          w_pp = w_pp - Hess(j,i)*v_pp(:,j)
        ENDDO ARNOLDI_ITERATION
        Hess(i+1,i) = norm_p(w_pp);
        v_pp(:,i+1) = w_pp/Hess(i+1,i);

        !-----
        ! Apply Givens Rotations
        !-----
        GIVENS_ROTATIONS:DO j = MAX(1,n-1),(i-1)
          temp1       =  cs(j)*Hess(j,i) + sn(j)*Hess(j+1,i)
          Hess(j+1,i) = -sn(j)*Hess(j,i) + cs(j)*Hess(j+1,i)
          Hess(j,i)   =  temp1;
        ENDDO GIVENS_ROTATIONS
        temp2 = DSQRT(Hess(i,i)*Hess(i,i) + Hess(i+1,i)*Hess(i+1,i))
        cs_k = Hess(i,i)/temp2
        sn_k = Hess(i+1,i)/temp2
        Hess(i,i) = cs_k*Hess(i,i) + sn_k*Hess(i+1,i);
        Hess(i+1,i) = zero;
        cs(i) = cs_k;
        sn(i) = sn_k;

        !-----
        ! Update residual vector
        !-----
        beta(i+1) = -sn_k*beta(i)
        beta(i)   =  cs_k*beta(i)
        error     =  DABS(beta(i+1))/b_norm;
        k = i;
 IF(numpe==1)WRITE(43,*) iters, error

        CALL STOP_COND(IsConverged, error, ltol)

        IF(IsConverged /= 0) EXIT GMRES_ITERATIONS;
      ENDDO GMRES_ITERATIONS

      !-----
      ! Construct solution
      !-----
      ALLOCATE(HessInv(k,k), ytemp(k))
      HessInv = zero;
      ytemp   = zero;
      CALL Invert2(Hess(1:k,1:k), HessInv, k)
      ytemp = MATMUL(HessInv,beta(1:k))
      IF(option == 0)  Hess = zero;
      x_pp = x_pp + MATMUL(z_pp(:,1:k), ytemp(1:k));

      DEALLOCATE(HessInv, ytemp)
      pmul  = zero;
      qmul  = zero;
      r0_pp = zero;
      v_pp  = zero;
      w_pp  = zero;
      sn    = zero;
      cs    = zero;
      beta  = zero;
      CALL STOP_COND(IsConverged, error, ltol)
      IF(IsConverged /= 0) EXIT RESTART_ITERATIONS;
    ENDDO RESTART_ITERATIONS
    x_vec = x_pp

    !-----
    ! Calculate the error norm
    !-----
    b_norm   = norm_p(b_vec)
    CALL PARAMATVEC(A_mat,x_pp,r0_pp,pmul,qmul,NodalMask &
                   ,nel_pp,nn_pp,neqs_pp,ntots,nod,nodof)
    !PARAMATVEC(A_mat, x_vec, r0_pp)
    r0_pp = b_vec - r0_pp;
    errnorm1 = norm_p(r0_pp)
    error = errnorm1/b_norm
 IF(numpe==1)CLOSE(43)

    !-----
    ! Deallocate arrays and exit
    !-----
    DEALLOCATE(pmul, qmul, x_pp)
    DEALLOCATE(r0_pp, Hess, sn, cs)
    DEALLOCATE(w_pp, v_pp, z_pp, beta)
    RETURN
  END SUBROUTINE SolveLinearSystem_GMRESR
!-------------------------------
!-------------------------------
!-------------------------------
END MODULE Parallel_FEA_LinearSolvers
