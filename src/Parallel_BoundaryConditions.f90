MODULE Parallel_BoundaryConditions
  USE precision;      USE global_variables;
  USE maths;          USE new_library;
  USE MP_INTERFACE;   USE gather_scatter;
  CONTAINS
!-------------------------------
! Apply Dirchlet boundary conditions by matrix elimination
!-------------------------------
SUBROUTINE DIRICHLET_BC(stork_pp, rmul_pp, g_g_pp, val_pp, ntots)
   IMPLICIT NONE
   INTEGER                  :: i;
   INTEGER,   INTENT(IN)    :: ntots
   INTEGER,   INTENT(IN)    :: g_g_pp(ntots)
   REAL(iwp), INTENT(IN)    :: val_pp(ntots)
   REAL(iwp), INTENT(INOUT) :: stork_pp(ntots,ntots),rmul_pp(ntots)
   REAL(iwp), PARAMETER     :: zero = 0._iwp, one = 1._iwp;

   DOFS:DO i = 1,ntots
     IF(g_g_pp(i) == 0)THEN
       stork_pp(i,:) = zero;
       stork_pp(i,i) = one;
       rmul_pp(i)    = stork_pp(i,i)*val_pp(i);
     ENDIF
   ENDDO DOFS
   RETURN
ENDSUBROUTINE DIRICHLET_BC


!-------------------------------
! Apply 3-D No-rotation BC for XY-plane with radius smoothing
!-------------------------------
SUBROUTINE SMOOTHED_NO_ROTATION_BC(Km_mat, Rm_vec, utemp, dutemp, gg_coord, gg_pp &
                                     , centre, rad0, disp_MAP, ntots, ndim, nod, nel_pp)
  IMPLICIT NONE
  INTEGER                  :: Iel, Inode, I, J, K, L;
  INTEGER,   INTENT(IN)    :: ntots, ndim, nod, nel_pp;
  INTEGER,   INTENT(IN)    :: disp_MAP(nod*ndim), gg_pp(nod*ndim,nel_pp);
  REAL(iwp), INTENT(IN)    :: rad0, centre(ndim), dutemp(ntots,nel_pp);
  REAL(iwp), INTENT(IN)    :: gg_coord(nod,ndim,nel_pp), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT) :: Km_mat(ntots,ntots,nel_pp), Rm_vec(ntots,nel_pp);
  REAL(iwp), PARAMETER     :: zero=0._iwp, one=1._iwp, two = 2._iwp;
  REAL(iwp), PARAMETER     :: tol=1.0E-03_iwp;
  REAL(iwp)                :: cx, cy, ux, uy, dux, duy;
  REAL(iwp)                :: CR2_avg, r_accum, N_nodes;
  REAL(iwp)                :: rx, ry, rx0, ry0, rx_rot, ry_rot;
  REAL(iwp)                :: Jxx, Jxy, Jyx, Jyy
  REAL(iwp)                :: Kx_vec(ntots), Ky_vec(ntots);
  REAL(iwp)                :: Kx0_vec(ntots), Ky0_vec(ntots);


  !
  ! Calculate the average radius
  !
  r_accum = zero;
  N_nodes = zero;
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);

        IF(DABS(cx*cx + cy*cy - rad0*rad0) < tol)THEN
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1); 
          J  = disp_MAP(K+2);
          ux = utemp(I,Iel)
          uy = utemp(J,Iel)
          dux = dutemp(I,Iel)
          duy = dutemp(J,Iel)
          r_accum = r_accum + (cx + ux + dux)*(cx + ux + dux) + (cy + uy + duy)*(cy + uy + duy);
          N_nodes = N_nodes + one;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(r_accum,r_accum,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  CALL MPI_ALLREDUCE(N_nodes,N_nodes,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  CR2_avg = r_accum/N_nodes;


  !
  ! Apply the Boundary conditions
  !
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,Iel) == 0)THEN
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);

        IF(DABS(cx*cx + cy*cy - rad0*rad0) < tol)THEN
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1); 
          J  = disp_MAP(K+2);

          dux = dutemp(I,Iel)
          duy = dutemp(J,Iel)

          ux = utemp(I,Iel);
          uy = utemp(J,Iel);

          rx0 = Rm_vec(I,Iel);
          ry0 = Rm_vec(J,Iel);

          Kx0_vec(:) = Km_mat(I,:,Iel)
          Ky0_vec(:) = Km_mat(J,:,Iel)

          rx_rot = cx*cx*CR2_avg - ( (cx+ux+dux)*(cx+ux+dux)*(cx*cx + cy*cy) )
          ry_rot = cy*cy*CR2_avg - ( (cy+uy+duy)*(cy+uy+duy)*(cx*cx + cy*cy) ) 

          rx  = (one + rx0)*rx_rot + rx0;
          ry  = (one + ry0)*ry_rot + ry0;

          Jxx = two*( (cx*cx/N_nodes)- (cx*cx + cy*cy) )*(cx + ux + dux)
          Jxy = two*(cx*cx/N_nodes)*(cy + uy + duy)
          Jyx = two*(cy*cy/N_nodes)*(cx + ux + dux)
          Jyy = two*( (cy*cy/N_nodes)- (cx*cx + cy*cy) )*(cy + uy + duy)

          Kx_vec    = Kx0_vec*(rx_rot + one);
          Kx_vec(I) = Kx_vec(I) + (one + rx0)*Jxx
          Kx_vec(J) = Kx_vec(J) + (one + rx0)*Jxy

          Ky_vec     = Ky0_vec*(ry_rot + one);
          Ky_vec(I) = Ky_vec(I) + (one + ry0)*Jyx
          Ky_vec(J) = Ky_vec(J) + (one + ry0)*Jyy

          Rm_vec(I,Iel)   = rx;
          Rm_vec(J,Iel)   = ry;
          Km_mat(I,:,Iel) = Kx_vec
          Km_mat(J,:,Iel) = Ky_vec
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  RETURN
ENDSUBROUTINE SMOOTHED_NO_ROTATION_BC

!-------------------------------
! Apply 3-D No-rotation BC for XY-plane with radius smoothing
!-------------------------------
SUBROUTINE RSMOOTHED_NO_ROTATION_BC_XY(Km_mat, Rm_vec, utemp, gg_coord, gg_pp &
                                     , centre, rad0, disp_MAP, ntots, ndim, nod, nel_pp)
  IMPLICIT NONE
  INTEGER                  :: Iel, Inode, I, J, K, L;
  INTEGER,   INTENT(IN)    :: ntots, ndim, nod, nel_pp;
  INTEGER,   INTENT(IN)    :: disp_MAP(nod*ndim), gg_pp(nod*ndim,nel_pp);
  REAL(iwp), INTENT(IN)    :: rad0, centre(ndim);
  REAL(iwp), INTENT(IN)    :: gg_coord(nod,ndim,nel_pp), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT) :: Km_mat(ntots,ntots,nel_pp), Rm_vec(ntots,nel_pp);
  REAL(iwp), PARAMETER     :: zero=0._iwp, one=1._iwp, two = 2._iwp;
  REAL(iwp), PARAMETER     :: three = 3._iwp, four = 4._iwp, six = 6._iwp;
  REAL(iwp), PARAMETER     :: tol=1.0E-03_iwp;
  REAL(iwp)                :: cx, cy, ux, uy;
  REAL(iwp)                :: Bx, By, r, r_avg, G, A, N, sqrG, sqrA, sqrN;
  REAL(iwp)                :: rx, ry, Jxx, Jxy, Jyy;
  REAL(iwp)                :: F1, F2, F3, Ftemp1, Ftemp2, Ftemp3;

  !
  ! Calculate the average radius
  !
  N = zero;
  A = zero;
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        r = DSQRT(cx*cx + cy*cy);

        IF(DABS(r - rad0) < tol)THEN
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1); 
          J  = disp_MAP(K+2);

          Bx = gg_coord(Inode,1,Iel) + utemp(I,Iel) - centre(1);
          By = gg_coord(Inode,2,Iel) + utemp(J,Iel) - centre(2);

          A = A + (Bx*Bx + By*By)
          N = N + one;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(A,A,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  CALL MPI_ALLREDUCE(N,N,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  sqrN  = DSQRT(N);
  sqrA  = DSQRT(A);
  r_avg = sqrA/sqrN;

  !
  ! Apply the Boundary conditions
  !
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,Iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        r = DSQRT(cx*cx + cy*cy);

        IF(DABS(r - rad0) < tol)THEN
          !
          ! Map displacements and residuals
          !
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1);
          J  = disp_MAP(K+2);
          ux = utemp(I,Iel);
          uy = utemp(J,Iel);

          !
          ! Build No-XY-rotation BC
          !
          rx =  -four*cy*(cx*uy - cy*ux)*(cx*uy - cy*ux)*(cx*uy - cy*ux);
          ry =   four*cx*(cx*uy - cy*ux)*(cx*uy - cy*ux)*(cx*uy - cy*ux);

          Jxx =  two*six*cy*cy*(cx*uy - cy*ux)*(cx*uy - cy*ux);
          Jxy = -two*six*cy*cx*(cx*uy - cy*ux)*(cx*uy - cy*ux);
          Jyy =  two*six*cx*cx*(cx*uy - cy*ux)*(cx*uy - cy*ux);

          Rm_vec(I,Iel)   = Rm_vec(I,Iel)   + rx;
          Km_mat(I,I,Iel) = Km_mat(I,I,Iel) + Jxx;
          Km_mat(I,J,Iel) = Km_mat(I,J,Iel) + Jxy;

          Rm_vec(J,Iel)   = Rm_vec(J,Iel)   + ry;
          Km_mat(J,I,Iel) = Km_mat(J,I,Iel) + Jxy;
          Km_mat(J,J,Iel) = Km_mat(J,J,Iel) + Jyy;

          !
          ! Build Average-Radius-smoothing BC
          !
          Bx   = gg_coord(Inode,1,Iel) + ux - centre(1);
          By   = gg_coord(Inode,2,Iel) + uy - centre(2);
          G    = Bx*Bx + By*By;
          sqrG = DSQRT(G);

          Ftemp1 =  (r_avg - r);
          Ftemp2 =  (one/(sqrA*sqrN) - one/(sqrG));
          Ftemp3 = -(one/(A*sqrA*sqrN) - one/(G*sqrG));
 
          F1 = two*Ftemp2*Ftemp2;
          F2 = two*Ftemp3*Ftemp1;
          F3 = two*Ftemp1*Ftemp2;

          rx = F1*F2*Bx;
          ry = F1*F2*By;

          Jxx = (F1 + F2)*Bx*Bx + F3;
          Jxy = (F1 + F2)*Bx*By;
          Jyy = (F1 + F2)*By*By + F3;

          Rm_vec(I,Iel)   = Rm_vec(I,Iel)   + rx;
          Km_mat(I,I,Iel) = Km_mat(I,I,Iel) + Jxx;
          Km_mat(I,J,Iel) = Km_mat(I,J,Iel) + Jxy;

          Rm_vec(J,Iel)   = Rm_vec(J,Iel)   + ry;
          Km_mat(J,I,Iel) = Km_mat(J,I,Iel) + Jxy;
          Km_mat(J,J,Iel) = Km_mat(J,J,Iel) + Jyy;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE RSMOOTHED_NO_ROTATION_BC_XY


SUBROUTINE RSMOOTHED_NO_ROTATION_BC_XY_RES(Rm_vec, utemp, gg_coord, gg_pp, centre, rad0 &
                                         , disp_MAP, ntots, ndim, nod, nel_pp)
  IMPLICIT NONE
  INTEGER                  :: Iel, Inode, I, J, K, L;
  INTEGER,   INTENT(IN)    :: ntots, ndim, nod, nel_pp;
  INTEGER,   INTENT(IN)    :: disp_MAP(nod*ndim), gg_pp(nod*ndim,nel_pp);
  REAL(iwp), INTENT(IN)    :: rad0, centre(ndim);
  REAL(iwp), INTENT(IN)    :: gg_coord(nod,ndim,nel_pp), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT) :: Rm_vec(ntots,nel_pp);
  REAL(iwp), PARAMETER     :: zero=0._iwp, one=1._iwp, two = 2._iwp;
  REAL(iwp), PARAMETER     :: three = 3._iwp, four = 4._iwp, six = 6._iwp;
  REAL(iwp), PARAMETER     :: tol=1.0E-03_iwp;
  REAL(iwp)                :: cx, cy, ux, uy;
  REAL(iwp)                :: Bx, By, r, r_avg, G, A, N, sqrG, sqrA, sqrN;
  REAL(iwp)                :: rx, ry, F1, F2, Ftemp1, Ftemp2, Ftemp3;

  !
  ! Calculate the average radius
  !
  N = zero;
  A = zero;
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        r = DSQRT(cx*cx + cy*cy);

        IF(DABS(r - rad0) < tol)THEN
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1); 
          J  = disp_MAP(K+2);

          Bx = gg_coord(Inode,1,Iel) + utemp(I,Iel) - centre(1);
          By = gg_coord(Inode,2,Iel) + utemp(J,Iel) - centre(2);

          A = A + (Bx*Bx + By*By)
          N = N + one;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(A,A,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  CALL MPI_ALLREDUCE(N,N,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  sqrN  = DSQRT(N);
  sqrA  = DSQRT(A);
  r_avg = sqrA/sqrN;

  !
  ! Apply the Boundary conditions
  !
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,Iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        r = DSQRT(cx*cx + cy*cy);

        IF(DABS(r - rad0) < tol)THEN
          !
          ! Map displacements and residuals
          !
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1);
          J  = disp_MAP(K+2);
          ux = utemp(I,Iel);
          uy = utemp(J,Iel);

          !
          ! Build No-XY-rotation BC
          !
          rx =  -four*cy*(cx*uy - cy*ux)*(cx*uy - cy*ux)*(cx*uy - cy*ux);
          ry =   four*cx*(cx*uy - cy*ux)*(cx*uy - cy*ux)*(cx*uy - cy*ux);

          Rm_vec(I,Iel)   = Rm_vec(I,Iel)   + rx;
          Rm_vec(J,Iel)   = Rm_vec(J,Iel)   + ry;

          !
          ! Build Average-Radius-smoothing BC
          !
          Bx   = gg_coord(Inode,1,Iel) + ux - centre(1);
          By   = gg_coord(Inode,2,Iel) + uy - centre(2);
          G    = Bx*Bx + By*By;
          sqrG = DSQRT(G);

          Ftemp1 =  (r_avg - r);
          Ftemp2 =  (one/(sqrA*sqrN) - one/(sqrG));
          Ftemp3 = -(one/(A*sqrA*sqrN) - one/(G*sqrG));

          F1 = two*Ftemp2*Ftemp2;
          F2 = two*Ftemp3*Ftemp1;

          rx = F1*F2*Bx;
          ry = F1*F2*By;

          Rm_vec(I,Iel)   = Rm_vec(I,Iel)   + rx;
          Rm_vec(J,Iel)   = Rm_vec(J,Iel)   + ry;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE RSMOOTHED_NO_ROTATION_BC_XY_RES



SUBROUTINE RSMOOTHED_NO_ROTATION_BC_XY_JAC(Km_mat, utemp, gg_coord, gg_pp &
                                     , centre, rad0, disp_MAP, ntots, ndim, nod, nel_pp)
  IMPLICIT NONE
  INTEGER                  :: Iel, Inode, I, J, K, L;
  INTEGER,   INTENT(IN)    :: ntots, ndim, nod, nel_pp;
  INTEGER,   INTENT(IN)    :: disp_MAP(nod*ndim), gg_pp(nod*ndim,nel_pp);
  REAL(iwp), INTENT(IN)    :: rad0, centre(ndim);
  REAL(iwp), INTENT(IN)    :: gg_coord(nod,ndim,nel_pp), utemp(ntots,nel_pp);
  REAL(iwp), INTENT(INOUT) :: Km_mat(ntots,ntots,nel_pp);
  REAL(iwp), PARAMETER     :: zero=0._iwp, one=1._iwp, two = 2._iwp;
  REAL(iwp), PARAMETER     :: three = 3._iwp, four = 4._iwp, six = 6._iwp;
  REAL(iwp), PARAMETER     :: tol=1.0E-03_iwp;
  REAL(iwp)                :: cx, cy, ux, uy;
  REAL(iwp)                :: Bx, By, r, r_avg, G, A, N, sqrG, sqrA, sqrN;
  REAL(iwp)                :: rx, ry, Jxx, Jxy, Jyy;
  REAL(iwp)                :: F1, F2, F3, Ftemp1, Ftemp2, Ftemp3;

  !
  ! Calculate the average radius
  !
  N = zero;
  A = zero;
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        r = DSQRT(cx*cx + cy*cy);

        IF(DABS(r - rad0) < tol)THEN
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1); 
          J  = disp_MAP(K+2);

          Bx = gg_coord(Inode,1,Iel) + utemp(I,Iel) - centre(1);
          By = gg_coord(Inode,2,Iel) + utemp(J,Iel) - centre(2);

          A = A + (Bx*Bx + By*By)
          N = N + one;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  CALL MPI_ALLREDUCE(A,A,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  CALL MPI_ALLREDUCE(N,N,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier);
  sqrN  = DSQRT(N);
  sqrA  = DSQRT(A);
  r_avg = sqrA/sqrN;

  !
  ! Apply the Boundary conditions
  !
  DO Iel = 1,nel_pp
    DO Inode = 1,nod
      L = Inode*ndim;
      IF(gg_pp(L,Iel) == 0)THEN
        !
        ! Find position relative to centre of rotation
        !
        cx = gg_coord(Inode,1,Iel) - centre(1);
        cy = gg_coord(Inode,2,Iel) - centre(2);
        r = DSQRT(cx*cx + cy*cy);

        IF(DABS(r - rad0) < tol)THEN
          !
          ! Map displacements and residuals
          !
          K  = ndim*(Inode-1);
          I  = disp_MAP(K+1);
          J  = disp_MAP(K+2);
          ux = utemp(I,Iel);
          uy = utemp(J,Iel);

          !
          ! Build No-XY-rotation BC
          !
          Jxx =  two*six*cy*cy*(cx*uy - cy*ux)*(cx*uy - cy*ux);
          Jxy = -two*six*cy*cx*(cx*uy - cy*ux)*(cx*uy - cy*ux);
          Jyy =  two*six*cx*cx*(cx*uy - cy*ux)*(cx*uy - cy*ux);

          Km_mat(I,I,Iel) = Km_mat(I,I,Iel) + Jxx;
          Km_mat(I,J,Iel) = Km_mat(I,J,Iel) + Jxy;
          Km_mat(J,I,Iel) = Km_mat(J,I,Iel) + Jxy;
          Km_mat(J,J,Iel) = Km_mat(J,J,Iel) + Jyy;

          !
          ! Build Average-Radius-smoothing BC
          !
          Bx   = gg_coord(Inode,1,Iel) + ux - centre(1);
          By   = gg_coord(Inode,2,Iel) + uy - centre(2);
          G    = Bx*Bx + By*By;
          sqrG = DSQRT(G);

          Ftemp1 =  (r_avg - r);
          Ftemp2 =  (one/(sqrA*sqrN) - one/(sqrG));
          Ftemp3 = -(one/(A*sqrA*sqrN) - one/(G*sqrG));
 
          F1 = two*Ftemp2*Ftemp2;
          F2 = two*Ftemp3*Ftemp1;
          F3 = two*Ftemp1*Ftemp2;

          Jxx = (F1 + F2)*Bx*Bx + F3;
          Jxy = (F1 + F2)*Bx*By;
          Jyy = (F1 + F2)*By*By + F3;

          Km_mat(I,I,Iel) = Km_mat(I,I,Iel) + Jxx;
          Km_mat(I,J,Iel) = Km_mat(I,J,Iel) + Jxy;
          Km_mat(J,I,Iel) = Km_mat(J,I,Iel) + Jxy;
          Km_mat(J,J,Iel) = Km_mat(J,J,Iel) + Jyy;
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE RSMOOTHED_NO_ROTATION_BC_XY_JAC


!-------------------------------
! Detects local boundary Face elements and returns Face steering array
!-------------------------------
!Warning very slow algorithm
SUBROUTINE BoundaryFaceElement_detection(gg_Face, nloadedFace, gnum_pp, boundary_N &
                                       , boundaryID, element, ndim, nod, nFace     &
                                       , nodFace, nbnd, nel_pp, npess, numpes)
  IMPLICIT NONE
  INTEGER                      :: iel, i, j, k, l;
  INTEGER                      :: BoundaryNods, BoundaryElements;
  CHARACTER(len=15)            :: face_element
  CHARACTER(len=15), INTENT(IN):: element
  INTEGER,   INTENT(IN)        :: ndim, nod, nFace, nodFace, nbnd, nel_pp;
  INTEGER,   INTENT(IN)        :: gnum_pp(nod,nel_pp), npess, numpes;
  INTEGER,   INTENT(IN)        :: boundary_N(nbnd), boundaryID;
  INTEGER,   INTENT(INOUT)     :: gg_Face(nFace+1,nel_pp), nloadedFace;
  INTEGER, ALLOCATABLE         :: BoundaryElementPos(:), FaceNodeMap(:,:);
  INTEGER, ALLOCATABLE         :: boundary_N2(:,:), gg_pp(:,:);
  INTEGER, ALLOCATABLE         :: gFace_temp(:,:), gnum_temp(:,:);
  LOGICAL                      :: FaceTest1 = .TRUE., FaceTest2 = .TRUE.;

  !-----
  ! Initialisation
  !-----
  ALLOCATE(FaceNodeMap(nFace,nodFace), gg_pp(nod,nel_pp))
  CALL Element_FaceMASK(FaceNodeMap, element, face_element, nface, nodFace, ndim)
  gg_pp = 0;
  
  !-----
  ! Flag all local boundary nodes
  !-----
  ALLOCATE(boundary_N2(nbnd,2))
  boundary_N2(:,1) = boundary_N(:);
  boundary_N2(:,2) = 0;
  ELEMENTS: DO iel = 1,nel_pp
    CALL FIND_G3(gnum_pp(:,iel),gg_pp(:,iel),boundary_N2)
  ENDDO ELEMENTS
  DEALLOCATE(boundary_N2)

  !
  ! First, second and third Pass are reduce the search space
  !

  !-----
  ! First-Pass approximating number of Elements with boundary faces
  !-----
  BoundaryElements = 0;
  DO iel = 1,nel_pp
    BoundaryNods = 0;
    DO i = 1,nod
      IF(gg_pp(i,iel)==0) BoundaryNods = BoundaryNods + 1;
    ENDDO
    IF(BoundaryNods >= nodFace) BoundaryElements = BoundaryElements + 1;
  ENDDO


  !-----
  ! Second pass Storing approximate Elements with boundary faces
  !-----
  IF(BoundaryElements > 0)THEN
    ALLOCATE(BoundaryElementPos(BoundaryElements));
    BoundaryElementPos = 0;
    k = 0;
    DO iel = 1,nel_pp
      BoundaryNods = 0;
      DO i = 1,nod
        IF(gg_pp(i,iel)==0) BoundaryNods = BoundaryNods + 1;
      ENDDO
      IF(BoundaryNods >= nodFace) THEN
        k = k + 1;
        BoundaryElementPos(k) = iel;
      ENDIF
    ENDDO
  ENDIF

  !-----
  ! Third Pass finding the element faces and setting them 
  !-----
  gg_Face = 0;
  IF(BoundaryElements > 0)THEN
    DO i = 1,BoundaryElements
      iel = BoundaryElementPos(i);
      DO j = 1,nFace
        FaceTest1 = .TRUE.
        DO k = 1,nodFace
          FaceTest2 = (gg_pp(FaceNodeMap(j,k),iel)==0);
          FaceTest1 = FaceTest2 .AND. FaceTest1;
        ENDDO
        IF(FaceTest1)THEN
          gg_Face(j+1,iel) = boundaryID;
          gg_Face(1,iel) = gg_Face(1,iel) + 1;
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  !-----
  ! Fourth pass calculate the number of loaded faces
  !-----
  nloadedFace = 0;
  DO iel = 1,nels_pp
    nloadedFace = nloadedFace + gg_Face(1,iel);
  ENDDO
  IF(BoundaryElements > 0) DEALLOCATE(BoundaryElementPos);
  DEALLOCATE(FaceNodeMap, gg_pp)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ier);
  RETURN
END SUBROUTINE BoundaryFaceElement_detection

!-------------------------------
! Calculates normal vector from tangents
!-------------------------------
SUBROUTINE NORMAL_VECTOR(normal, tangents, ndim, nvectors)
  IMPLICIT NONE
  INTEGER                 :: i, j, k, l;
  INTEGER,   INTENT(IN)   :: ndim, nvectors;
  REAL(iwp), INTENT(IN)   :: tangents(nvectors,ndim);
  REAL(iwp), INTENT(INOUT):: normal(ndim);

  !Calculate the normal vector
  normal = 0._iwp;
  SELECT CASE(ndim)
    CASE(1)
      normal(1) =  tangents(1,1);
    CASE(2)
      normal(1) =  tangents(1,2);
      normal(2) = -tangents(1,1);
    CASE(3)
      normal(1) =  tangents(1,2)*tangents(2,3) - tangents(1,3)*tangents(2,2);
      normal(2) =  tangents(1,3)*tangents(2,1) - tangents(1,1)*tangents(2,3);
      normal(3) =  tangents(1,1)*tangents(2,2) - tangents(1,2)*tangents(2,1);
    CASE DEFAULT
      WRITE(*,*) "ERROR IN NUMBER OF DIMENSIONS"
      WRITE(*,*) "only 1 to 3-D acceptable"
  END SELECT
  RETURN
END SUBROUTINE NORMAL_VECTOR

!-------------------------------
! Calculates normal vector derivative from tangents and SF-derivative
!-------------------------------
SUBROUTINE NORMAL_VECTOR_DER(NormalDer,tangents,Der,ndim,ntots,nvectors)
  IMPLICIT NONE
  INTEGER                 :: i, j, k, l;
  INTEGER,   INTENT(IN)   :: ndim, ntots, nvectors;
  REAL(iwp), INTENT(IN)   :: tangents(nvectors,ndim);
  REAL(iwp), INTENT(IN)   :: Der(nvectors,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: NormalDer(ndim,ntots);

  !Calculate the normal vector
  NormalDer = 0._iwp;
  SELECT CASE(ndim)
    CASE(1)
      NormalDer(1,:) =  tangents(1,1)*Der(1,1,:);
    CASE(2)
      NormalDer(1,:) =  tangents(1,2)*Der(1,1,:);
      NormalDer(2,:) = -tangents(1,1)*Der(1,2,:);
    CASE(3)
      NormalDer(1,:)  = NormalDer(1,:) + tangents(1,2)*Der(2,3,:);
      NormalDer(1,:)  = NormalDer(1,:) + Der(1,2,:)*tangents(2,3);
      NormalDer(1,:)  = NormalDer(1,:) - tangents(1,3)*Der(2,2,:);
      NormalDer(1,:)  = NormalDer(1,:) - Der(1,3,:)*tangents(2,2);

      NormalDer(2,:) =  NormalDer(2,:) + tangents(1,3)*Der(2,1,:);
      NormalDer(2,:) =  NormalDer(2,:) + Der(1,3,:)*tangents(2,1);
      NormalDer(2,:) =  NormalDer(2,:) - tangents(1,1)*Der(2,3,:);
      NormalDer(2,:) =  NormalDer(2,:) - Der(1,1,:)*tangents(2,3);

      NormalDer(3,:) =  NormalDer(3,:) + tangents(1,1)*Der(2,2,:);
      NormalDer(3,:) =  NormalDer(3,:) + Der(1,1,:)*tangents(2,2);
      NormalDer(3,:) =  NormalDer(3,:) - tangents(1,2)*Der(2,1,:);
      NormalDer(3,:) =  NormalDer(3,:) - Der(1,2,:)*tangents(2,1);
    CASE DEFAULT
      WRITE(*,*) "ERROR IN NUMBER OF DIMENSIONS"
      WRITE(*,*) "only 1 to 3-D acceptable"
  END SELECT
  RETURN
END SUBROUTINE NORMAL_VECTOR_DER

!-------------------------------
! Calculates normal vector to face
!-------------------------------
SUBROUTINE FACE_NORMAL_VECTOR(normal,Vcoords,Fcoords,tangents,Facefun &
                             ,nod,nodFace,ndim,nvectors)
  IMPLICIT NONE
  INTEGER                  :: i;
  INTEGER,   INTENT(IN)    :: ndim, nvectors, nod, nodFace;
  REAL(iwp), INTENT(IN)    :: tangents(nvectors,ndim), Facefun(nod);
  REAL(iwp), INTENT(IN)    :: Fcoords(nodFace,ndim), Vcoords(nod,ndim);
  REAL(iwp), INTENT(INOUT) :: normal(ndim);
  REAL(iwp), PARAMETER     :: one = 1._iwp, zero = 0._iwp;
  REAL(iwp)                :: VolCentroid(ndim),FaceCentroid(ndim)
  REAL(iwp)                :: FaceNormal(ndim), subprod

  DO I = 1,ndim
    VolCentroid(I)  = SUM(Vcoords(:,I))/REAL(nod,iwp)
    FaceCentroid(I) = DOT_PRODUCT(Facefun,Fcoords(:,I))
  ENDDO
  CALL NORMAL_VECTOR(normal, tangents, ndim, nvectors)
  FaceNormal = FaceCentroid-VolCentroid;
  subprod = DOT_PRODUCT(FaceNormal,normal)
  IF(subprod < zero) normal = -normal
  RETURN
END SUBROUTINE FACE_NORMAL_VECTOR

!-------------------------------
! Calculates normal vector to face
!-------------------------------
SUBROUTINE FACE_NORMAL_VECTOR_DER(NormalDer,Vcoords,Fcoords,Der,tangents,Facefun &
                             ,nod,nodFace,ntots,ndim,nvectors)
  IMPLICIT NONE
  INTEGER                 :: i;
  INTEGER,   INTENT(IN)   :: ndim, nvectors, nod, nodFace, ntots;
  REAL(iwp), INTENT(IN)   :: tangents(nvectors,ndim), Facefun(nod);
  REAL(iwp), INTENT(IN)   :: Fcoords(nodFace,ndim), Vcoords(nod,ndim);
  REAL(iwp), INTENT(IN)   :: Der(nvectors,ndim,ntots);
  REAL(iwp), INTENT(INOUT):: NormalDer(ndim,ntots);
  REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;
  REAL(iwp)               :: VolCentroid(ndim),FaceCentroid(ndim)
  REAL(iwp)               :: FaceNormal(ndim), normal(ndim), subprod

  DO I = 1,ndim
    VolCentroid(I)  = SUM(Vcoords(:,I))/REAL(nod,iwp)
    FaceCentroid(I) = DOT_PRODUCT(Facefun,Fcoords(:,I))
  ENDDO
  CALL NORMAL_VECTOR_DER(NormalDer,tangents,Der,ndim,ntots,nvectors)
  CALL NORMAL_VECTOR(normal, tangents, ndim, nvectors)
  FaceNormal = FaceCentroid-VolCentroid;
  subprod = DOT_PRODUCT(FaceNormal,normal)
  IF(subprod < zero) NormalDer = -NormalDer;
  RETURN
END SUBROUTINE FACE_NORMAL_VECTOR_DER

!-------------------------------
! Returns the Element face/boundary element mapping
!-------------------------------
SUBROUTINE Element_FaceMASK(EF_MAP, element, face_element, nface, nod_face, ndim)
  IMPLICIT NONE
  INTEGER                          :: i, j;
  INTEGER,           INTENT(IN)    :: nface, nod_face, ndim
  CHARACTER(len=15), INTENT(IN)    :: element
  CHARACTER(len=15), INTENT(INOUT) :: face_element
  INTEGER,           INTENT(OUT)   :: EF_MAP(nface,nod_face)

  SELECT CASE(ndim)
    CASE(2)
      face_element = 'line';
      SELECT CASE(element)
        CASE('quadrilateral')
          SELECT CASE(nod_face)
            CASE(2)
              EF_MAP(1,:) = (/1,2/);
              EF_MAP(2,:) = (/2,3/);
              EF_MAP(3,:) = (/3,4/);
              EF_MAP(4,:) = (/4,1/);
            CASE(3)
              EF_MAP(1,:) = (/1,2,3/);
              EF_MAP(2,:) = (/3,4,5/);
              EF_MAP(3,:) = (/5,6,7/);
              EF_MAP(4,:) = (/7,8,1/);
          END SELECT
        CASE('triangle')
          SELECT CASE(nod_face)
            CASE(2)
              EF_MAP(1,:) = (/1,2/);
              EF_MAP(2,:) = (/2,3/);
              EF_MAP(3,:) = (/3,1/);
            CASE(3)
              EF_MAP(1,:) = (/1,2,3/);
              EF_MAP(2,:) = (/3,4,5/);
              EF_MAP(3,:) = (/5,6,1/);
            CASE(4)
              EF_MAP(1,:) = (/1,2,3,4/);
              EF_MAP(2,:) = (/4,5,6,7/);
              EF_MAP(3,:) = (/7,8,9,1/);
            CASE(5)
              EF_MAP(1,:) = (/1,2,3,4,5/);
              EF_MAP(2,:) = (/5,6,7,8,9/);
              EF_MAP(3,:) = (/9,10,11,12,1/);
          END SELECT
      ENDSELECT
    CASE(3)
      SELECT CASE(element)
        CASE('hexahedron')
          face_element = 'quadrilateral';
          SELECT CASE(nod_face)
            CASE(4)
              EF_MAP(1,:) = (/1,5,6,2/);
              EF_MAP(2,:) = (/2,6,7,3/);
              EF_MAP(3,:) = (/3,7,8,4/);
              EF_MAP(4,:) = (/4,8,5,1/);
              EF_MAP(5,:) = (/1,2,3,4/);
              EF_MAP(6,:) = (/5,8,7,6/);
            CASE(8)
              EF_MAP(1,:) = (/1,17,5,13,6,18,2,9/);
              EF_MAP(2,:) = (/2,18,6,14,7,19,3,10/);
              EF_MAP(3,:) = (/3,19,7,15,8,20,4,11/);
              EF_MAP(4,:) = (/4,20,8,16,5,17,1,12/);
              EF_MAP(5,:) = (/1,9,2,10,3,11,4,12/);
              EF_MAP(6,:) = (/5,16,8,15,7,14,6,13/);
            CASE(9)
              EF_MAP(1,:) = (/1,17,5,13,6,18,2,9,23/);
              EF_MAP(2,:) = (/2,18,6,14,7,19,3,10,22/);
              EF_MAP(3,:) = (/3,19,7,15,8,20,4,11,24/);
              EF_MAP(4,:) = (/4,20,8,16,5,17,1,12,21/);
              EF_MAP(5,:) = (/1,9,2,10,3,11,4,12,25/);
              EF_MAP(6,:) = (/5,16,8,15,7,14,6,13,26/);
          END SELECT
        CASE('tetrahedron')
          face_element = 'triangle';
          SELECT CASE(nod_face)
            CASE(3)
              EF_MAP(1,:)  = (/1,2,3/);
              EF_MAP(2,:)  = (/1,4,2/);
              EF_MAP(3,:)  = (/2,4,3/);
              EF_MAP(4,:)  = (/3,4,1/);
            CASE(6)
              EF_MAP(1,:)  = (/1,5,2,6,3,7/);
              EF_MAP(2,:)  = (/1,8,4,9,2,5/);
              EF_MAP(3,:)  = (/2,9,4,10,3,6/);
              EF_MAP(4,:)  = (/3,10,4,8,1,7/);
          END SELECT
        CASE DEFAULT
          WRITE(*,*) "ERROR Invalid elements"
          EF_MAP = 1;
      ENDSELECT
  END SELECT
  RETURN
END SUBROUTINE Element_FaceMASK

!-------------------------------
!-------------------------------
!-------------------------------
ENDMODULE Parallel_BoundaryConditions









