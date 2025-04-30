

    J11 = J11 + MATMUL(TRANSPOSE(dE),MATMUL(C_ijkl,dE))
    DO s =1,nst
      R1  = R1  + S_ij(s)*dE(s,:)
      J11 = J11 + S_ij(s)*d2E(s,:,:)
    ENDDO

    CALL SHAPE_FUN(funP,points2,Ig);
    funnyP(1,:) = funP(abaqustosg)
    Cdef = MATMUL(TRANSPOSE(Fedef),Fedef);
    CALL INVERT2(Cdef,CdefInv,ndim);
    J3 = determinant(Fedef)
    pressure = DOT_PRODUCT(funP(abaqustosg),ptemp)
    DO s =1,nst
      CALL VOIGHT_ITERATOR(s,i,j,nst);
      R1  = R1  + pressure*CdefInv(i,j)*dE(s,:)
      J11 = J11 + pressure*CdefInv(i,j)*d2E(s,:,:)
      DO m = 1,ndofU;
        J12(m,:)=J12(m,:) + CdefInv(i,j)*funnyP(1,:)*dE(s,m)
        DO n = 1,ndofU; DO t =1,nst
          CALL VOIGHT_ITERATOR(t,k,l,nst);
          J11(m,n)=J11(m,n)-2._iwp*pressure*CdefInv(i,k)*CdefInv(j,l) &
                                  *dE(s,m)*dE(t,n)
        ENDDO; ENDDO
      ENDDO
    ENDDO

  Y     = 1.00_iwp;
  nu    = 0.49999_iwp;
  mu    = (4.6_iwp/2.2_iwp)/(Y/(2._iwp+2._iwp*nu))
  lmbda = Y*nu/((1+nu)*(1._iwp-2._iwp*nu))


    R2  = R2  + (DLOG(J3) - (pressure/lmbda))*funnyP(1,:)
    J22 = (-J3/lmbda)*MATMUL(TRANSPOSE(funnyP(:,:)),funnyP(:,:))


  !-----
  ! Integrate the residuals
  ! and the Jacobians
  !-----
  DO M = 1,ntots
    GAUSS_PTS2: DO Ig = 1,nip
      DO I = 1,nst
        Rm(M) = Rm(M) + dE(I,M,Ig)*S_ij(I,Ig)*det(Ig)*weights(Ig);
        Km(M,:) = Km(M,:) + S_ij(I,Ig)*d2E(I,:,M,Ig)*det(Ig)*weights(Ig);
        DO J = 1,nst
          Km(M,:) = Km(M,:) + dE(I,M,Ig)*C_ijkl(I,J,Ig)*dE(J,:,Ig)*det(Ig)*weights(Ig);
        ENDDO
      ENDDO
    END DO GAUSS_PTS2
  ENDDO


  !
  ! Residual and Jacobian eliminated Pressure DOFs
  !
  DO i = (ndofU+ndofP+1),ntots
    Rm(I)   = zero;
    Km(I,:) = zero;
    Km(I,I) = one;
  ENDDO
