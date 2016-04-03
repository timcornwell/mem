      SUBROUTINE FRPRMN(P,N,FTOL,ITER,FRET)
      PARAMETER (NMAX=1000,EPS=1.E-10)
      INTEGER ITMAX
      DIMENSION P(N),G(NMAX),H(NMAX),XI(NMAX)
      FP=FUNC(P)
      CALL DFUNC(P,XI)
      DO 11 J=1,N
        G(J)=-XI(J)
        H(J)=G(J)
        XI(J)=H(J)
11    CONTINUE
      ITMAX=ITER
      DO 14 ITS=1,ITMAX
        ITER=ITS
        CALL LINMIN(P,XI,N,FRET)
        IF(2.*ABS(FRET-FP).LE.FTOL*(ABS(FRET)+ABS(FP)+EPS))RETURN
        FP=FUNC(P)
        CALL DFUNC(P,XI)
        GG=0.
        DGG=0.
        DO 12 J=1,N
          GG=GG+G(J)**2
C          DGG=DGG+XI(J)**2
          DGG=DGG+(XI(J)+G(J))*XI(J)
12      CONTINUE
        IF(GG.LE.0.0)RETURN
        GAM=DGG/GG
        DO 13 J=1,N
          G(J)=-XI(J)
          H(J)=G(J)+GAM*H(J)
          XI(J)=H(J)
13      CONTINUE
14    CONTINUE
      PAUSE 'FRPR maximum iterations exceeded'
      RETURN
      END
C
      FUNCTION F1DIM(X)
      PARAMETER (NMAX=100)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      DIMENSION XT(NMAX)
      DO 11 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
11    CONTINUE
      F1DIM=FUNC(XT)
      RETURN
      END
C
      FUNCTION DF1DIM(X)
      PARAMETER (NMAX=100)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      DIMENSION XT(NMAX),DF(NMAX)
      DO 11 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
11    CONTINUE
      CALL DFUNC(XT,DF)
      DF1DIM=0.
      DO 12 J=1,NCOM
        DF1DIM=DF1DIM+DF(J)*XICOM(J)
12    CONTINUE
      RETURN
      END
C
      SUBROUTINE LINMIN(P,XI,N,FRET)
      PARAMETER (NMAX=100,TOL=1.E-6)
      EXTERNAL F1DIM, DF1DIM
      DIMENSION P(N),XI(N)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      NCOM=N
      DO 11 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
11    CONTINUE
      AX=0.
      XX=1.
      BX=2.
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=DBRENT(AX,XX,BX,F1DIM,DF1DIM,TOL,XMIN)
      DO 12 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
12    CONTINUE
      RETURN
      END
C
      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20)
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
1     IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 1
      ENDIF
      RETURN
      END
C
      FUNCTION DBRENT(AX,BX,CX,F,DF,TOL,XMIN)
      PARAMETER (ITMAX=100,ZEPS=1.0E-10)
      LOGICAL OK1,OK2
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DX=DF(X)
      DV=DX
      DW=DX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          D1=2.*(B-A)
          D2=D1
          IF(DW.NE.DX) D1=(W-X)*DX/(DX-DW)
          IF(DV.NE.DX) D2=(V-X)*DX/(DX-DV)
          U1=X+D1
          U2=X+D2
          OK1=((A-U1)*(U1-B).GT.0.).AND.(DX*D1.LE.0.)
          OK2=((A-U2)*(U2-B).GT.0.).AND.(DX*D2.LE.0.)
          OLDE=E
          E=D
          IF(.NOT.(OK1.OR.OK2))THEN
            GO TO 1
          ELSE IF (OK1.AND.OK2)THEN
            IF(ABS(D1).LT.ABS(D2))THEN
              D=D1
            ELSE
              D=D2
            ENDIF
          ELSE IF (OK1)THEN
            D=D1
          ELSE
            D=D2
          ENDIF
          IF(ABS(D).GT.ABS(0.5*OLDE))GO TO 1
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(DX.GE.0.) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=0.5*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
          FU=F(U)
        ELSE
          U=X+SIGN(TOL1,D)
          FU=F(U)
          IF(FU.GT.FX)GO TO 3
        ENDIF
        DU=DF(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          DV=DW
          W=X
          FW=FX
          DW=DX
          X=U
          FX=FU
          DX=DU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. W.EQ.X) THEN
            V=W
            FV=FW
            DV=DW
            W=U
            FW=FU
            DW=DU
          ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
            V=U
            FV=FU
            DV=DU
          ENDIF
        ENDIF
11    CONTINUE
      PAUSE 'DBRENT exceeded maximum iterations.'
3     XMIN=X
      DBRENT=FX
      RETURN
      END
