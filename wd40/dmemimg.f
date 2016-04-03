      PROGRAM DMEMIMG
C
C Driver for MEM routines. Solves circular convolution.
C
      INTEGER MAXIMSZ
      PARAMETER(MAXIMSZ=100)
C
      INTEGER J, NITER, IMSZ, DATASZ
      REAL TOL, MODEL(MAXIMSZ), DOBS(MAXIMSZ)
      REAL MEM(MAXIMSZ), DEF(MAXIMSZ), SIGMA(MAXIMSZ), LAMBDA(MAXIMSZ)
C
      REAL MM(6*MAXIMSZ)
      COMMON /MEMCOM/ MM
C
      WRITE(*,*) 'Maximum entropy deconvolution program'
C
      TOL = 1E-20
      DATASZ = 64
      IMSZ = DATASZ
      NITER=100
C
C Initialize model data
C
      DO 1 I=1,IMSZ
         MODEL(I)=0.0
    1 CONTINUE
      MODEL(33)=900.0
      DO 2 I=34,48
         MODEL(I)=20.0
    2 CONTINUE
      MODEL(49)=200.0
C
C Find total flux
C
      FLUX = 0.0
      DO 3 I=1,IMSZ
         FLUX=FLUX+MODEL(I)
 3    CONTINUE
C
C Initialize SIGMA
C
      DO 4 J=1,DATASZ
         SIGMA(J)=1.0E-6
 4    CONTINUE
C
C Find predicted data
C
      CALL MMAX (IMSZ, DATASZ, MODEL, DOBS)
C
C Set up default image
C
      DO 9 I = 1,IMSZ
         DEF(I) = FLUX/FLOAT(IMSZ)
 9    CONTINUE
C
C Actually do iterations
C
      CALL MMSOLVE (DATASZ, IMSZ, TOL, NITER, FLUX, FLOAT(DATASZ), 
     $   DEF, DOBS, SIGMA, MEM, LAMBDA)
C
      WRITE(*,1000)
      WRITE(*,1010)
      DO 50 I = 1, IMSZ
         WRITE(6,1020) I, MODEL(I), DEF(I), DOBS(I), MEM(I),
     $      MEM(I)-MODEL(I)
 50   CONTINUE
C
 1000 FORMAT (1X,'PixelNo',7X,'Model',5X,'a priori',3X,
     &        'DOBS',5X,'Image',6X,'Residuals')
 1010 FORMAT ('---------------------------------------------',
     &         '--------------------')
 1020 FORMAT (2X,I2,10X,F5.1,6X,F5.1,5X,F5.1,4X,F6.2,4X,F8.3)
C
      END
C----------------------------------------------------------------------
      SUBROUTINE MMAX(IMSZ, DATASZ, IN, OUT)
C
C Returns A * X i.e. OUT(J)=SUM(I) [A(I,J) * IN(I)]
C
      INTEGER DATASZ, IMSZ
      REAL IN(*), OUT(*)
      CALL CONVOLVEARRAY (DATASZ, IN, OUT)
C
      END
C----------------------------------------------------------------------
      SUBROUTINE MMATD(DATASZ, IMSZ, IN, OUT)
C
C Returns A^T * D i.e. OUT(I)=SUM(J) [A(J,I) * IN(J)]
C
      INTEGER DATASZ, IMSZ
      REAL IN(*), OUT(*)
      CALL CONVOLVEARRAY (DATASZ, IN, OUT)
C
      END
C---------------------------------------------------------
      SUBROUTINE CONVOLVEARRAY(N,S,T)
C
C  Convolve circularly 
C
      REAL T(*),S(*),TEMP(512)
      INTEGER N,J
      TEMP(1)=(S(N-1)+S(N)+S(1)+S(2)+S(3))
      TEMP(2)=(S(N)+S(1)+S(2)+S(3)+S(4))
      DO 1  J=3,N-2
    1 TEMP(J)=(S(J-2)+S(J-1)+S(J)+S(J+1)+S(J+2))
      TEMP(N-1)=(S(N-3)+S(N-2)+S(N-1)+S(N)+S(1))
      TEMP(N)=(S(N-2)+S(N-1)+S(N)+S(1)+S(2))
      DO 2 J=1,N
    2 T(J) = 0.2*TEMP(J)
      END
