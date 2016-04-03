C@
C This set of MEM routines solves A*X=D+E where:
C    X is an unknown Image
C    D is some measured Data
C    E is an error term
C    A is a linear operator i.e. a matrix in this discrete representation
C
C The solution is that which:
C    Maximizes the entropy
C       H=-SUM(I) [X(I)*LN(X(I)/M(I))]
C       where M(I) is a default image
C    Fits the data such that
C       SUM(J) [(Dobs(J)-DPredicted(J))/SIGMA(J)]**2 = OMEGA
C    Has total flux
C       SUM(I) X(I) = FLUX
C
C The algorithm used is that due to Wilczek and Drapatz. The dual
C problem of finding the Lagrange multipliers is used instead of the
C primal approach of finding X(I) directly. It is easy to show that the
C MEM image is given by:
C
C       X(I) = FLUX * EXP(-SUM(J) [ A(J,I) * LAMBDA(J)] ) / Z
C
C where the partition function is:
C
C   Z = SUM(I) EXP(-SUM(J) [ A(J,I) * LAMBDA(J)] )
C
C This version calls the Numerical Recipes Conjugate Gradients routine.
C
      SUBROUTINE MMSOLVE (DSZ, ISZ, TOL, NITER, IFLUX, OMEGA, DEFI, 
     $   DOBSI, SIGMAI, MEMF, LAMBDA)
C
C Main routine for MEM solution.
C
      REAL TOL, MEMF(*), DEFI(*), DOBSI(*), SIGMAI(*), LAMBDA(*),
     $   IFLUX, OMEGA
      INTEGER DSZ, ISZ, NITER
C
      REAL FLUX
      INTEGER DATASZ, IMSZ, MEM, DEF, DOBS, DPRED, SIGMA
      COMMON /MEMCOMI/ FLUX, DATASZ, IMSZ, MEM, DEF, DOBS, DPRED,
     $   SIGMA
      REAL MM(1)
      COMMON /MEMCOM/ MM
C
      REAL MU
C
C  Set up pointers
C
      DATASZ = DSZ
      IMSZ = ISZ
      MEM = 1
      DEF = MEM + IMSZ
      DOBS = DEF + IMSZ + 1
      DPRED = DOBS + DATASZ + 1
      SIGMA = DPRED + DATASZ + 1 
C
      FLUX = IFLUX
      DO 1 J = 1, DATASZ
         MM(DOBS+J-1) = DOBSI(J)
         MM(SIGMA+J-1) = SIGMAI(J)
 1    CONTINUE
C
      MU = 0.0
      IF(OMEGA.GT.0.0) THEN
         MM(DOBS+DATASZ) = OMEGA
         DO 100 J = 1, DATASZ
            MU = MU + (MM(DOBS+J-1)*MM(SIGMA+J-1))**2
 100     CONTINUE
         MU = SQRT(MU/(4.0*OMEGA))
      ENDIF
      LAMBDA(DATASZ+1) = MU
C
      DO 10 I = 1, IMSZ
         MM(DEF+I-1) = DEFI(I)
 10   CONTINUE
C
      CALL FRPRMN (LAMBDA, DATASZ, TOL, NITER, Z)
C
      DO 20 I = 1, IMSZ
         MEMF(I) = MM(MEM+I-1) 
 20   CONTINUE
C
      END
C
      SUBROUTINE F (LAMBDA, Z)
      REAL LAMBDA(*), Z
      INTEGER I
      REAL MAXDR, LAMTRM, NORM
      DATA MAXDR /1E38/
      REAL FLUX
      INTEGER DATASZ, IMSZ, MEM, DEF, DOBS, DPRED, SIGMA
      COMMON /MEMCOMI/ FLUX, DATASZ, IMSZ, MEM, DEF, DOBS, DPRED,
     $   SIGMA
      REAL MM(1)
      COMMON /MEMCOM/ MM
C
C Find A^T * Lambda
C
      CALL MMATD (DATASZ, IMSZ, LAMBDA, MM(MEM))
C
C Find new Image and partition function
C
      LAMTRM = LOG(MAXDR)
      NORM = 0.0
      DO 10 I = 1, IMSZ
         MM(MEM+I-1) = MM(DEF+I-1) * EXP(MIN(-MM(MEM+I-1),LAMTRM))
         NORM = NORM + MM(MEM+I-1)
  10  CONTINUE
      DO 20 I = 1, IMSZ
         MM(MEM+I-1) = FLUX * MM(MEM+I-1) / NORM
  20  CONTINUE
C
      Z = FLUX*LOG(NORM)
C
      END
C
      REAL FUNCTION FUNC(LAMBDA)
      REAL LAMBDA(*)
      INTEGER J
      REAL Z, MU, OMEGA
C
      REAL FLUX
      INTEGER DATASZ, IMSZ, MEM, DEF, DOBS, DPRED, SIGMA
      COMMON /MEMCOMI/ FLUX, DATASZ, IMSZ, MEM, DEF, DOBS, DPRED,
     $   SIGMA
      REAL MM(1)
      COMMON /MEMCOM/ MM
C
C Find Image for this set of Lagrange Multipliers
C
      CALL F (LAMBDA, Z)
      DO 30 J = 1, DATASZ
         Z = Z + LAMBDA(J) * MM(DOBS+J-1)
   30 CONTINUE
      MU = LAMBDA(DATASZ+1)
      IF(MU.GT.0.0) THEN
         OMEGA = MM(DOBS+DATASZ)
         Z = Z + MU * OMEGA 
         DO 40 J = 1,DATASZ
            Z = Z + (LAMBDA(J)*MM(SIGMA+J-1))**2 / (4.0*MU)
   40    CONTINUE
      ENDIF
      FUNC = Z
      RETURN
      END
C
      SUBROUTINE DFUNC (LAMBDA, GRADZ)
      REAL LAMBDA(*), GRADZ(*)
C
      REAL FLUX
      INTEGER DATASZ, IMSZ, MEM, DEF, DOBS, DPRED, SIGMA
      COMMON /MEMCOMI/ FLUX, DATASZ, IMSZ, MEM, DEF, DOBS, DPRED,
     $   SIGMA
      REAL MM(1)
      COMMON /MEMCOM/ MM
C
      REAL MU, Z
      INTEGER J
C
C Find current Image
C
      CALL F (LAMBDA, Z)
C
C Find predicted data for this image using A*X
C
      CALL MMAX (IMSZ, DATASZ, MM(MEM), MM(DPRED))
C
C Now calculate gradient
C
      DO 10 J = 1, DATASZ
         GRADZ(J) = (MM(DOBS+J-1) - MM(DPRED+J-1))
 10   CONTINUE
C
      MU=LAMBDA(DATASZ+1)
      IF(MU.GT.0.0) THEN
         GRADZ(DATASZ+1)=0.0
         DO 20 J = 1, DATASZ
            GRADZ(DATASZ+1)=GRADZ(DATASZ+1)
     $         + (LAMBDA(J)*MM(SIGMA+J-1))**2
            GRADZ(J) = GRADZ(J)
     $         + LAMBDA(J)*MM(SIGMA+J-1)**2/(2.0*MU) 
 20      CONTINUE
      ENDIF
C
      END
C
