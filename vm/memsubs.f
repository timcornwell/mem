C----------------------------------------------------------------------
C
C  Maximum Entropy Deconvolution Program
C
C     T.J. Cornwell
C     National Radio Astronomy Observatory
C     P.O. Box 0,
C     Socorro, New Mexico, 87801, USA.
C
C (C) Copyright
C     Associated Universities, Inc.
C     1990
C
C----------------------------------------------------------------------
      SUBROUTINE MEMONEIT (DOINIT, IMSZ, MEM, DEF, DEFLEV, STP, GCH, WT,
     1   Q, TFLUX, TCHISQ, TOL, FLUX, CHISQ, ENTRPY, NRMGRD, CALGCH)
C
C Does one MEM iteration. In this version, the pixels must be contiguous
C in memory. The routine CALGCH is called to evaluate both chi-squared
C and the gradient of chi-squared with respect to the image. You should
C call this routine repeatedly until FLUX and TFLUX, CHISQ and TCHISQ
C agree and NRMGRD is reasonably small e.g. < 0.05 say. TOL can be
C about 0.3 for reasonable speed in most cases. For difficult cases set
C it to be much smaller. Q is an estimate of the diagonal element of 
C grad-grad chi-squared, and should be initially set to some reasonable
C value. IF DEFLEV is greater than zero then it is taken as a constant
C default level and the default image is ignored.
C     This routine performs the MEM algorithm of Cornwell and Evans
C with some modifications suggested by Bob Sault.
C
C       DOINIT	LOG	input	Initialize ?
C	IMSZ	INT	input	Number of pixels
C	MEM	REAL(*)	input	MEM image
C	DEF	REAL(*)	input	default image
C	DEFLEV	REAL	input	Constant default level
C	STP	REAL(*)	input	Name of current step image
C	GCH	REAL(*)	input	Name of Grad chi-squared image
C	WT	REAL(*)	input	Weights (1/sigma**2)
C	Q	REAL	i/o	Current estimate of diagonal
C	TFLUX	REAL	input	Required Total flux
C	TCHISQ	REAL	input	Required chi-squared
C	TOL	REAL	input	Tolerance for solution
C	FLUX	REAL	output	Actual current flux
C	CHISQ	REAL	output	Actual current chi-squared
C	ENTRPY	REAL	output	Entropy of image
C	NRMGRD	REAL	output	Normalized gradient
C	CALGCH	EXT	input	Routine to evaluate chi-squared and 
C				the gradient:
C				CALL CALGCH(IMSZ, MEM, WT, GCH, CHISQ)
C----------------------------------------------------------------------
      LOGICAL 		DOINIT
      INTEGER		IMSZ
      REAL		MEM(*), DEF(*), STP(*), GCH(*), WT(*)
Cf2py intent(in) imsz
Cf2py intent(inout) mem
Cf2py intent(in) def
Cf2py intent(in) deflev
Cf2py intent(inout) stp
Cf2py intent(in) gch
Cf2py intent(in) wt
Cf2py depend(imsz) mem
Cf2py depend(imsz) def
Cf2py depend(imsz) stp
Cf2py depend(imsz) gsh
Cf2py depend(imsz) wt
      REAL		DEFLEV, Q, TFLUX, TCHISQ, TOL, FLUX, CHISQ
      REAL		ENTRPY, NRMGRD
Cf2py intent(inout) q
Cf2py intent(in) tflux
Cf2py intent(in) tchisq
Cf2py intent(in) tol
Cf2py intent(out) flux
Cf2py intent(out) chisq
Cf2py intent(out) entrpy
Cf2py intent(out) nrmgrd
C
      REAL		SCALE, SCALEM, EPS, PONE, PZERO, IMLO, 
     1			IMMAX, IMMIN
      REAL		ALPHA, BETA
      DOUBLE PRECISION	GDG(4,4)
      SAVE		ALPHA, BETA
C======================================================================
C
C If this the first call
C
      IF (DOINIT) THEN
         ALPHA = 0.0
         BETA = 0.0
         FLUX = 0.0
         CALL CALGCH (IMSZ, MEM, WT, GCH, CHISQ)
         CALL MEMCHALB (IMSZ, MEM, DEF, DEFLEV, GCH, WT, ALPHA, BETA,
     $      Q, FLUX, TFLUX, CHISQ, TCHISQ, TOL, NRMGRD, GDG)
      END IF
C
C Find Step to take
C
      CALL MEMCLSTP (IMSZ, MEM, DEF, DEFLEV, STP, GCH, WT, ALPHA, BETA,
     $   Q, FLUX, IMMAX, IMMIN, NRMGRD, GDG, PZERO)
C
C Limit the step to less than the tolerance
C
      SCALE = 1.0
      SCALEM = 1.0
      IF (NRMGRD.GT.0.0) SCALEM = TOL/NRMGRD
      SCALE = MIN(1.0, SCALEM)
C
C Take step
C
      IMLO = 0.1*IMMIN + 1E-8
      CALL MEMTAKST (IMSZ, MEM, STP, 1.0, SCALE, IMLO, MEM)
C
C  Now calculate the gradient after the step to see if it was a good idea.
C  The value of the gradient will be used to revise the step length.
C  First, calculate the residuals for this new image
C
      CALL CALGCH (IMSZ, MEM, WT, GCH, CHISQ)
C
C  Calculate dot product of gradient and step
C
       CALL MEMGDS  (IMSZ, MEM, DEF, DEFLEV, STP, GCH, ALPHA, BETA,
     $   PONE)
C
C  This is the optimum step
C
       EPS = 1.0
       IF (PZERO.NE.PONE) EPS = PZERO/(PZERO-PONE)
       IF (SCALE.NE.0.0) EPS = MIN(EPS, SCALEM/SCALE)
       IF (EPS.LE.0.0) THEN
          EPS = 1.0
       END IF
C
C  now step to the optimum point, clipping the step to prevent excessive
C  change on any one iteration
C
      IF (ABS(EPS-1.0).GT.TOL) THEN
         CALL MEMTAKST (IMSZ, MEM, STP, 1.0, SCALE*(EPS-1.0), IMLO, MEM)
         CALL CALGCH (IMSZ, MEM, WT, GCH, CHISQ)
      END IF
C
C  calculate entropy and flux of new image
C
      CALL MEMCALHF (IMSZ, MEM, DEF, DEFLEV, FLUX, ENTRPY)
C
C Re-adjust the estimate of the beam volume.
C
      Q = Q*(1/MAX(0.5,MIN(2.0,EPS))+1.0)/2.0
C
C Change alpha and beta
C
      CALL MEMCHALB (IMSZ, MEM, DEF, DEFLEV, GCH, WT, ALPHA, BETA,
     $   Q, FLUX, TFLUX, CHISQ, TCHISQ, TOL, NRMGRD, GDG)

      END
C***********************************************************************
C
C The following are internal routines and should not be called directly
C
C***********************************************************************
C
      SUBROUTINE MEMCALHF(IMSZ, IN1, IN2, DEFLEV, FLUX, ENTRPY)
C
C Calculate entropy, flux, etc
C
C
C IMSZ       INT  I  ARRAY SIZE
C IN1(*)     REAL I  ARRAY 1
C IN2(*)     REAL I  ARRAY 2
C DEFLEV     REAL I  CONSTANT DEFAULT LEVEL
C FLUX       REAL O  FLUX OF IN1
C ENTRPY     REAL O  ENTRPY OF IN1 RELATIVE TO IN2
C-------------------------------------------------------------------------
C
      INTEGER 	IMSZ
      REAL	IN1(*), IN2(*), DEFLEV
      REAL	FLUX, ENTRPY
C
      INTEGER 	IMCNTR
C======================================================================
C
      ENTRPY = 0.0
      FLUX = 0.0
C
      IF (DEFLEV.GT.0.0) THEN
         DO 10 IMCNTR = 1,IMSZ
           FLUX = FLUX + IN1(IMCNTR)
           ENTRPY = ENTRPY - IN1(IMCNTR)*LOG(IN1(IMCNTR))
  10     CONTINUE
         ENTRPY = ENTRPY + FLUX*LOG(DEFLEV)
      ELSE
         DO 20 IMCNTR = 1,IMSZ
           FLUX = FLUX + IN1(IMCNTR)
           ENTRPY = ENTRPY - IN1(IMCNTR)*LOG(IN1(IMCNTR)/IN2(IMCNTR))
  20     CONTINUE
      END IF
C
      IF (FLUX.GT.0.0) THEN
         ENTRPY = ENTRPY/FLUX + LOG(FLOAT(IMSZ))
      ELSE
         ENTRPY = 0.0
      END IF
C
      END
      SUBROUTINE MEMCHALB (IMSZ, MEM, DEF, DEFLEV, GCH, WT, ALPHA, BETA,
     $   Q, FLUX, TFLUX, CHISQ, TCHISQ, TOL, NRMGRD, GDG)
C
C Find changes in alpha and beta
C
C  IMSZ       INT    I    IMAGE SIZE
C  MEM(*)     REAL   I    CURRENT SOLUTION
C  DEF(*)     REAL   I    DEFAULT IMAGE
C  DEFLEV     REAL   I    CONSTANT DEFAULT LEVEL
C  GCH(*)     REAL   I    GRAD CHI-SQUARED
C  WT(*)      REAL   I    Weights
C  ALPHA      REAL   I    LAGRANGE MULTIPLIER FOR CHISQ
C  BETA       REAL   I    LAGRANGE MULTIPLIER FOR FLUX
C  Q          REAL   IO   DIAG. VALUE GRAD GRAD CHISQ
C  TOL        REAL   I    TOLERANCE FOR SOLUTIION
C  NRMGRD     REAL   O    NORMALISED GRADIENT
C  GDG(4,4)   REAL   O    GRADIENT DOT PRODUCTS
C  GDSTEP     REAL   O    CURRENT GRADIENT . STEP
C-------------------------------------------------------------------------
C
      INTEGER 	IMSZ
      REAL	MEM(*), DEF(*), DEFLEV, GCH(*), WT(*), ALPHA, BETA, Q
      REAL	TOL, TFLUX, FLUX, CHISQ, TCHISQ, NRMGRD
      DOUBLE PRECISION	GDG(4,4)
C
      INTEGER H,C,F,J
      PARAMETER (H=1)
      PARAMETER (C=2)
      PARAMETER (F=3)
      PARAMETER (J=4)
C
      INTEGER 	IMCNTR
      INTEGER 	AXIS1,AXIS2
C
      REAL	RHESS, GRAD(J), GGC, LENGTH, LOGDFLEV
C======================================================================
C
C Initialise
C
       DO 1 AXIS2 = H,J
       DO 2 AXIS1 = H,J
          GDG(AXIS1,AXIS2) = 0.0
   2   CONTINUE
   1   CONTINUE
       GGC = ALPHA*Q
C
C Now do entire image
C
       IF (DEFLEV.GT.0.0) THEN
          LOGDFLEV = LOG(DEFLEV)
          DO 10 IMCNTR = 1,IMSZ
            RHESS = MEM(IMCNTR)/(1+GGC*WT(IMCNTR)*MEM(IMCNTR))
            GRAD(H) = -LOG(MEM(IMCNTR)) + LOGDFLEV
            GRAD(C) = GCH(IMCNTR)
            GDG(H,H) = GDG(H,H)+ GRAD(H)*RHESS*GRAD(H)
            GDG(H,C) = GDG(H,C)+ GRAD(H)*RHESS*GRAD(C)
            GDG(H,F) = GDG(H,F)+ GRAD(H)*RHESS
            GDG(C,C) = GDG(C,C)+ GRAD(C)*RHESS*GRAD(C)
            GDG(C,F) = GDG(C,F)+ GRAD(C)*RHESS
            GDG(F,F) = GDG(F,F)+ RHESS
  10     CONTINUE
      ELSE
         DO 20 IMCNTR = 1,IMSZ
            RHESS = MEM(IMCNTR)/(1+GGC*WT(IMCNTR)*MEM(IMCNTR))
            GRAD(H) = -LOG(MEM(IMCNTR)/DEF(IMCNTR))
            GRAD(C) = GCH(IMCNTR)
            GDG(H,H) = GDG(H,H)+ GRAD(H)*RHESS*GRAD(H)
            GDG(H,C) = GDG(H,C)+ GRAD(H)*RHESS*GRAD(C)
            GDG(H,F) = GDG(H,F)+ GRAD(H)*RHESS
            GDG(C,C) = GDG(C,C)+ GRAD(C)*RHESS*GRAD(C)
            GDG(C,F) = GDG(C,F)+ GRAD(C)*RHESS
            GDG(F,F) = GDG(F,F)+ RHESS
  20     CONTINUE
      END IF
C
C Evaluate the gradient, normalised appropriately. 
C
      GDG (H,J) = GDG(H,H) - ALPHA * GDG(H,C) - BETA * GDG (H,F)
      GDG (C,J) = GDG(H,C) - ALPHA * GDG(C,C) - BETA * GDG (C,F)
      GDG (F,J) = GDG(H,F) - ALPHA * GDG(C,F) - BETA * GDG (F,F)
      GDG (J,J) = GDG(H,H) + ALPHA**2 *GDG(C,C) + BETA**2 * GDG(F,F)
     1   + 2 * ALPHA * BETA * GDG (C,F) - 2 * ALPHA * GDG (H,C)
     2   - 2 * BETA * GDG (H,F)      
      LENGTH = GDG(H,H) + ALPHA**2 * GDG(C,C) + BETA**2 * GDG(F,F)
C
      IF ((ALPHA.EQ.0.0).AND.(BETA.EQ.0.0)) THEN
         LENGTH = GDG(F,F)
      END IF
      NRMGRD = GDG(J,J)/LENGTH
      IF (ALPHA.EQ.0.0) NRMGRD = 0.0
C
      IF (NRMGRD.LE.TOL) THEN
         CALL MEMUPDAT (GDG, TOL, TFLUX, TCHISQ, FLUX, CHISQ, LENGTH,
     $      ALPHA, BETA)
      ELSE 
         CALL MEMINIAB (GDG, ALPHA, BETA, TFLUX)
      END IF
C
      END
      SUBROUTINE MEMCLSTP(IMSZ, MEM, DEF, DEFLEV, STP, GCH, WT, ALPHA, 
     1   BETA, Q, FLUX, IMMAX, IMMIN, NRMGRD, GDG, GDSTEP)
C
C Calculate the next change in the mem image from the usual
C Newton-raphson formula. the inverse of the hessian is approximated by
C only taking the inverse of the diagonal elements
C
C  IMSZ       INT    I    IMAGE SIZE
C  MEM(*)     REAL   I    CURRENT SOLUTION
C  DEF(*)     REAL   I    DEFAULT IMAGE
C  DEFLEV     REAL   I    CONSTANT DEFAULT LEVEL
C  GCH(*)     REAL   I    Grad Chisquared
C  WT(*)      REAL   I    Weightsp
C  STP(*)     REAL   IO   STEP IMAGE
C  ALPHA      REAL   I    LAGRANGE MULTIPLIER FOR CHISQ
C  BETA       REAL   I    LAGRANGE MULTIPLIER FOR FLUX
C  FLUX       REAL   O    ACTUAL FLUX
C  Q          REAL   IO   DIAG. VALUE GRAD GRAD CHISQ
C  IMMAX      REAL   O    IMAGE MAX.
C  IMMIN      REAL   O    IMAGE MIN.
C  NRMGRD     REAL   O    NORMALISED GRADIENT
C  GDG(4,4)   DBLE   O    GRADIENT DOT PRODUCTS
C  GDSTEP     REAL   O    CURRENT GRADIENT . STEP
C
C-------------------------------------------------------------------------
C
      INTEGER 	IMSZ
      REAL	MEM(*), DEF(*), DEFLEV, GCH(*), WT(*), STP(*)
      REAL	IMMIN, IMMAX, NRMGRD, GDSTEP, FLUX, Q, ALPHA, BETA
      DOUBLE PRECISION	GDG(4,4)
C
      INTEGER H,C,F,J
      PARAMETER (H=1)
      PARAMETER (C=2)
      PARAMETER (F=3)
      PARAMETER (J=4)
C
      INTEGER 	IMCNTR
      INTEGER 	AXIS1,AXIS2
C
      REAL	RHESS,GRAD(J),GGC, LENGTH, LOGDFLEV
C======================================================================
C
C Initialise
C
       IMMAX = -1E20
       IMMIN =  1E20
       FLUX = 0.0
       DO 1 AXIS2 = H,J
       DO 2 AXIS1 = H,J
          GDG(AXIS1,AXIS2) = 0.0
   2   CONTINUE
   1   CONTINUE
       GGC = ALPHA*Q
C
C Now do entire image
C
       IF (DEFLEV.GT.0.0) THEN
          LOGDFLEV = LOG(DEFLEV)
          DO 10 IMCNTR = 1,IMSZ
            FLUX = FLUX + MEM(IMCNTR)
            IMMIN = MIN (IMMIN, MEM(IMCNTR))
            IMMAX = MAX (IMMAX, MEM(IMCNTR))
            RHESS = MEM(IMCNTR)/(1+GGC*WT(IMCNTR)*MEM(IMCNTR))
            GRAD(H) = -LOG(MEM(IMCNTR)) + LOGDFLEV
            GRAD(C) = GCH(IMCNTR)
            GRAD(J) = GRAD(H) - ALPHA*GRAD(C) - BETA
            STP (IMCNTR) = RHESS*GRAD(J)
            GDG(H,H) = GDG(H,H)+ GRAD(H)*RHESS*GRAD(H)
            GDG(H,C) = GDG(H,C)+ GRAD(H)*RHESS*GRAD(C)
            GDG(H,F) = GDG(H,F)+ GRAD(H)*RHESS
            GDG(C,C) = GDG(C,C)+ GRAD(C)*RHESS*GRAD(C)
            GDG(C,F) = GDG(C,F)+ GRAD(C)*RHESS
            GDG(F,F) = GDG(F,F)+ RHESS
  10     CONTINUE
      ELSE
         DO 20 IMCNTR = 1,IMSZ
            FLUX = FLUX + MEM(IMCNTR)
            IMMIN = MIN (IMMIN, MEM(IMCNTR))
            IMMAX = MAX (IMMAX, MEM(IMCNTR))
            RHESS = MEM(IMCNTR)/(1+GGC*WT(IMCNTR)*MEM(IMCNTR))
            GRAD(H) = -LOG(MEM(IMCNTR)/DEF(IMCNTR))
            GRAD(C) = GCH(IMCNTR)
            GRAD(J) = GRAD(H) - ALPHA*GRAD(C) - BETA
            STP (IMCNTR) = RHESS*GRAD(J)
            GDG(H,H) = GDG(H,H)+ GRAD(H)*RHESS*GRAD(H)
            GDG(H,C) = GDG(H,C)+ GRAD(H)*RHESS*GRAD(C)
            GDG(H,F) = GDG(H,F)+ GRAD(H)*RHESS
            GDG(C,C) = GDG(C,C)+ GRAD(C)*RHESS*GRAD(C)
            GDG(C,F) = GDG(C,F)+ GRAD(C)*RHESS
            GDG(F,F) = GDG(F,F)+ RHESS
  20     CONTINUE
      END IF
C
C Evaluate the gradient, normalised appropriately. 
C
      GDG (H,J) = GDG(H,H) - ALPHA * GDG(H,C) - BETA * GDG (H,F)
      GDG (C,J) = GDG(H,C) - ALPHA * GDG(C,C) - BETA * GDG (C,F)
      GDG (F,J) = GDG(H,F) - ALPHA * GDG(C,F) - BETA * GDG (F,F)
      GDG (J,J) = GDG(H,H) + ALPHA**2 *GDG(C,C) + BETA**2 * GDG(F,F)
     1   + 2 * ALPHA * BETA * GDG (C,F) - 2 * ALPHA * GDG (H,C)
     2   - 2 * BETA * GDG (H,F)      
      GDSTEP = GDG(J,J)
      LENGTH = GDG(H,H) + ALPHA**2 * GDG(C,C) + BETA**2 * GDG(F,F)
      IF (LENGTH.LE.0.0) THEN
         LENGTH = GDG(F,F)
      END IF
      NRMGRD = GDG(J,J)/LENGTH
C
      END
      SUBROUTINE MEMGDS (IMSZ, MEM, DEF, DEFLEV, STP, GCH, ALPHA, 
     1   BETA, GDSTEP)
C
C Evaluate the dot product of the current image and the step
C
C  IMSZ       INT    I    IMAGE SIZE
C  MEM(*)     REAL   I    CURRENT SOLUTION
C  DEF(*)     REAL   I    DEFAULT IMAGE
C  GCH(*)     REAL   I    GRAD CHI-SQUARED
C  STP(*)     REAL   I    STEP IMAGE
C  ALPHA      REAL   IO   LAGRANGE MULTIPLIER FOR CHISQ
C  BETA       REAL   IO   LAGRANGE MULTIPLIER FOR FLUX
C  GDSTEP     REAL   O    CURRENT GRAD . STEP
C-------------------------------------------------------------------------
C
      INTEGER 	IMSZ
      REAL	MEM(*), DEF(*), DEFLEV, GCH(*), STP(*), ALPHA, BETA
      REAL	GDSTEP
C
      INTEGER 	IMCNTR
      REAL	LOGDFLEV
C======================================================================
C
C Initialise
C
      GDSTEP = 0.0
C
C Now do entire image
C
       IF (DEFLEV.GT.0.0) THEN
          LOGDFLEV = LOG(DEFLEV)
          DO 10 IMCNTR = 1,IMSZ
            GDSTEP = GDSTEP + STP(IMCNTR) * 
     1        (-LOG(MEM(IMCNTR)) + LOGDFLEV - ALPHA*GCH(IMCNTR) 
     2         - BETA)
  10     CONTINUE
      ELSE
         DO 20 IMCNTR = 1,IMSZ
            GDSTEP = GDSTEP + STP(IMCNTR) * 
     1         (-LOG(MEM(IMCNTR)/DEF(IMCNTR)) - ALPHA*GCH(IMCNTR) 
     2          - BETA)
  20     CONTINUE
      END IF
C
      END
      SUBROUTINE MEMINIAB (GDG, ALPHA, BETA, TFLUX)
C
C Change alpha as much as required and allowed. the changes
C In alpha are such that the resulting normalised gradient does
C Not exceed the tolerance
C
C
C GDG(4, 4)    DBLE  I  GRADIENT DOT PRODUCTS
C ALPHA        REAL  O  LAGRANGE MULTIPLIER FOR FLUX
C BETA         REAL  O  LAGRANGE MULTIPLIER FOR CHISQ
C-------------------------------------------------------------------------
C
      REAL	ALPHA, BETA, TFLUX
      DOUBLE PRECISION	GDG(4, 4)
C
      INTEGER H, C, F, J
      PARAMETER (H=1)
      PARAMETER (C=2)
      PARAMETER (F=3)
      PARAMETER (J=4)
      DOUBLE PRECISION	DET
C======================================================================

      IF (TFLUX.LE.0.0) THEN
         ALPHA = MAX (0.0, SNGL(GDG(H, C)/GDG(C, C)))
         BETA = 0.0
      ELSE
         DET = GDG(C, C)*GDG(F, F) - GDG(C, F)*GDG(C, F)
         ALPHA = (GDG(F, F)*GDG(H, C)-GDG(C, F)*GDG(H, F))/DET
         BETA  = (GDG(C, C)*GDG(H, F)-GDG(C, F)*GDG(H, C))/DET
      END IF
C
      END
      SUBROUTINE MEMTAKST (N, IN1, IN2, WT1, WT2, IMLO, OUT)
C
C Add two arrays together with different weights, also clip.
C
C
C  N     INT  I  ARRAY SIZE
C  IN1      REAL I  ARRAY 1
C  IN2      REAL I  ARRAY 2
C  WT1      REAL I  WEIGHT FOR ARRAY 1
C  WT2      REAL I  WEIGHT FOR ARRAY 2
C  IMLO     REAL I  MINIMUM ALLOWED FOR OUT
C  OUT      REAL O  OUTPUT ARRAY
C
C-------------------------------------------------------------------------
C
      INTEGER 	N, I
      REAL	IN1(N), IN2(N), OUT(N)
      REAL	WT1, WT2, IMLO
C======================================================================
C
      DO 10 I = 1,N
        OUT(I) = MAX ( (WT1 * IN1(I) + WT2 * IN2(I)), IMLO)
  10  CONTINUE
C
      END
      SUBROUTINE MEMUPDAT (GDG, TOL, TFLUX, TCHISQ, FLUX, CHISQ, 
     &   LENGTH, ALPHA, BETA)
C
C  Change alpha as much as required and allowed. the changes
C  in alpha are such that the resulting normalised gradient does
C  not exceed the tolerance
C
C GDG(4,4)     DBLE  I  GRADIENT DOT PRODUCTS
C TOL          REAL  I  TOLERANCE
C TFLUX        REAL  I  TARGET FLUX
C TCHISQ       REAL  I  TARGET CHISQ
C FLUX         REAL  I  ACTUAL FLUX
C CHISQ        REAL  I  ACTUAL CHISQ
C LENGTH       REAL  I  LENGTH
C ALPHA        REAL  O  LAGRANGE MULTIPLIER FOR FLUX
C BETA         REAL  O  LAGRANGE MULTIPLIER FOR CHISQ
C-------------------------------------------------------------------------
C
      REAL	TOL, TFLUX, TCHISQ, FLUX, CHISQ, ALPHA, BETA, 
     1		LENGTH
      DOUBLE PRECISION	GDG(4,4)
C
      DOUBLE PRECISION	DALPHA, DAMIN, DAMAX, DBETA, DBMAX, DBMIN, 
     $   DCHISQ, DFLUX, A, B, DET
C
      INTEGER H, C, F, J
      PARAMETER (H=1)
      PARAMETER (C=2)
      PARAMETER (F=3)
      PARAMETER (J=4)
C======================================================================
      A = GDG(C,J)/GDG(C,C)
      B = A**2 - (GDG(J, J) - TOL*LENGTH)/ GDG(C, C)
      IF (B.GT.0.0) THEN
         B = SQRT(B)
         DAMAX = (A + B)
         DAMIN = (A - B)
      ELSE
         DAMAX = 0.0
         DAMIN = 0.0
      END IF
C
      IF (TFLUX.LE.0.0) THEN
         DCHISQ = CHISQ - TCHISQ + GDG (C, J)
         DALPHA = DCHISQ/GDG(C, C)
C
         DALPHA = MAX(DAMIN, MIN(DAMAX, DALPHA))
         ALPHA = MAX (0.0, ALPHA + SNGL(DALPHA))
      ELSE
         A = GDG(F,J)/GDG(F,F)
         B = A**2 - (GDG(J, J) - TOL*LENGTH)/GDG(F,F)
         IF (B.GT.0.0) THEN
            B = SQRT(B)
            DBMAX = (A + B)
            DBMIN = (A - B)
         ELSE
            DBMAX = 0.0
            DBMIN = 0.0
         END IF
C
         DCHISQ = CHISQ - TCHISQ + GDG(C, J)
         DFLUX  = FLUX -  TFLUX +  GDG(F, J)
         DET = GDG(C, C)*GDG(F, F) - GDG(F, C)**2
         DALPHA = (GDG(F, F)*DCHISQ - GDG(C, F)*DFLUX)/DET
         DBETA =  (GDG(C, C)*DFLUX  - GDG(C, F)*DCHISQ)/DET
C
         DALPHA = MAX(DAMIN, MIN(DAMAX, DALPHA))
         ALPHA = MAX (0.0, ALPHA + SNGL(DALPHA))
         DBETA = MAX(DBMIN, MIN(DBMAX, DBETA))
         BETA = BETA + DBETA
C
      END IF
C
      END


