*DECK LMPAR
      SUBROUTINE LMPAR (N, R, LDR, IPVT, DIAG, QTB, DELTA, PAR, X,
     +   SIGMA, WA1, WA2)
C***BEGIN PROLOGUE  LMPAR
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SNLS1 and SNLS1E
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (LMPAR-S, DMPAR-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N matrix A, an N by N nonsingular DIAGONAL
C     matrix D, an M-vector B, and a positive number DELTA,
C     the problem is to determine a value for the parameter
C     PAR such that if X solves the system
C
C           A*X = B ,     SQRT(PAR)*D*X = 0 ,
C
C     in the least squares sense, and DXNORM is the Euclidean
C     norm of D*X, then either PAR is zero and
C
C           (DXNORM-DELTA) .LE. 0.1*DELTA ,
C
C     or PAR is positive and
C
C           ABS(DXNORM-DELTA) .LE. 0.1*DELTA .
C
C     This subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     QR factorization, with column pivoting, of A. That is, if
C     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
C     columns, and R is an upper triangular matrix with diagonal
C     elements of nonincreasing magnitude, then LMPAR expects
C     the full upper triangle of R, the permutation matrix P,
C     and the first N components of (Q TRANSPOSE)*B. On output
C     LMPAR also provides an upper triangular matrix S such that
C
C            T   T                   T
C           P *(A *A + PAR*D*D)*P = S *S .
C
C     S is employed within LMPAR and may be of separate interest.
C
C     Only a few iterations are generally needed for convergence
C     of the algorithm. If, however, the limit of 10 iterations
C     is reached, then the output PAR will contain the best
C     value obtained so far.
C
C     The subroutine statement is
C
C       SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA,
C                        WA1,WA2)
C
C     where
C
C       N is a positive integer input variable set to the order of R.
C
C       R is an N by N array. On input the full upper triangle
C         must contain the full upper triangle of the matrix R.
C         On output the full upper triangle is unaltered, and the
C         strict lower triangle contains the strict upper triangle
C         (transposed) of the upper triangular matrix S.
C
C       LDR is a positive integer input variable not less than N
C         which specifies the leading dimension of the array R.
C
C       IPVT is an integer input array of length N which defines the
C         permutation matrix P such that A*P = Q*R. Column J of P
C         is column IPVT(J) of the identity matrix.
C
C       DIAG is an input array of length N which must contain the
C         diagonal elements of the matrix D.
C
C       QTB is an input array of length N which must contain the first
C         N elements of the vector (Q TRANSPOSE)*B.
C
C       DELTA is a positive input variable which specifies an upper
C         bound on the Euclidean norm of D*X.
C
C       PAR is a nonnegative variable. On input PAR contains an
C         initial estimate of the Levenberg-Marquardt parameter.
C         On output PAR contains the final estimate.
C
C       X is an output array of length N which contains the least
C         squares solution of the system A*X = B, SQRT(PAR)*D*X = 0,
C         for the output PAR.
C
C       SIGMA is an output array of length N which contains the
C         diagonal elements of the upper triangular matrix S.
C
C       WA1 and WA2 are work arrays of length N.
C
C***SEE ALSO  SNLS1, SNLS1E
C***ROUTINES CALLED  ENORM, QRSOLV, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  LMPAR
      INTEGER N,LDR
      INTEGER IPVT(*)
      REAL DELTA,PAR
      REAL R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA1(*),WA2(*)
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      REAL DXNORM,DWARF,FP,GNORM,PARC,PARL,PARU,P1,P001,SUM,TEMP,ZERO
      REAL R1MACH,ENORM
      SAVE P1, P001, ZERO
      DATA P1,P001,ZERO /1.0E-1,1.0E-3,0.0E0/
C***FIRST EXECUTABLE STATEMENT  LMPAR
      DWARF = R1MACH(1)
C
C     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
C     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
C
      NSING = N
      DO 10 J = 1, N
         WA1(J) = QTB(J)
         IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA1(J) = ZERO
   10    CONTINUE
      IF (NSING .LT. 1) GO TO 50
      DO 40 K = 1, NSING
         J = NSING - K + 1
         WA1(J) = WA1(J)/R(J,J)
         TEMP = WA1(J)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 30
         DO 20 I = 1, JM1
            WA1(I) = WA1(I) - R(I,J)*TEMP
   20       CONTINUE
   30    CONTINUE
   40    CONTINUE
   50 CONTINUE
      DO 60 J = 1, N
         L = IPVT(J)
         X(L) = WA1(J)
   60    CONTINUE
C
C     INITIALIZE THE ITERATION COUNTER.
C     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
C     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
C
      ITER = 0
      DO 70 J = 1, N
         WA2(J) = DIAG(J)*X(J)
   70    CONTINUE
      DXNORM = ENORM(N,WA2)
      FP = DXNORM - DELTA
      IF (FP .LE. P1*DELTA) GO TO 220
C
C     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
C     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
C     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
C
      PARL = ZERO
      IF (NSING .LT. N) GO TO 120
      DO 80 J = 1, N
         L = IPVT(J)
         WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
   80    CONTINUE
      DO 110 J = 1, N
         SUM = ZERO
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 100
         DO 90 I = 1, JM1
            SUM = SUM + R(I,J)*WA1(I)
   90       CONTINUE
  100    CONTINUE
         WA1(J) = (WA1(J) - SUM)/R(J,J)
  110    CONTINUE
      TEMP = ENORM(N,WA1)
      PARL = ((FP/DELTA)/TEMP)/TEMP
  120 CONTINUE
C
C     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
C
      DO 140 J = 1, N
         SUM = ZERO
         DO 130 I = 1, J
            SUM = SUM + R(I,J)*QTB(I)
  130       CONTINUE
         L = IPVT(J)
         WA1(J) = SUM/DIAG(L)
  140    CONTINUE
      GNORM = ENORM(N,WA1)
      PARU = GNORM/DELTA
      IF (PARU .EQ. ZERO) PARU = DWARF/MIN(DELTA,P1)
C
C     IF THE INPUT PAR LIES OUTSIDE OF THE INTERVAL (PARL,PARU),
C     SET PAR TO THE CLOSER ENDPOINT.
C
      PAR = MAX(PAR,PARL)
      PAR = MIN(PAR,PARU)
      IF (PAR .EQ. ZERO) PAR = GNORM/DXNORM
C
C     BEGINNING OF AN ITERATION.
C
  150 CONTINUE
         ITER = ITER + 1
C
C        EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
C
         IF (PAR .EQ. ZERO) PAR = MAX(DWARF,P001*PARU)
         TEMP = SQRT(PAR)
         DO 160 J = 1, N
            WA1(J) = TEMP*DIAG(J)
  160       CONTINUE
         CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SIGMA,WA2)
         DO 170 J = 1, N
            WA2(J) = DIAG(J)*X(J)
  170       CONTINUE
         DXNORM = ENORM(N,WA2)
         TEMP = FP
         FP = DXNORM - DELTA
C
C        IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
C        OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
C        IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
C
         IF (ABS(FP) .LE. P1*DELTA
     1       .OR. PARL .EQ. ZERO .AND. FP .LE. TEMP
     2            .AND. TEMP .LT. ZERO .OR. ITER .EQ. 10) GO TO 220
C
C        COMPUTE THE NEWTON CORRECTION.
C
         DO 180 J = 1, N
            L = IPVT(J)
            WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
  180       CONTINUE
         DO 210 J = 1, N
            WA1(J) = WA1(J)/SIGMA(J)
            TEMP = WA1(J)
            JP1 = J + 1
            IF (N .LT. JP1) GO TO 200
            DO 190 I = JP1, N
               WA1(I) = WA1(I) - R(I,J)*TEMP
  190          CONTINUE
  200       CONTINUE
  210       CONTINUE
         TEMP = ENORM(N,WA1)
         PARC = ((FP/DELTA)/TEMP)/TEMP
C
C        DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
C
         IF (FP .GT. ZERO) PARL = MAX(PARL,PAR)
         IF (FP .LT. ZERO) PARU = MIN(PARU,PAR)
C
C        COMPUTE AN IMPROVED ESTIMATE FOR PAR.
C
         PAR = MAX(PARL,PAR+PARC)
C
C        END OF AN ITERATION.
C
         GO TO 150
  220 CONTINUE
C
C     TERMINATION.
C
      IF (ITER .EQ. 0) PAR = ZERO
      RETURN
C
C     LAST CARD OF SUBROUTINE LMPAR.
C
      END
