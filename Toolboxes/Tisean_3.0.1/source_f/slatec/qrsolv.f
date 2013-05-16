*DECK QRSOLV
      SUBROUTINE QRSOLV (N, R, LDR, IPVT, DIAG, QTB, X, SIGMA, WA)
C***BEGIN PROLOGUE  QRSOLV
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SNLS1 and SNLS1E
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QRSOLV-S, DQRSLV-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N matrix A, an N by N diagonal matrix D,
C     and an M-vector B, the problem is to determine an X which
C     solves the system
C
C           A*X = B ,     D*X = 0 ,
C
C     in the least squares sense.
C
C     This subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     QR factorization, with column pivoting, of A. That is, if
C     A*P = Q*R, where P is a permutation matrix, Q has orthogonal
C     columns, and R is an upper triangular matrix with diagonal
C     elements of nonincreasing magnitude, then QRSOLV expects
C     the full upper triangle of R, the permutation matrix P,
C     and the first N components of (Q TRANSPOSE)*B. The system
C     A*X = B, D*X = 0, is then equivalent to
C
C                  T       T
C           R*Z = Q *B ,  P *D*P*Z = 0 ,
C
C     where X = P*Z. If this system does not have full rank,
C     then a least squares solution is obtained. On output QRSOLV
C     also provides an upper triangular matrix S such that
C
C            T   T               T
C           P *(A *A + D*D)*P = S *S .
C
C     S is computed within QRSOLV and may be of separate interest.
C
C     The subroutine statement is
C
C       SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA)
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
C       X is an output array of length N which contains the least
C         squares solution of the system A*X = B, D*X = 0.
C
C       SIGMA is an output array of length N which contains the
C         diagonal elements of the upper triangular matrix S.
C
C       WA is a work array of length N.
C
C***SEE ALSO  SNLS1, SNLS1E
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  QRSOLV
      INTEGER N,LDR
      INTEGER IPVT(*)
      REAL R(LDR,*),DIAG(*),QTB(*),X(*),SIGMA(*),WA(*)
      INTEGER I,J,JP1,K,KP1,L,NSING
      REAL COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
      SAVE P5, P25, ZERO
      DATA P5,P25,ZERO /5.0E-1,2.5E-1,0.0E0/
C***FIRST EXECUTABLE STATEMENT  QRSOLV
      DO 20 J = 1, N
         DO 10 I = J, N
            R(I,J) = R(J,I)
   10       CONTINUE
         X(J) = R(J,J)
         WA(J) = QTB(J)
   20    CONTINUE
C
C     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
C
      DO 100 J = 1, N
C
C        PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
C        DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
C
         L = IPVT(J)
         IF (DIAG(L) .EQ. ZERO) GO TO 90
         DO 30 K = J, N
            SIGMA(K) = ZERO
   30       CONTINUE
         SIGMA(J) = DIAG(L)
C
C        THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
C        MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
C        BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
C
         QTBPJ = ZERO
         DO 80 K = J, N
C
C           DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C           APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
C
            IF (SIGMA(K) .EQ. ZERO) GO TO 70
            IF (ABS(R(K,K)) .GE. ABS(SIGMA(K))) GO TO 40
               COTAN = R(K,K)/SIGMA(K)
               SIN = P5/SQRT(P25+P25*COTAN**2)
               COS = SIN*COTAN
               GO TO 50
   40       CONTINUE
               TAN = SIGMA(K)/R(K,K)
               COS = P5/SQRT(P25+P25*TAN**2)
               SIN = COS*TAN
   50       CONTINUE
C
C           COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
C           THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
C
            R(K,K) = COS*R(K,K) + SIN*SIGMA(K)
            TEMP = COS*WA(K) + SIN*QTBPJ
            QTBPJ = -SIN*WA(K) + COS*QTBPJ
            WA(K) = TEMP
C
C           ACCUMULATE THE TRANSFORMATION IN THE ROW OF S.
C
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 70
            DO 60 I = KP1, N
               TEMP = COS*R(I,K) + SIN*SIGMA(I)
               SIGMA(I) = -SIN*R(I,K) + COS*SIGMA(I)
               R(I,K) = TEMP
   60          CONTINUE
   70       CONTINUE
   80       CONTINUE
   90    CONTINUE
C
C        STORE THE DIAGONAL ELEMENT OF S AND RESTORE
C        THE CORRESPONDING DIAGONAL ELEMENT OF R.
C
         SIGMA(J) = R(J,J)
         R(J,J) = X(J)
  100    CONTINUE
C
C     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
C     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
C
      NSING = N
      DO 110 J = 1, N
         IF (SIGMA(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
         IF (NSING .LT. N) WA(J) = ZERO
  110    CONTINUE
      IF (NSING .LT. 1) GO TO 150
      DO 140 K = 1, NSING
         J = NSING - K + 1
         SUM = ZERO
         JP1 = J + 1
         IF (NSING .LT. JP1) GO TO 130
         DO 120 I = JP1, NSING
            SUM = SUM + R(I,J)*WA(I)
  120       CONTINUE
  130    CONTINUE
         WA(J) = (WA(J) - SUM)/SIGMA(J)
  140    CONTINUE
  150 CONTINUE
C
C     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
C
      DO 160 J = 1, N
         L = IPVT(J)
         X(L) = WA(J)
  160    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE QRSOLV.
C
      END
