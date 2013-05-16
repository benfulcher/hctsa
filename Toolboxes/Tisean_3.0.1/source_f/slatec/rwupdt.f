*DECK RWUPDT
      SUBROUTINE RWUPDT (N, R, LDR, W, B, ALPHA, COS, SIN)
C***BEGIN PROLOGUE  RWUPDT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SNLS1 and SNLS1E
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (RWUPDT-S, DWUPDT-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an N by N upper triangular matrix R, this subroutine
C     computes the QR decomposition of the matrix formed when a row
C     is added to R. If the row is specified by the vector W, then
C     RWUPDT determines an orthogonal matrix Q such that when the
C     N+1 by N matrix composed of R augmented by W is premultiplied
C     by (Q TRANSPOSE), the resulting matrix is upper trapezoidal.
C     The orthogonal matrix Q is the product of N transformations
C
C           G(1)*G(2)* ... *G(N)
C
C     where G(I) is a Givens rotation in the (I,N+1) plane which
C     eliminates elements in the I-th plane. RWUPDT also
C     computes the product (Q TRANSPOSE)*C where C is the
C     (N+1)-vector (b,alpha). Q itself is not accumulated, rather
C     the information to recover the G rotations is supplied.
C
C     The subroutine statement is
C
C       SUBROUTINE RWUPDT(N,R,LDR,W,B,ALPHA,COS,SIN)
C
C     where
C
C       N is a positive integer input variable set to the order of R.
C
C       R is an N by N array. On input the upper triangular part of
C         R must contain the matrix to be updated. On output R
C         contains the updated triangular matrix.
C
C       LDR is a positive integer input variable not less than N
C         which specifies the leading dimension of the array R.
C
C       W is an input array of length N which must contain the row
C         vector to be added to R.
C
C       B is an array of length N. On input B must contain the
C         first N elements of the vector C. On output B contains
C         the first N elements of the vector (Q TRANSPOSE)*C.
C
C       ALPHA is a variable. On input ALPHA must contain the
C         (N+1)-st element of the vector C. On output ALPHA contains
C         the (N+1)-st element of the vector (Q TRANSPOSE)*C.
C
C       COS is an output array of length N which contains the
C         cosines of the transforming Givens rotations.
C
C       SIN is an output array of length N which contains the
C         sines of the transforming Givens rotations.
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
C***END PROLOGUE  RWUPDT
      INTEGER N,LDR
      REAL ALPHA
      REAL R(LDR,*),W(*),B(*),COS(*),SIN(*)
      INTEGER I,J,JM1
      REAL COTAN,ONE,P5,P25,ROWJ,TAN,TEMP,ZERO
      SAVE ONE, P5, P25, ZERO
      DATA ONE,P5,P25,ZERO /1.0E0,5.0E-1,2.5E-1,0.0E0/
C***FIRST EXECUTABLE STATEMENT  RWUPDT
      DO 60 J = 1, N
         ROWJ = W(J)
         JM1 = J - 1
C
C        APPLY THE PREVIOUS TRANSFORMATIONS TO
C        R(I,J), I=1,2,...,J-1, AND TO W(J).
C
         IF (JM1 .LT. 1) GO TO 20
         DO 10 I = 1, JM1
            TEMP = COS(I)*R(I,J) + SIN(I)*ROWJ
            ROWJ = -SIN(I)*R(I,J) + COS(I)*ROWJ
            R(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES W(J).
C
         COS(J) = ONE
         SIN(J) = ZERO
         IF (ROWJ .EQ. ZERO) GO TO 50
         IF (ABS(R(J,J)) .GE. ABS(ROWJ)) GO TO 30
            COTAN = R(J,J)/ROWJ
            SIN(J) = P5/SQRT(P25+P25*COTAN**2)
            COS(J) = SIN(J)*COTAN
            GO TO 40
   30    CONTINUE
            TAN = ROWJ/R(J,J)
            COS(J) = P5/SQRT(P25+P25*TAN**2)
            SIN(J) = COS(J)*TAN
   40    CONTINUE
C
C        APPLY THE CURRENT TRANSFORMATION TO R(J,J), B(J), AND ALPHA.
C
         R(J,J) = COS(J)*R(J,J) + SIN(J)*ROWJ
         TEMP = COS(J)*B(J) + SIN(J)*ALPHA
         ALPHA = -SIN(J)*B(J) + COS(J)*ALPHA
         B(J) = TEMP
   50    CONTINUE
   60    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE RWUPDT.
C
      END
