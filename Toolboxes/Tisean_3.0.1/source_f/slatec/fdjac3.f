*DECK FDJAC3
      SUBROUTINE FDJAC3 (FCN, M, N, X, FVEC, FJAC, LDFJAC, IFLAG,
     +   EPSFCN, WA)
C***BEGIN PROLOGUE  FDJAC3
C***SUBSIDIARY
C***PURPOSE  Subsidiary to SNLS1 and SNLS1E
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (FDJAC3-S, DFDJC3-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine computes a forward-difference approximation
C     to the M by N Jacobian matrix associated with a specified
C     problem of M functions in N variables.
C
C     The subroutine statement is
C
C       SUBROUTINE FDJAC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
C
C     where
C
C       FCN is the name of the user-supplied subroutine which
C         calculates the functions. FCN must be declared
C         in an external statement in the user calling
C         program, and should be written as follows.
C
C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C         INTEGER LDFJAC,M,N,IFLAG
C         REAL X(N),FVEC(M),FJAC(LDFJAC,N)
C         ----------
C         When IFLAG.EQ.1 calculate the functions at X and
C         return this vector in FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless
C         the user wants to terminate execution of FDJAC3.
C         In this case set IFLAG to a negative integer.
C
C       M is a positive integer input variable set to the number
C         of functions.
C
C       N is a positive integer input variable set to the number
C         of variables. N must not exceed M.
C
C       X is an input array of length N.
C
C       FVEC is an input array of length M which must contain the
C         functions evaluated at X.
C
C       FJAC is an output M by N array which contains the
C         approximation to the Jacobian matrix evaluated at X.
C
C       LDFJAC is a positive integer input variable not less than M
C         which specifies the leading dimension of the array FJAC.
C
C       IFLAG is an integer variable which can be used to terminate
C         THE EXECUTION OF FDJAC3. See description of FCN.
C
C       EPSFCN is an input variable used in determining a suitable
C         step length for the forward-difference approximation. This
C         approximation assumes that the relative errors in the
C         functions are of the order of EPSFCN. If EPSFCN is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
C
C       WA is a work array of length M.
C
C***SEE ALSO  SNLS1, SNLS1E
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  FDJAC3
      INTEGER M,N,LDFJAC,IFLAG
      REAL EPSFCN
      REAL X(*),FVEC(*),FJAC(LDFJAC,*),WA(*)
      INTEGER I,J
      REAL EPS,EPSMCH,H,TEMP,ZERO
      REAL R1MACH
      SAVE ZERO
      DATA ZERO /0.0E0/
C***FIRST EXECUTABLE STATEMENT  FDJAC3
      EPSMCH = R1MACH(4)
C
      EPS = SQRT(MAX(EPSFCN,EPSMCH))
C      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES
C           ARE TO BE RETURNED BY FCN.
      IFLAG = 1
      DO 20 J = 1, N
         TEMP = X(J)
         H = EPS*ABS(TEMP)
         IF (H .EQ. ZERO) H = EPS
         X(J) = TEMP + H
         CALL FCN(IFLAG,M,N,X,WA,FJAC,LDFJAC)
         IF (IFLAG .LT. 0) GO TO 30
         X(J) = TEMP
         DO 10 I = 1, M
            FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE FDJAC3.
C
      END
