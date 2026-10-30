*DECK CHKDER
      SUBROUTINE CHKDER (M, N, X, FVEC, FJAC, LDFJAC, XP, FVECP, MODE,
     +   ERR)
C***BEGIN PROLOGUE  CHKDER
C***PURPOSE  Check the gradients of M nonlinear functions in N
C            variables, evaluated at a point X, for consistency
C            with the functions themselves.
C***LIBRARY   SLATEC
C***CATEGORY  F3, G4C
C***TYPE      SINGLE PRECISION (CHKDER-S, DCKDER-D)
C***KEYWORDS  GRADIENTS, JACOBIAN, MINPACK, NONLINEAR
C***AUTHOR  Hiebert, K. L. (SNLA)
C***DESCRIPTION
C
C   This subroutine is a companion routine to SNLS1,SNLS1E,SNSQ,and
C   SNSQE which may be used to check the calculation of the Jacobian.
C
C     SUBROUTINE CHKDER
C
C     This subroutine checks the gradients of M nonlinear functions
C     in N variables, evaluated at a point X, for consistency with
C     the functions themselves. The user must call CKDER twice,
C     first with MODE = 1 and then with MODE = 2.
C
C     MODE = 1. On input, X must contain the point of evaluation.
C               On output, XP is set to a neighboring point.
C
C     MODE = 2. On input, FVEC must contain the functions and the
C                         rows of FJAC must contain the gradients
C                         of the respective functions each evaluated
C                         at X, and FVECP must contain the functions
C                         evaluated at XP.
C               On output, ERR contains measures of correctness of
C                          the respective gradients.
C
C     The subroutine does not perform reliably if cancellation or
C     rounding errors cause a severe loss of significance in the
C     evaluation of a function. Therefore, none of the components
C     of X should be unusually small (in particular, zero) or any
C     other value which may cause loss of significance.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE CHKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of functions.
C
C       N is a positive integer input variable set to the number
C         of variables.
C
C       X is an input array of length N.
C
C       FVEC is an array of length M. On input when MODE = 2,
C         FVEC must contain the functions evaluated at X.
C
C       FJAC is an M by N array. On input when MODE = 2,
C         the rows of FJAC must contain the gradients of
C         the respective functions evaluated at X.
C
C       LDFJAC is a positive integer input parameter not less than M
C         which specifies the leading dimension of the array FJAC.
C
C       XP is an array of length N. On output when MODE = 1,
C         XP is set to a neighboring point of X.
C
C       FVECP is an array of length M. On input when MODE = 2,
C         FVECP must contain the functions evaluated at XP.
C
C       MODE is an integer input variable set to 1 on the first call
C         and 2 on the second. Other values of MODE are equivalent
C         to MODE = 1.
C
C       ERR is an array of length M. On output when MODE = 2,
C         ERR contains measures of correctness of the respective
C         gradients. If there is no severe loss of significance,
C         then if ERR(I) is 1.0 the I-th gradient is correct,
C         while if ERR(I) is 0.0 the I-th gradient is incorrect.
C         For values of ERR between 0.0 and 1.0, the categorization
C         is less certain. In general, a value of ERR(I) greater
C         than 0.5 indicates that the I-th gradient is probably
C         correct, while a value of ERR(I) less than 0.5 indicates
C         that the I-th gradient is probably incorrect.
C
C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
C                 tions. In Numerical Methods for Nonlinear Algebraic
C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
C                 1988.
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  CHKDER
      INTEGER M,N,LDFJAC,MODE
      REAL X(*),FVEC(*),FJAC(LDFJAC,*),XP(*),FVECP(*),ERR(*)
      INTEGER I,J
      REAL EPS,EPSF,EPSLOG,EPSMCH,FACTOR,ONE,TEMP,ZERO
      REAL R1MACH
      SAVE FACTOR, ONE, ZERO
C
      DATA FACTOR,ONE,ZERO /1.0E2,1.0E0,0.0E0/
C***FIRST EXECUTABLE STATEMENT  CHKDER
      EPSMCH = R1MACH(4)
C
      EPS = SQRT(EPSMCH)
C
      IF (MODE .EQ. 2) GO TO 20
C
C        MODE = 1.
C
         DO 10 J = 1, N
            TEMP = EPS*ABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = EPS
            XP(J) = X(J) + TEMP
   10       CONTINUE
         GO TO 70
   20 CONTINUE
C
C        MODE = 2.
C
         EPSF = FACTOR*EPSMCH
         EPSLOG = LOG10(EPS)
         DO 30 I = 1, M
            ERR(I) = ZERO
   30       CONTINUE
         DO 50 J = 1, N
            TEMP = ABS(X(J))
            IF (TEMP .EQ. ZERO) TEMP = ONE
            DO 40 I = 1, M
               ERR(I) = ERR(I) + TEMP*FJAC(I,J)
   40          CONTINUE
   50       CONTINUE
         DO 60 I = 1, M
            TEMP = ONE
            IF (FVEC(I) .NE. ZERO .AND. FVECP(I) .NE. ZERO
     1          .AND. ABS(FVECP(I)-FVEC(I)) .GE. EPSF*ABS(FVEC(I)))
     2         TEMP = EPS*ABS((FVECP(I)-FVEC(I))/EPS-ERR(I))
     3                /(ABS(FVEC(I)) + ABS(FVECP(I)))
            ERR(I) = ONE
            IF (TEMP .GT. EPSMCH .AND. TEMP .LT. EPS)
     1         ERR(I) = (LOG10(TEMP) - EPSLOG)/EPSLOG
            IF (TEMP .GE. EPS) ERR(I) = ZERO
   60       CONTINUE
   70 CONTINUE
C
      RETURN
C
C     LAST CARD OF SUBROUTINE CHKDER.
C
      END
