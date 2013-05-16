*DECK RFFTI1
      SUBROUTINE RFFTI1 (N, WA, IFAC)
C***BEGIN PROLOGUE  RFFTI1
C***PURPOSE  Initialize a real and an integer work array for RFFTF1 and
C            RFFTB1.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTI1-S, CFFTI1-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   Subroutine RFFTI1 initializes the work arrays WA and IFAC which are
C   used in both RFFTF1 and RFFTB1.  The prime factorization of N and a
C   tabulation of the trigonometric functions are computed and stored in
C   IFAC and WA, respectively.
C
C   Input Argument
C
C   N       the length of the sequence to be transformed.
C
C   Output Arguments
C
C   WA      a real work array which must be dimensioned at least N.
C
C   IFAC    an integer work array which must be dimensioned at least 15.
C
C   The same work arrays can be used for both RFFTF1 and RFFTB1 as long
C   as N remains unchanged.  Different WA and IFAC arrays are required
C   for different values of N.  The contents of WA and IFAC must not be
C   changed between calls of RFFTF1 or RFFTB1.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL, and
C           (c) changing definition of variable TPI by using
C               FORTRAN intrinsic functions instead of DATA
C               statements.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from subsidiary to user-callable.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTI1
      DIMENSION WA(*), IFAC(*), NTRYH(4)
      SAVE NTRYH
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C***FIRST EXECUTABLE STATEMENT  RFFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.*ATAN(1.)
      ARGH = TPI/N
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = LD*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
