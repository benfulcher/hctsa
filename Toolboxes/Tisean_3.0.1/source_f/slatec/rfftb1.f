*DECK RFFTB1
      SUBROUTINE RFFTB1 (N, C, CH, WA, IFAC)
C***BEGIN PROLOGUE  RFFTB1
C***PURPOSE  Compute the backward fast Fourier transform of a real
C            coefficient array.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTB1-S, CFFTB1-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   Subroutine RFFTB1 computes the real periodic sequence from its
C   Fourier coefficients (Fourier synthesis).  The transform is defined
C   below at output parameter C.
C
C   The arrays WA and IFAC which are used by subroutine RFFTB1 must be
C   initialized by calling subroutine RFFTI1.
C
C   Input Arguments
C
C   N       the length of the array R to be transformed.  The method
C           is most efficient when N is a product of small primes.
C           N may change so long as different work arrays are provided.
C
C   C       a real array of length N which contains the sequence
C           to be transformed.
C
C   CH      a real work array of length at least N.
C
C   WA      a real work array which must be dimensioned at least N.
C
C   IFAC    an integer work array which must be dimensioned at least 15.
C
C           The WA and IFAC arrays must be initialized by calling
C           subroutine RFFTI1, and different WA and IFAC arrays must be
C           used for each different value of N.  This initialization
C           does not have to be repeated so long as N remains unchanged.
C           Thus subsequent transforms can be obtained faster than the
C           first.  The same WA and IFAC arrays can be used by RFFTF1
C           and RFFTB1.
C
C   Output Argument
C
C   C       For N even and for I = 1,...,N
C
C                C(I) = C(1)+(-1)**(I-1)*C(N)
C
C                     plus the sum from K=2 to K=N/2 of
C
C                      2.*C(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                     -2.*C(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C           For N odd and for I = 1,...,N
C
C                C(I) = C(1) plus the sum from K=2 to K=(N+1)/2 of
C
C                     2.*C(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
C
C                    -2.*C(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
C
C   Notes:  This transform is unnormalized since a call of RFFTF1
C           followed by a call of RFFTB1 will multiply the input
C           sequence by N.
C
C           WA and IFAC contain initialization calculations which must
C           not be destroyed between calls of subroutine RFFTF1 or
C           RFFTB1.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RADB2, RADB3, RADB4, RADB5, RADBG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from subsidiary to user-callable.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTB1
      DIMENSION CH(*), C(*), WA(*), IFAC(*)
C***FIRST EXECUTABLE STATEMENT  RFFTB1
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL RADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END
