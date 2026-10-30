*DECK RFFTF1
      SUBROUTINE RFFTF1 (N, C, CH, WA, IFAC)
C***BEGIN PROLOGUE  RFFTF1
C***PURPOSE  Compute the forward transform of a real, periodic sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTF1-S, CFFTF1-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   Subroutine RFFTF1 computes the Fourier coefficients of a real
C   periodic sequence (Fourier analysis).  The transform is defined
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
C   C       C(1) = the sum from I=1 to I=N of R(I)
C
C           If N is even set L = N/2; if N is odd set L = (N+1)/2
C
C             then for K = 2,...,L
C
C                C(2*K-2) = the sum from I = 1 to I = N of
C
C                     C(I)*COS((K-1)*(I-1)*2*PI/N)
C
C                C(2*K-1) = the sum from I = 1 to I = N of
C
C                    -C(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C           If N is even
C
C                C(N) = the sum from I = 1 to I = N of
C
C                     (-1)**(I-1)*C(I)
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
C***ROUTINES CALLED  RADF2, RADF3, RADF4, RADF5, RADFG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from subsidiary to user-callable.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTF1
      DIMENSION CH(*), C(*), WA(*), IFAC(*)
C***FIRST EXECUTABLE STATEMENT  RFFTF1
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL RADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END
