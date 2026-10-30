*DECK RGAUSS
      FUNCTION RGAUSS (XMEAN, SD)
C***BEGIN PROLOGUE  RGAUSS
C***PURPOSE  Generate a normally distributed (Gaussian) random number.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  L6A14
C***TYPE      SINGLE PRECISION (RGAUSS-S)
C***KEYWORDS  FNLIB, GAUSSIAN, NORMAL, RANDOM NUMBER, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Generate a normally distributed random number, i.e., generate random
C numbers with a Gaussian distribution.  These random numbers are not
C exceptionally good -- especially in the tails of the distribution,
C but this implementation is simple and suitable for most applications.
C See R. W. Hamming, Numerical Methods for Scientists and Engineers,
C McGraw-Hill, 1962, pages 34 and 389.
C
C             Input Arguments --
C XMEAN  the mean of the Guassian distribution.
C SD     the standard deviation of the Guassian function
C          EXP (-1/2 * (X-XMEAN)**2 / SD**2)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  RAND
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   910819  Added EXTERNAL statement for RAND due to problem on IBM
C           RS 6000.  (WRB)
C***END PROLOGUE  RGAUSS
      EXTERNAL RAND
C***FIRST EXECUTABLE STATEMENT  RGAUSS
      RGAUSS = -6.0
      DO 10 I=1,12
        RGAUSS = RGAUSS + RAND(0.0)
 10   CONTINUE
C
      RGAUSS = XMEAN + SD*RGAUSS
C
      RETURN
      END
