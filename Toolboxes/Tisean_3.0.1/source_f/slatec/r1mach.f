      FUNCTION r1MACH (I)
c this is not the original one from slatec 
      dimension const(5)
c small:
      DATA const(1) / 1.18E-38      /
c large:
      DATA const(2) / 3.40E+38      /
c diff:
      DATA const(3) / 0.595E-07     /
      DATA const(4) / 1.19E-07      /
c log10:
      DATA const(5) / 0.30102999566 /

C***FIRST EXECUTABLE STATEMENT  R1MACH
C
      R1MACH = const(I)
      RETURN
C
      END


