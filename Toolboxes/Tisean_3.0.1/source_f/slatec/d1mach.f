      DOUBLE PRECISION FUNCTION D1MACH (I)
c this is not the original one from slatec 
      double precision const(5)
c small:
      DATA const(1) / 2.23D-308  /
c large:
      DATA const(2) / 1.79D+308  /
c diff:
      DATA const(3) / 1.11D-16   /
      DATA const(4) / 2.22D-16   /
c log10:
      DATA const(5) / 0.301029995663981195D0 / 

C***FIRST EXECUTABLE STATEMENT  D1MACH
C
      D1MACH = const(I)
      RETURN
C
      END
