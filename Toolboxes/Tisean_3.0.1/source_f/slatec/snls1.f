*DECK SNLS1
      SUBROUTINE SNLS1 (FCN, IOPT, M, N, X, FVEC, FJAC, LDFJAC, FTOL,
     +   XTOL, GTOL, MAXFEV, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO,
     +   NFEV, NJEV, IPVT, QTF, WA1, WA2, WA3, WA4)
C***BEGIN PROLOGUE  SNLS1
C***PURPOSE  Minimize the sum of the squares of M nonlinear functions
C            in N variables by a modification of the Levenberg-Marquardt
C            algorithm.
C***LIBRARY   SLATEC
C***CATEGORY  K1B1A1, K1B1A2
C***TYPE      SINGLE PRECISION (SNLS1-S, DNLS1-D)
C***KEYWORDS  LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
C             NONLINEAR LEAST SQUARES
C***AUTHOR  Hiebert, K. L., (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of SNLS1 is to minimize the sum of the squares of M
C       nonlinear functions in N variables by a modification of the
C       Levenberg-Marquardt algorithm.  The user must provide a subrou-
C       tine which calculates the functions.  The user has the option
C       of how the Jacobian will be supplied.  The user can supply the
C       full Jacobian, or the rows of the Jacobian (to avoid storing
C       the full Jacobian), or let the code approximate the Jacobian by
C       forward-differencing.   This code is the combination of the
C       MINPACK codes (Argonne) LMDER, LMDIF, and LMSTR.
C
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE SNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
C      *                 GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO
C      *                 ,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
C       INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
C       INTEGER IPVT(N)
C       REAL FTOL,XTOL,GTOL,EPSFCN,FACTOR
C       REAL X(N),FVEC(M),FJAC(LDFJAC,N),DIAG(N),QTF(N),
C      *     WA1(N),WA2(N),WA3(N),WA4(M)
C
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to SNLS1 and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from SNLS1.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  If the user wants to supply the Jacobian
C         (IOPT=2 or 3), then FCN must be written to calculate the
C         Jacobian, as well as the functions.  See the explanation
C         of the IOPT argument below.
C         If the user wants the iterates printed (NPRINT positive), then
C         FCN must do the printing.  See the explanation of NPRINT
C         below.  FCN must be declared in an EXTERNAL statement in the
C         calling program and should be written as follows.
C
C
C         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C         INTEGER IFLAG,LDFJAC,M,N
C         REAL X(N),FVEC(M)
C         ----------
C         FJAC and LDFJAC may be ignored     , if IOPT=1.
C         REAL FJAC(LDFJAC,N)                , if IOPT=2.
C         REAL FJAC(N)                       , if IOPT=3.
C         ----------
C           If IFLAG=0, the values in X and FVEC are available
C           for printing.  See the explanation of NPRINT below.
C           IFLAG will never be zero unless NPRINT is positive.
C           The values of X and FVEC must not be changed.
C         RETURN
C         ----------
C           If IFLAG=1, calculate the functions at X and return
C           this vector in FVEC.
C         RETURN
C         ----------
C           If IFLAG=2, calculate the full Jacobian at X and return
C           this matrix in FJAC.  Note that IFLAG will never be 2 unless
C           IOPT=2.  FVEC contains the function values at X and must
C           not be altered.  FJAC(I,J) must be set to the derivative
C           of FVEC(I) with respect to X(J).
C         RETURN
C         ----------
C           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
C           and return this vector in FJAC.  Note that IFLAG will
C           never be 3 unless IOPT=3.  FVEC contains the function
C           values at X and must not be altered.  FJAC(J) must be
C           set to the derivative of FVEC(LDFJAC) with respect to X(J).
C         RETURN
C         ----------
C         END
C
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of SNLS1.  In this case, set
C         IFLAG to a negative integer.
C
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=2 or 3, then the user must supply the
C         Jacobian, as well as the function values, through the
C         subroutine FCN.  If IOPT=2, the user supplies the full
C         Jacobian with one call to FCN.  If IOPT=3, the user supplies
C         one row of the Jacobian with each call.  (In this manner,
C         storage can be saved because the full Jacobian is not stored.)
C         If IOPT=1, the code will approximate the Jacobian by forward
C         differencing.
C
C       M is a positive integer input variable set to the number of
C         functions.
C
C       N is a positive integer input variable set to the number of
C         variables.  N must not exceed M.
C
C       X is an array of length N.  On input, X must contain an initial
C         estimate of the solution vector.  On output, X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length M which contains the functions
C         evaluated at the output X.
C
C       FJAC is an output array.  For IOPT=1 and 2, FJAC is an M by N
C         array.  For IOPT=3, FJAC is an N by N array.  The upper N by N
C         submatrix of FJAC contains an upper triangular matrix R with
C         diagonal elements of nonincreasing magnitude such that
C
C                T     T           T
C               P *(JAC *JAC)*P = R *R,
C
C         where P is a permutation matrix and JAC is the final calcu-
C         lated Jacobian.  Column J of P is column IPVT(J) (see below)
C         of the identity matrix.  The lower part of FJAC contains
C         information generated during the computation of R.
C
C       LDFJAC is a positive integer input variable which specifies
C         the leading dimension of the array FJAC.  For IOPT=1 and 2,
C         LDFJAC must not be less than M.  For IOPT=3, LDFJAC must not
C         be less than N.
C
C       FTOL is a non-negative input variable.  Termination occurs when
C         both the actual and predicted relative reductions in the sum
C         of squares are at most FTOL.  Therefore, FTOL measures the
C         relative error desired in the sum of squares.  Section 4 con-
C         tains more details about FTOL.
C
C       XTOL is a non-negative input variable.  Termination occurs when
C         the relative error between two consecutive iterates is at most
C         XTOL.  Therefore, XTOL measures the relative error desired in
C         the approximate solution.  Section 4 contains more details
C         about XTOL.
C
C       GTOL is a non-negative input variable.  Termination occurs when
C         the cosine of the angle between FVEC and any column of the
C         Jacobian is at most GTOL in absolute value.  Therefore, GTOL
C         measures the orthogonality desired between the function vector
C         and the columns of the Jacobian.  Section 4 contains more
C         details about GTOL.
C
C       MAXFEV is a positive integer input variable.  Termination occurs
C         when the number of calls to FCN to evaluate the functions
C         has reached MAXFEV.
C
C       EPSFCN is an input variable used in determining a suitable step
C         for the forward-difference approximation.  This approximation
C         assumes that the relative errors in the functions are of the
C         order of EPSFCN.  If EPSFCN is less than the machine preci-
C         sion, it is assumed that the relative errors in the functions
C         are of the order of the machine precision.  If IOPT=2 or 3,
C         then EPSFCN can be ignored (treat it as a dummy argument).
C
C       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
C         internally set.  If MODE = 2, DIAG must contain positive
C         entries that serve as implicit (multiplicative) scale factors
C         for the variables.
C
C       MODE is an integer input variable.  If MODE = 1, the variables
C         will be scaled internally.  If MODE = 2, the scaling is speci-
C         fied by the input DIAG.  Other values of MODE are equivalent
C         to MODE = 1.
C
C       FACTOR is a positive input variable used in determining the ini-
C         tial step bound.  This bound is set to the product of FACTOR
C         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
C         itself.  In most cases FACTOR should lie in the interval
C         (.1,100.).  100. is a generally recommended value.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. Appropriate
C         print statements must be added to FCN (see example) and
C         FVEC should not be altered.  If NPRINT is not positive, no
C         special calls to FCN with IFLAG = 0 are made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  improper input parameters.
C
C         INFO = 1  both actual and predicted relative reductions in the
C                   sum of squares are at most FTOL.
C
C         INFO = 2  relative error between two consecutive iterates is
C                   at most XTOL.
C
C         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
C
C         INFO = 4  the cosine of the angle between FVEC and any column
C                   of the Jacobian is at most GTOL in absolute value.
C
C         INFO = 5  number of calls to FCN for function evaluation
C                   has reached MAXFEV.
C
C         INFO = 6  FTOL is too small.  No further reduction in the sum
C                   of squares is possible.
C
C         INFO = 7  XTOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 8  GTOL is too small.  FVEC is orthogonal to the
C                   columns of the Jacobian to machine precision.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       NFEV is an integer output variable set to the number of calls to
C         FCN for function evaluation.
C
C       NJEV is an integer output variable set to the number of
C         evaluations of the full Jacobian.  If IOPT=2, only one call to
C         FCN is required for each evaluation of the full Jacobian.
C         If IOPT=3, the M calls to FCN are required.
C         If IOPT=1, then NJEV is set to zero.
C
C       IPVT is an integer output array of length N.  IPVT defines a
C         permutation matrix P such that JAC*P = Q*R, where JAC is the
C         final calculated Jacobian, Q is orthogonal (not stored), and R
C         is upper triangular with diagonal elements of nonincreasing
C         magnitude.  Column J of P is column IPVT(J) of the identity
C         matrix.
C
C       QTF is an output array of length N which contains the first N
C         elements of the vector (Q transpose)*FVEC.
C
C       WA1, WA2, and WA3 are work arrays of length N.
C
C       WA4 is a work array of length M.
C
C
C 4. Successful Completion.
C
C       The accuracy of SNLS1 is controlled by the convergence parame-
C       ters FTOL, XTOL, and GTOL.  These parameters are used in tests
C       which make three types of comparisons between the approximation
C       X and a solution XSOL.  SNLS1 terminates when any of the tests
C       is satisfied.  If any of the convergence parameters is less than
C       the machine precision (as defined by the function R1MACH(4)),
C       then SNLS1 only attempts to satisfy the test defined by the
C       machine precision.  Further progress is not usually possible.
C
C       The tests assume that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently.  If these conditions
C       are not satisfied, then SNLS1 may incorrectly indicate conver-
C       gence.  If the Jacobian is coded correctly or IOPT=1,
C       then the validity of the answer can be checked, for example, by
C       rerunning SNLS1 with tighter tolerances.
C
C       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
C         of a vector Z, then this test attempts to guarantee that
C
C               ENORM(FVEC) .LE. (1+FTOL)*ENORM(FVECS),
C
C         where FVECS denotes the functions evaluated at XSOL.  If this
C         condition is satisfied with FTOL = 10**(-K), then the final
C         residual norm ENORM(FVEC) has K significant decimal digits and
C         INFO is set to 1 (or to 3 if the second test is also satis-
C         fied).  Unless high precision solutions are required, the
C         recommended value for FTOL is the square root of the machine
C         precision.
C
C       Second Convergence Test.  If D is the diagonal matrix whose
C         entries are defined by the array DIAG, then this test attempts
C         to guarantee that
C
C               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL).
C
C         If this condition is satisfied with XTOL = 10**(-K), then the
C         larger components of D*X have K significant decimal digits and
C         INFO is set to 2 (or to 3 if the first test is also satis-
C         fied).  There is a danger that the smaller components of D*X
C         may have large relative errors, but if MODE = 1, then the
C         accuracy of the components of X is usually related to their
C         sensitivity.  Unless high precision solutions are required,
C         the recommended value for XTOL is the square root of the
C         machine precision.
C
C       Third Convergence Test.  This test is satisfied when the cosine
C         of the angle between FVEC and any column of the Jacobian at X
C         is at most GTOL in absolute value.  There is no clear rela-
C         tionship between this test and the accuracy of SNLS1, and
C         furthermore, the test is equally well satisfied at other crit-
C         ical points, namely maximizers and saddle points.  Therefore,
C         termination caused by this test (INFO = 4) should be examined
C         carefully.  The recommended value for GTOL is zero.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of SNLS1 can be due to improper input
C       parameters, arithmetic interrupts, or an excessive number of
C       function evaluations.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1
C         or IOPT .GT. 3, or N .LE. 0, or M .LT. N, or for IOPT=1 or 2
C         LDFJAC .LT. M, or for IOPT=3 LDFJAC .LT. N, or FTOL .LT. 0.E0,
C         or XTOL .LT. 0.E0, or GTOL .LT. 0.E0, or MAXFEV .LE. 0, or
C         FACTOR .LE. 0.E0.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by SNLS1.  In this
C         case, it may be possible to remedy the situation by rerunning
C         SNLS1 with a smaller value of FACTOR.
C
C       Excessive Number of Function Evaluations.  A reasonable value
C         for MAXFEV is 100*(N+1) for IOPT=2 or 3 and 200*(N+1) for
C         IOPT=1.  If the number of calls to FCN reaches MAXFEV, then
C         this indicates that the routine is converging very slowly
C         as measured by the progress of FVEC, and INFO is set to 5.
C         In this case, it may be helpful to restart SNLS1 with MODE
C         set to 1.
C
C
C 6. Characteristics of the Algorithm.
C
C       SNLS1 is a modification of the Levenberg-Marquardt algorithm.
C       Two of its main characteristics involve the proper use of
C       implicitly scaled variables (if MODE = 1) and an optimal choice
C       for the correction.  The use of implicitly scaled variables
C       achieves scale invariance of SNLS1 and limits the size of the
C       correction in any direction where the functions are changing
C       rapidly.  The optimal choice of the correction guarantees (under
C       reasonable conditions) global convergence from starting points
C       far from the solution and a fast rate of convergence for
C       problems with small residuals.
C
C       Timing.  The time required by SNLS1 to solve a given problem
C         depends on M and N, the behavior of the functions, the accu-
C         racy requested, and the starting point.  The number of arith-
C         metic operations needed by SNLS1 is about N**3 to process each
C         evaluation of the functions (call to FCN) and to process each
C         evaluation of the Jacobian it takes M*N**2 for IOPT=2 (one
C         call to FCN), M*N**2 for IOPT=1 (N calls to FCN) and
C         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
C         can be evaluated quickly, the timing of SNLS1 will be
C         strongly influenced by the time spent in FCN.
C
C       Storage.  SNLS1 requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
C         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
C         locations and N integer storage locations, in addition to
C         the storage required by the program.  There are no internally
C         declared storage arrays.
C
C *Long Description:
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), and X(3)
C       which provide the best fit (in the least squares sense) of
C
C             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
C
C       to the data
C
C             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
C                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
C
C       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
C       I-th component of FVEC is thus defined by
C
C             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
C
C       **********
C
C       PROGRAM TEST
C C
C C     Driver for SNLS1 example.
C C
C       INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,
C      *        NWRITE
C       INTEGER IPVT(3)
C       REAL FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN
C       REAL X(3),FVEC(15),FJAC(15,3),DIAG(3),QTF(3),
C      *     WA1(3),WA2(3),WA3(3),WA4(15)
C       REAL ENORM,R1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 1
C       M = 15
C       N = 3
C C
C C     The following starting values provide a rough fit.
C C
C       X(1) = 1.E0
C       X(2) = 1.E0
C       X(3) = 1.E0
C C
C       LDFJAC = 15
C C
C C     Set FTOL and XTOL to the square root of the machine precision
C C     and GTOL to zero.  Unless high precision solutions are
C C     required, these are the recommended settings.
C C
C       FTOL = SQRT(R1MACH(4))
C       XTOL = SQRT(R1MACH(4))
C       GTOL = 0.E0
C C
C       MAXFEV = 400
C       EPSFCN = 0.0
C       MODE = 1
C       FACTOR = 1.E2
C       NPRINT = 0
C C
C       CALL SNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
C      *           GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
C      *           INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
C       FNORM = ENORM(M,FVEC)
C       WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
C      *        5X,' NUMBER OF JACOBIAN EVALUATIONS',I10 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
C       END
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
C C     This is the form of the FCN routine if IOPT=1,
C C     that is, if the user does not calculate the Jacobian.
C       INTEGER M,N,IFLAG
C       REAL X(N),FVEC(M)
C       INTEGER I
C       REAL TMP1,TMP2,TMP3,TMP4
C       REAL Y(15)
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C       END
C
C
C       Results obtained with different compilers or machines
C       may be slightly different.
C
C       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
C
C       NUMBER OF FUNCTION EVALUATIONS        25
C
C       NUMBER OF JACOBIAN EVALUATIONS         0
C
C       EXIT PARAMETER                         1
C
C       FINAL APPROXIMATE SOLUTION
C
C        0.8241058E-01  0.1133037E+01  0.2343695E+01
C
C
C       For IOPT=2, FCN would be modified as follows to also
C       calculate the full Jacobian when IFLAG=2.
C
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C C
C C     This is the form of the FCN routine if IOPT=2,
C C     that is, if the user calculates the full Jacobian.
C C
C       INTEGER LDFJAC,M,N,IFLAG
C       REAL X(N),FVEC(M)
C       REAL FJAC(LDFJAC,N)
C       INTEGER I
C       REAL TMP1,TMP2,TMP3,TMP4
C       REAL Y(15)
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       IF(IFLAG.NE.1) GO TO 20
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C C
C C     Below, calculate the full Jacobian.
C C
C    20    CONTINUE
C C
C       DO 30 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
C          FJAC(I,1) = -1.E0
C          FJAC(I,2) = TMP1*TMP2/TMP4
C          FJAC(I,3) = TMP1*TMP3/TMP4
C    30    CONTINUE
C       RETURN
C       END
C
C
C       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
C         LDFJAC would be set to 3, and FCN would be written as
C         follows to calculate a row of the Jacobian when IFLAG=3.
C
C       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
C C     This is the form of the FCN routine if IOPT=3,
C C     that is, if the user calculates the Jacobian row by row.
C       INTEGER M,N,IFLAG
C       REAL X(N),FVEC(M)
C       REAL FJAC(N)
C       INTEGER I
C       REAL TMP1,TMP2,TMP3,TMP4
C       REAL Y(15)
C       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
C      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
C      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
C      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     Insert print statements here when NPRINT is positive.
C C
C       RETURN
C     5 CONTINUE
C       IF( IFLAG.NE.1) GO TO 20
C       DO 10 I = 1, M
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
C    10    CONTINUE
C       RETURN
C C
C C     Below, calculate the LDFJAC-th row of the Jacobian.
C C
C    20 CONTINUE
C
C       I = LDFJAC
C          TMP1 = I
C          TMP2 = 16 - I
C          TMP3 = TMP1
C          IF (I .GT. 8) TMP3 = TMP2
C          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
C          FJAC(1) = -1.E0
C          FJAC(2) = TMP1*TMP2/TMP4
C          FJAC(3) = TMP1*TMP3/TMP4
C       RETURN
C       END
C
C***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
C                 implementation and theory.  In Numerical Analysis
C                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
C                 Watson, Editor), Lecture Notes in Mathematics 630,
C                 Springer-Verlag, 1978.
C***ROUTINES CALLED  CHKDER, ENORM, FDJAC3, LMPAR, QRFAC, R1MACH,
C                    RWUPDT, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SNLS1
      INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
      INTEGER IJUNK,NROW,IPVT(*)
      REAL FTOL,XTOL,GTOL,FACTOR,EPSFCN
      REAL X(*),FVEC(*),FJAC(LDFJAC,*),DIAG(*),QTF(*),WA1(*),WA2(*),
     1     WA3(*),WA4(*)
      LOGICAL SING
      EXTERNAL FCN
      INTEGER I,IFLAG,ITER,J,L,MODECH
      REAL ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM,ONE,PAR,
     1     PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO,SUM,TEMP,TEMP1,
     2     TEMP2,XNORM,ZERO
      REAL R1MACH,ENORM,ERR,CHKLIM
      CHARACTER*8 XERN1
      CHARACTER*16 XERN3
C
      SAVE CHKLIM, ONE, P1, P5, P25, P75, P0001, ZERO
      DATA CHKLIM/.1E0/
      DATA ONE,P1,P5,P25,P75,P0001,ZERO
     1     /1.0E0,1.0E-1,5.0E-1,2.5E-1,7.5E-1,1.0E-4,0.0E0/
C
C***FIRST EXECUTABLE STATEMENT  SNLS1
      EPSMCH = R1MACH(4)
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
      NJEV = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (IOPT .LT. 1 .OR. IOPT .GT. 3 .OR. N .LE. 0 .OR.
     1    M .LT. N .OR. LDFJAC .LT. N .OR. FTOL .LT. ZERO
     2    .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO
     3    .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 300
      IF (IOPT .LT. 3 .AND. LDFJAC .LT. M) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
C
C     EVALUATE THE FUNCTION AT THE STARTING POINT
C     AND CALCULATE ITS NORM.
C
      IFLAG = 1
      IJUNK = 1
      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = ENORM(M,FVEC)
C
C     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
C
      PAR = ZERO
      ITER = 1
C
C     BEGINNING OF THE OUTER LOOP.
C
   30 CONTINUE
C
C        IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
         IF (NPRINT .LE. 0) GO TO 40
         IFLAG = 0
         IF (MOD(ITER-1,NPRINT) .EQ. 0)
     1      CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
         IF (IFLAG .LT. 0) GO TO 300
   40    CONTINUE
C
C        CALCULATE THE JACOBIAN MATRIX.
C
      IF (IOPT .EQ. 3) GO TO 475
C
C     STORE THE FULL JACOBIAN USING M*N STORAGE
C
      IF (IOPT .EQ. 1) GO TO 410
C
C     THE USER SUPPLIES THE JACOBIAN
C
         IFLAG = 2
         CALL FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
         NJEV = NJEV + 1
C
C             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
C
         IF (ITER .LE. 1) THEN
            IF (IFLAG .LT. 0) GO TO 300
C
C           GET THE INCREMENTED X-VALUES INTO WA1(*).
C
            MODECH = 1
            CALL CHKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
C
C           EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,LDFJAC)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
            DO 350 I = 1, M
               MODECH = 2
               CALL CHKDER(1,N,X,FVEC(I),FJAC(I,1),LDFJAC,WA1,
     1              WA4(I),MODECH,ERR)
               IF (ERR .LT. CHKLIM) THEN
                  WRITE (XERN1, '(I8)') I
                  WRITE (XERN3, '(1PE15.6)') ERR
                  CALL XERMSG ('SLATEC', 'SNLS1', 'DERIVATIVE OF ' //
     *               'FUNCTION ' // XERN1 // ' MAY BE WRONG, ERR = ' //
     *               XERN3 // ' TOO CLOSE TO 0.', 7, 0)
               ENDIF
  350       CONTINUE
         ENDIF
C
         GO TO 420
C
C     THE CODE APPROXIMATES THE JACOBIAN
C
410      IFLAG = 1
         CALL FDJAC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
         NFEV = NFEV + N
  420    IF (IFLAG .LT. 0) GO TO 300
C
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
         CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
C
C        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
C        QTF.
C
         DO 430 I = 1, M
            WA4(I) = FVEC(I)
  430         CONTINUE
         DO 470 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 460
            SUM = ZERO
            DO 440 I = J, M
               SUM = SUM + FJAC(I,J)*WA4(I)
  440          CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 450 I = J, M
               WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  450          CONTINUE
  460       CONTINUE
            FJAC(J,J) = WA1(J)
            QTF(J) = WA4(J)
  470       CONTINUE
         GO TO 560
C
C        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE.
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX
C        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY
C        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST
C        N COMPONENTS IN QTF.
C
  475    DO 490 J = 1, N
            QTF(J) = ZERO
            DO 480 I = 1, N
               FJAC(I,J) = ZERO
  480          CONTINUE
  490        CONTINUE
         DO 500 I = 1, M
            NROW = I
            IFLAG = 3
            CALL FCN(IFLAG,M,N,X,FVEC,WA3,NROW)
            IF (IFLAG .LT. 0) GO TO 300
C
C            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN.
C
            IF(ITER .GT. 1) GO TO 498
C
C            GET THE INCREMENTED X-VALUES INTO WA1(*).
C
            MODECH = 1
            CALL CHKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
C
C            EVALUATE AT INCREMENTED VALUES, IF NOT ALREADY EVALUATED.
C
            IF(I .NE. 1) GO TO 495
C
C            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(*).
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA1,WA4,FJAC,NROW)
            NFEV = NFEV + 1
            IF(IFLAG .LT. 0) GO TO 300
495         CONTINUE
            MODECH = 2
            CALL CHKDER(1,N,X,FVEC(I),WA3,1,WA1,WA4(I),MODECH,ERR)
            IF (ERR .LT. CHKLIM) THEN
               WRITE (XERN1, '(I8)') I
               WRITE (XERN3, '(1PE15.6)') ERR
               CALL XERMSG ('SLATEC', 'SNLS1', 'DERIVATIVE OF FUNCTION '
     *            // XERN1 // ' MAY BE WRONG, ERR = ' // XERN3 //
     *            ' TOO CLOSE TO 0.', 7, 0)
            ENDIF
498         CONTINUE
C
            TEMP = FVEC(I)
            CALL RWUPDT(N,FJAC,LDFJAC,WA3,QTF,TEMP,WA1,WA2)
  500       CONTINUE
         NJEV = NJEV + 1
C
C        IF THE JACOBIAN IS RANK DEFICIENT, CALL QRFAC TO
C        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF.
C
         SING = .FALSE.
         DO 510 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) SING = .TRUE.
            IPVT(J) = J
            WA2(J) = ENORM(J,FJAC(1,J))
  510       CONTINUE
         IF (.NOT.SING) GO TO 560
         CALL QRFAC(N,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
         DO 550 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 540
            SUM = ZERO
            DO 520 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
  520         CONTINUE
            TEMP = -SUM/FJAC(J,J)
            DO 530 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  530          CONTINUE
  540       CONTINUE
            FJAC(J,J) = WA1(J)
  550       CONTINUE
  560    CONTINUE
C
C        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
         IF (ITER .NE. 1) GO TO 80
         IF (MODE .EQ. 2) GO TO 60
         DO 50 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
C
C        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
C        AND INITIALIZE THE STEP BOUND DELTA.
C
         DO 70 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
         XNORM = ENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   80    CONTINUE
C
C        COMPUTE THE NORM OF THE SCALED GRADIENT.
C
         GNORM = ZERO
         IF (FNORM .EQ. ZERO) GO TO 170
         DO 160 J = 1, N
            L = IPVT(J)
            IF (WA2(L) .EQ. ZERO) GO TO 150
            SUM = ZERO
            DO 140 I = 1, J
               SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
            GNORM = MAX(GNORM,ABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
C
C        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
C
         IF (GNORM .LE. GTOL) INFO = 4
         IF (INFO .NE. 0) GO TO 300
C
C        RESCALE IF NECESSARY.
C
         IF (MODE .EQ. 2) GO TO 190
         DO 180 J = 1, N
            DIAG(J) = MAX(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
C
C        BEGINNING OF THE INNER LOOP.
C
  200    CONTINUE
C
C           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
C
            CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2,
     1                 WA3,WA4)
C
C           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
            DO 210 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
            PNORM = ENORM(N,WA3)
C
C           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
            IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)
C
C           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
            IFLAG = 1
            CALL FCN(IFLAG,M,N,WA2,WA4,FJAC,IJUNK)
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = ENORM(M,WA4)
C
C           COMPUTE THE SCALED ACTUAL REDUCTION.
C
            ACTRED = -ONE
            IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C           COMPUTE THE SCALED PREDICTED REDUCTION AND
C           THE SCALED DIRECTIONAL DERIVATIVE.
C
            DO 230 J = 1, N
               WA3(J) = ZERO
               L = IPVT(J)
               TEMP = WA1(L)
               DO 220 I = 1, J
                  WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
            TEMP1 = ENORM(N,WA3)/FNORM
            TEMP2 = (SQRT(PAR)*PNORM)/FNORM
            PRERED = TEMP1**2 + TEMP2**2/P5
            DIRDER = -(TEMP1**2 + TEMP2**2)
C
C           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C           REDUCTION.
C
            RATIO = ZERO
            IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
C
C           UPDATE THE STEP BOUND.
C
            IF (RATIO .GT. P25) GO TO 240
               IF (ACTRED .GE. ZERO) TEMP = P5
               IF (ACTRED .LT. ZERO)
     1            TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
               IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
               DELTA = TEMP*MIN(DELTA,PNORM/P1)
               PAR = PAR/TEMP
               GO TO 260
  240       CONTINUE
               IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 250
               DELTA = PNORM/P5
               PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
C
C           TEST FOR SUCCESSFUL ITERATION.
C
            IF (RATIO .LT. P0001) GO TO 290
C
C           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
            DO 270 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
            DO 280 I = 1, M
               FVEC(I) = WA4(I)
  280          CONTINUE
            XNORM = ENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  290       CONTINUE
C
C           TESTS FOR CONVERGENCE.
C
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     1          .AND. P5*RATIO .LE. ONE) INFO = 1
            IF (DELTA .LE. XTOL*XNORM) INFO = 2
            IF (ABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     1          .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
            IF (INFO .NE. 0) GO TO 300
C
C           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
            IF (NFEV .GE. MAXFEV) INFO = 5
            IF (ABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH
     1          .AND. P5*RATIO .LE. ONE) INFO = 6
            IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
            IF (GNORM .LE. EPSMCH) INFO = 8
            IF (INFO .NE. 0) GO TO 300
C
C           END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
C
            IF (RATIO .LT. P0001) GO TO 200
C
C        END OF THE OUTER LOOP.
C
         GO TO 30
  300 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'SNLS1',
     +   'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'SNLS1',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 4) CALL XERMSG ('SLATEC', 'SNLS1',
     +   'THIRD CONVERGENCE CONDITION, CHECK RESULTS BEFORE ACCEPTING.',
     +   1, 1)
      IF (INFO .EQ. 5) CALL XERMSG ('SLATEC', 'SNLS1',
     +   'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
      IF (INFO .GE. 6) CALL XERMSG ('SLATEC', 'SNLS1',
     +   'TOLERANCES TOO SMALL, NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE SNLS1.
C
      END
