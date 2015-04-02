/* solve_chol - solve a linear system A*X = B using the Cholesky factorization
 * of A (square, symmetric and positive definite) using LAPACK/DPOTRS.
 *
 * Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2014-02-13.
 *
 * Modifications from 2014-02-09 as suggested by Todd Small:
 *  - types of q,m,n changed from long to integer, q error check had wrong sign
 */

#include "mex.h"
#include <string.h>

#ifdef MEX_INFORMATION_VERSION             /* now we are compiling for Matlab */
  #if defined(_WIN64)
    #define integer  long long
  #else
    #define integer  long
  #endif
  #define doublereal double
#else                                      /* now we are compiling for Octave */
  #ifdef __APPLE__
    #include <Accelerate/Accelerate.h>
    typedef __CLPK_integer    integer;
    typedef __CLPK_doublereal doublereal;
  #else
    typedef int    integer;
    typedef double doublereal;
  #endif
#endif

#if !defined(_WIN32) || !defined(MEX_INFORMATION_VERSION) /* not Win32/Matlab */
  #define dpotrs dpotrs_
#endif

extern integer dpotrs_(char *, integer *, integer *, doublereal *, integer *,
                                            doublereal *, integer *, integer *);  

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{                                        /* Matlab call: X = solve_chol(A, B) */
  doublereal *C;
  integer n, m, q;

  if (nrhs != 2 || nlhs > 1)                                   /* check input */
    mexErrMsgTxt("Usage: X = solve_chol(A, B)");
  n = mxGetM(prhs[0]);                    /* number of rows in inputs A and B */
  if (n != mxGetN(prhs[0]))
    mexErrMsgTxt("First argument matrix must be square.");
  if (n != mxGetM(prhs[1]))
    mexErrMsgTxt("Both argument matrices must have the same number of rows.");
  m = mxGetN(prhs[1]);                  /* number of colums in second input B */

  plhs[0] = mxCreateDoubleMatrix(n, m, mxREAL);         /* space for output X */
  C = mxGetPr(plhs[0]);

  if (n == 0) return;             /* if argument was empty matrix, do no more */
  memcpy( C, mxGetPr(prhs[1]), m*n*mxGetElementSize(plhs[0]) );  /* copy data */
  dpotrs("U", &n, &m, mxGetPr(prhs[0]), &n, C, &n, &q);       /* solve system */
  if (q < 0) mexErrMsgTxt("Error: illegal input to solve_chol");
}
