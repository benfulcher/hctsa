/* solve_chol - solve a linear system A*X = B using the cholesky factorization
   of A (where A is square, symmetric and positive definite.

   Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2010-12-16. */

#include "mex.h"
#include <math.h>
#include <string.h>

extern int dpotrs_(char *, long *, long *, double *, long *, double *, long *, long *);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *C;
  long n, m, q;

  if (nrhs != 2 || nlhs > 1)                              /* check the input */
    mexErrMsgTxt("Usage: X = solve_chol(R, B)");
  n = mxGetN(prhs[0]);
  if (n != mxGetM(prhs[0]))
    mexErrMsgTxt("Error: First argument matrix must be square");
  if (n != mxGetM(prhs[1]))
    mexErrMsgTxt("Error: First and second argument matrices must have same number of rows");
  m = mxGetN(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(n, m, mxREAL);  /* allocate space for output */
  C = mxGetPr(plhs[0]);

  if (n==0) return;               /* if argument was empty matrix, do no more */
    memcpy( C, mxGetPr(prhs[1]), n*m*sizeof(double) );/* copy argument matrix */
    dpotrs_("U", &n, &m, mxGetPr(prhs[0]), &n, C, &n, &q);    /* solve system */
  if (q > 0)
    mexErrMsgTxt("Error: illegal input to solve_chol");
}

