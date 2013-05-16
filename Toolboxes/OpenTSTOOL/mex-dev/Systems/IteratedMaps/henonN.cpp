#include <cmath>
#include "mex.h"

#ifdef MATLAB_MEX_FILE
#undef malloc
#undef realloc
#undef free
#define malloc mxMalloc
#define realloc mxRealloc
#define free mxFree
#define printf mexPrintf
#endif

// mex henonN.cpp -O

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{			
	/* check input args */
	if (nrhs < 3)
	{
		mexErrMsgTxt("HenonN : N, number of samples to compute and parameter vector must be given");
		return;
	}
	const long N = (long) *(double *)mxGetPr(prhs[0]);
	
	if (N < 2) {
		mexErrMsgTxt("N must be at least 2");
		return;
	}		

	/* handle matrix I/O */
	
	const long length = (long) *(double *)mxGetPr(prhs[1]);
	
	if (length < 1) {
		mexErrMsgTxt("requested number of samples must be positive");
		return;
	}		

	const double* params = (double *)mxGetPr(prhs[2]);

	if (mxGetM(prhs[2])*mxGetN(prhs[2]) < 2 + N) {
		mexErrMsgTxt("Need parameter vector of shape [a,b,x0,x1, ...]");
		return;
	}	
	
	const double a = params[0];
	const double b = params[1];

	double* x = new double[N];
	double* xn = new double[N];
	
	for (long n=0; n < N; n++)
		x[n] = params[2+n];
		
	plhs[0] = mxCreateDoubleMatrix(length, N, mxREAL);
    double* out = (double *) mxGetPr(plhs[0]);
	
	for (long i = 0; i < length; i++) {
		long n;
		
		//xn[0] = 1 + a * x[N-2] * x[N-2] + b * x[N-1];
		xn[0] = a - x[N-2] * x[N-2] - b * x[N-1];
		
		for (n=1; n < N; n++)
			xn[n] = x[n-1];
			
		for (n=0; n < N; n++) {
			out[i + length*n] = x[n] = xn[n];
		}
	}
	
	delete[] x;
	delete[] xn;
}	

