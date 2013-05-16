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


void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{			
	/* check input args */
	if (nrhs < 2)
	{
		mexErrMsgTxt("baker : number of samples to compute and parameter vector must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long length = (long) *(double *)mxGetPr(prhs[0]);
	
	if (length < 1) {
		mexErrMsgTxt("requested number of samples must be positive");
		return;
	}		

	const double* params = (double *)mxGetPr(prhs[1]);

	if (mxGetM(prhs[1])*mxGetN(prhs[1]) < 5) {
		mexErrMsgTxt("Need parameter vector of length 5 (eta,l1,l2,x0,y0)");
		return;
	}	
	
	const double eta = params[0];
	const double l1 = params[1];
	const double l2 = params[2];

	double xn = params[3];
	double yn = params[4];
	
	plhs[0] = mxCreateDoubleMatrix(length, 2, mxREAL);
    double* out = (double *) mxGetPr(plhs[0]);
	
	for (long i = 0; i < length; i++) {
		double xn1, yn1;
		
		if (yn < eta) {
			xn1 = l1 * xn;
			yn1 = yn / eta;
		} else {
			xn1 = 0.5 + l2*xn;
			yn1 = (yn-eta)/(1-eta);
		}
		
		out[i] = xn = xn1;
		out[length+i] = yn = yn1;
	}
}	

