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

using namespace std;

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{			
	/* check input args */
	if (nrhs < 2)
	{
		mexErrMsgTxt("Henon : number of samples to compute and parameter vector must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long length = (long) *(double *)mxGetPr(prhs[0]);
	
	if (length < 1) {
		mexErrMsgTxt("requested number of samples must be positive");
		return;
	}		

	const double* params = (double *)mxGetPr(prhs[1]);

	if (mxGetM(prhs[1])*mxGetN(prhs[1]) < 4) {
		mexErrMsgTxt("Need parameter vector of length 4 (a,b,x0,y0)");
		return;
	}	
	
	const double a = params[0];
	const double b = params[1];

	double xn = params[2];
	double yn = params[3];
	
	plhs[0] = mxCreateDoubleMatrix(length, 2, mxREAL);
    double* out = (double *) mxGetPr(plhs[0]);
	
	for (long i = 0; i < length; i++) {
		const double xn1 = 1 + a*xn*xn + b * yn;
		const double yn1 = xn;
		
		if (fabs(xn1) > 1000000) {
			mexWarnMsgTxt("Values exceed threshold of 1000000, limiting values");
			out[i] = xn = 1000000;
		}
		else {
			out[i] = xn = xn1;
		}
		out[length+i] = yn = yn1;
	}
}	

