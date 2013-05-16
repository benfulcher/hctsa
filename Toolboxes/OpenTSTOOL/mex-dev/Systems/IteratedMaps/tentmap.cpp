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

double g(const double y, const double h, const double e)
{
	const double t = 2.0 * h;
	const double p = pow(t, e);

	return (pow(y + t, e) - p) / (pow(1+t, e) - p);
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{			
	/* check input args */
	if (nrhs < 2)
	{
		mexErrMsgTxt("Tentmap : number of samples to compute and parameter vector must be given");
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
		mexErrMsgTxt("Need parameter vector of length 4 (h, e, s, x0)");
		return;
	}	
	
	const double h = params[0];
	const double e = params[1];
	const double s = params[2];

	double x0 = params[3];
	
	plhs[0] = mxCreateDoubleMatrix(length, 1, mxREAL);
    double* out = (double *) mxGetPr(plhs[0]);
	
	
	
	for (long i = 0; i < length; i++) {
		const double x = s * (1 - g(fabs(2 * x0 - 1), h, e));  		
		out[i] = x0 = x;
	}
}	

