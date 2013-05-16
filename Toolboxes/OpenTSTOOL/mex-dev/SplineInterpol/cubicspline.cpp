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

#include "cubic_spline.h"

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{       
	/* check input args */

	if (nrhs < 3)
	{
		mexErrMsgTxt("Cubic spline : YY = MYSPLINE(X,Y,XX)");
		return;
	}
	
	/* handle matrix I/O */
	
	const double* x = (double *)mxGetPr(prhs[0]);
	const long Nx = mxGetM(prhs[0])*mxGetN(prhs[0]); 	
	
	const double* y = (double *)mxGetPr(prhs[1]);
	const long Ny = mxGetM(prhs[1])*mxGetN(prhs[1]); 		
	
	const double* xx = (double *)mxGetPr(prhs[2]);
	const long Nxx = mxGetM(prhs[2])*mxGetN(prhs[2]); 		
	
	if (Nx < 2) {
		mexErrMsgTxt("There should be at least two data points.");
		return;
	}
	
	if (Ny < Nx) {
		mexErrMsgTxt("Abscissa and ordinate vector should be of the same length.");
		return;
	}
	
	plhs[0] = mxCreateDoubleMatrix(1, Nxx, mxREAL);      
    double* const out = (double *) mxGetPr(plhs[0]);   

 	cubic_spline<const double*> spline(Nx, 0, 0, (const double*) x, (const double*) y);
		
	for (long i = 0; i < Nxx; i++) {
		out[i] = spline.eval(xx[i]);	
	}
}
