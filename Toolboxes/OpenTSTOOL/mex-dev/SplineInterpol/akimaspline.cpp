#include <cmath>
#include "mextools/mextools.h"
#include "akimaspline.h"

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{       
	/* check input args */

	if (nrhs < 3)
	{
		mexErrMsgTxt("Akima interpolation : YY = AKIMA(X,Y,XX)");
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
	
	plhs[0] = mxCreateDoubleMatrix(Nxx, 1, mxREAL);      
    double* const out = (double *) mxGetPr(plhs[0]);   

 	akima_interpolator<const double*> akima(Nx, (const double*) x, (const double*) y);
		
	if (akima.geterr() == 1)	
	{
		mexErrMsgTxt("Need more input data points");
		return;
	} else if (akima.geterr() == 2) {
		mexErrMsgTxt("Data point must be given with ascending x values");
		return;		
	} else if (akima.geterr()) {
		mexErrMsgTxt("Error computing Akima coefficients");
		return;		
	}

	for (long i = 0; i < Nxx; i++) {
		out[i] = akima.eval(xx[i]);	
	}
}
