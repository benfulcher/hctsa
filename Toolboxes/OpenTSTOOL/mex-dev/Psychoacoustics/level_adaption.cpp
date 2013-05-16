// C.Merkwirth,J.Wichard DPI Goettingen 1998 
// adaptive level control
// Pueschels Nachregelketten mit Muenkners Dynamik-Kompression zur automatischen
// Aussteuerung von Signalen

#include <cmath>
#include "mex.h"

static double max(double a,double b) {return ((a >= b) ? a : b);}
static double min(double a,double b) {return ((a <= b) ? a : b);}

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
	if (nrhs < 4)
	{
		mexErrMsgTxt("adaptive level  : input samples, vector of time constants, upper dynamic limit and quietness threshold must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N = mxGetM(prhs[0]);	
	const long channels = mxGetN(prhs[0]);
	const double* const p  = (double *)mxGetPr(prhs[0]);
	
	const double* const time_constants  = (double *)mxGetPr(prhs[1]);
	const long number_dividers = (long) max(mxGetM(prhs[1]), mxGetN(prhs[1]));
	
	if (number_dividers < 1) {
		mexErrMsgTxt("Vector of time constants must be given");
		return;			
	}	
	
	for (long i=0; i < number_dividers; i++) {
		if (time_constants[i] < 0) {
			mexErrMsgTxt("Time constants must be positive");
			return;
		}
	}
	
	const double dynamic_limit = (double) *(double *)mxGetPr(prhs[2]) - 1;
	
	if (dynamic_limit < 0) {
		mexErrMsgTxt("Dynamic limit must be greater 1");
		return;			
	}	
	
	const double threshold = (double) *(double *)mxGetPr(prhs[3]);
	
	if (threshold <= 0) {
		mexErrMsgTxt("Threshold must be positive");
		return;			
	}		
	
	plhs[0] = mxCreateDoubleMatrix(N, channels, mxREAL);
	double* out = (double *) mxGetPr(plhs[0]);		

	double* const dividers = (double*) malloc(number_dividers * sizeof(double));

	for (long j=0; j < channels; j++) {
		const double* column = p + j * N;
		
		for (long k=0; k < number_dividers; k++) dividers[k] = 1;
		
		for (long i=0; i < N; i++) {
			double value = column[i];
			
			for (long k = 0; k < number_dividers; k++) {
				value /= dividers[k];
				const double tmp = fabs(value);
				if (value >= 0) 
					value = (tmp > 1) ? dynamic_limit*(2.0/(1.0+exp(-2.0 * (tmp-1)/dynamic_limit))-1)+1 : tmp;
				else
					value = (tmp > 1) ? -dynamic_limit*(2.0/(1.0+exp(-2.0 * (tmp-1)/dynamic_limit))-1)+1 : -tmp;
				
				dividers[k] = (max(fabs(value), threshold) + dividers[k]*(time_constants[k]-1))/time_constants[k];
			}
			*(out++) = value;
		}		
	}							
	
	free(dividers);
}	

