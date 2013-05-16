// C.Merkwirth,J.Wichard DPI Goettingen 1998 
// simple moving average filter, works along columns, even for multidimensional
// matrices

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

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{		
	/* check input args */
	if (nrhs < 2)
	{
		mexErrMsgTxt("moving average filter : input data and length (or window vector) must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N = mxGetM(prhs[0]);	
	const long channels = mxGetN(prhs[0]);
	const double* p  = (double *)mxGetPr(prhs[0]);
	
	if (max(mxGetM(prhs[1]), mxGetN(prhs[1])) <= 1) {	// length is given as scalar paramter
		const long len = (long) *(double *)mxGetPr(prhs[1]);
		
		if (len < 2) {
			mexErrMsgTxt("Average length too short");
			return;
		}	

		const long zeros_to_pad = len-1;
		const long begin_pad = zeros_to_pad/2;
		const long end_pad =  zeros_to_pad - begin_pad;

		plhs[0] = mxCreateDoubleMatrix(N, channels, mxREAL);
		double* out = (double *) mxGetPr(plhs[0]);		

		for (long j=0; j < channels; j++) {
			const double* column = p + j * N;

			for (long i=0; i < N; i++) {
				double sum = 0;
				for (long k = 0; k < len; k++) {
					const long abs_index = i+k;
					if ((abs_index >= begin_pad) && (abs_index < begin_pad + N))
						sum += column[abs_index - begin_pad];
				}
				*(out++) = sum/len;
			}		
		}				
	}	
	else {		// a window vector is given, so length is implicitly given
                const long len = (long) max(mxGetM(prhs[1]), mxGetN(prhs[1])); 
		const double* weights  = (double *)mxGetPr(prhs[1]);
		double total_weight = 0;
		
		for (long k = 0; k < len; k++) total_weight += weights[k];
		
		if (total_weight == 0) {
			mexErrMsgTxt("Mean of weight vector is zero");
			return;
		}	
		
		const long zeros_to_pad = len-1;
		const long begin_pad = zeros_to_pad/2;
		const long end_pad =  zeros_to_pad - begin_pad;

		plhs[0] = mxCreateDoubleMatrix(N, channels, mxREAL);
		double* out = (double *) mxGetPr(plhs[0]);		

		for (long j=0; j < channels; j++) {
			const double* column = p + j * N;

			for (long i=0; i < N; i++) {
				double sum = 0;
				for (long k = 0; k < len; k++) {
					const long abs_index = i+k;
					if ((abs_index >= begin_pad) && (abs_index < begin_pad + N))
						sum += column[abs_index - begin_pad] * weights[k];
				}
				*(out++) = sum/total_weight;
			}		
		}							
	}
}	

