#include <cmath>

#include "mextools/mextools.h"
#include "mextools/Utilities.h"
#include "mextools/Utilities.cpp"

static uint32 seed = 667297391;  // just a seed

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{		
	/* check input args */
	if (nrhs < 1)
	{
		mexErrMsgTxt("Mersenne Twister Random Number Generator : Number of random numbers [0,1) must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N = (long) *(double *)mxGetPr(prhs[0]);
	
	plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
	double* x = (double *) mxGetPr(plhs[0]);
	
	RNG rng(seed);
	
	for (long i= 0; i < N; i++) x[i] = rng.dRand();
	
	seed = rng.lRand() | (uint32) 1;

#ifdef VERBOSE	
	printf("New seed : %u \n", seed);
#endif 
}	




