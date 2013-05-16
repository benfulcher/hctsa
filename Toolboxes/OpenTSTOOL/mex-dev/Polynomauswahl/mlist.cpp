// Compute list of monoms of degree maxgrad or less
// mex -I.. mlist.cpp -O

#include "mextools/mextools.h"

static long maxlen;
static long* a;
static long counter;
static int count_flag;
static double* out;

void g(const long pos, const long rem_limit)
{
    if (pos == maxlen) {
			if (count_flag) {
				counter++;
			} else {
            	for (long i=0; i < maxlen;i++) {
                	*(out++) =  a[i];
				}
			}
    } else {
    	for (long limit=0; limit <= rem_limit; limit++) {
            	a[pos] = limit;
            	g(pos+1, rem_limit - limit);
    	}
	}
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	/* check input args */
	if (nrhs < 2)
	{
		mexErrMsgTxt("monomlist : Input arguments D and maxgrad must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long D = (long) *((double *)mxGetPr(prhs[0]));
	const long maxgrad = (long) *((double *)mxGetPr(prhs[1]));

	if ((D < 1) || (maxgrad <0)) {
		mexErrMsgTxt("Wrong parameters given");
		return;
	}
	
	maxlen = D;
    a = new long[D];

	count_flag = 1;
	counter = 0;
	g(0, maxgrad);	// first do a run without output to see how many polynoms we will get

	count_flag = 0;
	plhs[0] = mxCreateDoubleMatrix(D, counter, mxREAL);
	out = (double *) mxGetPr(plhs[0]);		
	
    g(0, maxgrad);

    delete[] a;
}



