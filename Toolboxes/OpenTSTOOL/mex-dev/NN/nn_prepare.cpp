// This mex-function just does the preprocessing for a nearest neighbor searcher
// The data structure that is computed during the preprocessing phase can then
// be used to do several calls to the mex-function nn_search without doing the 
// time-consuming preprocessing again (of course only when the input data set
// has not changed)

#include "mex.h"
#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

template<class Searcher>
mxArray* nn_prepare(Searcher& searcher)
{
	if (searcher.geterr()) {
		mexErrMsgTxt("Error preparing searcher");
	}	 

	return searcher.store();
}	

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{		
	long minpoints = 64;
	char* metric = 0;
	
	/* check input args */
	if (nrhs < 1)
	{
		mexErrMsgTxt("Prepare nearest neighbor search : data set of points (row vectors) must be given, cluster minpoints and or metric is optional");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N 		= mxGetM(prhs[0]);	
	const long dim  	= mxGetN(prhs[0]);
	const double* p 	= (double *)mxGetPr(prhs[0]);
	
	if (N < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}		

	if (dim < 1) {
		mexErrMsgTxt("Data points must be at least of dimension one");
		return;
	}	
	
	for (long i=1; i<nrhs; i++) {
		if (mxIsChar(prhs[i])) {
    		long buflen = (mxGetM(prhs[i]) * mxGetN(prhs[i])) + 1;
 			metric = (char*) mxMalloc(buflen);
        	mxGetString(prhs[i], metric, buflen); 
		}
		else if (mxIsDouble(prhs[i])) {
			minpoints = (long) *((double *)mxGetPr(prhs[i]));
		}
	}
		
	if (minpoints < 2) {
		mexErrMsgTxt("cluster minpoints must be greater two");
		return;	
	}
	
	if ((metric == 0) || (!strncmp("euclidian", metric, strlen(metric)))) {
		point_set<euclidian_distance> points(N,dim, p);	
		ATRIA<point_set<euclidian_distance> > searcher(points, 0, minpoints);	
		plhs[0] = nn_prepare(searcher);
		mxSetField(plhs[0], 0, "optional", mxCreateString("euclidian"));
	}
	else if (!strncmp("maximum", metric, strlen(metric))){
		point_set<maximum_distance> points(N,dim, p);	
		ATRIA<point_set<maximum_distance> > searcher(points, 0, minpoints);	
		plhs[0] = nn_prepare(searcher);
		mxSetField(plhs[0], 0, "optional", mxCreateString("maximum"));
	}
	else {
		mexErrMsgTxt("Unknown type of metric");
	}
	
	mxFree(metric);
}	

