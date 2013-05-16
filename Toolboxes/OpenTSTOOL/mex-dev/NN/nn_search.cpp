// Fast nearest neighbor searcher
// 
// Unlike the mex-function fnearneigh, this mex-function takes the preprocessing
// from the mex-function prepare_nnsearch. This can be used to avoid unnecessary
// re-computing when the input data set has not changed.

// Output of this program are the indices of the nearest neighbors.
// If a second output argument is specified, the distances to the neighbors are
// returned

#include <cmath>

// C++ Standard Template Library
// this files have to be located before(!!!) the #ifdef MATLAB_MEX_FILE sequence,
// otherwise the defines will break the STL
#include <vector>

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

template<class Searcher>
void nn_search(Searcher& searcher, double* const nn, double* const dists, const long R, const long N,
				const long dim, const double* p, const double* ref, const long NNR, const long past, const double epsilon,
				const int ref_or_direct)
{	
	if (searcher.geterr()) {
		mexErrMsgTxt("Error preparing searcher, maybe wrong preprocessing data were given or the point set has changed");
	}	 
	
	double*	const coord = (double*) malloc(dim * sizeof(double)); 	
	
	for (long n=0; n < R; n++) { /* iterate over all reference points */ 
		vector<neighbor> v;
		
		if (ref_or_direct) {
                  const long actual = (long) ref[n]-1;		/* Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...*/
			for (long k=0; k < dim; k++) coord[k] = p[actual+k*N];
			searcher.search_k_neighbors(v, NNR, coord, actual-past, actual+past, epsilon);
		} else {
			for (long k=0; k < dim; k++) coord[k] = ref[n+k*R];
			searcher.search_k_neighbors(v, NNR, coord, -1, -1, epsilon);
		}	
		
		for (long k = 0; k < v.size(); k++) { 	// v is the sorted vector of neighbors
				nn[n+k*R] = v[k].index() +1;	// convert indices back to Matlab (1..N) style 
				dists[n+k*R] = v[k].dist();
		}
	}
	
	free(coord);
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{	
	int ref_or_direct = 1;			// ref_or_direct = 1 means interpret second input argument as reference indices
									// = 0 means interpret second input argument as reference vectors
	
	long past = 0;
	double epsilon = 0;
	
	/* check input args */
	
	if (nrhs < 4)
	{
		mexErrMsgTxt("Fast nearest neighbour searcher : Data set of points (row vectors), preprocessing data, reference indices or reference points \nand number of nearest neighbours must be given, epsilon (default = 0) is optional");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N 		= mxGetM(prhs[0]);	
	const long dim  	= mxGetN(prhs[0]);
	const double* p 	= (double *)mxGetPr(prhs[0]);
	
	double* ref 	= (double *)mxGetPr(prhs[2]);
	long R; 										// number of query (reference) points
	
	const long NNR 	= (long) *((double *)mxGetPr(prhs[3]));
	
	if (N < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}		
	if (dim < 1) {
		mexErrMsgTxt("Data points must be at least of dimension one");
		return;
	}	
	if (NNR<1) {
		mexErrMsgTxt("At least one nearest neighbour must be requested");
		return;
	}	
	if ((mxGetN(prhs[2]) == 0) || (mxGetN(prhs[2]) == 0)) {
		mexErrMsgTxt("Wrong reference indices or reference points given");
		return;
	}	

	
	if (mxGetN(prhs[2]) == 1) {
		R = mxGetM(prhs[2]);
		ref_or_direct = 1;	
	} 
	else if ((mxGetM(prhs[2]) == 1) && (mxGetN(prhs[2]) != dim)) {
		R = mxGetN(prhs[2]);
		ref_or_direct = 1;	
	}
	else if (mxGetN(prhs[2]) == dim) {
		R = mxGetM(prhs[2]);
		ref_or_direct = 0;	
	} else  {
		mexErrMsgTxt("Cannot determine if second argument are reference indices or reference points");
		return;
	}	
	
/*******************************/
/*** Quick&Dirty patch, to manage onedimensional "direct" query-points. I. Wedekind, 16.9.2002 */
	if((mxGetN(prhs[2])==dim) && (ref_or_direct==1) && (nrhs<5))
	    {
		ref_or_direct=0;
	    }
/*******************************/

	if (R < 1) {
		mexErrMsgTxt("At least one reference index or point must be given");
		return;
	}	

	if (ref_or_direct) {		// interpret second argument as list of indices into the data set given as first argument 
		if (nrhs < 5)
		{
			mexErrMsgTxt("Fast nearest neighbour searcher : Data set of points (row vectors), preprocessing data, reference indices,\nnumber of nearest neighbours and past must be given");
			return;
		}
		past	= (long) *((double *)mxGetPr(prhs[4]));
		for (long i=0; i < R; i++) {
			if ((ref[i] < 1) || (ref[i]>N)) {
				mexErrMsgTxt("Reference indices out of range");
				return;
			}	
		}
		if ((N - (2*past)-1) < NNR)
		{
			mexErrMsgTxt("Fast nearest neighbour searcher : To many neighbors for each query point are requested");
			return;
		}
		if (nrhs > 5) 
			epsilon = (double) *((double *)mxGetPr(prhs[5]));		// support approximative queries		
	}  else {
		if (nrhs > 4) 
			epsilon = (double) *((double *)mxGetPr(prhs[4]));		// support approximative queries	
	}
	
	plhs[0] = mxCreateDoubleMatrix(R, NNR, mxREAL);
	double* nn = (double *) mxGetPr(plhs[0]);
	
	double* dists;
	
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(R, NNR, mxREAL);
		dists = (double *) mxGetPr(plhs[1]);
	} else
		dists = (double *) malloc(R*NNR * sizeof(double));
	
	char* metric = 0;
	
	if (mxIsChar(mxGetField(prhs[1], 0, "optional"))) {	
    	long buflen = (mxGetM(mxGetField(prhs[1], 0, "optional")) * mxGetN(mxGetField(prhs[1], 0, "optional"))) + 1;
 		metric = (char*) mxMalloc(buflen);
        mxGetString(mxGetField(prhs[1], 0, "optional"), metric, buflen); 
	}
	
	if ((metric == 0) || (!strncmp("euclidian", metric, strlen(metric)))) {
		point_set<euclidian_distance> points(N,dim, p);	
		ATRIA< point_set<euclidian_distance> > searcher(points, prhs[1]);	// this constructor used the data from the preprocessing that are given by the second input argument
#ifdef PROFILE
		cout << "read in atria structure" << endl;
#endif	
		nn_search(searcher, nn, dists, R, N, dim, p, ref, NNR, past, epsilon, ref_or_direct);
		if (nlhs > 2) {
			plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
			*((double *) mxGetPr(plhs[2])) = searcher.search_efficiency();
		}	
	}
	else if (!strncmp("maximum", metric, strlen(metric))){
		point_set<maximum_distance> points(N,dim, p);	
		ATRIA< point_set<maximum_distance> > searcher(points, prhs[1]);
#ifdef PROFILE
		cout << "read in atria structure" << endl;
#endif		
		nn_search(searcher, nn, dists, R, N, dim, p, ref, NNR, past, epsilon, ref_or_direct);
		if (nlhs > 2) {
			plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
			*((double *) mxGetPr(plhs[2])) = searcher.search_efficiency();
		}		
	}
	else {
		mexErrMsgTxt("Unknown type of metric used to create ATRIA structure");
	}
		
	mxFree(metric);		
			
	if (!(nlhs > 1)) free(dists);
}	

