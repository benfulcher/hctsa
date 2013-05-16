// Count nearest neighbors within distance eps
// from some given reference indices (which point to vectors inside the data set)
// or to explicitly given reference vectors

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

template<class Searcher>
void range_search(Searcher& searcher, const int nlhs, mxArray  *plhs[], const int nrhs, const mxArray  *prhs[], 
const long R, const double* range_vals, const long past, const double* ref, const long N, const long dim, const double* p, 
const int ref_or_direct, const long N_range)
{
	double*	const coord = new double[dim]; 
	plhs[0] = mxCreateDoubleMatrix((long) R, 1, mxREAL);
	double* out = (double *) mxGetPr(plhs[0]);	
	
	if (nlhs > 1) { 	// return number of neighbors and a cell array of vectors containing indices and distances to the neighbors
		plhs[1] = mxCreateCellMatrix(R, 2);
		
		for (long n=0; n < R; n++) { /* iterate over all reference points */ 
			long count;
			double range;
			vector<neighbor> v;
			
			if (N_range == R) {
				range = range_vals[n];
			} else {
				range = *range_vals;
			}

			if (ref_or_direct) {
				const long actual = (long) ref[n]-1;		/* Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...*/
				for (long k=0; k < dim; k++) coord[k] = p[actual+k*N];
				count = searcher.search_range(v, range, coord, actual-past, actual+past);
			} else {
				for (long k=0; k < dim; k++) coord[k] = ref[n+k*R];
				count = searcher.search_range(v, range, coord, -1, -1);
			}	
			
			mxArray* cellcontent1 = mxCreateDoubleMatrix(1, (long) (out[n] = count), mxREAL);	// vector of indices
			mxArray* cellcontent2 = mxCreateDoubleMatrix(1, (long) (out[n] = count), mxREAL);	// vector of distances
			
			double* indices = (double *)mxGetPr(cellcontent1);
			double* distances = (double *)mxGetPr(cellcontent2);
			
			for (long i=0; i < count; i++) {
				indices[i] = v[i].index()+1; 		// convert C indices back to Matlab indices
				distances[i] = v[i].dist();
			}
			
			mxSetCell(plhs[1], n, cellcontent1);
			mxSetCell(plhs[1], (long) n+R, cellcontent2);
		}		
	} else {			// only return one output argument : the number of neighbors for each reference point
		for (long n=0; n < R; n++) { /* iterate over all reference points */ 
			double range;

			if (N_range == R) {
				range = range_vals[n];
			} else {
				range = *range_vals;
			}

			if (ref_or_direct) {
				const long actual = (long) ref[n]-1;		/* Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...*/
				for (long k=0; k < dim; k++) coord[k] = p[actual+k*N];
				out[n] = searcher.count_range(range, coord, actual-past, actual+past);
			} else {
				for (long k=0; k < dim; k++) coord[k] = ref[n+k*R];
				out[n] = searcher.count_range(range, coord, -1, -1);
			}	
		}		
	}
		
	delete[] coord;
}


void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{	
	int ref_or_direct = 1;			// ref_or_direct = 1 means interpret third input argument as reference indices
									// = 0 means interpret second input argument as reference vectors
	
	long past = 0;
	
	/* check input args */
	if (nrhs < 4)
	{
		mexErrMsgTxt("Count nearest neighbor within distance range : Data set of points (row vectors), atria structure, reference indices or reference points \nand range must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N 		= mxGetM(prhs[0]);	
	const long dim  	= mxGetN(prhs[0]);
	const double* p 	= (double *)mxGetPr(prhs[0]);
	
	double* ref 	= (double *)mxGetPr(prhs[2]);
	long R; 										// number of query (reference) points
	
	const double* range = (double *)mxGetPr(prhs[3]);	// range can be one value for all reference points
	long N_range = mxGetM(prhs[3]) * mxGetN(prhs[3]);		// or one value for each reference point; thus allowed values for N_range : 1 or R
	
	if (N < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}		
	if (dim < 1) {
		mexErrMsgTxt("Data points must be at least of dimension one");
		return;
	}	
	if (range<0) {
		mexErrMsgTxt("Range must be positive");
		return;
	}	
	if ((mxGetN(prhs[1]) == 0) || (mxGetN(prhs[1]) == 0)) {
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

	if (R < 1) {
		mexErrMsgTxt("At least one reference index or point must be given");
		return;
	}	
	if ((N_range != R) && (N_range != 1)) {
		mexErrMsgTxt("One global range value or one range value for each reference points are allowed");
		return;
	}	

	if (ref_or_direct) {		// interpret second argument as list of indices into the data set given as first argument 
		if (nrhs < 5)
		{
			mexErrMsgTxt("Count nearest neighbor within distance range : Data set of points (row vectors), atria structure, reference indices or reference points \nrange and past must be given must be given");
			return;
		}
		past = (long) *((double *)mxGetPr(prhs[4]));
		for (long i=0; i < R; i++) {
			if ((ref[i] < 1) || (ref[i]>N)) {
				mexErrMsgTxt("Reference indices out of range");
				return;
			}	
		}	
	} 
	
	char* metric = 0;
	
	if (mxIsChar(mxGetField(prhs[1], 0, "optional"))) {	
    	long buflen = (mxGetM(mxGetField(prhs[1], 0, "optional")) * mxGetN(mxGetField(prhs[1], 0, "optional"))) + 1;
 		metric = (char*) mxMalloc(buflen);
        mxGetString(mxGetField(prhs[1], 0, "optional"), metric, buflen); 
	}
	
	if ((metric == 0) || (!strncmp("euclidian", metric, strlen(metric)))) {
		point_set<euclidian_distance> points(N,dim, p);	
		ATRIA< point_set<euclidian_distance> > searcher(points, prhs[1]);	// this constructor used the data from the preprocessing that are given by the second input argument
		range_search(searcher, nlhs, plhs, nrhs, prhs, R, range, past, ref, N, dim, p, ref_or_direct, N_range);
	}
	else if (!strncmp("maximum", metric, strlen(metric))){
		point_set<maximum_distance> points(N,dim, p);	
		ATRIA< point_set<maximum_distance> > searcher(points, prhs[1]);
		range_search(searcher, nlhs, plhs, nrhs, prhs, R, range, past, ref, N, dim, p, ref_or_direct, N_range);
	}
	else {
		mexErrMsgTxt("Unknown type of metric in ATRIA structure");
	}
	
	mxFree(metric);		
}	


