// largest lyapunov exponent for reconstructed time series data
// cmerk DPI Goettingen 1998

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

static inline long lmax(long a,long b) {return ((a >= b) ? a : b);}
static inline long lmin(long a,long b) {return ((a <= b) ? a : b);}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	int atria_preprocessing_given = 0;		// flag wheter preprocessing is already given on command line

	ATRIA< point_set<euclidian_distance> >* searcher = 0;
	// try to see if the first parameter given is an atria structure
	// If this is true, the order of input parameters is shifted by one

	// try to see if the first parameter given is an atria structure
	// If this is true, the order of input parameters is shifted by one
	if ((nrhs) && (mxIsStruct(prhs[0]))) {
		atria_preprocessing_given = 1;
		prhs++; 	// these two lines enable us to use the old argument parsing block without changing it
		nrhs--;
	}

	/* check input args */
	if (nrhs < 5)
	{
		mexErrMsgTxt("Largest lyapunov exponent : Data set of points (row vectors), reference indices, length of prediction, maximal number of neighbours and past must be given");
		return;
	}

	/* handle matrix I/O */

	const long N = mxGetM(prhs[0]);
	const long dim = mxGetN(prhs[0]);
	const double* p 	= (double *)mxGetPr(prhs[0]);

	const long R = mxGetM(prhs[1]) * mxGetN(prhs[1]);			// number of reference points
	const double* ref 	= (double *)mxGetPr(prhs[1]);

	const long length 	= (long) *((double *)mxGetPr(prhs[2])); 	// length of prediction
	const long NNR = (long) *((double *)mxGetPr(prhs[3]));			// number of nearest neighbors
	const long past	= (long) *((double *)mxGetPr(prhs[4])); 		// number of points to be excluded from search

	if (fabs((double) N) - fabs((double)length) < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}
	if (dim < 1) {
		mexErrMsgTxt("Data points must be at least of dimension one");
		return;
	}
	if (length <= 0) {
		mexErrMsgTxt("Prediction length must be greater zero");
		return;
	}
	if ((mxGetN(prhs[1]) == 0) || (mxGetM(prhs[1]) == 0)) {
		mexErrMsgTxt("Wrong reference indices given");
		return;
	}
	if ((NNR < 1) || (NNR > N)) {
		mexErrMsgTxt("Number of nearest neighbours must be within 1 and number of points in the data set");
		return;
	}

	// mexPrintf("Number of reference points : %d\n", R);

	if (R < 1) {
		mexErrMsgTxt("At least one reference index or point must be given");
		return;
	}

	for (long i=0; i < R; i++) { 	// check reference indices
		if ((ref[i] < 1) || (ref[i]>N-length)) {
			mexErrMsgTxt("Reference indices out of range");
			return;
		}
	}

#include "create_searcher.cpp"

	plhs[0] = mxCreateDoubleMatrix(length+1, 1, mxREAL);
	double* x = (double *) mxGetPr(plhs[0]);

	const double logarithm_of_2 = log(2.0);
	int issue_warning = 1;
	long count = 0; 	// count the number of pairs with distance > 0

	for (long step = 0; step <= length; step++) x[step] = 0;

	for (long n=0; n < R; n++) { 		/* main loop, iterate over all reference points */
		vector<neighbor>::iterator i;
		vector<neighbor> v;

		const long actual = (long) ref[n]-1;		/* Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...*/

		// try to find at least NNR valid pairs (reference point - neighbor)
		// a valid pair must fullfil several criteria :
		// 1) both must have a known future of at least lenght steps, so that we
		//    don't run out of the data set size when computing the divergence
		// 2) initial distance must be greater zero
		// 3) inital indices must differ by at least "past" steps
		for (long nnr=NNR; nnr < 4 * NNR; nnr++) {
			vector<neighbor> v_test;
			v.erase(v.begin(), v.end());

			searcher->search_k_neighbors(v_test, nnr, points.point_begin(actual), actual-past, actual+past);

			for (i = v_test.begin(); i < v_test.end(); i++) {
				if (((*i).dist() > 0) && ((*i).index() < N-length)) {
					v.push_back(*i);
				}
			}
			if (v.size() >= NNR)
				break;
		}

		for (i = v.begin(); i < v.end(); i++) { 	// v is the sorted vector of neighbors
			const long nnindex = (*i).index();
			const double nndist = (*i).dist();
			const double log_nndist = log(nndist);

			count++;

			for (long step = 0; step <= length; step++) {
				const double dist = points.distance(actual+step, nnindex+step);

				if (dist > 0) {
					x[step] += log(dist) - log_nndist;
				}
				else {
					if (issue_warning) {
						issue_warning = 0;	// don't issue this warning again
						mexPrintf("Largelyap warning : cannot take logarithm of distance zero \n");
					}
				}
			}
		}
	}

	if (count > 0)	{
		for (long step = 0; step <= length; step++) x[step] /= (count*logarithm_of_2);
	}

	delete searcher;
}
