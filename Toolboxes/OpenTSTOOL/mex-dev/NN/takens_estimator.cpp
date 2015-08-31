// takens estimator for the correlation dimension (D_2)
//
// cmerk DPI Goettingen 1998

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

template<class METRIC>
void compute(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[], const int atria_preprocessing_given,  METRIC dummy)
{
	long i,j,n,k;		/* loop variables */

	/* handle matrix I/O */

	const long N = mxGetM(prhs[0]);
	const long dim = mxGetN(prhs[0]);
	const double* p  	= (double *)mxGetPr(prhs[0]);

	const long R = max(mxGetM(prhs[1]), mxGetN(prhs[1]));
	const double* ref 	= (double *)mxGetPr(prhs[1]);

	const double  relative_range = (double) *((double *)mxGetPr(prhs[2]));
	const long past	= (long) *((double *)mxGetPr(prhs[3]));

	if (N < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}
	if (dim < 2) {
		mexErrMsgTxt("Data points must be at least of dimension two");
		return;
	}
	if ((mxGetN(prhs[1]) == 0) || (mxGetM(prhs[1]) == 0)) {
		mexErrMsgTxt("Wrong reference indices given");
		return;
	}

	if (R < 1) {
		mexErrMsgTxt("At least one reference index or point must be given");
		return;
	}

	for (i=0; i < R; i++) {
		if ((ref[i] < 1) || (ref[i]>N)) {
			mexErrMsgTxt("Reference indices out of range");
			return;
		}
	}

	point_set<METRIC> points(N,dim, p);
	ATRIA< point_set<METRIC> >* searcher = 0;

#ifdef MATLAB_MEX_FILE
	if (atria_preprocessing_given) {
		searcher = new ATRIA< point_set<METRIC> >(points, prhs[-1]);	// this constructor used the data from the preprocessing
		if (searcher->geterr()) {
			delete searcher;
			searcher = 0;
		}
	}
#endif

	if (searcher == 0) {
		searcher = new ATRIA< point_set<METRIC> >(points);
	}

	if (searcher->geterr()) {
		mexErrMsgTxt("Error preparing searcher");
		return;
	}

	double range;

	if (relative_range > 0)
		range = relative_range * searcher->data_set_radius();  // compute the maximal search radius using information about attractor size
	else
		range = fabs(relative_range); 	// if relative_range is negativ, use it's value (without sign) as maximal search radius

	// mexPrintf("Number of reference points            : %d\n", R);
	// mexPrintf("Upper bound for attractor size        : %f\n", 2 * searcher->data_set_radius());
	// mexPrintf("Maximal search radius                 : %f\n", range);

	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* out = (double *) mxGetPr(plhs[0]);

	unsigned long counter = 0;
	double sum = 0;

	for (n=0; n < R; n++) { 		/* iterate over all reference points */
          const long actual = (long) ref[n]-1;		/* Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...*/

		if (actual > past) {
			vector<neighbor> v;

			searcher->search_range(v, range, points.point_begin(actual), actual-past, N); // don't search points from [actual-past .. N-1]
			//overall_points += (actual-past);	// count the total number of points pairs that were at least theoretically tested

			if (v.size() > 0) {
				vector<neighbor>::iterator i;

				for (i = v.begin(); i < v.end(); i++) { // v is unsorted
					const double dist = (*i).dist();
					if (dist > 0) {
						sum += log(dist/range);
						counter++;
					}
				}
			}
		}
	}

	if (counter > 0) *out = -((double)counter)/sum;
	else *out = 0;

	delete searcher;
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	int atria_preprocessing_given = 0;		// flag wheter preprocessing is already given on command line

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
	if (nrhs < 4)
	{
		mexErrMsgTxt("Takens estimator : Data set of points (row vectors), reference indices, relative range (relative to attractor diameter) and past must be given");
		return;
	}

	if (atria_preprocessing_given) {
#ifdef MATLAB_MEX_FILE
		char* metric = 0;

		if (mxIsChar(mxGetField(prhs[-1], 0, "optional"))) {
    		long buflen = (mxGetM(mxGetField(prhs[-1], 0, "optional")) * mxGetN(mxGetField(prhs[-1], 0, "optional"))) + 1;
 			metric = (char*) mxMalloc(buflen);
        	mxGetString(mxGetField(prhs[-1], 0, "optional"), metric, buflen);
		}

		if ((metric == 0) || (!strncmp("euclidian", metric, strlen(metric)))) {
			euclidian_distance dummy;
			// mexPrintf("Using euclidian metric to calculated distances\n");
			compute(nlhs, plhs, nrhs, prhs, 1, dummy);
		}
		else if ((!strncmp("maximum", metric, strlen(metric)))) {
			maximum_distance dummy;
			// mexPrintf("Using maximum metric to calculated distances\n");
			compute(nlhs, plhs, nrhs, prhs, 1, dummy);
		}
		else
			printf("ATRIA preprocessing structure was not created using a supported metric; doing preprocessing again\n");

		mxFree(metric);
#endif
	} else {
		euclidian_distance dummy;
		// mexPrintf("Using euclidian metric to calculated distances\n");
		compute(nlhs, plhs, nrhs, prhs, 0, dummy);
	}


}
