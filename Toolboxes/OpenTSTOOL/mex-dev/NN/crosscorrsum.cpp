// Compute cross correlation sum versus radius
// Reference points are out of the second data set (data set B)
// Neighbors are searched in the firstdata set (data set A)
// cmerk DPI Goettingen 2000

// mex -I. -I.. crosscorrsum.cpp -O

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

int be_verbose=0; 	// used to suppress any kind of diagnostic output

inline long lmax(const long a, const long b) { return a >= b ? a : b; }
inline long lmin(const long a, const long b) { return a <= b ? a : b; }

template<class METRIC>
void compute(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[], const int atria_preprocessing_given,  METRIC dummy)
{
	long bins = 32;
	double start_dist = 0;
	double maximal_search_radius = 0;
	double scale_factor = 1;
	long opt_flag = 0;	// 0 => eucl.norm, 1 => max.norm, utm, 2 => eucl, silent 3 => max., silent

	/* handle matrix I/O */

	const long N = mxGetM(prhs[0]);
	const long dim = mxGetN(prhs[0]);
	const double* p  	= (double *)mxGetPr(prhs[0]);

	const long R = mxGetM(prhs[1]);
	const double* ref 	= (double *)mxGetPr(prhs[1]);

	const double relative_range = (double) *((double *)mxGetPr(prhs[2]));

	if (nrhs > 3) bins = (long) *((double *)mxGetPr(prhs[3]));

	if (nrhs > 4) opt_flag = (long) *((double *)mxGetPr(prhs[4]));

	if (N < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}
	if (dim < 1) {
		mexErrMsgTxt("Data points must be at least of dimension one");
		return;
	}
	if (relative_range <= 0) {
		mexErrMsgTxt("Relative range must be greater zero");
		return;
	}
	if ((mxGetN(prhs[1]) == 0) || (mxGetM(prhs[1]) == 0)) {
		mexErrMsgTxt("Wrong reference indices given");
		return;
	}
	if (bins < 2) {
		mexErrMsgTxt("Number of bins should be two");
		return;
	}
	if (R < 1) {
		mexErrMsgTxt("At least one reference point must be given");
		return;
	}
	if (mxGetN(prhs[1]) != dim) {
		mexErrMsgTxt("Dimension of reference points and data set points is not the same");
		return;
	}
	if ((opt_flag < 0) || (opt_flag > 3)) {
		mexErrMsgTxt("Flag must be out of 0...3");
		return;
	}

	point_set<METRIC> points(N,dim, p);
	point_set<METRIC> ref_points(R, dim, ref);

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

	if (mxGetM(prhs[2])*mxGetN(prhs[2]) == 2) {
		start_dist = (double) ((double *)mxGetPr(prhs[2]))[0];
		maximal_search_radius = (double) ((double *)mxGetPr(prhs[2]))[1];

		if (start_dist <= 0) {
			mexErrMsgTxt("Starting radius is zero or negativ");
			return;
		}
		if (maximal_search_radius <= start_dist) {
			mexErrMsgTxt("Maximal search radius must be greater than starting radius");
			return;
		}
	}
	else {
		maximal_search_radius = relative_range * searcher->data_set_radius();  // compute the maximal search radius using information about attractor size

		// try to determine an estimate for the minimum inter-point distance in the data set, but greater zero

		double*	const coord = new double[dim];

		for (long n=0; n < 128; n++) {
			const long actual = (long) (((double)rand() * (double) (N-2))/(double)RAND_MAX);
			vector<neighbor> v;

			for (long k=0; k < dim; k++) coord[k] = points.coordinate(actual,k);

			searcher->search_k_neighbors(v, 1, coord, actual, actual);	// search the nearest neighbor

			if (v.size() > 0) {
	    		if (start_dist == 0) start_dist = v[0].dist();

				if (v[0].dist() > 0) {
					if (v[0].dist() < start_dist) start_dist = v[0].dist();
				}
			}
		}

		if (start_dist <= 0) {	// first try to search again for a minimum inter-point distance greater zero
			for (long n=0; n < 512; n++) {
				const long actual = (long) (((double)rand() * (double) (N-2))/(double)RAND_MAX);
				vector<neighbor> v;

				for (long k=0; k < dim; k++) coord[k] = points.coordinate(actual,k);
				searcher->search_k_neighbors(v, 1, coord, actual, actual);	// search the nearest neighbor

				if (v.size() > 0) {
		    		if (start_dist == 0) start_dist = v[0].dist();

					if (v[0].dist() > 0) {
						if (v[0].dist() < start_dist) start_dist = v[0].dist();
					}
				}
			}
		}

		delete[] coord;

		if (start_dist <= 0) {	// give up if we cannot find an interpoint distance greater zero
			mexErrMsgTxt("Cannot find an interpoint distance greater zero, maybe ill-conditioned data set given");
			return;
		}
		if (maximal_search_radius <= start_dist) {
			mexErrMsgTxt("Maximal search radius must be greater than starting radius");
			return;
		}
	}

	scale_factor = pow(maximal_search_radius/start_dist, 1.0/(bins-1));

	if (be_verbose)	{
		mexPrintf("Number of reference points            : %d\n", R);
		mexPrintf("Upper bound for attractor size        : %f\n", 2 * searcher->data_set_radius());
		mexPrintf("Starting distance                     : %f\n", start_dist);
		mexPrintf("Maximal search radius                 : %f\n", maximal_search_radius);
		mexPrintf("Number of partitions used             : %d\n\n", bins);
	}

	plhs[0] = mxCreateDoubleMatrix(bins, 1, mxREAL);
	double* corrsums = (double *) mxGetPr(plhs[0]);
	double* dists;

	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(bins, 1, mxREAL);
		dists = (double *) mxGetPr(plhs[1]);
	} else
		dists = new double[bins];

	double x = start_dist;

	// initialize vectors
	// corrsum[i] counts the number of points/distances (real) smaller than dists[i]
	for  (long bin=0; bin < bins; bin++) {
		corrsums[bin] = 0;
		dists[bin] = x;
		x *= scale_factor;
	}

	long total_pairs = 0;	// count the total number of points pairs that were at least theoretically tested

	for (long n=0; n < R; n++) { 		/* iterate over all reference points */
		vector<neighbor> v;
		const long pairs = N;

		if (pairs <= 0)
			continue;

		total_pairs += pairs;

		searcher->search_range(v, maximal_search_radius, ref_points.point_begin(n), -1, -1);

		if (v.size() > 0) {
			for (vector<neighbor>::iterator i = v.begin(); i < v.end(); i++) { // v is unsorted (!!!)
				const double dist = (*i).dist();

				int mi, ni; 	// index of first and last (inclusive) element of array dists

				if (dist >= dists[ni = bins-1]) 	// => we know that dist <= dists[ni]
					continue;

				if (dist < dists[mi = 0]) { 		// => we know that dists[mi] < dist
					corrsums[0]++;
					continue;
				}

				do
                {
                    const int ki = (mi + ni) / 2;	// "mi < ki < ni" if "mi + 1 < ni"
                    if (dist < dists[ki])
						ni = ki;    	// "dists[mi] < dist <= dists[ni]" is still valid
					else
                        mi = ki;    	// "dists[mi] < dist <= dists[ni]" is still valid
                } while (ni > mi+1);

				// now, "ni == mi+1" and "dists[mi] < dist <= dists[ni]" is true,
				// so "dists[mi] < dist <= dists[mi+1]", so mi+1 is the desired bin where
				// to count this distance

				corrsums[ni]++;
			}
		}
	}

	if (total_pairs > 0) {
		unsigned long sum = 0;
		for  (long bin=0; bin < bins; bin++) {
			sum += (long) corrsums[bin];
			corrsums[bin] = ((double) sum) / (double) total_pairs;
		}
	}

	if (nlhs > 2) {
		plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
		*((double *) mxGetPr(plhs[2])) = (double) total_pairs;
	}

	if (!(nlhs > 1)) delete[] dists;

	delete searcher;
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	int atria_preprocessing_given = 0;		// flag wheter preprocessing is already given on command line
	long opt_flag = 0;	// 0 => eucl.norm, upper triangle matrix, 1 => max.norm, utm, 2 => eucl,full matrix, 3 => max.,full matrix

	be_verbose = 0; // Don't spit out feedback to user by default

	// try to see if the first parameter given is an atria structure
	// If this is true, the order of input parameters is shifted by one
	if ((nrhs > 0) && mxIsStruct(prhs[0])) {
		atria_preprocessing_given = 1;
		prhs++; 	// these two lines enable us to use the old argument parsing block without changing it
		nrhs--;
	}

	/* check input args */
	if (nrhs < 3)
	{
		mexErrMsgTxt("Correlation sum : Data set of points (row vectors), reference points and relative range (relative to attractor diameter) must be given, number of bins is optional");
		return;
	}

	if (nrhs > 4) opt_flag = (long) *((double *)mxGetPr(prhs[4]));

	if (opt_flag & (long) 2)
		be_verbose = 0;

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

			if (be_verbose)
				mexPrintf("Using euclidian metric to calculated distances\n");

			compute(nlhs, plhs, nrhs, prhs, 1, dummy);
		}
		else if ((!strncmp("maximum", metric, strlen(metric)))) {
			maximum_distance dummy;

			if (be_verbose)
				mexPrintf("Using maximum metric to calculated distances\n");

			compute(nlhs, plhs, nrhs, prhs, 1, dummy);
		}
		else
			printf("ATRIA preprocessing structure was not created using a supported metric; doing preprocessing again\n");

		mxFree(metric);
#endif
	} else {
		if (opt_flag & (long)1) {
			maximum_distance dummy;

			if (be_verbose)
				mexPrintf("Using maximum metric to calculated distances\n");
			compute(nlhs, plhs, nrhs, prhs, 0, dummy);
		} else {
			euclidian_distance dummy;
			if (be_verbose)
				mexPrintf("Using euclidian metric to calculated distances\n");

			compute(nlhs, plhs, nrhs, prhs, 0, dummy);
		}
	}
}
