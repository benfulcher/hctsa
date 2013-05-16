// crossprediction
// Search for functional dependencies between two data sets
//
//
// mex -I. -I.. -I./STL -I./NNSearcher crossprediction.cpp -O
//
// Origin : cmerk 1998
// Revision : cmerk Nov. 1999

#include "mex.h"

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

static double max(double a,double b) {return ((a >= b) ? a : b);}
static double min(double a,double b) {return ((a <= b) ? a : b);}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	int atria_preprocessing_given = 0;		// flag wheter preprocessing is already given on command line

	ATRIA< point_set<euclidian_distance> >* searcher = 0;

	// try to see if the first parameter given is an atria structure
	// If this is true, the order of input parameters is shifted by one
	if ((nrhs) && (mxIsStruct(prhs[0]))) {
		atria_preprocessing_given = 1;
		prhs++; 	// these two lines enable us to use the old argument parsing block without changing it
		nrhs--;
	}

	if (nrhs < 6)
	{
		mexErrMsgTxt("crossprediction : data set A of points (row vectors), data set B of points (row vectors), reference indices, taumax, number of NN and past must be given\n");
		return;
	}

	// Data set A is the data set in which the neighbor are searched
	// All reference indices must point into A

	const long N		= mxGetM(prhs[0]);
	const long dim  	= mxGetN(prhs[0]);
	const double* p  	= (double *)mxGetPr(prhs[0]);

	// Data set B is the data set in which the crossprediction coefficient is computed
	// Data set B should be at least of size A + taumax
	// The dimension of A and B need not be equal

	const long N2		= mxGetM(prhs[1]);
	const long dim2  	= mxGetN(prhs[1]);
	const double* p2  	= (double *)mxGetPr(prhs[1]);

	const double* ref 	= (double *)mxGetPr(prhs[2]);		// pointer to reference indices
	const long R = max(mxGetM(prhs[2]), mxGetN(prhs[2]));	// number of reference indices``

	const long taumax 	= (long) *((double *)mxGetPr(prhs[3]));
	const long NNR  	= (long) *((double *)mxGetPr(prhs[4]));
	const long past 	= (long) *((double *)mxGetPr(prhs[5]));

	if (N < 1) {
		mexErrMsgTxt("Data set A must contain at least two points (row vectors)");
		return;
	}
	if (taumax < 0) {
		mexErrMsgTxt("taumax may not be negative");
		return;
	}
	if (N2 < N + taumax) {
		mexErrMsgTxt("Data set B must contain at least as much points as data set A plus taumax");
		return;
	}
	if (NNR < 1) {
		mexErrMsgTxt("Number of nearest neighbors must be greater zero");
		return;
	}

	mexPrintf("Number of reference points : %d\n", R);

	if (R < 2) {
		mexErrMsgTxt("At least two reference indices must be given");
		return;
	}

	point_set<euclidian_distance> points(N, dim, p);
	point_set<euclidian_distance> points2(N2, dim2, p2);

	double average_distance = 0;		// "delta", compute an average distance in data set B
	long counter = 0;

	for (long i=0; i < R; i++) { 	// check reference indices
		if ((ref[i] < 1) || (ref[i]>N) ) {
			mexErrMsgTxt("Reference indices out of range");
			return;
		}
		average_distance += points2.distance((long) 0, (long) ref[i]-1);	// compute distances from reference points to point nr. 0
		counter++;
		if (i>0) {	// and compute distances among reference indices (which are hopefully randomly choosen) */
			average_distance += points2.distance((long) ref[i-1]-1, (long) ref[i]-1);
			average_distance += points2.distance((long) ref[0]-1, (long) ref[i]-1);
			counter+=2;
		}
	}

	average_distance /= counter;

	if (average_distance <= 0) {
		mexErrMsgTxt("Average interpoint-distance in data set B seems to be zero");
		return;
	}

	// Create searcher object
	// If the object is given as input argument, use it
	// Otherwise, create searcher by doing the preprocessing
	if (atria_preprocessing_given) {
		char* metric = 0;

		if (mxIsChar(mxGetField(prhs[-1], 0, "optional"))) {
    		long buflen = (mxGetM(mxGetField(prhs[-1], 0, "optional")) * mxGetN(mxGetField(prhs[-1], 0, "optional"))) + 1;
 			metric = (char*) mxMalloc(buflen);
        	mxGetString(mxGetField(prhs[-1], 0, "optional"), metric, buflen);
		}

		if ((metric == 0) || (!strncmp("euclidian", metric, strlen(metric)))) {
			searcher = new ATRIA< point_set<euclidian_distance> >(points, prhs[-1]);	// this constructor used the data from the preprocessing
			if (searcher->geterr()) {
				delete searcher;
				searcher = 0;
			}
		} else
			printf(" ATRIA preprocessing structure was not created using Euclidian metric; doing preprocessing again\n");

		mxFree(metric);
	}

	if (searcher == 0) {
		searcher = new ATRIA< point_set<euclidian_distance> >(points);
	}

	if (searcher->geterr()) {
		mexErrMsgTxt("Error preparing searcher");
		return;
	}
	// End of creating searcher object

	double*	const coord = (double*) malloc(dim * sizeof(double));

	plhs[0] = mxCreateDoubleMatrix(taumax+1, 1, mxREAL);
	double* x = (double *) mxGetPr(plhs[0]);	// uncorrected output values

	double* y;									// corrected output values

	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(taumax+1, 1, mxREAL);
		y = (double *) mxGetPr(plhs[1]);
	} else
		y = (double *) malloc(taumax+1 * sizeof(double));

	// initialize all variables

	for (long tau = 0; tau <= taumax; tau++) { x[tau] = y[tau] = 0; }

	long count_uncorrected = 0; 	// needed to normalize x
	double total_weight = 0;		// needed to normalize y

	for (long n=0; n < R; n++) { 		/* main loop, iterate over all reference points */
		vector<neighbor> v;

        // actual := actual reference point index
        const long actual = (long) ref[n]-1;       //  Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...

		for (long k=0; k < dim; k++) coord[k] = points.coordinate(actual,k);

		searcher->search_k_neighbors(v, NNR, coord, actual-past, actual+past);

		if (v.size() > 0) {				// v is the sorted vector of neighbors
			vector<neighbor>::iterator i;

			double sumA = 0;			// is used as a weight for the sum of distances in B
			double sumB = 0;

			count_uncorrected ++;

			// compute sums of distances for tau = 0
			for (i = v.begin(); i < v.end(); i++) {
				sumA += fabs((*i).dist());
				sumB += fabs(points2.distance((long) actual, (long) (*i).index()));
			}

			x[0] += sumB/v.size();

			if (sumA > 0) {
				total_weight += 1.0/sumA;
				y[0] += (sumB/sumA/v.size());
			}

			for (long tau = 1; tau <= taumax; tau++) {
				double sumB = 0;

				for (i = v.begin(); i < v.end(); i++) {
					sumB += fabs(points2.distance((long) actual+tau, (long) (*i).index()+tau));
				}

				x[tau] += sumB/v.size();
				if (sumA > 0) y[tau] += (sumB/sumA/v.size());
			}
		}
	}

	if (count_uncorrected*average_distance > 0)
		for (long tau = 0; tau <= taumax; tau++) x[tau] /= ((double) count_uncorrected * average_distance);

	if (average_distance*total_weight > 0)
		for (long tau = 0; tau <= taumax; tau++) y[tau] /= (total_weight * average_distance);


	free(coord);
	if (!(nlhs > 1)) free(y);

	delete searcher;
}

