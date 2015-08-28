// estimate minimum embedding dimension, using Cao's method

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

static double max(double a,double b) {return ((a >= b) ? a : b);}
static double min(double a,double b) {return ((a <= b) ? a : b);}
inline double squared(const double x) { return x*x; }

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	long i,j,n,k;		/* loop variables */

	/* check input args */
	if (nrhs < 3)
	{
		mexErrMsgTxt("Minimum embedding dimension : Data set of points (row vectors), reference indices and maximal number of neighbours must be given");
		return;
	}

	/* handle matrix I/O */

	const long N 		= mxGetM(prhs[0]); 		// number of points
	const long dim  	=  mxGetN(prhs[0]);		// maximal dimension
	const double* p 	= (double *)mxGetPr(prhs[0]);

	const double* ref 	= (double *)mxGetPr(prhs[1]);
	const long R 		= max(mxGetM(prhs[1]), mxGetN(prhs[1]));	// number of query (reference) points

	const long NNR  	= (long) *((double *)mxGetPr(prhs[2]));

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
	if ((NNR < 1) || (NNR > N)) {
		mexErrMsgTxt("Number of nearest neighbours must be within 1 and number of points in the data set");
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

	double*	const coord = (double*) malloc(dim * sizeof(double));

	// mexPrintf("Number of reference points            : %d\n", R);
	// mexPrintf("Maximal number of nearest neighbours  : %d\n", NNR);	

	plhs[0] = mxCreateDoubleMatrix(dim-1, 1, mxREAL);
	double* const Eout = (double *) mxGetPr(plhs[0]);

	double* Estar = 0;
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(dim-1, 1, mxREAL);
		Estar =	(double *) mxGetPr(plhs[1]);
	}
	else
		Estar = (double*) malloc((dim-1) * sizeof(double));

	for (long d=1; d < dim; d++) {
		double E = 0;
		double Est = 0;

		long counter = 0;

		point_set<euclidian_distance> points(N,d, p);	// this is possible because a matrix of size N by dim can be used as submatrix of size N by 1...dim
		ATRIA< point_set<euclidian_distance> > searcher(points);

		if (searcher.geterr()) {
			mexErrMsgTxt("Error preparing searcher");
			return;
		}

		for (n=0; n < R; n++) { 				/* iterate over all reference points */
			vector<neighbor> v;
			vector<neighbor>::iterator i;

			double x = 0;
			double y = 0;
			double z = 0;

			const long actual = (long) ref[n]-1;		/* Matlab to C means indices change from 1 to 0, 2 to 1, 3 to 2 ...*/

			for (k=0; k < d; k++) coord[k] = p[actual+k*N];

			searcher.search_k_neighbors(v, NNR, coord, actual, actual);

			for (i = v.begin(); i < v.end(); i++) { 	// v is the sorted vector of neighbors
				const long index = (*i).index();
				const double dist = (*i).dist();
				const double diff = p[index + d*N] - p[actual + d*N];

				y += dist;
				x += sqrt(squared(dist) + squared(diff));
				z += fabs(diff);
			}

			if (y > 0) {
				E += x/y;
				Est += z;
				counter++;
			}
		}

		Eout[d-1] = E/counter;
		Estar[d-1] = Est/counter;
	}

	free(coord);

	if (!(nlhs > 1)) free(Estar);

}
