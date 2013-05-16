// Compute histogram of return (recurrence) times

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "include.mex"

static double max(double a,double b) {return ((a >= b) ? a : b);}
static double min(double a,double b) {return ((a <= b) ? a : b);}
static inline long longabs(const long a) { return ((a<0) ? -a : a); }

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
		mexErrMsgTxt("Return time : Data set of points (row vectors), refernce indices, number of nearest neighbours, maxT and past must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N 		= mxGetM(prhs[0]);	
	const long dim  	= mxGetN(prhs[0]);
	const double* p 	= (double *)mxGetPr(prhs[0]);
	
	const long R = max(mxGetM(prhs[1]), mxGetN(prhs[1]));
	const double* ref 	= (double *)mxGetPr(prhs[1]);

	const long NNR 	= (long) *((double *)mxGetPr(prhs[2]));
	
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

	long maxT = (long) *((double *)mxGetPr(prhs[3])); 	// maximal return time, higher return times are neglected
	
	if (maxT < 2) {
		mexErrMsgTxt("maxT must be greater two");
		return;
	}	
	if (maxT >= N) maxT = N-1;	// this is the highest temporal distance that can be found in data set of length N

	const long past	= (long) *((double *)mxGetPr(prhs[4]));
	
	if (past < 0) {
		mexErrMsgTxt("past may not be negative");
		return;
	}	

	if (R < 1) {
		mexErrMsgTxt("At least one reference index or point must be given");
		return;
	}	

	for (long i=0; i < R; i++) {
		if ((ref[i] < 1) || (ref[i]>N)) {
			mexErrMsgTxt("Reference indices out of range");
			return;
		}	
	}

	//mexPrintf("%d %d %d %d %d\n", N, dim, NNR, maxT, past);	
						
#include "create_searcher.cpp"

	double*	const coord = (double*) malloc(dim * sizeof(double)); 
	
	plhs[0] = mxCreateDoubleMatrix(maxT, 1, mxREAL);
	double* hist = (double *) mxGetPr(plhs[0]);
	
	for (long n=0; n < maxT; n++) hist[n] = 0;	
	
	unsigned long counter = 0;
	
	for (long i=0; i < R; i++) {
		const long actual = (long) ref[i]-1;
		vector<neighbor> v;
		
		for (long k=0; k < dim; k++) coord[k] = p[actual+k*N];
		
		searcher->search_k_neighbors(v, NNR, coord, actual-past, actual+past);
	
		for (long k = 0; k < v.size(); k++) { 	// v is the sorted vector of neighbors
			const long T = longabs(v[k].index() - actual);	
			
			if (T == 0) {
				mexErrMsgTxt("Internal error : Return time of zero found");
				return;
			}	
			
			// hist[0] counts events T==1, hist[1] is for events T==2 etc.			
			if (T <= maxT) hist[T-1]++; 		
		}
	}
	
	// normalize values
	for (long n=1; n <= maxT; n++) hist[n-1] /= (2*NNR*(N-n));	
	
	free(coord);
	
	delete searcher;
}	



