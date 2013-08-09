// Box counting for a data set of points (row vectors) for
// different partition sizes (2, 3, 4, ...), (equidistant partitioning)
//
// Every value of input point set must be scaled to be within [0,1]
// Returns boxcounting, information and correlation dimension (D0, D1 and D2) (using log2)
// 
// A fast algorithm based on ternary search trees to store nonempty boxes is used
// C.Merkwirth,J.Wichard DPI Goettingen 1998 

// mex -I.. -I../.. boxcount.cpp -O

#include <cmath>

#include "mextools/mextools.h"
#include "ternary_search_tree.h"

#define PARTITIONMAX 16384
#define MY_EPS 1e-6
#define LOGARITHM_OF_2 0.69314718055995
#define MAXDIM 512

using namespace std;

inline double log_dualis(const double x) { return (log(x)/LOGARITHM_OF_2); }

inline double max(double a,double b) {return ((a >= b) ? a : b);}
inline double min(double a,double b) {return ((a <= b) ? a : b);}

class dim_estimator {
	protected:
		long* boxes; 	/* scaling of number of boxes for dimensions from 1 to dim */
		double* in;		/* scaling of information for dimensions from 1 to dim */
		double* co;		/* scaling of correlation for dimensions from 1 to dim */
		
		const long dim; 
		const long N;		
		long total_points;
	
	public:
		dim_estimator(const long n, const long Dim) : N(n), dim(Dim), total_points(0),
			boxes(0), in(0), co(0) { 
			boxes = new long[dim];
			in = new double[dim];
			co = new double[dim];

			for (long d=0; d < dim; d++) { 
				boxes[d] = 0;
				in[d] = co[d] = 0; 
			}
		}
		~dim_estimator() {
			delete[] co;
			delete[] in;
			delete[] boxes;
		};
		
		void operator()(const long mass, const int level) {
			// mass is the absolute frequency of points falling into that box at level level of the tree
			
			const double p_i = (((double)mass) / (double)N);		// relative frequency
	
		    total_points += mass;  			    		// Check we got all points
		    boxes[level]++;					    		// D0 = Boxen zaehlen (capacity dimension)
		    in[level] += (p_i * log_dualis(p_i));  		// D1 information dimension
		    co[level] += p_i * p_i;			    		// D2 correlation dimension				
		}

		long get_total_points() const { return total_points; }		// mainly used to check all points were inserted correctly
		double boxd(const long d) const { return  -log_dualis((double) boxes[d]); }
		double infod(const long d) const { return  in[d]; }
		double corrd(const long d) const { return  log_dualis(co[d]); }
};


void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	long r;
			
	/* check input args */
	if (nrhs < 2)
	{
		mexErrMsgTxt("Boxcount : Data set of points (row vectors) and number of partitions (per axis) must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N = mxGetM(prhs[0]);	
	const long dim = mxGetN(prhs[0]);
	const double* p  = (double *)mxGetPr(prhs[0]);

	const double* partitionsizes  = (double *)mxGetPr(prhs[1]);	// vector of partition sizes
	const long R = (long) max(mxGetM(prhs[1]), mxGetN(prhs[1]));		
	
	if (N < 1) {
		mexErrMsgTxt("Data set must consist of at least two points (row vectors)");
		return;
	}		
	if (dim < 1) {
		mexErrMsgTxt("Data points must be at least of dimension one");
		return;
	}	
	if (dim > MAXDIM) {
		mexErrMsgTxt("Maximal dimension exceeded");
		return;
	}	
	if (R < 1) {
		mexErrMsgTxt("At least one partition size must be given");
		return;
	}
	
	for (r=0; r < R; r++) {	// check partition sizes
          const long partitions = (long) partitionsizes[r];
		if (partitions < 2) {
			mexErrMsgTxt("Number of bins should be at least greater one");
			return;
		}	
		if (partitions > PARTITIONMAX) {
			mexErrMsgTxt("Number of bins too large");
			return;
		}
	}
	
	plhs[0] = mxCreateDoubleMatrix(R, dim, mxREAL);
	double* out = (double *) mxGetPr(plhs[0]);		// scaling of box counting dimension

	double* info;		//  scaling of information dimension
	double* corr;		//  scaling of correlation dimension
	
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(R, dim, mxREAL);
		info = (double *) mxGetPr(plhs[1]);
	} else
		info = (double *) malloc(R * dim * sizeof(double));

	if (nlhs > 2) {
		plhs[2] = mxCreateDoubleMatrix(R, dim, mxREAL);
		corr = (double *) mxGetPr(plhs[2]);
	} else
		corr = (double *) malloc(R * dim * sizeof(double));
	
	double minimum = p[0];
	double maximum = p[0];
		
	for (long i=1; i < N * dim; i++) {	
		if (minimum > p[i]) 
			minimum = p[i];
		if (maximum < p[i]) 
			maximum = p[i];	
	}

#ifdef VERBOSE				
	mexPrintf("Minimum and maximum of input data  : %lf and %lf \n", minimum, maximum);		
#endif
	
	if (minimum == maximum) {
		mexPrintf("Data values are indentical, nothing to compute\n");			
		for (long r=0; r < R; r++) { 		
			for (long d=0; d < dim; d++) {		
				out[r + d*R]  = 0;
				info[r + d*R] = 0;
				corr[r + d*R] = 0;
			}
		}
	} else {		// start of main loop over partition sizes
		for (long r=0; r < R; r++) { 							// we have R different partition sizes
			int key[MAXDIM];

			dim_estimator estimator(N, dim); 
			ternary_search_tree<int, dim_estimator> tree(dim);

			const long partitions = (long) partitionsizes[r];
			const double a = (((double) partitions)-MY_EPS) / (maximum - minimum);
			const double b = a * minimum;

#ifdef VERBOSE			
            // mexPrintf("Number of partitions per axis   : %d\n", partitions);
#endif
			// Tree fuellen
			for (long index=0; index < N; index++) {
				for (long d=0; d < dim; d++) {
					//const int x = floor((partitions-MY_EPS) * (p[index + d*N] - minimum) / (maximum- minimum));
                                  const int x = (int) floor(a*p[index + d*N] - b);
					
					if ((x < 0) || (x >= partitions)) {
						mexErrMsgTxt("Data values seem not to be in range [0,1]");
						return;
					}

					key[d] = x;	
				}
				if (tree.insert(key)) {
					mexErrMsgTxt("Ran out of memory");
					return;
				}
			}

			// Tree auszaehlen, dabei werden alle relevanten Groessen berechnet

			tree.traverse(estimator);

#ifdef VERBOSE
			mexPrintf("Total nodes allocated   : %d\n", tree.total_nodes());
			mexPrintf("Total memory allocated   : %d bytes\n", tree.total_nodes() * sizeof(Tnode<int>));
#endif

			if (estimator.get_total_points() != dim * N) {
				mexErrMsgTxt("Internal error, tree did not contain all points");
				return;
			}

			for (long d=0; d < dim; d++) {		
				out[r + d*R]  = estimator.boxd(d);
				info[r + d*R] = estimator.infod(d);
				corr[r + d*R] = estimator.corrd(d);
			}
		}
	}

	if (!(nlhs > 1)) free(info);
	if (!(nlhs > 2)) free(corr);
}	

