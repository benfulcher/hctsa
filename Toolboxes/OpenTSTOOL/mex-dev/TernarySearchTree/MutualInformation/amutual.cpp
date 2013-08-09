// Fast, but crude auto mutual information for a scalar
// timeseries for different timelags (0 .. maxtau)
//
// Input time series should be much longer than maximal
// timelag
// 
// Uses equidistant histogram boxes, so results are bad
// in a mathmatical sense
//
// A fast algorithm based on ternary search trees to store nonempty boxes is used
// C.Merkwirth,J.Wichard DPI Goettingen 1998 
//
// mex -I.. -I../.. amutual.cpp -O

#include <cmath>

#include "mextools/mextools.h"
#include "ternary_search_tree.h"

#define PARTITIONMAX 16384

using namespace std;

#define LOGARITHM_OF_2 0.69314718055995
inline double log_dualis(const double x) { return (log(x)/LOGARITHM_OF_2); }

inline double max(double a,double b) {return ((a >= b) ? a : b);}

class mt_ternary_tree : public ternary_search_tree<int>
{
	private:
	
	public:
		mt_ternary_tree() : ternary_search_tree<int>(2) {};
		~mt_ternary_tree() {};

		int insert2(const int* const key);
		
		template<class Evaluater>
		void traverse2(Evaluater& eval);
};

// return 0 on SUCCESS
int mt_ternary_tree::insert2(const int* const key)
{   
	int d;
	Tptr pp;
	
    Tptr* p = &root;
	long level = 0;			// level goes up to len-1 
	
    while(pp = *p) {		// as long as we encounter already exisiting nodes, we stay inside this while loop
		if ((d = key[level] - pp->splitkey) == 0) {		// go to next tree level
			pp->count++;
            p = &(pp->eqkid);
			if ((++level) == 2) return 0;
        } else if (d < 0) {
            p = &(pp->lokid);	/* move left in the current level */
		}
        else {
            p = &(pp->hikid);	/* move right in the current level */
		}
    }
    for (;;) {	/* once we find a node that is not allocated (==0), we must create every next node */
		if (freen-- == 0) {
			if (bufn == TST_BUFFERS) {
				//mexErrMsgTxt("Ran out of available buffers for tree nodes");
				return -1;		// FAILURE
			}			
			buf = new Tnode<int>[next_buf_size]; 
			freearr[bufn++] = buf;
			freen = next_buf_size-1;
			next_buf_size *= 2; 	// double size of the next buffer (this keeps overall number of allocations small)
		}
		*p = buf++;
        pp = *p;
        pp->splitkey = key[level];
		pp->count = 1;					// this node is newly created, so count is set to one 
		pp->level = level;
        pp->lokid = pp->eqkid = pp->hikid = 0;
        if ((++level) == 2) {
			pp->eqkid = (Tnode<int> *) key[0]; 	// dirty trick, for nodes of level 1 use eqkid as storage for key of level 0
			return 0;
		}	
        p = &(pp->eqkid);
    }
}

// traverse tree, execute the given function object on every node that is not empty
// used for mutual information computation
template<class Evaluater>
void mt_ternary_tree::traverse2(Evaluater& eval)
{   
	long buf_size = TST_BUFSIZE;		// the actual size of the buffer is doubling each iteration
	
	// traverse through all buffers that are completely filled 
	for (long i = 0; i < bufn; i++) {
		long number_nodes = buf_size;
		const Tptr b = (Tptr) freearr[i];
		
		if (i == bufn - 1)		// last buffer may not be completely filled
			number_nodes = buf_size-freen;
		
		buf_size *= 2;
		
		for (long j = 0; j < number_nodes; j++) {
			const Tptr p = b + j;	
			
			if (p->level) {		// ==> if (p->level == 1)
				eval(p->count, p->splitkey, (unsigned long) (p->eqkid));
			}
		}
	}
}

class mutinfo_estimator {
	protected:
		double mut_info;
		
		const double* marginal_histogram;
		
		const long N;			// total number of points			

	public:
		mutinfo_estimator(const long n, const double* Marginal_histogram) : mut_info(0), 
			marginal_histogram(Marginal_histogram), N(n)  {};
		~mutinfo_estimator() {};
		
		void operator()(const long mass, const int key, const int parent_key) {
			// mass is the absolute frequency of points falling into that box at level level of the tree
			
			const double pAB = (((double)mass) / (double)N);		// relative frequency

			const double pA = marginal_histogram[key];
			const double pB = marginal_histogram[parent_key];

		    mut_info += pAB * log_dualis(pAB/(pA * pB));
		}

		double result() const { return mut_info; }
};


void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{		
	/* check input args */
	if (nrhs < 3)
	{
		mexErrMsgTxt("amutual : Scalar time series, maximal timelag and number of partitions per axis must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long N = (long) max(mxGetM(prhs[0]), mxGetN(prhs[0]));	
	const double* p  = (double *)mxGetPr(prhs[0]);
 
    const long maxtau = (long) *((double *)mxGetPr(prhs[1]));

	const double* partitionsizes  = (double *)mxGetPr(prhs[2]);	// vector of partition sizes
	const long R = (long) max(mxGetM(prhs[2]), mxGetN(prhs[2]));		
	
	if (N-maxtau < 1) {
		mexErrMsgTxt("Length of time series much to short for given timelag");
		return;
	}	
	if (maxtau < 0) {
		mexErrMsgTxt("Maximal timelag may not be negative");
		return;
	}
	if (R < 1) {
		mexErrMsgTxt("At least one partition size must be given");
		return;
	}
	
	for (long r=0; r < R; r++) {	// check partition sizes
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

	plhs[0] = mxCreateDoubleMatrix(maxtau+1, R, mxREAL);
	double* out = (double *) mxGetPr(plhs[0]);		// scaling of box counting dimension
		
	int* ts = new int[N];  
		
	// start of main loop over partition sizes
	for (long r=0; r < R; r++) { 		// we have R different partition sizes
          const long partitions = (long) partitionsizes[r];
		
        // mexPrintf("Number of partitions per axis   : %d\n", partitions); // +BF commented out
		
		double* marginal_histogram = new double[partitions];
			
		for (long i=0; i < partitions; i++) marginal_histogram[i] = 0;	
			
		for (long index=0; index < N; index++) {
                  const long key = (long) floor(partitions * p[index]);
			
			if ((key < 0) || (key>=partitions)) 
				mexErrMsgTxt("Data values must be rang values in the range [0 1)");
				
			ts[index] = key;
			marginal_histogram[key]++;
		}

		for (long i=0; i < partitions; i++) marginal_histogram[i] /= N;

		for (long tau=0;  tau <= maxtau; tau++) {
			int key_vector[2];
			
			mutinfo_estimator estimator(N-maxtau, marginal_histogram);
			
			mt_ternary_tree tree;

			// Tree fuellen
			for (long index=0; index < N-maxtau; index++) {
				key_vector[0] = ts[index];
				key_vector[1] = ts[index+tau];
				tree.insert2(key_vector);
			}

			tree.traverse2(estimator);
				
			*(out++) = estimator.result();
		}
		
		delete[] marginal_histogram;
	}
	
	delete[] ts;
}	

