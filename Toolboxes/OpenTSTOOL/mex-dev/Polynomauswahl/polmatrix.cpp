// Compute list of monoms of degree maxgrad or less
// mex -I.. polmatrix.cpp -O

#include <cmath>
#include "mextools/mextools.h"

template<class T>
class interleaved_pointer	
{
	private:
		const T* ptr;
		const long increment;	
	public:
		typedef interleaved_pointer self;
		inline interleaved_pointer(const T* const p, const long inc) : ptr(p), increment(inc) {};
		inline T operator*() const { return *ptr; }
		inline self& operator++() {
    		ptr+=increment;
    		return *this;
  		}
  		inline self& operator--() {
    		ptr-=increment;
    		return *this;
  		}		
		bool operator==(const interleaved_pointer& x) {
			return ptr == x.ptr; 
		}
		bool operator!=(const interleaved_pointer& x) {
			return ptr != x.ptr; 
		}		
};	

class embedded_time_series_point_set {		
	protected:
		const long N;
		const long D;		// dimension
		const long DELAY;
		
		const double* const ts_ptr; 	
	public:
		typedef interleaved_pointer<double> point_iterator; // a smart pointer that iterates over the elements of one point in the embedded_time_series_point_set (points are row vectors)

		// n - Length of scalar time-series
		// dim - Embedding dimension
		// delay - Time delay in samples
		// ts - zero based vector containing time series data
		embedded_time_series_point_set(const long n, const long dim, const long delay, const double* ts) 
			: N(n-(dim-1)*delay), D(dim), DELAY(delay), ts_ptr(ts) {};
			
		~embedded_time_series_point_set() {};
		
		inline long dimension() const { return D; }; 
		inline long size() const { return N; };
				
		point_iterator point_begin(const long i) const { return point_iterator(ts_ptr + i + (D-1)*DELAY, -DELAY); }
		point_iterator point_end(const long i) const { return point_iterator(ts_ptr + i - DELAY, -DELAY); }	// past-the-end		
};

static long maxlen;
static long* a;
static long counter;
static int count_flag;
static double* out;

void g(const long pos, const long rem_limit)
{
    if (pos == maxlen) {
			if (count_flag) {
				counter++;
			} else {
            	for (long i=0; i < maxlen;i++) {
                	*(out++) =  a[i];
				}
			}
    } else {
    	for (long limit=0; limit <= rem_limit; limit++) {
            	a[pos] = limit;
            	g(pos+1, rem_limit - limit);
    	}
	}
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
	/* check input args */
	if (nrhs < 4)
	{
		mexErrMsgTxt("polmatrix : Input arguments D(imension), maxgrad, Delay and time-series must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long D = (long) *((double *)mxGetPr(prhs[0]));
	const long maxgrad = (long) *((double *)mxGetPr(prhs[1]));
	const long Delay = (long) *((double *)mxGetPr(prhs[2]));
	
	const double* ts = ((double *)mxGetPr(prhs[3]));
	const long n = mxGetM(prhs[3]) * mxGetN(prhs[3]);	// length of time-series
	
	if ((D < 1) || (maxgrad <0) || (Delay < 1)) {
		mexErrMsgTxt("Wrong parameters given");
		return;
	}
	
	embedded_time_series_point_set emb_vecs(n, D, Delay, ts); 
	const long N = emb_vecs.size(); 	// number of time-delay vectors
	
	if (N < 1) {
		mexErrMsgTxt("Wrong parameters given (too short time-series");
		return;
	}
	
	maxlen = D;
    a = new long[D];

	count_flag = 1;
	counter = 0;		// number of monoms
	g(0, maxgrad);		// first do a run without output to see how many polynoms we will get

	count_flag = 0;
	plhs[0] = mxCreateDoubleMatrix(D, counter, mxREAL);
	out = (double *) mxGetPr(plhs[0]);		
	const double* monoms = out;
	
    g(0, maxgrad);

	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(N, counter, mxREAL);
		double* polmat = (double *) mxGetPr(plhs[1]);		
		
		for (long i = 0; i < N; i++) {
			for (long j = 0; j < counter; j++) {
				const double* monom_exp = monoms + j * D;
				embedded_time_series_point_set::point_iterator vec = emb_vecs.point_begin(i);
				
				double x = 1;
				
				for (long d=0; d<D; d++) {
					x *= pow(*vec, *monom_exp);
					//cout << *vec << " " << *monom_exp << endl;
				
					++vec;
					++monom_exp;	
				}
			
				polmat[i + N*j] = x;
			}
		}	
	}

    delete[] a;
}



