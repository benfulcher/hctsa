#ifndef AKIMA_H
#define AKIMA_H

#include <cmath>

// Class for akima interpolation
// The akima coefficients are computed when the constructor is called.
// After construction, the coefficients cannot be altered in any way

template<class InputIterator>
class akima_interpolator
{
	private:
		const long n;		// number of input knots (= pairs of x and y coordinates)
		
		// knot data
		double* x;			// knot data is copied to x and y while construction
		double* y;			// so the user can savely forget about his original data vectors
						
		// spline coefficient
		double* b;
		double* c;
		double* d;

		int iflag;	// != 0 indicates error state 
		
	public:
		// be careful, X and Y must contain N elements
		akima_interpolator(const long N, InputIterator X, InputIterator Y);
		~akima_interpolator();
		
		double eval(const double u) const;		// evalute spline at u
		
		int geterr() const { return iflag; }
};

template<class InputIterator>
akima_interpolator<InputIterator>::akima_interpolator(const long N, InputIterator X, InputIterator Y)
	: n(N), x(0), y(0), b(0), c(0), d(0), iflag(0)
{
	long i;
	
	if (n < 2)
	{  /* no possible interpolation */
  		iflag = 1;
  		return;
	}

	// allocate memory for spline coefficients
	b = new double[n];
	c = new double[n];	
	d = new double[n];

	// allocate and copy vectors with knot data
	x = new double[n+4];	
	y = new double[n+4];

	// allocate vectors for Steigungsdata
	
	double* m = new double[n+4];	

	x += 2;		// shift pointer so that we can access x[-2], x[-1], x[0], ..., x[n], x[n+1]		
	y += 2;		// same for y
	m += 2; 	

	int ascend = 1;
	
	// copy input point data
	for (i=0; i < n; i++) {
		x[i] = *X++;
		y[i] = *Y++;
	}
	
	for (i = 1; i < n; ++i) if (x[i] <= x[i-1]) ascend = 0;
	
	if (!ascend)
	{
	   iflag = 2;	// data is not given with ascending x values
	   return;
	}

	// extrapolate 2 points left and two points right
	x[-1] = x[0] + x[1] - x[2];
	x[-2] = x[-1] + x[0] - x[1];

	x[n] = x[n-1] + x[n-2] - x[n-3];
	x[n+1] = x[n] + x[n-1] - x[n-2];

	for (i=0; i < n-1; i++) {
		m[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]);
	}
   
    y[-1] =(x[0]-x[-1])*(m[1]-2*m[0])+y[0];
    m[-1] =(y[0]-y[-1])/(x[0]-x[-1]);
	
    y[-2] =(x[-1]-x[-2])*(m[0]-2*m[-1])+y[-1];
    m[-2] =(y[-1]-y[-2])/(x[-1]-x[-2]);

    y[n] = (2*m[n-2]-m[n-3])*(x[n]-x[n-1])+y[n-1];
	m[n-1] = (y[n]-y[n-1])/(x[n]-x[n-1]);
	
    y[n+1] = (2*m[n-1]-m[n-2])*(x[n+1]-x[n])+y[n];	 
    m[n] = (y[n+1]-y[n])/(x[n+1]-x[n]);
	
	for (i=0; i < n; i++) {
	     const double term1 = std::fabs(m[i+1] - m[i]);
	     const double term2 = std::fabs(m[i-1] - m[i-2]);
		const double t1_t2 = term1 + term2;
		
		if (t1_t2 > 0) {
			b[i] = (term1 * m[i-1] + term2 * m[i]) / t1_t2;
		} else {
			b[i] = 0;
		}
	}
	
	// now calcualte the coefficients for the n-1 cubic polynomials
	for (i=0; i < n-1; i++) {		
	  	c[i] =(3*m[i]-2*b[i]-b[i+1]) / (x[i+1]-x[i]);
      	d[i] =(b[i]+b[i+1]-2*m[i]) / ((x[i+1]-x[i]) * (x[i+1]-x[i]));
	}
	
	m -= 2;
	delete[] m;
} 

template<class InputIterator>
akima_interpolator<InputIterator>::~akima_interpolator(){
	x -= 2;
	y -= 2;

	delete[] x;
	delete[] y;
	
	delete[] b;
	delete[] c;
	delete[] d;
}

template<class InputIterator>
double akima_interpolator<InputIterator>::eval(const double u) const
{
	static long last;
	long  i, j, k;
	double w;
	
	i = last;
	
	if (i >= n-1) i = 0;
	if (i < 0)  i = 0;

	// perform a binary search
	if ((x[i] > u) || (x[i+1] < u))
	{  
		i = 0;
		j = n-1;
		do
		{
			k = (i + j) / 2;         // split the domain to search 
			if (u < x[k])  j = k;    // move the upper bound 
			if (u >= x[k]) i = k;    // move the lower bound 
		} while (j > i+1);           //  there are no more segments to search 
	}
	
	last = i;

	// printf("%ld\n", i);

	w = u - x[i];
	w = y[i] + w * (b[i] + w * (c[i] + w * d[i]));
	
	return w;
}


#endif



