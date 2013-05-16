#include <cmath>
#include "mextools/mextools.h"

// Generalized dimension estimation (Renyi dimensions)
//
// Method given by Piet Schram, Willem van der Water "Generalized Dimensions from
// Nearest Neighbor Information"
// 
// We compute moments of the distances to the k-th nearest neighbor for a sequence of k (4..128)
// and compare these to theoretical values
//
// Robust estimation is used instead of mean square error fitting
//
// cmerk 1999

using namespace std;

// mex -I. -I.. gendimest.cpp -O

// determine machine precision of not defined elsewhere
#ifndef DBL_EPSILON
double DBL_EPSILON()
{
	double halfu = 0.5;

 	while(((double)(1.0 + halfu)) > 1.0) {
		halfu *= 0.5;
 	}

 	return 2 * halfu;
}	 
#endif

// This may look like C code, but it is really -*- C++ -*-
/*
 ************************************************************************
 *
 *                        Numerical Math Package
 *
 *                  Brent's one-dimensional minimizer 
 *
 *           finds a local minimum of a single argument function
 *                        over a given range
 *
 * Input
 *      double fminbr(ax,bx,f,tol)
 *      const double ax                 a and b, a < b, specify the interval
 *      const double bx                 the minimum is to be sought in
 *      UnivariateFunctor f             The function under consideration
 *      const double tol                Acceptable tolerance for the minimum
 *                                      location. It is an optional parameter
 *                                      with default value DBL_EPSILON
 *
 * Output
 *      Fminbr returns an estimate to the location of the minimum
 *      with accuracy 3*SQRT_EPSILON*abs(x) + tol.
 *      The procedure can only determine a local minimum, which coincides with
 *      the global one if and only if the function under investigation is
 *      unimodular.
 *      If a function being examined possesses no local minimum within
 *      the given interval, Fminbr returns either the left or the right end
 *      point of the interval, wherever the function value is smaller.
 *
 * Algorithm
 *      G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *      computations. M., Mir, 1980, p.202 of the Russian edition
 *
 * The function makes use of a "gold section" procedure combined with
 * a parabolic interpolation.
 * At each step the code operates three abscissae - x,v, and w.
 *      x - the last and the best approximation to the minimum location,
 *              i.e. f(x) <= f(a) or/and f(x) <= f(b)
 *          (if the function f has a local minimum in (a,b), then both
 *           conditions are met after one or two steps).
 *      v,w are previous approximations to the location of the minimum.
 *      They may coincide with a, b, or x (although the algorithm tries
 *      to make all u, v, and w distinct). 
 * Points x, v, and w are used to construct an interpolating parabola,
 * whose minimum is regarded as a new approximation to the minimum
 * of the function, provided the parabola's minimum falls within [a,b]
 * and reduces the current interval [a,b] to a larger extent than the
 * gold section procedure does.
 * When f(x) has a positive second derivative at the point of minimum
 * (which does not coincide with a or b) the procedure converges
 * superlinearly at a rate of about 1.324
 *
 * $Id: gendimest.cpp,v 1.3 2008/09/18 12:51:46 engster Exp $
 *
 ************************************************************************
 */

template<class UnivariateFunctor>
double fminbr(                          // An estimate to the min location
        const double ax,                // Specify the interval the minimum
        const double bx,                // to be sought in
        UnivariateFunctor& f,            // Function under investigation
        const double tol)               // Acceptable tolerance
{
  static const double r = (3-sqrt(5.0))/2;      // The golden section ratio
  
#ifndef DBL_EPSILON  
	static const double sqrt_eps = sqrt(DBL_EPSILON());
#else
	static const double sqrt_eps = sqrt(DBL_EPSILON);	
#endif
 
  double a = ax, b = bx;                // Current interval
  double v = a + r*(b-a);               // First step - always gold section
  double fv = f(v);
  double x = v;                         // the last and the best approximation
  double fx = fv;
  double w = v;                         // a previous approx to the min
  double fw = fv;

  while(1)               				// Main iteration loop
  {
    const double range = b-a;           // Interval where the minimum
                                        // is searched in
    const double midpoint = (a+b)/2;
    const double tol_act =              // The effective tolerance
                sqrt_eps*fabs(x) + tol/3;

       

    if( 2*fabs(x-midpoint) + range <= 4*tol_act ) {
		return x;                         // Acceptable approximation is found
	}
                                        // Compute a new step with the gold
                                        // section
    double new_step = r * ( x < midpoint ? b-x : a-x );


                        // Decide on the interpolation  
    if(fabs(x-w) >= tol_act  )          // If x and w are distinct
    {                                   // interpolatiom may be tried
      register double p;                // Interpolation step is calcula-
      register double q;                // ted as p/q; division operation
                                        // is delayed until last moment
      register double t;

      t = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*t;
      q = 2*(q-t);

      if( q > 0 )                       // Formulas above computed new_step
        p = -p;                         // = p/q with a wrong sign (on purpose).
      else                              // Correct this, but in such a way so
        q = -q;                         // that q would be positive

      if(fabs(p) < fabs(new_step*q) &&   // If x+p/q falls in [a,b] and is not
         p > q*(a-x+2*tol_act) &&       // too close to a and b, and isn't
         p < q*(b-x-2*tol_act)  )       // too large, it is accepted
           new_step = p/q;
                                        // If p/q is too large then the
                                        // gold section procedure would
                                        // reduce [a,b] to larger extent
    }

    if( fabs(new_step) < tol_act )       // Adjust the step to be not less
      new_step =  new_step > 0 ?        // than tolerance
        tol_act : -tol_act;

                                // Obtain the next approximation to min
                                // and reduce the encompassing interval
    register const double t = x + new_step;  // Tentative point for the min
    register const double ft = f(t);
    if( ft <= fx )
    {                                   // t is a better approximation
      ( t < x ? b : a ) = x;            // Reduce the interval so that
                                        // t would fall within it
      v = w;  w = x;  x = t;            // Assign the best approx to x
      fv=fw;  fw=fx;  fx=ft;
    }
    else                                // x remains the better approx
    {
      ( t < x ? a : b ) = t;            // Reduce the interval encompassing x
      
      if( ft <= fw || w==x )
      {
        v = w;  w = t;
        fv=fw;  fw=ft;
      }
      else if( ft<=fv || v==x || v==w )
        v = t, fv = ft;
    }
  }             // ===== End of loop =====

}

static double max(double a,double b) {return ((a >= b) ? a : b);}
static double min(double a,double b) {return ((a <= b) ? a : b);}

// compute log(gamma(xx))
static double gammaln(const double xx)
{
    const double cof[6] = { 76.18009172947146, -86.50532032941677,
            				24.01409824083091,-1.231739572450155,
            				0.1208650973866179e-2,-0.5395239384953e-5};
   
    double y = xx;
	const double x = xx;
	
    double tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
	
    double ser=1.000000000190015;
    
	for (int j=0;j<=5;j++) 
		ser += cof[j]/++y;
    
	return -tmp+log(2.5066282746310005*ser/x);
}

// compute first derivative of log(gamma(xx))
static double dloggamma(const double xx)
{
 	const double h = 0.001;
	return (gammaln(xx+h)-gammaln(xx-h)) / (2.0 * h);
}

// we don't employ mean square error, instead a robust log(1 + 0.5 * d^2) is used
double log_error(const long first, const long last, mvector& ref, mvector& y , const double a = 1)
{
	if (last-first+1 > 0) {
		double error = 0;

		for (long i = first; i <= last; i++)  {
			const double e = (y(i) - ref(i))/a;
			error += log(1 + 0.5 * e * e);
		}	

		return error;
	} 
	
	return 0;	
}

class Error_Function {
	public:
		mvector& y; 		// vector of reference values
		const long kmin;
		const long kmax;
		const double g;

		mvector x;
		mvector& z;			// computed function values

		struct Scale_Function {
			Error_Function& P;
			Scale_Function(Error_Function& p) : P(p) {};			
			double operator()(const double a) {
				for (long k = P.kmin; k <= P.kmax; k++)  {
					P.x(k) = a * P.z(k);
				}
				return log_error(P.kmin, P.kmax, P.y, P.x); 
			}
		};
				
	public:
		Error_Function(const long k_min, const long k_max, mvector& ref, mvector& fit, const double G)
			: y(ref), kmin(k_min), kmax(k_max), g(G), z(fit), x(k_max) {};
		~Error_Function() {};
		double operator()(const double D);		
};

double Error_Function::operator()(const double D) { 
	// First compute the expected curve F(k, gamma,D) = (GAMMA(k+gamma/D)/GAMMA(k))^(1/gamma)
	// for a given value of the dimension D (gamma denotes the moment that was computed)
	// It is not necessary to compute the exact values of the above denoted function,
	// just the increase of F with k is important : F(k+1) = F(k) * (k+g/D)/k

	if (g==0) {
		for (long k = kmin; k <= kmax; k++)  {
			z(k) = exp(dloggamma(k)/D);
		}				
	} else {
		double start = 1;
		z(kmin) = 1;

		for (long k = kmin; k <= kmax-1; k++)  {
			z(k+1) = pow(start *= (k+g/D)/k, 1.0/g);
		}		
	}

	// Adjust the absolute scaling of curve F for best match against the measured curve
	Scale_Function scale(*this);
	
	const double a = fminbr(0, 1e6, scale, 0.00001);
	
	for (long k = kmin; k <= kmax; k++)  {
		z(k) *= a;
	}	
	
	// Return the residual error (that should only depend on D)
	return log_error(kmin, kmax, y, z);		
}

void compute_moments(const double* distances, const long NNR, const long R, const long NR_gammas,
					 const double* gammas, double* mom)
{
	const double* ptr = distances;
	
	for (long k=0; k < NNR; k++) {
		long i;
		for (i = 0; i < NR_gammas; i++) 
			mom[k + i * NNR] = 0;

		for (long j=0; j < R; j++) {
			for (long i = 0; i < NR_gammas; i++) {
				const double g = gammas[i];
				if (g == 0) 
					mom[k + i * NNR] += log(*ptr);
				else
					mom[k + i * NNR] += pow(*ptr, g);
			}
			ptr++;
		}
		for (i = 0; i < NR_gammas; i++) {
			const double g = gammas[i];
			if (g == 0) {
				mom[k + i * NNR] = exp(mom[k + i * NNR]/R);
			} else {
				mom[k + i * NNR] = pow(mom[k + i * NNR] / R, 1.0/g);		
			}
		}
	}
}

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{				
	/* check input args */
	if (nrhs < 5)
	{
		mexErrMsgTxt("gendimest : distances, gammas, kmin_low, kmin_high and kmax must be given");
		return;
	}
	
	/* handle matrix I/O */
	
	const long R = mxGetM(prhs[0]);		// number of reference points
	const long NNR = mxGetN(prhs[0]);		// number of distances (to its nearest neighbors) for each reference point
	
	const double* const distances = (double *)mxGetPr(prhs[0]);	// matrix of distances to nearest neighbors
	
	const double* gammas = (double *)mxGetPr(prhs[1]);	// gamma values
	const long NR_gammas = mxGetM(prhs[1]) * mxGetN(prhs[1]);		// number of gamma values
	
	const long kmin_lower = (long) *(double *)mxGetPr(prhs[2]);
	const long kmin_upper = (long) *(double *)mxGetPr(prhs[3]);
	const long kmax = (long) *(double *)mxGetPr(prhs[4]);
	
	if (kmin_lower < 1) {
		mexErrMsgTxt("kmin must be positive");
		return;
	}	
	if (kmax <= kmin_upper) {
		mexErrMsgTxt("kmax must be greater kmin_upper");
		return;
	}	
	if (kmin_lower > kmin_upper) {
		mexErrMsgTxt("kmin_upper must be greater or equal to kmin_lower");
		return;
	}	
	if (NNR < kmax) {
		mexErrMsgTxt("kmax may not exceed length of v");
		return;
	}		

	mxArray* moments = mxCreateDoubleMatrix(NNR, NR_gammas, mxREAL);
	double* mom = (double *) mxGetPr(moments);

	mxArray* fits = mxCreateDoubleMatrix(NNR, NR_gammas, mxREAL);

	// compute curve from given distance values

	compute_moments(distances, NNR, R, NR_gammas, gammas, mom);

	// printf("Computed moments of neighbor distances\n");

	plhs[0] = mxCreateDoubleMatrix(NR_gammas, kmin_upper-kmin_lower+1, mxREAL);
    
	for (long kmin = kmin_lower; kmin <= kmin_upper; kmin++) {
		double* out = (double *) mxGetPr(plhs[0]) + (kmin-kmin_lower) * NR_gammas;

		for (long i = 0; i < NR_gammas; i++) { 
			const double g = gammas[i];
			const double DMin = max(0.05, -g / kmin);
			const double DMax = 128;

			if (DMin >= DMax) {
				mexErrMsgTxt("kmin to small");
				return;
			}		
		
			mvector v(moments, i+1);	// vector to hold the moments of nearest neighbor distances
			mvector z(fits, i+1);		// vector to hold the fitted curve
		
			Error_Function err(kmin, kmax, v, z, g);

			// the dimension is estimated by minimizing the mean square error
			// of the difference between the theoretical curve and the curve
			// that is computed using the distance information from nearest 
			// neighbors
			const double d = fminbr(DMin, DMax, err, 0.0001);
			
			out[i] = d;	// return estimated dimension
		}
	}
	if (nlhs > 1) { 		// return computed moments of distances
		plhs[1] = moments;
    }
	
	if (nlhs > 2) { 		// also return theoretical expected curve that fits best
		plhs[2] = fits;
	}
}	

