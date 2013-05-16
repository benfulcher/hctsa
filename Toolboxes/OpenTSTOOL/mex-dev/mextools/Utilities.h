#ifndef MY_UTILITIES
#define MY_UTILITIES

#include <math.h>

typedef unsigned long uint32;

#define MT_BUFSZ       (624)                 // length of state vector
// RANGE_LE1 puts a uint32 in the range [0,1] (0<=x<=1)
#define RANGE_LE1 2.3283064370807974e-10

// RANGE_LT1 puts a uint32 in the range [0,1) (0<=x<1)
#define RANGE_LT1 2.3283064370807974e-10

class RNG 
{
	private:
    	uint32   state[MT_BUFSZ+1];  // state vector + 1 extra to not violate ANSI C
    	uint32   *next;          	// next random value is computed from here
    	int      left;      		// can *next++ this many times before reloading

    	uint32 reloadMT();

	public:
    	void Seed(const uint32);
    	uint32 lRand();
    	inline double dRand() {
		
			return ( (double)lRand() * RANGE_LT1 ); 	// reals: [0,1)-interval
		}

    	RNG(const uint32 seed = 615460891);
};

#undef RANGE_LE1 
#undef RANGE_LT1 
#undef MT_BUFSZ

class My_Utilities {	// some utility functions
	private:
		RNG rng;	// random number generator
	
	public :
		
		// define a function that returns a random indices out of 0,1,2, ... ,N-1 
		uint32 randindex(const unsigned long N) { 
			double X = rng.dRand();
			uint32 R = (uint32) floor(X * ((double)N)); 
			return R;
		} 
			
		static double dmax(const double a, const double b) {return ((a >= b) ? a : b);}
		static double dmin(const double a, const double b) {return ((a <= b) ? a : b);}
		static long lmax(const long a, const long b) {return ((a >= b) ? a : b);}
		static long lmin(const long a, const long b) {return ((a <= b) ? a : b);}
		static double squared(const double a) { return a*a; }
		template<class T>
    	static void swap(T* a, const long i, const long j) { T tmp = a[i]; a[i] = a[j]; a[j] = tmp; }
};


#endif
