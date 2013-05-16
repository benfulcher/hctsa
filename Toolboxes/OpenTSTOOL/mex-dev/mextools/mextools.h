#ifndef MLABVECTOR_H
#define MLABVECTOR_H

#include <sys/types.h>
//#include <math.h>
//#include <stdio.h>
//#include <iostream.h>
//#include <math>
//#include <stdio>
#include <iostream>

// Determine if this C++ code is compiled to a mex-file
// Otherwise, attach code that will produce a stand-alone executable (using mex2main.*)
#ifdef MATLAB_MEX_FILE
	#include "mex.h"
	#undef malloc
	#undef realloc
	#undef free
	#define malloc mxMalloc
	#define realloc mxRealloc
	#define free mxFree
	#define printf mexPrintf
#else
	#include "mex2main.h"
	#include "mex2main.cpp"
#endif

// declare main class names 
class mvector;
class mmatrix;

// Simple double vector class
// The base used for indexing (zero in C, one in Fortran) is adjustable
// No automatic resizing, no bound checking, just a lightweight vector class
template<class T>
class onebased_vector {
	protected:
		T* v;
		int owner;
		const long base_index;
			
	public:
		onebased_vector(T* const V, const long Base_index = 1) : v(V-Base_index), owner(0), base_index(Base_index) {};
		onebased_vector(const long L, const long Base_index = 1) : v(0), owner(1), base_index(Base_index)
		{
			if (L > 0) {
				v = new T[L];
				v -= base_index;
			}
		}	
		~onebased_vector() { 
			if (owner) {
				v += base_index; 
				delete[] v; 
			}
		}
		
		T operator()(const long i) const { return v[i]; }
		T& operator()(const long i) { return v[i]; }
};

// Simple vector class for use in mex-files, one based indexing, so the C++ code will look very
// similar to the m-file code
//
// Can be used to access input arguments : mvector vec(prhs[0]);
// or to create output arguments : mvector out(plhs[0] = mxCreateDoubleMatrix(24, 1, mxREAL));
//
// for (long i=1; i <= 24; i++) out(i) = 0.0;

class mvector : public onebased_vector<double> {
	protected:
		const long Length;
		
		mvector(const mvector& a); 	 // give dummy copy constructor to prevent copying 
		mvector(double* V, const long length) :
			onebased_vector<double>((double*) V), Length(length) {};
				
	public:
		mvector(const mxArray* const a) : 
			onebased_vector<double>((double*) mxGetPr(a)), Length(mxGetM(a) * mxGetN(a)) {};
			
		// create a vector from a column (1..N) of a matlab M by N matrix
		mvector(const mxArray* a, const long col) : 
			onebased_vector<double>(((double*) mxGetPr(a))+ (col-1)*mxGetM(a)), Length(mxGetM(a)) {};

		mvector(const long L) : onebased_vector<double>(L), Length(L) {};
		long length() const { return Length; }
		
		typedef double* iterator;
		iterator begin() { return &(v[1]); }
		iterator end() { return &(v[Length+1]); }	// past the end (STL style)
		
		friend class mmatrix;
};

// Simple N by M double matrix with one based indexing for rows and for columns
class mmatrix : protected onebased_vector<double> {
	protected:
		const long M;
		const long N;
		
		mmatrix(const mmatrix& a);		// give dummy copy constructor to prevent copying 
		
	public:
		mmatrix(const mxArray* const a) :
			onebased_vector<double>((double*) mxGetPr(a), mxGetM(a)+1), M(mxGetM(a)), N(mxGetN(a)) {};
		mmatrix(const long m, const long n) : 
			onebased_vector<double>(m*n, m+1), M(m), N(n) {};	
			
		long getM() const { return M; }
		long getN() const { return N; }
		
		double operator()(const long i, const long j) const { return v[i+j*M]; }
		double& operator()(const long i, const long j) { return v[i+j*M]; }
};

#endif

