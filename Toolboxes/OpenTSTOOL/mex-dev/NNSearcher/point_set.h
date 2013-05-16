#ifndef POINT_SET_H
#define POINT_SET_H

#include <cmath>

// This file gives an example class for the implementation of a point_set which can
// be used by the nearest neighbor algorithm
// This particular implementation can be used for use with Matlab mex-files, where
// a point set is given as a Fortran style matrix (column major). The coordinates of
// one point is given by one row of this matrix

// interleaved_pointer	: a smart pointer that iterates with an interleave of #inc elements over an array
// May be used to iterate over rows or columns of a dense matrix of type T
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
		inline T operator[](const long i) const {
			return ptr[i*increment];
		}		
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

template<class METRIC>
class point_set_base {
	protected:
		const long N;				// number of points			
		point_set_base<METRIC>(const long n) : N(n) {};
	
	public:
		~point_set_base<METRIC>() {};
		inline long size() const { return N; };
		typedef METRIC Metric;
}; 

// The next class models a set of points in some space. Points can be accessed by an index, ranging from
// zero to N-1. The implementation of the points is not important as long as the class provides the
// possibility to calculate distances between two points of the point set and the distance between a point
// from the data set and an externally given point.
// Class point_set is parametrized by the METRIC that is used top compute distances.
// METRIC must be a class having an operator(). For possible implementations of a METRIC,
// see file "metric.h" in this directory
// This point_set implementation is the standard tstool way of handling matlab matrix data
template <class METRIC>
class point_set : public point_set_base<METRIC> {			// this implementation is for matlab matrix data
     // needed since gcc 3.4
     using point_set_base<METRIC>::N;
	protected:
		const long D;		// dimension
			
		const double* const matrix_ptr; 	// points are stored as row vectors of a N by D fortran style matrix
		const METRIC Distance;		// a function object that calculates distances		
	public:
		point_set(const long n, const long d, const double* mat) : point_set_base<METRIC>(n), D(d), matrix_ptr(mat), Distance() {};
		point_set(const long n, const long d, const double* mat, const METRIC& metr) : point_set_base<METRIC>(n), D(d), matrix_ptr(mat), 
			Distance(metr) {};	
			
		~point_set() {};
		inline long dimension() const { return D; }; 
				
		typedef interleaved_pointer<double> row_iterator; // a smart pointer that iterates over the elements of one point in the point_set (points are row vectors)
		
		row_iterator point_begin(const long n) const { return row_iterator(matrix_ptr + n, N); }
		row_iterator point_end(const long n) const { return row_iterator(matrix_ptr + n + N*D, N); }	// past-the-end	

		// indices may vary between 0 and N-1, or 0 and D-1
		// vec2 must be a double vector of length D, vec1[0] ... vec1[D-1]	
		double coordinate(const long n, const long d) const { return matrix_ptr[n + N*d]; };

		template<class ForwardIterator>
		inline double distance(const long index1, ForwardIterator vec2) const
		{
			return Distance(point_begin(index1), point_end(index1), vec2); 
		}

#ifdef PARTIAL_SEARCH
		template<class ForwardIterator>
		inline double distance(const long index1, ForwardIterator vec2, const double thresh) const
		{
			return Distance(point_begin(index1), point_end(index1), vec2, thresh); 
		}
#endif

		inline double distance(const long index1, const long index2) const
		{		
			return Distance(point_begin(index1), point_end(index1), point_begin(index2)); 
		}
			
		template<class ForwardIterator>
		inline void add(ForwardIterator vec, const long index) const
		{
			for (register long d=0; d < D; d++) vec[d] += matrix_ptr[index + N*d];
		}
};

// The next class is a models a set of points which is given by the (time-delay) embedding of a
// scalar time series. The time delay vectors are created "on the fly", they are not stored
// in memory. This slow down things a little bit, but keeps memory consumption really low.
// Even if we do not explicitly store the delay vectors in a matrix, we assume that delay vectors are
// row vectors.
// Attention : The elemtents of one vector are stored in opposite order compared to the order of the
// input time series. This means that the first component of a vector is the newest (which means having highest
// index) value of the section of the time-series.
template <class METRIC>
class embedded_time_series_point_set : public point_set_base<METRIC> {		
	protected:
		const long D;		// dimension
		
		const long DELAY;
		
		double* const ts_ptr; 	
		
		const METRIC Distance;		// a function object that calculates distances
					
	public:
		typedef interleaved_pointer<double> point_iterator; // a smart pointer that iterates over the elements of one point in the embedded_time_series_point_set (points are row vectors)

		// n - Length of scalar time-series
		// dim - Embedding dimension
		// delay - Time delay in samples
		// ts - zero based vector containing time series data
		embedded_time_series_point_set(const long n, const long dim, const long delay, double* ts) 
			: point_set_base<METRIC>(n-(dim-1)*delay), D(dim), DELAY(delay), ts_ptr(ts), Distance() {};
		embedded_time_series_point_set(const long n, const long dim, const long delay, double* ts, const METRIC& metr) 
			: point_set_base<METRIC>(n-(dim-1)*delay), D(dim), DELAY(delay), ts_ptr(ts), Distance(metr) {};	
			
		~embedded_time_series_point_set() {};
		
		inline long dimension() const { return D; }; 
				
		point_iterator point_begin(const long i) const { return point_iterator(ts_ptr + i + (D-1)*DELAY, -DELAY); }
		point_iterator point_end(const long i) const { return point_iterator(ts_ptr + i - DELAY, -DELAY); }	// past-the-end	

		// indices may vary between 0 and N-1, or 0 and D-1
		// vec2 must be a double vector of length D, vec1[0] ... vec1[D-1]	
		//double coordinate(const long n, const long d) const { return ts_ptr[n + d * DELAY]; };
	
		template<class ForwardIterator>
		inline double distance(const long index1, ForwardIterator vec2) const
		{
			return Distance(point_begin(index1), point_end(index1), vec2); 
		}

#ifdef PARTIAL_SEARCH
		template<class ForwardIterator>
		inline double distance(const long index1, ForwardIterator vec2, const double thresh) const
		{
			return Distance(point_begin(index1), point_end(index1), vec2, thresh); 
		}
#endif
		
		inline double distance(const long index1, const long index2) const
		{		
			return Distance(point_begin(index1), point_end(index1), point_begin(index2)); 
		}			
};


// calculate distances on the surface of a sphere of radius 1
// points must be given as tupels (phi, lambda) where phi is the latitude
// and lamda is the longitude
class spherical_distance {
	public: 
		spherical_distance() {};
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2) const 
		{
			const double phi1 = *first1;
			const double lamda1 = *(++first1);	
			const double phi2 = *first2;
			const double lamda2 = *(++first2);
			const double cose = sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2)*cos(lamda1-lamda2);			
  			return acos(cose);
		}
};

class spherical_point_set : public point_set<spherical_distance>
{
	public:
		spherical_point_set(const long n, const double* mat) :
			point_set<spherical_distance> (n, 2, mat) {};
};

// define a point_set that is transposed to the standard tstool notation, that means here points are row vectors
template <class METRIC>
class C_point_set  : public point_set_base<METRIC> {	// this implementation is for C sytle matrix data
	protected:
		const long D;		// dimension
		const double* const matrix_ptr; 	// points are stored as row vectors of a N by D C style matrix
		const METRIC Distance;		// a function object that calculates distances		
	public:
		C_point_set(const long n, const long d, const double* mat) : point_set_base<METRIC>(n), D(d), matrix_ptr(mat), Distance() {};
		C_point_set(const long n, const long d, const double* mat, const METRIC& metr) : point_set_base<METRIC>(n), D(d), matrix_ptr(mat), 
			Distance(metr) {};		
		~C_point_set() {};
		
		inline long dimension() const { return D; }; 

		typedef const double* row_iterator; // pointer that iterates over the elements of one point in the C_point_set (points are row vectors)
		
		row_iterator point_begin(const long n) const { return matrix_ptr + n*D; }
		row_iterator point_end(const long n) const { return matrix_ptr + (n+1)*D; }	// past-the-end	
	
		template<class ForwardIterator>
		inline double distance(const long index1, ForwardIterator vec2) const
		{
			return Distance(point_begin(index1), point_end(index1), vec2); 
		}

#ifdef PARTIAL_SEARCH
		template<class ForwardIterator>
		inline double distance(const long index1, ForwardIterator vec2, const double thresh) const
		{
			return Distance(point_begin(index1), point_end(index1), vec2, thresh); 
		}
#endif

		inline double distance(const long index1, const long index2) const
		{		
			return Distance(point_begin(index1), point_end(index1), point_begin(index2)); 
		}		
		
		double coordinate(const long n, const long d) const { return matrix_ptr[n*D + d]; };
};


#endif
