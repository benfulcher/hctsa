#ifndef METRIC_H
#define METRIC_H

#include "nn_aux.h"

// This header file defines templated function objects (functors) for distance calculations. They work on any reasonable
// container class by using forward iterators
// For good performance, no bound checking is done, so it's the programmer's resposibility to have everything properly
// allocated.
// 
// To support PARTIAL DISTANCE CALCULATION during the search phase (not the preprocessing phase!), operator() is overloaded.
// operator() with three arguments is the standard distance calculation and returns the exact distance.
// operator() with four arguments is the tresholded distance calculation with terminates as soon as the partial distance
// exceeds the given threshold value (but otherwise returns the exact distance). Carefull implementation can speed up
// the search considerably. When threshold is exceeded, INFINITY is returned to prevent the search functions to consider
// this value as valid distance.

class euclidian_distance {
	public:
		euclidian_distance() {}; 
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2) const 
		{
			const double y = (*first1 - *first2);
			double dist = y*y;
  			for (++first1, ++first2 ; first1 != last1; ++first1, ++first2) {
				const double x = (*first1 - *first2);
				dist += x * x;
			}
  			return sqrt(dist);
		}
		
		// support partial search : if partial distance exceeds thresh, stop computing of distance 
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2, const double thresh) const 
		{
			const double t = thresh * thresh;
			const double y = (*first1 - *first2);
			double dist = y*y;
			
			if (dist > t) 
				return INFINITY;

  			for (++first1, ++first2 ; first1 != last1; ++first1, ++first2) {
				const double x = (*first1 - *first2);
				dist += x * x;
				if (dist > t) 
					return INFINITY;
			}
  			return sqrt(dist);
		}		
};	

// This distance may be used with brute force searcher (where no assumptions on the metric are done).
// Do not use this with ATRIA !!!
class squared_euclidian_distance {
	public:
		squared_euclidian_distance() {}; 
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2) const 
		{
			const double y = (*first1 - *first2);
			double dist = y*y;
  			for (++first1, ++first2 ; first1 != last1; ++first1, ++first2) {
				const double x = (*first1 - *first2);
				dist += x * x;
			}
  			return dist;
		}
		
		// support partial search : if partial distance exceeds thresh, stop computing of distance 
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2, const double thresh) const 
		{
			const double y = (*first1 - *first2);
			double dist = y*y;
			
			if (dist > thresh) 
				return INFINITY;

  			for (++first1, ++first2 ; first1 != last1; ++first1, ++first2) {
				const double x = (*first1 - *first2);
				dist += x * x;
				if (dist > thresh) 
					return INFINITY;
			}
  			return dist;
		}		
};	

class exp_weighted_euclidian_distance {
	private:
		const double lamda; 	// 0 < lamda <= 1
	public:	
		exp_weighted_euclidian_distance() : lamda(1) {};
		exp_weighted_euclidian_distance(const double l) : lamda(l) {};
		 
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2) const 
		{
			const double y = (*first1 - *first2);
			double dist = y*y;
  			double weight = lamda;
			for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
				const double x = (*first1 - *first2);
				dist += weight * x * x;
				weight *= lamda;
			}
  			return sqrt(dist);
		}	
		// support partial search
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2, const double thresh) const 		
		{
			const double t = thresh * thresh;
			const double y = (*first1 - *first2);
			double dist = y*y;
			
			if (dist > t) 
				return INFINITY;
		
  			double weight = lamda;
			for (++first1, ++first2; first1 != last1; ++first1, ++first2) {
				const double x = (*first1 - *first2);
				dist += weight * x * x;
				if (dist > t) {
					return INFINITY;
				}
				weight *= lamda;		
			}
  			return sqrt(dist);
		}				
};	

class maximum_distance {
	public: 
		maximum_distance() {};
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2) const 
		{
			double dist = fabs(*first1 - *first2);
  			for (++first1, ++first2 ; first1 != last1; ++first1, ++first2) {
				const double x = fabs(*first1 - *first2);
				if (x > dist) dist = x;
			}
  			return dist;
		}
		// support partial search
		template <class ForwardIterator1, class ForwardIterator2>
		double operator() (ForwardIterator1 first1, ForwardIterator1 last1, ForwardIterator2 first2, const double thresh) const 
		{
			double dist = fabs(*first1 - *first2);
			if (dist > thresh) {
				return INFINITY;
			}
  			for (++first1, ++first2 ; first1 != last1; ++first1, ++first2) {
				const double x = fabs(*first1 - *first2);
				if (x > dist) {
					if (x > thresh) {
						return INFINITY;
					}
					dist = x;
				}
			}
  			return dist;
		}					
};	

#endif

