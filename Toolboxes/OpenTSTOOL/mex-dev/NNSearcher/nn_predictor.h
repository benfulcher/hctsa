#ifndef LOCAL_NN_PREDICTION_H
#define LOCAL_NN_PREDICTION_H

// Definitions for nearest neighbor searcher
#include "nn_aux.h"

struct constant_weight {
	double operator()(double x) { return 1.0; }; 	// return (x <= 1) ? 1.0  : 0.0; }
};

struct inverse_weight { 
	double operator()(double x) { 
		return 1.0/(1 + x);
	}  
};

struct gaussian_weight { 
	double operator()(double x) { 
		return exp(-x*x);
	}
};

struct linear_weight { 
	double operator()(double x) { 
		return 1 - x;
	}
};

// McNames biweight-function (1 - (r/rmax)^2)^2	  
struct biquadratic_weight { 
	double operator()(double x) { 
		const double y = 1 - x * x;
		return (y * y);				
	}
};

struct tricubic_weight { 
	double operator()(double x) { 
		const double y = 1 - x * x * x;
		return (y * y * y);				
	}
};

template<class ImageAccessor, class WeightingFunction>
struct local_approx {	
	ImageAccessor& getValue;
	WeightingFunction& weight;
	
	local_approx(ImageAccessor& gV, WeightingFunction& wg) : getValue(gV), weight(wg) {}
	~local_approx() {}

	double operator() (const vector<neighbor>& v) const {	
		const double r_max = v.back().dist();
		
		double totalWeight = 0;
		double y = 0;
		
		for (vector<neighbor>::const_iterator i=v.begin(); i != v.end(); ++i) {
			const double w  = weight(i->dist()/r_max);
			y += w * getValue(i->index());
			totalWeight += w;
		}
	
		if (totalWeight == 0) {
			cerr << "nn_predictor.h : local_approx : total weight is zero" << endl;
			return 0;
		} else {
			return y / totalWeight;
		}
	}
};

template<class NNSearcher>
class nn_predictor
{
	protected :
		int err;  		// error state, == 0 means OK, every other value is a failure
		NNSearcher& searcher;
		
	public:
		nn_predictor(NNSearcher& s); 
		~nn_predictor() {}		
		inline int geterr() const { return err; }
		
		template<class ForwardIterator, class LocalEstimator> 
		double predict(const long NNR, ForwardIterator preimage, LocalEstimator& loc_est, const long first = -1, const long last = -1);
		
};	

template<class NNSearcher>
nn_predictor<NNSearcher>::nn_predictor(NNSearcher& s) : err(s.geterr()), searcher(s)
{}

template<class NNSearcher>
template<class ForwardIterator, class LocalEstimator> 
double nn_predictor<NNSearcher>::predict(const long NNR, ForwardIterator preimage, LocalEstimator& loc_est, const long first, const long last)
{
	vector<neighbor> v;
	searcher.search_k_neighbors(v, NNR, preimage, first, last);
	
	if (v.size() != NNR)  {
		err = 1;
		return 0;
	} 
		
	return loc_est(v);
}


#endif

