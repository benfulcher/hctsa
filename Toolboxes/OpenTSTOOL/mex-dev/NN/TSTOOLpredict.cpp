// 
// predict : predict len next points of an time-series of vector points,
//           using a local constant model
// 
//   input arguments :
//    ts - scalar time series 
//	  dim - embedding dimension
//	  delay - time delay in samples
//    len    - length of prediction in samples
//    NNR    - number of nearest neighbours to use for each prediction
//    mode (optional) - method of computing predicted values
//    exp_weighted_factor (optional) - parameter for the euclidian weighted distance calculation, may be choosen out of ]0,1] (see McNames Leuven,Belgium 1998)
// 
//  cmerk 1998
//
// mex -I. -I.. predict.cpp -O

#include "mextools/mextools.h"

// this includes the code for the nearest neighbor searcher and the prediction routines
#include "NNSearcher/point_set.h"
#include "NNSearcher/nn_predictor.h"
#include "include.mex"

inline double uniform_weight(const double r) {
	if (r <= 1) 
		return 1;
	else 
		return 0;
}

inline double biquadratic_weight(const double r) {
	const double y = 1 - r * r;
	return (y * y);					// McNames biweight-function (1 - (r/rmax)^2)^2
}

// Ben Fulcher, 2015-03-04 -- added due to compiling errors from nn_predictor.h
inline double linear_weight(const double r) {
	return 1 - r;
}

inline double tricubic_weight(const double r) {
	const double y = 1 - r * r * r;
	return (y * y * y);
}


template<class Searcher, class Point_set>
void cross_validate(Searcher& searcher, const Point_set& px, const long trajectory_length, const int mode, const long NNR,
	const long start_index, const long length)
{
	typename Searcher::point_set::sample points_actual = searcher.get_point_set().actual_sample();
	typename Point_set::sample px_actual = px.actual_sample(); 
	
	for (long t = start_index; t < start_index + length; t++) {
		vector<neighbor> v;
		double pred = 0;	
		double weight = 0;
		
		searcher.search_k_neighbors(v, NNR + 1, px.point_begin(t), t-trajectory_length, t+trajectory_length);
		
		for (long k=0; k < NNR; k++) {
			double w;
		
			const double r = v[k].dist() / v[NNR].dist();
		
			switch(mode/2) {
				case 0  :
					w = uniform_weight(r);
                    break;
				case 1  : 
                    linear_weight lw;
					w = lw(r);
					break;
				case 2  :
					w = biquadratic_weight(r);
					break;
				case 3  :
					w = tricubic_weight(r);
					break;
				default :
					w = 1;
					break;						
			}

			weight += w;

			if (mode & 1) {	// intergrated prediction
				pred += w * (points_actual[v[k].index()+1] - points_actual[v[k].index()]);
			} else {
				pred += w * points_actual[v[k].index()+1];
			}
		}

		pred /= (weight);

		if (mode & 1)	// integrated prediction
			pred += px_actual[t]; 

		// Feed prediction back to px (which is a copy of the original time-series), so
		// that we can do an iterated prediction run
		px_actual[t+1] = pred;		
	}
}



void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{		
	/* check input args */
	
	if (nrhs < 3)
	{
		mexErrMsgTxt("Local constant prediction : Data set of points (row vectors), length, number of NNR must be given, stepsize and mode (0-3) is optional\n");
		return;
	}
	
	/* handle matrix and parameter I/O */
	
	long mode = 0;		
	long stepsize = 1;
	
	const long N = mxGetM(prhs[0]);	
	const long dim = mxGetN(prhs[0]);
	const double* const p = (double *)mxGetPr(prhs[0]);
	
	const long length 	= (long) *((double *)mxGetPr(prhs[1]));
	const long NNR  	= (long) *((double *)mxGetPr(prhs[2]));
		
	if (nrhs > 3) 	
		stepsize	= (long) *((double *)mxGetPr(prhs[3]));
	
	if (nrhs > 4) 
		mode  = (long) *((double *)mxGetPr(prhs[4])); 	// averaging mode
	
	if (length < 1) {
		mexErrMsgTxt("Length must be greater zero");
		return;
	}		 
	if (NNR < 1) {
		mexErrMsgTxt("Number of nearest neighbors must be greater zero");
		return;
	}		
	if (stepsize < 1) {
		mexErrMsgTxt("Stepsize must be greater zero");
		return;
	}		
	
	point_set<euclidian_distance> points(N, dim, p);
		
	mexPrintf("Length of input data set        : %d\n", N);
	mexPrintf("Dimension                       : %d\n", dim);
	mexPrintf("Length of prediction            : %d\n", length);
	mexPrintf("Number of nearest neighbors     : %d\n", NNR);
	mexPrintf("Stepsize                        : %d\n", stepsize);	
		
	plhs[0] = mxCreateDoubleMatrix(length, dim, mxREAL);
	double* const x = (double *) mxGetPr(plhs[0]);
	
	ATRIA<point_set<euclidian_distance> > searcher(points, 1);	
			
	if (searcher.geterr()) {	
		mexErrMsgTxt("Error preparing searcher : Inconsistent parameters given ?");
		return;
	}	 			
			
	double*	const coord = new double[dim];
				
	// start prediction with last value of the data set			
	for (long d=0; d < dim; d++) coord[d] = points.coordinate(N-1, d);  
				
	switch(mode) {
		case 1 :
		{
			mexPrintf("Local constant prediction using direct weighting\n");
		
			for (long t = 0; t < length; t++) {
					vector<neighbor> v;
					searcher.search_k_neighbors(v, NNR+1, coord, -1, -1);
					
					for (long d=0; d < dim; d++)
						coord[d] = 0;
					
					double weight = 0;
					
					for (long k=0; k < NNR; k++) {
						const double r = v[k].dist() / v[NNR].dist();
						const double w = biquadratic_weight(r);
						weight+=w;
						for (long d=0; d < dim; d++)
							coord[d] += w * points.coordinate(v[k].index()+stepsize, d);	
					}
					for (long d=0; d < dim; d++)
						x[t+d*length] = (coord[d] /= weight);	
			}
			break;
		}	
		case 2 :
		{
			mexPrintf("Local constant prediction using integrated mean\n");
					
			for (long t = 0; t < length; t++) {
					vector<neighbor> v;
					searcher.search_k_neighbors(v, NNR, coord, -1, -1);
					
					for (long k=0; k < NNR; k++) {
						for (long d=0; d < dim; d++)
							coord[d] += (points.coordinate(v[k].index()+stepsize, d) - points.coordinate(v[k].index(), d))/NNR;	
					}
					for (long d=0; d < dim; d++)
						x[t+d*length] = coord[d];
			}
			break;
		}
		case 3 :
		{
			mexPrintf("Local constant prediction using integrated weighting\n");
			
			for (long t = 0; t < length; t++) {
					vector<neighbor> v;
					searcher.search_k_neighbors(v, NNR+1, coord, -1, -1);
					
					double weight = 0;
					for (long k=0; k < NNR; k++) {
						const double r = v[k].dist() / v[NNR].dist();
						weight +=  biquadratic_weight(r);
					}
						
					for (long k=0; k < NNR; k++) {
						const double r = v[k].dist() / v[NNR].dist();
						const double w = biquadratic_weight(r);
						for (long d=0; d < dim; d++)
							coord[d] += w * (points.coordinate(v[k].index()+stepsize, d) - points.coordinate(v[k].index(), d)) / weight;	
					}
					for (long d=0; d < dim; d++)
						x[t+d*length] = coord[d];
			}
			break;
		}		
		default : 
		{
			mexPrintf("Local constant prediction using direct mean\n");
					
			for (long t = 0; t < length; t++) {
					vector<neighbor> v;
					searcher.search_k_neighbors(v, NNR, coord, -1, -1);
					
					for (long d=0; d < dim; d++)
						coord[d] = 0;					
					
					for (long k=0; k < NNR; k++) {
						for (long d=0; d < dim; d++)
							coord[d] += points.coordinate(v[k].index()+stepsize, d);	
					}
					for (long d=0; d < dim; d++)
						x[t+d*length] = (coord[d] /= NNR);
			}
			break;
		}	
	}
	
	delete[] coord;
}	


