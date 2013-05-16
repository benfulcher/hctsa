#include <cmath>
#include <stack>
#include <vector>
#include <algorithm>

#include "mextools/mextools.h"
#include "mextools/Utilities.h"
#include "mextools/Utilities.cpp"

using namespace std;

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{	
	/* check input args */
	if (nrhs < 3)
	{
		mexErrMsgTxt("randref : start,end,Nref must be given");
		return;
	}
	
	/* handle matrix I/O */

	const long start = (long) *((double *)mxGetPr(prhs[0]));
	const long end	= (long) *((double *)mxGetPr(prhs[1])); 
	const long Nref = (long) *((double *)mxGetPr(prhs[2]));
	const long N = end-start+1;
	
	if (N < 1) {
		mexErrMsgTxt("start or end are not valid");
		return;
	}		
	if (Nref < 0) {
		if (Nref == -1) {
			plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
			double* x = (double *) mxGetPr(plhs[0]);
		
			for (long i = start; i <= end; i++) {
				*(x++) = (double) i;
			}
		} else {
			mexErrMsgTxt("Number of reference indices must be positive");
		}
		return;
	}	
	if (N <= Nref) {
		plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
		double* x = (double *) mxGetPr(plhs[0]);
			
		for (long i=start; i <= end; i++)
			*(x++) = i;
		
		return;
	}
	
	plhs[0] = mxCreateDoubleMatrix(Nref, 1, mxREAL);
	double* x = (double *) mxGetPr(plhs[0]);
		
	vector<double> y;
	static My_Utilities utils;
	
	for (long i=start; i <= end; i++)
		y.push_back(i);
	
	//cout << y.size() << endl;
	
	for (long i=N-1; i >=N-Nref; i--) {
		const long j = utils.randindex(i);
		std::swap(y[i],y[j]);
	}
	
	copy(y.end()-Nref, y.end(), x);
	
	sort(x,x+Nref);
}



