#include "mex.h"

extern void M_wrapper_double( int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray*prhs[], 
    double (*f) (const double*, const int), int normalize );

extern void M_wrapper_int( int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray*prhs[], 
    int (*f) (const double*, const int), int normalize );