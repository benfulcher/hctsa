#include "mex.h"

#include "CO_AutoCorr.h"
#include "DN_HistogramMode_10.h"
#include "DN_HistogramMode_5.h"
#include "DN_OutlierInclude.h"
#include "FC_LocalSimple.h"
#include "IN_AutoMutualInfoStats.h"
#include "MD_hrv.h"
#include "PD_PeriodicityWang.h"
#include "SB_BinaryStats.h"
#include "SB_CoarseGrain.h"
#include "SB_MotifThree.h"
#include "SB_TransitionMatrix.h"
#include "SC_FluctAnal.h"
#include "SP_Summaries.h"
#include "butterworth.h"
#include "fft.h"
#include "helper_functions.h"
#include "histcounts.h"
#include "splinefit.h"
#include "stats.h"

void M_wrapper_double( int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray*prhs[], 
    double (*f) (const double*, const int), int normalize )
     
{ 
    double *inMatrix;       /* 1xN input matrix */
    int ncols;
    double *outMatrix; /* output matrix */
    
    // check inputs
    if(nrhs != 1) {
        mexErrMsgIdAndTxt("catch22:nrhs",
                          "One input required.");
    }
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("catch22:nlhs",
                          "One output required.");
    }
    if( !mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("catch22:notDouble",
            "Input vector must be type double.");
    }
    if(mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("catch22:notRowVector",
                          "Input must be a row vector.");
    }
    
    // get input
    inMatrix = mxGetPr(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    // set output
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    outMatrix = mxGetPr(plhs[0]);
    
    // calculate result
    if (normalize){
        
        double * y_zscored = malloc(ncols * sizeof * y_zscored);
        zscore_norm2(inMatrix, ncols, y_zscored);

        outMatrix[0] = f(y_zscored, ncols);

        free(y_zscored);
    } 
    else {
        outMatrix[0] = f(inMatrix, ncols);
    }   
    
    return;
    
}

void M_wrapper_int( int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray*prhs[], 
    int (*f) (const double*, const int), int normalize )
     
{ 
    double *inMatrix;       /* 1xN input matrix */
    int ncols;
    int *outMatrix; /* output matrix */
    
    // check inputs
    if(nrhs != 1) {
        mexErrMsgIdAndTxt("catch22:nrhs",
                          "One input required.");
    }
    if(nlhs > 1) {
        mexErrMsgIdAndTxt("catch22:nlhs",
                          "One output required.");
    }
    if( !mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("catch22:notDouble",
            "Input vector must be type double.");
    }
    if(mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("catch22:notRowVector",
                          "Input must be a row vector.");
    }
    
    // get input
    inMatrix = mxGetPr(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    // set output
    plhs[0] = mxCreateNumericMatrix(1,1, mxINT32_CLASS, mxREAL);
    outMatrix = (int*)mxGetData(plhs[0]);
    
    // calculate result
    if (normalize){
        
        double * y_zscored = malloc(ncols * sizeof * y_zscored);
        zscore_norm2(inMatrix, ncols, y_zscored);

        outMatrix[0] = f(y_zscored, ncols);

        free(y_zscored);
    } 
    else {
        outMatrix[0] = f(inMatrix, ncols);
    }   
    
    return;
    
}