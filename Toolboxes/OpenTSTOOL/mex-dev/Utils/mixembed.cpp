/***********************************************************************\
 **
 **    Title:         mixed state embedding
 **
 **    File:          mixembed.cpp
 **
 **    Contents:      main program
 **                   
 **    Author:        Kevin Bube - Drittes Physikal. Insitut Uni Goettingen
 **
 **    Version 0.1:   4.10.2001
 **            0.11:  5.10.2001  introduced type vector
 **            0.3:   12.11.2001 added more embedding options
 **            0.5:   24.5.2002  integrated into TSTOOL-CVS
 **
 **    Bugs:   Under-/Overflow checks are missing
 **
 **    This is the implementation of a matlab mex-file. Its purpose is
 **    to embed time series into a mixed state vector.
 **    The invokation in matlab is mixembedalg (X, Y). For details type
 **    "help mixembed" at matlab prompt.
 **
 **    Copyright 1997-2002 DPI Goettingen
 **    License http://www.physik3.gwdg.de/tstool/gpl.txt
 **
 \***********************************************************************/

#include "mex.h"
#include "matrix.h"

typedef struct{
    double* data;
    int rows, cols;
}matrix;


/*------------------------------------------------------------------------
 *
 * internal function 'embed'
 *
 * This function embeds time series into a mixed state vector.
 *
 * parameter:     ts - mxn matrix of n time series of length m
 *                dims - 1xn matrix of embedding dimensions. One for each
 *                       time series
 *                lags - 1xn matrix of time lags (in samples)
 *                shifts - 1xn-1 matrix which holds the amount of samples
 *                         which the time series is shifted relatively to the
 *                         first one
 *                startindex - pointer to an integer variable which holds the
 *                             index of the time series value, which is the
 *                             first value of the first embedding vector
 *                             (this is used mainly for debugging purposes)
 *
 * return value:  pointer to Matlab array of embedding vector
 *
 * preconditions: The dimensions of the time series must be the same. It
 *                Has to be checked if the values for the dimensions etc.
 *                are sane.
 *
 * side effects:  -
 *-----------------------------------------------------------------------*/
static mxArray *embed (matrix ts, matrix dims, matrix lags,
		       matrix shifts, int *startindex);

/*------------------------------------------------------------------------
 *
 * internal function 'max'
 *
 * This function returns the index of the maximum value of a vector. If
 * there are more elements with the maximum value, the smallest index is
 * returned.
 *
 * parameter:     vector - pointer to double vector which should be searched
 *                length - integer value of the vector length
 *
 * return value:  integer value of the index of the maximum element
 *
 * preconditions: -
 *
 * side effects:  -
 *-----------------------------------------------------------------------*/
static int max (const int *vector, int length);

/*------------------------------------------------------------------------
 *
 * internal function 'get_p_max_index'
 *
 * This function returns the first index of the first time series
 * embedding vector.
 *
 * parameter:     dims - 1xn dimension matrix
 *                lags - 1xn lags matrix
 *
 * return value:  integer value of the maximum index
 *
 * preconditions: The arrays must be of the same lengths. The values of the
 *                dimensions must be nonegaive. The ones of the time lags
 *                must be positive.
 *
 * side effects:  -
 *-----------------------------------------------------------------------*/
static int get_p_max_index (matrix dims, matrix lags);



/**************************************************************************
 *
 * main function
 *
 *************************************************************************/

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i;
    matrix ts, dims, lags, shifts;
    int startindex;
    int lags_need_free = false, shifts_need_free = false;

    /* shifts needs init
     */
    shifts.data = NULL;
    shifts.rows = 0;
    shifts.cols = 0;
    
    /*-----------------------------------------------------------------------
     * some simple sanity checks and data formating
     *---------------------------------------------------------------------*/

    switch (nrhs){
      case 4:
       if (!(mxIsDouble (prhs[3]) &&
	     !mxIsComplex (prhs[3]))){
	   mexErrMsgTxt ("corrupt input arg 4");
       }

       shifts.data = mxGetPr (prhs[3]);
       shifts.rows = mxGetM (prhs[3]);
       shifts.cols = mxGetN (prhs[3]);
       
       /* fall through here */
      case 3:
       if (!(mxIsDouble (prhs[2]) &&
	     !mxIsComplex (prhs[2]))){
	   mexErrMsgTxt ("corrupt input arg 3");
       }

       lags.data = mxGetPr (prhs[2]);
       lags.rows = mxGetM (prhs[2]);
       lags.cols = mxGetN (prhs[2]);

       /* fall through here */
      case 2:
       if (!(mxIsDouble (prhs[1]) &&
	     !mxIsComplex (prhs[1]) &&
	     mxIsDouble (prhs[0]) &&
	     !mxIsComplex (prhs[0]))){
	   mexErrMsgTxt ("corrupt input arg 1 or 2");
       }

       dims.data = mxGetPr (prhs[1]);
       dims.rows = mxGetM (prhs[1]);
       dims.cols = mxGetN (prhs[1]);
       ts.data = mxGetPr (prhs[0]);
       ts.rows = mxGetM (prhs[0]);
       ts.cols = mxGetN (prhs[0]);

       break;
      default:
       mexErrMsgTxt ("wrong number of input arguments");
    }


    /*
     * set default values if required
     */
    switch (nrhs){
      case 4:
       /* nothing to be done */
       break;
       
      case 3:
       if (ts.cols > 1){
	   /* we don't need this stuff in case we have only one time series
	    */
	   shifts.data = (double*) mxMalloc ((ts.cols-1) * sizeof (double));
	   if (shifts.data == NULL)
	       mexErrMsgTxt ("mexFunction(): error while allocating memory");
	   
	   for (i = 0; i < ts.cols-1; i++)
	       shifts.data[i] = 0.0;

	   shifts.rows = 1;
	   shifts.cols = ts.cols-1;
	   shifts_need_free = true;
       }
       break;

      case 2:
       if (ts.cols > 1){
	   /* this would lead to an error if we had only one time series
	    */
	   shifts.data = (double*) mxMalloc ((ts.cols-1) * sizeof (double));
	   if (shifts.data == NULL)
	       mexErrMsgTxt ("mexFunction(): error while allocating memory");
	   
	   for (i = 0; i < ts.cols-1; i++)
	       shifts.data[i] = 0.0;
	   
	   shifts.rows = 1;
	   shifts.cols = ts.cols-1;
	   shifts_need_free = true;
       }
	   
       lags.data = (double*) mxMalloc (ts.cols * sizeof (double));
       if (lags.data == NULL)
	   mexErrMsgTxt ("mexFunction(): error while allocating memory");

       for (i = 0; i < ts.cols; i++)
	   lags.data[i] = 1.0;

       lags.rows = 1;
       lags.cols = ts.cols;
       lags_need_free = true;
       break;

      default:
       /* this should not happen */
       mexErrMsgTxt ("internal error! please call the men in black");
    }


    /*
     * check dimensions of the matrices we were given
     */
    if (lags.rows != 1 ||
	lags.cols != ts.cols ||
	dims.rows != 1 ||
	dims.cols != ts.cols){
	mexErrMsgTxt ("wrong matrix dimensions");
    }

    if (ts.cols > 1){
	/* skip this if we have only one time series
	 */
	if (shifts.rows != 1 ||
	    shifts.cols != ts.cols-1){
	    mexErrMsgTxt ("wrong matrix dimensions");
	}
    }

    
    /*----------------------------------------------------------------------
     * do the deed and format output
     *--------------------------------------------------------------------*/
    plhs[0] = embed (ts, dims, lags, shifts, &startindex);
    plhs[1] = mxCreateDoubleMatrix (1, 1, mxREAL);
    (mxGetPr (plhs[1]))[0] = (double) startindex;

    /*
     * clean up
     */
    if (lags_need_free)
	mxFree (lags.data);
    if (shifts_need_free)
	mxFree (shifts.data);
}


static mxArray *embed (matrix ts, matrix dims, matrix lags,
		       matrix shifts, int *startindex)
{
    int i, j, k=0, l=0;
    matrix p;                                 /* embedding vector */
    mxArray *mx_p;
    int xindex;
    
    
    /* get first index of time that should be predicted */
    xindex = get_p_max_index (dims, lags);

    /* determine size of Matlab array and create it */
    /* We simply take the biggest shift when calculating p.rows.
       This is not the optimal way as some time steps may be dropped
       but in long time series this should not matter. */
    p.rows = ts.rows - xindex - max ((int *)shifts.data, shifts.cols);
    for (i = 0, p.cols = 0; i < ts.cols; i++)
	p.cols += ((int) dims.data[i]);

    mx_p = mxCreateDoubleMatrix (p.rows, p.cols, mxREAL);
    if (mx_p == NULL)
	mexErrMsgTxt ("embed(): error while allocating memory");
    p.data = mxGetPr (mx_p);

    /* finally it is time to get our hands dirty and fill the bucket */
    for (i = 0; i < ts.cols; i++){
	for (j = 0; j < dims.data[i]; j++){

	    /* there is no shift for the first time series */
	    if (i == 0){
              for (k = xindex - ((int)lags.data[i]*j) - 1;
		     k < ts.rows - lags.data[i]*j - 1;
		     k++){
		    p.data[l] = ts.data[k + i*ts.rows];
		    l++;
		}
	    }
	    else{
              for (k = xindex - ((int)lags.data[i]*j) - 1 + (int)shifts.data[i-1];
		     k < ts.rows - lags.data[i]*j - 1 + shifts.data[i-1];
		     k++){
		    p.data[l] = ts.data[k + i*ts.rows];
		    l++;
		}
	    }
	}
    }

    *startindex = xindex;
    return mx_p;
}


static int get_p_max_index (matrix dims, matrix lags)
{
    int *prod;
    int i;
    
    /* determine partial max. indizes */
    prod = (int*) mxMalloc (dims.cols * sizeof (int));
    if (prod == NULL)
	mexErrMsgTxt ("get_p_max_index(): error while allocating memory");

    for (i = 0; i < dims.cols; i++)
      prod[i] = ((int)dims.data[i] - 1) * (int)lags.data[i] + 1;

    i = max (prod, dims.cols);
    i = prod[i];
    
    mxFree (prod);
    return i;
}


static int max (const int *vector, int length)
{
    int index_max=0, i;

    for (i = 0; i < length; i++){
	if (vector[i] > vector[index_max])
	    index_max = i;
    }
    
    return index_max;
}
