/*=================================================================
 *
 * NEAREST.C	.MEX file corresponding to NEAREST.M
 *	        returns the index vector of nearest neighbours to 
 *              each embedded point in x
 *
 * The calling syntax is:
 *
 *		[ind] = nearest(X,tau)
 *
 *
 * This is a MEX-file for MATLAB.  
 *=================================================================*/
/* $Revision: 1.5 $ */
#include <math.h>
#include "mex.h"

void nearest(double	*x,
	     double	*ind,
	     double     *avect, /*normalisation factor*/
	     int tau,
	     int m,
	     int n)
/* for each column of x find the index of the RMS closest other column, excluding 
   those closer than tau to that column*/
{
  int i,j,k,oj,oi=0;
  int closest;
  double dist,diff,bestdist;
  for (i=0; i<n; i++) /*for each column of the matrix x*/{
    closest=0;
    bestdist=1e32;
    oj = 0;
    for (j=0; j<n; j++) /*compare to every other column*/{

      if (abs(i-j)>tau) /*the exclusion zone*/ {
	/*dist=rms(x[:,i]-x[:,j])*/
  	dist=0;

	for (k=0; k<m; k++) /* calculate rms */ {
	  diff= (*(x+oi+k)) - (*(x+oj+k));
	  dist += diff*diff*(*(avect+k));
	} /*of for k */

	if (dist<=bestdist) /*is this the best?*/ {
	  bestdist=dist;
	  closest=j;
	} /*of if (dist<bestdist)*/
      } /* of if (abs(i-j)>tau) */
      
      oj += m; /* =j*m */

    } /* of for j */

    *(ind+i)=closest+1;
    oi += m; /* =i*m */

  } /* of for i */
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *x,*ind,*avect;
    double dtau;
    int i,tau,mrows,ncols,mvect,nvect; 
    bool needfree=0;
    
    /* Check for proper number of arguments */
    
    if (nrhs>3) { 
	mexErrMsgTxt("Too many input arguments."); 
    } else if (nrhs == 0) {
        mexErrMsgTxt("Insufficient input arguments.");
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /*Assign a pointer to the input matrix*/
    x = mxGetPr(prhs[0]);

    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */     
    mrows = mxGetM(prhs[0]); 
    ncols = mxGetN(prhs[0]);
    
    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(1, ncols, mxREAL); 

    /* Get the size of the "forbidden zone" --- the second input arg. */
    if (nrhs>=2){
      dtau=mxGetScalar(prhs[1]);
      if (floor(dtau)!=dtau) {
	mexWarnMsgTxt("Second input should be an integer.");
      }
      tau=(int)floor(dtau);
    } else {
      tau=0;
    }

    /*norm weighting vector --- the third input arg */
    if (nrhs==3) {
      avect=mxGetPr(prhs[2]);
      mvect = mxGetM(prhs[2]); 
      nvect = mxGetN(prhs[2]);
      if (mvect*nvect<mrows) {
	mexErrMsgTxt("Third argument too short.");
      }
    } else {
      avect= (double *)malloc(sizeof(double)*mrows);
      for (i=0; i<mrows; i++) {
	*(avect+i) =1;
      }
      needfree=1;
    }


    /* Assign pointer to the ouput matrix */ 
    ind = mxGetPr(plhs[0]);
    
        
    /* Do the actual computations in a subroutine */
    nearest(x,ind,avect,tau,mrows,ncols); 
    
    /*free memory*/
    if (needfree) {
      free((double*)avect);
    }

    return;
    
}


