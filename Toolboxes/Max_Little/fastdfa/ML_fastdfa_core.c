/* Performs fast detrended fluctuation analysis on a nonstationary input signal.

   Useage:
   Inputs
    x          - input signal: must be a row vector
   Optional inputs:
    intervals  - List of sample interval widths at each scale
                 (If not specified, then a binary subdivision is constructed)

   Outputs:
    intervals  - List of sample interval widths at each scale
    flucts     - List of fluctuation amplitudes at each scale

   (c) 2006 Max Little. If you use this code, please cite:
   M. Little, P. McSharry, I. Moroz, S. Roberts (2006),
   Nonlinear, biophysically-informed speech pathology detection
   in Proceedings of ICASSP 2006, IEEE Publishers: Toulouse, France.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Real variable type */
#define  REAL        double

/* Input parameters */
#define  X_IN        prhs[0]
#define  INTVS_IN    prhs[1]

/* Output parameters */
#define  INTVS_OUT   plhs[0]
#define  FN_OUT      plhs[1]

/* function definition */
#define  SYNTAX   "[intervals, flucts] = fastdfa_core(x, intervals)"

/* Calculate accumulated sum signal */
static void cumulativeSum
(
   unsigned long  elements,
   REAL           *x,
   REAL           *y
)
{
   unsigned int i;
   REAL accum = 0.0f;
   
   for (i = 0; i < elements; i ++)
   {
      accum += x[i];
      y[i]   = accum;
   }
}

/* Calculate intervals if not specified */
static void calculateIntervals
(
   unsigned long  elements,
   unsigned long  scales,
   unsigned long  *intervals
)
{
   unsigned long  subdivs;
   REAL           idx_inc;
   long           scale;

   for (scale = scales - 1; scale >= 0; scale --)
   {
      subdivs   = 1 << scale;
      idx_inc   = (REAL)elements/(REAL)subdivs;
      intervals[scale] = (unsigned long)(idx_inc + 0.5f);
   }
}


/* Measure the fluctuations at each scale */
static void dfa
(
   REAL           *x,
   unsigned long  elements,
   unsigned long  *intervals,
   REAL           *flucts,
   unsigned long  scales
)
{
   unsigned long  idx, i, start, end, iwidth, accum_idx;
   long  scale;

   REAL  Sy, Sxy;                   /* y and x-y components of normal equations */
   REAL  Sx, Sxx;                   /* x-component of normal equations */
   REAL  a, b;                      /* Straight-line fit parameters */
   REAL  *trend;                    /* Trend vector */
   REAL  diff, accum, delta;        /* Intermediate variables */

   trend   = (REAL *)mxCalloc(elements, sizeof(REAL));

   /* Calculate successive accumulated data points for Sy and Sxy at each scale
      Since the sample indices are linear, there are simple closed forms for Sxx and Sx. */
   for (scale = scales - 1; scale >= 0; scale --)
   {
      /* Sx/Sxy accumulation over each interval */
      for (accum_idx = 0, idx = 0; idx < elements; idx += intervals[scale], accum_idx ++)
      {
         /* Find left and right-hand ends of interval */
         start  = idx;
         end    = idx + intervals[scale] - 1;

         /* We'll have to miss out an interval smaller than can fit at the end of the sequence */
         if (end >= elements)
         {
            /* Any left-over elements not accounted for are treated the same as the input vector */
            for (i = start; i < elements; i ++)
            {
               trend[i] = x[i];
            }
            break;
         }
         iwidth = end - start + 1;

         /* Accumulate Sy and Sxy of straight-line fit */
         Sy  = 0.0f;
         Sxy = 0.0f;
         for (i = start; i <= end; i ++)
         {
            Sy  += x[i];
            Sxy += x[i] * (REAL)i;
         }

         /* Calculate closed form values for Sx and Sxx */
         Sx     = ((REAL)end + (REAL)start) * (REAL)iwidth / 2.0;
         Sxx    = (REAL)iwidth * (2 * (REAL)end * (REAL)end + 2 * (REAL)start * (REAL)start +
                                  2 * (REAL)start * (REAL)end + (REAL)end - (REAL)start) / 6.0;
         delta  = (REAL)iwidth * Sxx - (Sx * Sx);

         /* Solve normal equations for straight-line fit */
         b      = (Sy * Sxx - Sx * Sxy) / delta;
         a      = ((REAL)iwidth * Sxy - Sx * Sy) / delta;

         /* Calculate straight-line fit (the trend) */
         for (i = start; i <= end; i ++)
         {
            trend[i] = a * (REAL)i + b;
         }
      }

      /* Calculate fluctuation at this scale */
      accum = 0.0f;
      for (i = 0; i < elements; i ++)
      {
         diff   = x[i] - trend[i];
         accum += diff * diff;
      }
      flucts[scale] = sqrt(accum / (REAL)elements);
   }

   /* Clean up */
   mxFree(trend);
}



/* Main entry point */
/* lhs - output parameters */
/* rhs - input parameters */
void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
#if !defined(V4_COMPAT)
    const mxArray *prhs[]         /* array of pointers to input arguments */
#else
    mxArray *prhs[]         /* array of pointers to input arguments */
#endif
)
{
   unsigned long  samples, dimensions, elements, entries;
   unsigned long  N_scales, *intervals, i;
   REAL           *flucts;
   
   REAL     *intvs_out;       /* vector of intervals at each scale */
   REAL     *fluct_out;       /* vector of fluctuations at each scale */
   REAL     *x_in;            /* input vector */
   REAL     *y_in;            /* cumulative sum of input vector */

   /* Check for proper number of arguments */
   if (((nrhs != 1) && (nrhs != 2)) || (nlhs != 2))
   {
      mexErrMsgTxt("Incorrect number of parameters.\nSyntax: "SYNTAX);
   }

   /* Get size of input sequence */
   samples    = mxGetM(X_IN);
   dimensions = mxGetN(X_IN);
   elements   = samples * dimensions;

   if (dimensions != 1)
   {
      mexErrMsgTxt("Input sequence must be a vector.");
   }

   if (!mxIsDouble(X_IN) || mxIsComplex(X_IN))
   {
      mexErrMsgTxt("Input sequence must be floating-point real.");
   }

   /* Get pointer access to real part of input sequence */
   x_in = mxGetPr(X_IN);

   /* Calculate accumulated sum signal */
   y_in = (REAL *)mxCalloc(elements, sizeof(REAL));
   cumulativeSum(elements, x_in, y_in);

   /* Get size of interval input vector if supplied */
   if (nrhs == 2)
   {
      entries    = mxGetM(INTVS_IN);
      dimensions = mxGetN(INTVS_IN);
      N_scales   = entries * dimensions;

      if (dimensions != 1)
      {
         mexErrMsgTxt("Intervals must be a vector.");
      }

      if (mxIsComplex(X_IN))
      {
         mexErrMsgTxt("Intervals must be real.");
      }

      if (N_scales < 2)
      {
         mexErrMsgTxt("Number of intervals must be greater than one.");
      }

      /* Get pointer access to real part of interval input vector */
      x_in = mxGetPr(INTVS_IN);
      intervals = (unsigned long *)mxCalloc(N_scales, sizeof(unsigned long));
      for (i = 0; i < N_scales; i ++)
      {
         intervals[i] = (unsigned long)x_in[i];
      }

      /* Check that the interval widths are sensible */
      for (i = 0; i < N_scales; i ++)
      {
         if ((intervals[i] > elements) || (intervals[i] < 3))
         {
            mexErrMsgTxt("Invalid interval size: must be between size of sequence x and 3.");
         }
      }
   }
   else
   {
      N_scales  = (unsigned long)(log10(elements)/log10(2.0));
      if (((REAL)(1 << (N_scales - 1))) > ((REAL)elements/2.5f))
      {
         N_scales --;
      }
      intervals = (unsigned long *)mxCalloc(N_scales, sizeof(unsigned long));
      calculateIntervals(elements, N_scales, intervals);
   }

   /* Calculate number of scales needed */
   flucts = (REAL *)mxCalloc(N_scales, sizeof(REAL));

   /* Measure the fluctuations at each scale */
   dfa(y_in, elements, intervals, flucts, N_scales);

   /* Create output vectors, get pointer access */
   INTVS_OUT = mxCreateDoubleMatrix(N_scales, 1, mxREAL);
   FN_OUT    = mxCreateDoubleMatrix(N_scales, 1, mxREAL);
   intvs_out = mxGetPr(INTVS_OUT);
   fluct_out = mxGetPr(FN_OUT);
   for (i = 0; i < N_scales; i ++)
   {
      intvs_out[i] = intervals[i];
      fluct_out[i] = flucts[i];
   }

   /* Release allocated memory. */
   mxFree(y_in);
   mxFree(intervals);
   mxFree(flucts);

   return;
}
