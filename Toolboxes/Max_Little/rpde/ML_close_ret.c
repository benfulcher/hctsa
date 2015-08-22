/* Close returns code by M. Little (c) 2006 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Real variable type */
#define  REAL        double

/* Input parameters */
#define  X_IN        prhs[0]
#define  EMBEDDIM_IN prhs[1]
#define  EMBEDDEL_IN prhs[2]
#define  ETA_IN      prhs[3]

/* Output parameters */
#define  CRS_OUT     plhs[0]

/* function definition */
#define  SYNTAX      "close_returns = close_ret(x, embed_dim, embed_delay, eta)"


/* Create embedded version of given sequence */
static void embedSeries
(
   unsigned long embedDims,      /* Number of dimensions to embed */
   unsigned long embedDelay,     /* The embedding delay */
   unsigned long embedElements,  /* Number of embedded points in embedded sequence */
   REAL          *x,             /* Input sequence */
   REAL          *y              /* (populated) Embedded output sequence */
)
{
   unsigned int i, d, inputDelay;

   for (d = 0; d < embedDims; d ++)
   {
      inputDelay = (embedDims - d - 1) * embedDelay;
      for (i = 0; i < embedElements; i ++)
      {
         y[i * embedDims + d] = x[i + inputDelay];
      }
   }
}


/* Search for first close returns in the embedded sequence */
static void findCloseReturns
(
   REAL           *x,               /* Embedded input sequence */
   REAL           eta,              /* Close return distance */
   unsigned long  embedElements,    /* Number of embedded points */
   unsigned long  embedDims,        /* Number of embedding dimensions */
   unsigned long   *closeRets        /* Close return time histogram */
)
{
   REAL  eta2 = eta * eta;
   REAL  diff, dist2;
   unsigned long  i, j, d, timeDiff, etaFlag;

   for (i = 0; i < embedElements; i ++)
   {
      closeRets[i] = 0;
   }

   for (i = 0; i < embedElements; i ++)
   {
      j = i + 1;
      etaFlag = 0;
      while ((j < embedElements) && !etaFlag)
      {
         dist2 = 0.0f;
         for (d = 0; d < embedDims; d ++)
         {
            diff   = x[i * embedDims + d] - x[j * embedDims + d];
            dist2 += diff * diff;
         }

         if (dist2 > eta2)
         {
            etaFlag = 1;
         }

         j ++;
      }

      etaFlag = 0;
      while ((j < embedElements) && !etaFlag)
      {
         dist2 = 0.0f;
         for (d = 0; d < embedDims; d ++)
         {
            diff   = x[i * embedDims + d] - x[j * embedDims + d];
            dist2 += diff * diff;
         }

         if (dist2 <= eta2)
         {
            timeDiff = j - i;
            closeRets[timeDiff] ++;
            etaFlag = 1;
         }

         j ++;
      }
   }
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
   long           rows, columns, vectorElements, embedElements, i;
   unsigned long  embedDims, embedDelay;
   REAL           etaIn;
   
   REAL           *CRSOut;          /* Output vector of close return counts */
   REAL           *sequenceIn;      /* Input vector */
   unsigned long  *closeRets;       /* Close return counts */
   REAL           *embedSequence;   /* Embedded input vector */

   /* Check for proper number of arguments */
   if ((nrhs != 4) || (nlhs != 1))
   {
      mexErrMsgTxt("Incorrect number of parameters.\nSyntax: "SYNTAX);
   }

   /* Checks on input sequence vector */
   rows           = mxGetM(X_IN);
   columns        = mxGetN(X_IN);
   vectorElements = columns * rows;
   if (columns != 1)
   {
      mexErrMsgTxt("Input sequence must be a row vector.");
   }
   if (!mxIsDouble(X_IN) || mxIsComplex(X_IN))
   {
      mexErrMsgTxt("Input sequence must be floating-point real.");
   }
   sequenceIn = mxGetPr(X_IN);

   /* Checks on close return distance */
   if (!mxIsDouble(ETA_IN) || mxIsComplex(ETA_IN))
   {
      mexErrMsgTxt("Close return distance eta must be floating-point real.");
   }
   etaIn = *mxGetPr(ETA_IN);

   /* Checks on embedding dimension */
   if (!mxIsNumeric(EMBEDDIM_IN) || mxIsComplex(EMBEDDIM_IN))
   {
      mexErrMsgTxt("Embedding dimension must be an integer.");
   }
   embedDims = (unsigned long)*mxGetPr(EMBEDDIM_IN);

   /* Checks on embedding delay */
   if (!mxIsNumeric(EMBEDDEL_IN) || mxIsComplex(EMBEDDEL_IN))
   {
      mexErrMsgTxt("Embedding delay must be an integer.");
   }
   embedDelay = (unsigned long)*mxGetPr(EMBEDDEL_IN);

   /* Create embedded version of input sequence */
   embedElements = vectorElements - ((embedDims - 1) * embedDelay);
   embedSequence = (REAL *)mxCalloc(embedElements * embedDims, sizeof(REAL));
   embedSeries(embedDims, embedDelay, embedElements, sequenceIn, embedSequence);

   /* Find close returns */
   closeRets = (unsigned long *)mxCalloc(embedElements, sizeof(unsigned long));
   findCloseReturns(embedSequence, etaIn, embedElements, embedDims, closeRets);

   /* Create output vectors, get pointer access */
   CRS_OUT = mxCreateDoubleMatrix(embedElements, 1, mxREAL);
   CRSOut  = mxGetPr(CRS_OUT);
   for (i = 0; i < embedElements; i ++)
   {
      CRSOut[i] = closeRets[i];
   }

   /* Release allocated memory. */
   mxFree(embedSequence);
   mxFree(closeRets);

   return;
}
