/*==========================================================
 * sampen_mex.c -
 *
 * Converted to mex file for Matlab by Ben Fulcher
 *
 * This is a MEX-file for MATLAB.
 *
 *==========================================================
 file: sampen.c	Doug Lake	2 August 2002
 			Last revised:	1 November 2004 (by george@mit.edu) 1.2
 -------------------------------------------------------------------------------
 sampen: calculate Sample Entropy
 Copyright (C) 2002-2004 Doug Lake

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 Place - Suite 330, Boston, MA 02111-1307, USA.  You may also view the agreement
 at http://www.fsf.org/copyleft/gpl.html.

 You may contact the author via electronic mail (dlake@virginia.edu).  For
 updates to this software, please visit PhysioNet (http://www.physionet.org/).

 _______________________________________________________________________________

 Revision history:
   1.0 (2 August 2002, Doug Lake)	Original version
   1.1 (6 January 2004, George Moody)	Removed limits on input series length
   1.2 (1 November 2004, George Moody)	Merged bug fixes from DL (normalize
 					by standard deviation, detect and
 					avoid divide by zero); changed code to
 					use double precision, to avoid loss of
 					precision for small m and large N

 Compile this program using any standard C compiler, linking with the standard C
 math library.  For example, if your compiler is gcc, use:
     gcc -o sampen -O sampen.c -lm

 For brief instructions, use the '-h' option:
     sampen -h

 Additional information is available at:
     http://www.physionet.org/physiotools/sampen/.

 */

#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* sampen2 calculates an estimate of sample entropy and the variance of the
   estimate. */
void sampen2(double *y, int mm, double r, int n)
{
    double *p = NULL;
    double *v1 = NULL, *v2 = NULL, *s1 = NULL, dv;
    int *R1 = NULL, *R2 = NULL, *F2 = NULL, *F1 = NULL, *F = NULL, FF;
    int *run = NULL, *run1 = NULL;
    double *A = NULL, *B = NULL;
    double *K = NULL, *n1 = NULL, *n2 = NULL;
    int MM;
    int m, m1, i, j, nj, jj, d, d2, i1, i2, dd;
    int nm1, nm2, nm3, nm4;
    double y1;

    mm++;
    MM = 2 * mm;

    if ((run = (int *) mxCalloc(n, sizeof(int))) == NULL)
	exit(1);
    if ((run1 = (int *) mxCalloc(n, sizeof(int))) == NULL)
	exit(1);
    if ((R1 = (int *) mxCalloc(n * MM, sizeof(int))) == NULL)
	exit(1);
    if ((R2 = (int *) mxCalloc(n * MM, sizeof(int))) == NULL)
	exit(1);
    if ((F = (int *) mxCalloc(n * MM, sizeof(int))) == NULL)
	exit(1);
    if ((F1 = (int *) mxCalloc(n * mm, sizeof(int))) == NULL)
	exit(1);
    if ((F2 = (int *) mxCalloc(n * mm, sizeof(int))) == NULL)
	exit(1);
    if ((K = (double *) mxCalloc((mm + 1) * mm, sizeof(double))) == NULL)
	exit(1);
    if ((A = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((B = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((p = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((v1 = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((v2 = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((s1 = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((n1 = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);
    if ((n2 = (double *) mxCalloc(mm, sizeof(double))) == NULL)
	exit(1);

    for (i = 0; i < n - 1; i++) {
	nj = n - i - 1;
	y1 = y[i];
	for (jj = 0; jj < nj; jj++) {
	    j = jj + i + 1;
	    if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) {
		run[jj] = run1[jj] + 1;
		m1 = (mm < run[jj]) ? mm : run[jj];
		for (m = 0; m < m1; m++) {
		    A[m]++;
		    if (j < n - 1)
			B[m]++;
		    F1[i + m * n]++;
		    F[i + n * m]++;
		    F[j + n * m]++;
		}
	    }
	    else
		run[jj] = 0;
	}			/* for jj */

	for (j = 0; j < MM; j++) {
	    run1[j] = run[j];
	    R1[i + n * j] = run[j];

	}
	if (nj > MM - 1)
	    for (j = MM; j < nj; j++)
		run1[j] = run[j];
    }				/* for i */

    for (i = 1; i < MM; i++)
	for (j = 0; j < i - 1; j++)
	    R2[i + n * j] = R1[i - j - 1 + n * j];
    for (i = MM; i < n; i++)
	for (j = 0; j < MM; j++)
	    R2[i + n * j] = R1[i - j - 1 + n * j];
    for (i = 0; i < n; i++)
	for (m = 0; m < mm; m++) {
	    FF = F[i + n * m];
	    F2[i + n * m] = FF - F1[i + n * m];
	    K[(mm + 1) * m] += FF * (FF - 1);
	}

    for (m = mm - 1; m > 0; m--)
	B[m] = B[m - 1];
    B[0] = (double) n *(n - 1) / 2;
    for (m = 0; m < mm; m++) {
	p[m] = (double) A[m] / B[m];
	v2[m] = p[m] * (1 - p[m]) / B[m];
    }
    dd = 1;
    for (m = 0; m < mm; m++) {
	d2 = m + 1 < mm - 1 ? m + 1 : mm - 1;
	for (d = 0; d < d2 + 1; d++) {
	    for (i1 = d + 1; i1 < n; i1++) {
		i2 = i1 - d - 1;
		nm1 = F1[i1 + n * m];
		nm3 = F1[i2 + n * m];
		nm2 = F2[i1 + n * m];
		nm4 = F2[i2 + n * m];
		for (j = 0; j < (dd - 1); j++) {
		    if (R1[i1 + n * j] >= m + 1)
			nm1--;
		    if (R2[i1 + n * j] >= m + 1)
			nm4--;
		}
		for (j = 0; j < 2 * (d + 1); j++)
		    if (R2[i1 + n * j] >= m + 1)
			nm2--;
		for (j = 0; j < (2 * d + 1); j++)
		    if (R1[i2 + n * j] >= m + 1)
			nm3--;
		K[d + 1 + (mm + 1) * m] +=
		    (double) 2 *(nm1 + nm2) * (nm3 + nm4);
	    }
	}
    }

    n1[0] = (double) n *(n - 1) * (n - 2);
    for (m = 0; m < mm - 1; m++)
	for (j = 0; j < m + 2; j++)
	    n1[m + 1] += K[j + (mm + 1) * m];
    for (m = 0; m < mm; m++) {
	for (j = 0; j < m + 1; j++)
	    n2[m] += K[j + (mm + 1) * m];
    }

    for (m = 0; m < mm; m++) {
	v1[m] = v2[m];
	dv = (n2[m] - n1[m] * p[m] * p[m]) / (B[m] * B[m]);
	if (dv > 0)
	    v1[m] += dv;
	s1[m] = (double) sqrt((double) (v1[m]));
    }

    for (m = 0; m < mm; m++) {
	if (p[m] == 0)
	    printf("No matches! SampEn((%d,%g,%d) = Inf"
		   " (standard deviation = Inf)!\n", m, r, n);
	else
	    printf("SampEn(%d,%g,%d) = %lf (standard deviation = %lf)\n",
		   m, r, n, -log(p[m]), s1[m]);
    }

    mxFree(A);
    mxFree(B);
    mxFree(p);
    mxFree(run);
    mxFree(run1);
    mxFree(s1);
    mxFree(K);
    mxFree(n1);
    mxFree(R1);
    mxFree(R2);
    mxFree(v1);
    mxFree(v2);
    mxFree(F);
    mxFree(F1);
    mxFree(F2);
}

/* sampen() calculates an estimate of sample entropy but does NOT calculate
   the variance of the estimate */
void sampen(double *y, int M, double r, int n, double *sampEnt)
{
    double *p = NULL;
    double *e = NULL;
    long *run = NULL, *lastrun = NULL, N;
    double *A = NULL, *B = NULL;
    int M1, j, nj, jj, m;
    int i;
    double y1;

    /* Allocate memory: */
    M++;
    if ((run = (long *) mxCalloc(n, sizeof(long))) == NULL)
	exit(1);
    if ((lastrun = (long *) mxCalloc(n, sizeof(long))) == NULL)
	exit(1);
    if ((A = (double *) mxCalloc(M, sizeof(double))) == NULL)
	exit(1);
    if ((B = (double *) mxCalloc(M, sizeof(double))) == NULL)
	exit(1);
    if ((p = (double *) mxCalloc(M, sizeof(double))) == NULL)
	exit(1);

    /* start running */
    for (i = 0; i < n - 1; i++) {
	nj = n - i - 1;
	y1 = y[i];
	for (jj = 0; jj < nj; jj++) {
	    j = jj + i + 1;
	    if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) {
		run[jj] = lastrun[jj] + 1;
		M1 = M < run[jj] ? M : run[jj];
		for (m = 0; m < M1; m++) {
		    A[m]++;
		    if (j < n - 1)
			B[m]++;
		}
	    }
	    else
		run[jj] = 0;
	}			/* for jj */
	for (j = 0; j < nj; j++)
	    lastrun[j] = run[j];
    }				/* for i */

    N = (long) (n * (n - 1) / 2);
    p[0] = A[0] / N;
    sampEnt[0] = -log(p[0]);
    /* printf("SampEn(0,%g,%d) = %lf\n", r, n, sampEnt[0]); */

    for (m = 1; m < M-1; m++) {
	p[m] = A[m] / B[m - 1];
	if (p[m] == 0)
    {
        sampEnt[m] = 0;
        /* printf("No matches! SampEn((%d,%g,%d) = Inf!\n", m, r, n); */
    }
	else
    {
        sampEnt[m] = -log(p[m]);
	    /* printf("SampEn(%d,%g,%d) = %lf\n", m, r, n, sampEnt[m]);
        printf("A = %lf\n",A[m]));
        printf("B = %lf\n",B[m])); */
    }
    }

    mxFree(A);
    mxFree(B);
    mxFree(p);
    mxFree(run);
    mxFree(lastrun);
}

static char *help_strings[] = {
 "usage: %s [OPTIONS ...] [TEXT-FILE]\n",
 "where OPTIONS may include:",
 " -h    print this usage summary",
 " -m M  set the maximum epoch length to M (default: 2)",
 " -n    normalize such that the mean of the input is 0 and the sample",
 "       variance is 1",
 " -r R  set the tolerance to R (default: 0.2)",
 " -v    output an estimate of the standard deviation of each SampEn",
 "TEXT-FILE should contain the time series (a column of decimal numbers);",
 "if omitted, sampen reads its standard input.  The output contains values of",
 "SampEn(m,r,N) where m is the epoch length, r is the tolerance, and N is the",
 "length of the input series.",
 NULL
};

void help()
{
    int i;

    (void) fprintf(stderr, help_strings[0], "sampen");
    for (i = 1; help_strings[i] != NULL; i++)
	(void) fprintf(stderr, "%s\n", help_strings[i]);
}

/* The gateway function */
void mexFunction(int nlhs,              /* number of expected outputs */
                 mxArray *plhs[],       /* array of pointers to output arguments */
                 int nrhs,              /* number of inputs */
                 const mxArray *prhs[]) /* array of pointers to input arguments */
{
    double *x_in;               /* Input time-series vector */
    int numSamples, i;             /* number of samples in the time series */
    double m;                   /* embedding dimension */
    double r;                   /* threshold value */
    double *sampEnt;            /* The sample entropies computed from the algorithm */
    double *sampEnt_out;     /* outputs */

    /* Check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("hctsa:sampen:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("hctsa:sampen:nlhs","One output required.");
    }

    /* check that the first input argument is double */
    if( !mxIsDouble(prhs[0]) ||
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("hctsa:sampen:notScalar","Input must be a vector of doubles.");
    }

    /* check that the second input argument is double */
    if( !mxIsDouble(prhs[1]) ||
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("hctsa:sampen:notScalar","Second input must be a double.");
    }

    /* check that the third input argument is double */
    if( !mxIsDouble(prhs[2]) ||
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("hctsa:sampen:notScalar","Third input must be a double.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[0])!=1) {
        mexErrMsgIdAndTxt("hctsa:sampen:notRowVector","Input must be a row vector.");
    }

    /* create a pointer to the real data in the input vector, x_in  */
    x_in = mxGetPr(prhs[0]);

    /* Value of embedding dimension, m */
    m = mxGetScalar(prhs[1]);

    /* Pointer to threshold value, r */
    r = mxGetScalar(prhs[2]);

    /* Number of samples in the time series */
    numSamples = mxGetN(prhs[0]);

    /* Calculate number of scales needed */
    sampEnt = (double *)mxCalloc(m+1, sizeof(double));

    /* mexPrintf("\nRunning SampEn(%f,%f) on a time series with %d samples\n",m,r,numSamples); */

    /* call the computational routine */
    sampen(x_in, m, r, numSamples, sampEnt);

    /* assign output: */
    /* mexPrintf("\nSampEn(%f,%f) completed! Assigning output....\n",m,r); */

    plhs[0] = mxCreateDoubleMatrix(m+1, 1, mxREAL);
    sampEnt_out = mxGetPr(plhs[0]);

    for (i = 0; i < m; i ++)
    {
       sampEnt_out[i] = sampEnt[i];
       /*    mexPrintf("\n %u SampEn(%u) = %f.\n",i,i,sampEnt[i]); */
    }

    /* Release allocated memory. */
    mxFree(sampEnt);

    return;
}
