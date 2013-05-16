/* Matlab MEX implementation of Kalafut-Visscher step-finding algorithm.

   Usage:
   Input:
    x          - Input signal: must be a row vector

   Output:
    intervals  - List of step points: zeros should be ignored

   (c) Max Little, 2010. Adapted from code written by B. Kalafut et al.
   The algorithm is described in B. Kalafut, K. Visscher,
   "An objective, model-independent method for detection of non-uniform
   steps in noisy signals", Comp. Phys. Comm., 179(2008), 716-723.
   If you use this code for your research, please cite:
   "Steps and bumps: precision extraction of discrete states of molecular
   machines using physically-based, high-throughput time series analysis"
   Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Input parameters */
#define  DATA_IN     prhs[0]

/* Output parameters */
#define  STEPS_OUT   plhs[0]

/*Improved step detector, take 3*/

typedef struct dl_list
{
    struct dl_list *next;
    struct dl_list *prior;

    double *data;
    int size;
    double var;
}node;

double variance(double *, int);
double mean(double *, int);
double SIC(double, int, int);
node * insert(node *,node *);

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
    double *data;            /* input vector */
    double *output;          /* output steps vector */

    double SIC_curr,SIC_opt,grand_var,var1,var2,var_prior,var_post;
    int steploc,c;
    int numsteps=0;
    int flag;
    node *currnode,*newnode,*optnode,*currnode2;
    node *base=(node *)mxCalloc(1,sizeof(node));

    int size = mxGetM(DATA_IN);  /* get it's size */
    data = mxGetPr(DATA_IN);    /* get pointer access to input vector */
    STEPS_OUT = (double *)mxCreateDoubleMatrix(size, 1, mxREAL);
    output = mxGetPr(STEPS_OUT);

    base->data=data; //initializing base node
    base->prior=NULL;
    base->next=NULL;
    base->size=size;
    base->var=variance(data,size);

    SIC_opt=SIC(base->var,base->size,2);//initialize optimum SIC

    for(;;) //empty loop control.  We'll terminate with a "break" below
    {
        currnode=base;
        flag=0;
        while(currnode!=NULL) //break loop if there's no such thing as a next node
        {
            var_prior=0.;//sum up n*variance for previous nodes
            currnode2=currnode->prior;
            while(currnode2!=NULL)
            {
                var_prior+=currnode2->size*currnode2->var;
                currnode2=currnode2->prior;
            }

            var_post=0.; //now sum it up for following nodes
            currnode2=currnode->next;
            while(currnode2!=NULL)
            {
                var_post+=currnode2->size*currnode2->var;
                currnode2=currnode2->next;
            }
            for(c=0;c<(currnode->size-1);c++)//test all possible change points in the current node
            {
                var1=variance(currnode->data,c+1);
                var2=variance(currnode->data+c+1,currnode->size-(c+1));

                grand_var=(c+1)*var1+(currnode->size-(c+1))*var2;
                grand_var+=var_prior+var_post;
                grand_var=grand_var/(size);
                SIC_curr=SIC(grand_var,size,3+numsteps);
                if(SIC_curr<SIC_opt)
                {
                    SIC_opt=SIC_curr;
                    steploc=c;
                    flag=1;
                    optnode=currnode;
                }
            }
            currnode=currnode->next; //move on to next node
        }
        if (flag==1)
        {
                newnode=insert(optnode,optnode->next); //we break "optnode" in two at the change point.  new node is right half
                newnode->data=optnode->data+steploc+1; //give the new node its data pointer
                newnode->size=optnode->size-(steploc+1); //set new sizes
                optnode->size=steploc+1;
                optnode->var=variance(optnode->data,optnode->size); //set new variances
                newnode->var=variance(newnode->data,newnode->size);
                numsteps++;//we've added a step
        } else
            break;
    }

    steploc=0;
    for(c=0,currnode=base;currnode!=NULL;c++,currnode=currnode->next) //we read out the data
    {
        steploc+=currnode->size;
        output[c]=steploc;
    }

    /*Now we free dynamically allocated memory for good measure.  No need to do this if calling from Labview.*/
    for(currnode=base;;)
    {
        if(currnode->next!=NULL)
            currnode=currnode->next;
        else
            break;
    }
    while(currnode!=base)
    {
        currnode=currnode->prior;
        mxFree(currnode->next);
    }
    mxFree(base);

    return;
}


/*Allocates a new linked-list node and inserts it between prev and post*/
node * insert(node *prev, node *post)
{
    node *newnode=(node *)mxCalloc(1,sizeof(node));
    if (prev!=NULL) prev->next=newnode;
    if (post!=NULL) post->prior=newnode;
    newnode->prior=prev;
    newnode->next=post;
    return newnode;
}

/*compute (MLE of) mean of 1-D array of double*/
double mean(double *data, int size)
{
    double sum=0.;
    int c;
    for(c=0;c<size;c++)
        sum+=data[c];
    return sum/size;
}

/*compute MLE of variance of 1-D array of double*/
double variance(double *data, int size)
{
    double avg=mean(data,size);
    int c;
    double sum=0;
    for(c=0;c<size;c++)
        sum+=(data[c]-avg)*(data[c]-avg);
    return sum/size;
}

/*
 * Computes SIC for model with mean changepoints
 * one parameter for each mean and one for the variance
 * constant terms are dropped.
 *
 */
double SIC(double variance,int size,int params)
{
    return (double)params*log(size)+(double)size*log(variance);
}
