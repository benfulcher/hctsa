/* build surrogates */
/* Copyright (C) Andreas Schmitz (1997) */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include "rand.h"
#include "dft.h"
#include "rank.h"

/* permute data (do p exchanges) */ 
void perm(double *data, long n, long p)
{
   long i,change1,change2;
   double help;

   for(i=0;i<p;i++) {
      change1=(long)((rnd_1279()/(double)ULONG_MAX)*n);
      change2=(long)((rnd_1279()/(double)ULONG_MAX)*n);
      help=data[change1];
      data[change1]=data[change2];
      data[change2]=help;
   }
}


/* with "Wiener filter" and re-indexing */
long pol_surr(long nmax, double *y, double *z, double *arr_sort, 
   double *ampl_org, long iters)
{
   double *phase, *y0;
   long n_ampl,it,i,*index;

   n_ampl=nmax/2+1;
   
   phase=(double *)malloc(n_ampl*sizeof(double));
   y0=(double *)malloc(nmax*sizeof(double));
   index=(long *)malloc(nmax*sizeof(long));
   if (phase==NULL || y0==NULL || index==NULL) {
      fprintf(stderr,"no memory\n");
      exit(5);
   }

   /* begin with random data */
   memcpy((void *)z,(void *)y,nmax*sizeof(double));
   perm(y,nmax,10*nmax);

   /* Polish; z has exact spectrum, y has exact amplitudes */
   for(it=0;it<iters;it++) {
      memcpy((void *)y0,(void *)y,nmax*sizeof(double));
      dft (y,nmax,NULL    ,phase);      /* "Wiener filter" */
      idft(y,nmax,ampl_org,phase);      /* spectrum right */
      memcpy((void *)z,(void *)y,nmax*sizeof(double));
      indexx(nmax,y,index); /* rescale */
      for (i=0;i<nmax;i++) y[index[i]]=arr_sort[i];
      if(memcmp((void *)y0,(void *)y,nmax*sizeof(double))==0) break;
   }

   free(index);
   free(phase);
   return it;
}
