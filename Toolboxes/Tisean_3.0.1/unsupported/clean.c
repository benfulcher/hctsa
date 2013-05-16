#include <math.h>
#include <stdlib.h>
#include <stddef.h>

#include "clean.h"
#include "eigen.h"

extern struct list *top, *akt, *postakt;

double **matrix(long nr, long nc)
{
   long i;
   double **m;
   
   m=(double **)malloc(nr*sizeof(double*));
   for(i=0; i<nr; i++)
      m[i]=(double *)malloc(nc*sizeof(double));
   return m;
}

long neigh(double *y, long n, long d, long m, long kmax, long npmax, double eps,
   long *nlist)
{
   long nn, i, nfound;

   for(nfound=0, nn=n-1; nn >= (m-1)*d && nfound < kmax; nn--){
      for(i=0; i<m; i++)
         if(fabs((double)(y[n-i*d]-y[nn-i*d]))>=eps) break;
      if(i==m) nlist[nfound++]=nn;   /* fill list of neighbours */
      if(npmax && nn < n-npmax) break;
   }
   return nfound;
}

long train(long n, double *y, long d, long m, long nq, double **cm, 
          long kmin, long kmax, long npmax, double eps, double *r)
{
   long i,j,np,nfound,iq;  
   double s;
   static double **c=NULL;
   static long *nlist=NULL;
   struct list *sec;

   if(c==NULL) c=matrix(m,m);
   if(nlist==NULL) nlist=(long *)malloc(kmax*sizeof(long));
   nfound=neigh(y,n,d,m,kmax,npmax,eps,nlist);

   for(i=0; i<m; i++){
      for(s=y[n-(m-i-1)*d], np=0; np<nfound; np++)
         s+=y[nlist[np]-(m-i-1)*d];
      cm[n][i]=s/(nfound+1); 
   }

   if (nfound<kmin) return 0;

   sec=malloc(sizeof(*sec));
   sec->tcm=(double *)malloc(m*sizeof(double));
   sec->tse=matrix(m,m);
   sec->n=n;
   for(i=0; i<m; i++){
      for(s=0, np=0; np<nfound; np++)
         s+=cm[nlist[np]][i];
      sec->tcm[i]=2*cm[n][i]-s/nfound;
   }
   for(i=0; i<m; i++){
      for(j=i; j<m; j++){      /*  compute covariance matrix */
         for(s=0, np=0; np<nfound; np++){
            s+=(y[nlist[np]-(m-i-1)*d]-sec->tcm[i])
               *(y[nlist[np]-(m-j-1)*d]-sec->tcm[j]);
   	    }
         c[j][i]=c[i][j]=r[i]*r[j]*s/nfound;
      }
   }                  
   eigen(c,m,nq);       /* find eigenvectors (increasing)*/
   for(i=0; i<m; i++)
      for(iq=0; iq<nq; iq++)
         sec->tse[iq][i]=c[i][iq];
   sec->pre=top;
   top=sec;
   return 1;
}

void clean(long n, double *y, double *yc, long d, long m, long nq, 
   long kmin, long kmax, long npmax, double eps, double delta, double *r,
   long stdp, long *stdpm)
{
   long i, j, iq, flag=1, stcnt;
   double s;
   static double **cm=NULL;

   cm=(double **)realloc(cm, (n+1)*sizeof(double*));
   cm[n]=(double *)malloc(m*sizeof(double));

   for(akt=top, stcnt=0; akt!=NULL && stcnt<stdp;
                postakt=akt, akt=akt->pre, stcnt++){
      if(npmax && akt->n<n-npmax){ /* old, remove from stack */
         if (akt==top)
            top = top->pre;
	 else
            postakt->pre = akt->pre;
      } else {
         for(i=0;i<m;i++) 
            if(fabs((double)(y[n-i*d]-akt->tcm[i]))>=delta) break;
         if(i==m) break; /* representative found */
      }
   }
   *stdpm+=stcnt;
   if (akt==NULL || stcnt==stdp)
      flag=train(n,y,d,m,nq,cm,kmin,kmax,npmax,eps,r);
   else {
      if (akt!=top){
         postakt->pre = akt->pre;
         akt->pre = top; 
         top = akt;
      }
      for(j=0;j<m;j++) cm[n][j]=top->tcm[j];
   }

   if (flag) {
      for(i=0; i<m; i++){
         for(s=0, iq=0; iq<nq; iq++)
            for(j=0; j<m; j++) 
               s+=(y[n-(m-j-1)*d]- top->tcm[j])*
                  top->tse[iq][i]* top->tse[iq][j]*r[j];
         yc[n-(m-i-1)*d]-=(y[n-(m-i-1)*d]-top->tcm[i] -s/r[i])/r[i];
      }
   }
}



