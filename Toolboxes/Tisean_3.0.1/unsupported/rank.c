/* box assisted sorting/ranking utilities    */
/* Copyright (C) T. Schreiber (1998)         */

#include <stdlib.h>
#include <stdio.h>
#define MAX(a,b) ((a)>(b)? (a) : (b))
#define MIN(a,b) ((a)<(b)? (a) : (b))
#define EMPTY (-1)

void minmax(long nmax, double *x, double *xmin, double *xmax)
{
   long i;

   for(*xmin=x[0], *xmax=x[0], i=0; i<nmax; i++){
       *xmax=MAX(*xmax,x[i]);
       *xmin=MIN(*xmin,x[i]);
   }
}

void rank(long nmax, double *x, long *list)
{
   long *jptr, nl, n, i, ip, ipp;
   double xmin, xmax, sc, xn;

   nl=nmax/2;
   jptr=(long *)malloc(nl*sizeof(long));
   for(n=0; n<nl; n++) 
      jptr[n]=EMPTY;
   minmax(nmax, x, &xmin, &xmax);
   if(xmax==xmin){
      for(n=0; n<nmax; n++)
         list[n]=n;
      return;
   }
   sc=(nl-1)/(xmax-xmin);

   for(n=0; n<nmax; n++){
      xn=x[n];
      i=(long)((xn-xmin)*sc);
      ip=jptr[i];
      if( (ip==EMPTY) || (xn<=x[ip]) ) {
         jptr[i]=n;
      } else {
	do{
           ipp=ip;
	   ip=list[ip];
        }
	while( (ip!=EMPTY) && (xn>x[ip]) );
	list[ipp]=n;
      }
      list[n]=ip;
   }

   for(n=0, i=0; i<nl; i++){
      ip=jptr[i];
      while(ip!=EMPTY){
         ipp=ip;
         ip=list[ip];
         list[ipp]=++n;
      }
   }

   for(n=0; n<nmax; n++) list[n]-=1;
   free(jptr);
}

void rank2index(long nmax, long *list)
{
   long n, ib, im, it;

   for(n=0; n<nmax; n++)
      list[n]=-list[n];

   for(n=0; n<nmax; n++){
      if(list[n]<0){
         ib=n;
         im=-list[n];
         while(1){
            it=-list[im];
            list[im]=ib;
	    if(it==n){
               list[n]=im;
               break;
            } 
            ib=im;
            im=it;
         }
      }
   }
}

void indexx(long nmax, double *x, long *list)
{
   rank(nmax,x,list);
   rank2index(nmax,list);
}

void rank2sort(long nmax, double *x, long *list)
{
   long n, ib, it;
   double hb, ht;

   for(n=0; n<nmax; n++)
      list[n]=-list[n];

   for(n=0; n<nmax; n++){
      if(list[n]<0){
         ib=n;
         hb=x[n];
         while(1){
            it=-list[ib];
            list[ib]=it;
	    ht=x[it];
            x[it]=hb;
	    if(it==n) break;
            ib=it;
            hb=ht;
          }
       }
   }
}

void sort(long nmax, double *x)
{
   long *list;

   list=(long *)malloc(nmax*sizeof(long));
   rank(nmax,x,list);
   rank2sort(nmax,x,list);
   free(list);
}


