/* Functions for dft and idft (uses sing.c) */
/* Copyright (C) Andreas Schmitz (1997) */

#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <errno.h>

#include "sing.h"
#include "primes.h"

/* give me a data-array, I'll give you amplitudes and phases */
void dft(double *data, long count, double *ampl, double *phase)
{
   double *realpt,*imagpt;
   double norm;
   long halfcount,i,j;

   /* count is even */
   if (count%2 == 0) {
      halfcount=count/2;
      realpt=(double *)malloc(2*(halfcount+1)*sizeof(double));
      imagpt=(double *)(realpt+(halfcount+1));
      
      for( i=0,j=0; i<halfcount; i++ ) {
	 realpt[i] = data[j++];
	 imagpt[i] = data[j++];
      }
      realpt[halfcount] = imagpt[halfcount] = 0.0;

      sing(realpt, imagpt, halfcount, -1 );
      realtr( realpt, imagpt, halfcount, -1 );

      norm=2.0*count;
      for(i=1;i<halfcount;i++) {
	 if (ampl) ampl[i]=sqrt(realpt[i]*realpt[i]+imagpt[i]*imagpt[i])/norm;
	 if (phase) phase[i]=-atan2(imagpt[i],realpt[i]);
      }
      if (ampl) {
	 ampl[0]=realpt[0]/norm;
	 ampl[halfcount]=realpt[halfcount]/norm;
      }
      if (phase) {
	 phase[0]=0.0;
	 phase[halfcount]=0.0;
      }
      
      free(realpt);
   }
   /* count is odd */
   else {
      halfcount=count/2;
      realpt=(double *)malloc(2*(count+1)*sizeof(double));
      imagpt=(double *)(realpt+(count+1));

      for( i=0; i<count; i++ ) {
	 realpt[i] = data[i];
	 imagpt[i] = 0.0;
      }

      sing(realpt, imagpt, count, -1 );

      norm=count;
      for(i=1;i<=halfcount;i++) {
	 if (ampl) ampl[i]=sqrt(realpt[i]*realpt[i]+imagpt[i]*imagpt[i])/norm;
	 if (phase) phase[i]=-atan2(imagpt[i],realpt[i]);
      }
      if (ampl) ampl[0]=realpt[0]/norm;
      if (phase) phase[0]=0.0;

      free(realpt);
   }
}

/* give me amplitudes and phases, I'll give you real data */
void idft(double *data, long count, double *ampl, double *phase)
{
   double *realpt,*imagpt;
   long halfcount,i,j;

   /* count is even */
   if (count%2 == 0) {
      halfcount=count/2;
      realpt=(double *)malloc(2*(halfcount+1)*sizeof(double));
      imagpt=(double *)(realpt+(halfcount+1));

      realpt[0]=ampl[0];
      imagpt[0]=0.0;
      for(i=1;i<halfcount;i++) {
	 realpt[i]= ampl[i]*cos(phase[i]);
	 imagpt[i]=-ampl[i]*sin(phase[i]);
      }
      realpt[halfcount]=ampl[halfcount];
      imagpt[halfcount]=0.0;
      
      realtr( realpt, imagpt, halfcount, 1 );
      sing( realpt, imagpt, halfcount, 1 );

      for( i=0, j=0; i<halfcount; i++ ) {
	 data[j++] = realpt[i];
	 data[j++] = imagpt[i];
      }
      free(realpt);
   }
   /* count is odd */
   else {
      halfcount=count/2;
      realpt=(double *)malloc(2*(count+1)*sizeof(double));
      imagpt=(double *)(realpt+(count+1));

      realpt[0]=ampl[0];
      imagpt[0]=0.0;
      for(i=1;i<=halfcount;i++) {
	 realpt[count-i]=  realpt[i]= ampl[i]*cos(phase[i]);
	 imagpt[count-i]=-(imagpt[i]=-ampl[i]*sin(phase[i]));
      }

      sing( realpt, imagpt, count, 1 );

      for( i=0; i<count; i++ ) {
	 data[i] = realpt[i];
      }
      free(realpt);
   }
}

/* give back a lower value, so that product of primefactors <= 20000 */
/* highest allowed primefactor is 17 */
long strip_primes(long n)
{
   long k,i,product,prime,lastprime,table_nr=7;   /* 7th prime is 17 */
   
   if (n<=primetable[table_nr-1]) return n;

   for (;n>=primetable[table_nr-1];n--) {
      k=n;
      product=1;
      lastprime=0;
      for (i=0;i<table_nr;) {
	 prime=primetable[i];
	 if (k%prime == 0) {
	    k/=prime;
	    if (prime!=lastprime) {
	       product*=prime;
	       lastprime=prime;
	    }
	 }
	 else i++;
      }
      if (k==1 && product<20000) break;
   }

   return n;
}
