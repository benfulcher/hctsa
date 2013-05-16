/* Calculate third order autocorrelation */

/****************************************************************/
/* this source is freely distributable                          */
/* copyright (C) Andreas Schmitz (1997)                         */
/* Wuppertal University, Germany                                */
/****************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "arguments.h"

#define MAXNAMELEN 512      /* maximal lenght of file name */
#define DEF_POINTS 10000000 /* default number of points */
#define DEF_DELAY1 1         /* default delay */

void usage(char *progname)
{
   fprintf(stderr,"\nUsage: %s [-d# -D# -l# -x# -c# -h] file(s)\n",progname);
   fprintf(stderr,"\t-d <delay one (1)>\n");
   fprintf(stderr,"\t-D <delay two (2x delay one)>\n");
   fprintf(stderr,"\t-l <number of points (whole file)>\n");
   fprintf(stderr,"\t-x <number of values to be skipped (0)>\n");
   fprintf(stderr,"\t-c <column to be read (1 or file,#)>\n");
   fprintf(stderr,"\t-h <show this message>\n\n");
   exit(0);
}

double third_auto(double *y,unsigned long n,long d1,long d2)
{
   double autocorr=0.0,mean=0.0,skew=0.0,ym;
   unsigned long i,maxd;

   maxd=d2;
   if (d1>d2) maxd=d1;
   
   for (i=maxd;i<n;i++) mean+=y[i];
   mean/=(double)(n-maxd);

   for (i=maxd;i<n;i++) {
      ym=y[i]-mean;
      skew+=ym*ym*ym;
      autocorr+=ym*(y[i-d1]-mean)*(y[i-d2]-mean);
   }

   return (autocorr/skew);
}

int main(int argc,char *argv[])
{
   double *y=NULL;
   char file[MAXNAMELEN]; 
   long  col, delay1, delay2, ifi, n, nmax, nexcl;
   FILE *fpin=stdin;

   whatido(argc,argv,"three point autocorrelation",argv[0]);
   nmax=lcan(argc,argv,'l',DEF_POINTS);
   nexcl=lcan(argc,argv,'x',0);
   col=lcan(argc,argv,'c',0);
   delay1=lcan(argc,argv,'d',DEF_DELAY1);
   delay2=lcan(argc,argv,'D',2*delay1);

   /* process all files */
   for (ifi=1; ifi<=nstrings(); ifi++) {

      /* open file */
      col=getfile(file,nthstring(argc,argv,ifi),col);
      if(file[0]!='-'){
	 fpin=fopen(file,"r");
	 if(fpin==NULL) {
	    fprintf(stderr,"cannot open file %s\n",file);
	    return 5;
	 }
	 putcol(file,col);
      }
      else
         fpin=stdin;

      /* read file */
      for(n=0; n<nexcl; n++) if(readdum(fpin)==EOF) return 4;
      for(n=0; n<nmax; n++){
	 y=(double *)realloc(y, (size_t)((n+1)*sizeof(double)));
	 if(readval(fpin,y+n,col)==EOF) break;
      }

      if ((n<=delay1)||(n<=delay2)) {
	 fprintf(stderr,
		 "delay too large or time series too short for file %s\n",
		 file);
	 free(y);
	 return 5;
      };

      printf("%lf %s\n",third_auto(y,n,delay1,delay2),file);
   }

   free(y);
   return 0;
}
