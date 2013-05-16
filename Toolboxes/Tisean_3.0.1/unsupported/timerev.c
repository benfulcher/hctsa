/* Calculate statistics for time reversibility */

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
#define DEF_DELAY 1         /* default delay */

void usage(char *progname)
{
   fprintf(stderr,"\nUsage: %s [-d# -l# -x# -c# -h] file(s)\n",progname);
   fprintf(stderr,"\t-d <delay (1)>\n");
   fprintf(stderr,"\t-l <number of points (whole file)>\n");
   fprintf(stderr,"\t-x <number of values to be skipped (0)>\n");
   fprintf(stderr,"\t-c <column to be read (1 or file,#)>\n");
   fprintf(stderr,"\t-h <show this message>\n\n");
   exit(0);
}

double t_rev(double *y,unsigned long nmax,long del)
{
   double squared=0.0,cubed=0.0,diff;
   unsigned long i;

   for (i=0;i<nmax-del;i++) {
      diff=y[i+del]-y[i];
      squared+=diff*diff;
      cubed+=diff*diff*diff;
   }

   return (cubed/squared);
}

int main(int argc,char *argv[])
{
   double *y=NULL;
   char file[MAXNAMELEN]; 
   long  col, delay, n, nmax, nexcl, ifi;
   FILE *fpin=stdin;

   whatido(argc,argv,"time reversal asymmetry", argv[0]);
   nmax=lcan(argc,argv,'l',DEF_POINTS);
   nexcl=lcan(argc,argv,'x',0);
   col=lcan(argc,argv,'c',0);
   delay=lcan(argc,argv,'d',DEF_DELAY);

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

      if (n<=delay) {
	 fprintf(stderr,
		 "delay too large or time series too short for file %s\n",
		 file);
	 free(y);
	 return 5;
      };

      printf("%lf %s\n",t_rev(y,n,delay),file);
      fclose(fpin);
   }
	 
   free(y);
   return 0;
}
