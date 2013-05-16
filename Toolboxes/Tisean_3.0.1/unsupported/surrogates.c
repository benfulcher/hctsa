/* Create Surrogate data */

/****************************************************************/
/* this source is freely distributable                          */
/* copyright (C) Andreas Schmitz (1997)                         */
/* Wuppertal University, Germany                                */
/****************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include "arguments.h"
#include "surro.h"
#include "dft.h"
#include "rank.h"
#include "rand.h"

#define MAXNAMELEN 512      /* maximal lenght of file name */
#define DEF_POINTS 10000000 /* default number of points */
#define DEF_ITER   100000   /* default iterations for polishing */
#define DEF_NUMB   1        /* default number of surrogates */

void usage(char *progname)
{
   fprintf(stderr,"\nUsage: %s [-n# -i# -I# -l# -x# -c# -h] file\n",progname);
   fprintf(stderr,"\t-n <number of surrogates (1)>\n");
   fprintf(stderr,"\t-i <number of iterations (until no change)>\n");
   fprintf(stderr,"\t-I <seed for random numbers (0 = call to time())>\n");
   fprintf(stderr,"\t-l <number of points (whole file)>\n");
   fprintf(stderr,"\t-x <number of values to be skipped (0)>\n");
   fprintf(stderr,"\t-c <column to be read (1 or file,#)>\n");
   fprintf(stderr,"\t-h <show this message>\n\n");
   exit(0);
}

int main(int argc,char *argv[])
{
   double *y=NULL, *z=NULL, *ampl_org, *arr_sort, sy, sy2, syz;
   char file[MAXNAMELEN], *outfile; 
   long number, col, iter, it, nmax, nexcl, n, i, seed=67654748, filen, iget;
   FILE *fpin=stdin, *fpout=stdout;

   whatido(argc,argv,"create iterative surrogate data",argv[0]);
   nmax=lcan(argc,argv,'l',DEF_POINTS);
   nexcl=lcan(argc,argv,'x',0);
   col=lcan(argc,argv,'c',0);
   number=lcan(argc,argv,'n',DEF_NUMB);
   iter=lcan(argc,argv,'i',DEF_ITER);
   seed=-lcan(argc,argv,'I',seed);
   if(seed==0) seed=time(NULL);
   rnd_init((unsigned long)seed);
   outfile=getout(argc,argv,&iget);

   col=getfile(file,nthstring(argc,argv,1),col);
   if(file[0]=='-')
      file[0]='\0';
   else {
      fprintf(stderr,"trying to open file %s\n", file);
      fpin=fopen(file,"r");
      if(fpin==NULL) {
	fprintf(stderr,"cannot open file %s\n",file);
	return 5;
      }
   }
   if(number>1) iget=1; 
   if(iget){
      if(file[0]=='\0') strcpy(file,"stdin");
      putcol(file,col);
      strcat(file, "_surr");
   }
   else
      file[0]='\0';

   if(outfile[0]!='\0') strcpy(file,outfile);
   filen=(long)strlen(file);

   for(n=0; n<nexcl; n++) if(readdum(fpin)==EOF) return 4;
   for(n=0; n<nmax; n++){
      y=(double *)realloc(y, (size_t)((n+1)*sizeof(double)));
      if(readval(fpin,y+n,col)==EOF) break;
   }
   fprintf(stderr,"read %ld values\n",n);

   nmax=strip_primes(n);   /* reduce filelength -> low prime factors */
   if (n!=nmax) fprintf(stderr,"shortened to length: %ld\n",nmax);

   arr_sort=(double *)malloc((size_t)(nmax*sizeof(double)));
   ampl_org=(double *)malloc((size_t)((nmax/2+1)*sizeof(double)));
   z=(double *)malloc(nmax*sizeof(double));
   if (arr_sort==NULL || ampl_org==NULL || z==NULL) {
      fprintf(stderr,"no memory\n");
      exit(5);
   }

   dft(y,nmax,ampl_org,NULL); /* original amplitudes */

   for (i=0; i<nmax; i++)    /* make a sorted copy of the org data */
      arr_sort[i]=y[i];   
   sort(nmax,arr_sort);
      
   for (i=0; i<number; i++) {
      if(file[0]!='\0'){
         if(number>1) sprintf(file+filen,"_%3.3ld",i+1);
         fprintf(stderr,"%s ",file);
         fpout=fopen(file,"w");
         if(fpout==NULL) {
	   fprintf(stderr,"cannot open file %s\n",file);
	   return 5;
         }
      }
      it=pol_surr(nmax,y,z,arr_sort,ampl_org,iter);
      for(n=0, sy=0, sy2=0, syz=0; n<nmax; n++){
	 sy+=y[n];
         sy2+=y[n]*y[n];
         syz+=(y[n]-z[n])*(y[n]-z[n]);
      }
      sy=sqrt((sy2/nmax)-(sy/nmax)*(sy/nmax));
      syz=sqrt(syz/nmax);
      fprintf(stderr,"(%ld iterations, relative discrepancy  %lf)\n", 
         it,syz/sy);
      for(n=0; n<nmax; n++)
         fprintf(fpout,"%lf %lf\n",y[n],z[n]);
      fclose(fpout);
   }
   
   free(y);
   free(arr_sort);

   return 0;
}
