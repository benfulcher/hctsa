/* nonlinear noise reduction */
/* Copyright (C) Marcus Richter and Thomas Schreiber (1997) */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stddef.h>

#define MAXNAMELEN 512      /* maximal lenght of file name */
#define DEF_STACK  500      /* default stack depth */
#define DEF_POINTS 10000000 /* default number of points */

#include "arguments.h"
#include "eigen.h"
#include "clean.h"

struct list *top, *akt, *postakt;

void usage(char *progname)
{
   fprintf(stderr, 
     "\nUsage: %s -m# -q# -r# -K# [-d# -k# -T# -w# -s# -l# -x# -c# -h] file\n",
      progname);
   fprintf(stderr,"\t-m <embedding dimension>\n");
   fprintf(stderr,"\t-q <dimension of manifold projected on>\n");
   fprintf(stderr,"\t-r <diameter of neighbourhood>\n");
   fprintf(stderr,"\t-K <maximal number of neighbours>\n");
   fprintf(stderr,"\t-d <delay (1)>\n");
   fprintf(stderr,"\t-k <minimal number of neighbours (m)>\n");
   fprintf(stderr,"\t-T <maximal time in the past considered as neighbours (all)>\n");
   fprintf(stderr,"\t-w <range represented by each centre (e/4)>\n");
   fprintf(stderr,"\t-s <stack depth (%d)>\n",DEF_STACK);
   fprintf(stderr,"\t-l <number of points (whole file)>\n");
   fprintf(stderr,"\t-x <number of values to be skipped (0)>\n");
   fprintf(stderr,"\t-c <column to be read (1 or file,#)>\n");
   fprintf(stderr,"\t-h <show this message>\n\n");
   exit(0);
}

int main(int argc, char *argv[])
{
   double *y=NULL, *yc=NULL, eps, delta, *r;
   long col, m, d=1, nmax, nexcl, i, n, kmin,nq,kmax,npmax=0,stdpm=0,stdp,iget;
   FILE *fpin=stdin, *fpout=stdout;
   char file[MAXNAMELEN], *outfile;

   whatido(argc,argv,"nonlinear noise reduction in a data stream",argv[0]);
   nmax=lcan(argc,argv,'l',DEF_POINTS);
   nexcl=lcan(argc,argv,'x',0);
   col=lcan(argc,argv,'c',0);
   m=lmust(argc,argv,'m');
   nq=lmust(argc,argv,'q');
   eps=fmust(argc,argv,'r');
   d=lcan(argc,argv,'d',d);
   kmax=lmust(argc,argv,'K');
   kmin=lcan(argc,argv,'k',m);
   npmax=lcan(argc,argv,'T',npmax);
   delta=fcan(argc,argv,'w',eps/4);
   stdp=lcan(argc,argv,'s',DEF_STACK);
   outfile=getout(argc,argv,&iget);

   col=getfile(file,nthstring(argc,argv,1),col);

   if(file[0]=='-') 
      file[0]='\0';
   else
      fpin=fopen(file,"r");

   if(iget){
      if(file[0]=='\0') strcpy(file,"stdin");
      putcol(file,col);
      strcat(file, "_c" );
   }
   if(outfile[0]!='\0') strcpy(file,outfile);

   if(iget || outfile[0]!='\0'){
      fpout=fopen(file,"w");
      if(fpout==NULL) {
        fprintf(stderr,"cannot open file %s\n",file);
        return 5;
      }
      fprintf(stderr,"writing to %s\n",file);
   }

   top=NULL;
   r=(double *)malloc(m*sizeof(double));
   for(i=0; i<m; i++)                /* $r=1/\sqrt{P},\;{\rm tr}(1/r)=1$ */
      r[i]= (i==0 || i==m-1)? (2*SMALL+m-2)/SMALL: 2*SMALL+m-2;

   for(n=0; n<nexcl; n++) if(readdum(fpin)==EOF) return 4;
   for(n=0; n<nmax; n++){
      y=(double *)realloc(y, (n+1)*sizeof(double));
      yc=(double *)realloc(yc, (n+1)*sizeof(double));
      if(readval(fpin,y+n,col)==EOF) break;
      yc[n]=y[n];
      if(n>=(m-1)*d){
	clean(n,y,yc,d,m,nq,kmin,kmax,npmax,eps,delta,r,stdp,&stdpm); 
        fprintf(fpout,"%lf %lf\n",yc[n-(m-1)*d],y[n-(m-1)*d]-yc[n-(m-1)*d]);
      }
   }
   nmax=n;
   for(n=nmax-(m-1)*d; n<nmax;n++)
      fprintf(fpout,"%lf %lf\n",yc[n],y[n]-yc[n]);

   for(n=0, eps=0; n<nmax; n++)
      eps+=(yc[n]-y[n])*(yc[n]-y[n]);
   if(n>0) {
      fprintf(stderr,"rms correction: %lf\n", sqrt(eps/nmax));
      fprintf(stderr,"average stack usage: %lf\n", (double)stdpm/nmax);
   }
   fclose(fpout);
}



