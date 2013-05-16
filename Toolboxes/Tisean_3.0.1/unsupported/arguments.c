/* read command line options */
/* Copyright (C) Thomas Schreiber (1997) */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#define MAX(a,b) ((a)>(b)? (a) : (b))

void usage(char *);

long *largs=NULL;
long nargs;

void argdel(int argc, long i)
{
   if(largs==NULL){
      largs=calloc(argc, sizeof(long));
      nargs=argc-1;
   }
   if(i==0 || i>=argc) return;
   if(largs[i]) return;
   largs[i]=1;
   nargs-=1;
}

long nstrings()
{
   return nargs==0 ? 1: nargs;
}

char *nthstring(int argc, char *argv[], long n)
{
   long iv, i;

   for(iv=0, i=1; i<argc; i++){
      if (!largs[i]) iv++;
      if(iv==n) break;
   }
   return i>=argc? "-": argv[i];
}

long lopt(int argc, char *argv[], char c)
{
   long i;

   for (i=1; i<argc; i++) 
      if (argv[i][0]=='-' && argv[i][1]==c){
         argdel(argc,i);
         return 1;
      }
   return 0;
}  

long lmust(int argc, char *argv[], char c)
{
   long lval, i;

   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-' && argv[i][1]==c) {
         argdel(argc,i);
         if ((long)strlen(argv[i])==2) {
	    if (++i < argc)
	       if(sscanf(argv[i],"%ld",&lval)==1){
                  argdel(argc,i);
                  return lval;
	       }
         }
         else 
	    if(sscanf(argv[i]+2,"%ld",&lval)==1) return lval;
         usage(argv[0]);
      }
   }
   usage(argv[0]);
}  

double fmust(int argc, char *argv[], char c)
{
   long i;
   double fval;

   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-' && argv[i][1]==c) {
         argdel(argc,i);
         if ((long)strlen(argv[i])==2) {
	    if (++i < argc)
	       if(sscanf(argv[i],"%lf",&fval)==1){
                  argdel(argc,i);
                  return fval;
	       }
         }
         else 
	    if(sscanf(argv[i]+2,"%lf",&fval)==1) return fval;
         usage(argv[0]);
      }
   }
   usage(argv[0]);
}  

long lcan(int argc, char *argv[], char c, long ldef)
{
   long lval, i;

   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-' && argv[i][1]==c) {
         argdel(argc,i);
         if ((long)strlen(argv[i])==2) {
	    if (++i < argc)
	       if(sscanf(argv[i],"%ld",&lval)==1){
                  argdel(argc,i);
                  return lval;
	       }
         }
         else 
	    if(sscanf(argv[i]+2,"%ld",&lval)==1) return lval;
         return ldef;
      }
   }
   return ldef;
}  

char *scan(int argc, char *argv[], char c, char* sdef)
{
   long i;

   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-' && argv[i][1]==c) {
         argdel(argc,i);
         if ((long)strlen(argv[i])==2) {
            if (++i < argc) {
               if(argv[i][0]=='-') return sdef;
               argdel(argc,i);
               return argv[i];
	    }
         }
         else 
	    return argv[i][2]=='-'? sdef: argv[i]+2;
         return sdef;
      }
   }
   return sdef;
}  

double fcan(int argc, char *argv[], char c, double fdef)
{
   long i;
   double fval;

   for (i=1; i<argc; i++) {
      if (argv[i][0]=='-' && argv[i][1]==c) {
         argdel(argc,i);
         if ((long)strlen(argv[i])==2) {
	    if (++i < argc)
	       if(sscanf(argv[i],"%lf",&fval)==1){
                  argdel(argc,i);
                  return fval;
	       }
         }
         else 
	    if(sscanf(argv[i]+2,"%lf",&fval)==1) return fval;
         usage(argv[0]);
      }
   }
   return fdef;
}  

static long line_count;

long getfile(char *file, char *file_raw, long col)
{
   char *com;

   line_count=0;
   strcpy(file,file_raw);
   if(col==0){
      com=strrchr(file,(int)',');

      if(com!=NULL){
         if(strrchr(com,(int)'_')==NULL) 
            if(sscanf(com+1,"%ld",&col)==1) *com='\0';
      }
   }
   if(col>0) fprintf(stderr,"reading from column %ld\n",col);
   return col;
}

void putcol(char *file, long col)
{
   if(col>0) sprintf(file+strlen(file),",%ld",col);
}

long readline(FILE *fp, double *y, long col)
{
   long i, ic, icc;
   char c;

   if(col==0) col=1;
   for(i=1; i<=col; i++){
      ic=fscanf(fp,"%lf",y);
      if(ic==EOF) return ic;
      if(ic!=1) break;
   }
   line_count++;
   do
      icc=fscanf(fp,"%c",&c);
   while(c!='\n' && icc!=EOF);
   return ic;
}

long readdum(FILE *fp)
{
   long ic;
   char c;

   do
      ic=fscanf(fp,"%c",&c);
   while(c!='\n' && ic!=EOF);
   line_count++;
   return ic;
}

long readval(FILE *fp, double *y, long col)
{
   long ic;

   do{
      ic=readline(fp,y,col);
      if(ic!=EOF && ic!=1) 
         fprintf(stderr,"data in line %ld ignored\n", line_count);
   }
   while(ic!=EOF && ic!=1);
   return ic;
}

void whatido(int argc, char *argv[], char *text, char *progname)
{
   fprintf(stderr,"%s: %s\n", progname, text);
   if(lopt(argc,argv,'h')) usage(progname);
   argdel(argc,0);
}

char *getout(int argc, char *argv[], long *iget)
{
   char *fout;

   *iget=0;
   fout=scan(argc, argv,'o',"");
   if(fout[0]!='\0') return fout;
   *iget=lopt(argc,argv,'o');
   return fout;
}
