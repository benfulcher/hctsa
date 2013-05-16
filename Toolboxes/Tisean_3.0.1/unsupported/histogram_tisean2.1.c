/*Author: Rainer Hegger. Last modified Nov 14, 2000*/
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Makes a histogram of the data"

unsigned long length=ULONG_MAX;
unsigned long base=50;
unsigned long exclude=0;
unsigned int column=1;
unsigned int verbosity=0xff;
double size;
char my_stdout=1,gotsize=0;
char *outfile=NULL;
char *infile=NULL;

double *series;
double average,var;
double min,max;
unsigned int bsize=1024,bsize1=1023;
long *box,*list;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
	  " possible datafile.\nIf no datafile is given stdin is read. "
	  " Just - also means stdin\n");
  fprintf(stderr,"\t-l length of file [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-c column to read [default %d]\n",column);
  fprintf(stderr,"\t-b # of bases [default %ld]\n",base);
  fprintf(stderr,"\t-r size of boxes [default (data interval)/50]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.dat ;"
	  " If no -o is given: stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **str)
{
  char *out;

  if ((out=check_option(str,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(str,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(str,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(str,n,'b','u')) != NULL)
    sscanf(out,"%lu",&base);
  if ((out=check_option(str,n,'r','f')) != NULL) {
    gotsize=1;
    sscanf(out,"%lf",&size);
  }
  if ((out=check_option(str,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(str,n,'o','o')) != NULL) {
    my_stdout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void put_in_boxes(void)
{
  unsigned int i,x;
  
  for (i=0;i<bsize;i++)
    box[i]= -1;

  for (i=0;i<length;i++) {
    x=(unsigned int)(series[i]/size)&bsize1;
    printf("%ld\n",x);
    list[i]=box[x];
    box[x]=i;
  }
}

unsigned long neighbors(double x)
{
  unsigned long found=0;
  int i,i1;
  long element;
  static double size2;
  
  size2=size/2.0;
  
  i=(int)((x+size/2.)/size)&bsize1;
  for (i1=i-1;i1<=i+1;i1++) {
    element=box[i1&bsize1];
    while (element != -1) {
      if (fabs(x-series[element]) <= size2)
	found++;
      element=list[element];
    }
  }
  return found;
}

int main(int argc,char **argv)
{
  char stdi=0;
  unsigned long i;
  double x,norm,step;
  FILE *fout;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,&column,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,1));
      strcpy(outfile,infile);
      strcat(outfile,".his");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,1));
      strcpy(outfile,"stdin.his");
    }
  }
  if (!my_stdout)
    test_outfile(outfile);

  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&average,&var);
  rescale_data(series,length,&min,&max);
  
  
  if (base > 0) {
    if (!gotsize)
      size=1.0/(double)base;
    else {
      size /= max;
    }
  }

  check_alloc(box=(long*)malloc(sizeof(long)*bsize));
  check_alloc(list=(long*)malloc(sizeof(long)*length));


  if (base > 0)
    put_in_boxes();
  
  step=1.0;
  if (base > 0)
    step=size;//1.0/base;
  norm=1.0/length/size/max;

  if (!my_stdout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(fout,"#interval of data: [%e:%e]\n",min,max+min);
    fprintf(fout,"#average= %e\n",average);
    fprintf(fout,"#standard deviation= %e\n",var);
    for (i=0;i<base;i++) {
      x=(double)(i*step);
      fprintf(fout,"%e %e\n",x*max+min,(double)(neighbors(x))*norm);
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#interval of data: [%e:%e]\n",min,max+min);
    fprintf(stdout,"#average= %e\n",average);
    fprintf(stdout,"#standard deviation= %e\n",var);
    for (i=0;i<base;i++) {
      x=(double)(i*step);
      fprintf(stdout,"%e %e\n",x*max+min,(double)(neighbors(x))*norm);
      fflush(stdout);
    }
  }
  return 0;
}
