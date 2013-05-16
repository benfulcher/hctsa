/*Author: Rainer Hegger. Last modified Apr 28, 1999*/
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsa.h"

unsigned long length=ULONG_MAX;
unsigned long base=50;
unsigned long exclude=0;
unsigned int column=1;
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
  fprintf(stderr,"\nMakes a histogram of the data\n\n");
  fprintf(stderr," Usage: " << progname << " [options]\n");
  fprintf(stderr," options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
	  " possible datafile.\nIf no datafile is given stdin is read. "
	  " Just - also means stdin\n");
  fprintf(stderr,"\t-l length of file [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-c column to read [default is %d]\n",column);
  fprintf(stderr,"\t-b # of bases [default %d]\n",base);
  fprintf(stderr,"\t-e size of boxes (in rescaled ([0:1]) units) "
	  " [default 2/bases]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.dat ;"
	  " If no -o is given: stdout]\n");
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
  if ((out=check_option(str,n,'e','f')) != NULL) {
    gotsize=1;
    sscanf(out,"%lf",&size);
  }
  if ((out=check_option(str,n,'o','o')) != NULL) {
    my_stdout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
  if (size > 0.5) {
    fprintf(stderr,"Size must be smaller 0.5");
    exit(127);
  }
}

void put_in_boxes(void)
{
  unsigned int i,x;
  
  for (i=0;i<bsize;i++)
    box[i]= -1;

  for (i=0;i<length;i++) {
    x=(unsigned int)(series[i]/size)&bsize1;
    list[i]=box[x];
    box[x]=i;
  }
}

unsigned long neighbors(double x)
{
  unsigned long found=0;
  int i,i1;
  long element;
  static double size2=size/2.0;
  
  i=(int)(x/size)&bsize1;
  for (i1=i-1;i1<=i+1;i1++) {
    element=box[i1&bsize1];
    while (element != -1) {
      if (fabs(x-series[element]) < size2)
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
  
  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
  infile=search_datafile(argc,argv,&column);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL)
    if (!stdi) {
      outfile=new char[strlen(infile)+5];
      strcpy(outfile,infile);
      strcat(outfile,".his");
    }
    else {
      outfile=new char[10];
      strcpy(outfile,"stdin.his");
    }
  if (!my_stdout)
    test_outfile(outfile);

  if (base > 0)
    if (!gotsize)
      size=1.0/base;

  series=(double*)get_series(infile,&length,exclude,column);
  variance(series,length,&average,&var);
  rescale_data(series,length,&min,&max);

  box=new long[bsize];
  list=new long[length];


  if (base > 0)
    put_in_boxes();
  
  double step=1.0;
  if (base > 0)
    step=1.0/base;
  double x;
  double norm=(double)base/length/max;

  if (!my_stdout) {
    ofstream fout(outfile,ios::out);
    fout << "#interval of data: [" << min << ":" << max+min << "]\n";
    fout << "#average= " << average << '\n';
    fout << "#standard deviation= " << var << '\n';
    for (i=0;i<base;i++) {
      x=(double)(i*step);
      fout << x*max+min << " " << double(neighbors(x))*norm << '\n';
    }
    fout.close();
  }
  else {
    cout << "#interval of data: [" << min << ":" << max+min << "]\n";
    cout << "#average= " << average << '\n';
    cout << "#standard deviation= " << var << '\n';
    for (i=0;i<base;i++) {
      x=(double)(i*step);
      cout << x*max+min << " " << double(neighbors(x))*norm << '\n';
      cout.flush();
    }
  }
  return 0;
}
