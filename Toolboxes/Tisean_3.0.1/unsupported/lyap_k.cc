//Author: Rainer Hegger. Last modified: Apr 23, 1999
#include <fstream.h>
#include <strstream.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsa.h"

#define BOX 128
const unsigned int ibox=BOX-1;

unsigned long length=ULONG_MAX;
unsigned long exclude=0;
unsigned long reference=ULONG_MAX;
unsigned int maxdim=2;
unsigned int mindim=2;
unsigned int delay=1;
unsigned int column=1;
unsigned int epscount=5;
unsigned int maxiter=50;
unsigned int window=0;
double epsmin=1.e-3,epsmax=1.e-2;
char eps0set=0,eps1set=0;
char *outfile=NULL;
char *infile=NULL;

double *series,**lyap;
long box[BOX][BOX],*liste,**lfound,*found,**count;
double max,min;

void show_options(char *progname)
{
  cerr << "Estimate the largest Lyapunov exponent using the Kantz algorithm\n";
  cerr << "\nUsage: " << progname << " [options]\n";
  cerr << "Options:\n";
  cerr << "Everything not being a valid option will be interpreted as a"
    " possible datafile.\nIf no datafile is given stdin is read. Just - also"
    " means stdin\n";
  cerr << "  -l # of data [default: whole file]\n";
  cerr << "  -x # of lines to be ignored [default: 0]\n";
  cerr << "  -c column to read [default: 1]\n";
  cerr << "  -M maxdim [default: 2]\n";
  cerr << "  -m mindim [default: 2]\n";
  cerr << "  -d delay [default: 1]\n";
  cerr << "  -r mineps [default: (data interval)/1000]\n";
  cerr << "  -R maxeps [default: (data interval)/100]\n";
  cerr << "  -# # of eps [default: 5]\n";
  cerr << "  -n # of reference points [default: # of data]\n";
  cerr << "  -s # of iterations [default: 50]\n";
  cerr << "  -t time window [default: 0]\n";
  cerr << "  -o outfile [default: 'datafile'.lyap]\n";
  cerr << "  -h show these options\n";
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
  if ((out=check_option(str,n,'M','u')) != NULL)
    sscanf(out,"%u",&maxdim);
  if ((out=check_option(str,n,'m','u')) != NULL)
    sscanf(out,"%u",&mindim);
  if ((out=check_option(str,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(str,n,'r','f')) != NULL) {
    eps0set=1;
    sscanf(out,"%lf",&epsmin);
  }
  if ((out=check_option(str,n,'R','f')) != NULL) {
    eps1set=1;
    sscanf(out,"%lf",&epsmax);
  }
  if ((out=check_option(str,n,'#','u')) != NULL)
    sscanf(out,"%u",&epscount);
  if ((out=check_option(str,n,'n','u')) != NULL)
    sscanf(out,"%lu",&reference);
  if ((out=check_option(str,n,'s','u')) != NULL)
    sscanf(out,"%u",&maxiter);
  if ((out=check_option(str,n,'t','u')) != NULL)
    sscanf(out,"%u",&window);
  if ((out=check_option(str,n,'o','o')) != NULL)
    if (strlen(out) > 0)
      outfile=out;
}

void put_in_boxes(double eps)
{
  unsigned long i;
  long j,k;
  static unsigned long blength=length-(maxdim-1)*delay-maxiter;

  for (i=0;i<BOX;i++)
    for (j=0;j<BOX;j++)
      box[i][j]= -1;

  for (i=0;i<blength;i++) {
    j=(long)(series[i]/eps)&ibox;
    k=(long)(series[i+delay]/eps)&ibox;
    liste[i]=box[j][k];
    box[j][k]=i;
  }
}

void lfind_neighbors(long act,double eps)
{
  unsigned int hi,k,k1;
  long i,j,i1,i2,j1,element,lwindow;
  double dx,eps2=sqr(eps);

  lwindow=(long)window;
  for (hi=0;hi<maxdim-1;hi++)
    found[hi]=0;
  i=(long)(series[act]/eps)&ibox;
  j=(long)(series[act+delay]/eps)&ibox;
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&ibox;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&ibox];
      while (element != -1) {
	if ((element < (act-lwindow)) || (element > (act+lwindow))) {
	  dx=sqr(series[act]-series[element]);
	  if (dx <= eps2) {
	    for (k=1;k<maxdim;k++) {
	      k1=k*delay;
	      dx += sqr(series[act+k1]-series[element+k1]);
	      if (dx <= eps2) {
		k1=k-1;
		lfound[k1][found[k1]]=element;
		found[k1]++;
	      }
	      else
		break;
	    }
	  }
	}
	element=liste[element];
      }
    }
  }
}

void iterate_points(long act)
{
  double **lfactor;
  double *dx;
  unsigned int i,j,l,l1;
  long k,element,**lcount;
  
  lfactor=new double*[maxdim-1];
  lcount=new long*[maxdim-1];
  for (i=0;i<maxdim-1;i++) {
    lfactor[i]=new double[maxiter+1];
    lcount[i]=new long[maxiter+1];
  }
  dx=new double[maxiter+1];

  for (i=0;i<=maxiter;i++)
    for (j=0;j<maxdim-1;j++) {
      lfactor[j][i]=0.0;
      lcount[j][i]=0;
    }

  for (j=mindim-2;j<maxdim-1;j++) {
    for (k=0;k<found[j];k++) {
      element=lfound[j][k];
      for (i=0;i<=maxiter;i++)
	dx[i]=sqr(series[act+i]-series[element+i]);
      for (l=1;l<j+2;l++) {
	l1=l*delay;
	for (i=0;i<=maxiter;i++)
	  dx[i] += sqr(series[act+i+l1]-series[element+l1+i]);
      }
      for (i=0;i<=maxiter;i++)
	if (dx[i] > 0.0){
	  lcount[j][i]++;
	  lfactor[j][i] += dx[i];
	}
    }
  }
  for (i=mindim-2;i<maxdim-1;i++)
    for (j=0;j<=maxiter;j++)
      if (lcount[i][j]) {
	count[i][j]++;
	lyap[i][j] += log(lfactor[i][j]/lcount[i][j])/2.0;
      }
  
  for (i=0;i<maxdim-1;i++){
    delete[] lfactor[i];
    delete[] lcount[i];
  }
  delete[] lcount;
  delete[] lfactor;
  delete[] dx;
}

int main(int argc,char **argv)
{
  char stdi=0;
  double eps_fak;
  double epsilon;
  unsigned int i,j,l;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
  infile=search_datafile(argc,argv,&column);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL)
    if (!stdi) {
      outfile=new char[strlen(infile)+6];
      sprintf(outfile,"%s.lyap",infile);
    }
    else {
      outfile=new char[11];
      sprintf(outfile,"stdin.lyap");
    }
  test_outfile(outfile);

  series=get_series(infile,&length,exclude,column);
  rescale_data(series,length,&min,&max);

  if (eps0set)
    epsmin /= max;
  if (eps1set)
    epsmax /= max;

  if (epsmin >= epsmax) {
    epsmax=epsmin;
    epscount=1;
  }
  
  if (reference > (length-maxiter-(maxdim-1)*delay))
    reference=length-maxiter-(maxdim-1)*delay;
  if (maxdim < 2)
    maxdim=2;
  if (mindim < 2)
    mindim=2;
  if (mindim > maxdim)
    maxdim=mindim;
  
  liste=new long[length];
  found=new long[maxdim-1];
  lfound=new long*[maxdim-1];
  for (i=0;i<maxdim-1;i++)
    lfound[i]=new long[length];
  count=new long*[maxdim-1];
  for (i=0;i<maxdim-1;i++)
    count[i]=new long[maxiter+1];
  lyap=new double*[maxdim-1];
  for (i=0;i<maxdim-1;i++)
    lyap[i]=new double[maxiter+1];

  if (epscount == 1)
    eps_fak=1.0;
  else
    eps_fak=pow(epsmax/epsmin,1.0/(double)(epscount-1));

  ofstream fout(outfile);
  for (l=0;l<epscount;l++) {
    epsilon=epsmin*pow(eps_fak,(double)l);
    for (i=0;i<maxdim-1;i++)
      for (j=0;j<=maxiter;j++) {
	count[i][j]=0;
	lyap[i][j]=0.0;
      }
    put_in_boxes(epsilon);
    for (i=0;i<reference;i++) {
      lfind_neighbors(i,epsilon);
      iterate_points(i);
    }
    cerr.flags(ios::scientific);
    cerr << "epsilon= " << epsilon*max << '\n';
    for (i=mindim-2;i<maxdim-1;i++) {
      fout.flags(ios::scientific);
      fout << "#epsilon= " << epsilon*max << " dim= " << i+2 << '\n';
      for (j=0;j<=maxiter;j++)
	if (count[i][j])
	  fout << j << " " << lyap[i][j]/count[i][j] << " " 
	       << count[i][j] <<'\n';
      fout << '\n';
    }
    fout.flush();
  }
  fout.close();
  return 0;
}
