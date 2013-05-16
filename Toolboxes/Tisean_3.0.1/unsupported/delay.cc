//Author: Rainer Hegger. Last modified Mar 20, 1999
#include <fstream.h>
#include <strstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tsa.h"

unsigned long length=ULONG_MAX;
unsigned long exclude=0;
unsigned int column=1;
int delay=1;
int dim=2;
char *outfile=NULL;
char *infile=NULL;
char stdo=1;

double *series;

void show_options(char *progname)
{
  cerr << "\nCreates a delay repesentation of the data\n\n";
  cerr << "\nUsage: " << progname << " [options]\n";
  cerr << "Options:\n";
  cerr << "Everything not being a valid option will be interpreted as a"
    " possible datafile.\nIf no datafile is given stdin is read. Just - also"
    " means stdin\n";
  cerr << "\t-l # of data [default: whole file]\n";
  cerr << "\t-x # of rows to ignore [default: 0]\n";
  cerr << "\t-c column to read [default: 1]\n";
  cerr << "\t-m dimension [default: 2]\n";
  cerr << "\t-d delay [default: 1]\n";
  cerr << "\t-o output file [default: 'datafile'.del, without -o: stdout]\n";
  cerr << "\t-h show these options\n";
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
  if ((out=check_option(str,n,'m','u')) != NULL)
    sscanf(out,"%u",&dim);
  if ((out=check_option(str,n,'d','u')) != NULL)
    sscanf(out,"%u",&delay);
  if ((out=check_option(str,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stin=0;
  unsigned long i;
  int j;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);

  infile=search_datafile(argc,argv,&column);
  if (infile == NULL)
    stin=1;

  if (outfile == NULL) 
    if (!stin) {
      outfile=new char[strlen(infile)+5];
      strcpy(outfile,infile);
      strcat(outfile,".del");
    }
    else {
      outfile=new char[10];
      strcpy(outfile,"stdin.del");
    }
  if (!stdo)
    test_outfile(outfile);

  if (delay < 1)
    delay=1;
  if (dim < 2)
    dim=2;
  if (column < 1)
    column=1;
  
  series=get_series(infile,&length,exclude,column);
  
  if (stdo) {
    for (i=(dim-1)*delay;i<length;i++) {
      for (j=0;j<dim-1;j++)
	cout << series[i-(dim-1-j)*delay] << " ";
      cout << series[i] << '\n';
    }
  }
  else {
    ofstream fout(outfile);
    for (i=(dim-1)*delay;i<length;i++) {
      for (j=0;j<dim-1;j++)
	fout << series[i-(dim-1-j)*delay] << " ";
      fout << series[i] << '\n';
    }
    fout.close();
  }
  return 0;
}
