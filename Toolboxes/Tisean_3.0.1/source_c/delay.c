/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Rainer Hegger. Last modified (rewritten in C) Aug 22, 2004*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include "routines/tsa.h"

#define WID_STR "Produces delay vectors"


unsigned long length=ULONG_MAX;
unsigned long exclude=0;
unsigned int verbosity=0xff;
int delay=1;
unsigned int indim=1,embdim=2;
char *column=NULL,*format=NULL,*multidelay=NULL;
char *outfile=NULL;
char *infile=NULL;
char dimset=0,formatset=0,embset=0,mdelayset=0,delayset=0;
char stdo=1;

double **series;
unsigned int *formatlist,*delaylist;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr,"\nUsage: %s [options]\n",progname);
  fprintf(stderr,"Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
	  " possible datafile.\nIf no datafile is given stdin is read."
	  " Just - also means stdin\n");
  fprintf(stderr,"\t-l # of data [default: whole file]\n");
  fprintf(stderr,"\t-x # of rows to ignore [default: 0]\n");
  fprintf(stderr,"\t-M num. of columns to read [default: %u]\n",indim);
  fprintf(stderr,"\t-c columns to read [default: 1,...,M]\n");
  fprintf(stderr,"\t-m dimension [default: 2]\n");
  fprintf(stderr,"\t-F format of the delay vector (see man page)\n");
  fprintf(stderr,"\t-d delay [default: 1]\n");
  fprintf(stderr,"\t-D multi delay list (see man page)\n");
  fprintf(stderr,"\t-V verbosity level [default: 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-o output file [default: 'datafile'.del, "
	  "without -o: stdout]\n");
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
  if ((out=check_option(str,n,'c','s')) != NULL)
    column=out;
  if ((out=check_option(str,n,'M','u')) != NULL) {
    sscanf(out,"%u",&indim);
    dimset=1;
  }
  if ((out=check_option(str,n,'F','s')) != NULL) {
    format=out;
    formatset=1;
  }
  if ((out=check_option(str,n,'m','u')) != NULL) {
    sscanf(out,"%u",&embdim);
    embset=1;
  }
  if ((out=check_option(str,n,'d','u')) != NULL) {
    sscanf(out,"%u",&delay);
    delayset=1;
  }
  if ((out=check_option(str,n,'D','s')) != NULL) {
    multidelay=out;
    mdelayset=1;
  }
  if ((out=check_option(str,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(str,n,'o','o')) != NULL) {
    stdo=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

void create_format_list(void)
{
  unsigned int i=0,num=0,sum=0;

  while (format[i]) {
    if (!(isdigit(format[i])) && !(format[i] == ',')) {
      fprintf(stderr,"Wrong format of -F parameter. Exiting!\n");
      exit(DELAY_WRONG_FORMAT_F);
    }
    i++;
  }

  i=0;
  while (format[i]) {
    if (format[i++] == ',')
      num++;
  }

  check_alloc(formatlist=(unsigned int*)malloc(sizeof(int)*(num+1)));
  for (i=0;i<=num;i++) {
    sscanf(format,"%d",&formatlist[i]);
    if (i<num) {
      while ((*format) != ',')
	format++;
    }
    format++;
  }

  if (dimset && ((num+1) != indim)) {
    fprintf(stderr,"Number of dimensions in -F is not equal to -M. Exiting!\n");
    exit(DELAY_DIM_NOT_EQUAL_F_M);
  }

  for (i=0;i<=num;i++)
    sum += formatlist[i];
  if (embset && (sum != embdim)) {
    fprintf(stderr,"The dimensions given in -m and -F are not equal!"
	    " Exiting\n");
    exit(DELAY_DIM_NOT_EQUAL_F_m);
  }
  if (!dimset)
    indim=num+1;
  if (!embset)
    embdim=sum;
}

void create_delay_list(void)
{
  unsigned int i=0,num=0;

  while (multidelay[i]) {
    if (!(isdigit(multidelay[i])) && !(multidelay[i] == ',')) {
      fprintf(stderr,"Wrong format of -D parameter. Exiting!\n");
      exit(DELAY_WRONG_FORMAT_D);
    }
    i++;
  }

  i=0;
  while (multidelay[i]) {
    if (multidelay[i++] == ',')
      num++;
  }

  check_alloc(delaylist=(unsigned int*)malloc(sizeof(int)*(num+1)));
  for (i=0;i<=num;i++) {
    sscanf(multidelay,"%d",&delaylist[i]);
    if (i<num) {
      while ((*multidelay) != ',')
	multidelay++;
    }
    multidelay++;
  }

  if ((num+1) != (embdim-indim)) {
    fprintf(stderr,"Wrong number of delays. See man page. Exiting!\n");
    exit(DELAY_WRONG_NUM_D);
  }
}

int main(int argc,char **argv)
{
  char stin=0;
  unsigned long i;
  int j,k;
  unsigned int alldim,maxemb,emb,rundel,delsum,runmdel;
  unsigned int *inddelay;
  FILE *fout;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif


  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stin=1;

  if (outfile == NULL) {
    if (!stin) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,1));
      strcpy(outfile,infile);
      strcat(outfile,".del");
    }
    else {
      check_alloc(outfile=(char*)calloc(10,1));
      strcpy(outfile,"stdin.del");
    }
  }
  if (!stdo)
    test_outfile(outfile);

  if (delayset && mdelayset) {
    fprintf(stderr,"-d and -D can't be used simultaneously. Exiting!\n");
    exit(DELAY_INCONS_d_D);
  }

  if (delay < 1) {
    fprintf(stderr,"Delay has to be larger than 0. Exiting!\n");
    exit(DELAY_SMALL_ZERO);
  }

  if (!formatset && (embdim%indim)) {
    fprintf(stderr,"Inconsistent -m and -M. Please set -F\n");
    exit(DELAY_INCONS_m_M);
  }
  if (formatset) {
    create_format_list();
  }
  else {
    check_alloc(formatlist=(unsigned int*)malloc(sizeof(int)*indim));
    for (i=0;i<indim;i++) {
      formatlist[i]=embdim/indim;
    }
  }

  alldim=0;
  for (i=0;i<indim;i++)
    alldim += formatlist[i];

  if (mdelayset) {
    create_delay_list();
  }

  check_alloc(inddelay=(unsigned int*)malloc(sizeof(int)*alldim));

  rundel=0;
  if (!mdelayset) {
    for (i=0;i<indim;i++) {
      delsum=0;
      inddelay[rundel++]=delsum;
      for (j=1;j<formatlist[i];j++) {
	delsum += delay;
	inddelay[rundel++]=delsum;
      }
    }
  }
  else {
    runmdel=0;
    for (i=0;i<indim;i++) {
      delsum=0;
      inddelay[rundel++]=delsum;
      for (j=1;j<formatlist[i];j++) {
	delsum += delaylist[runmdel++];
	inddelay[rundel++]=delsum;
      }
    }
  }

  maxemb=0;
  for (i=0;i<alldim;i++)
    maxemb=(maxemb<inddelay[i])?inddelay[i]:maxemb;

  if (column == NULL) {
    series=get_multi_series(infile,&length,exclude,&indim,"",dimset,verbosity);
  }
  else {
    series=get_multi_series(infile,&length,exclude,&indim,column,dimset,
			    verbosity);
  }

  if (stdo) {
    if (verbosity)
      fprintf(stderr,"Writing to stdout\n");
    for (i=maxemb;i<length;i++) {
      rundel=0;
      for (j=0;j<indim;j++) {
	emb=formatlist[j];
	for (k=0;k<emb;k++)
	  fprintf(stdout,"%e ",series[j][i-inddelay[rundel++]]);
      }
      fprintf(stdout,"\n");
    }
  }
  else {
    fout=fopen(outfile,"w");
    if (verbosity)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=maxemb;i<length;i++) {
      for (j=0;j<indim;j++) {
	rundel=0;
	emb=formatlist[j];
	for (k=0;k<emb;k++)
	  fprintf(fout,"%e ",series[j][i-inddelay[rundel++]]);
      }
      fprintf(fout,"\n");
    }
    fclose(fout);
  }

  if (formatlist != NULL)
    free(formatlist);
  if (delaylist != NULL)
    free(delaylist);
  free(inddelay);

  return 0;
}
