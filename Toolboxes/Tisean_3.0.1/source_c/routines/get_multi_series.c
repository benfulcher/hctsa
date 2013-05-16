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
/*Author: Rainer Hegger Last modified: Sep 3, 1999 */
/*Note: Keep in mind that the first index runs the dimension,
  the second the time series index */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "tsa.h"
#include "tisean_cec.h"

#define SIZE_STEP 1000
extern void check_alloc(void*);

double **get_multi_series(char *name,unsigned long *l,unsigned long ex,
			  unsigned int *col,char *which,char colfix,
			  unsigned int verbosity)
{
  char *input,**format;
  int i,j;
  unsigned int *hcol,maxcol=0,colcount=0;
  unsigned long count,max_size=SIZE_STEP,hl,allcount;
  int input_size=INPUT_SIZE;
  double **x;
  FILE *fin;

  if (strlen(which) > 0) {
    colcount=1;
    for (i=0;i<strlen(which)-1;i++) {
      if (!isdigit((unsigned int)which[i]) && (which[i] != ',')) {
	fprintf(stderr,"Wrong format in the column string."
		" Has to be num,num,num,...,num\n");
	exit(GET_MULTI_SERIES_WRONG_TYPE_OF_C);
      }
      if (which[i] == ',') {
	colcount++;
	which[i]=' ';
      }
    }
    if (!isdigit((unsigned int)which[strlen(which)-1])) {
	fprintf(stderr,"Wrong format in the column string."
		" Has to be num,num,num,...,num\n");
	exit(GET_MULTI_SERIES_WRONG_TYPE_OF_C);
    }      
  }
  if (!colfix && (*col < colcount))
    *col=colcount;

  check_alloc(input=(char*)calloc((size_t)input_size,(size_t)1));
  check_alloc(hcol=(unsigned int*)malloc(sizeof(unsigned int)* *col));
  while ((int)(*which) && isspace((unsigned int)(*which)))
    which++;
  if (*which)
    for (i=0;i< *col-1;i++) {
      sscanf(which,"%u",&hcol[i]);
      if (hcol[i] > maxcol)
	maxcol=hcol[i];
      while ((int)(*which) && !isspace((unsigned int)(*which)))
	which++;
      while ((int)(*which) && isspace((unsigned int)(*which)))
	which++;
      if (!((int)(*which)))
	break;
    }
  else
    i= -1;
  
  if (*which)
    sscanf(which,"%u",&hcol[i]);
  else
    for (j=i+1;j< *col;j++)
      hcol[j]= ++maxcol;
  
  if (verbosity&VER_INPUT) {
    fprintf(stderr,"Using columns: ");
    for (i=0;i< *col;i++)
      fprintf(stderr,"%d ",hcol[i]);
    fprintf(stderr,"\n");
  }

  check_alloc(format=(char**)malloc(sizeof(char*)* *col));
  for (i=0;i< *col;i++) {
    check_alloc(format[i]=(char*)calloc((size_t)(4*hcol[i]),(size_t)1));
    strcpy(format[i],"");
    for (j=1;j<hcol[i];j++)
      strcat(format[i],"%*lf");
    strcat(format[i],"%lf");
  }
  free(hcol);
  
  check_alloc(x=(double**)malloc(sizeof(double*)* *col));
  for (i=0;i< *col;i++)
    check_alloc(x[i]=(double*)malloc(sizeof(double)*max_size));
  hl= *l;

  count=0;
  allcount=0;
  if (name == NULL) {
    for (i=0;i<ex;i++)
      if ((input=myfgets(input,&input_size,stdin,verbosity)) == NULL)
	break;
    while ((count < hl) && 
	   ((input=myfgets(input,&input_size,stdin,verbosity)) != NULL)) {
      if (count == max_size) {
	max_size += SIZE_STEP;
	for (i=0;i< *col;i++)
	  check_alloc(x[i]=(double*)realloc(x[i],sizeof(double)*max_size));
      }
      allcount++;
      for (i=0;i< *col;i++)
	if (sscanf(input,format[i],&x[i][count]) != 1) {
	  if (verbosity&VER_INPUT)
	    fprintf(stderr,"Line %lu ignored: %s",allcount,input);
	  break;
	}
      if (i == *col)
	count++;
    }
  }
  else {
    fin=fopen(name,"r");
    for (i=0;i<ex;i++)
      if ((input=myfgets(input,&input_size,fin,verbosity)) == NULL)
	break;
    while ((count < hl) && 
	   ((input=myfgets(input,&input_size,fin,verbosity)) != NULL)) {
      if (count == max_size) {
	max_size += SIZE_STEP;
	for (i=0;i< *col;i++)
	  check_alloc(x[i]=(double*)realloc(x[i],sizeof(double)*max_size));
      }
      allcount++;
      for (i=0;i< *col;i++)
	if (sscanf(input,format[i],&x[i][count]) != 1) {
	  if (verbosity&VER_INPUT)
	    fprintf(stderr,"Line %lu ignored: %s",allcount,input);
	  break;
	}
      if ((count == 0) && (i == *col) && (verbosity&VER_FIRST_LINE)) {
	fprintf(stderr,"get_multi_series: first data item(s) used:\n");
	for (i=0;i< *col;i++)
	  fprintf(stderr,"%lf ",x[i][0]);
	fprintf(stderr,"\n");
      }
      if (i == *col)
	count++;
    }
    fclose(fin);
  }
  
  for (i=0;i< *col;i++)
    free(format[i]);
  free(format);
  free(input);

  *l = count;  
  if (*l == 0) {
    fprintf(stderr,"0 lines read. It makes no sense to continue. Exiting!\n");
    exit(GET_MULTI_SERIES_NO_LINES);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Use %lu lines.\n",*l);
  }

  if (max_size > count)
    for (i=0;i< *col;i++) 
      check_alloc(x[i]=(double*)realloc(x[i],sizeof(double)*count));
  
  return x;
}
#undef SIZE_STEP
