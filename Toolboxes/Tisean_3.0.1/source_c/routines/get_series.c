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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsa.h"
#include "tisean_cec.h"

#define SIZE_STEP 1000
extern void check_alloc(void*);

double *get_series(char *name,unsigned long *l,unsigned long ex,
		unsigned int col,unsigned int verbosity)
{
  char *input,*format;
  int i;
  unsigned long count,allcount,max_size=SIZE_STEP,hl;
  int input_size=INPUT_SIZE;
  double *x;
  FILE *fin;
  
  check_alloc(input=(char*)calloc((size_t)input_size,(size_t)1));
  check_alloc(format=(char*)calloc((size_t)(4*col),(size_t)1));
  strcpy(format,"");
  for (i=1;i<col;i++)
    strcat(format,"%*lf");
  strcat(format,"%lf");
  
  check_alloc(x=(double*)malloc(sizeof(double)*max_size));
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
	check_alloc(x=(double*)realloc(x,sizeof(double)*max_size));
      }
      allcount++;
      if (sscanf(input,format,&x[count]) != 1) {
	if (verbosity&VER_INPUT)
	  fprintf(stderr,"Line %lu ignored: %s",allcount,input);
      }
      else
	count++;
      if ((verbosity&VER_FIRST_LINE) && (count == 0))
	fprintf(stderr,"get_series: first data item used:\n%lf\n",x[0]);
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
	check_alloc(x=(double*)realloc(x,sizeof(double)*max_size));
      }
      allcount++;
      if (sscanf(input,format,&x[count]) != 1) {
	if (verbosity&VER_INPUT)
	  fprintf(stderr,"Line %lu ignored: %s",allcount,input);
      }
      else
	count++;
    }
    fclose(fin);
  }
  free(input);
  
  *l = count;
  if (*l == 0) {
    fprintf(stderr,"0 lines read. It makes no sense to continue. Exiting!\n");
    exit(GET_SERIES_NO_LINES);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Use %lu lines.\n",*l);
  }
  if (max_size > count)
    check_alloc(x=(double*)realloc(x,sizeof(double)*count));
  
  return x;
}
#undef SIZE_STEP
