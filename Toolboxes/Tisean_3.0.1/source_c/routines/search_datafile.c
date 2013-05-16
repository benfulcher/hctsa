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
#include <ctype.h>
#include "tsa.h"

char check_col(char *col)
{
  int i;
  
  for (i=0;i<strlen(col);i++)
    if (!isdigit((unsigned int)col[i])) {
      fprintf(stderr,"Column must be a unsigned integer. Ignoring it!\n");
      return 0;
    }
  return 1;
}

char look_for_column(char *name,unsigned int *col)
{
  char *hcol,*hname;
  char vcol=0;
  int j,in;

  check_alloc(hname=(char*)calloc(strlen(name)+1,1));
  check_alloc(hcol=(char*)calloc(strlen(name)+1,1));
  j=0;
  while (*(name+j) != '\0') {
    if (*(name+j) == ',') {
      in=sscanf(name+j+1,"%s",hcol);
      if (in > 0)
	vcol=check_col(hcol);
      *(name+j)='\0';
      break;
    }
    *(hname+j)=*(name+j);
    j++;
  }
  *col=(unsigned int)atoi(hcol);
  free(hname);
  free(hcol);

  return vcol;
}

char* search_datafile(int n,char **names,unsigned int *col,
		      unsigned int verbosity)
{
  char valid=0,validcol=0;
  char *retname=NULL;
  int i;
  unsigned int hcol;
  FILE *test;

  for (i=n-1;i>0;i--) {
    if (names[i] != NULL) {
      valid=0;
      if (strcmp(names[i],"-")) {
	if (col != 0)
	  validcol=look_for_column(names[i],&hcol);
	test=fopen(names[i],"r");
	if (test == NULL) {
	  fprintf(stderr,"File %s not found!\n",names[i]);
	}
	else {
	  fclose(test);
	  if ((col != 0) && (validcol == 1))
	    *col=hcol;
	  if (col != 0) {
	    if (verbosity&VER_INPUT)
	      fprintf(stderr,"Using %s as datafile, reading column %u\n",
		      names[i],*col);
	  }
	  else {
	    if (verbosity&VER_INPUT)
	      fprintf(stderr,"Using %s as datafile!\n",names[i]);
	  }
	  check_alloc(retname=(char*)calloc(strlen(names[i])+1,(size_t)1));
	  strcpy(retname,names[i]);
	  names[i]=NULL;
	  return retname;
	}
      }
      else {
	valid=1;
	break;
      }
    }
  }

  if (valid == 1) {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Reading input from stdin!\n");
    return NULL;
  }
  
  if (verbosity&VER_INPUT) {
    if ((col != 0) && (validcol == 1))
      fprintf(stderr,"Reading input from stdin, using column %u!\n",*col);
    else
      fprintf(stderr,"Reading input from stdin!\n");
  }

  return NULL;
}
