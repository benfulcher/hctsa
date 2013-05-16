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
/*Author: Rainer Hegger Last modified: Sep 4, 1999 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tsa.h"

char* myfgets(char *str,int *size,FILE *fin,unsigned int verbosity)
{
  char *ret;
  char *hstr=NULL;
  char last;

  ret=fgets(str,*size,fin);
  if (ret == NULL)
    return NULL;

  last=str[strlen(str)-1];

  while (last != '\n') {
    *size += INPUT_SIZE;
    check_alloc(hstr=(char*)calloc((size_t)INPUT_SIZE,(size_t)1));
    check_alloc(str=realloc(str,(size_t)*size));
    ret=fgets(hstr,INPUT_SIZE,fin);
    strcat(str,hstr);
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Line in file too long. Increasing input size\n");
    last=str[strlen(str)-1];
    free(hstr);
  }

  if (ret == NULL)
    return NULL;
  else
    return str;
}
