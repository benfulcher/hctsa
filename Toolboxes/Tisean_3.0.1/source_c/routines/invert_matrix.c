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
/* Author: Rainer Hegger Last modified: Sep 5, 2004*/
/* Changes: 
 * Sep 5, 2004: added the extern check_alloc line
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern void check_alloc(void*);

double **invert_matrix(double **mat,unsigned int size)
{
  int i,j,k;
  double **hmat,**imat,*vec;
  extern void solvele(double**,double*,unsigned int);

  check_alloc(hmat=(double**)malloc(sizeof(double*)*size));
  for (i=0;i<size;i++) {
    check_alloc(hmat[i]=(double*)malloc(sizeof(double)*size));
  }

  check_alloc(imat=(double**)malloc(sizeof(double*)*size));
  for (i=0;i<size;i++) {
    check_alloc(imat[i]=(double*)malloc(sizeof(double)*size));
  }

  check_alloc(vec=(double*)malloc(sizeof(double)*size));
  
  for (i=0;i<size;i++) {
    for (j=0;j<size;j++) {
      vec[j]=(i==j)?1.0:0.0;
      for (k=0;k<size;k++)
	hmat[j][k]=mat[j][k];
    }
    solvele(hmat,vec,size);
    for (j=0;j<size;j++)
      imat[j][i]=vec[j];
  }
  
  free(vec);
  for (i=0;i<size;i++)
    free(hmat[i]);
  free(hmat);

  return imat;
}
