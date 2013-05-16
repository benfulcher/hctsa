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
/* Author: Rainer Hegger Last modified: Aug 14th, 1998 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tisean_cec.h"

void solvele(double **mat,double *vec,unsigned int n)
{
  double vswap,*mswap,*hvec,max,h,pivot,q;
  int i,j,k,maxi;

  for (i=0;i<n-1;i++) {
    max=fabs(mat[i][i]);
    maxi=i;
    for (j=i+1;j<n;j++)
      if ((h=fabs(mat[j][i])) > max) {
	max=h;
	maxi=j;
      }
    if (maxi != i) {
      mswap=mat[i];
      mat[i]=mat[maxi];
      mat[maxi]=mswap;
      vswap=vec[i];
      vec[i]=vec[maxi];
      vec[maxi]=vswap;
    }
    
    hvec=mat[i];
    pivot=hvec[i];
    if (fabs(pivot) == 0.0) {
      fprintf(stderr,"Singular matrix! Exiting!\n");
      exit(SOLVELE_SINGULAR_MATRIX);
    }
    for (j=i+1;j<n;j++) {
      q= -mat[j][i]/pivot;
      mat[j][i]=0.0;
      for (k=i+1;k<n;k++)
	mat[j][k] += q*hvec[k];
      vec[j] += q*vec[i];
    }
  }
  vec[n-1] /= mat[n-1][n-1];
  for (i=n-2;i>=0;i--) {
    hvec=mat[i];
    for (j=n-1;j>i;j--)
      vec[i] -= hvec[j]*vec[j];
    vec[i] /= hvec[i];
  }
}
