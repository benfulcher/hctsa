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
/*Author: Rainer Hegger Last modified: March 1st, 1998 */
#include <math.h>

unsigned long find_neighbors(double *s,long **box,long *list,double *x,
			     unsigned long l,unsigned int bs,unsigned int dim,
			     unsigned int del,double eps,unsigned long *flist)
{
  unsigned long nf=0;
  int i,i1,i2,j,j1,k,k1;
  int ib=bs-1;
  long element;
  double dx;
  
  k=(int)((dim-1)*del);
  i=(int)(x[-k]/eps)&ib;
  j=(int)(x[0]/eps)&ib;
  
  for (i1=i-1;i1<=i+1;i1++) {
    i2=i1&ib;
    for (j1=j-1;j1<=j+1;j1++) {
      element=box[i2][j1&ib];
      while (element != -1) {
	for (k=0;k<dim;k++) {
	  k1= -k*(int)del;
	  dx=fabs(x[k1]-s[element+k1]);
	  if (dx > eps)
	    break;
	}
	if (k == dim)
	  flist[nf++]=element;
	element=list[element];
      }
    }
  }
  return nf;
}
