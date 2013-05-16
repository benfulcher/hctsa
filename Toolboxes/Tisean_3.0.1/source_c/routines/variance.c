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
/*Author: Rainer Hegger Last modified: May 23th, 1998 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tisean_cec.h"

void variance(double *s,unsigned long l,double *av,double *var)
{
  unsigned long i;
  double h;
  
  *av= *var=0.0;

  for (i=0;i<l;i++) {
    h=s[i];
    *av += h;
    *var += h*h;
  }
  *av /= (double)l;
  *var=sqrt(fabs((*var)/(double)l-(*av)*(*av)));
  if (*var == 0.0) {
    fprintf(stderr,"Variance of the data is zero. Exiting!\n\n");
    exit(VARIANCE_VAR_EQ_ZERO);
  }
}

